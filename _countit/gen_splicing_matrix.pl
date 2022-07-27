#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;

use Quantas;

my $prog = basename ($0);
my $verbose = 0;

my $type = 'cass';
my $base = "";
my $suffix = "";

my $minCoverage = 10;
my $maxStd = 0.1;
my $naString = "";
my $average = 0;

my $id2gene2symbolFile = "";
my $printInfo = 0;

GetOptions ("t|type:s"=>\$type,
	"base:s"=>\$base,
	"suffix:s"=>\$suffix,
	"avg"=>\$average,
	"min-cov:i"=>\$minCoverage,
	"max-std:f"=>\$maxStd,
	"na-string:s"=>\$naString,
	"id2gene2symbol:s"=>\$id2gene2symbolFile,
	"print-info"=>\$printInfo,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate AS matrix\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " <in.conf> [string]: the first column is the dir or file name, and the second column is the group name\n";
	print " -base         [string] : base dir of input data\n";
	print " -type         [string] : AS type ([cass]|taca|alt5|alt3|mutx|iret|apa|snv|snv2ss|ss)\n";
	print " -suffix       [string] : suffix of output file to be appended to group name (default=none)\n";
	print " --avg                  : use average instead of sum\n";
	print " --min-cov     [int]    : min coverage ($minCoverage)\n";
	print " --max-std     [float]  : max standard deviation ($maxStd)\n";
	print " --na-string   [string] : na string (default:empty)\n";
	print " --id2gene2symbol [file]: mapping file of id to gene to symbol\n";
	print " --print-info           : print AS information columns\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outFile) = @ARGV;

if ($base ne '')
{
	Carp::croak "dir $base does not exist\n" unless -d $base;
}

print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readConfigFile ($configFile, $base, $type, $suffix);

print "done.\n" if $verbose;


print "loading mapping file of id to gene to symbol...\n" if $verbose;
my %id2gene2symbolHash;
if (-f $id2gene2symbolFile)
{
	my $fin;
	open ($fin, "<$id2gene2symbolFile") || Carp::croak "cannot open file $id2gene2symbolFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my ($id, $geneId, $symbol) = split (/\t/, $line);
		$id2gene2symbolHash{$id} = "$geneId//$symbol";
	}	

	close ($fin);
}
elsif ($id2gene2symbolFile ne '')
{
	Carp::croak "cannot open file $id2gene2symbolFile to read\n";
}


my $n = keys %id2gene2symbolHash;

print "$n mapping entries loaded\n" if $verbose;



print "loading and aggregating data of individual samples ...\n" if $verbose;

#my %sampleData;
my $ASInfo;
my $nAS = 0;
my $iter = 0;

my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;
my @groupData;

for (my $g = 0; $g < @groupNames; $g++)
{
	my $gName = $groupNames[$g];
	my $samples = $groups->{$gName}->{"samples"};
	my $nsamples = @$samples;

	foreach my $s (@$samples)
	{
		print "$iter: group=$gName, sample=$s\n" if $verbose;
		my $inputFile = $base ne '' ? "$base/$s" : $s;
        if (-d $inputFile)
        {
            $inputFile = "$inputFile/$type.count.txt";
        }

		my $sdata = readASDataFile ($inputFile, $type);
		$ASInfo = $sdata->{"ASInfo"};
		if ($nAS != 0)
		{
			Carp::croak "data inconsistency detected\n" if @$ASInfo != $nAS;
		}
		else
		{
			$nAS = @$ASInfo;
		}
		#$sampleData{$s} = $sdata->{"data"};
	
		my $data = $sdata->{'data'};
		for (my $i = 0; $i < $nAS; $i++)
        {
            my $d = $data->[$i];
            for (my $j = 0; $j < @$d; $j++)
            {
                $groupData[$g][$i][$j] += $d->[$j];
                $groupData[$g][$i][$j] /= $nsamples if $average;
            }
        }			

		$iter++;
	}
}

print "$iter samples, $nAS events loaded.\n" if $verbose;



my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

my $header = "";
if ($printInfo)
{
	$header = join ("\t", "#chrom", "chromStart", "chromEnd", "event_id", "score", "strand", "type", "isoforms");
}
else
{
	$header = "#event_id";
}

if (-f $id2gene2symbolFile)
{
	$header .= "\tNAME";
}


print $fout join ("\t", $header, @groupNames), "\n";

for (my $i = 0; $i < $nAS; $i++)
{
	my @out;
	for (my $g = 0; $g < @groupNames; $g++)
	{
		my $d = $groupData[$g][$i];

		my ($in, $ex);
		if ($type eq 'cass')
		{
			$in = $d->[3] + $d->[4];
			$ex = $d->[5]*2;
		}
		elsif ($type eq 'alt3' || $type eq 'alt5')
		{
			$in = $d->[3];
			$ex = $d->[4];
		}
		elsif ($type eq 'iret')
		{
			$in = $d->[0];
			$ex = $d->[1] * 2;
		}
		elsif ($type eq 'mutx')
		{
			$in = $d->[3] + $d->[4];
			$ex = $d->[6] + $d->[7];
		}
		elsif ($type eq 'taca')
		{
			$in = $d->[3];
			my $asId = $ASInfo->[$i][3];
			my @cols = split ("-", $asId);
        	my $nAltExon = $cols[2] eq 'sr' ? $cols[1] : $cols[2];
			
			$ex = $d->[4] * ($nAltExon+1);
		}
		elsif ($type eq 'apa' || $type eq 'snv' || $type eq 'snv2ss' || $type eq 'ss')
		{
			$in = $d->[0]; #site 1
			$ex = $d->[1]; #site 2
		}
		else
		{
			Carp::croak "incorrect AS type: $type\n";
		}

		my $n = $in + $ex;

		if ($n < $minCoverage)
		{
			$out[$g] = $naString;
			next;
		}

		my $p = $in / $n;
		my $std = sqrt ($p * (1-$p) / $n);
		if ($std > $maxStd)
		{
			$out[$g] = $naString;
		}
		else
		{
			$out[$g] = $p;
		}
	}

	my $gene2symbol = exists $id2gene2symbolHash{$ASInfo->[$i][3]} ? $id2gene2symbolHash{$ASInfo->[$i][3]} : "NA//NA";

	if (-f $id2gene2symbolFile)
	{
		if ($printInfo)
		{
			print $fout join ("\t", @{$ASInfo->[$i]}, $gene2symbol, @out), "\n";
		}
		else
		{
			print $fout join ("\t", $ASInfo->[$i][3], $gene2symbol, @out), "\n";
		}
	}
	else
	{
		if ($printInfo)
		{
			print $fout join ("\t", @{$ASInfo->[$i]}, @out), "\n";
		}
		else
		{
			print $fout join ("\t", $ASInfo->[$i][3], @out), "\n";
		}
	}
}


close ($fout);




