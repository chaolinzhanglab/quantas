#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;


my $prog = basename ($0);
my $verbose = 0;

my $type = 'cass';
my $base = "";

my $minCoverage = 10;
my $maxStd = 0.1;
my $naString = "";
my $average = 0;

my $id2gene2symbolFile = "";

GetOptions ("t|type:s"=>\$type,
	"base:s"=>\$base,
	"avg"=>\$average,
	"min-cov:i"=>\$minCoverage,
	"max-std:f"=>\$maxStd,
	"na-string:s"=>\$naString,
	"id2gene2symbol:s"=>\$id2gene2symbolFile,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate AS matrix\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " <in.conf> [string]: the first column is the dir or file name, and the second column is the group name\n";
	print " -base         [string] : base dir of input data\n";
	print " -type         [string] : AS type ($type)\n";
	print " --avg                  : use average instead of sum\n";
	print " --min-cov     [int]    : min coverage ($minCoverage)\n";
	print " --max-std     [float]  : max standard deviation ($maxStd)\n";
	print " --na-string   [string] : na string (default:empty)\n";
	print " --id2gene2symbol [file]: mapping file of id to gene to symbol\n";
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
my $groups = readConfigFile ($configFile, $base);

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
my $n = keys %id2gene2symbolHash;

print "$n mapping entries loaded\n" if $verbose;



print "loading data of individual samples ...\n" if $verbose;

my %sampleData;
my $ASInfo;
my $nAS = 0;
my $iter = 0;

my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;

foreach my $gName (@groupNames)
{
	my $samples = $groups->{$gName}->{"samples"};
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
		$sampleData{$s} = $sdata->{"data"};
		$iter++;
	}
}

print "$iter samples, $nAS events loaded.\n" if $verbose;


print "aggregating samples in the same group ...\n" if $verbose;


my @groupData;

for (my $g = 0; $g < @groupNames; $g++)
{
	my $gName = $groupNames[$g];
	my $samples = $groups->{$gName}->{"samples"};

	foreach my $s (@$samples)
	{
		print "sample=$s\n" if $verbose;
		my $data = $sampleData{$s};
		for (my $i = 0; $i < $nAS; $i++)
		{
			my $d = $data->[$i];

			my $nsamples = @$d;
			for (my $j = 0; $j < @$d; $j++)
            {
                $groupData[$g][$i][$j] += $d->[$j];
				$groupData[$g][$i][$j] /= $nsamples if $average;
            }
		}
	}
}


my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

if (-f $id2gene2symbolFile)
{
	print $fout join ("\t", "#event_id", "NAME", @groupNames), "\n";
}
else
{
	print $fout join ("\t", "#event_id", @groupNames), "\n";
}


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
        	my $nAltExon = $cols[2];
			$ex = $d->[4] * ($nAltExon+1);
		}
		elsif ($type eq 'apa')
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
		print $fout join ("\t", $ASInfo->[$i][3], $gene2symbol, @out), "\n";
	}
	else
	{
		print $fout join ("\t", $ASInfo->[$i][3], @out), "\n";
	}
}


close ($fout);




sub readConfigFile
{
	my ($configFile, $base) = @_;
	my $fin;
	open ($fin, "<$configFile") || Carp::croak "cannot open file $configFile to read\n";
	my $i = 0;
	my %groups;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		next if $line=~/^\#/;
		my ($sampleName, $groupName) = split (/\t/, $line);
		$groups{$groupName}->{"id"} = $i++ unless exists $groups{$groupName};
		push @{$groups{$groupName}->{"samples"}}, $sampleName;

		my $inputFile = $base ne '' ? "$base/$sampleName" : $sampleName;
		if (-d $inputFile)
		{
			$inputFile = "$inputFile/$type.count.txt";
		}

		Carp::croak "Input file $inputFile does not exist\n" unless -f $inputFile;
	}
	close ($fin);
	return \%groups;
}

sub readASDataFile
{
    my ($inputFile, $type) = @_;

    my $fin;
    my @data;
    my @ASInfo;
    open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";
    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line =~/^\s*$/;
        next if $line =~/^\#/;

        my @cols = split (/\t/, $line);
        my (@infoCols, @dataCols);

        if ($type eq 'cass' || $type eq 'iret' || $type eq 'mutx' || $type eq 'taca')
        {
            @infoCols = @cols[0..7];
            @dataCols = @cols[8..$#cols];
        }
        elsif ($type eq 'alt3' || $type eq 'alt5')
        {
            @infoCols = @cols[0..7];
            push @infoCols, $cols[10];
            @dataCols = @cols[8..9];
            push @dataCols, @cols[11..$#cols];
        }
		elsif ($type eq 'apa')
        {
            #polyA seq data
            @infoCols = @cols[0..7];
            push @infoCols, @cols[10..$#cols];
            @dataCols = @cols[8..9];
        }
        else
        {
            Carp::croak "incorrect AS type: $type\n";
        }

        push @ASInfo, \@infoCols;
        push @data, \@dataCols;
    }
    close ($fin);

    return {ASInfo=>\@ASInfo, data=>\@data};
}


