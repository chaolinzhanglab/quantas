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
my $pseudoCount = 1;
my $log2 = 0;
my $count = 0;

my $base = "";
my $suffix = "";
#my $method = "mean";  #sum

GetOptions (
	"base:s"=>\$base,
	"suffix:s"=>\$suffix,
	"pseudocount:f"=>\$pseudoCount,
	"log2"=>\$log2,
	"raw-count"=>\$count,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate expression matrix\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " -base         [string] : base dir of input data\n";
	print " -suffix       [string] : add suffix to the file names\n";
	print " -pseudocount  [float]  : pseudocount to be used ($pseudoCount)\n";
	print " -log2                  : report log2 transformed RPKM\n";
	print " --raw-count            : report the raw read count (will suppress -log2)\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outFile) = @ARGV;

print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readExprConfigFile ($configFile, $base, $suffix);

print "done.\n" if $verbose;

print "loading data of individual samples ...\n" if $verbose;

my %sampleData;
my $geneId;
my $n = 0;
my $iter = 0;
my $geneInfo;

my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;

foreach my $gName (@groupNames)
{
	my $samples = $groups->{$gName}->{"samples"};
	foreach my $s (@$samples)
	{
		print "$iter: group=$gName, sample=$s\n" if $verbose;
		my $inputFile = $s;
		$inputFile = "$base/$inputFile" if $base ne '';
		my $sdata = readExprDataFile ($inputFile, $pseudoCount);
		$geneInfo = $sdata->{"geneInfo"};
		if ($n != 0)
		{
			Carp::croak "data inconsistency detected\n" if @$geneInfo != $n;
		}
		else
		{
			$n = @$geneInfo;
		}
		$sampleData{$s} = $sdata->{"data"};
		$iter++;
	}
}

print "$iter samples, $n genes loaded.\n" if $verbose;


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
		for (my $i = 0; $i < $n; $i++)
		{
			my $d = $data->[$i];
			
			$groupData[$g][$i][0] += $d->[0]; #tagNum
			$groupData[$g][$i][1] += $d->[1]; #RPKM
		}
	}
}


my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

#print $fout join ("\t", "#gene_id", "gene_symbol", @groupNames), "\n";
print $fout join ("\t", "gene_id", "NAME", @groupNames), "\n";

for (my $i = 0; $i < $n; $i++)
{
	my @out;
	for (my $g = 0; $g < @groupNames; $g++)
	{
		my $gName = $groupNames[$g];
		my $d = $groupData[$g][$i];
		my $samples = $groups->{$gName}->{"samples"};
		my $k = @$samples;
		if ($count)
		{
			#raw read count
			$out[$g] = $d->[0]; 
		}
		else
		{
			if ($log2 && $d->[1] == 0)
			{
				$out[$g] = "NA";
			}
			else
			{
				$out[$g] = $log2 ? log ($d->[1]/$k) / log(2) : $d->[1]/$k;
			}
		}
	}

	#print $fout join ("\t", @{$geneInfo->[$i]}, @out), "\n";
	print $fout join ("\t", $geneInfo->[$i][0], $geneInfo->[$i][1], @out), "\n";
}


close ($fout);



