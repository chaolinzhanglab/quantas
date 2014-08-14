#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;


my $prog = basename ($0);
my $verbose = 0;
my $log2 = 0;

my $base = "";
#my $method = "mean";  #sum

GetOptions (
	"base:s"=>\$base,
	"log2"=>\$log2,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate expression matrix\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " -base         [string] : base dir of input data\n";
	print " -log2                  : report log2 transformed RPKM\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outFile) = @ARGV;

print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readConfigFile ($configFile, $base);

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
		
		my $sdata = readExprDataFile ($inputFile);
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

print $fout join ("\t", "#gene_id", "gene_symbol", "exon_len", @groupNames), "\n";

for (my $i = 0; $i < $n; $i++)
{
	my @out;
	for (my $g = 0; $g < @groupNames; $g++)
	{
		my $gName = $groupNames[$g];
		my $d = $groupData[$g][$i];
		my $samples = $groups->{$gName}->{"samples"};
		my $k = @$samples;	
		$out[$g] = $log2 ? log ($d->[1]/$k) / log(2) : $d->[1]/$k;
	}

	print $fout join ("\t", @{$geneInfo->[$i]}, @out), "\n";
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
	}
	close ($fin);
	return \%groups;
}

sub readExprDataFile
{
	my ($inputFile) = @_;
	
	my $fin;
	my @data;
	my @geneInfo;
	open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";
	
	my $totalTagNum = 0;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		
		my ($geneId, $symbol, $tagNum, $exonLen, $RPKM) = split (/\t/, $line);
		
		push @geneInfo, [$geneId, $symbol, $exonLen];
		push @data, [$tagNum,$RPKM];
		
		$totalTagNum += $tagNum;
	}
	close ($fin);

	#recalculate RPKM by adding pseudo count
	for (my $i = 0; $i < @geneInfo; $i++)
	{
		my $exonLen = $geneInfo[$i][2];
		my $tagNum = $data[$i][0];
		
		$tagNum = 1 if $tagNum == 0;
		
		$data[$i][1] = $tagNum * 1e9 / $exonLen / $totalTagNum;
	}

	return {geneInfo=>\@geneInfo, data=>\@data};
}


