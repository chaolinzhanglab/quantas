#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;


my $prog = basename ($0);
my $verbose = 0;

#my $method = "mean";  #sum

GetOptions (
	"v|verbose"=>\$verbose
);

if (@ARGV < 3)
{
	print "combine gene expression summary from technical replicates\n";
	print "Usage $prog [options] <in1.expr.txt> <in2.expr.txt> [...] <out.expr.txt>\n";
	print " -v                     : verbose\n";
	exit (1);
}

my @inFiles = @ARGV;
my $outFile = pop @inFiles;



print "loading data of individual samples ...\n" if $verbose;

my @sampleData;
my $geneId;
my $n = 0;
my $iter = 0;
my $geneInfo;

foreach my $inputFile (@inFiles)
{
	print "$iter: $inputFile\n" if $verbose;
	
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
	$sampleData[$iter] = $sdata->{"data"};
	$iter++;
}

print "$iter samples, $n genes loaded.\n" if $verbose;



print "combine samples ...\n" if $verbose;

my @groupData;
my $totalTagNum = 0;

$iter = 0;
foreach my $data (@sampleData)
{
	print "$iter ...\n" if $verbose;
	$iter++;

	for (my $i = 0; $i < $n; $i++)
	{
		my $d = $data->[$i];
			
		$groupData[$i][0] += $d->[0]; #tagNum
		#$groupData[$i][1] += 0; #RPKM
		$totalTagNum += $d->[0];
	}
}

print "recalculating RPKM\n" if $verbose;

#recalculate RPKM by adding pseudo count
for (my $i = 0; $i < @$geneInfo; $i++)
{
	my $exonLen = $geneInfo->[$i][2];
	my $tagNum = $groupData[$i][0];
		
	$tagNum = 1 if $tagNum == 0;
	$groupData[$i][1] = $tagNum * 1e9 / $exonLen / $totalTagNum;
}

my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

print $fout join ("\t", "#gene_id", "gene_symbol", "tag_num", "exon_len", "RPKM"), "\n";

for (my $i = 0; $i < $n; $i++)
{
	my ($geneId, $symbol, $exonLen) = @{$geneInfo->[$i]};
	my ($tagNum, $RPKM) = @{$groupData[$i]};

	print $fout join ("\t", $geneId, $symbol, $tagNum, $exonLen, $RPKM), "\n";
}

close ($fout);


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


