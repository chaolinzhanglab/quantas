#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;

use MyConfig;
use Quantas;

my $prog = basename ($0);
my $cmdDir = dirname($0);
my $type = "";
my $verbose = 0;

my $cache = getDefaultCache ($prog);

#my $method = "mean";  #sum

GetOptions (
	"type=s"=>\$type,
	"cache:s"=>\$cache,
	"v|verbose"=>\$verbose
);

if (@ARGV < 3)
{
	print "combine summaries from technical replicates\n";
	print "Usage $prog [options] <in1.txt> <in2.txt> [...] <out.txt>\n";
	print " -type  [string]        : [expr]|bed6|bedGraph|iret|cass|alt5|alt3|mutx|taca|ss\n";
	print " -cache [string]        : cache dir ($cache)\n";
	print " -v                     : verbose\n";
	exit (1);
}

my @inFiles = @ARGV;
my $outFile = pop @inFiles;


if ($type eq 'bedGraph')
{
	my $ret = system ("mkdir $cache");
	Carp::croak "cannot create dir $cache\n" unless $ret == 0;
}


print "loading data of individual samples ...\n" if $verbose;

my $sampleNum = @inFiles;
#my @sampleData;
my $geneId;
my $rowNum = 0;
my $iter = 0;
my $info;

my $verboseFlag = $verbose ? " -v" : "";
my $tmpBedFile = "$cache/tmp.bed";

my @groupData;
my $totalTagNum = 0;

foreach my $inputFile (@inFiles)
{
	print "$iter/$sampleNum: $inputFile\n" if $verbose;
	$iter++;

	if ($type eq 'expr')
	{

		my $sdata = readExprDataFile ($inputFile);
		my $data = $sdata->{'data'};
		$info = $sdata->{"info"};
		if ($rowNum == 0)
		{
			$rowNum = @$info;
		}
		else
		{
			Carp::croak "data inconsistency detected\n" if @$info != $rowNum;
		}

		for (my $i = 0; $i < $rowNum; $i++)
		{
			my $d = $data->[$i];
			
			$groupData[$i][0] += $d->[0]; #tagNum
			$totalTagNum += $d->[0];
		}

	}
	elsif ($type eq 'iret' || $type eq 'cass' || $type eq 'alt5' || $type eq 'alt3' || $type eq 'mutx' || $type eq 'taca' || $type eq 'ss')
	{
		my $sdata = readASDataFile ($inputFile, $type);	
		my $data = $sdata->{'data'};
		$info = $sdata->{'info'};
		if ($rowNum == 0)
        {
            $rowNum = @$info;
        }
        else
        {
            Carp::croak "data inconsistency detected\n" if @$info != $rowNum;
        }

      	for (my $i = 0; $i < $rowNum; $i++)
        {
            my $d = $data->[$i];
			#alt5: isoform1Tags    isoform2Tags    altSSDistance   exonTags    proximalJunctionTags    distalJunctionTags
			#cass: isoform1Tags    isoform2Tags    exonTags    inclusionJunction1Tags  inclusionJunction2Tags  skippingJunctionTags
			#iret: isoform1Tags    isoform2Tags    retainedIntron5SSTags   retainedIntron3SSTags   junctionTags
			my $n = @$d;
			map {$groupData[$i][$_] += $d->[$_]} (0 .. ($n-1));
        }

	}
	elsif ($type eq 'bed6')
	{
		my $sdata = readBed6DataFile ($inputFile);
		my $data = $sdata->{'data'};
		$info = $sdata->{"info"};
		if ($rowNum == 0)
        {
            $rowNum = @$info;
        }
        else
        {
            Carp::croak "data inconsistency detected\n" if @$info != $rowNum;
        }

        for (my $i = 0; $i < $rowNum; $i++)
        {
            my $d = $data->[$i];
            #score
            $groupData[$i] += $d;
        }

	}
	elsif ($type eq 'bedGraph')
	{
		#for bedgraph file, we assume one track in each file and no header line
		my $cmd = "awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\tt\\t\"\$4}' $inputFile >> $tmpBedFile";
		$cmd = "zcat $inputFile | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\tt\\t\"\$4}' >> $tmpBedFile" if $inputFile =~/\.gz$/;
		print $cmd, "\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed:$?\n" unless $ret == 0;
	}
	else
    {
        Carp::croak "incorrec type:$type\n";
    }
}

print "$iter samples, $rowNum rows loaded.\n" if $verbose;



if ($type eq 'expr')
{

	print "recalculating RPKM\n" if $verbose;

	#recalculate RPKM by adding pseudo count
	for (my $i = 0; $i < @$info; $i++)
	{
		my $exonLen = $info->[$i][2];
		my $tagNum = $groupData[$i][0];
		
		$tagNum = 1 if $tagNum == 0; #pseudo count
		$groupData[$i][1] = $tagNum * 1e9 / $exonLen / $totalTagNum;
	}
}
elsif ($type eq 'bedGraph')
{
	print "combine samples ...\n" if $verbose;
	
	my $cmd = "perl $cmdDir/tag2profile.pl $verboseFlag -big -weight -exact -of bedgraph $tmpBedFile $outFile";
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" unless $ret == 0;
	system("rm -rf $cache");
}

if ($type ne 'bedGraph')
{
	print "write output to $outFile...\n" if $verbose;
	my $fout;
	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

	#print header line
	if ($type eq 'expr')
	{
		print $fout join ("\t", "#gene_id", "gene_symbol", "tag_num", "exon_len", "RPKM"), "\n";
	}
	elsif ($type eq 'cass')
	{
		print $fout join ("\t", "#chrom", qw(chromStart  chromEnd    name    score   strand  type    isoformIDs  isoform1Tags    isoform2Tags exonTags    inclusionJunction1Tags  inclusionJunction2Tags  skippingJunctionTags)), "\n";
	}
	elsif ($type eq 'alt5' || $type eq 'alt3')
	{
		print $fout join ("\t", "#chrom", qw(chromStart  chromEnd    name    score   strand  type    isoformIDs  isoform1Tags    isoform2Tags altSSDistance   exonTags    proximalJunctionTags    distalJunctionTags)), "\n";
	}
	elsif ($type eq 'mutx')
	{
		print $fout join ("\t", "#chrom", qw(chromStart  chromEnd    name    score   strand  type    isoformIDs  isoform1Tags    isoform2Tags 5'ExonTags  5'ExonJunction1Tags 5'ExonJunction2Tags 3'ExonTags  3'ExonJunction1Tags 3'ExonJunction2Tags)), "\n";
	}
	elsif ($type eq 'taca')
	{
		print $fout join ("\t", "#chrom", qw(chromStart  chromEnd    name    score   strand  type    isoformIDs  isoform1Tags    isoform2Tags exonTags    inclusionJunctionTags   skippingJunctionTags)), "\n";
	}
	elsif ($type eq 'iret')
	{
		print $fout join ("\t", "#chrom", qw(chromStart  chromEnd    name    score   strand  type    isoformIDs  isoform1Tags    isoform2Tags    retainedIntron5SSTags   retainedIntron3SSTags   junctionTags)), "\n";
	}
	elsif ($type eq 'ss')
	{
		print $fout join ("\t", "#chrom", qw(chromStart  chromEnd    name    score   strand  type    isoformIDs  isoform1Tags    isoform2Tags)), "\n";
	}
	elsif ($type eq 'bed6')
	{
		#no header
	}
	else
	{
		Carp::croak "incorrec type:$type\n";
	}

	#print data
	for (my $i = 0; $i < $rowNum; $i++)
	{
		if ($type eq 'expr')
		{
			my ($geneId, $symbol, $exonLen) = @{$info->[$i]};
			my ($tagNum, $RPKM) = @{$groupData[$i]};
	
			print $fout join ("\t", $geneId, $symbol, $tagNum, $exonLen, $RPKM), "\n";
		}
		elsif ($type eq 'iret' || $type eq 'cass' || $type eq 'mutx' || $type eq 'taca' || $type eq 'ss')
		{
			print $fout join ("\t", @{$info->[$i]}, @{$groupData[$i]}), "\n";
		}
		elsif ($type eq 'alt5' || $type eq 'alt3')
		{
			my @i = @{$info->[$i]};
			my ($isoform1Tags, $isoform2Tags, $exonTags, $proximalJunctionTags, $distalJunctionTags)  = @{$groupData[$i]};
			my $altSSDistance = pop @i;
	
			print $fout join ("\t", @i, $isoform1Tags, $isoform2Tags, $altSSDistance, $exonTags, $proximalJunctionTags, $distalJunctionTags), "\n";
		}
		elsif ($type eq 'bed6')
		{
			my ($chrom, $chromStart, $chromEnd, $name, $strand) = @{$info->[$i]};
			my $score = $groupData[$i];
	
			print $fout join ("\t", $chrom, $chromStart, $chromEnd, $name, $score, $strand), "\n";
		}
		else
		{
	   	 	Carp::croak "incorrec type:$type\n";
		}
	}
	
	close ($fout);
}



