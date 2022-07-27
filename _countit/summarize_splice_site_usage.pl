#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

use MyConfig;
use Bed;
use Common;


my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $verbose = 0;
my $geneBedFile = "";
my $weight = 0;
my $big = 0;
my $separateStrand = 0;
my $anchor = 1; 

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

GetOptions (
		'big'=>\$big,
		'weight'=>\$weight,
		'ss'=>\$separateStrand,
		'anchor:i'=>\$anchor,
		'gene:s'=>\$geneBedFile,
		'c|cache:s'=>\$cache,
		'keep-cache'=>\$keepCache,
		'v|verbose'=>\$verbose);

if (@ARGV != 3)
{
	print "summarize usage of each splice sites\n";
	print "Usage: $prog [options] <intron.bed> <tag.bed> <summary.txt>\n";
	print " <intron.bed> -- bed file of introns to be considered\n";
	print " <tag.bed> -- bed file of all tags, gz file allowed\n";
	print "OPTIONS:\n";
	print " --gene   [file]: a gene bed file used to filter reads that are not bound by genes\n";
	print " -big           : the tag file is big\n";
	print " -weight        : weight tags according to score\n";
	print " --ss           : consider the two strands separately\n";
	print " --anchor [int] : minimum anchor to count overlapping tags ($anchor nt)\n";
	print " -c             : cache dir ($cache)\n";
	print " --keep-cache   : keep cache files\n";
	print " -v             : verbose\n";
	exit (1);
}

Carp::croak "anchor must be at least 1\n" unless $anchor >= 1;

system ("mkdir $cache") unless -d $cache;
my ($intronBedFile, $tagBedFile, $summaryFile) = @ARGV;


my $bigFlag = $big ? '-big' : '';
my $verboseFlag = $verbose ? '-v' : '';
my $ssFlag = $separateStrand ? '--ss' : '';
my $weightFlag = $weight ? '-weight' : '';


if ($geneBedFile ne '')
{
	Carp::croak "$geneBedFile does not exist\n" unless -f $geneBedFile;
	
	my $tagBedFile2 = "$cache/tag.bed";
	my $cmd = "perl $cmdDir/tagoverlap.pl $verboseFlag $bigFlag $ssFlag -region $geneBedFile --denom tag --complete-overlap --keep-tag-name --keep-score --non-redundant $tagBedFile $tagBedFile2";

	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

	$tagBedFile = $tagBedFile2;
}

my $tagJunctionBedFile = "$cache/tag.junction.bed";

my $cmd = "grep -v \"^track\" $tagBedFile | awk '{if(NF==12 && \$10>1) {print \$0}}' > $tagJunctionBedFile";
$cmd = "gunzip -c $tagBedFile | grep -v \"^track\" | awk '{if(NF==12 && \$10>1) {print \$0}}' > $tagJunctionBedFile" if $tagBedFile =~/\.gz$/;
$cmd = "bunzip2 -c $tagBedFile | grep -v \"^track\" | awk '{if(NF==12 && \$10>1) {print \$0}}' > $tagJunctionBedFile" if $tagBedFile =~/\.bz2$/;

print $cmd, "\n" if $verbose;
my $ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

my $cleanTagBedFile = "$cache/tag.clean.bed";

my $minTagLen = 2 * $anchor - 1;

my $anchor2 = $anchor - 1; #this is because the junciton nucleotide 
$cmd = "grep -v \"^track\" $tagBedFile | awk '{if(\$3-\$2> $minTagLen) {print \$1\"\\t\"\$2+$anchor2\"\\t\"\$3-$anchor2\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}}' > $cleanTagBedFile";
$cmd = "gunzip -c $tagBedFile | grep -v \"^track\" | awk '{if(\$3-\$2>= $minTagLen) {print \$1\"\\t\"\$2+$anchor2\"\\t\"\$3-$anchor2\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}}' > $cleanTagBedFile" if $tagBedFile =~/\.gz$/;
$cmd = "bunzip2 -c $tagBedFile | grep -v \"^track\" | awk '{if(\$3-\$2>= $minTagLen) {print \$1\"\\t\"\$2+$anchor2\"\\t\"\$3-$anchor2\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}}' > $cleanTagBedFile" if $tagBedFile =~/\.bz2$/;

print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;


######################################################
#count junction tags for each splice sites
######################################################

print "loading introns ...\n" if $verbose;

my (%junctionStartHash, %junctionEndHash);
my ($n, $nStart, $nEnd) = 0;
my $fin;
open ($fin, "<$intronBedFile") || Carp::croak "cannot open file $intronBedFile to read\n";
while (my $j = readNextBedLine ($fin))
{
    my $chrom = $j->{'chrom'};
    my $strand = $j->{'strand'};

    my $jStart = $j->{'chromStart'} - 1;
    my $jEnd = $j->{'chromEnd'} + 1;

	if (not exists $junctionStartHash{$chrom}{$strand}{$jStart})
	{
		$nStart++;
    	$junctionStartHash{$chrom}{$strand}{$jStart} = 0;
	}
	
	if (not exists $junctionEndHash{$chrom}{$strand}{$jEnd})
	{
		$nEnd++;
    	$junctionEndHash{$chrom}{$strand}{$jEnd} = 0;
	}
	$n++;
}
close ($fin);

print "$n introns loaded. detected $nStart junction starts, $nEnd junciton ends\n" if $verbose;


print "count junction tags for each intron ...\n" if $verbose;

open ($fin, "<$tagJunctionBedFile") || Carp::croak "cannot open file $tagJunctionBedFile to read\n";
my $iter = 0;
while (my $t = readNextBedLine ($fin))
{
	print "$iter ...\n" if $iter % 100000 == 0 && $verbose;
    $iter++;

	next unless exists $t->{'blockCount'} && $t->{'blockCount'} >= 2;
	my $chrom = $t->{'chrom'};
	my $chromStart = $t->{'chromStart'};
	my $chromEnd = $t->{'chromEnd'};
	my $strand = $t->{'strand'};
	my $score = $weight ? $t->{'score'} : 1;
	
	for (my $i = 0; $i < $t->{'blockCount'} - 1; $i++)
	{
		my $jStart = $chromStart + $t->{'blockStarts'}[$i] + $t->{'blockSizes'}[$i] - 1;
		my $jEnd = $chromStart + $t->{'blockStarts'}[$i+1];
		my @strands = $separateStrand ? ($strand) : ('+', '-');
		
		foreach my $s (@strands)
		{
			$junctionStartHash{$chrom}{$s}{$jStart} += $score if exists $junctionStartHash{$chrom}{$s}{$jStart} && $jStart - $chromStart + 1 >= $anchor;
			$junctionEndHash{$chrom}{$s}{$jEnd} += $score if exists $junctionEndHash{$chrom}{$s}{$jEnd} && $chromEnd - $jEnd + 1 >= $anchor;
		}
	}
}
close ($fin);

print "saving junction tag count for each splice site ...\n" if $verbose;
my $siteJunctionTagCountFile = "$cache/site.junction_tagcount.bed";

my $fout;
open ($fout, ">$siteJunctionTagCountFile") || Carp::croak "cannot open file $siteJunctionTagCountFile to write\n";

foreach my $chrom (sort keys %junctionStartHash)
{
    foreach my $strand (sort keys %{$junctionStartHash{$chrom}})
    {
		#print "processing $chrom, strand $strand ...\n" if $verbose;

		my @junctionStart = sort {$a <=> $b} keys %{$junctionStartHash{$chrom}{$strand}};
		foreach my $jStart (@junctionStart)
		{
			my $score = $junctionStartHash{$chrom}{$strand}{$jStart};
			my $type = $strand eq '+' ? '5ss' : '3ss';
			my $name = "$chrom:" . ($jStart+1) . "//$strand//$type";
			print $fout join ("\t", $chrom, $jStart, $jStart + 1, $name, $score, $strand), "\n";
		}

		my @junctionEnd = sort {$a <=> $b} keys %{$junctionEndHash{$chrom}{$strand}};
		foreach my $jEnd (@junctionEnd)
		{
			my $score = $junctionEndHash{$chrom}{$strand}{$jEnd};
			my $type = $strand eq '+' ? '3ss' : '5ss';
			my $name = "$chrom:" . ($jEnd+1) . "//$strand//$type";
			print $fout join ("\t", $chrom, $jEnd, $jEnd + 1, $name, $score, $strand), "\n";
		}
	}
}
close ($fout);

##########################################
#count overlapping tags for each splice sites
##########################################

print "count number of tags overlapping with each splice site ...\n" if $verbose;

my $siteOverlapTagCountFile = "$cache/site.overlap_tagcount.bed";

$cmd = "perl $cmdDir/tag2profile.pl $verboseFlag $bigFlag $weightFlag $ssFlag -c $cache/tmp_exon_tag_count -region $siteJunctionTagCountFile $cleanTagBedFile $siteOverlapTagCountFile";
print $cmd, "\n" if $verbose;

$ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

$cmd = "echo \"#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\ttype\tisoformIDs\tisoform1Tags\tisoform2Tags\" > $summaryFile";
print $cmd, "\n" if $verbose;

$ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

$cmd = "paste $siteJunctionTagCountFile $siteOverlapTagCountFile | awk '{n=\$11;if(n<\$5) {n=\$5}; s=0;if(n>0){s=\$5/n}; print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"s\"\\t\"\$6\"\\tsplice_site\\tSPLICE/UNSPLICE\\t\"\$5\"\\t\"n-\$5}' >> $summaryFile";
print $cmd, "\n" if $verbose;

$ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;


system ("rm -rf $cache") unless $keepCache;

