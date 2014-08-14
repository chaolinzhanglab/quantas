#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $verbose = 0;

my $separateStrand = 0;
my $isoformBedFile = "";
my $cleanIsoform = 0;
my $bigExonSize = 400;
my $outputDir = "./gapless_out";

my $big = 0;
my $printSingleton = 0;
my $usePreCalScore = 0;

my $sizeDistFile = "";

my $keepCache = 0;

GetOptions ("isoform=s"=>\$isoformBedFile,
	"clean-isoform"=>\$cleanIsoform,
	"E|big-exon-size:i"=>\$bigExonSize,
	"big"=>\$big,
	"use-pre-calc-score"=>\$usePreCalScore,
	"use-pre-calc-size-dist:s"=>\$sizeDistFile,
	"ss"=>\$separateStrand,
	"print-singleton"=>\$printSingleton,
	"o|outp-dir:s"=>\$outputDir,
	"keep-cache"=>\$keepCache,
	"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print "bayesian analysis of splicing isoform structure from paired-end mRNA-seq data\n";
	print "Usage: $prog [options] <read1.bed> <read2.bed>\n";
	print " -isoform [file]                : file name of splice isoform database\n";
	print " --use-pre-calc-score           : use pre-calculated isoform prior score\n";
	print " --clean-isoform                : clean isoforms (will supress --use-pre-calc-score)\n";
	print " -E        [int]                : big exon size (default=$bigExonSize)\n";
	print " --use-pre-calc-size-dist [file]: use pre-calculated size distribution file\n";
	print " -big                           : read number is big (i.e. over ~6M reads)\n";
	print " --ss                           : separate the two strands\n";
	print " --print-singleton              : print reads even when they are not found to be a legitimate pair\n";
	print " -o        [dir]                : output dir (default=$outputDir)\n"; 
	print " --keep-cache                   : keep cache files when the job is done\n";
	print " -v                             : verbose\n";
	exit (1);
}


my ($read1BedFile, $read2BedFile) = @ARGV;

if ($sizeDistFile ne '')
{
	Carp::croak "size distribution file does not exists\n" unless -f $sizeDistFile;
}

my $verboseFlag = $verbose ? '-v' : '';
my $bigFlag = $big ? "-big" : "";
my $separateStrandFlag = $separateStrand ? "--ss" : '';

system ("mkdir $outputDir") unless -d $outputDir;

my $tmpDir = "$outputDir/tmp";
system ("mkdir $tmpDir");

$usePreCalScore = 0 if $cleanIsoform;

my $bigExonBedFile = "$tmpDir/bigExon.bed";

if (not -f $isoformBedFile)
{
	Carp::croak "the isoform database file $isoformBedFile does not exist\n";
}

if ($cleanIsoform)
{
	print "clean isoforms ...\n" if $verbose;

	#replace the name so that they are simple and unique
	my $cleanIsoformBedFile = "$tmpDir/isoform.clean.bed";
	my $cmd = "awk 'BEGIN{i=0} {if(NF<10) {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"i\"\\t\"\$5\"\\t\"\$6} else {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"i\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12}; i=i+1}' $isoformBedFile > $cleanIsoformBedFile";
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
	#$cleanIsoformBedFile = $isoformBedFile;

	#remove redundancy
	my $cleanNrIsoformBedFile = "$tmpDir/isoform.clean.nr.bed";
	$cmd = "perl $cmdDir/bed2representative.pl $bigFlag $verboseFlag -s --same-block-count --same-intron -c $tmpDir/clean_bed_tmp_files $cleanIsoformBedFile $cleanNrIsoformBedFile";
	print $cmd, "\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

	$isoformBedFile = $cleanNrIsoformBedFile;
}

my $scoredIsoformBedFile = "$tmpDir/isoform.scored.bed";
if ($usePreCalScore)
{
	print "using pre-calculated isoform score ...\n" if $verbose;
	$scoredIsoformBedFile = $isoformBedFile;
}
else
{
	print "scoring isoforms using observed junction reads ...\n" if $verbose;
	
	print "retrieving junction reads ...\n" if $verbose;

	my $junctionRead1BedFile = "$tmpDir/read1.junction.bed";
	my $junctionRead2BedFile = "$tmpDir/read2.junction.bed";

	my $cmd = "awk '{if(\$10>1) {print \$0}}' $read1BedFile > $junctionRead1BedFile";
	my $ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

	$cmd = "awk '{if(\$10>1) {print \$0}}' $read2BedFile > $junctionRead2BedFile";
	$ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

	my $junctionReadBedFile = "$tmpDir/read.junction.bed";
	$cmd = "cat $junctionRead1BedFile $junctionRead2BedFile > $junctionReadBedFile";
	$ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

	$cmd = "perl $cmdDir/score_exon_trio.pl $verboseFlag $bigFlag -c $tmpDir/score_junction_tmp_files $isoformBedFile $junctionReadBedFile $scoredIsoformBedFile";
	print $cmd, "\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
}

if (-f $sizeDistFile)
{
	print "use provided insert size distribution file $sizeDistFile ...\n" if $verbose;
}
else
{
	print "estimating insert size distribution...\n" if $verbose;

	$sizeDistFile = "$outputDir/size_dist.txt";
	
	print "retrieving genomic reads...\n" if $verbose;

	my $genomicRead1BedFile = "$tmpDir/read1.genomic.bed";
	my $genomicRead2BedFile = "$tmpDir/read2.genomic.bed";

	my $cmd = "awk '{if(NF<=9 || \$10==1) {print \$0}}' $read1BedFile > $genomicRead1BedFile";
	my $ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

	$cmd = "awk '{if(NF<=9 || \$10==1) {print \$0}}' $read2BedFile > $genomicRead2BedFile";
	$ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

	my $bigExonBedFile = "$tmpDir/bigExon.bed";
	my $pairOnBigExonBedFile = "$tmpDir/pair.bigexon.bed";

	print "extracting exons with a size >= $bigExonSize to $bigExonBedFile...\n" if $verbose;
	
	$cmd = "awk '{if ((\$10==1 || NF <= 6) && \$3-\$2>=$bigExonSize && \$1 != \"chrM\") {print \$0}}' $isoformBedFile > $bigExonBedFile";
	$ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

	$cmd = "perl $cmdDir/combine_paired_end.pl $verboseFlag $bigFlag $separateStrandFlag -cache $tmpDir/exon_size_tmp_files -gene $bigExonBedFile $genomicRead1BedFile $genomicRead2BedFile $pairOnBigExonBedFile";
	print $cmd, "\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
	$cmd = "awk '{print \$3-\$2}' $pairOnBigExonBedFile | sort | uniq -c | awk '{print \$2\"\t\"\$1}' | sort -k 1 -n > $sizeDistFile";
	$ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
}

print "inferring structures of PE reads ...\n" if $verbose;

my $pairBedFile = "$outputDir/pair.gapless.bed";
my $printSingletonFlag = $printSingleton ? '--print-singleton' : '';
my $cmd = "perl $cmdDir/combine_paired_end.pl $verboseFlag $bigFlag $separateStrandFlag -cache $tmpDir/gapless_tmp_files -gene $scoredIsoformBedFile -use-prior -size-dist $sizeDistFile $printSingletonFlag $read1BedFile $read2BedFile $pairBedFile";
print $cmd, "\n" if $verbose;

my $ret = system ($cmd);
Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

system ("rm -rf $tmpDir") unless $keepCache;











