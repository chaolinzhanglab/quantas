#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Carp;
use Data::Dumper;
use File::Basename;

use MyConfig;

my $prog = basename ($0);
my $progDir = dirname ($0);

my $verbose = 0;

my $donorBedFile = "";
my $acceptorBedFile = "";
my $inBedFile = "";

my $ext = 8;
my $big = 0;

my $cache = getDefaultCache ($prog);
my $keepCache = 0;


GetOptions (
	"big"=>\$big,
	"map:s"=>\$inBedFile,
	"5:s"=>\$donorBedFile,
	"3:s"=>\$acceptorBedFile,
	"ext:i"=>\$ext,
	"c:s"=>\$cache,
	"keep-cache"=>\$keepCache,
	"v"=>\$verbose);


if (@ARGV != 2)
{
	print "get unspliced reads from sam file\n";
	print "Usage: $prog [options] <in.sam or in.bam> <out.sam>\n";
	print " -map  [string]: bed file of mapped reads\n";
	print " -big          : the bed file of mapped reads is big\n";
	print " -5    [string]: bed file of 5' splice sites\n";
	print " -3    [string]: bed file of 3' splice sites\n";
	print " -ext  [int]   : size on each side of the splice site required to cover ($ext)\n";
	print " -c    [string]: cache dir ($cache)\n";
	print " --keep-cache  : keep cache\n";
	print " -v            : verbose\n";
	exit (1);
} 

my ($inSamFile, $outSamFile) = @ARGV;

my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir $cache\n" if $ret != 0;


print "prepare splice sites ...\n" if $verbose;

my $spliceSiteBedFile = "$cache/ss.ext.bed";
my $verboseFlag = $verbose ? "-v" : "";

my $cmd = "perl $progDir/bedExt.pl $verboseFlag -l \"-$ext\" -r " . ($ext-1) ." $donorBedFile - > $spliceSiteBedFile";
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

$cmd = "perl $progDir/bedExt.pl $verboseFlag -l \"-" . ($ext-1) . "\" -r $ext $acceptorBedFile - >> $spliceSiteBedFile";
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;


print "get unspliced reads ...\n" if $verbose;

my $bigFlag = $big ? "-big" : "";
my $unsplicedBedFile = "$cache/unspliced.bed";
$cmd = "perl $progDir/tagoverlap.pl $bigFlag -c $cache/tagoverlap -region $spliceSiteBedFile -denom region --complete-overlap -d \"#\" --keep-score -v $inBedFile $unsplicedBedFile";
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

my $unsplicedIdPair = "$cache/unspliced.idpair";

$cmd = "cut -f 4 $unsplicedBedFile | awk -F \"#\" '{print \$1\"\\t\"\$2}' > $unsplicedIdPair";
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

my $tmpSamFile = "$cache/out.sam";

$cmd = "perl $progDir/removeRow.pl -r $inSamFile $unsplicedIdPair | perl $progDir/selectRow.pl - $unsplicedIdPair > $tmpSamFile";
if ($inSamFile =~/\.bam$/)
{
	$cmd = "samtools view $inSamFile | perl $progDir/removeRow.pl -r - $unsplicedIdPair | perl $progDir/selectRow.pl - $unsplicedIdPair > $tmpSamFile";
}
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

my $unsplicedId = "$cache/unspliced.id";
$cmd = "awk '{print \"SS:Z:\"\$2}' $unsplicedIdPair > $unsplicedId";
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

$cmd = "paste $tmpSamFile $unsplicedId > $outSamFile";
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

system ("rm -rf $cache") unless $keepCache;



