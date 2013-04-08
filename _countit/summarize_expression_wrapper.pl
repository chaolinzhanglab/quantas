#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Carp;

use MyConfig;
use Common;


my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $exonBedFile = "";
my $e2gFile = "";
my $separateStrand = 0;

my $locationFile = "";
my $dbkey = "";

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

my $big = 0;
my $weight = 0;

my $verbose = 0;


my @ARGV0 = @ARGV;




#priority to search index file:  1. given directly, 2 specified by a galaxy location file and dbkey 3. build from scratch from a fasta file

#when location of genome index, dbkey, and junction bed files are specified

GetOptions (
		"big"=>\$big,
		"exon:s"=>\$exonBedFile,
		"e2g:s"=>\$e2gFile,
		"loc:s"=>\$locationFile,
		"dbkey:s"=>\$dbkey,
		"ss"=>\$separateStrand,
		'weight'=>\$weight,
		"cache:s"=>\$cache,
		'keep-cache'=>\$keepCache,
		"v"=>\$verbose);

if (@ARGV !=2)
{
	print "Usage: $prog [options] <input.bed> <out.txt>\n";
	print " -big                 : big file\n";
	print " -exon       [file]   : exon bed file\n";
	print " -e2g        [file]   : exonid-to-geneid-to-genesymbol mapping\n";
	print " -loc        [file]   : location file name\n";
   	print " -dbkey      [string] : dbkey\n";
	print " -weight              : weight each tag by the score column\n";
	print " --ss                 : consider the two strands seperately\n";
	print " --cache     [string] : cache dir ($cache)\n";
	print " --keep-cache         : do not delete cache files after the job done\n";
	print " -v                   : verbose\n";
	exit (0);
}

sub getPathFromLocationFile
{
	my ($inFile, $dbkey) = @_;
	my %pathHash;
	my $fin;
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^#/;

		my ($db, $analysis, $path, $type) = split (/\s+/, $line);
		if ($db eq $dbkey && $analysis eq 'rnaseq')
		{
			$pathHash{$type} = $path;
		}
	}
	close ($fin);
	return \%pathHash;
}

print "CMD = $prog ", join (" ", @ARGV0), "\n" if $verbose;
my ($inBedFile, $outFile) = @ARGV;


my $verboseFlag = $verbose ? "-v" : "";
my $bigFlag = $big ? '-big' : '';
my $weightFlag = $weight ? "--weight" : "";

if ($exonBedFile eq "" && $e2gFile eq "")
{
	print "location file = $locationFile\n" if $verbose;
	print "dbkey = $dbkey\n" if $verbose;
	my $path = getPathFromLocationFile ($locationFile, $dbkey);
	
	$exonBedFile = $path->{'exon'} if exists $path->{'exon'};
	$e2gFile = $path->{'id2gene2symbol'} if exists $path->{'id2gene2symbol'};
}

Carp::croak "cannot find exon bed file: $exonBedFile\n" unless -f $exonBedFile;
Carp::croak "cannot find exon-to-gene-to-symbol mapping file: $e2gFile\n" unless -f $e2gFile;

$cache = getFullPath ($cache);
system ("mkdir $cache");



print "get count of tags on each core exon\n" if $verbose;
my $exonTagCountBedFile = "$cache/exon.tag.count.bed";
my $ssFlag = $separateStrand ? '--ss' : '';
my $cmd = "perl $cmdDir/tag2profile.pl $bigFlag $ssFlag $weightFlag $verboseFlag -region $exonBedFile -of bed  $inBedFile $exonTagCountBedFile";
print "$cmd\n" if $verbose;
my $ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

print "summarize gene expression\n" if $verbose;

$cmd = "perl $cmdDir/summarize_expression.pl $verboseFlag -e2g $e2gFile $exonTagCountBedFile $outFile";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD $cmd failed: $?\n" if $ret != 0;

system ("rm -rf $cache") unless $keepCache;


