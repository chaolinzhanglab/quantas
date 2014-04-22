#!/usr/bin/perl -w

use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

use MyConfig;

my $gene2symbolFile = "";
my $geneExtBedFile = "";

my $verbose = 0;

my $prog = basename ($0);
my $cmdDir = dirname ($0);
my $cache = getDefaultCache ($prog);

GetOptions (
		"gene-ext-file:s"=>\$geneExtBedFile,
		"gene2symbol:s"=>\$gene2symbolFile,
		"cache:s"=>\$cache,
		"v"=>\$verbose);


if (@ARGV != 3)
{
	print "assign polyA site to genes\n";
	print "$prog [options] <gene.meta.bed> <polyA.bed> <out.bed>\n";
	print " --gene-ext-file [string] : extended gene bed file\n";
	print " --gene2symbol            : gene2symbol file\n";
	print " -cache          [string] : cache dir ($cache)\n";
	print " -v                       : verbose\n";
	exit (1);
}


my ($geneBedFile, $polyABedFile, $outFile) = @ARGV;

my $ret = system("mkdir $cache");
Carp::croak "cannot create dir $cache\n" if $ret != 0;


#compare polyA site with genes

print "match polyA sites with genes ...\n" if $verbose;

my $polyA_vs_geneBedFile = "$cache/polyA_vs_gene.bed";

my $verboseFlag = $verbose ? "-v" : "";
my $cmd = "perl $cmdDir/tagoverlap.pl $verboseFlag -region $geneBedFile -ss $polyABedFile $polyA_vs_geneBedFile";
print "$cmd\n" if $verbose;

$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

my $polyA_vs_geneIdPairFile = "$cache/polyA_vs_gene.idpair";
$cmd = "cut -f 4 $polyA_vs_geneBedFile | awk -F \"//\" '{print \$1\"\\t\"\$2}' > $polyA_vs_geneIdPairFile";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

my $polyA_vs_geneIdPairSortFile = "$cache/polyA_vs_gene.idpair.sort";
$cmd = "perl $cmdDir/selectRow.pl -f 3 -s $polyA_vs_geneIdPairFile $polyABedFile > $polyA_vs_geneIdPairSortFile";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;


#compare polyA site with genes with extension in 3'UTR
my $polyA_vs_geneExtIdPairSortFile = "$cache/polyA_vs_gene.ext.idpair.sort";

if (-f $geneExtBedFile)
{
	print "match polyA sites with extended regions ...\n" if $verbose;
	my $polyA_vs_geneExtBedFile = "$cache/polyA_vs_gene.ext.bed";
	$cmd = "perl $cmdDir/tagoverlap.pl $verboseFlag -region $geneExtBedFile -ss $polyABedFile $polyA_vs_geneExtBedFile";
	print "$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

	my $polyA_vs_geneExtIdPairFile = "$cache/polyA_vs_gene.ext.idpair";
	$cmd = "cut -f 4 $polyA_vs_geneExtBedFile | awk -F \"//\" '{print \$1\"\\t\"\$2}' > $polyA_vs_geneExtIdPairFile";
	print "$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

	$cmd = "perl $cmdDir/selectRow.pl -f 3 -s $polyA_vs_geneExtIdPairFile $polyABedFile > $polyA_vs_geneExtIdPairSortFile";
	print "$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;
}
else
{
	my $cmd = "touch $polyA_vs_geneExtIdPairSortFile";
	print "$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;
}

#combine matched gene id: if no overlapping genes, find that withing ext 1k

print "combine matched gene ids ...\n" if $verbose;
$cmd = "awk '{print \$0\"\\t1\"}' $polyA_vs_geneIdPairSortFile > $cache/tmp1";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

$cmd = "awk '{print \$0\"\\t0\"}' $polyA_vs_geneExtIdPairSortFile > $cache/tmp2";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

$cmd = "cat $cache/tmp1 $cache/tmp2 > $cache/tmp";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

$cmd = "perl ~/scripts/uniqRow.pl -c max_num -id 0 -value 2 -v $cache/tmp $cache/tmp3";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

my $tmpOutFile = (-f $gene2symbolFile) ? "$cache/out" : $outFile;
$cmd = "cut -f 1-2 $cache/tmp3 | grep -v \"||\" > $tmpOutFile";
print "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;


if (-f $gene2symbolFile)
{

	print "get gene-to-symbol mapping ...\n" if $verbose;
	my $gene2symbolSortFile = "$cache/gene2symbol.sort";
	
	my $cmd = "perl ~/scripts/selectRow.pl -f 1 -p -pt \"\" -s $gene2symbolFile $tmpOutFile > $gene2symbolSortFile";
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;

	$cmd = "paste $tmpOutFile $gene2symbolSortFile | awk '{print \$1\"\\t\"\$2\"//\"\$4}' > $outFile";
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;
}

system ("rm -rf $cache");




