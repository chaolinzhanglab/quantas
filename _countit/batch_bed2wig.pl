#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;
use Quantas;

my $prog = basename ($0);
my $cmdDir = dirname ($0);
my $verbose = 0;

my $base = "";
my $suffix = "";
my $inputFormat = "bed";
my $weight = "";
my $cache = "";

#my $method = "mean";  #sum

GetOptions (
	"format:s"=>\$inputFormat,
	"weight"=>\$weight,
	"base:s"=>\$base,
	"suffix:s"=>\$suffix,
	"cache:s"=>\$cache,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate expression matrix\n";
	print "Usage $prog [options] <in.conf> <outdir>\n";
	print " -format       [string] : input format ([bed]|bedgraph), gz or bz2 files are accepted\n";
	print " -weight       [string] : weight individual tags (for bed files only; bedgraphs are always weighted)\n";
	print " -base         [string] : base dir of input data\n";
	print " -suffix       [string] : add suffix to the file names\n";
	print " -cache        [string] : cache dir\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outDir) = @ARGV;

my $verboseFlag = $verbose ? "-v" : "";

Carp::croak "incorrect input format: $inputFormat\n" if $inputFormat ne 'bed' && $inputFormat ne 'bedgraph';



print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readExprConfigFile ($configFile, $base, $suffix);

print "done.\n" if $verbose;


system ("mkdir $outDir");

print "generating wig files for each group...\n" if $verbose;
my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;

my $iter = 0;
foreach my $gName (@groupNames)
{
	my $samples = $groups->{$gName}->{"samples"};
	my $cmd = "cat ";

	my $f = $samples->[0];
	if ($f=~/\.gz$/i)
	{
		$cmd = "zcat ";
	}
	elsif ($f=~/\.bz2/i)
	{
		$cmd = "bzcat ";
	}

	foreach my $s (@$samples)
	{
		print "$iter: group=$gName, sample=$s\n" if $verbose;
		my $inputFile = $s;
		$inputFile = "$base/$inputFile" if $base ne '';
		
		$cmd .= " $inputFile";
	}

	if ($inputFormat eq 'bedgraph')
	{
		$cmd .= " | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\tt\\t\"\$4}'";
	}

	my $weightFlag = $inputFormat eq 'bedgraph' || $weight ? '--weight' : '';
	my $cacheFlag = $cache ne ''? "-c $cache" : '';
	$cmd .= " | perl $cmdDir/tag2profile.pl $verboseFlag $cacheFlag $weightFlag -big -exact -of bedgraph - $outDir/$gName.wig";

	print $cmd, "\n";
	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

	$iter++;
}


