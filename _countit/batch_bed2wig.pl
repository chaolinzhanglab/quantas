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
my $cmdDir = dirname ($0);
my $verbose = 0;

my $base = "";
my $suffix = "";
my $inputFormat = "bed";
my $weight = "";
my $qsub = 0;
my $queueName = "";
my $jobName = "job";

my $cache = getDefaultCache ($prog);
my $keepCache = 0;


#my $method = "mean";  #sum

GetOptions (
	"format:s"=>\$inputFormat,
	"weight"=>\$weight,
	"base:s"=>\$base,
	"suffix:s"=>\$suffix,
	"qsub"=>\$qsub,
	"queue:s"=>\$queueName,
	"job:s"=>\$jobName,
	"cache:s"=>\$cache,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate expression matrix\n";
	print "Usage $prog [options] <in.conf> <outdir>\n";
	print " --format       [string] : input format ([bed]|bedgraph), gz or bz2 files are accepted\n";
	print " --weight       [string] : weight individual tags (for bed files only; bedgraphs are always weighted)\n";
	print " --base         [string] : base dir of input data\n";
	print " --suffix       [string] : add suffix to the file names\n";
	print " --qsub                  : run jobs using sge\n";
	print " --queue        [string] : queque name\n";
	print " --job          [string] : job name\n";
	print " --cache        [string] : cache dir\n";
	print " --keep-cache            : keep cache dir\n";
	print " -v                      : verbose\n";
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

system ("mkdir $cache");

my $iter = 0;
my $scriptListFile = "$cache/script.list";
my $scriptDir = "$cache/scripts";

my $fout;
if ($qsub)
{
	open ($fout, ">$scriptListFile") || Carp::croak "cannot open $scriptListFile to write\n";
	system ("mkdir $scriptDir");	
}

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

	my $outputFile = "$outDir/$gName";
	$outputFile .= $suffix if $suffix ne "";
	
	if (@$samples == 1 && $inputFormat eq 'bedgraph')
	{
		#only one bedgraph file, we just make a copy
		my $inputFile = $samples->[0];
		$inputFile = "$base/$inputFile" if $base ne '';

		$cmd = "cp $inputFile $outputFile";
	}
	else
	{
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
		$cmd .= " | perl $cmdDir/tag2profile.pl $verboseFlag -c $cache/tag2profile_$iter $weightFlag -big -exact -of bedgraph - $outputFile";
	}

	print $cmd, "\n";
	
	if ($qsub)
	{
		my $scriptFile = "$scriptDir/script_$iter.sh";
		my $fout2;
		open ($fout2, ">$scriptFile") || Carp::croak "cannot open file $scriptFile to write\n";
		print $fout2 $cmd, "\n";
		close ($fout2);
		print $fout $scriptFile, "\n";
	}
	else
	{
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;
	}
	$iter++;
}

if ($qsub)
{
	close ($fout);

	my $cmd = "perl $cmdDir/batchQsub.pl $verboseFlag -c $cache/qsub -j $jobName -q $queueName --wait $scriptListFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);
}

system ("rm -rf $cache") unless $keepCache;

