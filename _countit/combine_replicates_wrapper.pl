#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;

use Quantas;
use MyConfig;

my $prog = basename ($0);
my $cmdDir = dirname ($0);
my $verbose = 0;
my $cache = getDefaultCache ($prog);

my $type = 'cass';
my $base = ""; #base dir of input files
my $suffix = ""; #suffix of output files

my $qsub = 0;
my $queues = "";
my $jobName = "combine_replicates_$$"; #a unique job name

GetOptions ("t|type:s"=>\$type,
	"base:s"=>\$base,
	"suffix:s"=>\$suffix,
	"qsub"=>\$qsub,
	"q:s"=>\$queues,
	"job-name:s"=>\$jobName,
	"cache:s"=>\$cache,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "combine replicates into multiple groups defined in a configuration file\n";
	print "Usage $prog [options] <in.conf> <outdir>\n";
	print " <in.conf> [string]: the first column is the dir or file name, and the second column is the group name\n";
	print " -base         [string] : base dir of input data\n";
	print " -type         [string] : [expr]|bed6|bedGraph|iret|cass|alt5|alt3|mutx|taca|ss\n";
	print " -suffix       [string] : suffix of output file to be appended to group name (default=none)\n";
	print " -qsub                  : run jobs on cluster\n";
	print " -q            [string] : queues to use\n";
	print " -job-name     [string] : job name\n";
	print " -cache        [string] : cache ($cache)\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outDir) = @ARGV;

if ($base ne '')
{
	Carp::croak "dir $base does not exist\n" unless -d $base;
}

my $verboseFlag = $verbose ? "-v" : "";

#for qsub
my $scriptDir = "$cache/scripts";
if ($qsub)
{
	my $ret = system ("mkdir $cache");
	Carp::croak "cannot mkdir $cache\n" unless $ret == 0;
	
	$ret = system ("mkdir $scriptDir");
	Carp::croak "cannot mkdir $scriptDir\n" unless $ret == 0;
}


print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readConfigFile ($configFile, $base, $type, $suffix);

print "done.\n" if $verbose;

my $ret = system ("mkdir $outDir");
Carp::croak "cannot create dir $outDir\n" unless $ret == 0;


print "combine individual samples in each group ...\n" if $verbose;

my $iter = 0;
my $giter = 0;

my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;
foreach my $gName (@groupNames)
{
	my $samples = $groups->{$gName}->{"samples"};	
	my $outFile = "$outDir/$gName" . $suffix;
	if (@$samples < 2)
	{
		warn "only one sample in group $gName\n";
		next;
	}

	my $jobCache = $qsub ? "$cache/job_$giter" . "_cache" : $cache;
	my $cmd = "perl $cmdDir/combine_replicates.pl $verboseFlag -cache $jobCache -type $type";
	
	foreach my $s (@$samples)
	{
		print "$iter: group=$gName, sample=$s\n" if $verbose;
		my $inputFile = $base ne '' ? "$base/$s" : $s;
        if (-d $inputFile)
        {
            $inputFile = "$inputFile/$type.count.txt";
        }

		$cmd .= " $inputFile";
		$iter++;
	}
	$cmd .= " $outFile";
	print $cmd, "\n" if $verbose;
	
	if ($qsub)
	{
		my $scriptFile = "$scriptDir/script_$giter.sh";
		my $fout;
		open ($fout, ">$scriptFile") || Carp::croak "cannot open file $scriptFile to write\n";
		print $fout "#!/bin/sh\n";
		print $fout "echo \"$cmd\"\n";
		print $fout $cmd, "\n";
		close ($fout);
	}
	else
	{
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed:$?\n" unless $ret == 0;
	}

	$giter++;
}

if ($qsub)
{
	my $scriptListFile = "$cache/scripts.list";
	system ("ls -rt $scriptDir/*.sh > $scriptListFile");
	Carp::croak "The script list file at $scriptListFile was not generated properly\n" unless -f $scriptListFile;

	my $qsubCache = "$cache/qsub";
	my $queueFlag = $queues ne '' ? "-q $queues" : "";
	my $cmd = "perl $cmdDir/batchQsub.pl -v --wait $queueFlag -c $qsubCache -j $jobName $scriptListFile";

	print "CMD=$cmd\n" if $verbose;
	system ($cmd);
}


