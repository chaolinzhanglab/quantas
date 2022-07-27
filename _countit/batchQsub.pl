#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Carp;
use File::Basename;

use MyConfig;
use SGE;


my $prog = basename ($0);


my $cache = MyConfig::getDefaultCache ($prog); 

my $jobName = $ENV{'LOGNAME'};
my $queueFile = "";
my $dir = "";
my $qsubOK = 0;
my $memory = "";
my $time = "";
my $cores = 1;

my $noSGE = 0;
my $waitForSGE = 0;
my $verbose = 0;

GetOptions ('c|cache:s'=>\$cache,
		'j|jobname:s'=>\$jobName,
		'm|memory:s'=>\$memory,
		't|time:s'=>\$time,
		'p:i'=>\$cores,
		'q|queues:s'=>\$queueFile,
		'd|dir:s'=>\$dir,
		'wait'=>\$waitForSGE,
		'no-sge'=>\$noSGE,
		'v'=>\$verbose);


if (@ARGV != 1)
{
	print "submit jobs via qsub\n";
	print "Usage: $prog [options] <script.list>\n";
	print " -c [string]: cache dir for qsub output [$cache]\n";
	print " -j [string]: job names [$jobName]\n";
	print " -p    [int]: number of cores for parallel execution ($cores)\n";
	print " -m [string]: memory (e.g. 2G, default is not to specify)\n";
	print " -t [string]: running time (e.g. 2:: for 2 hours, and :20: for 20 min)\n";
	print " -d [string]: dir to be attached to each script [$dir]\n";
	print " -q [string]: a list of queue to use or file name specifying queques to use\n";
	print " --wait     : wait until the jobs are finished before exiting\n";
	print " --no-sge   : do not submit jobs to SGE; run locally\n";
	print " -v         : verbose\n";
	exit(1);
}


die "$cache already exist\n" if (-d $cache);
die "$cache already exist\n" if (-f $cache);

my $queueStr = "";

if ($noSGE == 0)
{
	my $testQsub = `which qsub`;
	chomp $testQsub;
	$qsubOK = 1 if ($testQsub =~/\/qsub$/);
}

system ("mkdir $cache") if $qsubOK;


if (-f $queueFile)
{
	print "reading queue information from $queueFile ...\n" if $verbose;
	my $fin;
	open ($fin, "<$queueFile") || Carp::croak "can not open file $queueFile to read\n";
	my @queues =<$fin>;
	close ($fin);

	chomp @queues;
	$queueStr = join (",", @queues);
}
elsif ($queueFile ne '')
{
	#assuming this is a list of queues
	$queueStr = $queueFile
}

my $scriptListFile = $ARGV [0];

my $fin;
my @scriptList;
open ($fin, "<$scriptListFile") || Carp::croak "can not open file $scriptListFile to read\n";
while (my $f = <$fin>)
{
	chomp $f;
	next if $f=~/^\s*$/;
	
	if (length($dir) > 0 && (-d $dir))
	{
		$f = "$dir/$f";
	}

	push @scriptList, $f;
}

my @jobIds;
if ($qsubOK)
{
	my $batchScriptFile = "$cache/batch.sh";
	my $fout;
	open ($fout, ">$batchScriptFile") || Carp::croak "cannot open file $batchScriptFile to write\n";

	print $fout "#!/bin/bash\n";

	my $njobs = @scriptList;
	my $paramLine = "#\$ -t 1-$njobs";
	$paramLine .= " -pe smp $cores" if $cores > 1;
	$paramLine .= " -N $jobName -cwd -v TMPDIR=$cache -V -e $cache -o $cache";
	
	my $resource = "";
	$resource .= "mem=$memory," if $memory;
	$resource .= "time=$time," if $time;	
	chop $resource if $resource;

	$paramLine .= " -l \"$resource\"" if $resource;
	$paramLine .= " -q $queueStr" if $queueStr;

	print $fout $paramLine, "\n";
	
	print $fout "files=(", join (" ", @scriptList), ")\n";
	print $fout "f=\${files\[\$SGE_TASK_ID-1\]}\n";
	print $fout "bash \$f\n";
	close ($fout);

	#my $ret = `qsub $batchScriptFile`;
	#chomp $ret;
	my @lines = `qsub $batchScriptFile`;
	my $jobId = "";
	foreach my $ret (@lines)
	{
		#in case there are warning messages, ignore them
		chomp $ret;
		print "$ret", "\n" if $verbose;
		if ($ret=~/^Your job/g)
		{
			my @cols = split (/\s+/, $ret);
			$jobId = $cols[2];
		}
	}

	Carp::croak "invalid job id: $jobId\n" unless $jobId =~/\d+/;
	$jobId =~/^(\d+)\.(\d+)-(\d+):\d+$/;
	$jobId = $1;
	my $taskIdStart = $2;
	my $taskIdEnd = $3;

	push @jobIds, $jobId;

	print "job id=$jobId ($taskIdStart - $taskIdEnd) submited ...\n" if $verbose;
}
else
{
	foreach my $f (@scriptList)
	{
		print "sh $f> $f.out 2> $f.err\n";
		my $ret = system ("sh $f > /dev/null");
		
		#my $ret = system ("sh $f > $f.out 2> $f.err");
		Carp::croak "jobs failed: $!\n" if $ret != 0;
	}
}
close ($fin);

if ($qsubOK && $waitForSGE && @jobIds > 0)
{
	print "job ids=", join ("\t", @jobIds), "\n" if $verbose;
	waitUntilSGEJobsDone (\@jobIds, $verbose);
}

