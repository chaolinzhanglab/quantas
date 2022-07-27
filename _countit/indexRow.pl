#!/usr/bin/env perl
#

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

my $prog = basename($0);
my $keyIndex = -1;
my $skipLines = 0;

my $verbose = 0;

GetOptions ('key:i'=>\$keyIndex,
			'skip:i'=>\$skipLines,
			'v'=>\$verbose);

if (@ARGV != 2)
{
	print "index fasta file\n";
	print "Usage: $prog [options] <in.txt> <out.idx>\n";
	print " OPTIONS:\n";
	print " -key   [int]: column id used as key (-1 = no key, start from 0)\n";
	print " --skip [int]: number of lines to skip (default=0)\n";
	print " -v          : verbose\n";
	exit(1);
}


my ($inFile, $outIndexFile) = @ARGV;

my $fin;
open ($fin, "<$inFile") || Carp::croak "can not open file $inFile to read\n";

my $fout;
open ($fout, ">$outIndexFile") || Carp::croak "can not open file $outIndexFile to write\n";

print $fout "#file = $inFile\n";
my $pointer = 0;
my $iter = 0;
while (my $line = <$fin>)
{
	print "$iter ...\n" if $verbose && $iter % 100000 == 0;
	if ($iter >= $skipLines)
	{
		chomp $line;
		my $key = "";
	
		if ($keyIndex >= 0)
		{
			$key = "-";
			my @cols = split (/\t/, $line, $keyIndex + 2);
			$key = $cols[$keyIndex] if $keyIndex < @cols;
			print $fout join ("\t", $key, $pointer), "\n";
		}
		else
		{
			print $fout $pointer, "\n";
		}
	}
	$iter++;
	$pointer = tell ($fin);	
}
close ($fin);
close ($fout);


