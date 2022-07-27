#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Math::CDF qw(:all);

use Vcf;
use Common;


my $prog = basename ($0);
my $verbose = 0;

my $type = "sense"; #antisense
my $minDepth = 5; # >=20
my $minStrandBias = 0.95; #>=10

my $printAll = 0;

my $progDir = dirname ($0);

GetOptions (
	"d:i"=>\$minDepth,
	"s:f"=>\$minStrandBias,
	"all"=>\$printAll,
	"type:s"=>\$type,
	"v"=>\$verbose);

if (@ARGV != 1)
{
	print "determine strand by read count for strand-specific libraries\n";
	print "Usage: $prog [options] <in.vcf>\n";
	print " <in.vcf> : use - for STDIN\n";
    print " -d     [int]    : minimum depth to report a potential SNV ($minDepth)\n";
    print " -s     [float]  : minimum strand bias to call strand ($minStrandBias)\n";
	print " -all            : print all sites without filtering\n";
	print " -type  [string] : library type ([sense]|antisense)\n";
	print " -v              : verbose\n";
	exit (1);
}


my ($inVcfFile) = @ARGV;

my $fin;

if ($inVcfFile eq '-')
{
	$fin = *STDIN;
}
else
{
	open ($fin, "<$inVcfFile") || Carp::croak "cannot open file $inVcfFile to read\n";
}

my $iter = 0;

print generateVcfHeader (), "\n";

while (my $line =<$fin>)
{
	chomp $line;
	next if $line =~/^\#/;

	print STDERR "$iter ...\n" if $verbose && $iter % 10000 == 0;
	$iter++;
	my $snv = lineToVcf ($line);

	my $snvInfo = $snv->{'info'};
	Carp::croak "no DP4 info for SNV:", Dumper ($snv), "\n" unless exists $snvInfo->{'DP4'};
	my $dp4 = $snvInfo->{'DP4'};

	my ($refPos, $refNeg, $altPos, $altNeg) = split (/\,/, $dp4);	
	
	my $pos = $refPos + $altPos;
	my $neg = $refNeg + $altNeg;

	my $total = $pos + $neg;

	my $strand = ".";
	my $totalSense = $pos > $neg ? $pos : $neg;	

	if ($totalSense < $minDepth)
	{	
		$strand = ".";
	}
	elsif ($pos / $total >= $minStrandBias)	
	{
		$strand = $type eq 'sense' ? '+' : '-';
	}
	elsif ($neg / $total >= $minStrandBias)
	{
		$strand = $type eq 'sense' ? '-' : '+';
	}
	else
	{
		$strand =  '.';
	}

	if ($strand eq '.')
	{
		next unless $printAll;
	}
	
	$snvInfo->{'strand'} = $strand;
	print vcfToLine ($snv), "\n";
}

close ($fin) unless $inVcfFile eq '-';



