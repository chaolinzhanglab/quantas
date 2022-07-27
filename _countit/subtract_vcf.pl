#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Carp;
use Data::Dumper;
use File::Basename;

use Vcf;

my $prog = basename ($0);

my $vcfListFile = "";
my $minAltBase = 0;
my $minDepth = 0;
my $minMaf = 0;
my $fixAltBase = 0;
my $outFile = "";

my $verbose = 0;

GetOptions (
	"l:s"=>\$vcfListFile,
	"d:i"=>\$minDepth,
	"a:i"=>\$minAltBase,
	"r:f"=>\$minMaf,
	"x"=>\$fixAltBase,
	"o:s"=>\$outFile,
	"v"=>\$verbose);

if (@ARGV !=3)
{
	print "subtract read count in vcf2 from vcf1\n";
	print "Usage: $prog [options] <in1.vcf> <in2.vcf> <out.vcf>\n";
	print " [options\n";
	print " -v         : verbose\n";
	exit (1);
}

my ($inVcfFile1, $inVcfFile2, $outVcfFile) = @ARGV;

my %snvHash;


my ($fin1, $fin2, $fout);

open ($fin1, "<$inVcfFile1") || Carp::croak "cannot open file $inVcfFile1 to read\n";
open ($fin2, "<$inVcfFile2") || Carp::croak "cannot open file $inVcfFile2 to read\n";

if ($outVcfFile eq '-')
{
	$fout = *STDOUT;
}
else
{
	open ($fout, ">$outVcfFile") || Carp::croak "cannot open file $outVcfFile to write\n";
}

print $fout join ("\t", "#CHROM", "POS ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), "\n";

my $i = 0;

while (my $line = <$fin1>)
{
	chomp $line;
	#next if $line =~/^\s*$/;
	if ($line =~/^\#/)
	{
		#we assume two files have the same format
		$line=<$fin2>;
		next;
	}

		
	print STDERR "$i ...\n" if $verbose && $i % 100000 == 0;
	$i++;

	my $snv = lineToVcf ($line);

	my $line2 = <$fin2>; chomp $line2;
	my $snv2 = lineToVcf ($line2);

	Carp::croak "data inconsistency\n" unless $snv->{'id'} eq $snv2->{'id'};


	my $chrom = $snv->{'chrom'};
	my $position = $snv->{'position'};
	
	my $dp4 = $snv->{'info'}->{'DP4'};
	my $dp42 = $snv2->{'info'}->{'DP4'};

	my ($refPos, $refNeg, $altPos, $altNeg) = split (/\,/, $dp4);
	my ($refPos2, $refNeg2, $altPos2, $altNeg2) = split (/\,/, $dp42);

	$snv->{'info'}->{'DP4'} = join(",", $refPos - $refPos2, $refNeg - $refNeg2, $altPos - $altPos2, $altNeg - $altNeg2);
	
	$snv->{'info'}{'BC'} -= $snv2->{'info'}{'BC'} if exists $snv->{'info'}{'BC'} && exists $snv2->{'info'}{'BC'};
	
	print $fout vcfToLine ($snv), "\n";
}

close ($fin1);
close ($fin2);
close ($fout) if $outVcfFile ne '-';



