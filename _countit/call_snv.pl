#!/usr/bin/perl -w

use strict;
use Carp;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Math::CDF qw(:all);

use Vcf;
use Common;


my $prog = basename ($0);
my $verbose = 0;
my $errorRate = 0.01; # 0.01, per nucleotide error rate
my $minDepth = 5; # >=20
my $minAltBase = 2; #>=2
my $minAltBaseStrand = 1; #>=1 
my $minBC = 1;

my $minRefStrandBiasP = 0.01; #>0.01
my $maxRefStrandBiasFold = 2; #<=2
my $minAltStrandBiasP = 0.01; #>0.01
my $maxAltStrandBiasFold = 2; #<=2
my $maxAltP = 0.01;	#<0.01



my $progDir = dirname ($0);

GetOptions (
	"d:i"=>\$minDepth,
	"a:i"=>\$minAltBase,
	"as:i"=>\$minAltBaseStrand,
	'bc:i'=>\$minBC,
	"p:f"=>\$maxAltP,
	"pr:f"=>\$minRefStrandBiasP,
	"fr:f"=>\$maxRefStrandBiasFold,
	"pa:f"=>\$minAltStrandBiasP,
	"fa:f"=>\$maxAltStrandBiasFold,
	"e:f"=>\$errorRate,
	"v"=>\$verbose);

if (@ARGV != 1)
{
	print "call statistically significant snv\n";
	print "Usage: $prog [options] <in.vcf>\n";
	print " By default, the information will be added to the INFO field without filtering\n";
	print " <in.vcf> : use - for STDIN\n";
	print " -e     [float]  : error rate to test alternative alleles ($errorRate)\n";
    print " -d     [int]    : minimum depth to report a potential SNV ($minDepth)\n";
    print " -a     [int]    : minimum alt base to report a potential SNV ($minAltBase)\n";
	print " -as    [int]    : minimum alt base on each strand to report a potential SNV ($minAltBaseStrand)\n";
	print " -bc    [int]    : minimum biological complexity ($minBC)\n";
	print " -p     [float]  : p value of the alt base ($maxAltP)\n";
	print " -pr    [float]  : p value of the reference allele strand bias($minRefStrandBiasP)\n";
	print " -fr    [float]  : fold change of the reference allele strand bias ($maxRefStrandBiasFold)\n";
	print " -pa    [float]  : p value of the alternative allele strand bias ($minAltStrandBiasP)\n";
	print " -fa    [float]  : fold change of the alt allele strand bias ($maxAltStrandBiasFold)\n";
	print " -v              : verbose\n";
	exit (1);
}


Carp::croak "error rate must be positive\n" unless $errorRate > 0;


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
	
	my $refSum = $refPos+$refNeg;
	my $altSum = $altPos+$altNeg;
	my $total = $refSum + $altSum;

	next unless $total >= $minDepth && $altSum >= $minAltBase && $altPos >= $minAltBaseStrand && $altNeg >= $minAltBaseStrand && $snvInfo->{'BC'} >= $minBC;

	#these are two-tail test
	my $refStrandBiasP = binomTest ($refPos, $refSum, 0.5);
	next if $refStrandBiasP < $minRefStrandBiasP && max($refPos, $refNeg) / (min($refPos, $refNeg)+0.00001) > $maxRefStrandBiasFold;

	my $altStrandBiasP = binomTest ($altPos, $altSum, 0.5);
	next if $altStrandBiasP < $minAltStrandBiasP && max($altPos, $altNeg) / (min($altPos, $altNeg)+0.00001) > $maxAltStrandBiasFold;

	my $altP = 1- pbinom ($altSum-1, $total, $errorRate); #the probability of having $altSum or more out of $total
	$altP = 1 if $altP > 1;
	$altP = 0 if $altP <= 0;
	
	next unless $altP <= $maxAltP;

	my $score = 1000;
	$score = -log($altP) / log(10) * 10 if $altP > 0;
	$score = 0 if $score <= 0;
	
	$snvInfo->{'altP'} = sprintf ("%.4g", $altP);
	$snvInfo->{'refStrandP'} = sprintf ("%.4g", $refStrandBiasP);
	$snvInfo->{'altStrandP'} = sprintf ("%.4g", $altStrandBiasP);
	$snv->{'qual'} = sprintf ("%.4g", $score);
	
	$snv->{'id'} = $snv->{'chrom'} . ":" . ($snv->{'position'}+1) if $snv->{'id'} eq '.';
	
	print vcfToLine ($snv), "\n";
}

close ($fin) unless $inVcfFile eq '-';



