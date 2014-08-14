#!/usr/bin/perl -w

use strict;
use Carp;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Math::CDF qw(:all);


use Sequence;
use Vcf;
use Common;
use MyConfig;

my $prog = basename ($0);
my $verbose = 0;
my $progDir = dirname ($0);

my $knownSiteFile = "";
my $keepGeneInfo = 0;
my $keepScore = 0;
my $countStrand = "";


GetOptions (
	"l:s"=>\$knownSiteFile,
	"count-strand:s"=>\$countStrand,
	"keep-score"=>\$keepScore,
	"keep-gene"=>\$keepGeneInfo,
	"v"=>\$verbose);

if (@ARGV != 2)
{
	print "summarize snv extracted from RNASeq data\n";
	print "Usage: $prog [options] <in.vcf> <out.txt>\n";
	print " <in.vcf>                : use - for STDIN\n";
	print " -l      [string]        : a vcf file with list of known sites (optionally with strand info) to summarize\n";
	print " --count-strand [string] : only count reads on the specified strand ([+]|-)\n";
	print " --keep-score            : keep score in snv\n";
	print " --keep-gene             : keep gene information\n";
	print " -v                      : verbose\n";
	exit (1);
}


my ($inVcfFile, $outFile) = @ARGV;


if ($countStrand ne '')
{
	Carp::croak "incorrect parameter for --count-strand: $countStrand\n" unless $countStrand eq '+' || $countStrand eq '-';
}

my %knownSiteHash;
if ($knownSiteFile ne '')
{
	print STDERR "load known sites from $knownSiteFile ...\n" if $verbose;
	my $snvs = readVcfFile ($knownSiteFile, $verbose);
	my $n = @$snvs;
	print STDERR "$n sites loaded\n" if $verbose;

	my $i = 0;
	foreach my $s (@$snvs)
	{
		my $chromPos = $s->{"chrom"} . ":" . $s->{'position'};
		$s->{'info'}{'strand'} = $countStrand if $countStrand ne ''; #override the strand to count

		$knownSiteHash{$chromPos} = $s;
		$s->{'iter'} = $i++;
	}
	$snvs = [];
}


my $fin;
my $fout;
if ($inVcfFile eq '-')
{
	$fin = *STDIN;
}
else
{
	open ($fin, "<$inVcfFile") || Carp::croak "cannot open file $inVcfFile to read\n";
}

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

my $iter = 0;
my $header = join ("\t", "#chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "form", "altBaseCount", "refBaseCount", "refBase", "altBase", "refSense", "refAntisense", "altSense", "altAntisense");
$header = join ("\t", $header, "geneId", "symbol", "region") if $keepGeneInfo;

print $fout $header, "\n";


while (my $line =<$fin>)
{
	chomp $line;
	next if $line =~/^\#/;

	print "$iter ...\n" if $verbose && $iter % 10000 == 0;
	$iter++;
	my $snv = lineToVcf ($line);

	#Carp::croak Dumper ($snv), "\n";
	my $chrom = $snv->{'chrom'};
	my $position = $snv->{'position'}; #zero-based coordinate
	my $chromPos = "$chrom:$position";
	
	my $name = $snv->{'id'};
	my ($strand, $gene, $region);
	
	if ($knownSiteFile ne '')
	{
		#if a list of known sites provided, limit to those sites and use the the strand information provided in the list
		next unless exists $knownSiteHash{$chromPos};
		Carp::croak "in consistency in altBase: snv=", Dumper ($snv), "known sites=", Dumper ($knownSiteHash{$chromPos}), "\n" unless $snv->{'altBase'} eq $knownSiteHash{$chromPos}{'altBase'};
		Carp::croak "no strand information found:", Dumper ($knownSiteHash{$chromPos}), "\n" unless exists $knownSiteHash{$chromPos}{'info'}{'strand'};		

		$strand = $knownSiteHash{$chromPos}{'info'}{'strand'};
		$gene = $knownSiteHash{$chromPos}{'info'}{'gene'} if exists $knownSiteHash{$chromPos}{'info'}{'gene'};
		$region = $knownSiteHash{$chromPos}{'info'}{'region'} if exists $knownSiteHash{$chromPos}{'info'}{'region'};
		$name = $knownSiteHash{$chromPos}{'id'};
	}
	else
	{
		#all sites and get the strand information there
		Carp::croak "no strand information found:", Dumper ($snv), "\n" unless exists $snv->{'info'}{'strand'};
		$strand = $snv->{'info'}{'strand'};
		$gene = $snv->{'info'}->{'gene'} if exists $snv->{'info'}->{'gene'};
		$region = $snv->{'info'}->{'region'} if exists $snv->{'info'}->{'region'};
	}

	Carp::croak "incorrect strand info for site $name\n" unless $strand=~/[\+\-]/;

	my $snvInfo = $snv->{'info'};
	Carp::croak "no DP4 info for SNV $chromPos\n" unless exists $snvInfo->{'DP4'};
	my $dp4 = $snvInfo->{'DP4'};

	my ($refSense, $refAntisense, $altSense, $altAntisense) = split (/\,/, $dp4);	
	my ($refBase, $altBase) = ($snv->{'refBase'}, $snv->{'altBase'});
	
	if ($strand eq '-' || $strand eq '.-')
	{
		#Carp::croak "chromPos=$chromPos, ref=$refBase, alt=$altBase\n";
		$refBase = complement ($refBase);
		$altBase = complement ($altBase);
		($refAntisense, $refSense, $altAntisense, $altSense) = ($refSense, $refAntisense, $altSense, $altAntisense);
	}
	
	my $refSum = $countStrand ne '' ? $refSense : $refSense+$refAntisense;
	my $altSum = $countStrand ne '' ? $altSense : $altSense+$altAntisense;
	#note strand has already been flipped above

	my $total = $refSum + $altSum;
	my $score = $keepScore ? $snv->{'qual'} : $total;
	
	my $outLine = join ("\t", $chrom, $position, $position+1, $name, $score, $strand, "SNV", "ALT/REF", $altSum, $refSum,
            $refBase, $altBase, $refSense, $refAntisense, $altSense, $altAntisense);

	if ($keepGeneInfo)
	{
		my ($geneId, $symbol) = split ("//", $gene);
		$outLine = join ("\t", $outLine, $geneId, $symbol, $region);
	}
	
	if ($knownSiteFile eq '')
	{
		print $fout $outLine, "\n";
	}
	else
	{
		$knownSiteHash{$chromPos}->{'line'} = $outLine;
	}
}

if ($knownSiteFile ne '')
{
	foreach my $site (sort {$a->{'iter'} <=> $b->{'iter'}} values %knownSiteHash)
	{
		Carp::croak "site does not exist in file:", Dumper ($site), "\n" unless exists $site->{'line'};
		print $fout $site->{'line'}, "\n";
	}
}

close ($fout);
close ($fin) unless $inVcfFile eq '-';


