#!/usr/bin/env perl

use strict;
use warnings;


use File::Basename;
use Getopt::Long;
use Carp;


use Bed;
use MyConfig;


my $prog = basename ($0);
my $progDir = dirname ($0);
my $verbose = 0;


my $cache = MyConfig::getDefaultCache ($prog);

my $minTotal = 20;
my $minRatio = 0.2;



GetOptions ("v"=>\$verbose,
	'c:s'=>\$cache,
	'N:i'=>\$minTotal,
	'r:f'=>\$minRatio);

if (@ARGV != 3)
{
	print "define APAT events based on annotated exons and polyA sites detected in polyA seq\n";
	print " Usage: $prog [options] <site.bed> <exon.bed> <out.bed>\n";
	print " -c [string]: cache dir ($cache)\n";
	#print " -N [int]  : min total number of tags for the gene ($minTotal)\n";
	#print " -r [float]: min ratio of tags for the site among all sites in the gene ($minRatio)\n"; 
	print " -v        : verbose\n";
	exit (1);
}


my ($siteBedFile, $exonBedFile, $outBedFile) = @ARGV;

system ("mkdir $cache");

print "read polyA sites from $siteBedFile ...\n" if $verbose;

my $sites = readBedFile ($siteBedFile, $verbose);

my $n = @$sites;

print "$n sites loaded\n" if $verbose;


print "read exons from $exonBedFile ...\n" if $verbose;

my $exons = readBedFile ($exonBedFile, $verbose);

$n = @$exons;
print "$n exons loaded\n" if $verbose;


print "clean exons ...\n" if $verbose;

my %exonHash;
my %exonHashById;

foreach my $e (@$exons)
{
	my $strand = $e->{'strand'};
	my $coord = $strand eq '+' ? $e->{'chromStart'} : $e->{'chromEnd'};
	my $exonLen = $e->{'chromEnd'} - $e->{'chromStart'} + 1;
	my $exonId = $e->{'name'};

	Carp::croak "exonId not unique:$exonId\n" if exists $exonHashById {$exonId};
	$exonHashById {$exonId} = $e;

	if (exists $exonHash{$coord})
	{
		my $e0 = $exonHash{$coord};
		my $exonLen0 = $e0->{'chromEnd'} - $e0->{'chromStart'} + 1;

		if ($exonLen0 < $exonLen)
		{
			$e0->{'print'} = 0;
			$e->{'print'} = 1;
			$exonHash{$coord} = $e;
		}
	}
	else
	{
		$e->{'print'} = 1;
		$exonHash{$coord} = $e;
	}
}


my $cleanExonBedFile = "$cache/exon.clean.bed";
my $fout;

open ($fout, ">$cleanExonBedFile") || Carp::croak "cannot open file $cleanExonBedFile to write\n";
foreach my $e (@$exons)
{
	print $fout bedToLine ($e), "\n" if exists $e->{'print'} && $e->{'print'} == 1;
}
close ($fout);

my $cleanExonExtBedFile = "$cache/exon.clean.ext.bed";


my $cmd = "perl $progDir/bedExt.pl -r 1000 -v $cleanExonBedFile $cleanExonExtBedFile";
print "extend exons ...\n" if $verbose;
my $ret = system ($cmd);
Carp::croak "cmd=$cmd failed:$?\n" unless $ret == 0;


print "mapping each polyA site to exon ...\n" if $verbose;

my $siteToExonIdPairFile = "$cache/site2exon.idpair";
$cmd = "perl $progDir/polyA2gene.pl -v --gene-ext-file $cleanExonExtBedFile $cleanExonBedFile $siteBedFile $siteToExonIdPairFile";
$ret = system ($cmd);
Carp::croak "cmd=$cmd failed:$?\n" unless $ret == 0;


print "load site-2-exon mapping ...\n" if $verbose;

my $fin;
open ($fin, "<$siteToExonIdPairFile") || Carp::croak "cannot open file $siteToExonIdPairFile to read\n";

my %siteToGeneHash;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	my ($siteId, $exonId) = split ("\t", $line);
	
	$siteToGeneHash{$siteId} = $exonId;
}

print "reorganizing polyA sites...\n";

my %siteHash;
my $i = 0;
foreach my $s (@$sites)
{
	print "$i ...\n" if $verbose && $i % 10000 == 0;
	$i++;
	my $siteId = $s->{"name"};
	next unless exists $siteToGeneHash{$siteId};
	my $exonId = $siteToGeneHash{$siteId};
	push @{$siteHash{$exonId}->{'sites'}}, $s;
	$siteHash{$exonId}->{'N'} += $s->{'score'};
}

$n = keys %siteHash;

print "polyA sites are sorted into $n exons\n" if $verbose;


print "dump alt. polyA sites ...\n" if $verbose;

open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";
$i = 0;

foreach my $exonId (sort {$a cmp $b} keys %siteHash)
{

	print "$i ...\n" if $verbose && $i % 10000 == 0;
	$i++;

	#next unless $siteHash{$geneId}->{'N'} >= $minTotal;
	my @sitesOnExon = sort {$a->{'chromStart'} <=> $b->{'chromStart'}} @{$siteHash{$exonId}->{'sites'}};
	
	next unless @sitesOnExon > 1;

	my $strand = $sitesOnExon[0]->{'strand'};
	foreach my $s (@sitesOnExon)
	{
		if ($s->{'strand'} ne $strand)
		{
			print "inconsistency of strand for sites in $exonId\n";
			$strand = ".";
			last;
		}
	}

	next unless $strand eq '+' || $strand eq '-';
	
	my $e = $exonHashById{$exonId};
	@sitesOnExon = reverse (@sitesOnExon) if $strand eq '-';	
	
	my $iter = 0;
	foreach my $s (@sitesOnExon)
	{
		$s->{'name'} = $exonId . "[A$iter][". $s->{'score'} ."]";
		$iter++;
 
		my $s2 = $s;
		if ($strand eq '+')
		{
			$s2->{'chromStart'} = $e->{'chromStart'};
		}
		else
		{
			$s2->{'chromEnd'} = $e->{'chromEnd'};
		}
		print $fout bedToLine ($s2), "\n";
	}
}

close ($fout);

system ("rm -rf $cache");
