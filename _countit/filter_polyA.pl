#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;
use Carp;


use Bed;

my $prog = basename ($0);
my $verbose = 0;

my $minTotal = 20;
my $minRatio = 0.2;



GetOptions ("v"=>\$verbose,
	'N:i'=>\$minTotal,
	'r:f'=>\$minRatio);

if (@ARGV != 3)
{
	print "filter polyA sites from polyA seq\n";
	print " Usage: $prog [options] <in.bed> <site2gene> <out.bed>\n";
	print " -N [int]  : min total number of tags for the gene ($minTotal)\n";
	print " -r [float]: min ratio of tags for the site among all sites in the gene ($minRatio)\n"; 
	print " -v        : verbose\n";
	exit (1);
}


my ($inBedFile, $siteToGeneFile, $outBedFile) = @ARGV;


print "read polyA sites from $inBedFile ...\n" if $verbose;

my $sites = readBedFile ($inBedFile, $verbose);

my $n = @$sites;

print "$n sites loaded\n" if $verbose;

print "load site-2-gene mapping ...\n" if $verbose;

my $fin;
open ($fin, "<$siteToGeneFile") || Carp::croak "cannot open file $siteToGeneFile to read\n";

my %siteToGeneHash;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	my ($siteId, $geneId) = split ("\t", $line);
	
	$siteToGeneHash{$siteId} = $geneId;
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
	my $geneId = $siteToGeneHash{$siteId};
	push @{$siteHash{$geneId}->{'sites'}}, $s;
	$siteHash{$geneId}->{'N'} += $s->{'score'};
}

$n = keys %siteHash;

print "polyA sites are sorted into $n genes\n" if $verbose;

my $fout;

print "filter polyA sites ...\n" if $verbose;

open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";
$i = 0;

foreach my $geneId (sort {$a <=> $b} keys %siteHash)
{

	print "$i ...\n" if $verbose && $i % 10000 == 0;
	$i++;

	next unless $siteHash{$geneId}->{'N'} >= $minTotal;
	my $sites = $siteHash{$geneId}->{'sites'};
	
	foreach my $s (@$sites)
	{
		next unless $s->{'score'} / $siteHash{$geneId}->{'N'} >= $minRatio;
		#$s->{'name'} .= "//" . $geneId;

		print $fout bedToLine ($s), "\n";
	}
}

close ($fout);




