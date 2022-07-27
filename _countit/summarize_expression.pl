#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;
use Data::Dumper;

use Bed;


my $prog = basename ($0);
my $exon2geneFile = "";
my $verbose = "";

GetOptions ("e2g:s"=>\$exon2geneFile,
		"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print "estimate gene expression level\n";
	print "Usage: $prog [options] <exon.count.bed> <out.txt>\n";
	print " -e2g   [file]: exon2gene map file\n";
	print " -v           : verbose\n";
	exit (0);
}
my ($exonBedFile, $outFile) = @ARGV;

print "loading exons from $exonBedFile ...\n" if $verbose;

my $exons = readBedFile ($exonBedFile, $verbose);

my $nexon = @$exons;

print "$nexon exons loaded ...\n";

my %exon2geneHash;
if (-f $exon2geneFile)
{
	print "loading exon2gene mapping from $exon2geneFile ...\n";
	my $fin;
	open ($fin, "<$exon2geneFile") || Carp::croak "can not open file $exon2geneFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my ($exonId, $geneId, $geneSymbol) = split (/\t/, $line);
		next unless $geneId && length($geneId) > 0;
		$exon2geneHash{$exonId} = {geneid=> $geneId, symbol=>$geneSymbol};
	}
	close ($fin);
}


#also get the mapping from exon ids if available


print "counting tag numbers in each gene ...\n";

my %genes;
my $totalTagNumber = 0;
foreach my $e (@$exons)
{
	my $exonId = $e->{"name"};
	#get additional gene ids from exon id
	#this is not really necessary
	#if ($exonId =~/^\w{2}\-(\d+)\-/)
	#{
	#	$exon2geneHash{$exonId} = {geneid=>$1};
	#}

	next unless exists $exon2geneHash{$exonId};
	my $geneId = $exon2geneHash{$exonId}->{'geneid'};

	$genes{$geneId}->{"length"} += $e->{"chromEnd"} - $e->{"chromStart"} + 1;
	$genes{$geneId}->{"count"} += $e->{"score"};
	$genes{$geneId}->{'symbol'} = $exon2geneHash{$exonId}->{'symbol'};
	$totalTagNumber += $e->{'score'};
}

my $fout;
open ($fout, ">$outFile") || Carp::croak "can not open file $outFile to write\n";
print $fout "#", join ("\t", "gene_id", "gene_symbol", "tag_num", "exon_len", "RPKM"), "\n";

foreach my $geneId (sort {$a <=> $b} keys %genes)
{
	my $g = $genes{$geneId};
	#Carp::croak Dumper ($g), "\n";
	print $fout join ("\t", $geneId, $g->{'symbol'}, $g->{"count"}, $g->{"length"}, $g->{"count"}/$g->{"length"} * 1e9/$totalTagNumber), "\n";
}
close ($fout);






