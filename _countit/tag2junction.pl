#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Carp;

use Bio::SeqIO;

use MyConfig;
use Bed;
use Common;

my $prog = basename($0);

my $bigFile = 0;	#if yes, we need to use cache
my $separateStrand = 0;
my $weight = 0;	#the way to calculate the score in each cluster, without weight means the number of tags, other wise, the sum of scores
my $genomeDir = "";
my $organism = "mm9";
my $verbose = 0;
my $canonical = 0; #GT/AG, GC/AG, AT/AC
my $anchor = 0;	#min size of the terminal block

my $cache = getDefaultCache ($prog);

GetOptions (
			'big'=>\$bigFile,
			'weight'=>\$weight,
			'ss'=>\$separateStrand,
			'can'=>\$canonical,
			'a:i'=>\$anchor,
			'org:s'=>\$organism,
			'gd:s'=>\$genomeDir,
			'c|cache:s'=>\$cache,
			'v|verbose'=>\$verbose
			);

if (@ARGV != 2)
{
	print "identify unique junctions from junctionreads and infer strand\n";
	
	print "Usage: $prog [options] <in.bed> <out.bed>\n";
	print " -big             : set when the input file is big\n";
	print " -ss              : separate strand\n";
	print " -weight          : consider the weight of each tag\n";
	print " -can             : canonical splice sites only (GT/AG, GC/AG, AT/AC)\n";
	print " -a      [int]    : anchor size of the terminal block ($anchor)\n";
	print " -org    [string] : organism <[mm9]|mm8|mm6|mm5|rn4|rn3|hg17|hg18|sacCer1|dm2|ce2>\n";
	print " -gd     [string] : directory of fasta files [$genomeDir] (will override -org)\n";
	print " -c      [string] : cache dir ($cache)\n";
	print " -v               : verbose (on|[off])\n";
	exit (1);
}


if ($genomeDir ne '')
{
	Carp::croak "$genomeDir does not exists\n" unless -d $genomeDir;
	$genomeDir = getFullPath ($genomeDir);
}
elsif ($organism ne '')
{
	#if organism is specified, use that
	$genomeDir = getGenomeDir($organism);
	$genomeDir .="/fasta";
	Carp::croak "$genomeDir does not exists\n" unless -d $genomeDir;
}
else
{
	Carp::croak "No species specified\n";
}


my ($inBedFile, $outBedFile) = @ARGV;

Carp::croak "$cache already exists\n" if -d $cache;
system ("mkdir $cache");


my %tagCount;

print "reading tags from $inBedFile ...\n" if $verbose;
if ($bigFile)
{
	my $ret = splitBedFileByChrom ($inBedFile, $cache, $verbose);
	%tagCount = %$ret;
}
else
{
	my $tags = readBedFile ($inBedFile, $verbose);
	foreach my $t (@$tags)
	{
		my $chrom = $t->{"chrom"};
		push @{$tagCount{$chrom}}, $t;
	}
}

print "get tag count broken down into chromosomes ...\n" if $verbose;

foreach my $chrom (sort keys %tagCount)
{

	my $n = $tagCount{$chrom};
	$n = ref($n) eq 'HASH' ? $n = $n->{'n'} : @$n;
	print "$chrom : $n tags\n" if $verbose;
}


print "\n\nclustering tags ...\n" if $verbose;


my @strands = qw(+ -);


my $fout;
open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";


foreach my $chrom (sort keys %tagCount)
{
	#load tags
	my $tags;
	if ($bigFile)
	{
		my $tmpFile = $tagCount{$chrom}->{'f'};
		print "loading tags on chromsome $chrom from $tmpFile...\n" if $verbose;
		$tags = readBedFile ($tmpFile, $verbose);
	}
	else
	{
		$tags = $tagCount{$chrom};
	}

	my $n = @$tags;
	print "$n tags loaded on chromosome $chrom\n" if $verbose;
	Carp::croak "No block information\n" unless exists $tags->[0]->{'blockStarts'};

	
	#identify introns
	print "clustering $chrom ...\n" if $verbose;
	my %introns;
	foreach my $t (@$tags)
	{
		my $strand = $t->{'strand'};
		my $chrom = $t->{'chrom'};
		for (my $i = 0; $i < $t->{'blockCount'} -1; $i++)
		{
			my $upstreamExonStart = $t->{'chromStart'} + $t->{'blockStarts'}->[$i];
			my $intronStart = $t->{'chromStart'} + $t->{'blockStarts'}->[$i] + $t->{'blockSizes'}->[$i];
			my $intronEnd = $t->{'chromStart'} + $t->{'blockStarts'}->[$i+1] - 1;
			my $downstreamExonEnd = $intronEnd + $t->{'blockSizes'}->[$i+1];
			
			#if the terminal block is too short, skip the junction
			if ($i == 0)
			{
				next if $t->{'blockSizes'}->[$i] < $anchor;
			}
			
			if ($i == $t->{'blockCount'} - 2)
			{
				next if $t->{'blockSizes'}->[$i+1] < $anchor;
			}

			my $key = "$chrom:$intronStart:$intronEnd";
			$key .= ":$strand" if $separateStrand;
			my $score = $weight ? $t->{'score'} : 1;
			if (exists $introns{$key})
			{
				my $i = $introns{$key};
				$i->{'chromStart'} = $upstreamExonStart if $upstreamExonStart < $i->{'chromStart'};
				$i->{'chromEnd'} = $downstreamExonEnd if $downstreamExonEnd > $i->{'chromEnd'};
				$i->{'count'} += $score;
			}
			else
			{
				$introns{$key} = {chromStart=>$upstreamExonStart, chromEnd=>$downstreamExonEnd, count=>$score};
			}
		}
	}

	$n = keys %introns;
	print "$n introns identified on chrom $chrom\n" if $verbose;


	#output
	my $chromSeqStr = "";
	
	if ($canonical)
	{
		print "loading sequences of $chrom ...\n" if $verbose;
		my $chromFastaFile = "$genomeDir/$chrom.fa";
		my $seqIO = Bio::SeqIO->new (-file=>$chromFastaFile, format=>'Fasta');
		$chromSeqStr = $seqIO->next_seq()->seq();
	}

	my $iter = 0;
	foreach my $key (sort keys %introns)
	{
		my $i = $introns{$key};

		my ($chrom, $intronStart, $intronEnd, $strand) = split (":", $key);
		my $chromStart = $i->{'chromStart'};
		my $chromEnd = $i->{'chromEnd'};

		if ($canonical)
		{
			my $donor = substr ($chromSeqStr, $intronStart, 2);
			my $acceptor = substr ($chromSeqStr, $intronEnd - 1, 2);
			
			my $spliceSites = uc ($donor. $acceptor);
			Carp::croak "incorret sequences for splice site: $spliceSites\n" unless length($spliceSites) == 4;
			
			if ($separateStrand)
			{
				#get the sense strand of the splice sites
				$spliceSites = revcom ($spliceSites) if $strand eq '-';
				my $isCanonical = $spliceSites eq 'GTAG' || $spliceSites eq 'GCAG' || $spliceSites eq 'ATAC';
				if ($isCanonical == 0)
				{
					warn "non-canonical splice sites=$spliceSites, skipped: $key\n";
					next;
				}
			}
			else
			{
				my $isCanonical = $spliceSites eq 'GTAG' || $spliceSites eq 'GCAG' || $spliceSites eq 'ATAC';
				$isCanonical = $isCanonical || $spliceSites eq 'CTAC' || $spliceSites eq 'CTGC' || $spliceSites eq 'GTAT';
				if ($isCanonical == 0)
				{
					warn "non-canonical splice sites=$spliceSites, skipped: $key\n";
					next;
				}
				
				#try to infer the sense strand based on the splice sites
				$strand = ($spliceSites eq 'GTAG' || $spliceSites eq 'GCAG' || $spliceSites eq 'ATAC') ? '+' : '-';
			}
		}
		
		my $name = "$chrom-j-$iter";
		$iter++;

		$strand = '+' unless defined $strand; #in this case, the reads are not strand specific, and cannot be inferred

		print  $fout join ("\t", $chrom, $chromStart, $chromEnd + 1, sprintf("%s[%d]", $name, int($i->{'count'}+0.5)), $i->{'count'}, $strand,
			$chromStart, $chromEnd + 1, 0, 2, 
			join(",", $intronStart - $chromStart, $chromEnd - $intronEnd),
			join(",", 0, $intronEnd + 1 - $chromStart)), "\n";
	}
	print "$iter introns passed filtering on chrom $chrom\n" if $verbose;
}
close ($fout);

system ("rm -rf $cache");




