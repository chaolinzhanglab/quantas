#!/usr/bin/perl -w

use strict;


use Getopt::Long;
use Carp;
use Data::Dumper;
use File::Basename;

use Vcf;

my $prog = basename ($0);

my $vcfListFile = "";
my $minAltBase = 0;
my $minDepth = 0;
my $fixAltBase = 0;

my $verbose = 0;

GetOptions (
	"l:s"=>\$vcfListFile,
	"d:i"=>\$minDepth,
	"a:i"=>\$minAltBase,
	"x"=>\$fixAltBase,
	"v"=>\$verbose);

if (@ARGV < 1 && $vcfListFile eq '')
{
	print "extract SNVs from RNASeq bam files\n";
	print "Usage: $prog [options] <in1.vcf> [in2.vcf]\n";
	print "Note: in DP4, counts of REF is not accurate due to missing info in some samples\n";
	print " [options\n";
	print " -l [int]   : list of vcf files\n";
	print " -d [int]   : minimum depth ($minDepth)\n";
	print " -a [int]   : minimum alt base to report a potential SNV ($minAltBase)\n";
	print " -x         : do not reassign altBase\n";
	print " -v         : verbose\n";
	exit (1);
}

my @vcfFiles = @ARGV;

if (-f $vcfListFile)
{
	my $fin;
	open ($fin, "<$vcfListFile") || Carp::croak "cannot open $vcfListFile ...\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\#/;
		next if $line =~/^\s*$/;
		push @vcfFiles, $line;
	}
	close ($fin);
}

my $n = @vcfFiles;
print STDERR "$n vcf files to be processed\n" if $verbose;


my %snvHash;


my $iter = 0;
foreach my $vcfFile (@vcfFiles)
{
	print STDERR "$iter of $n: $vcfFile ...\n" if $verbose;
	$iter++;

	my $fin;
	open ($fin, "<$vcfFile") || Carp::croak "cannot open file $vcfFile to read\n";


	my $i = 0;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		
		print STDERR "$i ...\n" if $verbose && $i % 10000 == 0;
		$i++;

		my $snv = lineToVcf ($line);
		my $chrom = $snv->{'chrom'};
		my $position = $snv->{'position'};
		my $dp4 = $snv->{'info'}->{'DP4'};
		my ($refPos, $refNeg, $altPos, $altNeg) = split (/\,/, $dp4);
		my $refBase = $snv->{'refBase'};
		my $altBase = $snv->{'altBase'};
	
		if (exists $snvHash{$chrom} && $snvHash{$chrom}->{$position})
		{
			my $existingSnv = $snvHash{$chrom}->{$position};
			$existingSnv->{'info'}->{'BC'} += ($altPos + $altNeg > 0 ? 1 : 0);
			$existingSnv->{'+'}->{$refBase} += $refPos;
			$existingSnv->{'-'}->{$refBase} += $refNeg;
			$existingSnv->{'+'}->{$altBase} += $altPos;
			$existingSnv->{'-'}->{$altBase} += $altNeg;
			$existingSnv->{'id'} = $snv->{'id'};
		}
		else
		{
			map{$snv->{'+'}->{$_} = 0} qw(A C G T);
			map{$snv->{'-'}->{$_} = 0} qw(A C G T);
			$snv->{'+'}->{$refBase} += $refPos;
			$snv->{'-'}->{$refBase} += $refNeg;
			$snv->{'+'}->{$altBase} += $altPos;
			$snv->{'-'}->{$altBase} += $altNeg;
			$snv->{'info'}->{'BC'} = $altPos + $altNeg > 0 ? 1 : 0;
			$snvHash{$chrom}->{$position} = $snv;
		}
	}
}

srand (1); #used by assignAltBase

print join ("\t", "#CHROM", "POS ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), "\n";

foreach my $chrom (sort keys %snvHash)
{			
	my $snvChrom = $snvHash{$chrom};
	my $n = keys %$snvChrom;

	print STDERR "$chrom ($n snvs)...\n" if $verbose;
	
	my $iter = 0;
	foreach my $pos (sort {$a <=>$b} keys %$snvChrom)
	{
		print STDERR "$iter ...\n" if $verbose && $iter % 10000 == 0;
		$iter++;

		my $snv = $snvChrom->{$pos};

		my %readBaseHash = (
			'A'=>$snv->{'+'}->{'A'} + $snv->{'-'}->{'A'},
			'C'=>$snv->{'+'}->{'C'} + $snv->{'-'}->{'C'},
			'G'=>$snv->{'+'}->{'G'} + $snv->{'-'}->{'G'},
			'T'=>$snv->{'+'}->{'T'} + $snv->{'-'}->{'T'});

		my $refBase = $snv->{'refBase'};
    	my $altBase = $fixAltBase ? $snv->{'altBase'} : assignAltBase (\%readBaseHash, $refBase);
		
		my $altBaseSum = $readBaseHash{$altBase};
		my $total = $readBaseHash{$refBase} + $altBaseSum;
		next unless $altBaseSum >= $minAltBase && $total >= $minDepth;
		
		my $DP4 = join(",", $snv->{'+'}->{$refBase}, $snv->{'-'}->{$refBase}, $snv->{'+'}->{$altBase}, $snv->{'-'}->{$altBase});
		
		#note pos is zero-based coordinates
		$snv->{'altBase'} = $altBase;
		$snv->{'info'}->{'DP4'} = $DP4;
		print vcfToLine ($snv), "\n";
	}
}


