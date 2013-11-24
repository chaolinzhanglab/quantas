#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;
use Bio::SeqIO;

use Motif;
use Common;


my $prog = basename($0);
my $verbose = 0;

my $singleStrand = 0;
my $threshold = 0;
my $functionalDepth = 0;
my $baseCompStr = "";
my $motifType = "count";

GetOptions ('s|single-strand'=>\$singleStrand,
		'motif-type:s'=>\$motifType,
		'c:s'=>\$baseCompStr,
		't|threshold:f'=>\$threshold,
		'func-depth'=>\$functionalDepth,
		'v'=>\$verbose
);


if (@ARGV != 3)
{
	print "search motif sites\n";
	print "Usage: $prog [options] <input.fa> <motif.txt> <out.bed>\n";
	print "OPTION:\n";
	print " --motif-type [string]: ([count]|pwm)\n";
	print " -c           [string]: background base composition (for count matrix)\n";
	print " -t           [float] : motif threshold ($threshold)\n";
	print " --func-depth         : use functional depth rather than the PWM score\n";
	print " -s                   : search only the given strand\n";
	print " -v                   : verbose\n";
	exit (1);
}

my ($fastaFile, $motifFile, $outBedFile) = @ARGV;# = $ARGV[0];

my ($a, $c, $g, $t) = (0,0,0,0);
if ($baseCompStr ne '')
{
	($a, $c, $g, $t) = split (",", $baseCompStr);
}
elsif ($motifType eq 'count')
{
	print "calculate base composition ...\n" if $verbose;
	
	my $seqIO = Bio::SeqIO->new (-file =>$fastaFile, -format => 'Fasta');
	while (my $seq = $seqIO->next_seq())
	{
		my $seqStr = $seq->seq();

		my $seqId = $seq->id();

		print "$seqId ...\n" if $verbose;
		my $a2 = ($seqStr=~tr/aA//);
		my $c2 = ($seqStr=~tr/cC//);
		my $g2 = ($seqStr=~tr/gG//);
		my $t2 = ($seqStr=~tr/tT//);

		$a += $a2;
		$g += $g2;
		$c += $c2;
		$t += $t2;
	}
}

my $N=$a + $g + $c + $t;

my $baseComp;

if ($N > 0)
{
	$baseComp = $singleStrand? {A=>$a / $N, C=>$c/ $N, G=>$g/$N, T=>$t/$N} : {A=>($a+$t)/2/$N, C=> ($g+$c)/2/$N, G=> ($g+$c)/2/$N, T=>($a+$t)/2/$N};

	print "base comp: A=", $baseComp->{'A'}, ", C=", $baseComp->{'C'}, ", G=", $baseComp->{'G'}, ", T=", $baseComp->{'T'}, "\n" if $verbose;
}

print "read motif file $motifFile ...\n" if $verbose;

my $motif = readMotifFile ($motifFile);
my $matrix = getMatrix($motif->[0]);

$matrix = countToStormoMatrix ($matrix, $baseComp) if $motifType eq 'count';

print "search motif occurrence ...\n";

my $seqIO = Bio::SeqIO->new (-file =>$fastaFile, -format => 'Fasta');

my $motifWidth = @$matrix;

my $maxScore = getMaxMatrixScore ($matrix);
my $minScore = getMinMatrixScore ($matrix);


print "maxScore = $maxScore, minScore=$minScore\n" if $verbose;
my $matrixRC = revComMatrix($matrix);

my $fout;

open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile\n";

my $iter = 0;
my $j = 0;
while (my $seq = $seqIO->next_seq())
{
	my $seqStr = $seq->seq();
	my $seqLen = $seq->length ();
	
	my $chr = $seq->id ();
	$chr =~s/\.fa$//g;

	print "$j ...\n" if $j % 10000 == 0 && $verbose;
	$j++;

	for (my $i = 0; $i < $seqLen - $motifWidth + 1; $i++)
	{
		#print "$i ...\n" if $i % 100000 == 0 && $verbose;

		my $siteSeq = substr ($seqStr, $i, $motifWidth);
		
		my $score = getMatrixScore ($matrix, $siteSeq);
		$score = ($score - $minScore) / ($maxScore - $minScore) if $functionalDepth;

		if ($score >= $threshold)
		{
			print $fout join ("\t", $chr, $i, $i + $motifWidth, "s_" . $iter, $score, '+'), "\n";
			$iter++;
		}
		
		next if $singleStrand;

		my $scoreRC = getMatrixScore ($matrixRC, $siteSeq);
		$scoreRC = ($scoreRC  - $minScore) / ($maxScore - $minScore);

		if ($scoreRC >= $threshold)
		{
			print $fout join ("\t", $chr, $i, $i + $motifWidth, "s_" . $iter, $scoreRC, '-'), "\n";
			$iter++;
		}
	}
}

close ($fout);







