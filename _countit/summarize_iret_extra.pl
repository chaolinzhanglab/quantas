#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Carp;

use Bed;
use Common;

my $prog = basename ($0);
my $readLen = 1;
my $verbose = 0;

GetOptions ("read-len=i"=>\$readLen,
		"v"=>\$verbose);

if (@ARGV != 6)
{
	print "summarize extra metrics for iret\n";
	print "Usage: $prog [options] <intron.id.map.txt> <iret.summary.txt> <expr.summary.txt> <intron.count.txt> <exon.count.txt> <out.txt>\n";
	print "read-len [int]: read length\n";
	print " -v : verbose\n";
	exit (1);
}


my ($intronIdMapFile, $iretSummaryFile, $exprSummaryFile, $intronCountFile, $exonCountFile, $outFile) = @ARGV;


print "loading intron info ...\n" if $verbose;
my $fin;
open ($fin, "<$intronIdMapFile") || Carp::croak "cannot open file $intronIdMapFile to read\n";
my %intronHash;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\#/;
	my ($intronId, $geneId, $geneSymbol, $upExonId, $dnExonId, $overlapRefSeqInternalExon) = split ("\t", $line);
	$intronHash{$intronId} = {gid=>$geneId, symbol=>$geneSymbol, up=>$upExonId, dn=>$dnExonId, refseq=>$overlapRefSeqInternalExon};
}

close ($fin);

print "loading iret summary ...\n" if $verbose;
open ($fin, "<$iretSummaryFile") || Carp::croak "cannot open file $iretSummaryFile to read\n";
while (my $line = <$fin>)
{
	chomp $line;
    next if $line=~/^\#/;
	my @cols = split ("\t", $line);
	my $intronId = $cols[3];
	$intronHash{$intronId}{'iret'} = \@cols;
}
close ($fin);


print "loading expression ...\n" if $verbose;

open ($fin, "<$exprSummaryFile") || Carp::croak "cannot open file $exprSummaryFile to read\n";
my %exprHash;

my $totalTagNum = 0;
while (my $line = <$fin>)
{
    chomp $line;
    next if $line=~/^\#/;
    my ($gid, $symbol, $tagNum, $exonLen, $RPKM) = split ("\t", $line);
	$exprHash{$gid} = $RPKM; #{len=>$exonLen, n=>$tagNum, RPKM=>$RPKM};
	$totalTagNum += $tagNum;
}
close ($fin);

=remove
#recalculate rpkm
foreach my $gid (keys %exprHash)
{
	my $g = $exprHash{$gid};
	my $tagNum = $g->{'n'};
	$tagNum = 1 if $tagNum == 0;
	$g->{'RPKM'} = $tagNum * 1e9 / $exonLen / $totalTagNum;
}
=cut

print "loading exon count ...\n" if $verbose;

my %exonHash;
open ($fin, "<$exonCountFile") || Carp::croak "cannot open file $exonCountFile to read\n";
while (my $line = <$fin>)
{
    chomp $line;
    next if $line=~/^\#/;
	my $e = lineToBed ($line);
	my $ec = $e->{'score'};
	my $el = $e->{'chromEnd'} - $e->{'chromStart'} + 1;
	my $eRPKM = $ec * 1e9 / ($el+$readLen - 1) / $totalTagNum;
	$exonHash{$e->{'name'}} = $eRPKM;
}
close ($fin);


print "summarize iret metrics ...\n";
my $fout;
open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
open ($fin, "<$intronCountFile") || Carp::croak "cannot open file $intronCountFile to read\n";


print $fout join ("\t", "#chrom", "chromStart", "chromEnd", "name", "estScore", "strand",
           "gene", "upstream_exon_id", "downstream_exon_id", "overlap_refseq_internal_exon", "coverage", "imbalance", "PIR1f", "PIR1c", "PIR2f", "PIR2c",
            "intron_read_count", "retain_5ss_count", "retain_3ss_count", "spliced_count", "intron_RPKM", "upstream_exon_RPKM", "downstream_exon_RPKM", "gene_RPKM"), "\n";

while (my $line = <$fin>)
{
	chomp $line;
    next if $line=~/^\#/;
	my $i = lineToBed ($line);

	my $name = $i->{'name'};
	my $ic = $i->{'score'};
	my $il = $i->{'chromEnd'} - $i->{'chromStart'} + 1;
	
	my $iRPKM = $ic * 1e9 / ($il+$readLen - 1) / $totalTagNum;

	my $g = $intronHash{$name};
	my $gid = $g->{'gid'};
	my $symbol = $g->{'symbol'};
	my $ueid = $g->{'up'};
	my $deid = $g->{'dn'};
	my ($chrom, $chromStart, $chromEnd, $name2, $estScore, $strand, $type, $isoformIDs, $isoform1Tags, $isoform2Tags, $retain5SSTag, $retain3SSTag, $juncTag) = @{$g->{'iret'}};

	my $gRPKM = exists $exprHash{$gid} ? $exprHash{$gid} : 'NA';
	my $ueRPKM = exists $exonHash{$ueid} ? $exonHash{$ueid} : 'NA';
	my $deRPKM = exists $exonHash{$deid} ? $exonHash{$deid} : 'NA';

	my $coverage = min($retain5SSTag, $retain3SSTag) + $juncTag;
	my $imbalance = $retain5SSTag + $retain3SSTag > 0 ? abs ($retain5SSTag - $retain3SSTag) / ($retain5SSTag + $retain3SSTag) : 0;
	
	my $PIR1f = $isoform1Tags + 2 * $isoform2Tags > 0 ? $isoform1Tags/($isoform1Tags + 2 * $isoform2Tags) * 100 : 'NA';
	my $PIR1c = $coverage > 0 ? min($retain5SSTag, $retain3SSTag) / $coverage * 100 : 'NA';
	
	my $PIR2f = $ueRPKM + $deRPKM > 0 ? $iRPKM * 2 / ($ueRPKM + $deRPKM) * 100 : 'NA';
	my $PIR2c = max ($ueRPKM, $deRPKM) > 0 ? $iRPKM / max ($ueRPKM, $deRPKM) * 100 : 'NA';
	$PIR2c = max ($ueRPKM, $deRPKM, $gRPKM) > 0 ? $iRPKM / max ($ueRPKM, $deRPKM, $gRPKM) * 100 : 'NA' if $gRPKM ne 'NA';
	if ($PIR2f ne 'NA')
	{
		$PIR2f = 100 if $PIR2c > 100;
		$PIR2f = 100 if $PIR2c > 100;	
	}

	print $fout join ("\t", $chrom, $chromStart, $chromEnd, $name, $estScore, $strand, 
			"$gid//$symbol", $ueid, $deid, $g->{'refseq'}, $coverage, $imbalance, $PIR1f, $PIR1c, $PIR2f, $PIR2c, 
			$ic, $retain5SSTag, $retain3SSTag, $juncTag, $iRPKM, $ueRPKM, $deRPKM, $gRPKM), "\n";
}

close ($fin);








