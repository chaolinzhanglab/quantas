#!/usr/bin/perl -w


use Getopt::Long;
use File::Basename;
use Carp;

use Bed;

my $prog = basename ($0);

my $verbose = 0;
GetOptions ("v"=>\$verbose);

if (@ARGV != 3)
{
	print "summarize polyA sites\n";
	print "Usage: $prog [options] <sample.count.bed> <polyA2gene.map> <out.txt>\n";
	print " -v : verbose\n";
	exit (1);
}


my ($polyACountBedFile, $polyAToGeneFile, $outFile) = @ARGV;


print "loading polyA count file ...\n" if $verbose;

my $fin;
my %polyAHash;
open ($fin, "<$polyACountBedFile") || Carp::croak "cannot open file $polyACountBedFile ...\n";
my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
    next if $line =~/^\s*$/;


	print "$iter ...\n" if $verbose && $iter % 50000 == 0;
	$iter++;

	my $b = lineToBed ($line);

	my $polyAId = $b->{"name"};
	$polyAHash{$polyAId} = $b;
}
close ($fin);

print "$iter polyA sites loaded\n" if $verbose;

print "loading polyA to gene mapping ...\n" if $verbose;

my %polyAToGeneHash;
open ($fin, "<$polyAToGeneFile") || Carp::croak "cannot open file $polyAToGeneFile ...\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	my ($polyAId, $geneId) = split ("\t", $line);
	Carp::croak "polyA site $polyAId does not exist in the count file\n" unless exists $polyAHash{$polyAId};	

	push @{$polyAToGeneHash{$geneId}}, $polyAHash{$polyAId};
}
close ($fin);


my $fout;

print "output results ...\n" if $verbose;
open ($fout, ">$outFile") || Carp::croak "cannot open $outFile to write\n";
$iter = 0;
print $fout join ("\t", "#chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "isoform", "isoform2Count", "isoform1Count", "geneId", "site1Pos", "site2Pos"), "\n";
foreach my $geneId (sort keys %polyAToGeneHash)
{
	print "$iter ...\n" if $verbose && $iter % 1000 == 0;
	$iter++;

	my $polyAs = $polyAToGeneHash{$geneId};
	next unless @$polyAs > 1;

	my $strand = $polyAs->[0]->{"strand"};
	for (my $i = 1; $i < @$polyAs; $i++)
	{
		Carp::croak "different strands found for polyA sites in gene $geneId\n" if $strand ne $polyAs->[$i]->{'strand'};
	}

	my @sortedPolyAs = $strand eq '+' ? (sort {$a->{"chromStart"} <=> $b->{'chromStart'}} @$polyAs) : (sort {$b->{"chromStart"} <=> $a->{'chromStart'}} @$polyAs);

	for (my $i = 0; $i < @sortedPolyAs; $i++)
	{
		my $p = $sortedPolyAs[$i];
		for (my $j = $i+1; $j < @sortedPolyAs; $j++)
		{
			my $q = $sortedPolyAs[$j];
			my ($chromStart, $chromEnd) = sort {$a <=> $b} ($p->{'chromStart'}, $q->{'chromStart'});
			$chromStart -= 50;
			$chromEnd += 50;

			my $name = join ("//", $p->{'name'}, $q->{'name'});

			print $fout join ("\t", $p->{'chrom'}, $chromStart, $chromEnd + 1, $name, $p->{'score'} + $q->{'score'}, $p->{'strand'},
				"apa", "P5/P3", $p->{'score'}, $q->{'score'}, $geneId, $p->{'chromStart'}, $q->{'chromStart'}), "\n";
				#$geneId, $p->{'chromStart'}, $q->{'chromStart'}, $p->{'score'}, $q->{'score'}), "\n";
		}
	}
}

close ($fout);






