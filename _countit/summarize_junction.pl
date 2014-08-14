#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

use MyConfig;
use Bed;
use Common;


my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $verbose = 0;
my $weight = 0;
my $big = 0;
my $separateStrand = 0;

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

GetOptions (
		'big'=>\$big,
		'weight'=>\$weight,
		'ss'=>\$separateStrand,
		'c|cache:s'=>\$cache,
		'keep-cache'=>\$keepCache,
		'v|verbose'=>\$verbose);

if (@ARGV != 3)
{
	print "summarize the number of reads for each junction\n";
	print "Usage: $prog [options] <intron.bed> <tag.bed> <intron.count.bed>\n";
	print " <intron.bed> -- bed file of unique introns\n";
	print " <tag.bed>    -- bed file of all tags, gz file allowed\n";
	print "OPTIONS:\n";
	#print " -big           : the tag file is big\n";
	print " -weight        : weight tags according to score\n";
	print " --ss           : consider the two strands separately\n";
	print " -c             : cache dir ($cache)\n";
	#print " --keep-cache   : keep cache files\n";
	print " -v             : verbose\n";
	exit (1);
}


my ($intronBedFile, $tagBedFile, $outBedFile) = @ARGV;
	
my $bigFlag = $big ? '-big' : '';
my $verboseFlag = $verbose ? '-v' : '';
my $ssFlag = $separateStrand ? '--ss' : '';


my %intronHash;

print "loading introns from $intronBedFile ...\n" if $verbose;
my $introns = readBedFile ($intronBedFile, $verbose);

my $iter = 0;
foreach my $i (@$introns)
{
	my $chrom = $i->{'chrom'};
	my $chromStart = $i->{'chromStart'};
	my $chromEnd = $i->{'chromEnd'};
	my $strand = $i->{'strand'};

	my $intronKey = "$chrom:$chromStart:$chromEnd";
	$intronKey .= ":$strand" if $separateStrand;
	$i->{'iter'} = $iter++;
	$i->{'score'} = 0;
	$intronHash{$intronKey} = $i;	
}

print "$iter introns loaded\n" if $verbose;

print "counting junction reads ...\n" if $verbose;

my $fin;
open ($fin, "<$tagBedFile") || Carp::croak "cannot open file $tagBedFile to read\n"; 

if ($tagBedFile =~/\.gz$/)
{
	open ($fin, "gunzip -c $tagBedFile | ")||Carp::croak "cannot open file $tagBedFile to read\n";
}
else
{
	open ($fin, "<$tagBedFile") || Carp::croak "cannot open file $tagBedFile to read\n";
}



$iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	print "$iter ...\n" if $iter % 100000 == 0 && $verbose;
	$iter++;	

	my $t = lineToBed ($line);
	next unless exists $t->{'blockCount'} && $t->{'blockCount'} > 1;

	my $strand = $t->{'strand'};
	my $chrom = $t->{'chrom'};
 
	for (my $i = 0; $i < $t->{'blockCount'} -1; $i++)
	{
		my $intronStart = $t->{'chromStart'} + $t->{'blockStarts'}->[$i] + $t->{'blockSizes'}->[$i];
		my $intronEnd = $t->{'chromStart'} + $t->{'blockStarts'}->[$i+1] - 1;

		my $key = "$chrom:$intronStart:$intronEnd";
        $key .= ":$strand" if $separateStrand;
        next unless exists $intronHash{$key};
		
		my $score = $weight ? $t->{'score'} : 1;
		$intronHash{$key}{'score'} += $score;
	}
}

close ($fin);

print "dumping output ...\n" if $verbose;

my $fout;
open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";
foreach my $i (sort {$a->{'iter'} <=> $b->{'iter'}} values %intronHash)
{
	print $fout bedToLine ($i), "\n";
}
close ($fout);



