#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Carp;

use MyConfig;
use Bed;

my $prog = basename($0);

my $bigFile = 0;	#if yes, we need to use cache

my $sameStrand = 0;	#cluster regions are in the same strand for a match
my $keepScore = 0; #keep the original score of representative tag (for collapse mode)

#singleton means no overlap with other regions
#my $weight = 0;	#the way to calculate the score in each cluster, without weight means the number of tags, other wise, the sum of scores
my $sameBlockCount = 0; #require the same block count to combine
my $sameIntron = 0; #they have the same intron in overlapping region

my $debug = 0;
my $verbose = 0;

my $cache = getDefaultCache ($prog);

Carp::croak "$cache already exists\n" if -d $cache;
system ("mkdir $cache");

GetOptions (
			'big'=>\$bigFile,
			's|same-strand'=>\$sameStrand,
			'same-block-count'=>\$sameBlockCount,
			'same-intron'=>\$sameIntron,
			#'weight'=>\$weight,
			#'keep-score'=>\$keepScore,
			'c|cache:s'=>\$cache,
			'd|debug'=>\$debug,
			'v|verbose'=>\$verbose
			);


if (@ARGV != 2)
{
	print "remove a transript whose exons are covered by another transcript\n";
	
	print "Usage: $prog [options] <in.bed> <out.bed>\n";
	print " -big               : set when the input file is big\n";
	print " -s                 : same strand required [off]\n";
	print " --same-block-count : require the same block count to combine\n";
	print " --same-intron      : require the same introns in overlapping region\n";
	#print " -weight           : consider the weight of each tag\n";
	print " -c      [string]   : cache dir ($cache)\n";
	print " -d                 : debug (on|[off])\n";
	print " -v                 : verbose (on|[off])\n";
	exit (1);
}


my ($inBedFile, $outBedFile) = @ARGV;

my %tagCount;

print "reading transcripts from $inBedFile ...\n" if $verbose;
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

print "get transcripts count broken down into chromosomes ...\n" if $verbose;

foreach my $chrom (sort keys %tagCount)
{

	my $n = $tagCount{$chrom};
	$n = ref($n) eq 'HASH' ? $n = $n->{'n'} : @$n;

	print "$chrom : $n transcripts\n" if $verbose;
}


print "\n\nclustering transcripts ...\n" if $verbose;

my @strands = ('b');
@strands = qw(+ -) if $sameStrand;


my $fout;
open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";


foreach my $s (@strands)
{
	
	print "processing strand $s ...\n" if $verbose;
	foreach my $chrom (sort keys %tagCount)
	{
		my $tags;
		if ($bigFile)
		{
			my $tmpFile = $tagCount{$chrom}->{'f'};
			print "loading transcripts on chromsome $chrom from $tmpFile...\n" if $verbose;
			$tags = readBedFile ($tmpFile, $verbose);
		}
		else
		{
			$tags = $tagCount{$chrom};
		}

		my $n = @$tags;
		print "$n transcripts loaded on chromosome $chrom\n" if $verbose;
	
		Carp::croak "No strand specified\n" if $sameStrand && (not exists $tags->[0]->{"strand"});

		print "clustering $chrom ...\n" if $verbose;
		
		#my @tags4cluster;
		#my $i = 0;
		#foreach my $t (@$tags)
		#{
		#	next if $s ne 'b' && $t->{"strand"} ne $s;
		#	$t->{"idx"} = $i++;
		#	push @tags4cluster, $t;
		#}
		#
		#my $clusters = clusterRegions (\@tags4cluster, $s);

		#cluster if two transcripts overlap
		
		#definition: clusterRegions ($regionsOnChrom, $strand, $maxGap, $overlapFraction, $collapse); 
		my $clusters = clusterRegions ($tags, $s, 0, 0, 0);
		my $nc = @$clusters;
		print "\n\n$nc clusters found\n" if $debug;

		#strand has been handled in clustering, no worry any more

		my $iter = 0;
		foreach my $clust (@$clusters)
		{
			#sort according to the size of the region
			my @tagsInClust = sort {$b->{'chromEnd'} - $b->{'chromStart'} <=> $a->{'chromEnd'} - $a->{'chromStart'}}@$tags[@$clust];
			
			for (my $i = 0; $i < @tagsInClust; $i++)
			{
				my $tag1 = $tagsInClust[$i];
				next if exists $tag1->{'remove'};
				my $t1Key = regionToKey ($tag1);
				my $t1blockCount = exists $tag1->{'blockCount'} ? $tag1->{'blockCount'} : 1;
				for (my $j = $i + 1; $j < @tagsInClust; $j++)
				{
					my $tag2 = $tagsInClust[$j];
					my $t2blockCount = exists $tag2->{'blockCount'} ? $tag2->{'blockCount'} : 1;
					next if $sameBlockCount && $t1blockCount != $t2blockCount;
					next unless $tag1->{'chromStart'} <= $tag2->{'chromStart'} && $tag1->{'chromEnd'} >= $tag2->{'chromEnd'};
					
					if (not $sameIntron)
					{
						#combine if the exon of tag2 is completely overed by exons in tag1
						my $t3 = combineRegions ([$tag1, $tag2]);
						my $t3Key = regionToKey ($t3);
						$tag2->{'remove'} = 1 if $t1Key eq $t3Key;
					}
					else
					{
						#the new algorithm: combine if tag2 is a fragment of tag1
						my $t3 = segmentRegion ($tag1, $tag2->{'chromStart'}, $tag2->{'chromEnd'});
						next unless $t3;#they do not share exonic region

						my $t3Key = regionToKey ($t3);
						my $t2Key = regionToKey ($tag2);

						$tag2->{'remove'} = 1 if $t2Key eq $t3Key;
					}
				}
				print $fout bedToLine ($tag1), "\n";
			}
		} #cluster
	}#chrom
}#strand
close ($fout);

system ("rm -rf $cache");

sub regionToKey
{
	my ($r) = @_;

	my $key = $r->{"chrom"} . ":" . $r->{"chromStart"} . "-" . $r->{"chromEnd"};

	my %r2 = %$r;#make a copy
	bedToFull (\%r2) unless exists $r2{'blockCount'};
	$key .= ":" . join ("-", @{$r2{"blockStarts"}}, "//", @{$r2{"blockSizes"}});
	return $key;
}

