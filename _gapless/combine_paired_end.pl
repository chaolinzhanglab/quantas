#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use File::Basename;

use MyConfig;
use Bed;
use Common;

my $prog = basename ($0);
my $geneContigBedFile = "";
#my $separateStrand = 0;
my $libraryType = "unstranded";
# or 'stranded-sa': read1 is sense and read2 is antisense
# 'stranded-as' : read1 is antisense and read2 is sense


my $insertSizeDistFile = "";
my %insertSizeDist = ();
my $maxSize = 400;
my $usePrior = 0;

#my $ignoreJunctionFromInput = 0;

my $printSingleton = 0;
my $big = 0;

my $cache = getDefaultCache ($prog);

my $verbose = 0;
my $debug = 0;

GetOptions ("v|verbose"=>\$verbose,
	"big"=>\$big,
	"gene:s"=>\$geneContigBedFile,
	'size-dist:s'=>\$insertSizeDistFile,
	'use-prior'=>\$usePrior,
	#"ss|sep-strand"=>\$separateStrand,
	"library-type:s"=>\$libraryType,
	#"ignore-junction"=>\$ignoreJunctionFromInput,
	"print-singleton"=>\$printSingleton,
	"cache:s"=>\$cache);


#assume 5' and 3' reads have the same name
#different pairs have unique names
#


#notes for stranded vs unstranded library, Chaolin Zhang, June 16, 2014
#for unstranded library, we assume read1 and read2 are on the opposite strand
#in the output, successfully paired reads will adopt the strand of read1, singleton reads will adopt the original strand

#for stranded library, we also gave the sense strand, based on the library type. specifically
#for stranded-as: the strand of paired read or singleton read1 is flipped
#for stranded-sa: the strand of singleton read2 is flipped



if (@ARGV != 3)
{
	print "combine paired reads in the same gene and with proper location/orientation\n";
	print "Usage: $prog [options] <5end.bed> <3end.bed> <pair.bed>\n";
	print " -gene           [file]  : gene contig bed file to limit paired reads in the same gene\n";
	print " -use-prior              : use informative priors specified in score column of the gene file\n";
	print " -size-dist      [file]  : file specifying insert size distribution\n";
	#print " -ss              : separate strand\n";
	print " --library-type  [string]: library  type ([unstranded]|stranded-as|stranded-sa)\n";
	#print " --ignore-junction: ignore SE junction reads in the input\n";
	print " -big                    : big read file\n";
	print " -print-singleton        : print reads not in pairs\n";
	print " -cache                  : cache dir ($cache)\n";
	print " -v                      : verbose\n";
	exit (1);
}


print "use prior=$usePrior\n" if $verbose;
my ($fivePrimeBedFile, $threePrimeBedFile, $combineBedFile) = @ARGV;

if ($libraryType ne 'unstranded' && $libraryType ne 'stranded-sa' && $libraryType ne 'stranded-as')
{
	Carp::croak "incorrect library type: $libraryType\n";
}


Carp::croak "$cache already exist\n" if -d $cache;

my ($contigCache, $fivePrimeCache, $threePrimeCache) = ("$cache/contig", "$cache/fivePrime", "$cache/threePrime");

if ($big)
{
	system ("mkdir $cache");
	system ("mkdir $contigCache");
	system ("mkdir $fivePrimeCache");
	system ("mkdir $threePrimeCache");
}

if ($insertSizeDistFile ne '')
{
	print "reading size distribution of cDNA insert from $insertSizeDistFile ...\n" if $verbose;

	Carp::croak "cannot open file $insertSizeDistFile ...\n" unless -f $insertSizeDistFile;
	my $fin;
	open ($fin, "<$insertSizeDistFile") || Carp::croak "cannot open file $insertSizeDistFile to read\n";
	my $freqTotal = 0;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		my ($len, $freq) = split ("\t", $line);
		next if $len >= $maxSize;

		$insertSizeDist{$len} = $freq;
		$freqTotal += $freq;
	}
	close ($fin);

	foreach my $len ((1..$maxSize))
	{
		if (exists $insertSizeDist{$len})
		{
			$insertSizeDist{$len} /= $freqTotal;
		}
		else
		{
			$insertSizeDist{$len} = 1e-6;#pseudo count
		}
	}
}


################################################################################################
#
my $contigHash = {};
Carp::croak "$geneContigBedFile does not exist\n" unless -f $geneContigBedFile;

my $n = 0;

print "loading gene contigs from $geneContigBedFile...\n" if $verbose;
if ($big)
{
	$contigHash = splitBedFileByChrom ($geneContigBedFile, $contigCache, $verbose);
	#Carp::croak Dumper ($contigHash), "\n";

	foreach my $chrom (sort keys %$contigHash)
	{
		$n+= $contigHash->{$chrom}->{'n'};
	}

	Carp::croak "empty file: $geneContigBedFile\n" unless $n > 0;
	print "$n genes loaded\n" if $verbose;
}
else
{
	my $contigs = readBedFile ($geneContigBedFile, $verbose);
	my $n = @$contigs;
	Carp::croak "empty file: $geneContigBedFile\n" unless $n > 0;
	print "$n genes loaded\n" if $verbose;

	my $i = 0;
	foreach my $c (@$contigs)
	{
		my $chrom = $c->{"chrom"};
		#$c->{"score"} = 0;
		#$c->{"name"} = $c->{"chrom"} . ":" . ($c->{"chromStart"} + 1) . "-" . ($c->{"chromEnd"} + 1) unless exists $c->{"name"};
		
		#note we now change the name of isoforms
		#so that they are unqiue and short
		$c = AnnotationIO::bed2Full ($c) unless exists $c->{'blockCount'};
		$c->{'name'} = $i;
		$i++;

		#if ($separateStrand)
		if ($libraryType ne 'unstranded')
		{
			Carp::croak "no strand information in ", Dumper ($c), "\n" unless exists $c->{"strand"};
			my $strand = $c->{"strand"};
			push @{$contigHash->{$chrom}->{$strand}}, $c;
		}
		else
		{
			push @{$contigHash->{$chrom}->{'b'}}, $c;
		}
	}
}


##########################################################################################################
my $fivePrimeChromHash = {};
my $threePrimeChromHash = {};

print "reading 5' reads from $fivePrimeBedFile ...\n" if $verbose;
if ($big)
{
	$fivePrimeChromHash = splitBedFileByChrom ($fivePrimeBedFile, $fivePrimeCache, $verbose);
}
else
{
	my $fivePrimeReads = readBedFile ($fivePrimeBedFile, $verbose);

	foreach my $r (@$fivePrimeReads)
	{
		#organize according to chrom, to search overlapping genes
		my $chrom = $r->{"chrom"};
		push @{$fivePrimeChromHash->{$chrom}}, $r;
	}
}

print "get read count broken down into chromosomes ...\n" if $verbose;

foreach my $chrom (sort keys %$fivePrimeChromHash)
{
	my $n = ref ($fivePrimeChromHash->{$chrom}) eq 'ARRAY' ? @{$fivePrimeChromHash->{$chrom}} : $fivePrimeChromHash->{$chrom}->{"n"};
	print "$chrom : $n reads\n" if $verbose;
}

print "reading 3' reads from $threePrimeBedFile ...\n" if $verbose;
if ($big)
{
	$threePrimeChromHash = splitBedFileByChrom ($threePrimeBedFile, $threePrimeCache, $verbose);
}
else
{
	my $threePrimeReads = readBedFile ($threePrimeBedFile, $verbose);

	foreach my $r (@$threePrimeReads)
	{
		#organize according to chrom, to search overlapping genes
		my $chrom = $r->{"chrom"};
		push @{$threePrimeChromHash->{$chrom}}, $r;
	}
}

foreach my $chrom (sort keys %$threePrimeChromHash)
{
	my $n = ref ($threePrimeChromHash->{$chrom}) eq 'ARRAY' ? @{$threePrimeChromHash->{$chrom}} : $threePrimeChromHash->{$chrom}->{"n"};
	print "$chrom : $n reads\n" if $verbose;
}


#############################################################################

#get the number of chroms with reads
my %chromHash;

foreach my $chrom ((keys %$fivePrimeChromHash, keys %$threePrimeChromHash))
{
	$chromHash{$chrom} = 1;
}

my $fout;
open ($fout, ">$combineBedFile") || Carp::croak "can not open file $combineBedFile to write\n";

foreach my $chrom (sort keys %chromHash)
{
	print "processing chrom $chrom ...\n" if $verbose;

	my ($fivePrimeReads, $threePrimeReads) = ([], []);
	my $contigsOnChrom = {};

	#load the reads on the chrom
	if ($big)
	{
		if (exists $contigHash->{$chrom} && exists $contigHash->{$chrom}->{'f'})
		{
			print "loading contigs from ", $contigHash->{$chrom}->{'f'}, "\n" if $verbose;
			my $contigs = readBedFile ($contigHash->{$chrom}->{'f'}, $verbose);
			my $n = @$contigs;
			print "$n contigs loaded\n" if $verbose;

			my $i = 0;
			foreach my $c (@$contigs)
			{
				#note that we changed the name of the contig here
				$c = AnnotationIO::bed2Full ($c) unless exists $c->{'blockCount'};
				$c->{'name'} = $i;
				$i++;

				#my $s = $separateStrand ? $c->{'strand'} : 'b';
				my $s = $libraryType eq 'unstranded' ? 'b' : $c->{'strand'};
				push @{$contigsOnChrom->{$s}}, $c;
			}
		}

		if (exists $fivePrimeChromHash->{$chrom} && exists $fivePrimeChromHash->{$chrom}->{"f"})
		{
			print "loading 5' reads from ", $fivePrimeChromHash->{$chrom}->{"f"}, "\n" if $verbose;
			$fivePrimeReads = readBedFile ($fivePrimeChromHash->{$chrom}->{"f"}, $verbose);
		
			my $n = @$fivePrimeReads;
			print "$n 5' reads loaded\n" if $verbose;
		}

		if (exists $threePrimeChromHash->{$chrom} && exists $threePrimeChromHash->{$chrom}->{"f"})
		{
			print "loading 3' reads from ", $threePrimeChromHash->{$chrom}->{"f"}, "\n" if $verbose;
			$threePrimeReads = readBedFile ($threePrimeChromHash->{$chrom}->{"f"}, $verbose);

			my $n = @$threePrimeReads;
			print "$n 3' reads loaded\n" if $verbose;
		}
	}
	else
	{
		#contig names have been fixed at this point
		$contigsOnChrom = $contigHash->{$chrom} if exists $contigHash->{$chrom};
		$fivePrimeReads = $fivePrimeChromHash->{$chrom} if exists $fivePrimeChromHash->{$chrom};
		$threePrimeReads = $threePrimeChromHash->{$chrom} if exists $threePrimeChromHash->{$chrom};
	}
	

	#find overlapping genes for each of the reads
	my @strand = ('b');
	#@strand = qw (+ -) if $separateStrand;
	@strand = qw (+ -) if $libraryType ne 'unstranded';

	
	print "find overlapping genes for each reads on $chrom ...\n" if $verbose;

	foreach my $s (@strand)
	{
		print "processing strand $s ...\n" if $verbose;
		next unless exists $contigsOnChrom->{$s};
		
		my @tags = ();

		foreach my $t (@$fivePrimeReads)
		{
			$t->{'3end'} = $t->{"strand"} eq '+' ? $t->{'chromEnd'} : $t->{"chromStart"};
			$t->{'len'} = exists $t->{'blockCount'} ? sum($t->{'blockSizes'}) : ($t->{'chromEnd'} - $t->{'chromStart'} + 1);

			#push @tags, $t if $s eq 'b' || $t->{"strand"} eq $s;
			#changed by Chaolin Zhang June 16, 2014
			if ($libraryType eq 'unstranded')
			{
				push @tags, $t;
			}
			elsif ($libraryType eq 'stranded-sa')
			{
				#require read1 to be on the same strand of contig
				push @tags, $t if $t->{'strand'} eq $s;
			}
			elsif ($libraryType eq 'stranded-as')
			{
				#require read1 to be on the opposite strand of config
				push @tags, $t if $t->{'strand'} ne $s;
			}
			else
			{
				Carp::croak "incorrect library type: $libraryType\n";
			}
		}
		my $n = @tags;

		#we cluster 5' and 3' tags separately to save some memory
		print "clustering $n 5' tags to genes on strand $s of $chrom ...\n" if $verbose;
		findOverlapTags (\@tags, $contigsOnChrom->{$s}) if @tags > 0;
		#note that strand requirement is already fullfilled in the filtered tags, so we do not check strand in this function

		@tags = ();

		foreach my $t (@$threePrimeReads)
		{
			$t->{'3end'} = $t->{"strand"} eq '+' ? $t->{'chromEnd'} : $t->{"chromStart"};
			$t->{'len'} = exists $t->{'blockCount'} ? sum($t->{'blockSizes'}) : ($t->{'chromEnd'} - $t->{'chromStart'} + 1);
			
			#push @tags, $t if $s eq 'b' || $t->{"strand"} ne $s; #note that the strand is flipped here
			#changed by Chaolin Zhang June 16, 2014

            if ($libraryType eq 'unstranded')
            {
                push @tags, $t;
            }
            elsif ($libraryType eq 'stranded-sa')
            {
                #require read2 to be on the opposite strand of contig
                push @tags, $t if $t->{'strand'} ne $s;
            }
            elsif ($libraryType eq 'stranded-as')
            {
                #require read2 to be on the same strand of config
                push @tags, $t if $t->{'strand'} eq $s;
            }
            else
            {
                Carp::croak "incorrect library type: $libraryType\n";
            }
		}

		$n = @tags;
		print "clustering $n 3' tags to genes on strand $s of $chrom ...\n" if $verbose;
		findOverlapTags (\@tags, $contigsOnChrom->{$s}) if @tags > 0;
		#note that strand requirement is already fullfilled in the filtered tags, so we do not check strand in this function

		@tags = ();
	}


	print "hashing 5' reads on $chrom by name ...\n" if $verbose;

	my %fivePrimeHash;
	my $iter = 0;

	foreach my $r5 (@$fivePrimeReads)
	{
		$iter++;
		print "$iter ...\n" if $verbose && $iter % 10000 == 0;

		$fivePrimeHash{$r5->{"name"}} = $r5;
	}


	print "matching 5' and 3' read mates on $chrom ...\n" if $verbose;

	$iter = 0;
	#my $found = 0; #debug

	foreach my $r3 (@$threePrimeReads)
	{
		$iter++;
		print "$iter ...\n" if $verbose && $iter % 10000 == 0;

		my $name = $r3->{"name"};
		#$found = 1 if $debug && $name eq 'HWUSI-EAS1600:DKO2_250:11_30_09:6:103:165:493#0';

		if (not exists $fivePrimeHash{$name})
		{
			$r3->{"name"} .= "//3";
			$r3->{"score"} = 1;
			
			flipStrand ($r3) if $libraryType eq 'stranded-sa'; 
			#this is because we know read2 is on antisense strand

			print $fout bedToLine (bedToFull($r3)), "\n" if $printSingleton;
			next;
		}

		my $r5 = $fivePrimeHash{$name};
		
		#check consistency
		my $chromStart1 = $r5->{"chromStart"};
		my $chromEnd1 = $r5->{"chromEnd"};
		my $chromStart2 = $r3->{"chromStart"};
		my $chromEnd2 = $r3->{"chromEnd"};

		my $sameGene = 0;

		my @paths;

		foreach my $g (keys %{$r3->{"gene"}})
		{
			if (exists $r5->{"gene"}->{$g})
			{
				$sameGene++;
				push @paths, $r5->{"gene"}->{$g};
			}
		}

		if ($sameGene && $r5->{"chrom"} eq $r3->{"chrom"}
			&& $r5->{"strand"} ne $r3->{"strand"}
			&& (($r5->{"strand"} eq '+' && $chromStart1 <= $chromStart2 && $chromEnd1 <= $chromEnd2) 
					|| $r5->{"strand"} eq '-' && $chromStart2 <= $chromStart1 && $chromEnd2 <= $chromEnd1))
		{
			#the two reads are on the same chromosome and have proper orientation
			#so check if they are a legal pair more carefully

			$r5->{"match"} = 1;
			
=debug
			if ($debug && $found)
			{
				my $nPath = @paths;

				print "track name=reads\n";
				print AnnotationIO::printBedRegionToString ($r5), "\n";
				print AnnotationIO::printBedRegionToString ($r3), "\n";

				print "track name=path description=\"$nPath paths found between ", Common::min ($r5->{"3end"}, $r3->{"3end"}), " and ", Common::max ($r5->{"3end"}, $r3->{"3end"}), "\"\n";
				#print Dumper (\@paths), "\n";

				foreach my $p (@paths)
				{
					print AnnotationIO::printBedRegionToString ($p), "\n";
				}
			}
=cut
			if ((-f $insertSizeDistFile) && 
					(($r5->{'strand'} eq '+' && $r5->{'3end'} < $r3->{'3end'}) || ($r5->{'strand'} eq '-' && $r5->{'3end'} > $r3->{'3end'})))
			{
				#there is a gap, we do probablistic inference of all possibilities
				#and the size distribution was provided

				#get all unique paths between the two reads (consider only the 3'end of each read)
				#my $uniqPaths = getUniqPaths (\@paths, Common::min ($r5->{"3end"}, $r3->{"3end"}), Common::max ($r5->{"3end"}, $r3->{"3end"}), $found);
				my $uniqPaths = getUniqPaths2 (\@paths, min ($r5->{"3end"}, $r3->{"3end"}), max ($r5->{"3end"}, $r3->{"3end"}));


=debug			
				if ($debug && $found)
				{
					my $nPath = @$uniqPaths;
					print "track name=uniq_path description=$nPath\n";
					foreach my $p (@$uniqPaths)
					{
						print AnnotationIO::printBedRegionToString ($p), "\n";
					}
				}
=cut
				my $probs = inferPosterior ($uniqPaths, $r5->{"len"} + $r3->{"len"}, \%insertSizeDist, $usePrior);

				for (my $i = 0; $i < @$uniqPaths; $i++)
				{
					my $pair = combineRegions ([$r5, $r3, $uniqPaths->[$i]]);
					$pair->{"score"} = $probs->[$i];
					$pair->{"name"} .= "//p$i//" . sprintf("%.2f", $probs->[$i]);
					
					#$pair adopts the same strand as $r5
					flipStrand ($pair) if $libraryType eq 'stranded-as';

					print $fout bedToLine ($pair), "\n";
				}
			}
			else
			{
				#the two reads overlap, so no gap in between. we can combine the two reads directly
				#or the pair length distribution was not provided, so we can not do the inference

				my $pair = combineRegions ([$r5, $r3]);
				$pair->{'score'} = 1;
				$pair->{'name'} .="//p";

				flipStrand ($pair) if $libraryType eq 'stranded-as';
				print $fout bedToLine ($pair), "\n";

			}
			#my $rp = combineRegions ([$r5, $r3]);
			
			#$rp->{"name"} .= "//p";
			#print $fout AnnotationIO::printBedRegionToString ($rp), "\n";
		}
		else
		{
			#illegal pair, although both reads mapped successfully
			$r3->{"name"} .= "//3";
			$r3->{"score"} = 1;
			flipStrand ($r3) if $libraryType eq 'stranded-sa';

			print $fout bedToLine (bedToFull($r3)), "\n" if $printSingleton;
		}
		#Carp::croak "done\n" if $found && $debug;
	}

	if ($printSingleton)
	{
		foreach my $r5 (@$fivePrimeReads)
		{
			$r5->{"name"} .= "//5";
			$r5->{"score"} = 1;

			flipStrand ($r5) if $libraryType eq 'stranded-as';
			print $fout bedToLine (bedToFull($r5)), "\n"
			unless exists $r5->{"match"};
		}
	}
}

close ($fout);

system ("rm -rf $cache");


sub inferPosterior
{
	my ($paths, $sequencedEndsSize, $sizeDist, $usePrior) = @_;
	return [1] if (@$paths < 2); #only one possibility

	#Carp::croak Dumper ($paths), "\n";
	foreach my $p (@$paths)
	{
		$p->{'len'} = exists $p->{'blockSizes'} ? sum ($p->{'blockSizes'}) : ($p->{"chromEnd"} - $p->{"chromStart"} + 1);
		#this is the length between the 3' most nucleotide of a pair
		$p->{'len'} += $sequencedEndsSize - 2;
		#this is the length between the 5' most nucleotide of a pair (the complete cDNA insert)

		$p->{'len'} = $maxSize if $p->{'len'} > $maxSize;
	}

	my @dist;
   
	if ($usePrior)	
	{
		#the score of each path is the number of junction reads  per junction, 
		#reflecting the prior knowledge of the relative abundance
		#Note: if a path do not have junctions (e.g. intron retention), there might be some problems
		#but in this case, the sizes usually differ dramatically, so that the prior play a minor role

		@dist = map {$sizeDist->{$_->{'len'}}} @$paths;
		#print "noninformative: ", join ("\t", @dist), "\n";

		@dist = map {$sizeDist->{$_->{'len'}} * ($_->{'score'} + 1)} @$paths; #1 is the pseudo count
		#Carp::croak "informative: ", join ("\t", @dist), "\n";
	}
	else
	{
		#non-informative prior
		@dist = map {$sizeDist->{$_->{'len'}}} @$paths;
	}

	#print Dumper ($paths), "\n";
	#Carp::croak Dumper (\@dist), '\n';

	#assume non-informative priors
	my $denom = sum (\@dist);
	my @posteriors = map {$_/$denom} @dist; 

	return \@posteriors;
}

sub getUniqPaths2
{
	#my ($paths, $chromStart, $chromEnd, $printDetail) = @_;# for debug
	my ($paths, $chromStart, $chromEnd) = @_;
	#paths have full 12 column format

	my %uniqPaths;
	foreach my $p (@$paths)
	{
		my $blockCount = $p->{"blockCount"};

		my $segment = segmentRegion ($p, $chromStart, $chromEnd);
		next unless ref($segment) && $segment->{"chromStart"} == $chromStart && $segment->{"chromEnd"} == $chromEnd;
		#start and end must be on exons

		if (not exists ($segment->{"blockCount"}))
		{
			$segment = bedToFull($segment);
		}

		next unless $segment->{'blockCount'} == $blockCount;
		#Nov. 17, 2010
		#
		#TO DO: if the path has additional junctions beyond the range, we do not consider because the prior score is not accurate
		#in this case, we assume there is a better path that just span the junctions defined by the range
		#
		#this is true in most of cases in the current method to build the isoform database, with an exception for exons
		#exon pairs or trios tolerate terminal exons observed in EST/mRNAs
		#but only terminal exons observed from UCSC/RefSeq known genes were included
		#
		#this can reject some of the retained introns, but these are rare and exons were not scored anyway
		#

		#with the same chromStart and end, block sizes and starts are sufficient to distinguish diffrerent clusters


		my $key = join ("-", @{$segment->{"blockSizes"}},"//", @{$segment->{"blockStarts"}});

		#if redundant paths exists, we first choose the one with least junctions and then choose the shortest one
		#to minimize irrelevant region outside the interval
		#this may affect the prior

=debug
if ($debug && $printDetail)
{
		if (exists $uniqPaths{$key})
		{
			print "compare current uniq path and the new path\n";
			print "current=", AnnotationIO::printBedRegionToString ($uniqPaths{$key}->{'p'}), "\n";
			print "current block count=", $uniqPaths{$key}->{'p'}->{'blockCount'}, "\n";
			print "current length=", ( $uniqPaths{$key}->{'p'}->{"chromEnd"} - $uniqPaths{$key}->{'p'}->{"chromStart"}), "\n";
			print "new path=", AnnotationIO::printBedRegionToString ($p), "\n";
			print "new block count=", $blockCount, "\n";
			print "new length=", $p->{"chromEnd"} - $p->{"chromStart"}, "\n";

		}
		else
		{
			print "find new unique path\n";
			print AnnotationIO::printBedRegionToString ($p), "\n";
		}

		#Carp::croak Dumper ($uniqPaths{$key}), "\n";
		#print Dumper ($p), "\n";
}
=cut
		if ((not exists $uniqPaths{$key})
			|| $uniqPaths{$key}->{'p'}->{'blockCount'} > $blockCount 
			|| ($uniqPaths{$key}->{'p'}->{'blockCount'} == $blockCount && $uniqPaths{$key}->{'p'}->{"chromEnd"} - $uniqPaths{$key}->{'p'}->{"chromStart"} > $p->{"chromEnd"} - $p->{"chromStart"}))
		{
			#$uniqPaths{$key} = $segment;
			$uniqPaths{$key} = {s=>$segment,p=>$p};
			#print $p->{'name'}, " is used as new unique path\n" if $debug && $printDetail;
		}
	
	}
	#my @uniqPaths = values %uniqPaths;
	my @uniqPaths = map {$uniqPaths{$_}->{'s'}} keys %uniqPaths;

	return \@uniqPaths;	
}



#find overlapping genes for each tag
sub findOverlapTags
{
	my ($tagsUnSorted, $regionsUnSorted) = @_;

	#here we only need to consider the last position of the reads
	my @tags = sort {$a->{"3end"} <=> $b->{"3end"}} @$tagsUnSorted;
	#my $tags = \@tags;

	my @regions = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regionsUnSorted;
	#my $regions = \@regions;
	
	my $firstTagIdx = 0; #the first tag that overlap with or on the right of the current window

	foreach my $r (@regions)
	{
		
		my $chromStart = $r->{"chromStart"};
		my $chromEnd = $r->{"chromEnd"};
		
		while ($firstTagIdx < @tags && $tags[$firstTagIdx]->{"3end"} < $chromStart)
		{
			$firstTagIdx++;
		}
		
		my $i = $firstTagIdx;
		
		while ($i < @tags && $tags[$i]->{"3end"} <= $chromEnd)
		{
			if ($tags[$i]->{"3end"} >= $chromStart)
			{
				$tags[$i]->{"gene"}->{$r->{"name"}} = $r;
			}
			$i++;
		}
	}
	@tags = ();
	@regions = ();
}

sub flipStrand
{
	my $read = $_[0];
	Carp::croak "incorrect strand info:", Dumper ($read), "\n"
	unless exists $read->{'strand'} && ($read->{'strand'} eq '+' || $read->{'strand'} eq '-');

	if ($read->{'strand'} eq '+')
	{
		$read->{'strand'} = '-';
	}
	elsif ($read->{'strand'} eq '-')
	{
		$read->{'strand'} = '+';
	}
}


