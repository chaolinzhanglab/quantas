#!/usr/bin/env perl

#This program is to take the output of getExonIntron.pl to 
#get the relevant regions for alternative 5'/3' splice sites


use strict;
use warnings;

use File::Basename;
use Carp;
use Getopt::Long;
use Data::Dumper;

use Bed;
use Common;


my $prog = basename ($0);

my $inType = 'alt5';
my $minMinorEvi = 0;
my $minSpace = 0;
my $enumerate = 0; 
#print only the regions defined by the two most abundant isoforms
#otherwise print all that satisfy the filtering

my $verbose = 0;



GetOptions ('if:s'=>\$inType,
		'evi:i'=>\$minMinorEvi,
		'space:i'=>\$minSpace,
		'enum'=>\$enumerate,
		'v|verbose'=>\$verbose);

if (@ARGV != 3)
{
	print "get relevant regions for alt. 5'/3' ss\n";
	print "Usage: $prog [options] <altss.exon.bed> <out.constitutive.bed> <out.alt.bed>\n";
	print "OPTIONS:\n";
	print " -if     [string] : input type of altss.exon.bed ([alt5]|alt3)\n"; 
	print " -evi    [int]    : minimal number of supporting evidence for the minor isoform ($minMinorEvi)\n";
	print " -space  [int]    : minimal space between the most major and minor isoform ($minSpace)\n";
	print " -enum            : print all combinations\n";
	print " -v: verbose\n";
	exit (1);
}

my @types = qw(alt5 alt3);

if (checkParamType ($inType, \@types) eq '')
{
	Carp::croak "incorrect input type: $inType, must be one of ", join ("|", @types), "\n";
}


my ($ASBedFile, $constBedFile, $altBedFile) = @ARGV;

print "loading $inType from $ASBedFile ...\n" if $verbose;
my $ASExons = readBedFile ($ASBedFile, $verbose);

my %ASEHash;
foreach my $ase (@$ASExons)
{
	my $name = $ase->{"name"};
	my $blockCount = 1;
	$blockCount = $ase->{"blockCount"} if exists $ase->{"blockCount"};
	Carp::croak "$name: blockCount != 1\n" unless $blockCount == 1;
	$name =~/^(.*?)\[(.*?)\]/;

	my $caseId = $1;
	my $isoformId = $2;
	$ASEHash{$caseId}->{$isoformId} = $ase;
}

my $n = keys %ASEHash;

print "$n AS events loaded\n" if $verbose;


my (@constRegions, @altRegions);
my $iter = 0;

foreach my $caseId (sort keys %ASEHash)
{
	print "$iter: processing $caseId ...\n" if $verbose && $iter % 5000 == 0;
	$iter++;

	my $isoforms = $ASEHash{$caseId};
	my @isoformIds = keys %$isoforms;
	
	if (@isoformIds < 2)
	{
		warn "only one isoform for $caseId\n";
		next;
	}
	
	my $name = $isoforms->{$isoformIds[0]}->{"name"};
	
	$name =~/^.*?\[.*?\]\[(.*?)\]/;
	my $evi = $1;
	my @evis = split ("/", $evi);
	
	for (my $i = 0; $i < @evis; $i++)
	{
		my $isoformId = "A$i";
		$isoforms->{$isoformId}->{"idx"} = $i if exists $isoforms->{$isoformId};
		$isoforms->{$isoformId}->{"evi"} = $evis[$i] if exists $isoforms->{$isoformId};
	}

   	if ($enumerate)
	{
		@isoformIds = sort {$isoforms->{$a}->{"idx"} <=> $isoforms->{$b}->{'idx'}} keys %$isoforms;
	}
	else
	{
		@isoformIds = sort {$isoforms->{$b}->{'evi'} <=> $isoforms->{$a}->{'evi'}} keys %$isoforms; #sort according to abundances
		@isoformIds = @isoformIds[0..1];
	}

	for (my $i = 0; $i < @isoformIds; $i++)
	{
		my $isoformId1 = $isoformIds[$i];
		my $isoform1 = $isoforms->{$isoformId1};
		next if $isoform1->{"evi"} < $minMinorEvi;

		for (my $j= $i+1; $j < @isoformIds; $j++)
		{
			my $isoformId2 = $isoformIds[$j];
			my $isoform2 = $isoforms->{$isoformId2};
			next if $isoform1->{"evi"} < $minMinorEvi;
			
			#print "isoform1 : ", $isoformId1, " evi=", $isoform1->{"evi"}, "\n" if $verbose;
			#print "isoform2 : ", $isoformId2, " evi=", $isoform2->{"evi"}, "\n" if $verbose;
	
			my $space = 0;
			if (($inType eq 'alt5' && $isoform1->{"strand"} eq '+') || 
				($inType eq 'alt3' && $isoform1->{"strand"} eq '-'))
			{
				$space = ABS ($isoform1->{"chromEnd"} - $isoform2->{"chromEnd"});
			}
			else
			{
				$space = ABS ($isoform1->{"chromStart"} - $isoform2->{"chromStart"});
			}

			#print " space = $space\n" if $verbose;
			next if $space < $minSpace;
	
			my $name = $caseId . "[" . $isoformId1 . "/" . $isoformId2 . "][" . $isoform1->{"evi"} . "/" . $isoform2->{"evi"} . "]";
	
			my ($rc, $ra);
			if (($inType eq 'alt5' && $isoform1->{"strand"} eq '+') ||
				($inType eq 'alt3' && $isoform1->{"strand"} eq '-'))
			{
				$rc = {
					chrom => $isoform1->{"chrom"},
					chromStart => $isoform1->{"chromStart"},
					chromEnd => ($isoform1->{"chromEnd"} <= $isoform1->{"chromEnd"})? $isoform1->{"chromEnd"} : $isoform2->{"chromEnd"},
					name => $name,
					score => max($isoform1->{"score"}, $isoform2->{"score"}),
					strand => $isoform1->{"strand"}
				};
			
				$ra = {
					chrom => $isoform1->{"chrom"},
					chromStart => (($isoform1->{"chromEnd"} <= $isoform2->{"chromEnd"})? $isoform1->{"chromEnd"} : $isoform2->{"chromEnd"}) + 1,
					chromEnd => ($isoform1->{"chromEnd"} >= $isoform2->{"chromEnd"})? $isoform1->{"chromEnd"} : $isoform2->{"chromEnd"},
					name => $name,
					score => max ($isoform1->{"score"}, $isoform2->{"score"}),
					strand => $isoform1->{"strand"}
				};
			}
			else
			{
				$rc = {
					chrom => $isoform1->{"chrom"},
					chromEnd => $isoform1->{"chromEnd"},
					chromStart => ($isoform1->{"chromStart"} >= $isoform2->{"chromStart"})? $isoform1->{"chromStart"} : $isoform2->{"chromStart"},
					name => $name,
					score => max ($isoform1->{"score"},$isoform2->{"score"}),
					strand => $isoform1->{"strand"}
				};
		
				$ra = {
					chrom => $isoform1->{"chrom"},
					chromStart => ($isoform1->{"chromStart"} <= $isoform2->{"chromStart"})? $isoform1->{"chromStart"} : $isoform2->{"chromStart"},
					chromEnd => (($isoform1->{"chromStart"} >= $isoform2->{"chromStart"})? $isoform1->{"chromStart"} : $isoform2->{"chromStart"}) - 1,
					name => $name,
					score => $isoform1->{"score"},
					strand => $isoform1->{"strand"}
				};
			}
			push @constRegions, $rc;
			push @altRegions, $ra;
		}
	}
}

writeBedFile (\@constRegions, $constBedFile);
writeBedFile (\@altRegions, $altBedFile);


