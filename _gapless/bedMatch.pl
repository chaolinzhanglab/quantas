#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Getopt::Long;

use MyConfig;
use Bed;
use Common;

my $prog = basename($0);

my $cache = MyConfig::getDefaultCache($prog);

my $verbose = 0;
my $sameStrand = 0;
my $separateStrand = 0;

my $useType = 0;
my $diffThreshold = 25;
my $delimiter = "//";

my $useBlock = 0;
my $sameName = 0;

my $big = 0;

my $keepScore = 1;


GetOptions ("s"=>\$sameStrand, #obsolete, keep it here for backward compatibility
		"ss"=>\$separateStrand,
		'big'=>\$big,
		'same-name'=>\$sameName,
		't|use-type'=>\$useType,
		'b|use-block'=>\$useBlock,
		'd:i'=>\$diffThreshold,
		"dd:s"=>\$delimiter,
		'keep-score:i'=>\$keepScore,
		'cache:s'=>\$cache,
		"v|verbose"=>\$verbose);

if (@ARGV != 3)
{
	print "match exact regions in two bed files by coordinates and strand\n";
	
	print "Usage: $prog [options] <first.bed> <second.bed> <out.bed>\n";
	print " -ss              : same strand required [off]\n";
	print " -big             : big file\n";
	print " -same-name       : require the same name for a match\n";
	print " --use-block      : consider exon/intron structure\n";
	print " --use-type       : treat terminal or internal exons differently\n";
	print "                    (5,0,3 in the score column for 5' terminal, internal, 3' terminal, respectively)\n";
	print " -d          [int]: tolerance of length difference for terminal exons ($diffThreshold)\n";
	print " -keep-score [int]: keep score of the first or second bed file ([1]|2)\n";
	print " -dd      [string]: delimiter ($delimiter)\n";
	print " -cache   [string]: cache dir($cache)\n"; 
	print " -v               : verbose [off]\n";
	exit (1);
}


my ($firstBedFile, $secondBedFile, $outBedFile) = @ARGV;

$sameStrand = $sameStrand | $separateStrand;

if ($useType)
{
	Carp::croak "must specify strand to treat terminat/internal exons differently\n" if $sameStrand == 0;
}

if ($useBlock && $useType)
{
	Carp::croak "will not consider terminal exons if block information is used\n";
}


my $cacheSub1 = "$cache/bed1";
my $cacheSub2 = "$cache/bed2";

my $firstRegionsHash ={};
my $secondRegionsHash ={};

if ($big)
{
	system ("mkdir $cache");
	system ("mkdir $cacheSub1");
	system ("mkdir $cacheSub2");

	print "split $firstBedFile by chromosomes ...\n" if $verbose;
	$firstRegionsHash = splitBedFileByChrom ($firstBedFile, $cacheSub1, $verbose);
	
	print "split $secondBedFile by chromosomes ...\n" if $verbose;
	$secondRegionsHash = splitBedFileByChrom ($secondBedFile, $cacheSub2, $verbose);


}
else
{
	$firstRegionsHash->{'all'} = {f=>$firstBedFile};
	$secondRegionsHash->{'all'} = {f=>$secondBedFile};
}


unlink $outBedFile if -f $outBedFile; #make sure the file is cleared because we will use append later

foreach my $chrom (sort keys %$firstRegionsHash)
{
	next unless exists $secondRegionsHash->{$chrom};

	print "processing chromosome $chrom ...\n" if $verbose && $big;

	my $f1 = $firstRegionsHash->{$chrom}->{'f'};
	my $f2 = $secondRegionsHash->{$chrom}->{'f'};

	print "loading regions from $f1 ...\n" if $verbose;
	my $firstRegions = readBedFile ($f1, $verbose);

	my %regionHash;

	foreach my $r (@$firstRegions)
	{
		if ($sameStrand)
		{
			Carp::croak "must specify strand" unless exists $r->{"strand"};
		}
		my $key = regionToKey ($r, $sameStrand); 
		push @{$regionHash{$key}}, $r;
	}

	my $nregion1 = @$firstRegions;
	print "$nregion1 regions loaded\n" if $verbose;

	print "reading regions from $f2 ...\n" if $verbose;
	my $secondRegions = readBedFile ($f2, $verbose);
	my $nregion2 = @$secondRegions;

	print "$nregion2 regions loaded\n" if $verbose;

	my @outRegions;
	foreach my $r2 (@$secondRegions)
	{
		if ($sameStrand)
		{
			Carp::croak "must specify strand", Dumper ($r2),  unless exists $r2->{"strand"};
		}
		my $key = regionToKey ($r2, $sameStrand); 

		if (exists $regionHash{$key})
		{
			foreach my $r1 (@{$regionHash{$key}})
			{
				my $diff = 0;
				if ($useType)
				{
					#keep only if the length are very similar for terminal exons
					next if $r1->{"score"} != $r2->{"score"}; # must be the same type
					$diff = ABS ($r1->{"chromStart"} - $r2->{"chromStart"}) + ABS ($r1->{"chromEnd"} - $r2->{"chromEnd"});

					#$diff = $r1->{"score"} == 5 ? Common::ABS ($r1->{"chromStart"} - $r2->{"chromStart"}) : Common::ABS ($r1->{"chromEnd"} - $r2->{"chromEnd"}) 
					#if $r1->{"strand"} eq '+';

					#$diff = $r1->{"score"} == 5 ? Common::ABS ($r1->{"chromEnd"} - $r2->{"chromEnd"}) : Common::ABS ($r1->{"chromStart"} - $r2->{"chromStart"})
					#if $r1->{"strand"} eq '-';

					#print "r1=", $r1->{"name"}, ", r2=", $r2->{"name"}, ", diff=$diff\n";
					next unless $diff <= $diffThreshold;
				}

				my %r1New = %$r1;
				$r1New{"name"} = join ($delimiter, $r1->{"name"}, $r2->{"name"});
				$r1New{"name"} .= $delimiter.$diff;
				$r1New{"score"} = $r2->{"score"} if $keepScore == 2 && exists $r2->{"score"};

				push @outRegions, \%r1New;
			}
		}
	}

	my $nout = @outRegions;
	print "$nout matches found\n" if $verbose;
	writeBedFile (\@outRegions, $outBedFile, '', 'a');
}

system ("rm -rf $cache") if -d $cache;

sub regionToKey
{
	my ($r, $sameStrand) = @_;

	my $key = $r->{"chrom"} . ":" . $r->{"chromStart"} . "-" . $r->{"chromEnd"};
	
	if ($useType)
	{
		Carp::croak "must specify strand to treat terminat/internal exons differently\n" unless exists $r->{"strand"};
		if ($r->{"score"} == 5)
		{
			$key = $r->{"chrom"} . ":" . $r->{"chromEnd"} if $r->{"strand"} eq '+';
			$key = $r->{"chrom"} . ":" . $r->{"chromStart"} if $r->{"strand"} eq '-';
		}
		elsif ($r->{"score"} == 3)
		{
			$key = $r->{"chrom"} . ":" . $r->{"chromStart"} if $r->{"strand"} eq '+';
			$key = $r->{"chrom"} . ":" . $r->{"chromEnd"} if $r->{"strand"} eq '-';
		}
		elsif ($r->{"score"} != 0)
		{
			Carp::croak "score must be 5, 0 or 3 if --use-type is on\n";
		}
	}
	elsif ($useBlock)
	{
		Carp::croak "must specify strand to treat terminat/internal exons differently\n" unless exists $r->{"strand"};
		my %r2 = %$r;#make a copy
		bedToFull (\%r2) unless exists $r2{'blockCount'};

		$key .= ":" . join ("-", @{$r2{"blockStarts"}}, "//", @{$r2{"blockSizes"}});
		%r2=();
	}

	$key .= ":" . $r->{"name"} if $sameName;
	$key .= ":" . $r->{"strand"} if $sameStrand && exists $r->{"strand"};
	
	return $key;
}

