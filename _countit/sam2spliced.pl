#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Carp;
use File::Basename;
use Getopt::Long;

use Sam;
use Bed;

my $prog = basename ($0);
my $printUniqOnly = 0;

my $donorBedFile = "";
my $acceptorBedFile = "";

my $filter = "splice"; #unsplice

my $verbose = 0;

GetOptions (
	"u|uniq"=>\$printUniqOnly,
	"5:s"=>\$donorBedFile,
	"3:s"=>\$acceptorBedFile,
	"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print STDERR "filter sam files to get certain types of reads\n";
	print STDERR "Usage: $prog [options] <in.sam> <out.sam>\n";
	print STDERR " <in.sam> : gzip compressed input file with .gz extension is allowed\n";
	print STDERR " <out.sam>: output sam file\n";
	print STDERR " You can also use - to specify STDIN for input or STDOUT for output\n\n";
	print STDERR "options:\n";
	print STDERR " -5 [string]: a bed file of donor splice sites\n";
	print STDERR " -3 [string]: a bed file of acceptor splice sites\n";
	print STDERR " -u         : limit to uniquely mapped reads\n";
	print STDERR " -f [string]: filter\n";
	print STDERR " -v,--verbose:         verbose\n";
	exit (1);
}

my ($inSamFile, $outSamFile) = @ARGV;

Carp::croak "no splice sites specified\n" if $donorBedFile eq '' || $acceptorBedFile eq '';

my %spliceSiteHash;
if (-f $donorBedFile)
{
	print STDERR "load splice site file $donorBedFile ...\n" if $verbose;
	my $spliceSites = readBedFile ($donorBedFile, $verbose);
	foreach my $s (@$spliceSites)
	{
		my $ssKey = genKey ($s->{'chrom'}, $s->{'chromStart'}, $s->{'strand'});
		$spliceSiteHash{'donor'}{$ssKey} = $s;
	}
}

if (-f $acceptorBedFile)
{
	print STDERR "load splice site file $acceptorBedFile ...\n" if $verbose;
	my $spliceSites = readBedFile ($acceptorBedFile, $verbose);
	foreach my $s (@$spliceSites)
	{
		my $ssKey = genKey ($s->{'chrom'}, $s->{'chromStart'}, $s->{'strand'});
		$spliceSiteHash{'acceptor'}{$ssKey} = $s;
	}
}


my ($fin, $fout);

if ($inSamFile eq "-")
{
    $fin = *STDIN;
}
else
{
	if ($inSamFile =~/\.gz$/)
	{
		open ($fin, "gunzip -c $inSamFile | ") || Carp::croak "cannot open file $inSamFile to read\n";
	}
	elsif ($inSamFile =~/\.bz2$/)
    {
        open ($fin, "bunzip2 -c $inSamFile | ") || Carp::croak "cannot open file $inSamFile to read\n";
    }
	elsif ($inSamFile =~/\.bam$/)
	{
		open ($fin, "samtools view $inSamFile | ") || Carp::croak "cannot open file $inSamFile to read\n";
	}
	else
	{
    	open ($fin, "<$inSamFile") || Carp::croak "cannot open file $inSamFile to read\n";
	}
}


if ( $outSamFile eq "-")
{
     $fout = *STDOUT;
}
else
{
    open ($fout, ">$outSamFile") || Carp::croak "cannot open file $outSamFile to write\n";
}

my $iter = 0;

while (my $line = <$fin>)
{
	chomp $line;

	#head lines
	if ($line=~/^\s*$/)
	{
		print $fout $line, "\n";
		next;
	}
	elsif ($line =~/^\@/)
	{
		print $fout $line, "\n";
		next;
	}

	print STDERR "$iter ...\n" if $verbose && $iter % 50000 == 0;
	$iter++;

	my $sam = lineToSam ($line);
	my $bed = samToBed ($sam, 1);
	#use RNA strand

	next unless $bed; #no alignment

	if ($printUniqOnly)
	{
		next unless $sam->{"TAGS"}=~/XT:A:U/;
	}

	my $flagInfo = $bed->{"flagInfo"};	
	next if $flagInfo->{'query_nomap'};

	if ($bed->{'blockCount'} > 1)
	{	
		my $chrom = $bed->{'chrom'};
		my $chromStart = $bed->{'chromStart'};
		my $strand = $bed->{'strand'};

		for (my $i = 0; $i < $bed->{'blockCount'}-1; $i++)
		{
			my $intronStart = $chromStart + $bed->{'blockStarts'}[$i] + $bed->{'blockSizes'}[$i];
			my $intronEnd = $chromStart + $bed->{'blockStarts'}[$i+1] - 1;

			($intronStart, $intronEnd) = ($intronEnd, $intronStart) if $strand eq '-';			

			my $donorKey = genKey ($chrom, $intronStart, $strand);
			if (exists $spliceSiteHash{'donor'}{$donorKey})
			{
				my %sam2 = %$sam;
				$sam2{'TAGS'} .= "\tSS:Z:" .$spliceSiteHash{'donor'}{$donorKey}{'name'};
				print $fout samToLine (\%sam2), "\n";
			}
			
			my $acceptorKey = genKey ($chrom, $intronEnd, $strand);
			if (exists $spliceSiteHash{'acceptor'}{$acceptorKey})
			{
				my %sam2 = %$sam;
				$sam2{'TAGS'} .= "\tSS:Z:".$spliceSiteHash{'acceptor'}{$acceptorKey}{'name'};
				print $fout samToLine (\%sam2), "\n";
			}
			
		}
	}
}

print STDERR "Done! Totally $iter lines processed! \n" if $verbose;

close ($fin) if $inSamFile ne '-';
close ($fout) if $outSamFile ne '-';

sub genKey
{
    my ($chrom, $pos, $strand) = @_;
    return "$chrom:$pos:$strand";
}

