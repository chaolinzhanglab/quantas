#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Carp;
use File::Basename;

use MyConfig;
use Bed;


my $prog = basename ($0);

my $verbose = 0;
my $exonBedFile = "";
my $intronBedFile = "";

my $cache = getDefaultCache ($prog);

GetOptions ("e=s"=>\$exonBedFile,
			"j=s"=>\$intronBedFile,
			"c:s"=>\$cache,
			"v"=>\$verbose);

if (@ARGV != 4)
{
	print "associate exonic SNV to splice sites\n";
	print "Usage: $prog [options] <in.snv.bed> <out.donor.bed> <out.acceptor.bed> <snv2splicesite.txt>\n";
	print " -e     [string]: exon bed file\n";
	print " -j     [string]: junction bed file\n";
	print " -v             : verbose\n";
	exit (1);
}

my ($snvBedFile, $donorBedFile, $acceptorBedFile, $snv2spliceSiteFile) = @ARGV;


my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir $cache\n" if $ret != 0;

print "loading introns ...\n" if $verbose;

my $fin;

my %donorHash;
my %acceptorHash;
open ($fin, "<$intronBedFile") || Carp::croak "cannot open file $intronBedFile to read\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	next if $line =~/^track/;
	
	my $intron = lineToBed ($line);
	my $chromStart = $intron->{'chromStart'};
	my $chromEnd = $intron->{'chromEnd'};

	my $donor = {chrom=>$intron->{'chrom'}, 
			chromStart=>$intron->{'chromStart'},
			chromEnd=>$intron->{'chromStart'},
			name=>$intron->{'name'},
			score=>$intron->{'score'},
			strand=>$intron->{'strand'}};
	my $acceptor = {chrom=>$intron->{'chrom'}, 
			chromStart=>$intron->{'chromEnd'},
			chromEnd=>$intron->{'chromEnd'},
			name=>$intron->{'name'},
			score=>$intron->{'score'},
			strand=>$intron->{'strand'}};

	($donor, $acceptor) = ($acceptor, $donor) if $intron->{'strand'} eq '-';
	
	$donor->{'name'} .= "[5SS]";
	$acceptor->{'name'} .= "[3SS]";

	my $donorKey = genKey ($donor->{'chrom'}, $donor->{'chromStart'}, $donor->{'strand'});
	my $acceptorKey = genKey ($acceptor->{'chrom'}, $acceptor->{'chromStart'}, $acceptor->{'strand'});


	$donorHash{$donorKey} = $donor unless exists $donorHash{$donorKey};
	$acceptorHash{$acceptorKey} = $acceptor unless exists $acceptorHash{$acceptorKey};
}

close ($fin);



print "get overlapping exons for each SNV ...\n" if $verbose;
my $snv_vs_exonFile = "$cache/snv_vs_exon.bed";

my $verboseFlag = $verbose ? "-v" : "";

my $cmd = "perl ~/scripts/tagoverlap.pl $verboseFlag -ss -region $snvBedFile $exonBedFile $snv_vs_exonFile";
print $cmd, "\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;


print "loading exon annotations for SNVs ...\n" if $verbose;

open ($fin, "<$snv_vs_exonFile") || Carp::croak "cannot open file $snv_vs_exonFile\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	my $bed = lineToBed ($line);
	my $name = $bed->{'name'};
	my ($exonId, $snvId) = split ("//", $name);
	my $chrom = $bed->{'chrom'};
	my $exonStart = $bed->{'chromStart'};
	my $exonEnd = $bed->{'chromEnd'};
	my $strand = $bed->{'strand'};

	my $donorPos = $exonEnd + 1;
	my $acceptorPos = $exonStart - 1;

	if ($strand eq '-')
	{
		($donorPos, $acceptorPos) = ($acceptorPos, $donorPos);
	}

	my $donorKey = genKey ($chrom, $donorPos, $strand);
	my $acceptorKey = genKey ($chrom, $acceptorPos, $strand);

	push @{$donorHash{$donorKey}{'snv'}}, $snvId if exists $donorHash{$donorKey};
	push @{$acceptorHash{$acceptorKey}{'snv'}}, $snvId if exists $acceptorHash{$acceptorKey};
}
close ($fin);

print "writing donors associated with SNVs ...\n" if $verbose;

my ($fout, $fout2);

open ($fout2, ">$snv2spliceSiteFile") || Carp::croak "cannot open file $snv2spliceSiteFile to write\n";

open ($fout, ">$donorBedFile") || Carp::croak "cannot open file $donorBedFile to write\n";

foreach my $donorKey (sort keys %donorHash)
{
	next unless exists $donorHash{$donorKey}{'snv'};
	print $fout bedToLine ($donorHash{$donorKey}), "\n";
	foreach my $snvId (@{$donorHash{$donorKey}{'snv'}})
	{
		print $fout2 join ("\t", $snvId, $donorHash{$donorKey}{'name'}, "5'SS"), "\n";
	}
}

close ($fout);

open ($fout, ">$acceptorBedFile") || Carp::croak "cannot open file $acceptorBedFile to write\n";

foreach my $acceptorKey (sort keys %acceptorHash)
{
	next unless exists $acceptorHash{$acceptorKey}{'snv'};
	print $fout bedToLine ($acceptorHash{$acceptorKey}), "\n";

	foreach my $snvId (@{$acceptorHash{$acceptorKey}{'snv'}})
    {
        print $fout2 join ("\t", $snvId, $acceptorHash{$acceptorKey}{'name'}, "3'SS"), "\n";
    }
}

close ($fout);
close ($fout2);






sub genKey
{
	my ($chrom, $pos, $strand) = @_;
	return "$chrom:$pos:$strand";
}





