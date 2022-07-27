#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

use MyConfig;
use Bed;


my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $verbose = 0;
my $big = 0;

my $cache = getDefaultCache ($prog);

##
#TODO: strictly speaking, we would like to separate strand for stranded libraries
#however, since we focus on canonical exon junctions, we assume or reads aligned to proper junctions are should be counted, no matter which strand the read is sequenced




GetOptions (
		'c|cache:s'=>\$cache,
		'big'=>\$big,
		'v|verbose'=>\$verbose);

if (@ARGV != 3)
{
	print "score exon pairs/trios by junction reads\n";
	print "Usage: $prog [options] <trio.bed> <tag.junction.bed> <out.bed>\n";
	print " <trio.bed> -- bed file of exon/pairs/trios\n";
	print " <tag.junction.bed> -- junction reads bed file\n";
	print "OPTIONS:\n";
	print " -big           : big tag file\n";
	print " -c             : cache dir ($cache)\n";
	print " -v             : verbose\n";
	exit (0);
}


my ($exonTrioBedFile, $tagJunctionBedFile, $outBedFile) = @ARGV;

my $verboseFlag = $verbose ? "-v" : '';
my $bigFlag = $big ? "--big" : '';

system ("mkdir $cache") unless -d $cache;

######################################################
#handle junction tags
######################################################

print "extract introns from $exonTrioBedFile ...\n" if $verbose;
my $exonTrioIntronBedFile = "$cache/trio.intron.bed";
my $cmd = "perl $cmdDir/gene2ExonIntron.pl $verboseFlag --allow-duplicate-name -oi $exonTrioIntronBedFile $exonTrioBedFile";
my $ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

print "extracting introns from $tagJunctionBedFile ...\n" if $verbose;
my $tagIntronBedFile = "$cache/tag.intron.bed";
$cmd = "perl $cmdDir/gene2ExonIntron.pl $verboseFlag --allow-duplicate-name -oi $tagIntronBedFile -nid $tagJunctionBedFile";
$ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;


print "match tag introns and exon trio introns ...\n" if $verbose;
my $intronMatchFile = "$cache/trio.vs.tag.match.bed";
$cmd = "perl $cmdDir/bedMatch.pl -c $cache/bed_match_tmp $verboseFlag $bigFlag $exonTrioIntronBedFile $tagIntronBedFile $intronMatchFile";
$ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;


print "count the number of tags for each exon-trio intron ...\n" if $verbose;
my $exonTrioJunctionTagCountFile = "$cache/trio.junction.tagcount.txt";

#$cmd = "awk '{print \$4}' $intronMatchFile | awk -F \"//\" '{print \$1}' | sort | uniq -c | awk '{print \$2\"\\t\"\$1}' > $exonTrioJunctionTagCountFile"; #$ASJunctionTagCountFile";
#$ret = system ($cmd);
#Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;
#the lines above were replaced to improve efficiency

my $fin;
open ($fin, "<$intronMatchFile") || Carp::croak "cannot open $intronMatchFile ...\n";

my %etJunctionTagCountHash;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	my @cols = split ("\t", $line, 5);
	my $name = $cols[3];
	my ($intronId, $useless) = split ("//", $name);
	$etJunctionTagCountHash{$intronId}++;
}
close ($fin);

my $fout;
open ($fout, ">$exonTrioJunctionTagCountFile") || Carp::croak "cannot open file $exonTrioJunctionTagCountFile to write\n";

foreach my $intronId (keys %etJunctionTagCountHash)
{
	print $fout $intronId, "\t", $etJunctionTagCountHash{$intronId}, "\n";
}
close ($fout);


print "loading exon trios from $exonTrioBedFile ...\n" if $verbose;

my $exonTrios = readBedFile ($exonTrioBedFile, $verbose);

my $n = @$exonTrios;
print "$n exon trios loaded\n" if $verbose;

my %exonTrioHash;

foreach my $et (@$exonTrios)
{
	$et->{"score"} = 0;
	$et->{"i"} = exists $et->{"blockCount"} ?  $et->{"blockCount"} - 1: 0;
	$exonTrioHash{$et->{"name"}} = $et;
}


print "reading junction tag count file from $exonTrioJunctionTagCountFile ...\n" if $verbose;

open ($fin, "<$exonTrioJunctionTagCountFile") || Carp::croak "cannot open file $exonTrioJunctionTagCountFile to read\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	my ($intronId, $count) = split ("\t", $line);
	$intronId=~/^(.*?)\_\d+$/;
	my $etId = $1;

	my $et = $exonTrioHash{$etId};

	$et->{"score"} += $count / $et->{"i"} if exists $et->{"i"};
}
close ($fin);

print "dump results to bed file $outBedFile ...\n" if $verbose;

writeBedFile ($exonTrios, $outBedFile);

