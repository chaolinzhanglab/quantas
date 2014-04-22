use strict;

use Carp;
use Getopt::Long;
use File::Basename;
use List::Util qw(min);

my $prog = basename ($0);

my $geneCol = -1;
my $pValue = 0;
my $pValueCol = -1;
my $fc = 0;
my $fcCol = -1;
my $exprCol = "";
my @exprCols = ();
my $coverage = 0;
my $coverageCol = -1;
my $dI = 0;
my $dICol = -1;
my $header = 0;

my $skipFail = 0;
my $negative = 0;
my $medianOut = "";
my $updnOut = "";

my $verbose = 0;

my @medians = ();

GetOptions (
		'e|event-col=i'=>\$geneCol,
		'p|p-value=f'=>\$pValue,
		'q|p-value-col=i'=>\$pValueCol,
		'f|FC=f'=>\$fc,
		'g|FC-col=i'=>\$fcCol,
		'r|rpkm-cols=s'=>\$exprCol,
		'c|coverage=f'=>\$coverage,
		'd|coverage-col=i'=>\$coverageCol,
		'i|inclusion=f'=>\$dI,
		'j|inclusion-col=i'=>\$dICol,
		'h|header-lines=i'=>\$header,
		's|skip-fail'=>\$skipFail,
		'n|negative'=>\$negative,
		'm|median-file=s'=>\$medianOut,
		'u|updn-file=s'=>\$updnOut,
		'v|verbose'=>\$verbose
) or die "Error: command line arguments\n";

if (@ARGV != 2)
{
	print "Filter expression or alternative splicing data\n";
	print "\nUSAGE: $prog [options] <input-file> <output-file>\n";
	print "\nFILES (tab-delimted):\n";
	print " INPUT       : custom, tab-delimted file (the columns are specified in the options)\n";
	print " OUTPUT      : DEFAULT              : same as input file, with an additional column specifying whether event meets filtering criteria [1/0/[-1]] (see -n flag)\n";
	print "             : ALTERNATE [-s flag]  : same as input file, but only include events that meet filtering criteria\n";
	print "             : ADDITIONAL [-m flag] : additional file containing events that are expressed above median in at least one sample (for GO background)\n";
	print "             : ADDITIONAL [-u flag] : two files, each the same as output file but only including filtered events whose expression is up or down\n";
	
	print "\nOPTIONS:\n";
	print "\n";
	print " -h [int]    : number of header lines [default = 0]\n";

	print "\n Filtering - Expression and splicing\n";
	print " -e [int]    : event (gene/exon) ID column\n";
	print " -p [float]  : p-value (or FDR) below this value\n";
	print " -q [int]    : p-value column (0-based)\n";
	
	print "\n Filtering - Expression-specific\n";
	print " -f [float]  : fold-change / log fold-change above this value\n";
	print " -g [int]    : fold-change column (0-based)\n";
	print " -r [x,y,z]  : comma-delimited list of columns (0-based, no spaces) containing RPKM values; at least one required to be > median\n";
	# print " -s [int]    : expression column (0-based)\n";

	print "\n Filtering - Splicing-specific\n";
	print " -c [float]  : coverage above this value\n";
	print " -d [int]    : coverage column (0-based)\n";	
	print " -i [float]  : abs(dI) above this value\n";
	print " -j [int]    : dI column (0-based)\n";

	print "\n Printing and output\n";
	print " -s          : skip (do not print) lines that do not pass filtering criteria [default = off]\n";
	print " -n          : print -1 for negative FC/dI [default = off]\n";
	print " -m [str]    : file to output genes which are above median RPKM in at least one sample (requires -e)\n";
	print " -u [str]    : file to output filtered up/dn genes (requires -e). Filename should contain 'UPDN' which will be replaced with 'up' or 'dn'. Turns on -n flag.\n";
	exit (1);
}

$negative = 1 if ($updnOut ne "");

# die "Error: event column required\n" if ($geneCol < 0);
die "Error: p-value column required\n" if ($pValue && $pValueCol < 0);
print "Warning: p-value cutoff not given\n" if ($pValueCol > 0 && $pValue == 0);
die "Error: fold-change column required\n" if ($fc && $fcCol < 0);
print "Warning: fold-change cutoff not given\n" if ($fcCol > 0 && $fc == 0);
# die "Error: expression column required\n" if ($expr && $exprCol < 0);
die "Error: coverage column required\n" if ($coverage && $coverageCol < 0);
print "Warning: coverage cutoff not given\n" if ($coverageCol > 0 && $coverage == 0);
die "Error: dI column required\n" if ($dI && $dICol < 0);
print "Warning: dI cutoff not given\n" if ($dICol > 0 && $dI == 0);
die "Error: event ID column required\n" if ($medianOut ne "" && $geneCol < 0);
die "Error: event ID column required\n" if ($updnOut ne "" && $geneCol < 0);
die "Error: One of FC or dI column required\n" if ($negative && !($dICol < 0 xor $fcCol < 0));
# die "Error: One of FC or dI required\n" if ($updnOut && !($dI xor $fc));
die "Error: Please enter RPKM as a comma-delimited list, e.g. -r 7,8,9,10\n" if ($exprCol ne "" && $exprCol !~ /^(\d,?)+$/);

my ($infile, $outfile) = @ARGV;

if ($exprCol ne "")
{
	my @RPKMs;
	my @exprColSplits = split(/,/, $exprCol);

	# double check
	foreach my $exprColSplit (@exprColSplits)
	{
		push (@exprCols, $exprColSplit) if ($exprColSplit =~ /^\d+$/);
	}

	foreach my $exprCol (@exprCols)
	{
		push(@RPKMs, ())
	}
	
	open (INFILE, "<$infile") || Carp::croak "cannot open $infile\n";
	my $countLines = 0;
	while (my $line = <INFILE>)
	{
		$countLines++;
		next if ($countLines <= $header);

		chomp $line;
		
		my @splits = split(/\t/, $line);

		for (my $i = 0; $i < scalar(@exprCols); $i++)
		{
			push (@{$RPKMs[$i]}, $splits[$exprCols[$i]]);
		}
	}

	close(INFILE);

	for (my $i = 0; $i < scalar(@exprCols); $i++)
	{
		$medians[$i] = median (@{$RPKMs[$i]});
		print "Median $i: $medians[$i]\n" if $verbose;
	}
}

# add UPDN keyword if not present already
$updnOut =~ s/(\.\w+)?$/.UPDN$1/ if ($updnOut ne "" && $updnOut !~ /UPDN/i);
(my $upOut = $updnOut) =~ s/UPDN/up/i;
(my $dnOut = $updnOut) =~ s/UPDN/dn/i;
print "$updnOut\t$upOut\t$dnOut\n" if $verbose;
# die;


open (INFILE, "<$infile") || Carp::croak "cannot open $infile\n";
open (OUTFILE, ">$outfile") || Carp::croak "cannot open $outfile for writing\n";
(open (OUTMED, ">$medianOut" || Carp::croak "cannot open $medianOut for writing")) if $medianOut ne "";
(open (OUTUP, ">$upOut" || Carp::croak "cannot open $upOut for writing")) if $updnOut ne "";
(open (OUTDN, ">$dnOut" || Carp::croak "cannot open $dnOut for writing")) if $updnOut ne "";


my $countLines = 0;
while (my $line = <INFILE>)
{
	$countLines++;
	chomp $line;

	if ($countLines <= $header)
	{
		if ($skipFail)
		{
			$line = "$line\n";
		}
		else
		{
			$line = "$line\tFilter\n";
		}

		print OUTFILE $line;
		print OUTUP $line if $updnOut ne "";
		print OUTDN $line if $updnOut ne "";
		
		next;
	}
	
	my @splits = split(/\t/, $line);

	# check filters
	my $passFilter = 1;
	$passFilter = 0 if ($pValue > 0 && $splits[$pValueCol] >= $pValue);
	$passFilter = 0 if ($fc > 0 && abs($splits[$fcCol]) <= $fc);
	# $passFilter = 0 if ($exprCol > 0 && $splits[$exprCol] == 0);
	$passFilter = 0 if ($coverage > 0 && $splits[$coverageCol] <= $coverage);
	$passFilter = 0 if ($dI > 0 && abs($splits[$dICol]) <= $dI);

	# print "$passFilter\t";
	my $overMedian = 0;
	for (my $i = 0; $i < scalar(@exprCols); $i++)
	{
		$overMedian = 1 if $splits[$exprCols[$i]] > $medians[$i];
	}
	$passFilter = 0 if ($exprCol ne "" && $overMedian == 0);
	print OUTMED "$splits[$geneCol]\n" if ($overMedian == 1 && $medianOut ne "");

	# my $dir = 1;
	# $dir = -1 if (($dI > 0 && $splits[$dICol] < 0) || ($fc > 0 && $splits[$fcCol] < 0));

	$passFilter *= -1 if ($negative && $splits[$dICol] < 0);
	$passFilter *= -1 if ($negative && $splits[$fcCol] < 0);

	next if ($skipFail && !$passFilter);

	if ($skipFail)
	{
		$line = "$line\n";
	}
	else
	{
		$line = "$line\t$passFilter\n";
	}

	print OUTFILE $line;
	print OUTUP $line if ($updnOut ne "" && $passFilter > 0);
	print OUTDN $line if ($updnOut ne "" && $passFilter < 0);
}

close(INFILE);
close (OUTFILE);
close (OUTMED) if $medianOut ne "";
close (OUTUP) if $updnOut ne "";
close (OUTDN) if $updnOut ne "";


# http://stackoverflow.com/questions/5119034/using-perl-to-find-median-mode-standard-deviation

sub median
{
	my @values = @_;
	my $median;
	my $mid = int @values/2;
	my @sorted_values = sort by_number @values;
	if (@values % 2) {
	    $median = $sorted_values[ $mid ];
	} else {
	    $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
	} 
}

sub by_number {
    if ($a < $b){ -1 } elsif ($a > $b) { 1 } else { 0 }
}
