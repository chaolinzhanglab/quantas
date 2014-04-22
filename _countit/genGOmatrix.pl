use strict;

#TODO: parallelize the hypergeometric calculations

use Carp;
use Getopt::Long;
use File::Basename;
use List::Util qw(min);
# use Math::BigFloat;
# use Parallel::ForkManager;

my $prog = basename ($0);

my $geneCol = 0;
# my $pValue = 0;
# my $pValueCol = -1;
# my $fc = 0;
# my $fcCol = -1;
# my $expr = -1;
# my $exprCol = -1;
# my $coverage = 0;
# my $coverageCol = -1;
# my $dI = 0;
# my $dICol = -1;
my $header = 0;
my $conf = 0;
my $sampleIn = "";
my $bgFile = "";
my $log = 0;

# store calculated log-factorials
# my %logfactStore = (0=>Math::BigFloat->new("0"));
my %logfactStore = (0=>0);

# for (my $i = 1; $i < 1000; $i++)
# {
# 	print log($i) . "\n";
# 	print log($i) . "\n";
# 	print log($i) . "\n";
# 	sleep(1);
# }



# GetOptions('coordinates=f{2}' => \@coor, 'rgbcolor=i{3}' => \@color);

GetOptions (
		'e|event-col=i'=>\$geneCol,
		# 'p|p-value=f'=>\$pValue,
		# 'q|p-value-col=i'=>\$pValueCol,
		# 'f|FC=f'=>\$fc,
		# 'g|FC-col=i'=>\$fcCol,
		# 'r|expression=i'=>\$expr,
		# 'c|coverage=f'=>\$coverage,
		# 'd|coverage-col=i'=>\$coverageCol,
		# 'i|inclusion=f'=>\$dI,
		# 'j|inclusion-col=i'=>\$dICol,
		'h|header-lines=i'=>\$header,
		'c|conf'=>\$conf,
		's|sample-name=s'=>\$sampleIn,
		'l|log'=>\$log
) or die "Error: command line arguments\n";

if (@ARGV != 4)
{
	print "Generate gene ontology matrix\n";
	print "Usage: $prog [options] <GO-file> <background-file> <input-file [.conf, see -c flag]> <output-file>\n";
	print "Formats (tab-delimited, header optional):\n";
	print " GO-file     : first column=event ID, second column=GO term\n";
	print " background  : first column=event ID\n";
	print " input files : first column=event ID; pre-filtered (e.g. only overexpressed genes)\n";
	print " input conf  : first column=file, second column(optional)=sample name; sorted\n";
	print "OPTIONS:\n";
	print " -e [int]    : event (gene/exon) ID column (0-based) [0]\n";
	# print " -p [float]  : p-value (or FDR) below this value\n";
	# print " -q [int]    : p-value column (0-based)\n";
	# print " -f [float]  : fold-change / log fold-change above this value\n";
	# print " -g [int]    : fold-change column (0-based)\n";
	# print " -r [int]    : column specifying if expression (RPKM) is above median [spec. by 1/0]. Can also use this for arbitrary filtering.\n";
	# print " -s [int]    : expression column (0-based)\n";	
	# print " -c [float]  : coverage above this value\n";
	# print " -d [int]    : coverage column (0-based)\n";	
	# print " -i [float]  : abs(dI) above this value\n";
	# print " -j [int]    : dI column (0-based)\n";
	print " -h [int]    : number of header lines [default = 0]\n";
	print " -c          : if specified, treat input file parameter as a configuration file. Otherwise, treat as input file.\n";
	print " -s [string] : optional sample name when single input file is given (not conf file)\n";
	print " -l          : display log10 values [off]\n";
	exit (1);
}

# die "Error: event column required\n" if ($geneCol < 0);
# die "Error: p-value column required\n" if ($pValue && $pValueCol < 0);
# die "Error: fold-change column required\n" if ($fc && $fcCol < 0);
# die "Error: expression column required\n" if ($expr && $exprCol < 0);
# die "Error: coverage column required\n" if ($coverage && $coverageCol < 0);
# die "Error: dI column required\n" if ($dI && $dICol < 0);
die "Error: Cannot specify sample name for configuration file. Specify names in conf file instead.\n" if ($sampleIn && $conf);

my ($GOFile, $bgFile, $inputConf, $outfile) = @ARGV;

open (GOFILE, "<$GOFile") || Carp::croak "cannot open $GOFile\n";
open (BGFILE, "<$bgFile") || Carp::croak "cannot open $bgFile\n";
open (INCONF, "<$inputConf") || Carp::croak "cannot open $inputConf\n";
open (OUTFILE, ">$outfile") || Carp::croak "cannot open $outfile for writing\n";

my @infiles = ();
my @sampleNames = ();

# if a single sample is given
if ($conf==0)
{
	@infiles = ($inputConf);
	@sampleNames = ($sampleIn);
}
# otherwise go through samples in conf file
else
{
	while (my $line = <INCONF>)
	{
		chomp $line;
		my ($infile, $sampleName) = split(/\t/, $line);
		push (@infiles, $infile);
		push (@sampleNames, $sampleName);
		
	}
	close(INCONF);
}

# check if files exists and assign ID to sample names for sorting
for (my $i = 0; $i < scalar(@infiles); $i++)
{
	my $infile = $infiles[$i];
	my $sampleName = $sampleNames[$i];

	open (INFILE, "<$infile") || Carp::croak "cannot open $infile\n";
	close (INFILE);

	if ($sampleName eq "")
	{
		# remove path and extension
		$sampleName = basename ($infile);
		$sampleName =~ s/\..*$//;
	}

	$sampleName = sprintf("%03d", $i+1) . ":$sampleName";
	$sampleNames[$i] =  $sampleName;
}

my %GOterms;
my %GOgenes;

# my %allGenes;
my %bgGenes;

while (my $line = <BGFILE>)
{
	chomp $line;
	my ($gene) = split(/\t/, $line);

	$bgGenes{$gene} = 1;
}
close (BGFILE);
print scalar(keys %bgGenes) . " events in background.\n";

# go through GO file
print "Processing GO file...\n";
while (my $line = <GOFILE>)
{
	chomp $line;
	my ($gene, $term) = split(/\t/, $line);

	$GOterms{$term}{$gene} = 1 if $bgGenes{$gene};
	$GOgenes{$gene}{$term} = 1 if $bgGenes{$gene};
}
close (GOFILE);
print scalar(keys %GOterms) . " GO terms matching background found.\n";



my %results;
my %sampleGenes;

my $countSamples = 0;
# read and process input files
# foreach my $infile (@infiles)
for (my $i = 0; $i < scalar(@infiles); $i++)
{
	my $infile = $infiles[$i];
	my $sampleName = $sampleNames[$i];
	print "Processing $sampleName...\n";

	# my %input = {};
	open (INFILE, "<$infile") || Carp::croak "cannot open file $infile to read\n";
	my $countLines = 0;
	while (my $line = <INFILE>)
	{
		$countLines++;
		next if ($countLines <= $header);

		chomp $line;
		my @splits = split(/\t/, $line);

		# add genes to index
		my $gene = $splits[$geneCol];
		# $allGenes{$gene} = 1 if (!exists($allGenes{$gene}));

		# check filters
		# next if ($pValue > 0 && $splits[$pValueCol] >= $pValue);
		# next if ($fc > 0 && abs($splits[$fcCol]) <= $fc);
		# next if ($expr > 0 && $splits[$expr] == 0);
		# next if ($coverage > 0 && $splits[$coverageCol] <= $coverage);
		# next if ($dI > 0 && $splits[$dICol] <= $dI);

		# $input{$gene} = 1;
		$sampleGenes{$sampleName}{$gene} = 1;
		# print "$gene\n";
	}
	close(INFILE);
	# die;
}

# print scalar (keys %allGenes) . " unique genes found.\n";

print "Calculating hypergeometric scores...\n";
foreach my $sample (sort keys %sampleGenes)
{
	print "$sample\n";
	# calculate hypergeometric score for each GO term
	my $countGO = 0;
	foreach my $go (sort keys %GOterms)
	{
		$countGO++;
		print ".";
		print "$countGO GO terms\n" if ($countGO%100 == 0);
		$results{$go}{$sample} = 1;

		# Numbers needed for hypergeometric test
		## a = population size (total background genes)
		## b = subset size (filtered genes)
		## c = interesting in population (# GO genes of specific term)
		## d = interesting in subset (# GO genes of specific term in filtered)

		my $a = scalar(keys %bgGenes);		
		my $b = 0;
		my $c = 0;
		my $d = 0;
		foreach my $gene (keys %bgGenes)
		{
			$b++ if ($sampleGenes{$sample}{$gene});
			$c++ if ($GOterms{$go}{$gene});
			$d++ if ($sampleGenes{$sample}{$gene} && $GOgenes{$gene}{$go});
		}

		my @hyperParams = ($c, $a-$c, $b, $d);
		# my @hyperParams = (300, 700, 100, 40);
		my $hyper2 = phyper(@hyperParams);
		# my $hyper2 = logfact($a);



		### OLD USING PHYPER (in R) ###

		# Use phyper (R)
		# Fast function, but it's slow to call R everytime (~1 call per second)
		## phyper(cool balls drawn, cool balls in urn, uncool balls in urn, balls drawn, lower.tail?)
		## Note:
		### lower.tail: if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
		### -> we need to subtract 1 from d, since P[X>d-1] = P[X≥d]
		
		# my $hyper = `R --slave -q -e "cat(phyper($d-1, $c, $a-$c, $b, lower.tail=FALSE))"`;
		# print "$hyper $hyper2 " . ($hyper2-$hyper) . "\n";

		# print "$a\t$b\t$c\t$d\t$hyper2\n";

		$hyper2 = log($hyper2)/log(10) if ($log);

		$results{$go}{$sample} = $hyper2;
	}
	print "\n";
}

# OUTFILE already opened above
my $out = "GO term";
foreach my $sample (sort keys %sampleGenes)
{
	$sample =~ s/^\d+://;
	$out .= "\t$sample";
}
print OUTFILE "$out\n";

foreach my $go (keys %GOterms)
{
	my $out = "$go";
	foreach my $sample (sort keys %sampleGenes)
	{
		$out .= "\t$results{$go}{$sample}";
	}
	print OUTFILE "$out\n";
}

close (OUTFILE);

# http://www.perlmonks.org/?node_id=466599
# http://www.perlmonks.org/?node_id=856885

# hypergeometric probability
sub dhyper {
	# There are n "good" and m "bad" balls in an urn.
	# Pick N of them. The probability of i selections:
	# (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
	my ($n, $m, $N, $i) = @_;
		
	my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
	my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
	
	return exp($loghyp1 - $loghyp2);
}

# cumulative hypergeometric probability
sub phyper {
	# There are n "good" and m "bad" balls in an urn.
	# Pick N of them. The probability of i or more successful selections (up to N, unless n < N): 
	my ($n, $m, $N, $i) = @_;

	my $total = 0;
	for (my $j = $i; $j <= min($n,$N); $j++)
	{
		$total += dhyper($n,$m,$N,$j);
	}

	$total = 1 if $total > 0.99999;

	return $total;
}

# calculates log factorials
# recursive because we might as well store everything without needing to recalculate every one of the previous values (e.g. we calculate logfact(2) when getting logfact(3))
sub logfact {
	# my $num = Math::BigFloat->new($_[0]);
	my $num = $_[0];
	# don't recalculate if we already have it
	return $logfactStore{$num} if exists($logfactStore{$num});

	# recursive: ln (current * (current-1)!) = ln (current) + lnfact (current-1)
	# my $logfact = Math::BigFloat->new();
	my $logfact = log($num) + logfact($num-1);
	$logfactStore{$num} = $logfact;
	return $logfact;
}
