#!/usr/bin/env perl

use strict;
use warnings;


use Carp;
use Getopt::Long;
use File::Basename;
use List::Util qw(min);
use MyConfig;


=header CHANGE HISTORY

02/16/2016 CZ: added group X vs. nonX comparison. cleaned the code
02/13/2015 CZ: added options to allow specification of AS type and AS test

=cut


# use Parallel::ForkManager;

my $start_run = time();

my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $cache = getDefaultCache ($prog);

my $type = "expr";
my $comparison = "pairwise"; # or 'other'
my $goConf = "";
my $base = "";
my $outDir = "";

#common parameters
my $pValue = 0.01;

#AS specific parameters
my $id2gene2symbolFile = "";
my $ASType = "cass";
my $ASTest = "fisher";
my $ASSubsetFile = "";
my $coverage = 20;
my $dI = 0.2;

#expression specific parameters
my $fc = 1;

my $geneCol = -1;
my $header = 0;

my $skipFail = 0;
my $negative = 0;
my $medianOut = "";
my $updnOut = "";

my $verbose = 0;

my @medians = ();

GetOptions (
		'type:s'=>\$type,
		'comp:s'=>\$comparison,
#		'go:s'=>\$goConf,
		'base:s'=>\$base,
		'out-dir:s'=>\$outDir,
		'id2gene2symbol=s'=>\$id2gene2symbolFile,
		'p|p-value=f'=>\$pValue,
		'f|FC=f'=>\$fc,
		'c|coverage=f'=>\$coverage,
		'i|inclusion=f'=>\$dI,
		'AS-type:s'=>\$ASType,
		'AS-test:s'=>\$ASTest,
		'AS-subset:s'=>\$ASSubsetFile,
		'cache:s'=>\$cache,
		'v|verbose'=>\$verbose
) or die "Error: command line arguments\n";

if (@ARGV != 1)
{
	print "Perform pairwise differential expression or differential splicing analysis from countit data.\n";
	print "USAGE: $prog [options] <input.conf>\n";
	print "\n";
	print "FILES:\n";
	print " input.conf  [string] : Col1 = expr/AS file; Col2 = sample name; tab-delimted\n";
	print "\n";
	print "OPTIONS:\n";
	print " -type       [string] : experiment type ([expr]|as)\n";
	print " -comp       [string] : comparison ([pairwise]|other)\n";
	print " -base       [string] : base directory for files in input-conf\n";
	print " --out-dir   [string] : output dir (default=current directory)\n";
	print "\n";
	print "Splicing-specific options:\n";
	print " --id2gene2symbol [string]: mapping of AS id to gene id and symbol\n";
	print " --AS-type   [string] : AS type ([cass]|taca|alt5|alt3|mutx|iret|apat|apa|snv)\n";
	print " --AS-test   [string] : statistical test ([fisher]|chisq|g|glm)\n";
	print " --AS-subset [string] : list of AS events used for further filtering\n";
	print "\n";
	print "Filtering - Expression and splicing\n";
	print " -p          [float]  : p-value (or FDR) below this value [0.01]\n";
	print "\n";
	print "Filtering - Expression-specific\n";
	print " -f          [float]  : fold-change / log fold-change above this value [1]\n";
	print "\n";
	print "Filtering - Splicing-specific\n";
	print " -c          [float]  : coverage above this value [20]\n";
	print " -i          [float]  : abs(dI) above this value [0.2]\n";
	print "\n";
	exit (1);
}

my ($inConf) = @ARGV;

my $outBase = basename($inConf);
$outBase =~ s/\.\w+$//;

my $verboseFlag = $verbose ? '-v' : '';


#TODO: add other types of AS
Carp::croak "only cass is implemented for AS analysis now\n" if $ASType ne 'cass';


if ($outDir ne '')
{
	runBashCommand ("mkdir $outDir");
	#$outBase = "$outDir/$outBase";
}
else
{
	$outDir = ".";
}

my @infiles = ();
my @groups = ();

print "### Cache dir: $cache\n\n";

Carp::croak "$cache already exists\n" if -d $cache || -f $cache;
my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir at $cache\n" unless $ret == 0;

my $expression = 1;
$expression = 0 if ($type =~ /^as$/i);


my $groups = readConfigFile ($inConf, $base);
my $iter = 0;
my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;

foreach my $gName (@groupNames)
{
	my $samples = $groups->{$gName}->{"samples"};
	foreach my $s (@$samples)
	{
		print "$iter: group=$gName, sample=$s\n" if $verbose;
		my $inputFile = $s;
		$inputFile = "$base/$inputFile" if $base ne '';
		
		print join ("\t", $inputFile, $gName), "\n" if $verbose;		
		$iter++;
	}
}
print "$iter samples in total.\n" if $verbose;

my $totalSampleNum = $iter;

print "generating configuration files for each comparison...\n" if $verbose;
my @comparisonInfo;

for (my $i = 0; $i < scalar(@groupNames); $i++)
{
	my $group1 = $groupNames[$i];
	if ($comparison ne 'pairwise')
	{
		#a or non-a
		my $group2 = "Non-" . $group1;
		my $outBasePair = "$outBase.$group1-$group2";
		my $outConfPair = "$outDir/$outBasePair.conf";
		open (OUTCONF, ">$outConfPair") || Carp::croak "cannot open file $outConfPair to write\n";
		
		my $samples = $groups->{$group1}->{"samples"};
		foreach my $s (@$samples)
		{
			print OUTCONF "$s\t$group1\n";
		}

		for (my $j = 0; $j < scalar(@groupNames); $j++)
		{
			next if $i == $j;
			my $other = $groupNames[$j];
			my $samples = $groups->{$other}->{"samples"};
			foreach my $s (@$samples)
			{
				print OUTCONF "$s\t$group2\n";
			}
		}
		close (OUTCONF);
		push @comparisonInfo, {
			group1=> $group1,
			group2=> $group2,
			group1SampleNum => scalar (@$samples),
			group2SampleNum => $totalSampleNum - scalar (@$samples),
			conf=>$outConfPair};
	}
	else
	{
		#pairwise
		for (my $j = $i+1; $j < scalar(@groupNames); $j++)
		{		
			# make pairwise conf file
			my $group2 = $groupNames[$j];

			my $outBasePair = "$outBase.$group1-$group2";
			my $outConfPair = "$outDir/$outBasePair.conf";

			open (OUTCONF, ">$outConfPair") || Carp::croak "cannot open file $outConfPair to write\n";
			my $samples = $groups->{$group2}->{"samples"};
			my $group2SampleNum = @$samples;
			foreach my $s (@$samples)
			{
				print OUTCONF "$s\t$group2\n";
			}

			$samples = $groups->{$group1}->{"samples"};
			my $group1SampleNum = @$samples;
			foreach my $s (@$samples)
			{
				# print OUTCONF "$base/" if $base;
				print OUTCONF "$s\t$group1\n";
			}
			close (OUTCONF);
			push @comparisonInfo, {
				group1=> $group1,
				group2=> $group2,
				group1SampleNum=> $group1SampleNum,
				group2SampleNum=> $group2SampleNum,
				conf=>$outConfPair};
		}
	}
}

print "performing comparisons ...\n" if $verbose;
foreach my $comp (@comparisonInfo)
{
	my $spliceFile = "";
	my $spliceFile = ".gene" unless $expression;

	my $group1 = $comp->{'group1'};
	my $group2 = $comp->{'group2'};
	my $group1SampleNum = $comp->{'group1SampleNum'};
	my $group2SampleNum = $comp->{'group2SampleNum'};
	my $outConfPair = $comp->{'conf'};
	my $outBasePair = "$outBase.$group1-$group2";
	
	print "$group1 ($group1SampleNum samples) vs. $group2 ($group2SampleNum samples)\n" if $verbose;
	
	if ($expression)
	{

		
		my $disp = $group1SampleNum < 2 || $group2SampleNum < 2 ? "common" : "tagwise";
		print "Using $disp dispersion\n" if $verbose;
		
		print "Testing pairwise differential expression ($group1 vs $group2)...\n" if $verbose;
		
		runBashCommand ("perl $cmdDir/test_expr_diff.pl $verboseFlag -rpkm -base $base -disp $disp -MAplot $outDir/$outBasePair.pdf $outConfPair $outDir/$outBasePair.txt");

		# filter results
		runBashCommand ("perl $cmdDir/filterExprAS.pl -e 0 -p $pValue -q 5 -f $fc -g 2 -r 6,7 -h 1 -n -m $outDir/$outBasePair.medGene.txt -u $outDir/$outBasePair.updn.txt $outDir/$outBasePair.txt $outDir/$outBasePair.filter.txt");

		runBashCommand ("cut -f1 $outDir/$outBasePair.up.txt >> $cache/$outBase.diffExprGenes.all.txt");
		runBashCommand ("cut -f1 $outDir/$outBasePair.dn.txt >> $cache/$outBase.diffExprGenes.all.txt");
		runBashCommand ("cat $outDir/$outBasePair.medGene.txt >> $cache/$outBase.medGene.all.txt");
		# exec '/bin/sleep', $p or die;
	}
	else
	{
		print "Testing pairwise differential splicing...\n" if $verbose;
		my $ASTestFlag = $ASTest ne '' ? "-test $ASTest" : '';
		my $ASTypeFlag = $ASType ne '' ? "-type $ASType" : '';

		my $cmd = "perl $cmdDir/test_splicing_diff.pl $verboseFlag -base $base $ASTypeFlag $ASTestFlag --min-cov $coverage --id2gene2symbol $id2gene2symbolFile $outConfPair $outDir/$outBasePair.txt";
		#right now the id2gene2symbol file is mandatory, which needs to be fixed			

		runBashCommand($cmd);

		# filter results
		runBashCommand("perl $cmdDir/filterExprAS.pl -e 4 -p $pValue -q 14 -c $coverage -d 9 -i $dI -j 12 -n -u $outDir/$outBasePair.UPDN.txt -h 1 $outDir/$outBasePair.txt $outDir/$outBasePair.filter.txt");
		runBashCommand("perl $cmdDir/filterExprAS.pl -e 4 -c $coverage -d 9 -h 1 $outDir/$outBasePair.txt $cache/$outBasePair.mincov.txt");
		runBashCommand("awk '\$28==1' $cache/$outBasePair.mincov.txt | cut -f1 | cut -f1 -d '/' >> $cache/$outBase.minCov.gene.txt");
		#note that the column number might change if the id2gene2symbol is not provided, so the hardcoded column number needs to be fixed
		#$NF will get the last column

		# for i in `ls *mincov20*`; do echo $i; awk '$28==1' $i | cut -f1 | cut -f1 -d "/" >> allgenes.mincov20.txt; done

		if ($ASSubsetFile ne '' && (-f $ASSubsetFile))
		{
			runBashCommand ("mv $outDir/$outBasePair.up.txt $outDir/$outBasePair.up.complete.txt");
			$cmd = "perl $cmdDir/selectRow.pl -h -q 4 $outDir/$outBasePair.up.complete.txt $ASSubsetFile > $outDir/$outBasePair.up.txt";
			runBashCommand($cmd);

			runBashCommand ("mv $outDir/$outBasePair.dn.txt $outDir/$outBasePair.dn.complete.txt");
			$cmd = "perl $cmdDir/selectRow.pl -h -q 4 $outDir/$outBasePair.dn.complete.txt $ASSubsetFile > $outDir/$outBasePair.dn.txt";
			runBashCommand($cmd);
		}

		runBashCommand ("cut -f5 $outDir/$outBasePair.up.txt | tail -n +2 >> $cache/$outBase.diffExprExons.all.txt");
		runBashCommand ("cut -f5 $outDir/$outBasePair.dn.txt | tail -n +2 >> $cache/$outBase.diffExprExons.all.txt");

		runBashCommand ("cut -f1 $outDir/$outBasePair.up.txt | tail -n +2 | cut -f1 -d '/' >> $outDir/$outBasePair.up.gene.txt");
		runBashCommand ("cut -f1 $outDir/$outBasePair.dn.txt | tail -n +2 | cut -f1 -d '/' >> $outDir/$outBasePair.dn.gene.txt");
	}
}

if ($expression)
{
	# unique genes expressed above median in at least one condition
	runBashCommand("sort -u $cache/$outBase.diffExprGenes.all.txt > $outDir/$outBase.diffExprGenes.all.uniq.txt");
	runBashCommand("sort -u $cache/$outBase.medGene.all.txt > $outDir/$outBase.medGene.all.uniq.txt");
}
else
{
	# unique exons meeting minimum coverage in at least one condition
	runBashCommand("sort -u $cache/$outBase.minCov.gene.txt > $outDir/$outBase.minCov.gene.uniq.txt");

	# unique diff expr exons
	runBashCommand("sort -u $cache/$outBase.diffExprExons.all.txt > $outDir/$outBase.diffExprExons.all.uniq.txt");
}

print "generate summary file ...\n" if $verbose;

my $fout;
open ($fout, ">$outDir/$outBase.summary.txt") || Carp::croak "cannot open file $outDir/$outBase.summary.txt to write\n";

if ($comparison eq 'pairwise')
{
	my @summary;
	for (my $i = 0; $i < @groupNames; $i++)
	{	
		my $group1 = $groupNames[$i];
		for (my $j = 0; $j < @groupNames; $j++)
		{
			if ($i == $j)
			{
				$summary[$i][$j] = 0;
				next;
			}

			my $group2 = $groupNames[$j];
			
			my $resultFile = $j>$i ? "$outDir/$outBase.$group1-$group2.up.txt" : "$outDir/$outBase.$group2-$group1.dn.txt";
			my $cmd = "wc -l $resultFile";
			my $n = `$cmd`;
			chomp $n;
			$n=~/^(\d+)/;
			$n = $1 -1;
			$summary[$i][$j] = $n;
		}
	}
	
	print $fout join ("\t", "group", @groupNames), "\n";
	for (my $i = 0; $i < @groupNames; $i++)
	{
		print $fout join ("\t", $groupNames[$i], @{$summary[$i]}), "\n";
	}
}
else
{
	print $fout join ("\t", "comparison", "up", "down"), "\n";
	for (my $i = 0; $i < @groupNames; $i++)
    {
        my $group1 = $groupNames[$i];
		my $group2 = "Non-" . $group1;
		my $outBasePair .= "$outBase.$group1-$group2";
		
		my @comp = ("up", "dn");	
		my @stat;	
		foreach my $c (@comp)
		{
			my $resultFile = "$outDir/$outBasePair.$c.txt";
			my $cmd = "wc -l $resultFile";
            my $n = `$cmd`;
			chomp $n;
			$n=~/^(\d+)/;
			$n = $1 -1;
			push @stat, $n;
		}
		print $fout join ("\t", "$group1-$group2", @stat), "\n";
	}
}
close ($fout);



# Genes that were differentially expressed in at least one pair of conditions:

# 	rm diffExprGenes.all.txt
# 	for i in `ls *filter.txt`; do tail -n+2 $i | awk '{if ($9!=0) {print $1}}' >> diffExprGenes.all.txt; done
# 	sort -u diffExprGenes.all.txt > diffExprGenes.all.uniq.txt

# Expr matrix
	
# Filter gene expression matrix:

# 	# word boundaries for grep
# 	sed -re 's/([[:digit:]]+)/^\1\\b/' diffExprGenes.all.uniq.txt > diffExprGenes.all.uniq.grep.txt
# 	grep -f pairwise/diffExprGenes.all.uniq.grep.txt DevCortex.expr.ab.txt > DevCortex.expr.ab.filter.txt

# 	# header row
# 	cat <(head -1 DevCortex.expr.ab.txt) DevCortex.expr.ab.filter.txt > tmp.txt
# 	mv tmp.txt DevCortex.expr.ab.filter.txt

# Filtered correlation matrix

# 	dcAS <- read.delim("DevCortex.expr.ab.filter.txt")
# 	write.table(cor(dcAS[,4:21], use="complete.obs"), file="DevCortex.expr.ab.filter.cor.txt",sep="\t",col.names=TRUE,row.names=TRUE)


my $end_run = time();
my $run_time = $end_run - $start_run;
print "Runtime: $run_time seconds\n";



sub runBashCommand
{
	my ($cmd) = @_;
	# my $ret = system ($cmd);
	# c flag allows commands with process substitution to be run
	my @args = ( "bash", "-c", $cmd );
	my $ret = system(@args);
	Carp::croak "cmd=$cmd failed\n" unless $ret == 0;
	print "cmd=$cmd\n";
	return $ret;	
}




sub readConfigFile
{
	my ($configFile, $base) = @_;
	my $fin;
	open ($fin, "<$configFile") || Carp::croak "cannot open file $configFile to read\n";
	my $i = 0;
	my %groups;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		next if $line=~/^\#/;
		my ($sampleName, $groupName) = split (/\t/, $line);
		Carp::croak "ill-formed conf file $inConf\n" if ($groupName eq "" || $sampleName eq "");
		$groups{$groupName}->{"id"} = $i++ unless exists $groups{$groupName};
		push @{$groups{$groupName}->{"samples"}}, $sampleName;
	}
	close ($fin);
	return \%groups;
}

