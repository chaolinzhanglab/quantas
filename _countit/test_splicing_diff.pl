#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;

use Common;
use MyConfig;

my $prog = basename ($0);
my $verbose = 0;

my $type = 'cass';
my $base = "";

my $minCoverage = 20;
my $FDR = 1;
my $deltaI = 0;
my $filterOutput = 0;

my $id2gene2symbolFile = "";

my $cache = getDefaultCache ($prog);

GetOptions ("t|type:s"=>\$type,
	"base:s"=>\$base,
	"min-cov:i"=>\$minCoverage,
	"FDR:f"=>\$FDR,
	"dI:f"=>\$deltaI,
	"id2gene2symbol:s"=>\$id2gene2symbolFile,
	"filter-output"=>\$filterOutput,
	"c:s"=>\$cache,
	"v|verbose"=>\$verbose
);


if (@ARGV != 2)
{
	print "statistical test of differential splicing\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " -base            [string] : base dir of input data\n";
	print " -type            [string] : AS type ($type)\n";
	print " --min-cov        [int]    : min coverage to calculate FDR ($minCoverage)\n";
	print " --id2gene2symbol [file]   : mapping file of id to gene to symbol\n";
	print " --filter-output           : filter output by FDR and dI\n";
	print " -FDR             [float]  : threshold on FDR ($FDR)\n";
	print " -dI              [float]  : threshold on dI ($deltaI)\n";
	print " -c               [dir]    : cache dir($cache)\n";
	print " -v                        : verbose\n";
	exit (1);
}


my ($configFile, $outFile) = @ARGV;

if ($base ne '')
{
	Carp::croak "dir $base does not exist\n" unless -d $base;
}

Carp::croak "$cache already exists\n" if -d $cache || -f $cache;
my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir at $cache\n" unless $ret == 0;


print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readConfigFile ($configFile, $base);

print "done.\n" if $verbose;

Carp::croak "must have two groups\n" unless (keys %$groups) == 2;


print "loading mapping file of id to gene to symbol...\n" if $verbose;
my %id2gene2symbolHash;
if (-f $id2gene2symbolFile)
{
	my $fin;
	open ($fin, "<$id2gene2symbolFile") || Carp::croak "cannot open file $id2gene2symbolFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my ($id, $geneId, $symbol) = split (/\t/, $line);
		$id2gene2symbolHash{$id} = "$geneId//$symbol";
	}	

	close ($fin);
}
my $n = keys %id2gene2symbolHash;

print "$n mapping entries loaded\n" if $verbose;


print "loading data of individual samples ...\n" if $verbose;

my %sampleData;
my $ASInfo;
my $nAS = 0;
my $iter = 0;

my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;

foreach my $gName (@groupNames)
{
	my $samples = $groups->{$gName}->{"samples"};
	foreach my $s (@$samples)
	{
		print "$iter: group=$gName, sample=$s\n" if $verbose;
		my $inputFile = "$s/$type.count.txt";
		$inputFile = "$base/$inputFile" if $base ne '';
		my $sdata = readASDataFile ($inputFile, $type);
		$ASInfo = $sdata->{"ASInfo"};
		if ($nAS != 0)
		{
			Carp::croak "data inconsistency detected\n" if @$ASInfo != $nAS;
		}
		else
		{
			$nAS = @$ASInfo;
		}
		$sampleData{$s} = $sdata->{"data"};
		$iter++;
	}
}

print "$iter samples, $nAS events loaded.\n" if $verbose;


print "aggregating samples in the same group ...\n" if $verbose;


my @groupData;

for (my $g = 0; $g < @groupNames; $g++)
{
	my $gName = $groupNames[$g];
	my $samples = $groups->{$gName}->{"samples"};

	foreach my $s (@$samples)
	{
		print "sample=$s\n" if $verbose;
		my $data = $sampleData{$s};
		for (my $i = 0; $i < $nAS; $i++)
		{
			my $d = $data->[$i];
			
			for (my $j = 0; $j < @$d; $j++)
			{
				$groupData[$g][$i][$j] += $d->[$j];
			}
		}
	}
}

print "performing statistical analysis ...\n" if $verbose;

my $pvalues = fisherTest (\@groupData, $cache, $verbose);

my @statResult;

my $npass = 0; 
#number of events that passed the filtering by coverage

for (my $i = 0; $i < $nAS; $i++)
{
	my ($phi1, $phi2, $n1, $n2, $inc1, $ex1, $inc2, $ex2);
	if ($type eq 'cass')
	{
		#0					1			 2                  3                      4                         5
		#isoform1Tags    isoform2Tags    exonTags        inclusionJunction1Tags inclusionJunction2Tags  skippingJunctionTags

		$inc1 = $groupData[0][$i][3] + $groupData[0][$i][4];
		$ex1 = $groupData[0][$i][5];

		$inc2 = $groupData[1][$i][3] + $groupData[1][$i][4];
		$ex2 = $groupData[1][$i][5];	

		$n1 = $inc1 + 2 * $ex1;
		$n2 = $inc2 + 2 * $ex2;
	
	
		$phi1 = $n1 > 0 ? $inc1/$n1 : 'NA';
		$phi2 = $n2 > 0 ? $inc2/$n2 : 'NA';
	}
	elsif ($type eq 'alt3' || $type eq 'alt5')
	{
		#0               1                 2           3                     4              
		#isoform1Tags	isoform2Tags	exonTags	proximalJunctionTags	distalJunctionTags
		$inc1 = $groupData[0][$i][3];
		$ex1 = $groupData[0][$i][4];

		$inc2 = $groupData[1][$i][3];
		$ex2 = $groupData[1][$i][4];

		$n1 = $inc1 + $ex1;
		$n2 = $inc2 + $ex2;

		$phi1 = $n1 > 0 ? $inc1/$n1 : 'NA';
		$phi2 = $n2 > 0 ? $inc2/$n2 : 'NA';
	}
	elsif ($type eq 'iret')
	{
		#0              1                2                    3
		#isoform1Tags	isoform2Tags	retainedIntronTags	junctionTags
		
		$inc1 = $groupData[0][$i][0];
		$ex1 = $groupData[0][$i][1];
	
		$inc2 = $groupData[1][$i][0];
		$ex2 = $groupData[1][$i][1];

		$n1 = $inc1 + 2 * $ex1;
		$n2 = $inc2 + 2 * $ex2;

		$phi1 = $n1 > 0 ? $inc1/$n1 : 'NA';
		$phi2 = $n2 > 0 ? $inc2/$n2 : 'NA';
	}
	elsif ($type eq 'mutx')
	{
		#0				1				2			3					4					5			6					7
		#isoform1Tags	isoform2Tags	5'ExonTags	5'ExonJunction1Tags	5'ExonJunction2Tags	3'ExonTags	3'ExonJunction1Tags	3'ExonJunction2Tags
		$inc1 = $groupData[0][$i][3] + $groupData[0][$i][4];
		$ex1 = $groupData[0][$i][6] + $groupData[0][$i][7];
		
		$inc2 = $groupData[1][$i][3] + $groupData[1][$i][4];
		$ex2 = $groupData[1][$i][6] + $groupData[1][$i][7];
		
		$n1 = $inc1 + $ex1;
		$n2 = $inc2 + $ex2;
		
		$phi1 = $n1 > 0 ? $inc1/$n1 : 'NA';
		$phi2 = $n2 > 0 ? $inc2/$n2 : 'NA';
	}
	elsif ($type eq 'taca')
	{
		#0              1				2			3						4
		#isoform1Tags	isoform2Tags	exonTags	inclusionJunctionTags	skippingJunctionTags

		$inc1 = $groupData[0][$i][3];
		$ex1 = $groupData[0][$i][4];
		
		$inc2 = $groupData[1][$i][3];
		$ex2 = $groupData[1][$i][4];
		
		my $asId = $ASInfo->[$i][3];
		my @cols = split ("-", $asId);

		my $nAltExon = $cols[2];

		$n1 = $inc1 + $ex1 * ($nAltExon+1);
		$n2 = $inc2 + $ex2 * ($nAltExon+1);
		
		$phi1 = $n1 > 0 ? $inc1/$n1 : 'NA';
		$phi2 = $n2 > 0 ? $inc2/$n2 : 'NA';
	}
	else
	{
		Carp::croak "incorrect AS type: $type\n";
	}	

	my $dI = $n1>0 && $n2 > 0 ? $phi1 - $phi2 : 'NA';
	my $cov = min ($inc1+$ex1, $inc2+$ex2, $inc1+$inc2, $ex1+$ex2);

	push @statResult, [$cov, $phi1, $phi2, $dI, $pvalues->[$i]];
	$npass++ if $cov >= $minCoverage;
}

print "$npass AS events has coverage >=$minCoverage\n" if $verbose;

my @statResultSort = sort {($b->[0]>=$minCoverage? 1:0) <=> ($a->[0] >= $minCoverage? 1:0) || $a->[4] <=> $b->[4]} @statResult;
#Carp::croak Dumper (\@statResultSort);

my $i;
for ($i = 0; $i < $nAS; $i++)
{
	my $q = $statResultSort[$i][4] * $npass / ($i+1);
	$statResultSort[$i][5] = $q;
	last if $q >= 1;
}

for (; $i < $nAS; $i++)
{
	$statResultSort[$i][5] = 1;
}

print "writing output to $outFile ...\n" if $verbose;
my $fout;
open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";


my @header;

push @header, "gene" if -f $id2gene2symbolFile;
push @header, qw(chrom chromStart chromEnd name score strand type isoformIDs);
push @header, "altSSDistance" if $type eq 'alt5' || $type eq 'alt3';

push @header, "coverage";
push @header, "I_g1($groupNames[0])";
push @header, "I_g2($groupNames[1])";
push @header, qw(dI_g1_vs_g2 pvalue FDR);

if ($type eq 'cass')
{
	push @header, qw(isoform1Tags_g1 isoform2Tags_g1 exonTags_g1 inclusionJunction1Tags_g1 inclusionJunction2Tags_g1 skippingJunctionTags_g1);
	push @header, qw(isoform1Tags_g2 isoform2Tags_g2 exonTags_g2 inclusionJunction1Tags_g2 inclusionJunction2Tags_g2 skippingJunctionTags_g2);
}
elsif ($type eq 'alt5' || $type eq 'alt3')
{
	push @header, qw(isoform1Tags_g1 isoform2Tags_g1 exonTags_g1 proximalJunctionTags_g1 distalJunctionTags_g1);
	push @header, qw(isoform1Tags_g2 isoform2Tags_g2 exonTags_g2 proximalJunctionTags_g2 distalJunctionTags_g2);
}
elsif ($type eq 'iret')
{
	push @header, qw(isoform1Tags_g1 isoform2Tags_g1 retainedIntronTags_g1 junctionTags_g1);
	push @header, qw(isoform1Tags_g2 isoform2Tags_g2 retainedIntronTags_g2 junctionTags_g2);
}
elsif ($type eq 'mutx')
{
	push @header, qw(isoform1Tags_g1 isoform2Tags_g1 5'ExonTags_g1 5'ExonJunction1Tags_g1 5'ExonJunction2Tags_g1 3'ExonTags_g1 3'ExonJunction1Tags_g1 3'ExonJunction2Tags_g1);
	push @header, qw(isoform1Tags_g2 isoform2Tags_g2 5'ExonTags_g2 5'ExonJunction1Tags_g2 5'ExonJunction2Tags_g2 3'ExonTags_g2 3'ExonJunction1Tags_g2 3'ExonJunction2Tags_g2);
}
elsif ($type eq 'taca')
{
	push @header, qw(isoform1Tags_g1 isoform2Tags_g1 exonTags_g1 inclusionJuncctionTags_g1 skippingJunctionTags_g1);
	push @header, qw(isoform1Tags_g2 isoform2Tags_g2 exonTags_g2 inclusionJuncctionTags_g2 skippingJunctionTags_g2);
}
else
{
	Carp::croak "incorrect AS type: $type\n";
}

print $fout join ("\t", @header), "\n";

for (my $i = 0; $i < $nAS; $i++)
{
	my @out;

	my $gene2symbol = exists $id2gene2symbolHash{$ASInfo->[$i]->[3]} ? $id2gene2symbolHash{$ASInfo->[$i]->[3]} : "NA//NA";	
	push @out, $gene2symbol if -f $id2gene2symbolFile;

	push @out, @{$ASInfo->[$i]};
	push @out, @{$statResult[$i]};
	push @out, @{$groupData[0][$i]};
	push @out, @{$groupData[1][$i]};

	#                0    1      2      3     4             5
	#@statResult, [$cov, $phi1, $phi2, $dI, $pvalues->[$i], FDR];

	if ($filterOutput == 0 || ($statResult[$i][0] >= $minCoverage && ABS($statResult[$i][3]) >= $deltaI && $statResult[$i][5] <= $FDR))
	{
		print $fout join ("\t", @out), "\n";
	}
}

close ($fout);

system ("rm -rf $cache");


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
		$groups{$groupName}->{"id"} = $i++ unless exists $groups{$groupName};
		push @{$groups{$groupName}->{"samples"}}, $sampleName;

		my $inputFile = "$sampleName/$type.count.txt";
		$inputFile = "$base/$inputFile" if $base ne '';

		Carp::croak "Input file $inputFile does not exist\n" unless -f $inputFile;
	}
	close ($fin);
	return \%groups;
}

sub readASDataFile
{
	my ($inputFile, $type) = @_;
	
	my $fin;
	my @data;
	my @ASInfo;
	open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		
		my @cols = split (/\t/, $line);
		my (@infoCols, @dataCols);
		
		if ($type eq 'cass' || $type eq 'iret' || $type eq 'mutx' || $type eq 'taca')
		{
			@infoCols = @cols[0..7];
			@dataCols = @cols[8..$#cols];
		}
		elsif ($type eq 'alt3' || $type eq 'alt5')
		{
			@infoCols = @cols[0..7];
			push @infoCols, $cols[10];
			@dataCols = @cols[8..9];
			push @dataCols, @cols[11..$#cols];
		}
		else
		{
			Carp::croak "incorrect AS type: $type\n";
		}
		
		push @ASInfo, \@infoCols;
		push @data, \@dataCols;
	}
	close ($fin);

	return {ASInfo=>\@ASInfo, data=>\@data};
}


sub fisherTest
{
	my ($groupData, $cache, $verbose)  = @_;

	print "generating data and script for statistical test ...\n" if $verbose;

	my $fout;

	my $testDataFile = "$cache/data.txt";
	open ($fout, ">$testDataFile") || Carp::croak "cannot open file $testDataFile to write\n";
	for (my $i = 0; $i < $nAS; $i++)
	{
		my $in1 = $groupData->[0][$i][0];
		my $ex1 = $groupData->[0][$i][1];
		my $in2 = $groupData->[1][$i][0];
		my $ex2 = $groupData->[1][$i][1];

		print $fout join ("\t", $in1, $ex1, $in2, $ex2), "\n";
	}
	close ($fout);


	#write R script
	my $testOutFile = "$cache/p.txt";
	
	my $testScriptFile = "$cache/test.R";
	open ($fout, ">$testScriptFile") || Carp::croak "cannot open file $testScriptFile to write\n";

	print $fout <<EOF;

data <- read.table ("$testDataFile", sep="\\t", header=F);
n <- dim(data)[1];

p <- 1:n;
for (i in 1:n)
{
    if (i-as.integer(i/1000) * 1000 == 0)
    {
        cat (sprintf("%d\\n",i));
    }

    dat <- matrix(data=as.numeric(data[i,]), nrow=2)
    out <- fisher.test(dat);
    p[i] <- out\$p.value;
}

write.table (p, "$testOutFile", sep="\\t", quote=F, row.names=F, col.names=F);

EOF
	
	close ($fout);

	my $cmd = "R --no-save < $testScriptFile";
	$cmd = "$cmd > /dev/null" unless $verbose;
	
	print "$cmd\n" if $verbose;

	my $ret = system ($cmd);
	Carp::croak "cmd=$cmd failed: $?\n" if $ret != 0;


	#read test output
	print "reading test results from $testOutFile ...\n" if $verbose;
	my $fin;
	open ($fin, "<$testOutFile") || Carp::croak "cannot open file $testOutFile to read\n";
	my @p;
	my $line;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		push @p, $line;
	}
	close ($fin);

	return \@p;
}


