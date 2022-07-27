#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Math::CDF qw(:all);

use Bed;
use Sequence;
use Vcf;
use Common;

use MyConfig;

my $prog = basename ($0);
my $verbose = 0;
my $annotFile = "";
my $dbkey = "mm10";
my $typeStr = "";

my $knownSNPBedFile = "";
my $excludeSNP = 0;

my $progDir = dirname ($0);
my $big = 0;
my $guess = 0;
my $printAll = 0;

my $cache = MyConfig::getDefaultCache ($prog);
my $keepCache = 0;



GetOptions (
	"annot=s"=>\$annotFile,
	"dbkey:s"=>\$dbkey,
	"snp:s"=>\$knownSNPBedFile,
	"ex-snp"=>\$excludeSNP,
	"type:s"=>\$typeStr,
	"big"=>\$big,
	"guess"=>\$guess,
	"all"=>\$printAll,
	"c:s"=>\$cache,
	"keep-cache"=>\$keepCache,
	"v"=>\$verbose);

if (@ARGV != 1)
{
	print "annotate gene information\n";
	print "Usage: $prog [options] <in.vcf>\n";
	print " in.vcf: .gz file accepted\n";
	print " -annot [string] : configuration of gene annotation files\n";
	print " -dbkey [string] : dbkey ($dbkey)\n";
	print " -snp   [string] : BED file with coordinates of known snps\n";
	print " -guess          : make the best guess of strand if ambiguous\n";
	print " -type  [string] : type of SNV (on the sense strand to include or exclude\n";
	print "                 : examples: IN:A>G;C>T or EX:A>G;C>T\n";
	print " -all            : include those without overlapping genes\n";
	print " --ex-snp        : exclude known snp in output (together with -snp)\n";
	print " -c     [string] : cache($cache)\n";
	print " --keep-cache    : keep cache\n";  
	print " -v              : verbose\n";
	exit (1);
}

my ($inVcfFile) = @ARGV;


#allowed SNV type (on the sense strand)

my %allowedTypeHash;

foreach my $b1 (qw(A C G T))
{
	foreach my $b2 (qw(A C G T))
	{
		next if $b1 eq $b2;	
		$allowedTypeHash{"$b1>$b2"} = 1;
	}
}

if ($typeStr ne '')
{
	Carp::croak "incorrect format:$typeStr\n" unless $typeStr =~/^(.*?):(.*?)$/;
	my $action = $1;
	my $str = $2;
	if ($action eq 'IN')
	{
		map {$allowedTypeHash{$_} = 0} keys %allowedTypeHash;
	}
	elsif ($action ne 'EX')
	{
		Carp::croak "incorrect format:$typeStr\n" unless $typeStr =~/^(.*?):(.*?)$/;
	}

	my @types = split (";", $str);
	foreach my $t (@types)
	{
		$t = uc($t);
		Carp::croak "unrecognized type: $t in $typeStr\n" unless exists $allowedTypeHash{$t};
		$allowedTypeHash{$t} = $action eq 'IN' ? 1 : 0;
	}
}



my $ret = system ("mkdir $cache");
Carp::croak "cannot create $cache\n" unless $ret == 0;

my $snvStrandBedFile = "$cache/snv.strand.bed";

#pos strand
#my $cmd = "grep -P -v \"^\#\" $inVcfFile | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1\".+\\t0\\t+\"}' > $snvStrandBedFile";
my $cmd = "grep -v \"^\#\" $inVcfFile | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1\".+\\t0\\t+\"}' > $snvStrandBedFile";
if ($inVcfFile =~/\.gz$/)
{
	$cmd = "zcat $inVcfFile | grep -v \"^\#\" | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1\".+\\t0\\t+\"}' > $snvStrandBedFile";
}

print STDERR "$cmd\n" if $verbose;

$ret = system ($cmd);
Carp::croak "CMD=$cmd failed: $?\n" unless $ret == 0;



#neg strand
#$cmd = "grep -P -v \"^\#\" $inVcfFile | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1\".-\\t0\\t-\"}' >> $snvStrandBedFile";
$cmd = "grep -v \"^\#\" $inVcfFile | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1\".-\\t0\\t-\"}' >> $snvStrandBedFile";
if ($inVcfFile =~/\.gz$/)
{
	$cmd = "zcat $inVcfFile | grep -v \"^\#\" | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1\".-\\t0\\t-\"}' >> $snvStrandBedFile";
}

print STDERR "$cmd\n" if $verbose;
$ret = system ($cmd);
Carp::croak "CMD=$cmd failed: $?\n" unless $ret == 0;

#add gene annotations
print STDERR "reading gene annotations ...\n" if $verbose;

my $verboseFlag = $verbose ? "-v" : "";
my $bigFlag = $big ? "-big" : "";

my $snvAnnotFile = "$cache/snv.annot.txt";
$cmd = "$progDir/bed2annotation.pl $verboseFlag $bigFlag -conf $annotFile -dbkey $dbkey -ss -gene -region -c $cache/annot -v $snvStrandBedFile $snvAnnotFile 1>&2";
print STDERR "$cmd\n" if $verbose;

$ret = system ($cmd);

Carp::croak "CMD=$cmd failed: $?\n" unless $ret == 0;


print STDERR "load gene annotations ...\n" if $verbose;
my $fin;

open ($fin, "<$snvAnnotFile") || Carp::croak "cannot open file $snvAnnotFile to read\n";

my %snvAnnotInfo;
my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\#/;
	print STDERR "$iter ...\n" if $verbose && $iter % 10000 == 0;
	$iter++;
	my ($name, $geneId, $geneSymbol, $region) = split ("\t", $line);
	next unless $geneId ne '';
	my ($chromPos, $strand) = split(/\./, $name);
	$snvAnnotInfo{$chromPos}{'gene'}{$strand} = {geneId=>$geneId, symbol=>$geneSymbol, region=>$region, strand=>$strand};
}
close ($fin);



if (-f $knownSNPBedFile)
{

	#
	print STDERR "annotating known snps ...\n" if $verbose;
	my $snvBedFile = "$cache/snv.nostrand.bed";
	#$cmd = "grep -P -v \"^\#\" $inVcfFile | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1}' > $snvBedFile";
	$cmd = "grep -v \"^\#\" $inVcfFile | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1}' > $snvBedFile";
	if ($inVcfFile =~/\.gz$/)
	{
		$cmd = "zcat $inVcfFile | grep -v \"^\#\" | awk '{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$1\":\"\$2-1}' > $snvBedFile";
	}
	print STDERR "$cmd\n" if $verbose;

	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" unless $ret == 0;

	my $snv_vs_snpFile = "$cache/snv_vs_knownsnp.bed";
	$cmd = "perl $progDir/tagoverlap.pl $verboseFlag -c $cache/snp -big -region $knownSNPBedFile $snvBedFile $snv_vs_snpFile 1>&2";
	print STDERR "$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" unless $ret == 0;

	#
	print STDERR "loading snp comparison results...\n" if $verbose;
	my $iter = 0;
	open ($fin, "<$snv_vs_snpFile") || Carp::croak "cannot open file $snv_vs_snpFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		print STDERR "$iter ...\n" if $verbose && $iter % 1000 == 0;
		$iter++;

		my $site = lineToBed ($line);
		my $name = $site->{'name'};
		my ($chromPos, $snpName) = split ("//", $name);
		next unless exists $snvAnnotInfo{$chromPos};

		$snvAnnotInfo{$chromPos}->{'snp'} = $snpName;
	}

	close ($fin);
}

if ($inVcfFile =~/\.gz$/)
{
	open ($fin, "gunzip -c $inVcfFile | ")||Carp::croak "cannot open file $inVcfFile to read\n";
}
elsif ($inVcfFile =~/\.bz2$/)
{
    open ($fin, "bunzip2 -c $inVcfFile | ")||Carp::croak "cannot open file $inVcfFile to read\n";
}
else
{
	open ($fin, "<$inVcfFile") || Carp::croak "cannot open file $inVcfFile to read\n";
}


$iter = 0;

print generateVcfHeader (), "\n";

while (my $line =<$fin>)
{
	chomp $line;
	next if $line =~/^\#/;

	print STDERR "$iter ...\n" if $verbose && $iter % 10000 == 0;
	$iter++;
	my $snv = lineToVcf ($line);

	#Carp::croak Dumper ($snv), "\n";
	my $chrom = $snv->{'chrom'};
	my $position = $snv->{'position'}; #zero-based coordinate
	my $chromPos = "$chrom:$position";

	if (not exists $snvAnnotInfo{$chromPos})
	{
		#no gene annotation
		#skip

		print vcfToLine ($snv), "\n" if $printAll;
		next;
	}	

	#get the gene with smaller id if both strand has overlapping genes
	#
	my @geneInfo = sort {$a->{'geneId'} <=> $b->{'geneId'}} values %{$snvAnnotInfo{$chromPos}{'gene'}};
	
	if (@geneInfo > 1)
	{
		next unless $guess;	
	}

	my $gInfo = $geneInfo[0];
	#Carp::croak Dumper ($gInfo), "\n";

	my $snvInfo = $snv->{'info'};
	
	if (exists $snvAnnotInfo{$chromPos}->{'snp'})
	{
		next if $excludeSNP;
		$snvInfo->{'snp'} = $snvAnnotInfo{$chromPos}->{'snp'};
	}

	$snvInfo->{'strand'} = @geneInfo > 1 ? "." . $gInfo->{'strand'} : $gInfo->{'strand'};	
	$snvInfo->{'gene'} = join ("//", $gInfo->{'geneId'}, $gInfo->{'symbol'});
	$snvInfo->{'region'} = $gInfo->{'region'};

	my $refBase = $snv->{'refBase'};
	my $altBase = $snv->{'altBase'};

	if ($snvInfo->{'strand'} eq '-' || $snvInfo->{'strand'} eq '.-')
	{
		$refBase = complement ($refBase);
		$altBase = complement ($altBase);
	}
	
	my $t = "$refBase>$altBase";
	Carp::croak "unrecognized type $t in SNV:", Dumper ($snv), "\n" unless exists $allowedTypeHash{$t};
	
	print vcfToLine ($snv), "\n" if $allowedTypeHash {$t};
}
close ($fin);

system ("rm -rf $cache") unless $keepCache;

