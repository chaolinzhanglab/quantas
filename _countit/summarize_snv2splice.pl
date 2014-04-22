#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Carp;
use File::Basename;
use Getopt::Long;

use MyConfig;
use Vcf;
use Sam;
use Sequence;


my $prog = basename ($0);
my $progDir = dirname ($0);

my $refGenome = "";
my $trimEnd = 5;
my $baseQ= 13;

my $snv2spliceSiteFile = "";

my $verbose = 0;
my $cache = getDefaultCache ($prog);
my $keepCache = 0;
my $samHeaderFile = "";

GetOptions (
	"snv2ss:s"=>\$snv2spliceSiteFile,
	"header:s"=>\$samHeaderFile,
	"r:s"=>\$refGenome,
	"t:i"=>\$trimEnd,
	"Q:i"=>\$baseQ,
	"c:s"=>\$cache,
	"keep-cache"=>\$keepCache,
	"v"=>\$verbose);

if (@ARGV != 3)
{
	print "summarize association of snv with splice sites\n";
	print "Usage: $prog [options] <snv.vcf> <in.sam> <out.txt>\n";
	print " -snv2ss [string]: map of snv to splice sites\n";
	print " -header [string]: file with sam header lines\n";
	print " -r      [string]: reference genome indexed by samtools\n";

	print " -t      [int]   : the [int] nucleotides at the end will be ignored ($trimEnd)\n";
	print " -Q      [int]   : skip bases with baseQ/BAQ smaller than INT ($baseQ)\n";
	print " -c      [string]: cache dir ($cache)\n";
	print " --keep-cache    : keep cache dir\n";
	print " -v              : verbose\n";
	exit (1);
}

my ($snvVcfFile, $inSamFile, $outFile) = @ARGV;

my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir $cache:$?\n" unless $ret == 0;


print "loading snvs ...\n" if $verbose;

my $snvs = readVcfFile ($snvVcfFile, $verbose);

my %snvHash;
foreach my $s (@$snvs)
{
	my $id = $s->{'id'};
	$snvHash{$id} = $s;
}
my $n = @$snvs;
print "$n snvs loaded\n" if $verbose;




print "loading reads ...\n" if $verbose;
my $reads = readSamFile ($inSamFile, $verbose);

my %readHash;
foreach my $r (@$reads)
{
	my $qName = $r->{'QNAME'};
	my $tags = $r->{'TAGS'};

	my $spliceSiteId = "";
	if ($tags =~/\sSS:Z:(\S*?)$/)
	{
		$spliceSiteId = $1;
	}
	else
	{
		Carp::croak "cannot find splice site id in:", Dumper ($r), "\n";
	}
	
	push @{$readHash{$spliceSiteId}}, $r;
}

$n = @$reads;
print "$n spliced reads loaded\n" if $verbose;


print "counting alleles ...\n" if $verbose;
my $fin;

my %snv2ssHash;
open ($fin, "<$snv2spliceSiteFile") || Carp::croak "cannot open file $snv2spliceSiteFile to read\n";

my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

my $iter = 0;

print $fout join ("\t", "#chrom", "chromStart", "chromEnd", "snvID", "score", "strand", "type", "isoformIDs", "altBaseCount", "refBaseCount", "refBase", "altBase", "spliceSiteID", "refSense", "refAntisense", "altSense", "altAntisense"), "\n";

while (my $line =<$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	print "$iter ...\n" if $verbose;
	$iter++;

	my ($snvId, $ssId, $ssType) = split ("\t", $line);

	my ($refSense, $refAntisense, $altSense, $altAntisense) = (0, 0, 0, 0);	

	next unless exists $snvHash{$snvId};

	my $snv = $snvHash{$snvId};
	Carp::croak "no strand information found:", Dumper ($snv), "\n" unless exists $snv->{'info'}{'strand'};
	
	if (exists $readHash{$ssId})
	{
		my $ret = summarizeSNV2SS ($snvHash{$snvId}, $readHash{$ssId}, $cache);
		
		$refSense = $ret->{'refPos'};
		$refAntisense = $ret->{'refNeg'};
		$altSense = $ret->{'altPos'};
		$altAntisense = $ret->{'altNeg'};
	}
	
	my $refBase = $snv->{'refBase'};
	my $altBase = $snv->{'altBase'};
	
	my $strand = $snv->{'info'}{'strand'};
	if ($strand eq '-' || $strand eq '.-')
	{
		$refBase = complement ($refBase);
        $altBase = complement ($altBase);
        ($refAntisense, $refSense, $altAntisense, $altSense) = ($refSense, $refAntisense, $altSense, $altAntisense);
	}
	
	my $refBaseCount = $refSense + $refAntisense;
	my $altBaseCount = $altSense + $altAntisense;
	
	print $fout join ("\t", $snv->{'chrom'}, $snv->{'position'}, $snv->{'position'} + 1, $snvId, $refBaseCount+$altBaseCount, $snv->{'info'}->{'strand'}, 
		'SNV2SPLICE', 'ALT/REF', $altBaseCount, $refBaseCount, $refBase, $altBase, $ssId, $refSense, $refAntisense, $altSense, $altAntisense), "\n";
}

close ($fin);
close ($fout);


sub summarizeSNV2SS
{
	my ($snv, $reads, $cache) = @_;

	my $tmpVcfFile = "$cache/snv.tmp.vcf";
	my $tmpSamFile = "$cache/read.tmp.sam";
	my $tmpBamPrefix = "$cache/read.tmp";
	my $tmpBamFile = "$tmpBamPrefix.bam";

	writeVcfFile ([$snv], $tmpVcfFile);
	writeSamFile ($reads, $tmpSamFile);
	
	my $cmd = "cat $samHeaderFile $tmpSamFile | samtools view -uSh - | samtools sort - $tmpBamPrefix 2> /dev/null";
	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

	my $readLen = length($reads->[0]->{'SEQ'});
	my $outVcfFile = "$cache/snv.count.vcf";

	$cmd = "perl $progDir/bam2vcf.pl -l $readLen -r $refGenome -s $tmpVcfFile -x -t $trimEnd -Q $baseQ -c $cache/vcf $tmpBamFile > $outVcfFile 2> /dev/null";
	$ret = system ($cmd);
    Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;
	
	my $site = readVcfFile ($outVcfFile);
	$site = $site->[0];
	my $dp4 = $site->{'info'}->{'DP4'};
	
	my ($refPos, $refNeg, $altPos, $altNeg) = split (/\,/, $dp4);
	return {refPos=>$refPos, refNeg=>$refNeg, altPos=>$altPos, altNeg=>$altNeg};
}



