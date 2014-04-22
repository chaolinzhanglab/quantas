#!/usr/bin/perl -w

use strict;


use Getopt::Long;
use Carp;
use Data::Dumper;
use File::Basename;

use MyConfig;
use Vcf;

my $prog = basename ($0);

my $big = 0;
my $bamListFile = "";
my $uniq = 0;
my $refGenome = "";
my $trimEnd = 5;
my $baseQ = 13;
my $offset = 33;
my $readLen = 101;
my $minDepth = 1;
my $minAltBase = 1;

my $knownSiteFile = "";
my $useKnownAltBase = 0;

my $mpileup = 0;

my $cache = MyConfig::getDefaultCache ($prog);
my $keepCache = 0;

my $verbose = 0;

GetOptions (
	"b:s"=>\$bamListFile,
	"l=i"=>\$readLen,
	"r=s"=>\$refGenome,
	"s:s"=>\$knownSiteFile,
	"big"=>\$big,
	"x"=>\$useKnownAltBase,
	"u"=>\$uniq,
	"t:i"=>\$trimEnd,
	"Q:i"=>\$baseQ,
	"d:i"=>\$minDepth,
	"a:i"=>\$minAltBase,
	"m"=>\$mpileup,
	"c:s"=>\$cache,
	"keep-cache"=>\$keepCache,
	"v"=>\$verbose);

if (@ARGV < 1 && $bamListFile eq '')
{
	print "extract SNVs from RNASeq bam files\n";
	print "Usage: $prog [options] <in.bam>\n";
	print " -l [int]    : read length\n";
	print " -r          : reference genome indexed by samtools\n";
	print "[options]\n";
	#print " -b [string] : a file with a list of bam files\n";
	print " -s [string] : a vcf file with a list of known sites (will disable -d -a)\n";
	print " -big        : the known site file is big\n";
	print " -x          : count the specified altbase (effective only with -s)\n";
	print " -u          : consider only uniquely mapped reads\n";
	print " -t [int]    : the [int] nucleotides at the end will be ignored ($trimEnd)\n";
	print " -Q [int]    : skip bases with baseQ/BAQ smaller than INT ($baseQ)\n";
	print " -d [int]    : minimum depth to report a potential SNV ($minDepth)\n";
	print " -a [int]    : minimum alt base to report a potential SNV ($minAltBase)\n";
	print " -c          : cache dir ($cache)\n";
	print " --keep-cache: keep cache\n";
	#print " -m         : the input is the mpileup file (for debug only)\n";
	print " -v          : verbose\n";
	exit (1);
}

my @bamFiles = @ARGV;

Carp::croak "trimmed too much!!\n" if $readLen - $trimEnd * 2 <= 0;


srand (1);

my $ret = system ("mkdir $cache");
Carp::croak "cannot make dir $cache:$?\n" unless $ret == 0;


if (-f $bamListFile)
{
	my $fin;
	open ($fin, "<$bamListFile") || Carp::croak "cannot open $bamListFile ...\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\#/;
		next if $line =~/^\s*$/;
		push @bamFiles, $line;
	}
	close ($fin);
}

my $n = @bamFiles;
print STDERR "$n bam files to be processed\n" if $verbose;

my $uniqFlag = $uniq ? "-q 1" : "";


my %knownSiteHash;
my $knownSitePosFile = "";
$n = 0;

if ($knownSiteFile ne '')
{
	#pos for known sites that are fed into samtools mpileup
	$knownSitePosFile = "$cache/site.pos.txt";
	my $cmd = "cut -f 1-2 $knownSiteFile > $knownSitePosFile";
	print STDERR "$cmd\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" unless $ret == 0;

	#process known sites to be loaded in hash tables all at once or in batches
	if ($big)
	{
		print STDERR "split known sites from $knownSiteFile ...\n" if $verbose;
		my $ret = splitVcfFileByChrom ($knownSiteFile, $cache, $verbose);
		%knownSiteHash = %$ret;
	}
	else
	{
		print STDERR "loading known sites from $knownSiteFile ...\n" if $verbose;

		my $sites = readVcfFile ($knownSiteFile, $verbose);
		foreach my $snv (@$sites)
		{
			my $chrom = $snv->{'chrom'};
			push @{$knownSiteHash{$chrom}}, $snv;
		}
	}

	print STDERR "get site count broken down into chromosomes ...\n" if $verbose;
	foreach my $chrom (sort keys %knownSiteHash)
	{

    	my $n = $knownSiteHash{$chrom};
    	$n = ref($n) eq 'HASH' ? $n = $n->{'n'} : @$n;
    	print STDERR "$chrom : $n sites\n" if $verbose;
	}
}


#-q 1: consider only uniquely mapped reads
#-d  : max coverage
#-D  : output per sample DP
#-g  : output vcf
#-I  : skip indels

my $fin;


if ($mpileup)
{
	open ($fin, "<".$bamFiles[0]) || Carp::croak "cannot open file $bamFiles[0] to read\n";
}
else
{

	my $cmd = "samtools mpileup -O $uniqFlag -d10000000 -f $refGenome " . join (" ", @bamFiles);
	$cmd = "samtools mpileup -l $knownSitePosFile -O $uniqFlag -d10000000 -f $refGenome " . join (" ", @bamFiles) if (-f $knownSitePosFile);

	print STDERR $cmd, "\n" if $verbose;
	open ($fin, "$cmd | ") || Carp::croak "cannot pipe samtools mpileup\n";
}

print join ("\t", "#CHROM", "POS ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), "\n";

my $iter = 0;
my $prevChrom = "";

my $knownSiteOnChrom = {};

my $fout;

my $grandTotalNonRefBase = 0;
my $grandTotalAltBase = 0;
my $grandTotalRead = 0;

while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	
	print STDERR "$iter ...\n" if $verbose && $iter % 1000000 == 0;
	$iter++;

	my ($chrom, $pos, $refBase, $count, $readBaseStr, $qualStr, $readPosStr) = split (/\t/, $line);
	$grandTotalRead++ if $readBaseStr =~/\^/;
	

	if ((-f $knownSiteFile) && $chrom ne $prevChrom)
	{
		#file handles for tmp site output
		close ($fout) if $prevChrom ne '';
	
		$prevChrom = $chrom;
		$knownSiteOnChrom = {};
		
		my $f = "$cache/$chrom.new.vcf";
		open ($fout, ">$f") || Carp::croak "cannot open $f to write\n"; 

		#load data for the current chrom
		my $sites = [];
		if ($big)
		{
			my $f = "$cache/$chrom.vcf";
			if (-f $f)
			{
				print STDERR "load known sites on $chrom...\n" if $verbose;
				$sites = readVcfFile ($f, $verbose);
				my $n = @$sites;
				print STDERR "$n sites loaded\n" if $verbose;

			}
		}
		else
		{
			if (exists $knownSiteHash{$chrom})
			{
				$sites = $knownSiteHash{$chrom};
			}
		}
		
		#sort sites into hash tables, and discard sites with ambiguous refBase
		foreach my $snv (@$sites)
		{
			my $refBase = $snv->{'refBase'};

			if ($refBase !~/[ACGTacgt]/)
			{
				print STDERR "non [ACGT] in refBase, will skip:", Dumper ($snv), "\n";
				next;
			}
		
			my $position = $snv->{'position'}; #zero-based
			$knownSiteOnChrom->{$position} = $snv;
		}
	}

	if (-f $knownSiteFile)
	{
		#limit to known sites
		next unless exists $knownSiteHash{$chrom} && exists $knownSiteOnChrom->{$pos-1}; #$pos is 1-based
	}

	$refBase = uc($refBase);

	next unless $refBase =~/[ACGT]/;

	#$readBaseStr=~s/\^\S|\$//g; #remove start mark and mapping score and end mark of a read
	#$readBaseStr=~s/[\^\S|\$]//g; #remove start mark and mapping score and end mark of a read
	
	#Carp::croak "readBaseStr=$readBaseStr\n";
	
	if (not (-f $knownSiteFile))
    {
		next unless $readBaseStr =~/[ACGTacgt]/;
		#no variants at all, do not waste time
		#note that we will not skip if a list of known sites is provided, in which case, we will try to extract info for all sites
	}

	my $readBase = splitReadBaseStr($readBaseStr);
	
	my @qual = split (//, $qualStr);
	my @readPos = split (/\,/, $readPosStr);

	Carp::croak "incorrect parsing of readBase: from=$readBaseStr, to=", join("|", @$readBase), "\n"
	if @$readBase != $count;

	my %readBaseHash = (
		'+'=>{A=>0, C=>0, G=>0, T=>0},
		'-'=>{A=>0, C=>0, G=>0, T=>0});

	for (my $i = 0; $i < $count; $i++)
	{
		next unless ord ($qual[$i]) - $offset >= $baseQ  && $readPos[$i] > $trimEnd && $readPos[$i] <= $readLen - $trimEnd;
	
		my $b = $readBase->[$i];
		if ($b eq '.')
		{
			$readBaseHash{'+'}->{$refBase}++;
		}
		elsif ($b eq ',')
		{
			$readBaseHash{'-'}->{$refBase}++;
		}
		elsif ($b eq 'A' || $b eq 'C' || $b eq 'G' || $b eq 'T')
		{
			$readBaseHash{'+'}->{$b}++;
			$grandTotalNonRefBase++;
		}
		elsif ($b eq 'a' || $b eq 'c' || $b eq 'g' || $b eq 't')
		{
			$readBaseHash{'-'}->{uc($b)}++;
			$grandTotalNonRefBase++;
		}
	}
	
	$readBaseHash{'A'} = $readBaseHash{'+'}->{'A'} + $readBaseHash{'-'}->{'A'};
	$readBaseHash{'C'} = $readBaseHash{'+'}->{'C'} + $readBaseHash{'-'}->{'C'};
	$readBaseHash{'G'} = $readBaseHash{'+'}->{'G'} + $readBaseHash{'-'}->{'G'};
	$readBaseHash{'T'} = $readBaseHash{'+'}->{'T'} + $readBaseHash{'-'}->{'T'};

	my $altBase = "";
	if ((-f $knownSiteFile) && $useKnownAltBase && exists $knownSiteOnChrom->{$pos-1})
	{
		$altBase = $knownSiteOnChrom->{$pos-1}->{'altBase'};
	}
	else
	{
		$altBase = assignAltBase (\%readBaseHash, $refBase);
	}

	my $altBaseSum = $readBaseHash{'+'}->{$altBase} + $readBaseHash{'-'}->{$altBase};
	$grandTotalAltBase += $altBaseSum;

	my $depth = $readBaseHash{'+'}->{$refBase} + $readBaseHash{'-'}->{$refBase} + $readBaseHash{'+'}->{$altBase} + $readBaseHash{'-'}->{$altBase};
	
	if (not (-f $knownSiteFile))
	{
		next unless $readBaseHash{$altBase} > 0;
		next unless $depth >= $minDepth && $altBaseSum >= $minAltBase;
	}

	my $DP4=join (",", $readBaseHash{'+'}->{$refBase}, $readBaseHash{'-'}->{$refBase}, $readBaseHash{'+'}->{$altBase}, $readBaseHash{'-'}->{$altBase});
	if (-f $knownSiteFile)
	{
		#this is the temp output file
		#sites with zero coverage will not be here, and they will be recovered in a post processing step
		print $fout join ("\t", $chrom, $pos, '.', $refBase, $altBase, 0, '.', "DP4=$DP4"), "\n";
	}
	else
	{
		#print site
		print join ("\t", $chrom, $pos, '.', $refBase, $altBase, 0, '.', "DP4=$DP4"), "\n";
	}
}

my $nonRefError = $grandTotalRead > 0 ? $grandTotalNonRefBase / $grandTotalRead / ($readLen - $trimEnd * 2) : "NA";
my $altError = $grandTotalRead > 0 ? $grandTotalAltBase / $grandTotalRead / ($readLen - $trimEnd * 2) : "NA";

#if a list of known sites are provided, the error rate is not accurate, so we do not report
print STDERR "total mapped reads = $grandTotalRead, total non-ref base = $grandTotalNonRefBase, non-ref per base=$nonRefError, total alt base = $grandTotalAltBase, alt per base=$altError\n" unless -f $knownSiteFile; 

close ($fin);
close ($fout) if $prevChrom ne '';

if (-f $knownSiteFile)
{
	#output known sites
	print STDERR "output updated sites ...\n" if $verbose;

	foreach my $chrom (sort keys %knownSiteHash)
	{
		print STDERR "processing $chrom ...\n" if $verbose;

		my %newSiteHash;
		my $f = "$cache/$chrom.new.vcf";
		if (-f $f)
		{
			my $sites = readVcfFile ($f, $verbose);
			foreach my $snv (@$sites)
			{
				my $position = $snv->{'position'};
				#there will be never ambiguous refBase at this stage
				$newSiteHash{$position} = $snv;
			}
		}
		
		my $sites;
		if ($big)
		{
			my $f = "$cache/$chrom.vcf";
			$sites = readVcfFile ($f, $verbose);
		}
		else
		{
			$sites = $knownSiteHash{$chrom};
		}
		
		#output
		foreach my $snv (@$sites)
		{
			my $position = $snv->{'position'}; #zero-based
			my $refBase = $snv->{'refBase'};

			if ($refBase !~/[ACGTacgt]/)
			{
				print STDERR "non [ACGT] in refBase, will skip:", Dumper ($snv), "\n";
				next;
			}
			
			if (exists $newSiteHash{$position})
			{
				my $newSnv = $newSiteHash{$position};
				$snv->{'altBase'} = $newSnv->{'altBase'};
				$snv->{'info'}->{'DP4'} = $newSnv->{'info'}->{'DP4'};
			}
			else
			{
				$snv->{'info'}->{'DP4'} = join (",", 0,0,0,0);
			}
			my $DP4 = $snv->{'info'}->{'DP4'};
			print join ("\t", $snv->{'chrom'}, $snv->{'position'} + 1, $snv->{'id'}, $snv->{'refBase'}, $snv->{'altBase'}, $snv->{'qual'}, $snv->{'filter'}, "DP4=$DP4"), "\n";
		}
	}
}

system ("rm -rf $cache") unless $keepCache;


sub splitReadBaseStr
{
	my $str=$_[0];
	my @ret;
					# ins                       #del start                del
	#while ($str =~/([\.\,]\+[0-9]+[ACGTNacgtn]+|[\.\,]\-[0-9]+[ACGTNacgtn]+|[\*\.\,ACGTNacgtn\>\<])/g)
	#{
	#	push @ret, $1;
	#}
	
	my @chars = split (//, $str);
	for (my $i = 0; $i < @chars; $i++)
	{
		my $c = $chars[$i];
		if ($c =~/[\*\.\,ACGTNacgtn\>\<]/)
		{
			push @ret, $c;
		}
		elsif ($c eq '+' || $c eq '-')
		{
			#indel
			my $type = $c;

			my $len = "";
			#get the length of indel from the next base
			for ($i++; $i < @chars; $i++)
			{
				my $c = $chars[$i];
				last unless $c =~/[0-9]/;
				$len .= $c;
			}

			#now the first base of indel sequence
			my $startIdx = $i;
			my $endIdx = $i + $len - 1;
			
			my $indelStr = join ("", @chars[$startIdx..$endIdx]);
			#ignore indels for now
			
			#push @ret, join("", $type, $len, $indelStr);
			
			#jump to the last base of indel str
			$i = $endIdx;
		}
		elsif ($c eq '$')	
		{
			#end of mapping, skip this char
			next;
		}
		elsif ($c eq '^')
		{
			#start of mapping, skip this and next char
			$i++;
			next;
		}
	}
	return \@ret;
}

