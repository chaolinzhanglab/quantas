#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

use Bed;
use Common;
use Align;

my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $verbose = 0;

my $sam = 0;
my $uniq = 0;
my $nsplit = 0;
my $splitSize = 4000000; #4 million

#my $separateStrand = 0;
my $libraryType = "unstranded"; 
# or 'stranded-sa': read1 is sense and read2 is antisense
# 'stranded-as' : read1 is antisense and read2 is sense

my $isoformBedFile = "";
my $cleanIsoform = 0;
my $bigExonSize = 400;
my $outputDir = "./gapless_out";

my $big = 0;
my $printSingleton = 0;
my $usePreCalScore = 0;

my $sizeDistOnly = 0;
my $sizeDistFile = "";
my $keepCache = 0;

GetOptions (
	"sam"=>\$sam,
	"uniq"=>\$uniq,
	"isoform=s"=>\$isoformBedFile,
	"split-size:i"=>\$splitSize,
	"split:i"=>\$nsplit,
	"clean-isoform"=>\$cleanIsoform,
	"E|big-exon-size:i"=>\$bigExonSize,
	"big"=>\$big,
	"use-pre-calc-score"=>\$usePreCalScore,
	"use-pre-calc-size-dist:s"=>\$sizeDistFile,
	#"ss"=>\$separateStrand,
	"library-type:s"=>\$libraryType,
	"size-dist-only"=>\$sizeDistOnly,
	"print-singleton"=>\$printSingleton,
	"o|outp-dir:s"=>\$outputDir,
	"keep-cache"=>\$keepCache,
	"v|verbose"=>\$verbose);


if (@ARGV != 1 && @ARGV != 2)
{
	print "bayesian analysis of splicing isoform structure from paired-end mRNA-seq data\n";
	print "Usage1: $prog [options] <read1.bed> <read2.bed>\n";
	print "Usage2: $prog [options] <read1_and_2.bed>\n";
	print "Usage3: $prog [options] <read1.sam> <read2.sam>\n";
	print "Usage4: $prog [options] <read1_and_2.sam>\n";
	print " gzip compressed input files with .gz extension are allowed\n";
	print " -sam                             : input files are in sam format\n";
	print " -uniq                            : keep only unique mapping in sam2bed conversion\n";
	print " --split-size             [int]   : number of SE. reads in each split ($splitSize)\n";
	print " -split                   [int]   : number of splits (will override --split-size)\n";
	print " -isoform                 [file]  : file name of splice isoform database\n";
	print " --use-pre-calc-score             : use pre-calculated isoform prior score\n";
	print " --clean-isoform                  : clean isoforms (will supress --use-pre-calc-score\n";
	print " -E                       [int]   : big exon size (default=$bigExonSize)\n";
	print " --use-pre-calc-size-dist [file]  : use pre-calculated size distribution file\n";
	print " -big                             : read number in each split is big (i.e. over ~3M pairs)\n";
	print " --library-type           [string]: library  type ([unstranded]|stranded-as|stranded-sa)\n";
	#print " --ss                             : separate the two strands\n";
	print " --size-dist-only                 : only estimate fragment size distribution\n";
	print " --print-singleton                : print reads even when they are not found to be a legitimate pair\n";
	print " -o                       [dir]   : output dir (default=$outputDir)\n"; 
	print " --keep-cache                     : keep cache files when the job is done\n";
	print " -v                               : verbose\n";
	exit (1);
}


my ($read1InFile, $read2InFile) = @ARGV;
Carp::croak "$read1InFile does not exists\n" unless -f $read1InFile;
if (@ARGV == 2)
{
	Carp::croak "$read2InFile does not exists\n" unless -f $read2InFile;
}

Carp::croak "$outputDir already exists" if -d $outputDir;

if ($sizeDistFile ne '')
{
	Carp::croak "size distribution file does not exists\n" unless -f $sizeDistFile;
}

if ($libraryType ne 'unstranded' && $libraryType ne 'stranded-as' && $libraryType ne 'stranded-sa')
{
	Carp::croak "incorrect library type: $libraryType\n";
}


print `date` if $verbose;

my $verboseFlag = $verbose ? '-v' : '';
my $bigFlag = $big ? "-big" : "";
#my $separateStrandFlag = $separateStrand ? "--ss" : '';

system ("mkdir $outputDir") unless -d $outputDir;

my $tmpDir = "$outputDir/tmp";
system ("mkdir $tmpDir");


#prepare input files
my $read1SingletonBedFile = "$tmpDir/read1.singleton.bed";
my $read2SingletonBedFile = "$tmpDir/read2.singleton.bed";

my $read1PairedBedFile = "$tmpDir/read1.paired.bed";
my $read2PairedBedFile = "$tmpDir/read2.paired.bed";
	
my ($n1s, $n2s, $np) = (0, 0, 0);
	
if ($sam)
{
	($n1s, $n2s, $np) = prepareSamFiles ($read1InFile, $read2InFile, $read1PairedBedFile, $read2PairedBedFile, $read1SingletonBedFile, $read2SingletonBedFile, $uniq);
}
else
{
	($n1s, $n2s, $np) = prepareBedFiles ($read1InFile, $read2InFile, $read1PairedBedFile, $read2PairedBedFile, $read1SingletonBedFile, $read2SingletonBedFile, $tmpDir);
}

print "n1s=$n1s, n2s=$n2s, np=$np\n" if $verbose;
Carp::croak "No pair was found. Check read names carefully\n" if $np == 0;


#split input bed files into chunks
#reads in the same pair will be in the same chunk
print "split paired reads ...\n" if $verbose;

if ($nsplit == 0)
{
	$nsplit = int (($np * 2 + $n1s + $n2s) / $splitSize);
	$nsplit = 1 if $nsplit == 0;
}

my $np_split = int ($np / $nsplit); 
$np_split++ if $np % $np_split != 0;

my $n1s_split = int ($n1s / $nsplit); 
$n1s_split++ if $n1s > 0 && $n1s % $n1s_split != 0;

my $n2s_split = int ($n2s / $nsplit); 
$n2s_split++ if $n2s > 0 && $n2s % $n2s_split != 0;


print "np_split = $np_split, n1s_split=$n1s_split, n2s_split=$n2s_split\n" if $verbose;

my $splitDir = "$tmpDir/split";
system ("mkdir $splitDir");

my $splitRead1BedStem = "$splitDir/r1.bed";
my $splitRead2BedStem = "$splitDir/r2.bed";

print "split paired 5' reads ...\n" if $verbose;
splitFileByRow ($read1PairedBedFile, $splitRead1BedStem, $np_split, $verbose, 0);

print "split unpaired 5' reads ...\n" if $verbose;
splitFileByRow ($read1SingletonBedFile, $splitRead1BedStem, $n1s_split, $verbose, 1);

print "split paired 3' reads ...\n" if $verbose;
splitFileByRow ($read2PairedBedFile, $splitRead2BedStem, $np_split, $verbose, 0);

print "split unpaired 3' reads ...\n" if $verbose;
splitFileByRow ($read2SingletonBedFile, $splitRead2BedStem, $n2s_split, $verbose, 1);

unlink $read1PairedBedFile, $read1SingletonBedFile, $read2PairedBedFile, $read2SingletonBedFile unless $keepCache;


#score isoforms
$usePreCalScore = 0 if $cleanIsoform;

my $bigExonBedFile = "$tmpDir/bigExon.bed";

if (not -f $isoformBedFile)
{
	Carp::croak "the isoform database file $isoformBedFile does not exist\n";
}

if ($cleanIsoform)
{
	print "clean isoforms ...\n" if $verbose;

	#replace the name so that they are simple and unique
	my $cleanIsoformBedFile = "$tmpDir/isoform.clean.bed";
	my $cmd = "awk 'BEGIN{i=0} {if(NF<10) {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"i\"\\t\"\$5\"\\t\"\$6} else {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"i\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12}; i=i+1}' $isoformBedFile > $cleanIsoformBedFile";
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
	#$cleanIsoformBedFile = $isoformBedFile;

	#remove redundancy
	my $cleanNrIsoformBedFile = "$tmpDir/isoform.clean.nr.bed";
	$cmd = "perl $cmdDir/bed2representative.pl $bigFlag $verboseFlag -s --same-block-count --same-intron -c $tmpDir/clean_bed_tmp_files $cleanIsoformBedFile $cleanNrIsoformBedFile";
	print $cmd, "\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

	$isoformBedFile = $cleanNrIsoformBedFile;
}

#get size distribution
if (-f $sizeDistFile)
{
	print "use provided insert size distribution file $sizeDistFile ...\n" if $verbose;
}
else
{
	print "estimating insert size distribution...\n" if $verbose;

	$sizeDistFile = "$outputDir/size_dist.txt";
	Carp::croak "$sizeDistFile already exists\n" if -f $sizeDistFile;

	system ("touch $sizeDistFile");

	my $bigExonBedFile = "$tmpDir/bigExon.bed";
	print "extracting exons with a size >= $bigExonSize to $bigExonBedFile...\n" if $verbose;
	
	my $cmd = "awk '{if ((\$10==1 || NF <= 6) && \$3-\$2>=$bigExonSize && \$1 != \"chrM\") {print \$0}}' $isoformBedFile > $bigExonBedFile";
	print $cmd, "\n" if $verbose;

	my $ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

	
	my %sizeDistHash;
	for (my $s = 0; $s < $nsplit; $s++)
	{
		print "processing split $s ...\n" if $verbose;

		print "retrieving genomic reads...\n" if $verbose;

		my $genomicRead1BedFile = "$tmpDir/read1.genomic.bed";
		my $genomicRead2BedFile = "$tmpDir/read2.genomic.bed";

		my $cmd = "awk '{if(NF<=9 || \$10==1) {print \$0}}' $splitRead1BedStem.$s > $genomicRead1BedFile";
		my $ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

		$cmd =    "awk '{if(NF<=9 || \$10==1) {print \$0}}' $splitRead2BedStem.$s > $genomicRead2BedFile";
		$ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

		my $pairOnBigExonBedFile = "$tmpDir/pair.bigexon.bed";

		$cmd = "perl $cmdDir/combine_paired_end.pl $verboseFlag $bigFlag --library-type $libraryType -cache $tmpDir/exon_size_tmp_files -gene $bigExonBedFile $genomicRead1BedFile $genomicRead2BedFile $pairOnBigExonBedFile";
		print $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
		
		system ("rm -rf $tmpDir/exon_size_tmp_files");

		#$cmd = "awk '{print \$3-\$2}' $pairOnBigExonBedFile | sort | uniq -c | awk '{print \$2\"\t\"\$1}' | sort -k 1 -n >> $sizeDistFile";
		#$ret = system ($cmd);
		#Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
		
		my $fin;
		open ($fin, "<$pairOnBigExonBedFile") || Carp::croak "cannot open file $pairOnBigExonBedFile to read\n";
		while (my $line =<$fin>)
		{
			next if $line =~/^\s*$/;
			chomp $line;
			my @cols = split ("\t", $line, 4);
			my $size = $cols[2] - $cols[1];
			$sizeDistHash{$size}++;
		}
		close ($fin);
	}

	
	#$cmd = "perl $cmdDir/uniqRow.pl -c sum -id 0 -value 1 $verboseFlag $sizeDistFile - | sort -k 1,1n > $sizeDistFile.tmp";
	#print $cmd, "\n" if $verbose;
	#$ret = system ($cmd);
	#Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
	
	#$cmd = "mv $sizeDistFile.tmp $sizeDistFile";
	#$ret = system ($cmd);
	#Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
	
	my $fout;
	open ($fout, ">$sizeDistFile") || Carp::croak "cannot open file $sizeDistFile to write\n";
	foreach my $size (sort {$a <=> $b} keys %sizeDistHash)
	{
		print $fout $size, "\t", $sizeDistHash{$size}, "\n";
	}
	close ($fout);
}

if ($sizeDistOnly)
{
	system ("rm -rf $tmpDir") unless $keepCache;
	print `date` if $verbose;
	exit (0);
}

my $scoredIsoformBedFile = "$tmpDir/isoform.scored.bed";
if ($usePreCalScore)
{
	print "using pre-calculated isoform score ...\n" if $verbose;
	$scoredIsoformBedFile = $isoformBedFile;
}
else
{
	print "scoring isoforms using observed junction reads ...\n" if $verbose;
	print "loading isoforms from $isoformBedFile ...\n" if $verbose;

	my $isoforms = readBedFile ($isoformBedFile, $verbose);

	#initialization
	map {$_->{'score'} = 0} @$isoforms;

	for (my $s = 0; $s < $nsplit; $s++)
	{
		print "processing split $s of $nsplit ...\n";
		print "retrieving junction reads ...\n" if $verbose;

		my $junctionRead1BedFile = "$tmpDir/read1.junction.bed";
		my $junctionRead2BedFile = "$tmpDir/read2.junction.bed";

		my $cmd = "awk '{if(\$10>1) {print \$0}}' $splitRead1BedStem.$s > $junctionRead1BedFile";
		print $cmd, "\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

		$cmd =    "awk '{if(\$10>1) {print \$0}}' $splitRead2BedStem.$s > $junctionRead2BedFile";
		print $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

		my $junctionReadBedFile = "$tmpDir/read.junction.bed";
		$cmd = "cat $junctionRead1BedFile $junctionRead2BedFile > $junctionReadBedFile";
		$ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;

		#library type not considered now, to be fixed
		$cmd = "perl $cmdDir/score_exon_trio.pl $verboseFlag $bigFlag -c $tmpDir/score_junction_tmp_files $isoformBedFile $junctionReadBedFile $scoredIsoformBedFile";
		print $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
		
		system ("rm -rf $tmpDir/score_junction_tmp_files");

		print "reading scored isoforms using split $s ...\n" if $verbose;
		my $splitScoredIsoforms = readBedFile ($scoredIsoformBedFile, $verbose);
		for (my $i = 0; $i < @$isoforms; $i++)
	   	{
			$isoforms->[$i]->{'score'} += $splitScoredIsoforms->[$i]->{'score'};
		}
		unlink $junctionRead1BedFile, $junctionRead2BedFile, $junctionReadBedFile;
	}
	print "writing scored isoforms to $scoredIsoformBedFile ...\n" if $verbose;

	writeBedFile ($isoforms, $scoredIsoformBedFile);
}



print "inferring structures of PE reads ...\n" if $verbose;

my $pairBedFile = "$outputDir/pair.gapless.bed";
my $printSingletonFlag = $printSingleton ? '--print-singleton' : '';

Carp::croak "$pairBedFile already exists\n" if -f $pairBedFile;

for (my $s = 0; $s < $nsplit; $s++)
{
	print "processing split $s ...\n" if $verbose;
	my $pairBedFileTmp = "$splitDir/pair.gapless.bed.$s";
	my $cmd = "perl $cmdDir/combine_paired_end.pl $verboseFlag $bigFlag --library-type $libraryType -cache $tmpDir/gapless_tmp_files -gene $scoredIsoformBedFile -use-prior -size-dist $sizeDistFile $printSingletonFlag $splitRead1BedStem.$s $splitRead2BedStem.$s $pairBedFileTmp";
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "Command [$cmd] crashed:$?\n" unless $ret == 0;
	system ("rm -rf $tmpDir/gapless_tmp_files");
	system ("cat  $pairBedFileTmp >> $pairBedFile");
}

system ("rm -rf $tmpDir") unless $keepCache;
print `date` if $verbose;



=head2 prepareSamFiles
#write input reads from sam files to four bed files
#read1/2 singleton
#read1/2 paired (in the same order)

#if read1 and read2 are provided in two separate files, they must be in the same order and the read names in a pair must be the same
#if read1 and read2 are provided in the same file, the line of read1 must by followed by the line of read2

=cut

sub prepareSamFiles
{
	my ($read1SamFile, $read2SamFile, $read1PairedBedFile, $read2PairedBedFile, $read1SingletonBedFile, $read2SingletonBedFile, $uniq) = @_;
	$uniq = 0 unless $uniq;	

	print "parsing input sam file(s) ...\n" if $verbose;

	my ($fin1, $fin2);

	if ($read1SamFile =~/\.gz$/)
	{
		open ($fin1, "gunzip -c $read1SamFile |") || Carp::croak "cannot open file $read1SamFile to read\n";
	}
	elsif ($read1SamFile =~/\.bz2$/)
    {
        open ($fin1, "bunzip2 -c $read1SamFile |") || Carp::croak "cannot open file $read1SamFile to read\n";
    }
	else
	{
		open ($fin1, $read1SamFile) || Carp::croak "cannot open file $read1SamFile to read\n";
	}

	if (-f $read2SamFile)
	{
		if ($read2SamFile =~/\.gz$/)
		{
			open ($fin2, "gunzip -c $read2SamFile |") || Carp::croak "cannot open file $read2SamFile to read\n";
		}
		elsif ($read2SamFile =~/\.bz2$/)
        {
            open ($fin2, "bunzip2 -c $read2SamFile |") || Carp::croak "cannot open file $read2SamFile to read\n";
        }
		else
		{
			open ($fin2, $read2SamFile) || Carp::croak "cannot open file $read2SamFile to read\n";
		}
	}

	my ($fout_r1_s, $fout_r1_p, $fout_r2_s, $fout_r2_p);
	my ($n1s, $n2s, $np) = (0, 0, 0);
	
	open ($fout_r1_s, ">$read1SingletonBedFile") || Carp::croak "cannot open file $read1SingletonBedFile to read\n";
	open ($fout_r2_s, ">$read2SingletonBedFile") || Carp::croak "cannot open file $read2SingletonBedFile to read\n";
	open ($fout_r1_p, ">$read1PairedBedFile") || Carp::croak "cannot open file $read1PairedBedFile to read\n";
	open ($fout_r2_p, ">$read2PairedBedFile") || Carp::croak "cannot open file $read2PairedBedFile to read\n";


	my $i = 0;
	while (my $read1_line = <$fin1>)
	{
		chomp $read1_line;
		next if $read1_line=~/^\s*$/;
		next if $read1_line=~/^\@/;

		print "$i ...\n" if $verbose && $i % 100000 == 0;
		$i++;
	
		my $read2_line = "";
		if (-f $read2SamFile)
		{
			while ($read2_line =<$fin2>)
			{
				next if $read2_line=~/^\s*$/;
				next if $read2_line=~/^\@/;
				last;
			}
		}
		else
		{
			while ($read2_line =<$fin1>)
			{
				next if $read2_line=~/^\s*$/;
				next if $read2_line=~/^\@/;
				last;
			}
		}

		my $read1_sam = lineToSam ($read1_line);
		my $read2_sam = lineToSam ($read2_line);
		
		#remove suffix if any
		$read1_sam->{"QNAME"} =~s/\/\d$//;
		$read2_sam->{"QNAME"} =~s/\/\d$//;

		Carp::croak "read1 and read2 does not have the same name:\n read1=$read1_line\nread2=$read2_line\n" 
		if $read1_sam->{"QNAME"} ne $read2_sam->{"QNAME"};

		my $read1_bed = samToBed ($read1_sam);
		my $read2_bed = samToBed ($read2_sam);

		$read1_bed = 0 if $read1_bed && $read1_bed->{"flagInfo"}->{'query_nomap'};
		$read2_bed = 0 if $read2_bed && $read2_bed->{"flagInfo"}->{'query_nomap'};
		
		if ($uniq)
		{
			if ($read1_bed)
			{
				$read1_bed = 0 unless $read1_sam->{"TAGS"} =~/XT:A:U/;
			}
			if ($read2_bed)
			{
				$read2_bed = 0 unless $read2_sam->{"TAGS"} =~/XT:A:U/;
			}
		}
	
		if ($read1_bed && $read2_bed)
		{
			$np++;
	
			print $fout_r1_p bedToLine ($read1_bed), "\n";
			print $fout_r2_p bedToLine ($read2_bed), "\n";
		}
		elsif ($read1_bed)
		{
			$n1s++;
			print $fout_r1_s bedToLine ($read1_bed), "\n";
		}
		elsif ($read2_bed)
		{
			$n2s++;
			print $fout_r2_s bedToLine ($read2_bed), "\n";
		}
	}

	close ($fin1);
	close ($fin2) if -f $read2SamFile;

	close ($fout_r1_s);
	close ($fout_r2_s);
	close ($fout_r1_p);
	close ($fout_r2_p);
	return ($n1s, $n2s, $np);

}





=head2 prepareBedFiles

#write input reads to four bed files
#read1/2 singleton
#read1/2 paired (in the same order)
#to make sure paired reads are in the same order, the input bed files need to be joint by read name first, which will be very time consuming

=cut

sub prepareBedFiles
{
	my ($read1BedFile, $read2BedFile, $read1PairedBedFile, $read2PairedBedFile, $read1SingletonBedFile, $read2SingletonBedFile, $cache) = @_;
	
	print "sort the original input files ...\n" if $verbose;

	my $sortedRead1BedFile = "$tmpDir/read1.sort.bed";
	my $sortedRead2BedFile = "$tmpDir/read2.sort.bed";
	my $jointReadBedFile = "$tmpDir/read.join.txt";

	if (-f $read2BedFile)
	{
		#read1 and read2 are provided in separate files
		my $cmd = "grep -v \"^track\" $read1BedFile | awk '{print \$4\"\\tr1\\t\"\$0}' | sort -T $cache -k 1b,1 > $sortedRead1BedFile";
		$cmd = "gunzip -c $read1BedFile | grep -v \"^track\" | awk '{print \$4\"\\tr1\\t\"\$0}' | sort -T $cache -k 1b,1 > $sortedRead1BedFile" if $read1BedFile=~/\.gz$/;
		$cmd = "bunzip2 -c $read1BedFile | grep -v \"^track\" | awk '{print \$4\"\\tr1\\t\"\$0}' | sort -T $cache -k 1b,1 > $sortedRead1BedFile" if $read1BedFile=~/\.bz2$/;

		print $cmd, "\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

		$cmd =    "grep -v \"^track\" $read2BedFile | awk '{print \$4\"\\tr2\\t\"\$0}' | sort -T $cache -k 1b,1 > $sortedRead2BedFile";
		$cmd = "gunzip -c $read2BedFile | grep -v \"^track\" | awk '{print \$4\"\\tr2\\t\"\$0}' | sort -T $cache -k 1b,1 > $sortedRead2BedFile" if $read2BedFile=~/\.gz$/;		
		$cmd = "bunzip2 -c $read2BedFile | grep -v \"^track\" | awk '{print \$4\"\\tr2\\t\"\$0}' | sort -T $cache -k 1b,1 > $sortedRead2BedFile" if $read2BedFile=~/\.bz2$/;

		print $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
	}
	else
	{
		#read1 and 2 are provided in the same file and distinguished by /1 and /2
		my $cmd = "grep -v \"^track\" $read1BedFile | grep -P \"\\/1\\t\" | sed 's/\\/1\\t/\\t/g' | awk '{print \$4\"\\tr1\\t\"\$0}'| sort -T $cache -k 1b,1 > $sortedRead1BedFile";
		$cmd = "gunzip -c $read1BedFile | grep -v \"^track\" | grep -P \"\\/1\\t\" | sed 's/\\/1\\t/\\t/g' | awk '{print \$4\"\\tr1\\t\"\$0}'| sort -T $cache -k 1b,1 > $sortedRead1BedFile"
		if $read1BedFile=~/\.gz$/;
		
		$cmd = "bunzip2 -c $read1BedFile | grep -v \"^track\" | grep -P \"\\/1\\t\" | sed 's/\\/1\\t/\\t/g' | awk '{print \$4\"\\tr1\\t\"\$0}'| sort -T $cache -k 1b,1 > $sortedRead1BedFile"
		if $read1BedFile=~/\.bz2$/;
	
		print $cmd, "\n" if $verbose;
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
	
		$cmd =    "grep -v \"^track\" $read1BedFile | grep -P \"\\/2\\t\" | sed 's/\\/2\\t/\\t/g' | awk '{print \$4\"\\tr2\\t\"\$0}'| sort -T $cache -k 1b,1 > $sortedRead2BedFile";
		$cmd = "gunzip -c $read1BedFile | grep -v \"^track\" | grep -P \"\\/2\\t\" | sed 's/\\/2\\t/\\t/g' | awk '{print \$4\"\\tr2\\t\"\$0}'| sort -T $cache -k 1b,1 > $sortedRead2BedFile"
		if $read1BedFile=~/\.gz$/;

		$cmd = "bunzip2 -c $read1BedFile | grep -v \"^track\" | grep -P \"\\/2\\t\" | sed 's/\\/2\\t/\\t/g' | awk '{print \$4\"\\tr2\\t\"\$0}'| sort -T $cache -k 1b,1 > $sortedRead2BedFile"
		if $read1BedFile=~/\.bz2$/;

		print $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
	}


	my $cmd = "join -j 1 -a 1 -a 2 $sortedRead1BedFile $sortedRead2BedFile > $jointReadBedFile";
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

	unlink $sortedRead1BedFile, $sortedRead2BedFile;

	print "parsing joint output file ...\n" if $verbose;

	my $fin;
	open ($fin, "<$jointReadBedFile") || Carp::croak "cannot open file $jointReadBedFile to read\n";

	my ($fout_r1_s, $fout_r1_p, $fout_r2_s, $fout_r2_p);
	my ($n1s, $n2s, $np) = (0, 0, 0);
	
	open ($fout_r1_s, ">$read1SingletonBedFile") || Carp::croak "cannot open file $read1SingletonBedFile to read\n";
	open ($fout_r2_s, ">$read2SingletonBedFile") || Carp::croak "cannot open file $read2SingletonBedFile to read\n";
	open ($fout_r1_p, ">$read1PairedBedFile") || Carp::croak "cannot open file $read1PairedBedFile to read\n";
	open ($fout_r2_p, ">$read2PairedBedFile") || Carp::croak "cannot open file $read2PairedBedFile to read\n";

	my $i = 0;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
	
		print "$i ...\n" if $verbose && $i % 100000 == 0;
		$i++;
	
		#print $line, "\n";
	
		my $read1_line = "";
		my $read2_line = "";
		my @cols = split (/\s+r1\s+/, $line);
		Carp::croak "inconsistency found in $line\n" if @cols > 2;
	
		if (@cols == 1)
		{
			#no read1, must have read2
			my $junk;
			($junk, $read2_line) = split (/\s+r2\s+/, $line);
		}
		else
		{
			#there is read1 in the line
			my $line = $cols[1];
			($read1_line, $read2_line) = split (/\s+r2\s+/, $line);
			$read2_line = "" unless $read2_line;
		}
		
		$read1_line =~s/\s+/\t/g;
		$read2_line =~s/\s+/\t/g;
	
		if ($read1_line ne '' && $read2_line ne '')
		{
			$np++;
	
			print $fout_r1_p $read1_line, "\n";
			print $fout_r2_p $read2_line, "\n";
		}
		elsif ($read1_line ne '')
		{
			$n1s++;
			print $fout_r1_s $read1_line, "\n";
		}
		elsif ($read2_line ne '')
		{
			$n2s++;
			print $fout_r2_s $read2_line, "\n";
		}
		else
		{
			Carp::croak "inconsistency found in line $line\n";
		}
	}

	close ($fin);
	close ($fout_r1_s);
	close ($fout_r2_s);
	close ($fout_r1_p);
	close ($fout_r2_p);

	unlink $jointReadBedFile;

	return ($n1s, $n2s, $np);
}






