#!/usr/bin/perl -w


use strict;
use Getopt::Long;

use Carp;
use File::Basename;
use Data::Dumper;

use MyConfig;

my $prog = basename ($0);
my $verbose = 0;

my @ARGV0 = @ARGV;


my $big = 0;
my $confFile = "";
my $dbkey = "";
my $weight = 0;
my $weightAvg = 0;
my $separateStrand = 0;

my $cache = getDefaultCache ($prog);
my $keepCache = 0;
my $combineOutput = 0;

my $get_cass = 0;
my $get_taca = 0;
my $get_alt5 = 0;
my $get_alt3 = 0;
my $get_mutx = 0;
my $get_iret = 0;
my $get_alts = 0;
my $get_altt = 0;

GetOptions (
	'conf=s'=>\$confFile,
	'dbkey=s'=>\$dbkey,
	'big'=>\$big,
	'weight'=>\$weight,
	'weight-avg'=>\$weightAvg,
	'ss'=>\$separateStrand,
	'cass'=>\$get_cass,
	'taca'=>\$get_taca,
	'alt5'=>\$get_alt5,
	'alt3'=>\$get_alt3,
	'mutx'=>\$get_mutx,
	'iret'=>\$get_iret,
	'alts'=>\$get_alts,
	'altt'=>\$get_altt,
	'c:s'=>\$cache,
	'keep-cache'=>\$keepCache,
	'v'=>\$verbose);


if (@ARGV != 2)
{
	print "Count reads on each isoform of AS events\n";

	print "Usage1: $prog [options] <in.bed> <out.dir>\n";
	print "Usage2: $prog [options] --all-in-one <in.bed> <out.txt>\n"; 
	print " -conf  [file]  : configuration file with input datasets\n";
	print " -dbkey [string]: genome build name (hg18|mm9)\n";
	print " -big           : big file\n";
	print " -weight        : weight tags using the score column\n";
	print " --weight-avg   : calculate weighted average instead of weighted sum\n";
	print " --ss           : consider the two strands separately\n";
	print " -cass          : annotate cass\n";
	print " -taca          : annotate taca\n";
	print " -alt5          : annotate alt5\n";
	print " -alt3          : annotate alt3\n";
	print " -mutx          : annotate mutx\n";
	print " -iret          : annotate iret\n";
	print " -alts          : annotate alts\n";
	print " -altt          : annotate altt\n";
	print " -c     [string]: cache dir ($cache) \n";
	print " --all-in-one   : combine all output in one file\n";
	print " --keep-cache   : keep cache when the job is done\n";
	print " -v             : verbose\n";
	exit (0);
}

print "CMD = $prog ", join (' ', @ARGV0), "\n" if $verbose;

my %analyses;

my $verboseFlag = $verbose ? '-v' : '';
my $bigFlag = $big ? '-big' : '';
my $weightFlag = $weight ? '-weight' : '';
my $weightAvgFlag = $weightAvg ? '-mean' : '';
my $ssFlag = $separateStrand ? '--ss' : '';

my $keepCacheFlag = $keepCache ? '--keep-cache' : '';

$analyses{'cass'} = 1 if $get_cass;
$analyses{'taca'} = 2 if $get_taca;
$analyses{'alt5'} = 3 if $get_alt5;
$analyses{'alt3'} = 4 if $get_alt3;
$analyses{'mutx'} = 5 if $get_mutx;
$analyses{'iret'} = 6 if $get_iret;
$analyses{'alts'} = 7 if $get_alts;
$analyses{'altt'} = 8 if $get_altt;


Carp::croak "no AS type chosen\n" if (keys %analyses) == 0;


my ($inBedFile, $outDir) = @ARGV;

my $outFile = "";
$outFile = $outDir if $combineOutput;


#unlink $outFile if -f $outFile;
system ("mkdir $cache");

if ($combineOutput == 0)
{
	my $ret = system("mkdir $outDir");
	Carp::croak "cannot create dir $outDir: $?\n" if $ret != 0;
}

print "separate genomic and junction reads ...\n" if $verbose;

#my $genomicTagBedFile = "$cache/tag.genome.bed";
#my $junctionTagBedFile = "$cache/tag.junction.bed";

#my $cmd = "grep -v \"^track\" $inBedFile | awk '{if(NF<=9 || \$10==1) {print \$0}}' > $genomicTagBedFile";
#my $ret = system ($cmd);
#print "CMD $cmd failed: $?\n" if $ret != 0;


#$cmd = "grep -v \"^track\" $inBedFile | awk '{if(NF==12 && \$10>1) {print \$0}}' > $junctionTagBedFile";
#$ret = system ($cmd);
#print "CMD $cmd failed: $?\n" if $ret != 0;

if ($combineOutput)
{
	my $header = "#" . join ("\t", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "isoformIDs", "isoform1Tags", "isoform2Tags", "extraInfoInOtherColumns");
	system ("echo \"$header\" > $outFile");
}

my $cmdDir = dirname ($0);
my $locationInfo = getLocationInfo ($confFile, $dbkey, "as");
foreach my $asType (sort {$analyses{$a} <=> $analyses{$b}} keys %analyses)
{
	print "analyze $asType ...\n" if $verbose;
	Carp::croak "annotation of AS type $asType does not exist\n" unless exists $locationInfo->{$asType} && -f  $locationInfo->{$asType};

	my $asBedFile =  $locationInfo->{$asType};
	my $tmpOutFile = "$cache/$asType.count.txt";

	my $typeCache = "$cache/cache_$asType";

	my $cmd = "perl $cmdDir/summarize_splicing.pl -type $asType $bigFlag $verboseFlag $keepCacheFlag -c $typeCache $ssFlag $weightFlag $weightAvgFlag $asBedFile $inBedFile $tmpOutFile";
	print "$cmd\n" if $verbose;
	my $ret = system ($cmd);
	print "CMD $cmd failed: $?\n" if $ret != 0;
	
	$cmd = $combineOutput ? "cat $tmpOutFile | grep -v \"^#\" >> $outFile" : "mv $tmpOutFile $outDir/";
	print "$cmd\n" if $verbose;
	$ret = system ($cmd);
	print "CMD $cmd failed: $?\n" if $ret != 0;
}
system ("rm -rf $cache") unless $keepCache;



sub getLocationInfo
{
	my ($conf, $dbkey, $analysis) = @_;

	my $fin;
	
	open ($fin, "<$conf") || Carp::croak "cannot open file $conf\n";
	
	my %ret;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;

		my ($db, $ana, $path, $type) = split (/\s+/, $line);
		$type = $ana unless $type;

		if ($db eq $dbkey && $ana eq $analysis)
		{
			Carp::croak "$path does not exist\n" unless -f $path;
			$ret{$type} = $path;
		}
	}
	close ($fin);
	return \%ret;
}

