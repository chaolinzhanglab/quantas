#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Data::Dumper;
use File::Basename;
use Getopt::Long;

use MyConfig;
use Common;
use Bed;


my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $configFile = "";
my $asType = "cass";
my $separateStrand = 0;

my $minCoverage = 10;
my $maxStd = 0.1;
my $naString = "";
my $id2gene2symbolFile = "";

my $verbose = 0;

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

GetOptions (
	"conf:s"=>\$configFile,
	"type:s"=>\$asType,
	"ss"=>\$separateStrand,
    "min-cov:i"=>\$minCoverage,
    "max-std:f"=>\$maxStd,
    "na-string:s"=>\$naString,
	"id2gene2symbol:s"=>\$id2gene2symbolFile,
	"c:s"=>\$cache,
	"keep-cache"=>\$keepCache,
	"v"=>\$verbose);


if (@ARGV != 4)
{
	print "generate exon inclusion matrix from GTEx junction read count matrix\n";
	print "$prog [options] <as.bed> <jcn.bed> <jcn_count.matrix.txt> <out.matrix.txt>\n";
	print " -conf         [string]: the first column is the dir or file name, and the second column is the group name\n";
	print " -type         [string] : AS type ([cass]|taca|alt5|alt3|mutx)\n";
	print " -ss                    : separate strand when matching junctions\n";
    print " --min-cov     [int]    : min coverage ($minCoverage)\n";
    print " --max-std     [float]  : max standard deviation ($maxStd)\n";
    print " --na-string   [string] : na string (default:empty)\n";
	print " --id2gene2symbol [file]: mapping file of id to gene to symbol\n";	
	print " -c            [string] : cache dir ($cache)\n";
	print " --keep-cache           : keep cache dir when done\n";
	print " -v                     : verbose\n";
	exit (1);
}


my ($ASBedFile, $junctionBedFile, $junctionMatrixFile, $outFile) = @ARGV;
my $verboseFlag = $verbose ? '-v' : '';

my $ret = system ("mkdir $cache");
Carp::croak "cannot open cache dir $cache:$?\n" unless $ret == 0; 


my $newJunctionMatrixFile = "$cache/jcn_count.matrix.txt";

if ($configFile ne '' && (-f $configFile))
{
	print "combine replicates in junction matrix ...\n" if $verbose;
	combineJunctionCountMatrix ($junctionMatrixFile, $configFile, $newJunctionMatrixFile);
}
else
{
	my $tmp = getFullPath ($junctionMatrixFile);
	my $cmd = "ln -s $tmp $newJunctionMatrixFile";
	print "$cmd\n" if $verbose;

	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed:$?\n" unless $ret == 0;
}


print "index junction count matrix file $newJunctionMatrixFile...\n" if $verbose;

my $junctionMatrixIndexFile = "$cache/junction.matrix.idx";

#index the junction count data using the junction id as key
#Also, we skip 2 lines because the junction matrix files from GTex has two extra lines. 
#The headline with sample names are not skipped here.  The key for the headline is 'junction_id'
#need to be fixed to allow more flexibility later

my $cmd = "perl $cmdDir/indexRow.pl $verboseFlag -key 0 --skip 2 $newJunctionMatrixFile $junctionMatrixIndexFile";
print "$cmd\n" if $verbose;

$ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" unless $ret == 0;


print "loading junction count index...\n" if $verbose;

my %junctionMatrixIndexHash;
my $fin;

open ($fin, "<$junctionMatrixIndexFile") || Carp::croak "cannot open file $junctionMatrixIndexFile to read\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	next if $line =~/^#/;

	my ($junctionId, $pointer) = split ("\t", $line);
	$junctionMatrixIndexHash {$junctionId} = $pointer;
}
close ($fin);


print "loading junction bed file $junctionBedFile ...\n" if $verbose;

my %junctionHash;

open ($fin, "<$junctionBedFile") || Carp::croak "cannot open file $junctionBedFile to read\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	
	my $junction = lineToBed ($line);	
	
	my $junctionId = $junction->{"name"};
	Carp::croak "junction $junctionId is not found in the junction count matrix file\n" unless exists $junctionMatrixIndexHash{$junctionId};

	my $junctionKey = $junction->{'chrom'} . ':' . $junction->{'chromStart'} . '-' . $junction->{'chromEnd'};		
	if ($separateStrand)
	{
		Carp::croak "no strand information in junction $junctionId: ", Dumper ($junction), "\n" unless exists $junction->{'strand'};
		$junctionKey .= ":" . $junction->{'strand'};
	}

	$junctionHash{$junctionKey} = {id=>$junctionId, pointer=>$junctionMatrixIndexHash{$junctionId}};
}

close ($fin);


################################


print "extract introns from $ASBedFile ...\n" if $verbose;
my $ASIntronBedFile = "$cache/as.intron.bed";

$cmd = "perl $cmdDir/gene2ExonIntron.pl $verboseFlag -oi $ASIntronBedFile $ASBedFile";
print $cmd, "\n" if $verbose;

$ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;


#organize introns for each AS events
my $ASIntrons = readBedFile ($ASIntronBedFile, $verbose);
my %ASHash;

foreach my $junction (@$ASIntrons)
{
	my $name = $junction->{'name'};
	my $junctionKey = $junction->{'chrom'} . ':' . $junction->{'chromStart'} . '-' . $junction->{'chromEnd'};
	if ($separateStrand)
	{
		Carp::croak "no strand information : ", Dumper ($junction), "\n" unless exists $junction->{'strand'};
		$junctionKey .= ":" . $junction->{'strand'};
	}	

	my ($asId, $isoformId, $evi, $junctionId);
	
	if ($name =~/^(.*?)\[(.*?)\]\[(.*?)\].*?\_(\d+)$/)
	{
		#to accomodate the pattern of new mutx ids
		($asId, $isoformId, $evi, $junctionId) = ($1, $2, $3, $4);
	}
	else
	{
		Carp::croak "incorrect format: $name\n";
	}


	if ($asType eq 'cass')
	{
		#$name = join ("", $asId, "[INC/SKIP][", $evi, "]");
		if ($isoformId eq 'INC')
		{
			$ASHash{$asId}->{"INC/SKIP"}->{"junction"}->{"INC"}->{$junctionId} = $junctionKey;
		}
		else
		{
			$ASHash{$asId}->{"INC/SKIP"}->{"junction"}->{"SKIP"} = $junctionKey;
		}
	}
	elsif ($asType eq 'taca')
	{
		if ($isoformId eq 'INC')
		{
			push @{$ASHash{$asId}->{"INC/SKIP"}->{"junction"}->{"INC"}}, $junctionKey;
		}
		else
		{
			$ASHash{$asId}->{"INC/SKIP"}->{"junction"}->{"SKIP"} = $junctionKey;
		}
	}
	elsif ($asType eq 'alt5' || $asType eq 'alt3')
	{
		my @evis = split (/\//, $evi);
		$isoformId =~/(\d+)$/;
		my $isoformIdx = $1;
		
		#enumerate all possible junctions with other isoforms	
		#the smaller isoform idx is always on the left
		for (my $i = 0; $i < $isoformIdx; $i++)
		{
			my $isoformIdPair = "A$i/A$isoformIdx";
			$ASHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId} = $junctionKey;
		}
		for (my $i = $isoformIdx + 1; $i < @evis; $i++)
		{
			my $isoformIdPair = "A$isoformIdx/A$i";
			$ASHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId} = $junctionKey;
		}
	}
	elsif ($asType eq 'mutx')
	{
		my @evis = split (/\//, $evi);
		$isoformId =~/(\d+)$/;
		my $isoformIdx = $1;
	
		#enumerate all possible junctions with other isoforms
		#the smaller isoform idx is always on the left
	
		for (my $i = 0; $i < $isoformIdx; $i++)
		{
			my $isoformIdPair = "M$i/M$isoformIdx";
			$ASHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId}->{$junctionId} = $junctionKey;
		}
	
		for (my $i = $isoformIdx + 1; $i < @evis; $i++)
		{
			my $isoformIdPair = "M$isoformIdx/M$i";
			$ASHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId}->{$junctionId} = $junctionKey;
		}
	}
	else
	{
		Carp::croak "wrong AS type: $asType\n";
	}
}


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
elsif ($id2gene2symbolFile ne '')
{
    Carp::croak "cannot open file $id2gene2symbolFile to read\n";
}


my $n = keys %id2gene2symbolHash;

print "$n mapping entries loaded\n" if $verbose;



############################################################
#output
############################################################


print "quantify exon inclusion ...\n" if $verbose;

my $fout;
open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

open ($fin, "<$newJunctionMatrixFile") || Carp::croak "cannot open file $newJunctionMatrixFile to read\n";
my $pointer = exists $junctionMatrixIndexHash{'junction_id'} ? $junctionMatrixIndexHash{'junction_id'} : $junctionMatrixIndexHash{'Name'};

seek ($fin, $pointer, 0); #go to the point
my $line = <$fin>;
chomp $line;
my @sampleNames = split ("\t", $line);
shift @sampleNames; shift @sampleNames; #remove the first two columns


if (-f $id2gene2symbolFile)
{
	print $fout join ("\t", "#event_id", "NAME", @sampleNames), "\n";
}
else
{
	print $fout join ("\t", "#event_id", @sampleNames), "\n";
}


my $ASEvents = readBedFile ($ASBedFile, $verbose);

$n = @$ASEvents;
print "$n entries loaded\n" if $verbose;

my %ASEoutput;

my $iter = 0;
foreach my $e (@$ASEvents)
{

	print "$iter ...\n" if $verbose && $iter % 1000 == 0;
	$iter++;

	my $name = $e->{"name"};

	my ($asId, $isoformId) = ("", "");
	
	if ($name =~/^(.*?)\[(.*?)\]/)
	{
		$asId = $1;
		$isoformId = $2;
	}
	else
	{
		Carp::croak "incorrect name format: $name\n";
	}

	#$asId =~/^(\w\w)/;
	#my $asType = $1;

	next if (($asType eq 'cass' || $asType eq 'taca') && $isoformId eq 'SKIP');

	next if exists $ASEoutput{$asId};	#already dumpped

	my $ASEs = $ASHash{$asId};

	my @psi;
	foreach my $isoformIdPair (sort keys %$ASEs)
	{
		my $ase = $ASEs->{$isoformIdPair};

		if ($asType eq 'cass')
		{
			
			my $inc1JunctionKey = exists $ase->{"junction"} && exists $ase->{"junction"}->{"INC"} && exists $ase->{"junction"}->{"INC"}->{0} ? $ase->{"junction"}->{"INC"}->{0} : "";
			my $inc2JunctionKey = exists $ase->{"junction"} && exists $ase->{"junction"}->{"INC"} && exists $ase->{"junction"}->{"INC"}->{1} ? $ase->{"junction"}->{"INC"}->{1} : "";
			my $skipJunctionKey = exists $ase->{"junction"} && exists $ase->{"junction"}->{"SKIP"} ? $ase->{"junction"}->{"SKIP"} : "";
	
			next if $inc1JunctionKey eq '' || $inc2JunctionKey eq '' || $skipJunctionKey eq '';
			next unless exists $junctionHash{$inc1JunctionKey} && exists $junctionHash{$inc2JunctionKey} && exists $junctionHash{$skipJunctionKey};

			my $inc1Junction = $junctionHash{$inc1JunctionKey};
			my $inc2Junction = $junctionHash{$inc2JunctionKey};
			my $skipJunction = $junctionHash{$skipJunctionKey};

			seek ($fin, $inc1Junction->{'pointer'}, 0);
			$line = <$fin>;
			my @inc1JunctionCount = split ("\t", $line);
			shift @inc1JunctionCount; shift @inc1JunctionCount;

			seek ($fin, $inc2Junction->{'pointer'}, 0);
			$line = <$fin>;
			my @inc2JunctionCount = split ("\t", $line);
			shift @inc2JunctionCount; shift @inc2JunctionCount;

			seek ($fin, $skipJunction->{'pointer'}, 0);
			$line = <$fin>;
			my @skipJunctionCount = split ("\t", $line);
			shift @skipJunctionCount; shift @skipJunctionCount;

			for (my $i = 0; $i < @inc1JunctionCount; $i++)
			{
				my $in = ($inc1JunctionCount[$i] + $inc2JunctionCount[$i]);
				my $ex = $skipJunctionCount[$i] * 2;

	        	my $n = $in + $ex;
				if ($n < $minCoverage)
       			{
            		$psi[$i] = $naString;
            		next;
        		}

 				my $p = $in / $n;
        		my $std = sqrt ($p * (1-$p) / $n);
        		if ($std > $maxStd)
        		{
            		$psi[$i] = $naString;
        		}
        		else
        		{
            		$psi[$i] = $p;
        		}
			}
		}
		elsif ($asType eq 'taca')
		{

			my $incJunctionKeys = exists $ase->{"junction"} && exists $ase->{"junction"}->{"INC"} ? $ase->{"junction"}->{"INC"} : "";
			my $skipJunctionKey = exists $ase->{"junction"} && exists $ase->{"junction"}->{"SKIP"} ? $ase->{"junction"}->{"SKIP"} : "";
	
			next if $skipJunctionKey eq '' || $incJunctionKeys eq '';

			my $good = 1;
			foreach my $incJunctionKey (@$incJunctionKeys)
			{
				if (not exists $junctionHash{$incJunctionKey})
				{
					$good = 0;
					last;
				}
			}
			next unless $good;
			next unless exists $junctionHash{$skipJunctionKey};

			my $numIncJunctions = @$incJunctionKeys;
			my @incJunctionCount;

			#inclusion junctions
			foreach my $incJunctionKey (@$incJunctionKeys)
            {
				my $incJunction = $junctionHash{$incJunctionKey};

				seek ($fin, $incJunction->{'pointer'}, 0);
				$line = <$fin>;
				my @iCount = split ("\t", $line);
				shift @iCount; shift @iCount;
				for (my $i = 0; $i < @iCount; $i++)
				{
					$incJunctionCount[$i] += $iCount[$i];
				}
			}

			my $skipJunction = $junctionHash{$skipJunctionKey};
			seek ($fin, $skipJunction->{'pointer'}, 0);
			$line = <$fin>;
			my @skipJunctionCount = split ("\t", $line);
			shift @skipJunctionCount; shift @skipJunctionCount;

			for (my $i = 0; $i < @incJunctionCount; $i++)
			{
				my $in = $incJunctionCount[$i];
				my $ex = $skipJunctionCount[$i] * $numIncJunctions;
	        	
				my $n = $in + $ex;
				if ($n < $minCoverage)
       			{
            		$psi[$i] = $naString;
            		next;
        		}

 				my $p = $in / $n;
        		my $std = sqrt ($p * (1-$p) / $n);
        		if ($std > $maxStd)
        		{
            		$psi[$i] = $naString;
        		}
        		else
        		{
            		$psi[$i] = $p;
        		}
			}
		}
		elsif ($asType eq 'alt5' || $asType eq 'alt3')
		{
			my ($isoformId1, $isoformId2) = split (/\//, $isoformIdPair);
			$isoformId1=~/(\d+)$/;
			my $isoformIdx1 = $1;

			$isoformId2=~/(\d+)$/;
			my $isoformIdx2 = $1;
			
			my $proximalIsoformId = ($isoformIdx1 > $isoformIdx2) ? $isoformId1 : $isoformId2;
			my $distalIsoformId = ($isoformIdx1 < $isoformIdx2) ? $isoformId1 : $isoformId2;
			
			($proximalIsoformId, $distalIsoformId) = ($distalIsoformId, $proximalIsoformId) if $asType eq 'alt3';

			my $proximalJunctionKey = exists $ase->{"junction"} && exists $ase->{"junction"}->{$proximalIsoformId} ? $ase->{"junction"}->{$proximalIsoformId} : "";
			my $distalJunctionKey = exists $ase->{"junction"} && exists $ase->{"junction"}->{$distalIsoformId} ? $ase->{"junction"}->{$distalIsoformId} : "";
	
			next if $proximalJunctionKey eq '' || $distalJunctionKey eq '';

			next unless exists $junctionHash{$proximalJunctionKey} && exists $junctionHash{$distalJunctionKey};

			my $proximalJunction = $junctionHash{$proximalJunctionKey};
			seek ($fin, $proximalJunction->{'pointer'}, 0);
			$line = <$fin>;
			my @proximalJunctionCount = split ("\t", $line);
			shift @proximalJunctionCount; shift @proximalJunctionCount;

			my $distalJunction = $junctionHash{$distalJunctionKey};
			seek ($fin, $distalJunction->{'pointer'}, 0);
			$line = <$fin>;
			my @distalJunctionCount = split ("\t", $line);
			shift @distalJunctionCount; shift @distalJunctionCount;

			for (my $i = 0; $i < @proximalJunctionCount; $i++)
			{
				my $in = $proximalJunctionCount[$i];
				my $ex = $distalJunctionCount[$i];
				
				my $n = $in + $ex;
				if ($n < $minCoverage)
       			{
            		$psi[$i] = $naString;
            		next;
        		}

 				my $p = $in / $n;
        		my $std = sqrt ($p * (1-$p) / $n);
        		if ($std > $maxStd)
        		{
            		$psi[$i] = $naString;
        		}
        		else
        		{
            		$psi[$i] = $p;
        		}
			}
		}
		elsif ($asType eq 'mutx')
		{
			my ($isoformId1, $isoformId2) = split (/\//, $isoformIdPair);

			$isoformId1=~/(\d+)$/;
			my $isoformIdx1 = $1;

			$isoformId2=~/(\d+)$/;
			my $isoformIdx2 = $1;

			my $leftIsoformId = ($isoformIdx1 < $isoformIdx2) ? $isoformId1 : $isoformId2;
			my $rightIsoformId = ($isoformIdx1 > $isoformIdx2) ? $isoformId1 : $isoformId2;

			my $leftExonJunction1Key = exists $ase->{"junction"} && exists $ase->{"junction"}->{$leftIsoformId} && exists $ase->{"junction"}->{$leftIsoformId}->{0} ? $ase->{"junction"}->{$leftIsoformId}->{0} : "";
			my $leftExonJunction2Key = exists $ase->{"junction"} && exists $ase->{"junction"}->{$leftIsoformId} && exists $ase->{"junction"}->{$leftIsoformId}->{1} ? $ase->{"junction"}->{$leftIsoformId}->{1} : "";
			
			my $rightExonJunction1Key = exists $ase->{"junction"} && exists $ase->{"junction"}->{$rightIsoformId} && exists $ase->{"junction"}->{$rightIsoformId}->{0} ? $ase->{"junction"}->{$rightIsoformId}->{0} : "";
			my $rightExonJunction2Key = exists $ase->{"junction"} && exists $ase->{"junction"}->{$rightIsoformId} && exists $ase->{"junction"}->{$rightIsoformId}->{1} ? $ase->{"junction"}->{$rightIsoformId}->{1} : "";
			
			next if $leftExonJunction1Key eq '' || $leftExonJunction2Key eq '' || $rightExonJunction1Key eq '' || $rightExonJunction2Key eq '';
			next unless exists $junctionHash{$leftExonJunction1Key} && exists $junctionHash{$leftExonJunction2Key} && exists $junctionHash{$rightExonJunction1Key} && exists $junctionHash{$rightExonJunction2Key};

			my $leftExonJunction1 = $junctionHash{$leftExonJunction1Key};
			seek ($fin, $leftExonJunction1->{'pointer'}, 0);
			$line = <$fin>;
			my @leftExonJunction1Count = split ("\t", $line);
			shift @leftExonJunction1Count; shift @leftExonJunction1Count;

			my $leftExonJunction2 = $junctionHash{$leftExonJunction2Key};
			seek ($fin, $leftExonJunction2->{'pointer'}, 0);
			$line = <$fin>;
			my @leftExonJunction2Count = split ("\t", $line);
			shift @leftExonJunction2Count; shift @leftExonJunction2Count;

			my $rightExonJunction1 = $junctionHash{$rightExonJunction1Key};
			seek ($fin, $rightExonJunction1->{'pointer'}, 0);
			$line = <$fin>;
			my @rightExonJunction1Count = split ("\t", $line);
			shift @rightExonJunction1Count; shift @rightExonJunction1Count;

			my $rightExonJunction2 = $junctionHash{$rightExonJunction2Key};
			seek ($fin, $rightExonJunction2->{'pointer'}, 0);
			$line = <$fin>;
			my @rightExonJunction2Count = split ("\t", $line);
			shift @rightExonJunction2Count; shift @rightExonJunction2Count;

			my @psi;
			for (my $i = 0; $i < @leftExonJunction1Count; $i++)
			{
				my $in = $leftExonJunction1Count[$i] + $leftExonJunction2Count[$i];
				my $ex = $rightExonJunction1Count[$i] + $rightExonJunction2Count[$i];
	
				my $n = $in + $ex;
				if ($n < $minCoverage)
       			{
            		$psi[$i] = $naString;
            		next;
        		}

 				my $p = $in / $n;
        		my $std = sqrt ($p * (1-$p) / $n);
        		if ($std > $maxStd)
        		{
            		$psi[$i] = $naString;
        		}
        		else
        		{
            		$psi[$i] = $p;
        		}
			}
		}

		my $gene2symbol = exists $id2gene2symbolHash{$e->{'name'}} ? $id2gene2symbolHash{$e->{'name'}} : "NA//NA";
		if (-f $id2gene2symbolFile)
    	{
			print $fout join ("\t", $e->{"name"}, $gene2symbol, @psi), "\n";
		}
		else
		{
			print $fout join ("\t", $e->{"name"}, @psi), "\n";
		}
	}	
	$ASEoutput{$asId} = 1;
}

close ($fin);
close ($fout);

system ("rm -rf $cache") unless $keepCache;


sub readConfigFile
{
    my $configFile = $_[0];
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
    }
    close ($fin);
    return \%groups;
}


sub combineJunctionCountMatrix
{

	my ($junctionMatrixFile, $configFile, $newJunctionMatrixFile) = @_;

    print "loading configuration file from $configFile ...\n" if $verbose;
    Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
    my $groups = readConfigFile ($configFile);

	my ($fin, $fout);

	open ($fin, "<$junctionMatrixFile") || Carp::croak "cannot open file $junctionMatrixFile to read\n";
	open ($fout, ">$newJunctionMatrixFile") || Carp::croak "cannot open file $newJunctionMatrixFile to write\n";

	my $line = <$fin>; chomp $line;
	print $fout $line, "\n";	

	$line = <$fin>;
	chomp $line;
	my ($junctionNumber, $sampleNumber) = split (/\s/, $line);
	my $groupNumber = keys %$groups;	

	print $fout join ("\t", $junctionNumber, $groupNumber), "\n";

	$line = <$fin>; #line with sample names
	chomp $line;

	my @cols = split (/\s+/, $line);
	shift @cols; shift @cols;

	Carp::croak "incorrect number of samples detected\n"	if (@cols != $sampleNumber);
	
	my %samples = map {$cols[$_] => $_} (0..($sampleNumber-1));

	#Carp::croak Dumper (\%samples), "\n";

	my @groupNames = sort {$groups->{$a}->{'id'} <=> $groups->{$b}->{'id'}} keys %$groups;
	#Carp::croak Dumper (\@groupNames), "\n";

	print $fout join ("\t", "Name", "Description", @groupNames), "\n";
	my $jiter = 0;

	print "combining junction counts ...\n" if $verbose;

	while ($line =<$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;

		print "$jiter ...\n" if $verbose && $jiter % 500 == 0;
		$jiter++;

		my @cols = split (/\s+/, $line);
		my $junctionId = shift @cols;
		my $transcriptId = shift @cols;

		Carp::croak "incorrect number of columns detected\n" if $#cols +1 != keys %samples;

		my @counts;
		my $giter = 0;

		foreach my $groupId (@groupNames)
		{
			my $groupSamples = $groups->{$groupId}->{'samples'};
			
			foreach my $s (@$groupSamples)
			{
				my $idx = $samples{$s};
				$counts[$giter] += $cols[$idx];
			}
			$giter++;
		}

		print $fout join ("\t", $junctionId, $transcriptId, @counts), "\n";
	}	

	close ($fin);
	close ($fout);
}


