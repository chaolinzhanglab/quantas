#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

use MyConfig;
use Bed;
use Common;


my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $asType = 'cass'; #alt5, alt3
my $verbose = 0;
my $weight = 0;
my $mean = 0;
my $big = 0;
my $separateStrand = 0;

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

GetOptions (
		'big'=>\$big,
		'type:s'=>\$asType,
		'weight'=>\$weight,
		'mean'=>\$mean,
		'ss'=>\$separateStrand,
		'c|cache:s'=>\$cache,
		'keep-cache'=>\$keepCache,
		'v|verbose'=>\$verbose);

if (@ARGV != 3)
{
	print "summarize the number of reads for each isoform\n";
	print "Usage: $prog [options] <as.bed> <tag.bed> <summary.txt>\n";
	print " <as.bed> -- bed file of AS events\n";
	print " <tag.bed> -- bed file of all tags, gz file allowed\n";
	print "OPTIONS:\n";
	print " -big           : the tag file is big\n";
	print " -type [string] : AS type ([cass]|taca|alt5|alt3|mutx|iret|apat|alts|altt)\n";
	print " -weight        : weight tags according to score\n";
	print " -mean          : find means instead of sum (default off)\n";
	print " --ss           : consider the two strands separately\n";
	print " -c             : cache dir ($cache)\n";
	print " --keep-cache   : keep cache files\n";
	print " -v             : verbose\n";
	exit (1);
}


system ("mkdir $cache") unless -d $cache;
my ($ASEventBedFile, $tagBedFile, $summaryFile) = @ARGV;
	
my $tagJunctionBedFile = "$cache/tag.junction.bed";

if ($asType ne 'altt' && $asType ne 'alts' && $asType ne 'apat') #alt start or alt termination
{
	my $cmd = "grep -v \"^track\" $tagBedFile | awk '{if(NF==12 && \$10>1) {print \$0}}' > $tagJunctionBedFile";
	$cmd = "gunzip -c $tagBedFile | grep -v \"^track\" | awk '{if(NF==12 && \$10>1) {print \$0}}' > $tagJunctionBedFile" if $tagBedFile =~/\.gz$/;
	$cmd = "bunzip2 -c $tagBedFile | grep -v \"^track\" | awk '{if(NF==12 && \$10>1) {print \$0}}' > $tagJunctionBedFile" if $tagBedFile =~/\.bz2$/;

	my $ret = system ($cmd);
	print "CMD $cmd failed: $?\n" if $ret != 0;
}

$weight = 1 if $mean;

my $summarizeMethod = $mean ? 'mean' : 'sum';
my $summaryFunc = $mean ? \&mean : \&sum;
my $bigFlag = $big ? '-big' : '';
my $verboseFlag = $verbose ? '-v' : '';
my $ssFlag = $separateStrand ? '--ss' : '';

##########################################
#handle exonic tags in the AS region
##########################################

print "extract alternative exonic regions from $ASEventBedFile ...\n" if $verbose;
my $ASExonBedFile = "$cache/as.exon.bed"; 

if ($asType eq 'cass' || $asType eq 'mutx' || $asType eq 'taca') 
{
	my $cmd = "perl $cmdDir/gene2ExonIntron.pl -v -internal -nid -oe $ASExonBedFile $ASEventBedFile";
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	
	Carp::croak	"CMD crashed: $cmd, $?\n" unless $ret == 0;
}
elsif ($asType eq 'alt5' || $asType eq 'alt3')
{
	my $twoExonBedFile = "$cache/twoexon.bed";
	my $cmd = "perl $cmdDir/gene2ExonIntron.pl -v -oe $twoExonBedFile $ASEventBedFile";
	print $cmd, "\n" if $verbose;

	my $ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;
	
	my $exonIdx = 0;
	$exonIdx = 1 if $asType eq 'alt3';
	my $oneExonBedFile = "$cache/oneexon.bed";

	#$cmd = "grep -P \"_$exonIdx\\t\" $twoExonBedFile > $oneExonBedFile";
	$cmd = "grep \"_$exonIdx\" $twoExonBedFile > $oneExonBedFile";
	$ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

	my $constBedFile = "$cache/const.exon.bed";
	$cmd = "perl $cmdDir/getASRegion.pl -if $asType -enum -v $oneExonBedFile $constBedFile $ASExonBedFile";
	$ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

	unlink $twoExonBedFile, $oneExonBedFile, $constBedFile;
}
elsif ($asType eq 'iret')
{
	#get the retained intron
	my $cmd = "perl $cmdDir/gene2ExonIntron.pl -v -internal -nid -oi $ASExonBedFile.tmp $ASEventBedFile";
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

	#$cmd = "perl ~/scripts/bedExt.pl -n up -l \"-5\" -r 5 $ASExonBedFile $ASExonBedFile.5SS";
	$cmd = "perl $cmdDir/bedExt.pl -n up -l \"5\" -r 5 $ASExonBedFile.tmp - | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\":5ss\\t\"\$5\"\\t\"\$6}' > $ASExonBedFile";	#2012-06-12
	print $cmd, "\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;
	
	#$cmd = "perl ~/scripts/bedExt.pl -n down -l \"-5\" -r 5 $ASExonBedFile $ASExonBedFile.3SS";
	$cmd = "perl $cmdDir/bedExt.pl -n down -l \"-5\" -r \"-5\" $ASExonBedFile.tmp - |awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\":3ss\\t\"\$5\"\\t\"\$6}' >> $ASExonBedFile";  #2012-06-12
	
	print $cmd, "\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

	system ("rm -rf $ASExonBedFile.tmp");
	
	#changes in 2012-06-12
	#we count the tags that overlaps with +5 relative to the 5'ss and -5 relative to the 3'ss as evidence of retained introns
	#this is to avoid mapping errors in which a junction read is mapped as an exon body read, when the overlap with one side is very short
	
	#11/20/2014, now we output reads overlap with 5'/3' end of the intron separately
}
elsif ($asType eq 'apat')
{
	my $cmd = "cp $ASEventBedFile $ASExonBedFile";
	print $cmd, "\n" if $verbose;
    my $ret = system ($cmd);
    Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;
}

elsif ($asType eq 'alts' || $asType eq 'altt')
{
	my $cmd = "perl $cmdDir/bedExt.pl -n up -l \"5\" -r 5 $ASEventBedFile $ASExonBedFile";
	$cmd = "perl $cmdDir/bedExt.pl -n down -l \"-5\" -r \"-5\" $ASEventBedFile $ASExonBedFile" if $asType eq 'altt';
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;
}


print "count number of tags on each AS exon ...\n" if $verbose;
my $ASExonTagCountFile = "$cache/as.exon.tagcount.bed";
my $weightFlag = $weight ? '-weight' : '';
$weightFlag = '-weight-avg' if $mean;

my $cmd = "perl $cmdDir/tag2profile.pl $verboseFlag $bigFlag $weightFlag $ssFlag -c $cache/tmp_exon_tag_count -region $ASExonBedFile $tagBedFile $ASExonTagCountFile";
print $cmd, "\n" if $verbose;

my $ret = system ($cmd);
Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

if ($asType eq 'taca')
{
	#get the total or average of multiple exons
	print "sum up tags on AS exons ...\n" if $verbose;
	my $cmd = "perl $cmdDir/uniqRow.pl -c $summarizeMethod -id 3 -value 4 -v $ASExonTagCountFile $ASExonTagCountFile"; 
	my $ret = system ($cmd);
    Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;	
}


######################################################
#handle junction tags
######################################################

my $ASJunctionTagCountFile = "$cache/as.junction.tagcount.bed";
if ($asType ne 'alts' && $asType ne 'altt' && $asType ne 'apat')
{

	print "extract introns from $ASEventBedFile ...\n" if $verbose;
	my $ASIntronBedFile = "$cache/as.intron.bed";
	$cmd = "perl $cmdDir/gene2ExonIntron.pl -v -oi $ASIntronBedFile $ASEventBedFile";
	print $cmd, "\n" if $verbose;

	$ret = system ($cmd);
	Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;

	#deal with mean
	my $meanFlag = $mean ? "-mean" : "";
	
	print "count junction tags for each intron ...\n" if $verbose;
	$cmd = "perl $cmdDir/summarize_junction.pl $verboseFlag $ssFlag $weightFlag $meanFlag $ASIntronBedFile $tagJunctionBedFile $ASJunctionTagCountFile";
	print $cmd, "\n" if $verbose;
	$ret = system ($cmd);
    Carp::croak "CMD crashed: $cmd, $?\n" unless $ret == 0;
}


########################################################################
#put counts of exonic tags and junctoin tags together into a hash table
########################################################################

my %tagCountHash;
print "reading exon tag count from $ASExonTagCountFile ...\n" if $verbose;

my $exonTagCount = readBedFile ($ASExonTagCountFile, $verbose);

my $n = @$exonTagCount;

print "$n entries loaded\n" if $verbose;

foreach my $e (@$exonTagCount)
{
	my $name = $e->{"name"};
	my $count = $e->{"score"};

	my ($asId, $isoformId, $evi) = ("", "", "");

	if ($name =~/^(.*?)\[(.*?)\]\[(.*?)\]/)
	{
		($asId, $isoformId, $evi) = ($1, $2, $3);
		$isoformId = "SPLICE" if $asId eq 'iret';
	}
	elsif ($asType eq 'iret')
	{
		#11/20/2014, cz
		#make an exception for intron retention to allow more flexible naming
		$asId = $name;
		$asId =~s/:[5|3]ss$//g;

		$isoformId = "SPLICE";
		$evi = 0;
	}
	else
	{
		Carp::croak "incorrect in format: $name\n";
	}

	if ($asType eq 'cass' || $asType eq 'taca')
	{
		my $isoformIdPair = "INC/SKIP";
		$tagCountHash{$asId}->{$isoformIdPair}->{"exon"} = $count;
	}
	elsif ($asType eq 'alt5' || $asType eq 'alt3')
	{
		my $isoformIdPair = $isoformId;
		$tagCountHash{$asId}->{$isoformIdPair}->{"exon"} = $count;
		$tagCountHash{$asId}->{$isoformIdPair}->{"space"} = $e->{'chromEnd'} - $e->{'chromStart'} + 1; #length of the alternative region
	}
	elsif ($asType eq 'mutx')
	{
		$isoformId=~/(\d+)$/;
		my $isoformIdx = $1;
		my @evis = split (/\//, $evi);

		#idpair always from 5' / 3'
		for (my $i = 0; $i < $isoformIdx; $i++)
		{
			#some of the pairs might not actually exists when M$i is not a canonical isoform
			#these pairs will be eliminated when output
			my $isoformIdPair = "M$i/M$isoformIdx";
			$tagCountHash{$asId}->{$isoformIdPair}->{"exon"}->{$isoformId} = $count;
		}

		for (my $i = $isoformIdx+1; $i < @evis; $i++)
		{
			#some of the pairs might not actually exists when M$i is not a canonical isoform
			#these pairs will be eliminated when output
			my $isoformIdPair = "M$isoformIdx/M$i";
			$tagCountHash{$asId}->{$isoformIdPair}->{"exon"}->{$isoformId} = $count;
		}
	}
	elsif ($asType eq 'iret')
	{
		my $isoformIdPair = "SPLICE/RET";
		
		if ($name =~/:5ss$/)
		{
			$tagCountHash{$asId}->{$isoformIdPair}->{"exon"}->{'5ss'} = $count;
		}
		elsif ($name =~/:3ss$/)
		{
			$tagCountHash{$asId}->{$isoformIdPair}->{"exon"}->{'3ss'} = $count;
		}
		else
		{
			Carp::croak "incorrect format: $name\n";
		}
	}
	elsif ($asType eq 'apat' || $asType eq 'alts' || $asType eq 'altt')
	{
		$tagCountHash{$asId}->{$isoformId}->{"exon"} = $count;
		$tagCountHash{$asId}->{$isoformId}->{'size'} = $e->{'chromEnd'} - $e->{'chromStart'} + 1 if $asType eq 'apat';
	}

	#$name = join ("", $asId, "[", $isoformId, "][", $evi, "]");
	#$tagCountHash{$asId}->{$name} = $count;
	
	#$tagCountHash{$asId}->{$isoformIdPair}->{"exon"} = $count;
}


#pair apat, alternative start or termination
if ($asType eq 'apat' || $asType eq 'alts' || $asType eq 'altt')
{
	foreach my $asId (keys %tagCountHash)
	{
		my $isoforms = $tagCountHash{$asId};
		my %isoformPairs;
		for (my $i = 0; $i < keys %$isoforms; $i++)
		{	
			next unless exists $isoforms->{"A$i"};
			for (my $j = $i+1; $j < keys %$isoforms; $j++)
			{
				next unless exists $isoforms->{"A$j"};
				$isoformPairs{"A$i/A$j"}->{"exon"}->{"A$i"} = $isoforms->{"A$i"}->{"exon"};
				$isoformPairs{"A$i/A$j"}->{"exon"}->{"A$j"} = $isoforms->{"A$j"}->{"exon"};
				
				if ($asType eq 'apat')
				{
					$isoformPairs{"A$i/A$j"}->{"size"}->{"A$i"} = $isoforms->{"A$i"}->{"size"};
					$isoformPairs{"A$i/A$j"}->{"size"}->{"A$j"} = $isoforms->{"A$j"}->{"size"};
				}
			}
		}
		$tagCountHash{$asId} = \%isoformPairs;
	}
}

#

#junction tag counts
if ($asType ne 'apat' && $asType ne 'alts' && $asType ne 'altt')
{
	my $junctionTagCount = readBedFile ($ASJunctionTagCountFile, $verbose);
	
	$n = @$junctionTagCount;
	
	print "$n entries loaded ...\n" if $verbose;
	
	foreach my $j (@$junctionTagCount)
	{
		my $name = $j->{"name"};
		my $count = $j->{"score"};

		#$name =~/^(.*?)\[(.*?)\]\[(.*?)\]\_(\d+)$/;
		
		my ($asId, $isoformId, $evi, $junctionId);
	
		if ($name =~/^(.*?)\[(.*?)\]\[(.*?)\].*?\_(\d+)$/)
		{
			#to accomodate the pattern of new mutx ids
			($asId, $isoformId, $evi, $junctionId) = ($1, $2, $3, $4);
			$isoformId = "SPLICE" if $asId eq 'iret';
		}
		elsif ($asType eq 'iret')
		{
			$asId = $name;
			$asId =~s/_\d+$//g;

			$isoformId = "SPLICE";
			$evi = 0;
			$junctionId = 0;
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
				$tagCountHash{$asId}->{"INC/SKIP"}->{"junction"}->{"INC"}->{$junctionId} = $count;
			}
			else
			{
				$tagCountHash{$asId}->{"INC/SKIP"}->{"junction"}->{"SKIP"} = $count;
			}
		}
		elsif ($asType eq 'taca')
		{
			if ($isoformId eq 'INC')
			{
				$tagCountHash{$asId}->{"INC/SKIP"}->{"junction"}->{"INC"} += $count;
				#report the tag number of individual junctions
				if (exists $tagCountHash{$asId}->{"INC/SKIP"}->{"junction"}->{"INC2"})
				{
					$tagCountHash{$asId}->{"INC/SKIP"}->{"junction"}->{"INC2"} .= "|$count";
				}
				else
				{
					$tagCountHash{$asId}->{"INC/SKIP"}->{"junction"}->{"INC2"} = $count;
				}
			}
			else
			{
				$tagCountHash{$asId}->{"INC/SKIP"}->{"junction"}->{"SKIP"} += $count;
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
				$tagCountHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId} = $count;
			}
			for (my $i = $isoformIdx + 1; $i < @evis; $i++)
			{
				my $isoformIdPair = "A$isoformIdx/A$i";
				$tagCountHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId} = $count;
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
				$tagCountHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId}->{$junctionId} = $count;
			}
	
			for (my $i = $isoformIdx + 1; $i < @evis; $i++)
			{
				my $isoformIdPair = "M$isoformIdx/M$i";
				$tagCountHash{$asId}->{$isoformIdPair}->{"junction"}->{$isoformId}->{$junctionId} = $count;
			}
		}
		elsif ($asType eq 'iret')
		{
			$tagCountHash{$asId}->{"SPLICE/RET"}->{"junction"} = $count;
		}
		#$tagCountHash{$asId}->{$name} = $count;
	}
}


############################################################
#output
############################################################


print "output summary file $summaryFile ...\n" if $verbose;

my $fout;
open ($fout, ">$summaryFile") || Carp::croak "cannot open file $summaryFile to write\n";

my @fixedColumnHeader = ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "isoformIDs", "isoform1Tags", "isoform2Tags");

if ($asType eq 'cass')
{
	print $fout "#", join ("\t", @fixedColumnHeader,
			"exonTags", "inclusionJunction1Tags", "inclusionJunction2Tags", "skippingJunctionTags"), "\n";
}
elsif ($asType eq 'taca')
{
	print $fout "#", join ("\t", @fixedColumnHeader,
			"exonTags", "inclusionJunctionTags", "skippingJunctionTags", "individualInclusionJunctionTags"), "\n";
}
elsif ($asType eq 'alt5' || $asType eq 'alt3')
{
	print $fout "#", join ("\t", @fixedColumnHeader, 
			"altSSDistance", "exonTags", "proximalJunctionTags", "distalJunctionTags"), "\n";
}
elsif ($asType eq 'mutx')
{
	print $fout "#", join ("\t", @fixedColumnHeader,
			"5'ExonTags", "5'ExonJunction1Tags", "5'ExonJunction2Tags", "3'ExonTags", "3'ExonJunction1Tags", "3'ExonJunction2Tags"), "\n";
}
elsif ($asType eq 'iret')
{
	print $fout "#", join ("\t", @fixedColumnHeader,
			"retainedIntron5SSTags", "retainedIntron3SSTags", "junctionTags"), "\n";
}
elsif ($asType eq 'apat')
{
	print $fout "#", join ("\t", @fixedColumnHeader,
			"commonSize", "altSize"), "\n";
}
elsif ($asType eq 'alts' || $asType eq 'altt')
{
	print $fout "#", join ("\t", @fixedColumnHeader), "\n";
}

my $ASEvents = readBedFile ($ASEventBedFile, $verbose);

$n = @$ASEvents;
print "$n entries loaded\n" if $verbose;

my %ASEoutput;

foreach my $e (@$ASEvents)
{
	my $name = $e->{"name"};

	my ($asId, $isoformId) = ("", "");
	
	if ($name =~/^(.*?)\[(.*?)\]/)
	{
		$asId = $1;
		$isoformId = $2;
		$isoformId = "SPLICE" if $asId eq 'iret';
	}
	elsif ($asType eq 'iret')
	{
		#11/20/2014, cz
		#make an exception for iret, since there is only one spliced isoform in the id

		$asId = $name;
		#$asId =~s/:[5|3]ss$//g;

		$isoformId = "SPLICE";
	}
	else
	{
		Carp::croak "incorrect name format: $name\n";
	}


	#$asId =~/^(\w\w)/;
	#my $asType = $1;

	next if (($asType eq 'cass' || $asType eq 'taca') && $isoformId eq 'SKIP');

	next if exists $ASEoutput{$asId};	#already dumpped

	my $ASEs = $tagCountHash{$asId};

	foreach my $isoformIdPair (sort keys %$ASEs)
	{
		my $ase = $ASEs->{$isoformIdPair};
		#next if exists $ASEoutput{$asId . "[" . $isoformIdPair . "]"};

		if ($asType eq 'cass')
		{
			my $exonTagCount = exists $ase->{"exon"} ? $ase->{"exon"} : 0;
			my $inc1JunctionTagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{"INC"} && exists $ase->{"junction"}->{"INC"}->{0} ? $ase->{"junction"}->{"INC"}->{0} : 0;
			my $inc2JunctionTagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{"INC"} && exists $ase->{"junction"}->{"INC"}->{1} ? $ase->{"junction"}->{"INC"}->{1} : 0;
			my $skipJunctionTagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{"SKIP"} ? $ase->{"junction"}->{"SKIP"} : 0;
	

			print $fout join ("\t", $e->{"chrom"}, $e->{"chromStart"}, $e->{"chromEnd"}, $e->{"name"}, $e->{"score"}, $e->{"strand"}, $asType, "INC/SKIP",
				#$summaryFunc->([$exonTagCount, $inc1JunctionTagCount, $inc2JunctionTagCount]), $skipJunctionTagCount,
				$exonTagCount, $skipJunctionTagCount, #07/22/2012, all tags are intesected with exons now, so we do not need to add junction tags explicitly
				$exonTagCount, $inc1JunctionTagCount, $inc2JunctionTagCount, $skipJunctionTagCount), "\n";
		}
		elsif ($asType eq 'taca')
		{
			my $exonTagCount = exists $ase->{"exon"} ? $ase->{"exon"} : 0;
			my $incJunctionTagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{"INC"} ? $ase->{"junction"}->{"INC"} : 0;
			my $incJunction2TagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{"INC2"} ? $ase->{"junction"}->{"INC2"} : "NA";

			my $skipJunctionTagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{"SKIP"} ? $ase->{"junction"}->{"SKIP"} : 0;
			print $fout join ("\t", $e->{"chrom"}, $e->{"chromStart"}, $e->{"chromEnd"}, $e->{"name"}, $e->{"score"}, $e->{"strand"}, $asType, "INC/SKIP",
				#$summaryFunc->([$exonTagCount, $incJunctionTagCount]), $skipJunctionTagCount,
				$exonTagCount, $skipJunctionTagCount, 
				$exonTagCount, $incJunctionTagCount, $skipJunctionTagCount, $incJunction2TagCount), "\n";
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

			my $exonTagCount = exists $ase->{"exon"} ? $ase->{"exon"} : 0;
			my $proximalJunctionTagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{$proximalIsoformId} ? $ase->{"junction"}->{$proximalIsoformId} : 0;
			my $distalJunctionTagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{$distalIsoformId} ? $ase->{"junction"}->{$distalIsoformId} : 0;
			my $space = exists $ase->{"space"} ? $ase->{"space"} : 0; #length of the alternative region
			
			print $fout join ("\t", $e->{"chrom"}, $e->{"chromStart"}, $e->{"chromEnd"}, $e->{"name"}, $e->{"score"}, $e->{"strand"}, $asType, "$proximalIsoformId/$distalIsoformId", 
					#$summaryFunc->([$exonTagCount, $proximalJunctionTagCount]), $distalJunctionTagCount,
					$exonTagCount, $distalJunctionTagCount,
					$space, $exonTagCount, $proximalJunctionTagCount, $distalJunctionTagCount), "\n";
		}
		elsif ($asType eq 'mutx')
		{
			my ($isoformId1, $isoformId2) = split (/\//, $isoformIdPair);
			next unless exists $ase->{"exon"}->{$isoformId1} && exists $ase->{"exon"}->{$isoformId2};
			#when a pair does not actually exist

			$isoformId1=~/(\d+)$/;
			my $isoformIdx1 = $1;

			$isoformId2=~/(\d+)$/;
			my $isoformIdx2 = $1;

			my $leftIsoformId = ($isoformIdx1 < $isoformIdx2) ? $isoformId1 : $isoformId2;
			my $rightIsoformId = ($isoformIdx1 > $isoformIdx2) ? $isoformId1 : $isoformId2;

			my $leftExonTagCount = exists $ase->{"exon"} && exists $ase->{"exon"}->{$leftIsoformId} ? $ase->{"exon"}->{$leftIsoformId} : 0;
			my $rightExonTagCount = exists $ase->{"exon"} && exists $ase->{"exon"}->{$rightIsoformId} ? $ase->{"exon"}->{$rightIsoformId} : 0;

			my $leftExonJunction1TagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{$leftIsoformId} && exists $ase->{"junction"}->{$leftIsoformId}->{0} ? $ase->{"junction"}->{$leftIsoformId}->{0} : 0;
			my $leftExonJunction2TagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{$leftIsoformId} && exists $ase->{"junction"}->{$leftIsoformId}->{1} ? $ase->{"junction"}->{$leftIsoformId}->{1} : 0;
			
			my $rightExonJunction1TagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{$rightIsoformId} && exists $ase->{"junction"}->{$rightIsoformId}->{0} ? $ase->{"junction"}->{$rightIsoformId}->{0} : 0;
			my $rightExonJunction2TagCount = exists $ase->{"junction"} && exists $ase->{"junction"}->{$rightIsoformId} && exists $ase->{"junction"}->{$rightIsoformId}->{1} ? $ase->{"junction"}->{$rightIsoformId}->{1} : 0;
			
			print $fout join ("\t", $e->{"chrom"}, $e->{"chromStart"}, $e->{"chromEnd"}, $e->{"name"}, $e->{"score"}, $e->{"strand"}, $asType, "$leftIsoformId/$rightIsoformId", 
				#$summaryFunc->([$leftExonTagCount, $leftExonJunction1TagCount, $leftExonJunction2TagCount]), $summaryFunc->([$rightExonTagCount, $rightExonJunction1TagCount, $rightExonJunction2TagCount]),	
				$leftExonTagCount, $rightExonTagCount,
				$leftExonTagCount, $leftExonJunction1TagCount, $leftExonJunction2TagCount, 
				$rightExonTagCount, $rightExonJunction1TagCount, $rightExonJunction2TagCount), "\n";
		}
		elsif ($asType eq 'iret')
		{
			my $intronTagCount5SS = exists $ase->{"exon"} && $ase->{'exon'}{'5ss'} ? $ase->{"exon"}{'5ss'} : 0;
			my $intronTagCount3SS = exists $ase->{"exon"} && $ase->{'exon'}{'3ss'} ? $ase->{"exon"}{'3ss'} : 0;
			my $intronTagCount = $intronTagCount5SS + $intronTagCount3SS;

			my $junctionTagCount = exists $ase->{"junction"} ? $ase->{"junction"} : 0;
			print $fout join ("\t", $e->{"chrom"}, $e->{"chromStart"}, $e->{"chromEnd"}, $e->{"name"}, $e->{"score"}, $e->{"strand"}, $asType, "RET/SPLICE",
				$intronTagCount, $junctionTagCount, 
				$intronTagCount5SS, $intronTagCount3SS, $junctionTagCount), "\n";
		}
		elsif ($asType eq 'apat')
		{
			my ($isoformId1, $isoformId2) = split (/\//, $isoformIdPair);
			print $fout join ("\t", $e->{"chrom"}, $e->{"chromStart"}, $e->{"chromEnd"}, $e->{"name"}, $e->{"score"}, $e->{"strand"}, $asType, $isoformIdPair,
				$ase->{"exon"}->{$isoformId1}, $ase->{"exon"}->{$isoformId2} - $ase->{"exon"}->{$isoformId1}, $ase->{"size"}->{$isoformId1}, $ase->{"size"}->{$isoformId2} - $ase->{"size"}->{$isoformId1}), "\n";
		}
		elsif ($asType eq 'alts' || $asType eq 'altt')
		{
			my ($isoformId1, $isoformId2) = split (/\//, $isoformIdPair);
			print $fout join ("\t", $e->{"chrom"}, $e->{"chromStart"}, $e->{"chromEnd"}, $e->{"name"}, $e->{"score"}, $e->{"strand"}, $asType, $isoformIdPair,
				$ase->{"exon"}->{$isoformId1}, $ase->{"exon"}->{$isoformId2}), "\n";
		}
	}
	$ASEoutput{$asId} = 1;
}

close ($fout);

system ("rm -rf $cache") unless $keepCache;

