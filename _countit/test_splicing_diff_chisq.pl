#!/usr/bin/perl -w

use strict;
use File::Basename;
use Carp;
use Getopt::Long;
use Math::CDF qw(:all); 
use Data::Dumper;


use Common;
use Bed;

my $prog = basename ($0);

my $groupStr = "";
my $outputFile = "";
my $verbose = 0;

my @ARGV0 = @ARGV;


GetOptions (
    'g=s'=>\$groupStr,
    'o=s'=>\$outputFile,
    'v'=> \$verbose);

if (@ARGV < 2)
{
    print "evaluating differential splicing\n";
    print "Usage: $prog [options] <in1.txt> <in2.txt> [in_n.txt]\n";
    print " -g [string]: group definition (e.g., 1,1,1,2,2,2)\n";
    print " -o [file]  : output file\n";
    print " -v         : verbose\n";
    exit (1);
}



my @inBedFiles = @ARGV;

print "CMD=", join(" ", $prog, @ARGV0), "\n" if $verbose;

my @group = split(/\,/, $groupStr);

Carp::croak "the number of input Bed files is inconsistent with group definition\n" if $#group != $#inBedFiles;

my $sampleNum = @group;
my %groupHash;

map {push @{$groupHash{$group[$_]}}, $_} (0 ..($sampleNum-1));

Carp::croak "group label must be 1 or 2\n" if @{$groupHash{1}} + @{$groupHash{2}} != $sampleNum;


print "reading tag numbers for each AS event ...\n" if $verbose;
my @data;
foreach my $inBedFile (@inBedFiles)
{
    print "reading $inBedFile ...\n" if $verbose;
    my $fin;
    my $iter = 0;
    open ($fin, "<$inBedFile") || Carp::croak "cannot open file $inBedFile to read\n";
    while (my $line = <$fin>)
    {
	chomp $line;
	next if $line =~/\#/;
	next if $line =~/^\s*$/;

	print "$iter ...\n" if $verbose && $iter % 5000 == 0;
	my ($chrom, $chromStart, $chromEnd, $name, $score, $strand, $type, $isoforms, $isoform1, $isoform2) = split (/\t/, $line);
	if ($#data >= $iter)
	{
	    my $as = $data[$iter];
	    Carp::croak "in consistency in $iter: $line\n" unless $name eq $as->{'name'};
	}
	else
	{
	    $data[$iter]= {
		chrom=>$chrom, chromStart=>$chromStart, chromEnd=>$chromEnd -1, name=>$name, score=>$score, strand=>$strand, 
		type=>$type, isoforms=>$isoforms};
	}
	push @{$data[$iter]->{'1'}}, $isoform1;
	push @{$data[$iter]->{'2'}}, $isoform2;

	$iter++;
    }
    close ($fin);
}

my $nevents = @data;

print "$sampleNum samples and $nevents events loaded\n" if $verbose;


print "evaluating differential splicing ...\n";

my $fout;
open ($fout, ">$outputFile") || Carp::croak "cannot open file $outputFile to write\n";

print $fout "#", join ("\t", 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 
				'type', 'isoforms', 'sample1/isoform1', 'sample1/isoform2', 'sample2/isoform1', 'sample2/isoform2', 'p'), "\n";

for (my $iter = 0; $iter < @data; $iter++)
{
    print "$iter ...\n" if $verbose && $iter % 5000 == 0;
    my $as = $data[$iter];
		#next unless $as->{'name'} eq 'CA-65945-54758-58020-58077-59557[INC][24/13]';

	my $isoform1 = $as->{'1'};
    my $isoform2 = $as->{'2'};

    my $group1 = $groupHash{'1'};
    my $group2 = $groupHash{'2'};

    #H0: no differential splicing
    my $chisqH0 = chisq ([$isoform1, $isoform2]);

	#print "chisqH0 = $chisqH0\n";

    #H1: there is differential splicing
    my $chisqH1 = 0;

    my @isoform1Sum;
    my @isoform2Sum;
    foreach my $groupId (1 ..2)
    {
		my $idx = $groupHash{$groupId};
		my @isoform1 = @{$isoform1}[@$idx];
		my @isoform2 = @{$isoform2}[@$idx];
	
		$isoform1Sum [$groupId-1] = Common::sum (\@isoform1);
		$isoform2Sum [$groupId-1] = Common::sum (\@isoform2);
	
		$chisqH1 += chisq([\@isoform1, \@isoform2]);

		#print "chisqH1 = $chisqH1\n";
    }
    
    my $p = 1;
    if ($chisqH1 > 1e-6)
    {
		#my $f = $chisqH0 / ($sampleNum -1) / ($chisqH1 / ($sampleNum -2));
		my $f = ($chisqH0 -$chisqH1) / ($chisqH1 / ($sampleNum -2));
		$f = 0 if $f < 0;
		#print "f=$f\n";

		$p = 1 - pf ($f, 1, $sampleNum -2);
		$p = 0 if $p < 0;

		#print "p=$p\n";
    }
    
    print $fout join ("\t", bedToLine ($as), $as->{'type'}, $as->{'isoforms'}, 
			$isoform1Sum[0], $isoform2Sum[0], $isoform1Sum[1], $isoform2Sum[1], $p), "\n";
    
}
close ($fout);

sub chisq
{
    my $dat = $_[0];
 
	#print Dumper ($dat), "\n";		
    my $nrow = @$dat;
    my $ncol = @{$dat->[0]};
    
    my @a; #rowSum / total
    my @b; #colSum / total
    
    for (my $i = 0; $i < $nrow; $i++)    
    {
		$a[$i] = Common::sum($dat->[$i]);
    }
    
    for (my $j = 0; $j < $ncol; $j++)
    {
		my @x = map {$dat->[$_]->[$j]} (0 .. ($nrow -1));
		$b[$j] = Common::sum(\@x);
    }

    my $total = Common::sum (\@a);
    return 0 if $total <= 0;

    @a = map {$_/$total} @a;
    @b = map {$_/$total} @b;

    my $chisq = 0;
    for (my $i = 0; $i < $nrow; $i++)
    {
		for (my $j = 0; $j < $ncol; $j++)
		{
		    my $dat_exp = $a[$i] * $b[$j] * $total;
	   		$chisq += ($dat->[$i][$j] - $dat_exp)**2 / $dat_exp if $dat_exp > 0;
		}
    }
	#print "chisq=$chisq\n";
    return $chisq;
}






