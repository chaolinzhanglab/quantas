#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;

use Scalar::Util qw(looks_like_number);
use Statistics::Basic qw(:all);

use MyConfig;

my $prog = basename ($0);
my $verbose = 0;
my $dispersion = "tagwise"; #"common"
my $MAPlotFile = "";
my $MAPlotFDR = 0.05;
my $MAPlotLogFC=1;

my $cache = getDefaultCache ($prog);
my $printRPKM = 0;
my $base = "";

#my $method = "mean";  #sum

GetOptions (
	"base:s"=>\$base,
	"disp:s"=>\$dispersion,
	"MAplot:s"=>\$MAPlotFile,
	"FDR:f"=>\$MAPlotFDR,
	"logFC:f"=>\$MAPlotLogFC,
	"rpkm"=>\$printRPKM,
	"c:s"=>\$cache,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "test differential expression\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " -base         [string] : base dir of input data\n";
	print " -disp         [string] : ([tagwise]|common|common:0.1)\n";
	print " -MAplot       [file]   : file name to write MAplot\n";
	print " -FDR          [float]  : FDR to highlight genes in MAplot ($MAPlotFDR)\n";
	print " -logFC        [float]  : log2FC to highlight genes in MAplot($MAPlotLogFC)\n";
	print " -rpkm                  : include average group RPKM values\n";
	print " -c            [dir]    : cache dir\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outFile) = @ARGV;


my $commonDisp = 0;

if ($dispersion =~/^common:(.*?)$/)
{
	$commonDisp = $1;
	$dispersion = "common";
}


Carp::croak "$cache already exists\n" if -d $cache || -f $cache;
my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir at $cache\n" unless $ret == 0;


print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readConfigFile ($configFile, $base);

print "done.\n" if $verbose;



my $configFile2 = "$cache/sample.conf";

my $fout;
open ($fout, ">$configFile2") || Carp::croak "cannot open file to write\n";

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
		
		print $fout join ("\t", $inputFile, $gName), "\n";		
		$iter++;
	}
}
close ($fout);
print "$iter samples in total.\n" if $verbose;




print "generating R scripts ...\n" if $verbose;


my $scriptFile = "$cache/script.R";
open ($fout, ">$scriptFile") || Carp::croak "cannot open file $scriptFile to write\n";

if ($dispersion eq 'common')
{
	print $fout <<EOF;

library (edgeR);

conf <- read.table ("$configFile2", sep="\\t", header=F);
files <- as.character (conf[,1]);
groups <- as.character (conf[,2]);

f <- files[1];

d <- read.table (f, sep="\\t", header=F);
gene.names <- d[,2];

d <- readDGE (files, columns=c(1,3), group=groups);

d\$common.dispersion <- $commonDisp;

if (d\$common.dispersion <= 0)
{
	if (min(table(groups))>=2)
	{
		d <- estimateCommonDisp(d);
	}
	else
	{
		d <- estimateGLMCommonDisp (d, method="deviance", robust=TRUE, subset = NULL);
	}
}

de.common <- exactTest(d, pair=c("$groupNames[1]", "$groupNames[0]"), );
p.adj <- p.adjust(de.common\$table\$PValue, method = "BH");
log2FC <- de.common\$table\$logFC;

EOF

if ($MAPlotFile)
{

print $fout <<EOF;
detags.diff.common <- rownames(de.common\$table[p.adj <$MAPlotFDR & abs(log2FC) > $MAPlotLogFC,])
pdf (file="$MAPlotFile")
plotSmear(de.common, pair=c("$groupNames[1]", "$groupNames[0]"), de.tags = detags.diff.common, main = "MA Plot Using Common Dispersion")
#abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)
dev.off()
EOF
}

	print $fout <<EOF;
gene.ids <- rownames(de.common\$table)

out <-cbind (gene.ids, gene.names, de.common\$table, p.adj);

write.table (out, "$outFile", sep="\\t", col.names=T, row.names=F, quote=F);

EOF

}
elsif ($dispersion eq 'tagwise')
{

	print $fout <<EOF;

library (edgeR);

conf <- read.table ("$configFile2", sep="\\t", header=F);
files <- as.character (conf[,1]);
groups <- as.character (conf[,2]);

f <- files[1];

d <- read.table (f, sep="\\t", header=F);
gene.names <- d[,2];

d <- readDGE (files, columns=c(1,3), group=groups);

d <- estimateCommonDisp(d)
#d <- estimateTagwiseDisp(d, prior.n = 10, grid.length = 500)
d <- estimateTagwiseDisp(d, prior.df=10, grid.length=500)


de.tagwise <- exactTest(d, pair=c("$groupNames[1]", "$groupNames[0]"))
p.adj <- p.adjust(de.tagwise\$table\$PValue, method = "BH");
log2FC <- de.tagwise\$table\$logFC;

EOF

if ($MAPlotFile)
{
	print $fout <<EOF;

detags.diff.tagwise <- rownames(de.tagwise\$table[p.adj <$MAPlotFDR & abs(log2FC) > $MAPlotLogFC,])
pdf (file="$MAPlotFile")
plotSmear(de.tagwise, pair=c("$groupNames[1]", "$groupNames[0]"), de.tags = detags.diff.tagwise, main = "MA Plot Using Tagwise Dispersion")
#abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)
dev.off()
EOF
}

	print	$fout <<EOF;

gene.ids <- rownames(de.tagwise\$table)

out <-cbind (gene.ids, gene.names, de.tagwise\$table, p.adj);

write.table (out, "$outFile", sep="\\t", col.names=T, row.names=F, quote=F);

EOF


}
else
{
	Carp::croak "wrong dispersion type: $dispersion\n";
}


close ($fout);

my $cmd = "R --no-save < $scriptFile";

$ret = system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless $ret == 0;

#system ("rm -rf $cache");

if ($printRPKM)
{	
	print "Calculating average RPKM...\n";

	# (my $outFileRPKM = $outFile) =~ s/.txt$/.rpkm.txt/;
	# $outFileRPKM .= ".rpkm" if ($outFileRPKM !~ m/.rpkm.txt$/);

	my $outFileRPKM = "$cache/$outFile.rpkm";

	my %RPKMs;

	foreach my $group (@groupNames)
	{
		my $samples = $groups->{$group}->{"samples"};
		foreach my $sample (@$samples)
		{
			my $inputFile = $sample;
			$inputFile = "$base/$inputFile" if $base ne '';
			open (INFILE, $inputFile) || Carp::croak "cannot open file $inputFile to read\n";
			my $totalTagNum = 0;
			while (my $line = <INFILE>)
			{
				chomp $line;
				next if $line=~/^\s*$/;
				next if $line=~/^\#/;
				my ($gene_id, $gene_symbol, $tag_num, $exon_length, $rpkm) = split(/\t/, $line);
				# $RPKMs{$gene_id}{$sample} = $rpkm if ($rpkm =~ m/^\d+\.?\d*[Ee]?-?\d*$/);
				#$RPKMs{$gene_id}{$sample} = $rpkm if (looks_like_number($rpkm));
				$RPKMs{$gene_id}{$sample} = [$tag_num, $exon_length];
				$totalTagNum += $tag_num;
			}
			close (INFILE);

			foreach my $gene_id (keys %RPKMs)
    		{
        		my $tagNum = $RPKMs{$gene_id}{$sample}[0];
				my $exonLen = $RPKMs{$gene_id}{$sample}[1];

        		$tagNum = 1 if $tagNum == 0;
        		$RPKMs{$gene_id}{$sample} = $tagNum * 1e9 / $exonLen / $totalTagNum;
    		}
		}

	}

	my %RPKMMeans;

	my %RPKMFilter; #whether the max of the groups is above the median
	foreach my $gene (keys %RPKMs)
	{
		$RPKMFilter{$gene} = 0; #note that we assume RPKM cannot be negative
		foreach my $group (@groupNames)
		{
			my $sum = 0;
			
			my $samples = $groups->{$group}->{"samples"};
			foreach my $sample (@$samples)
			{
				if (exists($RPKMs{$gene}{$sample}))
				{
					Carp::croak "negative RPKM, gene=$gene, sample=$sample, RPKM=", $RPKMs{$gene}{$sample}, "\n" if $RPKMs{$gene}{$sample} < 0;
					$sum += $RPKMs{$gene}{$sample};
				}
				else
				{
					Carp::carp "Missing value for gene $gene in sample $sample\n";
				}				
			}
			
			$RPKMMeans{$gene}{$group} = $sum / (scalar @$samples);
			$RPKMFilter{$gene} = $RPKMMeans{$gene}{$group} if $RPKMFilter{$gene} < $RPKMMeans{$gene}{$group};
		}
	}

	my $med = median (values %RPKMFilter);
	map {$RPKMFilter{$_} = $RPKMFilter{$_} > $med ? 1 : 0} keys %RPKMFilter;	


	open (INFILE, $outFile) || Carp::croak "cannot open $outFile to read\n";
	open (OUTFILE, ">$outFileRPKM") || Carp::croak "cannot open $outFileRPKM to write";

	my $count = -1;
	while (my $line = <INFILE>)
	{
		$count++;
		chomp $line;

		print OUTFILE $line;

		# header
		if ($count == 0)
		{
			foreach my $group (@groupNames)
			{
				print OUTFILE "\tRPKM($group)";
			}
			print OUTFILE "\tRPKMFilter";
		}
		else
		{
			my @splits = split(/\t/, $line);
			my $gene = $splits[0];

			foreach my $group (@groupNames)
			{
				print OUTFILE "\t$RPKMMeans{$gene}{$group}";
			}
			print OUTFILE "\t", $RPKMFilter{$gene};
		}

		print OUTFILE "\n";
	}
	close INFILE;
	close OUTFILE;

    system(("mv", $outFileRPKM, $outFile)) == 0 or die "Moving temp RPKM file failed: $?"
}

print "Done.\n";


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
	}
	close ($fin);
	return \%groups;
}


