#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;

use MyConfig;

my $prog = basename ($0);
my $verbose = 0;
my $dispersion = "tagwise"; #"common"
my $MAPlotFile = "";
my $cache = getDefaultCache ($prog);
my $base = "";

#my $method = "mean";  #sum

GetOptions (
	"base:s"=>\$base,
	"disp:s"=>\$dispersion,
	"MAplot:s"=>\$MAPlotFile,
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
detags.diff.common <- rownames(de.common\$table[p.adj <0.05 & abs(log2FC) > 1,])
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

detags.diff.tagwise <- rownames(de.tagwise\$table[p.adj <0.05 & abs(log2FC) > 1,])
pdf (file="$MAPlotFile")
plotSmear(d, pair=c("$groupNames[1]", "$groupNames[0]"), de.tags = detags.diff.tagwise, main = "MA Plot Using Tagwise Dispersion")
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


