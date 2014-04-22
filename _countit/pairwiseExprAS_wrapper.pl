use strict;

use Carp;
use Getopt::Long;
use File::Basename;
use List::Util qw(min);
use MyConfig;

# use Parallel::ForkManager;

my $start_run = time();

my $prog = basename ($0);
#my $cmdDir = "~/czlab_src";  
#for codes to be released to the public, directories cannot be hard coded
my $cmdDir = dirname ($0);

my $cache = getDefaultCache ($prog);
#we need to allow users to specify the cache dir, see below

my $type = "expr";
my $goConf = "";
my $base = "";
my $outDir = "";
my $id2gene2symbolFile = "";
#allow users to generate output to specified folder

my $geneCol = -1;
my $pValue = 0.01;
# my $pValueCol = -1;
my $fc = 1;
# my $fcCol = -1;
# my $exprCol = "";
# my @exprCols = ();
my $coverage = 20;
# my $coverageCol = -1;
my $dI = 0.2;
# my $dICol = -1;
my $header = 0;

my $skipFail = 0;
my $negative = 0;
my $medianOut = "";
my $updnOut = "";

my $verbose = 0;

my @medians = ();

GetOptions (
		'type:s'=>\$type,
		'go:s'=>\$goConf,
		'base:s'=>\$base,
		'out-dir:s'=>\$outDir,
		'id2gene2symbol=s'=>\$id2gene2symbolFile,
		# 'e|event-col=i'=>\$geneCol,
		'p|p-value=f'=>\$pValue,
		# 'q|p-value-col=i'=>\$pValueCol,
		'f|FC=f'=>\$fc,
		# 'g|FC-col=i'=>\$fcCol,
		# 'r|rpkm-cols=s'=>\$exprCol,
		'c|coverage=f'=>\$coverage,
		# 'd|coverage-col=i'=>\$coverageCol,
		'i|inclusion=f'=>\$dI,
		# 'j|inclusion-col=i'=>\$dICol,
		# 'h|header-lines=i'=>\$header,
		# 's|skip-fail'=>\$skipFail,
		# 'n|negative'=>\$negative,
		# 'm|median-file=s'=>\$medianOut,
		# 'u|updn-file=s'=>\$updnOut,
		'cache:s'=>\$cache,
		'v|verbose'=>\$verbose
) or die "Error: command line arguments\n";

if (@ARGV != 1)
{
	print "NEED TO ADD SCRIPTS DIR OPTION\n";
	print "Perform pairwise differential expression or differential splicing analysis from countit data.\n";
	print "\n";
	print "USAGE: $prog [options] <input.conf>\n";
	print "\n";
	print "FILES:\n";
	print " input.conf  [string] : Col1 = expr/AS file; Col2 = sample name; tab-delimted\n";
	# print " GO.conf     : first column=GO file, second column=short abbreviation to ID file (e.g. 'cc')\n";
	# print " GO file     : first column=event ID, second column=GO term\n";
	print "\n";
	print "OPTIONS:\n";
	print " -type       [string] : experiment type ([expr]|as)\n";
	# print " -go         : file containing list of GO-term files and abbreviations\n";
	print " -base       [string] : base directory for files in input-conf\n";
	print " --out-dir   [string] : output dir (default=current directory)\n";
	print "\n";
	print "Filtering - Expression and splicing\n";
	print " -p [float]  : p-value (or FDR) below this value [0.01]\n";
	print "\n";
	print "Filtering - Expression-specific\n";
	print " -f [float]  : fold-change / log fold-change above this value [1]\n";
	print "\n";
	print "Filtering - Splicing-specific\n";
	print " --id2gene2symbol [string]: mapping of AS id to gene id and symbol\n";
	print " -c [float]  : coverage above this value [20]\n";
	print " -i [float]  : abs(dI) above this value [0.2]\n";
	print "\n";
	exit (1);
}

my ($inConf) = @ARGV;

my $outBase = basename($inConf);
$outBase =~ s/\.\w+$//;

if ($outDir ne '')
{
	runBashCommand ("mkdir $outDir");
	#$outBase = "$outDir/$outBase";
}
else
{
	$outDir = ".";
}

my @infiles = ();
my @groups = ();

print "### Cache dir: $cache\n\n";

Carp::croak "$cache already exists\n" if -d $cache || -f $cache;
my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir at $cache\n" unless $ret == 0;

my $expression = 1;
$expression = 0 if ($type =~ /^as$/i);


=obsolete
my %gofiles = ();
if ($goConf ne "")
{
	open (GOIN, $goConf) || Carp::croak "Cannot open $goConf for reading\n";

	my $i = 0;
	while (my $line = <GOIN>)
	{
		$i++;
		chomp $line;
		my ($file, $key) = split(/\t/, $line);
		$key = $i if $key eq "";
		Carp::carp "Two GO files share same abbrevation... only keeping last instance.\n" if $gofiles{$key} ne "";
		Carp::croak "$file does not exist.\n" if !(-e $file);
		$gofiles{$key} = $file;
	}

	close (GOIN);
}
=cut

# open (INCONF, $inConf) || Carp::croak "cannot open $inConf\n";
# while (my $line = <INCONF>)
# {
# 	chomp $line;
# 	my ($infile, $group) = split(/\t/, $line);
# 	Carp::croak "ill-formed conf file $inConf\n" if ($group eq "" || $infile eq "");
# 	push (@infiles, $infile);
# 	push (@groups, $group);
	
# }
# close(INCONF);


my $groups = readConfigFile ($inConf, $base);
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
		
		print join ("\t", $inputFile, $gName), "\n";		
		$iter++;
	}
}
print "$iter samples in total.\n" if $verbose;








# my $pm = new Parallel::ForkManager(20); 

=obsolete
my $goConfOut = "$outBase.go.conf";
open (GOCONF, ">$cache/$goConfOut") || Carp::croak "cannot open file $cache/$goConfOut to write\n";
close (GOCONF);
=cut

for (my $i = 0; $i < scalar(@groupNames) - 1; $i++)
{
	for (my $j = $i+1; $j < scalar(@groupNames); $j++)
	{		
		# make pairwise conf file
		my $group1 = $groupNames[$i];
		my $group2 = $groupNames[$j];

		my $outBasePair .= "$outBase.$group1-$group2";
		my $outConfPair = "$outDir/$outBasePair.conf";
		# my $outConf = "$outBase.p${i}_${j}.conf";
		# print "$outConf\n";

		open (OUTCONF, ">$outConfPair") || Carp::croak "cannot open file $outConfPair to write\n";
		my $samples = $groups->{$group2}->{"samples"};
		foreach my $s (@$samples)
		{
			# print OUTCONF "$base/" if $base;
			print OUTCONF "$s\t$group2\n";
		}
		my $samples = $groups->{$group1}->{"samples"};
		foreach my $s (@$samples)
		{
			# print OUTCONF "$base/" if $base;
			print OUTCONF "$s\t$group1\n";
		}
		close (OUTCONF);

		my $spliceFile = "";
		my $spliceFile = ".gene" unless $expression;

		# open (GOCONF, ">>$cache/$goConfOut") || Carp::croak "cannot open file $cache/$goConfOut to write\n";
		# print GOCONF "$outBasePair.up$spliceFile.txt\t$group1 vs $group2 up\n";
		# print GOCONF "$outBasePair.dn$spliceFile.txt\t$group1 vs $group2 dn\n";
		# close (GOCONF);

		# $pm->start and next;

		if ($expression)
		{
			my $disp = "tagwise";
			print scalar($groups->{$group2}->{"samples"}) ."\n";
			$disp = "common" if (scalar(@{$groups->{$group1}->{"samples"}}) < 2 || scalar(@{$groups->{$group2}->{"samples"}} < 2));
			print "Using $disp dispersion\n";
		
			print "Testing pairwise differential expression ($group1 vs $group2)...\n";
			runBashCommand ("perl $cmdDir/test_expr_diff.pl -rpkm -base $base -disp $disp -MAplot $outDir/$outBasePair.pdf -v $outConfPair $outDir/$outBasePair.txt");

			# filter results
			runBashCommand ("perl $cmdDir/filterExprAS.pl -e 0 -p $pValue -q 5 -f $fc -g 2 -r 6,7 -h 1 -n -m $outDir/$outBasePair.medGene.txt -u $outDir/$outBasePair.updn.txt $outDir/$outBasePair.txt $outDir/$outBasePair.filter.txt");

			runBashCommand ("cut -f1 $outDir/$outBasePair.up.txt >> $cache/$outBase.diffExprGenes.all.txt");
			runBashCommand ("cut -f1 $outDir/$outBasePair.dn.txt >> $cache/$outBase.diffExprGenes.all.txt");
			runBashCommand ("cat $outDir/$outBasePair.medGene.txt >> $cache/$outBase.medGene.all.txt");
		    # exec '/bin/sleep', $p or die;
		}
		else
		{
			print "Testing pairwise differential splicing...\n";
			runBashCommand("perl $cmdDir/test_splicing_diff.pl -base $base -type cass --min-cov $coverage --id2gene2symbol $id2gene2symbolFile -v $outConfPair $outDir/$outBasePair.txt");
			#right now the id2gene2symbol file is mandatory, which needs to be fixed			

			# filter results
			runBashCommand("perl $cmdDir/filterExprAS.pl -e 4 -p $pValue -q 14 -c $coverage -d 9 -i $dI -j 12 -n -u $outDir/$outBasePair.UPDN.txt -h 1 $outDir/$outBasePair.txt $outDir/$outBasePair.filter.txt");
			runBashCommand("perl $cmdDir/filterExprAS.pl -e 4 -c $coverage -d 9 -h 1 $outDir/$outBasePair.txt $cache/$outBasePair.mincov.txt");
			runBashCommand("awk '\$28==1' $cache/$outBasePair.mincov.txt | cut -f1 | cut -f1 -d '/' >> $cache/$outBase.minCov.gene.txt");
			#note that the column number might change if the id2gene2symbol is not provided, so the hardcoded column number needs to be fixed
			#$NF will get the last column

			# for i in `ls *mincov20*`; do echo $i; awk '$28==1' $i | cut -f1 | cut -f1 -d "/" >> allgenes.mincov20.txt; done

			runBashCommand ("cut -f5 $outDir/$outBasePair.up.txt >> $cache/$outBase.diffExprExons.all.txt");
			runBashCommand ("cut -f5 $outDir/$outBasePair.dn.txt >> $cache/$outBase.diffExprExons.all.txt");

			runBashCommand ("cut -f1 $outDir/$outBasePair.up.txt | cut -f1 -d '/' >> $outDir/$outBasePair.up.gene.txt");
			runBashCommand ("cut -f1 $outDir/$outBasePair.dn.txt | cut -f1 -d '/' >> $outDir/$outBasePair.dn.gene.txt");

		}
		# $pm->finish;
	}
}

# $pm->wait_all_children;


if ($expression)
{
	# unique genes expressed above median in at least one condition
	runBashCommand("sort -u $cache/$outBase.diffExprGenes.all.txt > $outDir/$outBase.diffExprGenes.all.uniq.txt");
	runBashCommand("sort -u $cache/$outBase.medGene.all.txt > $outDir/$outBase.medGene.all.uniq.txt");

=too_much_stuff
	# generate expression matrix
	my $origMat = "$outBase.expr.mat.txt";
	runBashCommand("perl ~/czlab_src/countit/gen_expression_matrix.pl -base $base -v $inConf $origMat");

	# filter expression matrix
	my $filterMat = "$outBase.expr.mat.filter.txt";
	runBashCommand("sed -re 's/([[:digit:]]+)/" . '^\1\\\\b' . "/' $outBase.diffExprGenes.all.uniq.txt > " . "$cache/$outBase.diffExprGenes.all.uniq.grep.txt");
	runBashCommand("grep -f $cache/$outBase.diffExprGenes.all.uniq.grep.txt $origMat | cat > $cache/$filterMat");
	runBashCommand("cat <(head -1 $origMat) $cache/$filterMat > $filterMat");
=cut

	# GO matrix

	# if ($goConf ne "")
	# {
	# 	print "Performing GO analysis...\n";
	# 	runBashCommand("mv $cache/$goConfOut $goConfOut");
		

	# 	# my $pm = new Parallel::ForkManager(20); 

	# 	foreach my $key (keys %gofiles)
	# 	{
	# 		# $pm->start and next;
	# 		runBashCommand("perl ~/smw_scripts/DevCortex/genGOmatrix.pl -e 0 -h 1 -c $gofiles{$key} $outBase.medGene.all.uniq.txt $goConfOut $outBase.go.$key.expr.mat.txt");
	# 		# $pm->finish;
	# 	}

	# 	# $pm->wait_all_children;
	# }
}
else
{
	# unique exons meeting minimum coverage in at least one condition
	runBashCommand("sort -u $cache/$outBase.minCov.gene.txt > $outDir/$outBase.minCov.gene.uniq.txt");

	# unique diff expr exons
	runBashCommand("sort -u $cache/$outBase.diffExprExons.all.txt > $outDir/$outBase.diffExprExons.all.uniq.txt");

=too_much_stuff
	# generate expression matrix
	my $origMat = "$outBase.cass.mat.txt";
	runBashCommand("perl $cmdDir/countit/gen_splicing_matrix.pl --min-cov $coverage -v -type cass --max-std 0.1 --id2gene2symbol /ifs/data/c2b2/cz_lab/dbCASE/db_mm10/as_canonical/Mm.seq.all.AS.chrom.can.id2gene2symbol -base $base $inConf $origMat");
	

	# filter expression matrix

	my $filterMat = "$outBase.cass.mat.filter.txt";
	runBashCommand("sed -re 's/([[:digit:]]+)/" . '^\1\\\\b' . "/' $outBase.diffExprExons.all.uniq.txt >" . "$cache/$outBase.diffExprExons.all.uniq.grep.txt");
	runBashCommand("grep -Ff $outBase.diffExprExons.all.uniq.txt $origMat | cat > $cache/$filterMat");
	runBashCommand("cat <(head -1 $origMat) $cache/$filterMat > $filterMat");
=cut

	#again, it is not a good idea trying to do everything in a single script
	#also, hard coded dir is not acceptible in code for public release
	
	# GO matrix

	# if ($goConf ne "")
	# {
	# 	print "Performing GO analysis...\n";
	# 	runBashCommand("mv $cache/$goConfOut $goConfOut");
		

	# 	# my $pm = new Parallel::ForkManager(20); 

	# 	foreach my $key (keys %gofiles)
	# 	{
	# 		# $pm->start and next;
	# 		runBashCommand("perl ~/smw_scripts/DevCortex/genGOmatrix.pl -e 0 -h 1 -c $gofiles{$key} $outBase.minCov.gene.uniq.txt $goConfOut $outBase.go.$key.AS.mat.txt");
	# 		# $pm->finish;
	# 	}

	# 	# $pm->wait_all_children;
	# }
}


# Genes that were differentially expressed in at least one pair of conditions:

# 	rm diffExprGenes.all.txt
# 	for i in `ls *filter.txt`; do tail -n+2 $i | awk '{if ($9!=0) {print $1}}' >> diffExprGenes.all.txt; done
# 	sort -u diffExprGenes.all.txt > diffExprGenes.all.uniq.txt

# Expr matrix
	
# Filter gene expression matrix:

# 	# word boundaries for grep
# 	sed -re 's/([[:digit:]]+)/^\1\\b/' diffExprGenes.all.uniq.txt > diffExprGenes.all.uniq.grep.txt
# 	grep -f pairwise/diffExprGenes.all.uniq.grep.txt DevCortex.expr.ab.txt > DevCortex.expr.ab.filter.txt

# 	# header row
# 	cat <(head -1 DevCortex.expr.ab.txt) DevCortex.expr.ab.filter.txt > tmp.txt
# 	mv tmp.txt DevCortex.expr.ab.filter.txt

# Filtered correlation matrix

# 	dcAS <- read.delim("DevCortex.expr.ab.filter.txt")
# 	write.table(cor(dcAS[,4:21], use="complete.obs"), file="DevCortex.expr.ab.filter.cor.txt",sep="\t",col.names=TRUE,row.names=TRUE)


my $end_run = time();
my $run_time = $end_run - $start_run;
print "Runtime: $run_time seconds\n";



sub runBashCommand
{
	my ($cmd) = @_;
	# my $ret = system ($cmd);
	# c flag allows commands with process substitution to be run
	my @args = ( "bash", "-c", $cmd );
	my $ret = system(@args);
	Carp::croak "cmd=$cmd failed\n" unless $ret == 0;
	print "cmd=$cmd\n";
	return $ret;	
}




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
		Carp::croak "ill-formed conf file $inConf\n" if ($groupName eq "" || $sampleName eq "");
		$groups{$groupName}->{"id"} = $i++ unless exists $groups{$groupName};
		push @{$groups{$groupName}->{"samples"}}, $sampleName;
	}
	close ($fin);
	return \%groups;
}

