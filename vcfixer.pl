#!/usr/bin/perl
use strict ; use warnings;

# vcfixer   			# by M'Óscar 
my $version = "vcfixer_1.5.pl";
my $edited = "30.VII.2020";

############################

# Use this script to clean VCF files from missing values
# It works as dataset_fixer.sh and vcf_fixer
# For options and other information check the help information with "help" or "h" or so...
# vcf_missing.pl --help
# Usage: https://github.com/M-Osky/VCFixer



###########################################################################################################################
###############    Chagelog:                                                                       ########################
###############                                                                                    ########################
###############                                                                                    ########################
###############     Version 1.5 : 2020-07-30                                                       ########################
###############    Fixed some wrong messages when samples don't mach popmap                        ########################
###############    Fixed sed, now it will input Quality Metadata in inputed missing                ########################
###############    Metadata is complex and for now no editable from command line, only in script   ########################
###############    Check inputed metadata below default values for parameters                      ########################
###############    Trying to add a message when no missing are found when inputing                 ########################
###############                                                                                    ########################
###############    New vcfixer, equivalent to vcf_fixer version 1.4: 2020-05-19                    ########################
###############    Duplicated it in order  to fix bugs while vcf_fixer is running.                 ########################
###############    As of now vcf_fixer is outdated by this beta version                            ########################
###############    Major issue addressed:                                                          ########################
###############    The inputed genotypes were just the genotype "#/#" not quality information      ########################
###############    Some softwares read the quality information missing and took it as missing      ########################
###############    Now it's adding quality information to the inputed alleles                      ########################
###############    inputed quality is low, but should be enough to pass most software filters      ########################
###############    The metadata added is bellow default parameter values                           ########################
###############    Also: trying to fix the reported populations when a popmap is parsed  (389)     ########################
###############    BOTH THIS POINTS NEEDS WORK STILL                                               ########################
###############                                                                                    ########################
###############     Version 1.3 : 2020-03-12                                                       ########################
###############    Now checks if all samples from the file exist in popmap                         ########################
###############    fixing wrong population and missing counts (for the summary)                    ########################
###############    tried to fix some wrong population and missing counts that appeared             ########################
###############    in the log file. Debug 643; 649; 655; 1148, 1150; 1279; 1287;                   ########################
###############                                                                                    ########################
###############     Version 1.2 : 2019-10-28                                                       ########################
###############    debugging, trying to improve log file                                           ########################
###############                                                                                    ########################
###############     Version 1.1 : 2019-10-13                                                       ########################
###############    added popmap option as an alternative of reading pops from sample codes         ########################
###############    improved performance and trying to make it faster (spoiler alet: didn't)        ########################
###############                                                                                    ########################
###############     Version 1.0 : 2019-09-30                                                       ########################
###############    Now when there is multiple alleles with shared maximum frequency                ########################
###############    it will pick randomly one of those as a mode.                                   ########################
###############    Deleted some resundant chuncks of code. Improved output log files               ########################
###############    it will print the date and time it started and finished running                 ########################
###############                                                                                    ########################
###########################################################################################################################





#######################   PARAMETERS   #########################
# All of them can and should be set from the command line, check the help information.

my $poplength = 2;  		#how many characters long is the population code at the beginning of the sample code
#my $poplength = 3;  		#how many characters long is the population code at the beginning of the sample code

my $emptyrate = 0.8;           #ratio of missing from which samples will be considered empty
my $miss_samples = 0.3;            #ratio of missing from which samples must be deleted
my $miss_loci = 0.3;               #ratio of missing from which loci must be deleted

my $gralmiss = "pop";  		#how should the program replace the missing values?
#my $gralmiss = "global";  		#how should the program replace the missing values?
#my $gralmiss = "pop";  		#how should the program replace the missing values?
#my $gralmiss = "miss";  		#how should the program replace the missing values?
#my $gralmiss = "2/2";  		#how should the program replace the missing values?

my $misspop = "global";  		#how should the program replace when a SNP is missing in a entire population?
#my $misspop = "miss";  		#how should the program replace when a SNP is missing in a entire population?
#my $misspop = "global";  		#how should the program replace when a SNP is missing in a entire population?
#my $misspop = "2/2";  		#how should the program replace when a SNP is missing in a entire population?

my $minpop = 8;  		#minnimum number of samples per population in order to keep the population.

#output file name
my $tail = "not defined";  		#something to add at the end of the input name to generate the output name
#my $tail = "_something";  		#something to add at the end of the input name to generate the output name

my $outname = "not defined";  		#name or path for the output file
#my $outname = "podarcis";  		#name or path for the output file


# This should be fixed for all VCF produced by Populationns (Stacks)
#my $inputname = "populations.snps.vcf";  		# input file name, should be a string either alphanumeric or alphabetic.
my $inputname = "populations.snps.vcf";  		# input file name, should be a string either alphanumeric or alphabetic.

my $infocols = 9;  		# check how many columns of information has the VCF file before the first sample column
#my $infocols = 9;  		# check how many columns of information has the VCF file before the first sample column

my $sumpath = "summary_table_vcf_fixer.txt";  		#path to a file that will gather some of the details of ref_maps and vcf_fixer outputs
#my $sumpath = "summary_table.txt";  		#path to a file that will gather some of the details of ref_maps and vcf_fixer outputs
#my $sumpath = "no";  		#path to a file that will gather some of the details of ref_maps and vcf_fixer outputs

my $poplog = "ref_map.log";  		#name of the log file from populations run
#my $poplog = "populations.log";  		#name of the log file from populations run
#my $poplog = "ref_map.log";  		#name of the log file from populations run

my $popmap = "No default";  		# pop map to know the samples
#my $popmap = "none";
#my $popmap = "none";

my $quality = 0;		# input metadata about genotype quality



#################################################################################################################
###### quality metadata to add to inputed genotypes
#my $homoref = "\t0\/0:1:1,0:12:-0.05,-0.82,-1.42";
my $homoref = "\t0/0:1:1,0:12:-0.05,-0.82,-1.42";
#my $hetegen = "\t0\/1:1:1,1:12:-1.01,-0.05,-1.01";
my $hetegen = "\t0/1:1:1,1:12:-1.01,-0.05,-1.01";
#my $homoalt = "\t1\/1:1:0,1:12:-1.42,-0.82,-0.05";
my $homoalt = "\t1/1:1:0,1:12:-1.42,-0.82,-0.05";

#################################################################################################################




my $rawname = $inputname;
$rawname =~ s/^(.*?)\.vcf/$1/;

#Help!
my %arguments = map { $_ => 1 } @ARGV;
if(exists($arguments{"help"}) || exists($arguments{"--help"}) || exists($arguments{"-help"}) || exists($arguments{"h"}) || exists($arguments{"-h"}) || exists($arguments{"--h"})) {
	die "\n\n\t   $version   Help Information $edited\n\t--------------------------------------------------\n
	This program will read a VCF file, first will delete empty or monomorphic SNP, then empty samples, 
	then will delete loci with high missing rate, then samples with high missing rate,
	and finally  will input genotypes in the missing values left \"./.\".
	It is designed to work with the VCF file generated by the program \"populations\" (Stacks).
	Populations can be identified from sample names or from a tab separated popmap (Stacks format: one population and one sample per row).
	A log file, a list with deleted and kept SNPs, deleted and kept samples, and a new popmap will be saved saved in the output directory\n
	\n\tCommand line arguments and defaults:\n
	--input / --vcf           name (or path) of the VCF file. Default: $inputname\n
	--infocols                number of locus information columns before the first sample column. Default: $infocols\n
	--poplength               [int] how many characters of the sample name belong to the population code? Default: $poplength
	--popmap                  Alternatively provide a popmap to read which sample belongs to which population. $popmap\n
	--empty                   [float] missing rate from which a sample will be considered \"empty\" and deleted. Default: $emptyrate\n
	--miss_loci               [float] missing rate from which a loci should be deleted. Default: $miss_loci\n
	--miss_samples            [float] missing rate from which a sample should be deleted. Default: $miss_loci\n
	--minpop                  [int] minimum number of samples a population must have in order to keep it. Default: $minpop\n
	--gral_miss               How to replace the regular missing values? There are four options. Default: $gralmiss\n\t\t\t\t  \"pop\" to replace it with the population mode* (most frequent genotype)
	\t\t\t  \"global\" to replace the missing with the whole dataset mode*.
	\t\t\t  \"miss\" to leave it as missing. \"2/2\" or any other value to input that value.\n
	--pop_miss                What to input if a SNP is missing in an entire population? Three options available. Default: $misspop\n\t\t\t\t  \"global\" to input the global mode*, \"miss\" to keep them as missing,
	\t\t\t  \"2/2\", or \"5/5\", or any other value: to input a new genotype and remark its difference from the rest.\n
	--noquality               [flag] add this if you do not want $version to input quality/probability information in missing
	                          Some software can not process genotypes without quality metadata, by default will be also inputed,
	                          quality information is not editable from command line and should be handled carefully, values are:
	                          $homoref\n\t                          $hetegen\n\t                          $homoalt\n
	--summary                 path/name for a summary table that will gather some details of populations and vcf_fixer outputs.
	                          if \'--summary no\' the file will not be created. Default: $sumpath\n
	--poplog                  [optional] path or name of the populations (Stacks) log file. Default: $poplog
	                          if \'--poplog no\' it will not look for a logfile. Only used for summary table.
	                          If no path (only name) provided will look for file in vcf file location.\n 
	--out                     path or new name for the output file. By default will be \"input name\" + \"tail\" + \".vcf\"\n
	--tail                    Something to the end of the input file name to generate the output name.
	\t\t\t  If no tail provided, will add to the file name the sample number, SNP number, and missing handling:
	\t\t\t   1p: regular missing replaced with population mode
	\t\t\t   1g: regular missing replaced with global mode
	\t\t\t   1m: regular missing not handled - left as missing.
	\t\t\t   1x: regular missing replaced with custome input
	\t\t\t   0m: if not regular missing found in the file
	\t\t\t   2g: when a SNP is missing in the whole population, the global mode is input
	\t\t\t   2m: when missing in the whole population is left as missing
	\t\t\t   2x: custome input 
	\t\t\t   0p: there are not SNPs entirely missing in any population\n\n
	Command line call example:\n\tvcfixer.pl --input /usr/home/refmap/populations.snps.vcf --poplength 3 --gral_miss global --minpop 1 --summary no\n\n
		* When two or more alleles have the highest frequency, one of them will be picked randomly for each SNP\n\n\n";
}



################ PASSING ARGUMENTS

my $popstring = "not_defined";
use Getopt::Long;

GetOptions( "input=s" => \$inputname,    #   --input
            "vcf=s" => \$inputname,      #   --vcf
            "infocols=s" => \$infocols,      #   --infocols
            "empty=f" => \$emptyrate,      #   --empty
            "miss_loci=f" => \$miss_loci,      #   --miss_loci
            "miss_samples=f" => \$miss_samples,      #   --miss_samples
            "tail=s" => \$tail,      #   --tail
            "out=s" => \$outname,      #   --out
            "summary=s" => \$sumpath,      #   --summary
            "popmap=s" => \$popmap,      #   --popmap
            "poplength=i" => \$poplength,      #   --poplength
            "minpop=i" => \$minpop,      #   --minpop
            "poplog=s" => \$poplog,      #   --poplog
            "pop_miss=s" => \$misspop,      #   --pop_miss
			"noquality" => \$quality,        # --noquality
            "gral_miss=s" => \$gralmiss );   #   --gral_miss

my @populations = ();  		#needs to go







############### DIRECTORY PATH

use Cwd qw(cwd);
my $localdir = cwd;

my @directorypath = split('/' , $inputname);
my $pathlength = scalar @directorypath;
my $filename="nofilename";
my $subdirpath="nodirpath";
my $filepath="nofilepath";

if ($pathlength > 1) {
	$filename = $directorypath[-1];
	$filepath = $inputname;
	pop (@directorypath);
	$subdirpath = join ('/' , @directorypath);
}
elsif ($pathlength <= 1) {
	$filename = $inputname;
	$filepath = "$localdir" . "/" . "$inputname";
	$subdirpath = $localdir;
}

#print "\nPopulations to extract: $populations[0], $populations[1], $populations[2]\n";

#same with poplog
my $logname="nofilename";
my $subdirlog="nodirpath";
my $filelog="nofilepath";

if ($poplog ne "no") {
	my @logpath = split('/' , $poplog);
	my $loglength = scalar @logpath;
	if ($loglength > 1) {
		$logname = $logpath[-1];
		$filelog = $poplog;
		pop (@logpath);
		$subdirlog = join ('/' , @logpath);
	}
	elsif ($loglength <= 1) {
		$logname = $poplog;
		$subdirlog = $subdirpath;
		$filelog = "$subdirlog" . "/" . "$logname";
	}
}

print "\nReading $filename";
open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";
print "\n\nAnalizing loci information... \n";

##die "Log file: $logname, Log path: $subdirlog\nfull path log file: $filelog\n\n";

#use DateTime;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$mon = $mon + 1;
$year = $year + 1900;
my $datetime = "$year-$mon-$mday $hour:$min:$sec";

#my $datetime = "2019-mm-dd hh:mm:ss";			   # Seriously this is just temporary

my $mainparams = "--infocols => $infocols, --empty => $emptyrate, --miss_loci => $miss_loci, --miss_samples => $miss_samples\n--tail => $tail, --out => $outname, --summary => $sumpath, --poplength => $poplength\n--minpop => $minpop, --pop_miss => $misspop, --gral_miss => $gralmiss, --popmap => $popmap\n\n";

my @report=("\t$version log file", "\n", "input file name: $inputname", "Working directory: $localdir", "Executed on $datetime", " ", "Options:", "$mainparams");





my @introrows = ();  		#first first rows of metadata will be stored here
my @headers = ();  		#column headers will be stored here
my @locirows = ();  		#this will include all the non monomorphic loci information, one variable=one row

my $firstsample = $infocols;  		#first column with genotype information
my $infoloci = $infocols - 1;  		#last column with locus information
my $samplenum="none";  		# this variable will store the number of columns with SNP data
my $sampletot="none";  		# this variable will store the number of columns with sample tags
my $iniloci = 0;
my $goodloci = 0;
my $delete = 0;
my %uniquepops = ();
my $lastsample=0;
my @deletedloci = ("name\tScaffold\tposition\tID\tproblem");
my @locilist = ();
my @samplelist = ();




######## SAVE MOST COMMON GENOTYPE


while (<$VCFFILE>) {
	chomp;	#clean "end of line" symbols
	next if /^$/;  		#skip if blank
	next if /^\s*$/;  		#skip if only empty spaces
	my $line = $_;  		#save line
	my %better =();
	$line =~ s/\s+$//;  		#clean white tails in lines
	my @wholeline= split("\t", $line);  		#split columns as different elements of an array
	my $numcolumn = scalar @wholeline;
	$lastsample = $numcolumn -1;
	#if ($iniloci < 20) { print "\n$wholeline[2]: "; }
	if ($wholeline[0]=~ /^##.*?/) { 
		push (@introrows, $line);  		#save the first rows with metadata (they all start with "##")
		#print "Saving metadata $wholeline[0]\n";
	}
	elsif ($wholeline[0]=~ /^#.*?/) {
		@headers = @wholeline;
		foreach my $samplecols ($firstsample..$lastsample) {
			$sampletot = $numcolumn - $firstsample;
			my $popcode = substr ($wholeline[$samplecols], 0, $poplength);
			my $samplecode = $wholeline[$samplecols];
			push(@samplelist, $samplecode);
			$uniquepops{$popcode} = 1;
		}
	}
	else {
		$iniloci++;
		my $missingrate = 0;
		my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
		my @genotype = @SNP;  		#duplicate
		my @lociid = split (":", $wholeline[2]);
		my $scaffold = $wholeline[0];
		$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
		my $lociname = "$scaffold" . "_$lociid[0]";
		$samplenum = scalar @SNP;
		#print "\n$iniloci) sample tags: $sampletot -> genotypes at $lociname: $samplenum\n\n";
		if ($samplenum != $sampletot) {die "ERROR in VCF file! Discrepancy with number of samples\n\n$sampletot sample tags found: @samplelist\n\n$samplenum genotypes found for locus $lociname\n@SNP\n\nProgram will exit now. You should check your input file $inputname.\n\n";}
		my $last = $samplenum - 1;
		foreach my $count (0..$last) { 
			my @splitted = split(':', $genotype[$count]);  		#save only the genotype
			$genotype[$count] = $splitted[0];
		}
		#if ($iniloci < 20) {print "\n@genotype\n"; }
		foreach my $rowinfo (@genotype) {
			my $onegenotype = $rowinfo;
			
			#if ($onegenotype eq "./.") { next; }
			if (exists $better{$onegenotype}) { $better{$onegenotype} = $better{$onegenotype} +1; }
			else { $better{$onegenotype} = 1; }
			
		}
	
		my $missingnum = 0;
		if (exists $better{'./.'}) { $missingnum = $better{'./.'} ; }
		$missingrate = $missingnum / $samplenum;
		
		if (exists $better{'./.'}) { delete $better{'./.'} ; }
		#my $variability = %better;  		# Till now this worked, but is now giving back a fraction
		my $variability = scalar keys %better;
		#print "\nvariability = $variability\n";
		
		#die "\n\nDebugging\n\n";
		if($variability==1) {
			$delete++;
			#if ($iniloci < 20) { print "\n$iniloci) $lociname variability $variability". ":\n@genotype\nMonomorphic!\n"; }
			my $lociinfo = "$lociname\t$scaffold\t$wholeline[1]\t$wholeline[2]\tmonomorphic";
			push(@deletedloci, $lociinfo);
		}
		elsif ($missingrate >= $emptyrate) {
			$delete++;
			#if ($iniloci < 20) { print "\n$iniloci) $lociname variability $variability". ":\n@genotype\nToo mutch missing ($missingrate)\n"; }
			my $lociinfo = "$lociname\t$scaffold\t$wholeline[1]\t$wholeline[2]\tmissing rate $missingrate";
			push(@deletedloci, $lociinfo);
		}
		else {
			push (@locirows, $line);
			#my $iniloci = "$lociname\t$scaffold\t$wholeline[1]\t$wholeline[2]";
			#push(@locilist, $lociinfo);
			$goodloci++;
		}
		#$iniloci++;
	}
}

close $VCFFILE;

my $sampnum = $lastsample - $firstsample + 1;

my %popref = ();

if ($popmap ne "No default") {
	%uniquepops = ();
	open my $MAP, '<', $popmap or die "\nUnable to find or open $popmap: $!\n";
	while (<$MAP>) {
			my $line = $_;
			my $row = $line;
			$row =~ s/\s+$//;		#clean white tails in lines
			$row =~ s/^\s+//;	#clean white spaces at the beginning
			my @pair = split('\t', $row);
			$popref{$pair[0]} = $pair[1];
	}
	foreach my $key (keys %popref) {
		my $popcode = $popref{$key};
		$uniquepops{$popcode} = 1;
	}
}


my $popnum = scalar keys %uniquepops;		  #added "scalar keys" just in case
my @poplist = keys %uniquepops;  		#list of populations

if ($delete > 0) {print "Monomorphic or \"empty\" loci were deleted.\n$delete deleted from a total of $iniloci; $goodloci loci kept.\n" }
elsif ($delete == 0) { print "No monomorphic or empty loci found. "; }

if ($goodloci == 0) { die "\n\nWARNING: There are no loci left!\nCheck that the dataset and program options are corect.\nCheck the program help information if needed:  $version -h\n\n"; }

print "Now looking for \"empty\" samples... \n";

my $popindex = $popnum - 1;

#my @mostfreq = values %refgenotyped;
#print "\nGlobal averages: @mostfreq\n";


my $q=0;
#my $popmoda = "undefined";
my $fixedline = "undefined";
my @oldloci = @locirows;
my $snpnumber = scalar @locirows;
my $lastindex = $snpnumber -1;
my $pophandler = "nd";
my $gralhandler = "nd";
my @firstcols = ();
my @firstheaders = @headers[0..$infoloci];

#empty samples

#get info per sample
my $rownow = 1;
my %persampleind = ();
foreach my $row (@locirows) {
	#print "\n\n$row\n";
	my @wholeline= split("\t", $row);
	my @keepinfo = @wholeline[0..$infoloci];
	my $infoline = join("\t", @keepinfo);
	push (@firstcols, $infoline);
	#save the information of each sample in a different position of a hash
	#this is equivalent to translocate
	foreach my $samp ($firstsample..$lastsample) {
		if ($rownow == 1) { 
			my $samplename = $headers[$samp];
			$persampleind{$samp} = "$samplename\t$wholeline[$samp]";
			#print "\n $samp) $persampleind{$samp}\n";
			
		}
		elsif ($rownow <= $snpnumber) { $persampleind{$samp} = "$persampleind{$samp}\t$wholeline[$samp]"; }
	}
	$rownow++;
}
# now check the missingrate
my @deletedsamples = ();
my %goodsamples1 =();
my @samplecols1 = ();
my @samplenames = ();

#print "\none full: $persampleind{10}\n\n";

foreach my $samp ($firstsample..$lastsample) {
	my $locinum = 0;
	my $sampleline = $persampleind{$samp};
	my @wholeline= split("\t", $sampleline);
	my $sampleid = $wholeline[0];
	shift @wholeline;
	#print "ALL: @wholeline\n";
	my %eachsamp = ();
	$eachsamp{'./.'} = 0;
	foreach my $genot (@wholeline) {
		#print "$genot ";
		my @splitted = split(':', $genot);  		#save only the genotype
		my $onlygen = $splitted[0];
		if (exists $eachsamp{$onlygen}) { $eachsamp{$onlygen} = $eachsamp{$onlygen} +1; }
		else { $eachsamp{$onlygen} = 1; }
		$locinum++;
	}
	#print "\n\n";
	my $misscount = $eachsamp{'./.'};
	my $missrate = $misscount / $locinum;
	my $message = "$sampleid\tmissing rate $missrate";
	if ($missrate >= $emptyrate) { push (@deletedsamples, $message); }
	elsif ($missrate < $emptyrate) {
		$goodsamples1{$samp} = $sampleline;
		push (@samplecols1, $samp);  		#keep the unique sample col index
		push (@samplenames, $sampleid);
	}
}

@headers = (@firstheaders, @samplenames); #keep headers only for the samples that we have now
my $emptycount = scalar @deletedsamples;
if ($emptycount > 0) { print "\n$emptycount samples were considered \"empty\" and deleted (missing rate above or equal to $emptyrate).\n"; }
elsif ($emptycount == 0) { print "No samples deleted! " }

print "Now checking missing rate per loci... \n";

#untranspose (sort by locus)

my @goodloci1 = (); # all loci rows
my $indxmax = $snpnumber - 1;

foreach my $snp (0..$indxmax) {
	# get ready the loci info from each loci
	my $locusinf = $firstcols[$snp];
	#loop through info by sample
	foreach my $smp (@samplecols1) {
		my $sampleline = $goodsamples1{$smp};  		#get one sample line each time
		my @wholeline = split("\t", $sampleline);
		shift @wholeline;  		#delete the sample name
		$locusinf = "$locusinf\t$wholeline[$snp]";  		#add to the list of sample genotypes for the given locus
	}
	push (@goodloci1, $locusinf);  		#save all the info from the loci
}


#now check missing rate per loci and delete if needed
my $locinum1 = scalar @goodloci1;
my @goodloci2 =();
my $toomiss=0;
my $maxind = 0;
my @savedlist = ();
foreach my $lociline (@goodloci1) {
	my @wholeline= split("\t", $lociline);
	# get loci name
	#print "\n\n$lociline\n\n";
	my @lociid = split (":", $wholeline[2]);
	my $scaffold = $wholeline[0];
	$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
	my $lociname = "$scaffold" . "_$lociid[0]";
	# get number of loci
	my $maxcol = scalar @wholeline;
	$maxind = $maxcol - 1;
	#get only sample info
	my @genotypes = @wholeline[$firstsample..$maxind];
	my $samplecount = scalar @genotypes;
	my %eachlocus = ();
	$eachlocus{'./.'} = 0;
	foreach my $genot (@genotypes) {
		my @splitted = split(':', $genot);  		#save only the genotype
		my $onlygen = $splitted[0];
		if (exists $eachlocus{$onlygen}) { $eachlocus{$onlygen} = $eachlocus{$onlygen} +1; }
		else { $eachlocus{$onlygen} = 1; }
	}
	my $misscount = $eachlocus{'./.'};
	my $missrate = $misscount / $samplecount;
	if ($missrate >= $miss_loci) { 
		my $lociinfo = "$lociname\t$scaffold\t$wholeline[1]\t$wholeline[2]\tmissing rate $missrate";
		push(@deletedloci, $lociinfo);
		$toomiss++;
	}
	elsif ($missrate < $miss_loci) { 
		push (@goodloci2, $lociline); 
		push (@savedlist, $lociname);
	}
}
my $locinum2 = scalar @goodloci2;
if ($toomiss > 0) { print "$toomiss loci deleted with missing rate above $miss_loci\n"; }
elsif ($toomiss == 0) { print "No loci deleted! " }
#print "\n@savedlist\n";

if ($locinum2 ==0) {die "\n\nNO LOCI LEFT!!\n$version is done (and quite sad)\nMay be try different settings\n\n"; }
print "\nChecking missing rate per sample... \n";



#missing val per sample

#get info per sample
#my $newlastsample = $

@firstcols=();
$rownow = 1;
%persampleind = ();
foreach my $row (@goodloci2) {
	my @wholeline= split("\t", $row);
	#print "\n\n$row\n\n";
	my @keepinfo = @wholeline[0..$infoloci];
	my $infoline = join("\t", @keepinfo);
	#print "\n\n>> $infoline\n\n";
	push (@firstcols, $infoline);
	#save the information of each sample in a different position of a hash
	#this is equivalent to translocate
	foreach my $samp ($firstsample..$maxind) {
		if ($rownow == 1) { 
			my $samplename = $headers[$samp];
			$persampleind{$samp} = "$samplename\t$wholeline[$samp]";
			#print "\n $samp) $persampleind{$samp}\n";
			
		}
		elsif ($rownow <= $locinum2) { $persampleind{$samp} = "$persampleind{$samp}\t$wholeline[$samp]"; }
	}
	$rownow++;
}

#print "\n\n$persampleind{10}\n\n";

#foreach (@firstcols) {
#	print "\n\n > $_\n\n";
#}


#check missing
my %goodsamples2 =();
my @samplecols2 = ();
@samplenames =();
my $delsam=0;
foreach my $samp ($firstsample..$maxind) {
	my $locinum = 0;
	my $sampleline = $persampleind{$samp};
	my @wholeline= split("\t", $sampleline);
	my $sampleid = $wholeline[0];
	shift @wholeline;
	my %eachsamp = ();
	$eachsamp{'./.'} = 0;
	foreach my $genot (@wholeline) {
		my @splitted = split(':', $genot);  		#save only the genotype
		my $onlygen = $splitted[0];
		if (exists $eachsamp{$onlygen}) { $eachsamp{$onlygen} = $eachsamp{$onlygen} +1; }
		else { $eachsamp{$onlygen} = 1; }
		$locinum++;
	}
	my $misscount = $eachsamp{'./.'};
	my $missrate = $misscount / $locinum;
	if ($missrate >= $miss_samples) {
		my $message = "$sampleid\tmissing rate $missrate";
		push (@deletedsamples,$message);
		$delsam++;
	}
	elsif ($missrate < $emptyrate) {
		$goodsamples2{$samp} = $sampleline;
		#print "\n\n$sampleline\n\n";
		push (@samplecols2, $samp);  		#keep the unique sample col index
		push (@samplenames, $sampleid);
	}
}
if ($delsam > 0) { print "$delsam sample(s) deleted with missing rate above $miss_samples.\n"; }
elsif ($delsam == 0) { print "No samples deleted! " }

print "Now checking populations... \n";

#trim populations

my @delepops = ();
my @keepops = ();
my %popcount = ();
my @deletelist = ();

my $checksamp = scalar @samplecols2;
if ($checksamp == 0) { die "\n\nNO SAMPLES LEFT!!\n$version is done (and quite sad)\nMay be try different settings\n\n"; }


if ($popmap eq "No default") {
	foreach my $samp (@samplecols2) {
		my $sampline = $goodsamples2{$samp};
		my $poppart = substr ($sampline, 0, $poplength);
		if (exists $popcount{$poppart}) { $popcount{$poppart} = $popcount{$poppart} + 1; }
		else { $popcount{$poppart} = 1 }
		#print "$popcount{$poppart}\n";
	}
}
else {
	foreach my $samp (@samplecols2) {
		my $sampline = $goodsamples2{$samp};
		my @letsdothisagain = split ("\t", $sampline);
		if (exists $popref{$letsdothisagain[0]}) {
			my $gotpop = $popref{$letsdothisagain[0]};
			#print "sample $letsdothisagain[0] -> pop $gotpop.\n";
			if (exists $popcount{$gotpop}) { $popcount{$gotpop} = $popcount{$gotpop} + 1; }
			else { $popcount{$gotpop} = 1; }
		}
		else { push (@report, "WARNING!: Sample $letsdothisagain[0] does not exist in popmap and will not be analysed\n"); }
	}	
}


my $test = 0;
push (@report, "$popnum populations found: @poplist");
foreach my $popul (@poplist) {
	if (exists $popcount{$popul}) {
		if ($popcount{$popul} < $minpop) { push (@delepops, $popul); }
		elsif ($popcount{$popul} >= $minpop) { push(@keepops, $popul); }
	}
	else { push (@report, "Population $popul from popmap not found in the vcf file"); $test++; }
}
#if ($test > 0) { push (@report, ""); }
push (@report, "");
$test = 0;

my $smallpops = scalar @delepops;
my @samplecols3 = ();
@samplenames = ();
if ($smallpops > 0) {
	foreach my $shitpop (@delepops) {
		foreach my $samp (@samplecols2) {
			my $sampline = $goodsamples2{$samp};
			my @wholeline = split("\t", $sampline);
			my $sampleid = $wholeline[0];
			my $poppart = substr ($sampline, 0, $poplength);
			if ($popmap ne "No default") { $poppart = $popref{$sampleid}; }
			my $message = "$sampleid\tbelongs to $poppart with less than $minpop samples per population";
			if ($poppart eq $shitpop) { push(@deletedsamples, $message); }
			else {
				push (@samplecols3, $sampline);
				push (@samplenames, $sampleid);
				#print "\n\n >> $sampline\n\n";
			}
		}
	}
	if ($smallpops == 1) { print "Population @delepops had less than $minpop samples and was deleted.\n" }
	elsif ($smallpops > 1) { print "Populations @delepops, had less than $minpop samples and were deleted.\n" }
}
elsif ($smallpops == 0) { 
	print "All populations have enough samples.\n";
	foreach my $samp (@samplecols2) {
		my $sampline = $goodsamples2{$samp};
		my @wholeline = split("\t", $sampline);
		my $sampleid = $wholeline[0];
		push (@samplecols3, $sampline);
		push (@samplenames, $sampleid);
	}
}

@headers = (@firstheaders, @samplenames); #keep headers only for the samples that we have now

my $checkpops = scalar @keepops;
if ($checkpops == 0) { die "\n\nNO POPULATIONS LEFT!!\n$version is done (and quite sad)\nMay be try different settings\n\n"; }

print "\nAnalizing and inputing values on remaining missing... \n";

#translocate back to per-loci dataset

my @goodloci3 = (); # all loci rows
$indxmax = $locinum2 - 1;
my $finalsamp = scalar @samplecols3;
my $lastcol = $finalsamp - 1;

foreach my $snp (0..$indxmax) {
	# get ready the loci info from each loci
	my $locusinf = $firstcols[$snp];
	#print "\n\n$locusinf\n\n";
	#loop through info by sample
	foreach my $smp (0..$lastcol) {
		my $sampleline = $samplecols3[$smp];  		#get one sample line each time
		#print "\n\n > $sampleline\n\n";
		my @wholeline = split("\t", $sampleline);
		shift @wholeline;  		#delete the sample name
		$locusinf = "$locusinf\t$wholeline[$snp]";  		#add to the list of sample genotypes for the given locus
	}
	push (@goodloci3, $locusinf);  		#save all the info from the loci
	#print "\n\n$locusinf\n\n";
}

#####

#rename everything to their old names
my $locileft = scalar @goodloci3;
$lastindex = $locileft - 1;

my @checkline = split("\t", $goodloci3[0]);
my $colnum = scalar @checkline;
$lastsample = $colnum - 1;

@poplist = @keepops;
my $popleft = scalar @keepops;
$popindex = $popleft - 1;

@locirows = @goodloci3;

my $snpcount = 0;




##########################################
########### input missing  ###############
##########################################


$goodloci = 0;
my %refgenotyped = ();
my $warns = 0;


# replace missing

my $gralcount = 0;
my $popwisecount = 0;

if ($gralmiss eq "pop" && $misspop eq "global") {			   # pop - glob
	
	$gralhandler = "_1p";
	
	#compute global mode per SNP
	
	foreach my $linesnp (@locirows) {
		#my %freqgenot = ();
		my %freqgenot = ('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
		my $bienmoda = "none";
		my @wholeline = split ("\t", $linesnp);
		#print "\n\n$linesnp\n\n";
		my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
		my @genotype = @SNP;
		my @lociid = split (":", $wholeline[2]);
		my $scaffold = $wholeline[0];
		$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
		my $lociname = "$scaffold" . "_$lociid[0]";
		$samplenum = scalar @SNP;
		my $last = $samplenum - 1;
		foreach my $count (0..$last) { 
				my @splitted = split(':', $genotype[$count]);  		#save only the genotype
				$genotype[$count] = $splitted[0];
			}
		foreach my $onegenotype (@genotype) {
			$freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1;
			#if (exists $freqgenot{$onegenotype}) { $freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1; }
			#else { $freqgenot{$onegenotype} = 1; }
		}
		delete $freqgenot{'./.'};
		my $refal = $freqgenot{'0/0'};
		my $hetal = $freqgenot{'0/1'};
		my $altal = $freqgenot{'1/1'};
		
		if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
			my @highest = ('0/0', '0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
			my @highest = ('0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
			my @highest = ('0/0', '0/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
			my @highest = ('0/0', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		else { for my $bt (sort {$freqgenot{$a} <=> $freqgenot{$b} || $a cmp $b } (keys %freqgenot) ) { $bienmoda = $bt; } }
		
		$refgenotyped{$goodloci} = "$bienmoda";
		$goodloci++;
	}
	
	#now check the missing population-wise
	
	#gets samples per population
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $reportmisspop = 0;
		my $popwise = 0;
		my $modanomiss = "undefined";
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $printpop = "no need to input anything";
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		
		#if ($snpcount < 7) { print "@sampleindex\n"; }
		
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		
		
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			#my %popgenot=();
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			#print "\n\n@wholerow\n@sampleindex    ($population)\n\n";
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				#if ($snpcount < 7) { print "$allgenot[$num] "; }
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			#if ($snpcount < 7) { print "\n@allgenot\n"; }
						foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				#if ($snpcount < 7) { print "$onegenotype "; }
				#if (exists $popgenot{$onegenotype}) { $popgenot{$onegenotype} = $popgenot{$onegenotype} +1; }
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
				#else { $popgenot{$onegenotype} = 1; }
				#sort the hash acording to the values, times that each genotype appears
				#for $q (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $popmoda = $q; }
			}
			#sort the hash acording to the values, times that each genotype appears
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			
			my $refal = $popgenot{'0/0'};
			my $hetal = $popgenot{'0/1'};
			my $altal = $popgenot{'1/1'};
			
			if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
				my @highest = ('0/0', '0/1', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
				my @highest = ('0/1', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
				my @highest = ('0/0', '0/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
				my @highest = ('0/0', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			else { for my $nm (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $modanomiss = $nm; } }
			my $highestal = $popgenot{$modanomiss};
			#Replace the missing
			
			
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				
				$pophandler = "2g";
				my $inputmode = $refgenotyped{$snpcount};
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					
					#my $checkrow = $wholerow[$sampleposition];
					#my @matches = ($checkrow =~ /\.\/\./g);
					#my $getmatch = @matches;
					#$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$inputmode/g;
				}
				$popwisecount++;
				$popwise++;
				$printpop = "global mode was input";
				$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
				if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
				push (@report, "$shitlocus missing in all samples");
				$test++;
			}
			elsif ($missnum >= $highestal) {
				#$pophandler = "2g";
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$reportmisspop = $reportmisspop + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$modanomiss/g;
				}
				
				if ($warns == 0) {
					print "One or more SNPs have the same or higher frequency of missing (\"./.\") than the most frequent genotype. Check Log file.\n";
					$warns++;
				}
				if ($moretest == 0) { push (@report, "\nLoci at $population with same or higher missing (./.) than the most frequent genotype"); }
				push (@report, "$wholerow[0] $wholerow[1] ($wholerow[2]) missing in $missnum samples, population most frequent genotype \"$modanomiss\" present at $highestal samples was input.");
				$moretest++;
				$gralcount++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$modanomiss/g;
					$numgral = $getmatch;
					$reportmisspop = $reportmisspop + $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		push (@report, "\nSummary population $population: $reportmisspop general missing genotypes found (population mode was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, population mode was input.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by the population mode from each SNP\n\n";
}
elsif ($gralmiss eq "pop" && $misspop eq "miss") { #pop - miss
	$gralhandler = "_1p";
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $reportmisspop = 0;
		my $popwise = 0;
		my $modanomiss = "undefined";
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		
		#if ($snpcount < 7) { print "@sampleindex\n"; }
		
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		
		
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			#my %popgenot=();
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			#print "\n\n@wholerow\n@sampleindex    ($population)\n\n";
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				#if ($snpcount < 7) { print "$allgenot[$num] "; }
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			#if ($snpcount < 7) { print "\n@allgenot\n"; }
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				#if ($snpcount < 7) { print "$onegenotype "; }
				#if (exists $popgenot{$onegenotype}) { $popgenot{$onegenotype} = $popgenot{$onegenotype} +1; }
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
				#else { $popgenot{$onegenotype} = 1; }
				#sort the hash acording to the values, times that each genotype appears
				#for $q (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $popmoda = $q; }
			}
			#sort the hash acording to the values, times that each genotype appears
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			
			my $refal = $popgenot{'0/0'};
			my $hetal = $popgenot{'0/1'};
			my $altal = $popgenot{'1/1'};
			
			if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
				my @highest = ('0/0', '0/1', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
				my @highest = ('0/1', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
				my @highest = ('0/0', '0/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
				my @highest = ('0/0', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			else { for my $nm (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $modanomiss = $nm; } }
			my $highestal = $popgenot{$modanomiss};
			#Replace the missing
			
			
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2m";
				$popwisecount++;
				$popwise++;
				$printpop = "no need to input anything";
				$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
				if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
				push (@report, "$shitlocus missing in all samples");
				$test++;
			}
			elsif ($missnum >= $highestal) {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$reportmisspop = $reportmisspop + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$modanomiss/g;
				}
				
				if ($warns == 0) {
					print "One or more SNPs have the same or higher frequency of missing (\"./.\") than the most frequent genotype). Check Log file.\n";
					$warns++;
				}
				if ($moretest == 0) { push (@report, "\nLoci at $population with same or higher missing (./.) than the most frequent genotype"); }
				push (@report, "$wholerow[0] $wholerow[1] ($wholerow[2]) missing in $missnum samples, population most frequent genotype \"$modanomiss\" present at $highestal samples was input.");
				$moretest++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$modanomiss/g;
					$numgral = $getmatch;
					$reportmisspop = $reportmisspop + $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		push (@report, "\nSummary population $population: $reportmisspop general missing genotypes found (population mode was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, population mode was input.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by the population mode from each SNP\n\n";
}
elsif ($gralmiss eq "pop") { #pop - other
	$gralhandler = "_1p";
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest =0;
		my $numgral = 0;
		my $reportmisspop = 0;
		my $popwise = 0;
		my $modanomiss = "undefined";
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $printpop = "no needed to input anything.";
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		
		#if ($snpcount < 7) { print "@sampleindex\n"; }
		
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		
		
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			#my %popgenot=();
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			#print "\n\n@wholerow\n@sampleindex    ($population)\n\n";
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				#if ($snpcount < 7) { print "$allgenot[$num] "; }
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			#if ($snpcount < 7) { print "\n@allgenot\n"; }
						foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				#if ($snpcount < 7) { print "$onegenotype "; }
				#if (exists $popgenot{$onegenotype}) { $popgenot{$onegenotype} = $popgenot{$onegenotype} +1; }
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
				#else { $popgenot{$onegenotype} = 1; }
				#sort the hash acording to the values, times that each genotype appears
				#for $q (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $popmoda = $q; }
			}
			#sort the hash acording to the values, times that each genotype appears
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			
			my $refal = $popgenot{'0/0'};
			my $hetal = $popgenot{'0/1'};
			my $altal = $popgenot{'1/1'};
			
			if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
				my @highest = ('0/0', '0/1', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
				my @highest = ('0/1', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
				my @highest = ('0/0', '0/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
				my @highest = ('0/0', '1/1');
				my $sameal = scalar @highest;
				my $choosen = int(rand($sameal));
				$modanomiss = $highest[$choosen];
			}
			else { for my $nm (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $modanomiss = $nm; } }
			my $highestal = $popgenot{$modanomiss};
			#Replace the missing
			
			
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2x";
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					#$reportmisspop = $reportmisspop + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$misspop/g;
					$popwisecount++;
					$popwise++;
					$printpop = "\"$misspop\" was input.";
					$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
					if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
					push (@report, "$shitlocus missing in all samples.");
					$test++;
				}
			}
			elsif ($missnum >= $highestal) {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$reportmisspop = $reportmisspop + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$modanomiss/g;
				}
				
				if ($warns == 0) {
					print "One or more SNPs have the same or higher frequency of missing (\"./.\") than the most frequent genotype). Check Log file.\n";
					$warns++;
				}
				if ($moretest == 0) { push (@report, "\nLoci at $population with same or higher missing (./.) than the most frequent genotype"); }
				push (@report, "$wholerow[0] $wholerow[1] ($wholerow[2]) missing in $missnum samples, \"$misspop\" was input.");
				$moretest++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$modanomiss/g;
					$reportmisspop = $reportmisspop + $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $reportmisspop general missing genotypes found (population mode was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, population mode was input.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by the population mode from each SNP\n\n";
}
elsif ($gralmiss eq "global" && $misspop eq "global") {	# global - global
	$gralhandler = "_1g";
	
	# calculate global mode
	
	foreach my $linesnp (@locirows) {
		#my %freqgenot = ();
		my %freqgenot = ('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
		my $bienmoda = "none";
		my @wholeline = split ("\t", $linesnp);
		#print "\n\n$linesnp\n\n";
		my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
		my @genotype = @SNP;
		my @lociid = split (":", $wholeline[2]);
		my $scaffold = $wholeline[0];
		$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
		my $lociname = "$scaffold" . "_$lociid[0]";
		$samplenum = scalar @SNP;
		my $last = $samplenum - 1;
		foreach my $count (0..$last) { 
				my @splitted = split(':', $genotype[$count]);  		#save only the genotype
				$genotype[$count] = $splitted[0];
			}
		foreach my $onegenotype (@genotype) {
			$freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1;
			#if (exists $freqgenot{$onegenotype}) { $freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1; }
			#else { $freqgenot{$onegenotype} = 1; }
		}
		delete $freqgenot{'./.'};
		my $refal = $freqgenot{'0/0'};
		my $hetal = $freqgenot{'0/1'};
		my $altal = $freqgenot{'1/1'};
		
		if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
			my @highest = ('0/0', '0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
			my @highest = ('0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
			my @highest = ('0/0', '0/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
			my @highest = ('0/0', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		else { for my $bt (sort {$freqgenot{$a} <=> $freqgenot{$b} || $a cmp $b } (keys %freqgenot) ) { $bienmoda = $bt; } }
		
		$refgenotyped{$goodloci} = "$bienmoda";
		$goodloci++;
	}
	
	# check population-wise
	
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $numgral = 0;
		my $moretest = 0;
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			#my %popgenot=();
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
			}
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			
			#Replace the missing
			#print "$divers ";
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2g";
				my $inputmode = $refgenotyped{$snpcount};
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					
					$wholerow[$sampleposition]=~ s/\.\/\./$inputmode/g;
				}
				$popwisecount++;
				$popwise++;
				$printpop = "global mode was input";
				$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
				if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
				push (@report, "$shitlocus missing in all samples");
				$test++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$numgral = $getmatch;
					my $inputmode = $refgenotyped{$snpcount};
					$wholerow[$sampleposition]=~ s/\.\/\./$inputmode/g;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (global mode was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, global mode was input.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by the global mode from each SNP\n";
}
elsif ($gralmiss eq "global" && $misspop eq "miss") {	# global - miss
	$gralhandler = "_1g";
	
	# calculate global mode
	
	foreach my $linesnp (@locirows) {
		#my %freqgenot = ();
		my %freqgenot = ('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
		my $bienmoda = "none";
		my @wholeline = split ("\t", $linesnp);
		#print "\n\n$linesnp\n\n";
		my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
		my @genotype = @SNP;
		my @lociid = split (":", $wholeline[2]);
		my $scaffold = $wholeline[0];
		$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
		my $lociname = "$scaffold" . "_$lociid[0]";
		$samplenum = scalar @SNP;
		my $last = $samplenum - 1;
		foreach my $count (0..$last) { 
				my @splitted = split(':', $genotype[$count]);  		#save only the genotype
				$genotype[$count] = $splitted[0];
			}
		foreach my $onegenotype (@genotype) {
			$freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1;
			#if (exists $freqgenot{$onegenotype}) { $freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1; }
			#else { $freqgenot{$onegenotype} = 1; }
		}
		delete $freqgenot{'./.'};
		my $refal = $freqgenot{'0/0'};
		my $hetal = $freqgenot{'0/1'};
		my $altal = $freqgenot{'1/1'};
		
		if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
			my @highest = ('0/0', '0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
			my @highest = ('0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
			my @highest = ('0/0', '0/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
			my @highest = ('0/0', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		else { for my $bt (sort {$freqgenot{$a} <=> $freqgenot{$b} || $a cmp $b } (keys %freqgenot) ) { $bienmoda = $bt; } }
		
		$refgenotyped{$goodloci} = "$bienmoda";
		$goodloci++;
	}
	
	# population wise
	
	foreach my $pop (0..$popindex) {  			#each population
		my $moretest = 0;
		my $numgral = 0;
		$test = 0;
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			#my %popgenot=();
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
			}
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			#Replace the missing
			#print "$divers ";
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				
				$pophandler = "2m";
				$popwisecount++;
				$popwise++;
				$printpop = "no value was input";
				$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
				if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
				push (@report, "$shitlocus missing in all samples");
				$test++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$numgral = $getmatch;
					my $inputmode = $refgenotyped{$snpcount};
					$wholerow[$sampleposition]=~ s/\.\/\./$inputmode/g;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (global mode was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, global mode was input.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by the global mode from each SNP\n";
}
elsif ($gralmiss eq "global") {	# global - other
	$gralhandler = "_1g";
	
	# calculate global mode
	
	foreach my $linesnp (@locirows) {
		#my %freqgenot = ();
		my %freqgenot = ('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
		my $bienmoda = "none";
		my @wholeline = split ("\t", $linesnp);
		#print "\n\n$linesnp\n\n";
		my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
		my @genotype = @SNP;
		my @lociid = split (":", $wholeline[2]);
		my $scaffold = $wholeline[0];
		$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
		my $lociname = "$scaffold" . "_$lociid[0]";
		$samplenum = scalar @SNP;
		my $last = $samplenum - 1;
		foreach my $count (0..$last) { 
				my @splitted = split(':', $genotype[$count]);  		#save only the genotype
				$genotype[$count] = $splitted[0];
			}
		foreach my $onegenotype (@genotype) {
			$freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1;
			#if (exists $freqgenot{$onegenotype}) { $freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1; }
			#else { $freqgenot{$onegenotype} = 1; }
		}
		delete $freqgenot{'./.'};
		my $refal = $freqgenot{'0/0'};
		my $hetal = $freqgenot{'0/1'};
		my $altal = $freqgenot{'1/1'};
		
		if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
			my @highest = ('0/0', '0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
			my @highest = ('0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
			my @highest = ('0/0', '0/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
			my @highest = ('0/0', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		else { for my $bt (sort {$freqgenot{$a} <=> $freqgenot{$b} || $a cmp $b } (keys %freqgenot) ) { $bienmoda = $bt; } }
		
		$refgenotyped{$goodloci} = "$bienmoda";
		$goodloci++;
	}
	
	# population-wise
	
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			#my %popgenot=();
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
			}
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			#Replace the missing
			#print "$divers ";
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2x";
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					
					$wholerow[$sampleposition]=~ s/\.\/\./$misspop/g;
					$popwisecount++;
					$popwise++;
					$printpop = "\"$misspop\" was input";
					$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
					if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
					push (@report, "$shitlocus missing in all samples");
					$test++;
				}
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$numgral = $getmatch;
					my $inputmode = $refgenotyped{$snpcount};
					$wholerow[$sampleposition]=~ s/\.\/\./$inputmode/g;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (global mode was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, global mode was input.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by the global mode from each SNP\n";
}
elsif ($gralmiss eq "miss" && $misspop eq "global") { # miss - global
	$gralhandler = "_1m";
	
	# calculate global mode
	
	foreach my $linesnp (@locirows) {
		#my %freqgenot = ();
		my %freqgenot = ('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
		my $bienmoda = "none";
		my @wholeline = split ("\t", $linesnp);
		#print "\n\n$linesnp\n\n";
		my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
		my @genotype = @SNP;
		my @lociid = split (":", $wholeline[2]);
		my $scaffold = $wholeline[0];
		$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
		my $lociname = "$scaffold" . "_$lociid[0]";
		$samplenum = scalar @SNP;
		my $last = $samplenum - 1;
		foreach my $count (0..$last) { 
				my @splitted = split(':', $genotype[$count]);  		#save only the genotype
				$genotype[$count] = $splitted[0];
			}
		foreach my $onegenotype (@genotype) {
			$freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1;
			#if (exists $freqgenot{$onegenotype}) { $freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1; }
			#else { $freqgenot{$onegenotype} = 1; }
		}
		delete $freqgenot{'./.'};
		my $refal = $freqgenot{'0/0'};
		my $hetal = $freqgenot{'0/1'};
		my $altal = $freqgenot{'1/1'};
		
		if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
			my @highest = ('0/0', '0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
			my @highest = ('0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
			my @highest = ('0/0', '0/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
			my @highest = ('0/0', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		else { for my $bt (sort {$freqgenot{$a} <=> $freqgenot{$b} || $a cmp $b } (keys %freqgenot) ) { $bienmoda = $bt; } }
		
		$refgenotyped{$goodloci} = "$bienmoda";
		$goodloci++;
	}
	
	# population-wise
	
	
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $modanomiss = "undefined";
		
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				#if (exists $popgenot{$onegenotype}) { $popgenot{$onegenotype} = $popgenot{$onegenotype} +1; }
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
				#else { $popgenot{$onegenotype} = 1; }
				
				#sort the hash acording to the values, times that each genotype appears
				#for $q (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $popmoda = $q; }
			}
			
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			
			#Replace the missing
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
					$pophandler = "2g";
					my $inputmode = $refgenotyped{$snpcount};
					foreach my $id2 (0..$indexofindex) {
						my $sampleposition = $sampleindex[$id2];
						
						my $checkrow = $wholerow[$sampleposition];
						my @matches = ($checkrow =~ /\.\/\./g);
						my $getmatch = @matches;
						$gralcount = $gralcount + $getmatch;
						
						$wholerow[$sampleposition]=~ s/\.\/\./$inputmode/g;
					}
					$popwisecount++;
					$popwise++;
					$printpop = "global mode was input";
					$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
					if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
					push (@report, "$shitlocus missing in all samples");
					$test++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$numgral = $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (left as missing).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, they were left as missing.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were not replaced\n";
}
elsif ($gralmiss eq "miss" && $misspop eq "miss") { # miss - miss
	$gralhandler = "_1m";
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				#if (exists $popgenot{$onegenotype}) { $popgenot{$onegenotype} = $popgenot{$onegenotype} +1; }
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
				#else { $popgenot{$onegenotype} = 1; }
				
				#sort the hash acording to the values, times that each genotype appears
				#for $q (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $popmoda = $q; }
			}
			
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			
			#Replace the missing
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2m";
				$popwisecount++;
				$popwise++;
				$printpop = "no value was input";
				$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
				if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
				push (@report, "$shitlocus missing in all samples");
				$test++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$numgral = $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (left as missing).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, they were left as missing.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were not replaced\n";
}
elsif ($gralmiss eq "miss") { # miss - other
	$gralhandler = "_1m";
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		
		my $population = $poplist[$pop];
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot=('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			# save the genotypes and sort them by frequency of appearance
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				#if (exists $popgenot{$onegenotype}) { $popgenot{$onegenotype} = $popgenot{$onegenotype} +1; }
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
				#else { $popgenot{$onegenotype} = 1; }
				
				#sort the hash acording to the values, times that each genotype appears
				#for $q (sort {$popgenot{$a} <=> $popgenot{$b} || $a cmp $b } (keys %popgenot) ) { $popmoda = $q; }
			}
			
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			#Replace the missing
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2x";
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					
					$wholerow[$sampleposition]=~ s/\.\/\./$misspop/g;
					$popwisecount++;
					$popwise++;
					$printpop = "\"$misspop\" was input";
					$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
					if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
					push (@report, "$shitlocus missing in all samples");
					$test++;
				}
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$numgral = $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (left as missing).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, they were left as missing.\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were not replaced\n";
}
elsif ( $misspop eq "global") {
	$gralhandler = "_1x";
	
	# calculate global mode
	
	foreach my $linesnp (@locirows) {
		#my %freqgenot = ();
		my %freqgenot = ('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0);
		my $bienmoda = "none";
		my @wholeline = split ("\t", $linesnp);
		#print "\n\n$linesnp\n\n";
		my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
		my @genotype = @SNP;
		my @lociid = split (":", $wholeline[2]);
		my $scaffold = $wholeline[0];
		$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
		my $lociname = "$scaffold" . "_$lociid[0]";
		$samplenum = scalar @SNP;
		my $last = $samplenum - 1;
		foreach my $count (0..$last) { 
				my @splitted = split(':', $genotype[$count]);  		#save only the genotype
				$genotype[$count] = $splitted[0];
			}
		foreach my $onegenotype (@genotype) {
			$freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1;
			#if (exists $freqgenot{$onegenotype}) { $freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1; }
			#else { $freqgenot{$onegenotype} = 1; }
		}
		delete $freqgenot{'./.'};
		my $refal = $freqgenot{'0/0'};
		my $hetal = $freqgenot{'0/1'};
		my $altal = $freqgenot{'1/1'};
		
		if ($refal == $hetal && $refal == $altal && $hetal == $altal) {
			my @highest = ('0/0', '0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal < $hetal && $refal < $altal && $hetal == $altal) {
			my @highest = ('0/1', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal == $hetal && $refal > $altal && $hetal > $altal) {
			my @highest = ('0/0', '0/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		if ($refal > $hetal && $refal == $altal && $hetal > $altal) {
			my @highest = ('0/0', '1/1');
			my $sameal = scalar @highest;
			my $choosen = int(rand($sameal));
			$bienmoda = $highest[$choosen];
		}
		else { for my $bt (sort {$freqgenot{$a} <=> $freqgenot{$b} || $a cmp $b } (keys %freqgenot) ) { $bienmoda = $bt; } }
		
		$refgenotyped{$goodloci} = "$bienmoda";
		$goodloci++;
	}
	
	# population-wise
	
	
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $population = $poplist[$pop];
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot = ('./.' => "0", '0/0' => "0", '0/1' => "0", '1/1' => "0");
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
			}
			
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			#Replace the missing
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2g";
				my $inputmode = $refgenotyped{$snpcount};
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					
					$wholerow[$sampleposition]=~ s/\.\/\./$inputmode/g;
				}
				$popwisecount++;
				$popwise++;
				$printpop = "global mode was input";
				$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
				if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
				push (@report, "$shitlocus missing in all samples");
				$test++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$gralmiss/g;
					$numgral = $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (\"$gralmiss\" was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, they were replaced with \"$gralmiss\".\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by $gralmiss\n";
}
elsif ($misspop eq "miss") {
	$gralhandler = "_1x";
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $population = $poplist[$pop];
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot = ('./.' => "0", '0/0' => "0", '0/1' => "0", '1/1' => "0");
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
			}
			
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			#Replace the missing

			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2m";
				$popwisecount++;
				$popwise++;
				$printpop = "no value was input";
				$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
				if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
				push (@report, "$shitlocus missing in all samples");
				$test++;
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$gralmiss/g;
					$numgral = $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (\"$gralmiss\" was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, they were replaced with \"$gralmiss\".\n\n");
	print "\n $gralcount general missing values (found in few samples per population) were replaced by $gralmiss\n";
}
else {
	$gralhandler = "_1x";
	foreach my $pop (0..$popindex) {  			#each population
		$test = 0;
		my $moretest = 0;
		my $numgral = 0;
		my $population = $poplist[$pop];
		my $popmoda = "undefined";
		my $modanomiss = "undefined";
		
		my @sampleindex = ();
		my $popwise = 0;
		my $printpop = 0;
		my $shitlocus = 0;
		foreach my $samp ($firstsample..$lastsample) {  			#each column with samples
			my $poppart = "atlantis";
			if ($popmap eq "No default") { $poppart = substr ($headers[$samp], 0, $poplength); }  		#extract the population code from the sample name
			elsif ($popmap ne "No default") { $poppart = $popref{$headers[$samp]}; }  		#use popmap
			if ($population eq $poppart) { push (@sampleindex, $samp); }  		#save the positions of the samples that belong to the same population
		}
		my $numindex = scalar @sampleindex;
		my $indexofindex = $numindex - 1;
		# save only alleles from all the chunck of loci data
		foreach my $fila (0..$lastindex) {
			my %popgenot = ('./.' => "0", '0/0' => "0", '0/1' => "0", '1/1' => "0");
			my @allgenot = ();
			my $lineagain = $locirows[$fila];  		# each row
			my @wholerow = split("\t", $lineagain);  		# split by column
			my $num=0;
			foreach my $id (0..$indexofindex) {
				my $sampleposition = $sampleindex[$id];  		#each sample from that line (SNP) from that population
				my @holdallele = split (':', $wholerow[$sampleposition]);  		#split the alleles from the rest of the info
				$allgenot[$num] = $holdallele[0];  		#save the allele information
				$num++;
			}
			my $numgenot = scalar @allgenot;
			my $genotindex = $numgenot - 1;
			
			foreach my $genind (0..$genotindex) {
				my $onegenotype = $allgenot[$genind];
				$popgenot{$onegenotype} = $popgenot{$onegenotype} +1;
			}
			
			#my $divers = %popgenot;
			my $missnum = $popgenot{'./.'};
			delete $popgenot{'./.'};
			
			#Replace the missing
			if($missnum >= $numgenot) {
				# Replace missing values when missing in an entire population.
				$pophandler = "2x";
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$misspop/g;
					$popwisecount++;
					$popwise++;
					$printpop = "\"$misspop\" was input";
					$shitlocus = "$wholerow[0]" . ", $wholerow[1] ($wholerow[2])";
					if ($test == 0) { push (@report, "\nThere are SNPs missing in $population ($printpop)"); }
					push (@report, "$shitlocus missing in all samples");
					$test++;
				}
			}
			else {
				foreach my $id2 (0..$indexofindex) {
					my $sampleposition = $sampleindex[$id2];
					my $checkrow = $wholerow[$sampleposition];
					my @matches = ($checkrow =~ /\.\/\./g);
					my $getmatch = @matches;
					$gralcount = $gralcount + $getmatch;
					$wholerow[$sampleposition]=~ s/\.\/\./$gralmiss/g;
					$numgral = $getmatch;
				}
			}
			$fixedline = join ("\t", @wholerow);
			$locirows[$snpcount] = $fixedline;
			$snpcount++;
		}
		print "Population $population processed. ";
		push (@report, "\nSummary population $population: $numgral general missing genotypes found (\"$gralmiss\" was input).\n$test SNPs were missing in the whole population, $printpop.\n");
		if ($popwise == 1) { print "$popwise locus was missing in $population$printpop"; }
		elsif ($popwise > 1) { print "$popwise loci were missing in $population$printpop"; }
		print "\n";
		$snpcount = 0;
	}
	push (@report, "\n\nA total of $gralcount general missing genotypes found, they were replaced with \"$gralmiss\".\n\n");
	print "\n$gralcount general missing values (found in few samples per population) were replaced by $gralmiss\n";
}



if ($pophandler eq "nd") { $pophandler = "0p"; }
if ($gralcount == 0) { $gralhandler = "_0m"; }
if ($gralhandler eq "nd") { $gralhandler = "_0m"; }





print "\nAll missing were handled gracefully, now checking the data";




my @locifin = ();
#@locifin = @locirows;
my $locinum = 0;
# Doublecheck after inputing in missing
$delete=0;
my $moda = "none";
foreach my $line (@locirows) {
	my $mosthigh = "none";
	my %freqgenot = ('./.' => 0, '0/0' => 0, '0/1' => 0, '1/1' => 0, "$misspop" => 0, "$gralmiss" => 0 );
	my @wholeline= split("\t", $line);  		#split columns as different elements of an array
	#print "\n\n$line\n\n";
	my $numcolumn = scalar @wholeline;
	my @SNP = @wholeline[$firstsample..$lastsample]; #save the information about the loci
	my @genotype = @SNP;  		#duplicate
	
	my @lociid = split (":", $wholeline[2]);  		#extract the genotype
	my $scaffold = $wholeline[0];
	$scaffold =~ s/^Scaffold([0-9]*?)$/$1/;
	my $lociname = "$scaffold" . "_$lociid[0]";
	
	$samplenum = scalar @SNP;
	my $last = $samplenum - 1;
	
	foreach my $count (0..$last) { 
		my @splitted = split(':', $genotype[$count]);  		#save only the genotype
		$genotype[$count] = $splitted[0];
	}
	
	foreach my $data (0..$last) {
		my $onegenotype = $genotype[$data];
		$freqgenot{$onegenotype} = $freqgenot{$onegenotype} +1;
	}
	
	delete $freqgenot{'./.'};
	for my $nm (sort {$freqgenot{$a} <=> $freqgenot{$b} || $a cmp $b } (keys %freqgenot) ) { $mosthigh = $nm; }
	my $highestfreq = $freqgenot{$mosthigh};
	
	#my $variability = %freqgenot;
	#print "$variability ";
	if($highestfreq >= $samplenum) {
		$delete++;
		my $lociinfo = "$lociname\t$scaffold\t$wholeline[1]\t$wholeline[2]\tmonomorphic after fixing missing";
		push(@deletedloci, $lociinfo);
	}
	else {
		push (@locifin, $line);
		my $lociinfo = "$lociname\t$scaffold\t$wholeline[1]\t$wholeline[2]";
		push(@locilist, $lociinfo);
		$goodloci++;
	}
	$locinum++;
}

my $bestloci = scalar @locifin;


my $outfile = "undefined";
my @outpath = split('/', $outname);
my $outlength = scalar @outpath;
my $outdir="undefined";

if ($tail eq "not defined") {
	if ($bestloci > 999) {
		my $locishort = $bestloci;
		$locishort=~ s/^([0-9]*?)[0-9][0-9][0-9]$/$1k/;
		$tail = "_$samplenum" . "x$locishort$gralhandler$pophandler";
	}
	else {
		my $locishort = $bestloci;
		$tail = "_$samplenum" . "x$locishort$gralhandler$pophandler";
	}
}




my $namekept=0;

if ($outname eq "not defined") {
	my @namefull = split('\.' , $filename);
	pop (@namefull);
	$namekept = join('.', @namefull);		   # replaced '\.' with '.' because it was saving the file as "populations\.snps_tail.vcf"
	$namekept=~ s/(.*?)_raw.*?/$1/;
	$outdir = $subdirpath;
	$outfile = "$outdir" . "/$namekept$tail" . ".vcf";
}
elsif ($outlength > 1) {
	$namekept = $outpath[-1];
	$namekept=~ s/(.*?)\.vcf/$1/;
	$namekept=~ s/(.*?)\.VCF/$1/;
	$namekept=~ s/(.*?)_raw.*?/$1/;
	pop (@outpath);
	$outdir = join ('/', @outpath);
	$outfile = "$outdir" . "/$namekept$tail" . ".vcf";
}
elsif ($outlength <= 1) {
	$namekept = $outname;
	$namekept=~ s/(.*?)\.vcf/$1/;
	$namekept=~ s/(.*?)\.VCF/$1/;
	$namekept=~ s/(.*?)_raw.*?/$1/;
	$outdir = $subdirpath;
	$outfile = "$outdir" . "/$namekept$tail" . ".vcf";
}



my $headrow = join ("\t", @headers);
my $finallenght = scalar @headers;
my $indexcol = $finallenght - 1;
my @finalsamples = @headers[$firstsample..$indexcol];
my @neatvcf = (@introrows, $headrow, @locifin);


# save all files

open my $OUTVCF, '>', $outfile or die "\nUnable to create or save \"$outfile\": $!\n";
foreach (@neatvcf) {print $OUTVCF "$_\n";} # Print each entry in our array to the file
close $OUTVCF; 

#replace empty quality metadata from genotypes
if ($quality == 0) {
	print " and inputing quality and probability to new genotypes\n";
	#system("sed -i 's/\t0\/0:.:.:.:./$homoref/g' $outfile");
	system("sed -i 's#\t0/0:.:.:.:.#$homoref#g' $outfile");
	#system("sed -i 's/\t0\/1:.:.:.:./$hetegen/g' $outfile");
	system("sed -i 's#\t0/1:.:.:.:.#$hetegen#g' $outfile");
	#system("sed -i 's/\t1\/1:.:.:.:./$homoalt/g' $outfile");
	system("sed -i 's#\t1/1:.:.:.:.#$homoalt#g' $outfile");
}
else { print "\n"; }

if ($delete > 0) {print "$delete of $locinum loci were monomorphic after handling missings.\n";}
elsif ($delete == 0) { print "Datafile seems alright.\n"; }

print "Final dataset has $finalsamp samples and $bestloci SNPs\n";
push (@report, "\nFinal dataset has $finalsamp samples and $bestloci SNPs\n");


my $shitlocinum = scalar (@deletedloci);
my $deletedpath = "$outdir" . "/bad_loci";
open my $OUTBAD, '>', $deletedpath or die "\nUnable to create or save \"$deletedpath\": $!\n";
if ($shitlocinum >= 1) { foreach (@deletedloci) {print $OUTBAD "$_\n";} }	  # Print each entry in our array to the file
elsif ($shitlocinum == 0) { print $OUTBAD "All loci were good, none deleted!\n" }
close $OUTBAD; 

my $shitsamplenum = scalar (@deletedsamples);
my $deletedsamplepath = "$outdir" . "/bad_samples";
open my $OUTBS, '>', $deletedsamplepath or die "\nUnable to create or save \"$deletedsamplepath\": $!\n";
if ($shitsamplenum >= 1) { foreach (@deletedsamples) {print $OUTBS "$_\n";} } 	  # Print each entry in our array to the file
else { print $OUTBS "All samples were good, non deleted!" }
close $OUTBS; 

my $locipath = "$outdir" . "/good_loci";
open my $OUTGOOD, '>', $locipath or die "\nUnable to create or save \"$locipath\": $!\n";
foreach (@locilist) {print $OUTGOOD "$_\n";} # Print each entry in our array to the file
close $OUTGOOD; 

my $samppath = "$outdir" . "/good_samples";
open my $SAMPGOOD, '>', $samppath or die "\nUnable to create or save \"$samppath\": $!\n";
foreach (@finalsamples) {print $SAMPGOOD "$_\n";} # Print each entry in our array to the file
close $SAMPGOOD; 

my $mappath = "$outdir" . "/fixed_popmap";
open my $MAPGOOD, '>', $mappath or die "\nUnable to create or save \"$mappath\": $!\n";
foreach (@finalsamples) {
	my $sampname = $_;
	my $popname = substr ($sampname, 0, $poplength);
	if ($popmap ne "No default") { $popname = $popref{$sampname}; }
	print $MAPGOOD "$sampname\t$popname\n";
} # Print each entry in our array to the file
close $MAPGOOD; 


print "\n$namekept$tail" . ".vcf and list of loci and samples saved succesfully.\n";

my $popfin = scalar @keepops;

my $tabletop = "undefined";
my $tableinfo = "undefined";


# compose summary table

if ($sumpath ne "no") {
	print "Saving summary table as $sumpath... ";
	if ($poplog ne "no" && -e $filelog) {
		
		
		my $pattern3 = "private alleles";
		my $pattern4 = "Percent samples limit";
		my $pattern5 = "Locus Population limit";
		my $pattern6 = "Minor allele frequency";
		my $pattern7 = "defaultgrp";
		my $pattern9 = "Maximum observed heterozygosity cutoff:";
		my $pattern10 = "Working on [0-9]* samples";
		
		# Extract data from populations log file
		open my $output1, '<', $filelog or die "unable to open file '$filelog: $!";
		push (@report, "\n\n$logname file read to create the summary table.");
		my @privatealleles=();
		my $popini=0;
		my $r =0;
		my $p =0;
		my $m = 0;
		my $mh = 0;
		my $iniN = 0;
		my $date = 0;
		my $time = 0;
		my $outloop =0;

		while(my $line2=<$output1>) {
			if ($outloop ==0) {
				my $rundate = $line2;
				my @infoline= split(' ', $rundate);
				$date = $infoline[-2];
				$time = $infoline[-1];
			}
			if ($line2 =~ $pattern3) {
				my $private = $line2;
				my @privateline= split(' ', $private);
				my $pa = $privateline[-1];
				push (@privatealleles, $pa);
			}
			elsif ($line2 =~ $pattern4) {
				my $rflag = $line2;
				my @percentsamples= split(' ', $rflag);
				$r = $percentsamples[-1];
			}
			elsif ($line2 =~ $pattern5) {
				my $pflag = $line2;
				my @percentpop= split(' ', $pflag);
				$p = $percentpop[-1];
			}
			elsif ($line2 =~ $pattern10) {
				my $popmapinfo = $line2;
				my @popmapline= split(' ', $popmapinfo);
				$iniN = $popmapline[2];
			}
			elsif ($line2 =~ $pattern6) {
				my $mflag = $line2;
				my @minorallelefreq= split(' ', $mflag);
				$m = $minorallelefreq[-1];
			}
			elsif ($line2 =~ $pattern7) {
				my $poplist = $line2;
				my @populations= split(',', $poplist);
				$popini = scalar @populations;
			}
			elsif ($line2 =~ $pattern9) {
				my $heteroz = $line2;
				my @minmaxhet= split(' ', $heteroz);
				$mh = $minmaxhet[-1];
			}
			$outloop++;
		}
		close ($output1);
		
		
		my @allelenum = sort { $a <=> $b } @privatealleles;
		my $minim = $allelenum[0];
		my $maxim = $allelenum[-1];
		
		
		$tabletop = "Date\tTime\tPop_ini\tSample_ini\tm\tr\tp\tsnps_ini\tmax_het\tPamin\tPamax\tgral_miss\tmiss_popwise\tpop\tsamples\tsnps";
		$tableinfo = "$date\t$time\t$popini\t$iniN\t$m\t$r\t$p\t$iniloci\t$mh\t$minim\t$maxim\t$gralcount\t$popwisecount\t$popfin\t$finalsamp\t$bestloci";	
	}
	else {
		print "\nNo populations (Stacks) log file $filelog found.\n";
		push (@report, "\n\nNo ref_map or populations log file found. No \"$filelog\".\n\n");
		
		##############################################################################################  MAY NOT WORK TILL I GET HELP WITH THE PERL MODULE DATETIME
		
		my ($sec2,$min2,$hour2,$mday2,$mon2,$year2,$wday2,$yday2,$isdst2) = localtime(time);
		$mon2 = $mon2 + 1;
		$year2 = $year2 + 1900;
		my $time = "$year2-$mon2-$mday2 $hour2:$min2:$sec2";
		
		my $popini = $popnum;
		my $iniN = $sampnum;
		
		$tabletop = "Filepath\tDate_time\tPop_ini\tSample_ini\tsnps_ini\tgral_miss\tmiss_popwise\tpop\tsamples\tsnps";
		$tableinfo = "$filepath\t$time\t$popini\t$iniN\t$iniloci\t$gralcount\t$popwisecount\t$popfin\t$finalsamp\t$bestloci";	
		
	}
	
	my @pathsum = split('/', $sumpath);
	$pathlength = scalar @pathsum;
	#create directory if needed
	if ($pathlength > 1) {
		pop @pathsum;
		my $subdir = join ('/', @pathsum);
		unless(-e $subdir or mkdir $subdir) {die "Unable to create output directory $subdir\n"};
	}
	
	if (-e $sumpath) {
		open my $SUM, '>>', $sumpath or die "\nUnable to write or save \"$sumpath\": $!\n";
		print $SUM "$tableinfo\n";
		push (@report, "$sumpath found, will not overwrite headers. Headers used for this session:\n$tabletop\n\n");
		close ($SUM);
	}
		else {
		open my $SUM, '>', $sumpath or die "\nUnable to create or save \"$sumpath\": $!\n";
		print $SUM "$tabletop\n";
		print $SUM "$tableinfo\n";
		close ($SUM);
	}
	print "done!\n";
}

my $logfile = $version;
$logfile =~ s/(.*?)\.pl/$1\.log/;
my $warnings = scalar @report;


##############################################################################################  MAY NOT WORK TILL I GET HELP WITH THE PERL MODULE DATETIME

my ($sec3,$min3,$hour3,$mday3,$mon3,$year3,$wday3,$yday3,$isdst3) = localtime(time);
$mon3 = $mon3 + 1;
$year3 = $year3 + 1900;
my $ending = "$year3-$mon3-$mday3 $hour3:$min3:$sec3";


push (@report, "\n\nDone. Check bad_samples and bad_loci to see a list of the elements deleted and why were they deleted.\nCheck good_samples and good_loci for a list fo elements kept.\nfixed_popmap has the individual-population information for the remaining samples\n\nFinished $ending\n");

if ($warnings > 0) {
	my $reportway = "$outdir/$logfile";
	print "\nSaving $reportway\n\n";
	open my $LOG, '>', $reportway or die "\nUnable to create or save \"$reportway\": $!\n";
	foreach (@report) {print $LOG "$_\n";} # Print each entry in our array to the file
	close $LOG; 
}

print "\nALL DONE. $version finished!\n\n";
