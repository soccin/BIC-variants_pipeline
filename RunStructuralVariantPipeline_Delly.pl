#!/usr/bin/perl -w
##########StructuralVaraintFinder.pl########
#Author: Ronak Shah
#Date: 07/01/2013
#LastModified: 08/22/2013
#Version: 1.0
#Description: Get In the data and run the mapping and structural variant finding
#process.
#############
###Log:
##07/01/2013
#v1.0
use strict;
use Getopt::Long;
use IO::File;
use List::Util qw(max sum);
use Tie::IxHash;
use Vcf;
use File::Basename;
###use FindBin::Real qw(Bin);
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Schedule;
use Cluster;

my (
    $sampleFile,
    $stdNormal,
    $titleFile,
    $fof,
    $poolName,
    $projectName,
    $process,
    $datadir,
    $outdir,
    $mvFiles,
    $bwa,
    $baitIntervalFile,
    $targetIntervalFile,
    $targetIntervalBedFile,
    $refFile,
    $TMPDIR,
    $PICARD,
    $JAVA,
    $samtools,
    $fastqSource,
    $GATK,
    $PYTHON,
    $PYTHONPATH,
    $prog,
    $RefGeneFile,
    $RepeatRegionFile,
    $CancerCensusFile,
    $queue,
    $PERL,
    $DELLY,
    $QSUB,
    $standardNormalList,
    $MAPQ,
    $BASEQ,
    $barcodeFile,
    $adaptorFile,
    $TrimGalore,
    $nprocessors,
    $ZCAT,
    $GZIP,
    $FilterSV,
    $AnnotateSV,
    $DGvFile,
    $OverallSupportingReads,
    $OverallSupportingReadsHotspot,
    $SampleTumorSupportingReads,
    $SampleTumorSupportingReadsHotspot,
    $SampleNormalSupportingReads,
    $SampleNormalSupportingReadsHotspot,
    $OverallSupportingSplitReads,
    $OverallSupportingSplitReadsHotspot,
    $LengthOfSV,
    $OverallMapq,
    $OverallMapqHotspot,
    $SampleTumorGenotypeQualityFilter,
    $SampleTumorGenotypeQualityFilterHotspot,
    $SampleNormalGenotypeQualityFilter,
    $SampleNormalGenotypeQualityFilterHotspot,
    $HotspotFile,
    $DistanceBtwTumorNormalCTX,
    $dRANGER,
    $MCR,
    $pair_file,
    $bam_list,
    $scheduler,
    $priority_project,
    $priority_group
);

# This is variable is the path to the bin folder
my $binPath    = $Bin;
my $hg19DataPath = $binPath . "/data/";
my $configurationPath = $binPath . "/strvar/";   #path to configuration Files
my $HG19MAT = $binPath . "/data/RefSeq_hg19.mat";
my $config_file = $configurationPath . "template_dmp_sv.conf"; #default config file
my $ExcludeRegions = $hg19DataPath . "human.hg19.excl.tsv";

#--This variable holds the current time
my $now = time;

if (@ARGV < 1 or !GetOptions (
	'pre=s'                             => \$poolName,
	'config|c:s'                        => \$config_file,       
	'outputDirectory|o:s'               => \$outdir,
	'pair=s'                            => \$pair_file,
	'bam_list=s'                        => \$bam_list,
	'scheduler=s'                        => \$scheduler,
	'priority_project=s' => \$priority_project,
	'priority_group=s' => \$priority_group
))
	{
		Usage();
	}

#Check for Pool Name
if(!$poolName)
{
    die "Please specify the project name by -pre\n";  
}

if ($config_file) {
	print "The configration file in use is $config_file.\n";
}
else {
	print "Please specigy configuration file by -config\n";
	Usage();
	exit;
}

if ( !$outdir ) {
	print
"Please enter the directory where to write the data while/after being processed.See Usage\n";
	Usage();
	exit;
}
else {
	print "Results will be written in $outdir\n";
	`/bin/mkdir -p $outdir`;
	die "[ERROR]: Fail to create outdir: $outdir\n" if(!-d $outdir);
}

my %NormalTumorPair = ();
if(!$pair_file)
{
    die "Please specify the pairing file by using -pair PAIR_FILE\n";
}
else
{
    ReadNormalTumorPair();
}

my %SampleBamHash = ();
my %BamSampleHash = ();
if(!$bam_list)
{
    die "Please specify the bam list file by using -bam_list BAM_LIST_FILE\n";
}
else
{
    ReadBamList();
}


#Get Configration File details
my ($Version) = &GetConfiguration( $config_file, $outdir );
my $PrintConfigFile = "RunConfigration_StrVar.txt";
open( VERSION, ">", "$outdir/$PrintConfigFile" )
  || die "Cannot open $outdir/$PrintConfigFile, Reason:$!\n";

#Prin Version of tools and Files used
print VERSION "Tools|Files\tVersion\n";
while ( my ( $tools_files, $version ) = each(%$Version) ) {
	print VERSION "$tools_files\t$version\n";
}
close(VERSION);

#use default value to override config file
$TMPDIR = $outdir . "/tmp/";
$HotspotFile = $hg19DataPath . "v3clin_hg19_structuralvariants_geneInterval.txt";
$CancerCensusFile = $hg19DataPath . "cancer_gene_census.tsv";
$DGvFile = $hg19DataPath . "DGv_Annotation.tsv";
$RepeatRegionFile = $hg19DataPath . "repeatRegion.tsv";

`/bin/mkdir -p $TMPDIR`;
die "[ERROR]: Fail to create TMPDIR: $TMPDIR\n" if(!-d $TMPDIR);


# if(!$sampleFile)
# {
#     print "Sample sheet is not given in the configuration file. Looking into the raw data directory.\n";
#     $sampleFile = `ls $datadir/SampleSheet.csv`;
#     chomp($sampleFile);
    
#     if (!$sampleFile) {
#        print "Sample sheet could not be located in the raw data directory, pipeline requires a sample sheet to run. Please see usage.\n";
#         Usage();
#         exit(1);       
#     }
# }

# if (!$titleFile) {
#     print "Title file is not given in the configuration file. Looking into the raw data directory.\n";
#     $titleFile = `ls $datadir/title_file.txt`;
#     chomp($titleFile);
    
#     if (!$titleFile) {
#         print "Title file could not be located in the raw data directory, pipeline requires the title file to run. Please see usage.\n";
#         Usage();
#         exit(1);       
#     }
# }
# #Read Sample File
# my (
# 	$fcId,        $lane,    $sampleId, $sampleRef, $index,
# 	$description, $control, $recipe,   $operator,  $sampleProject
# ) = &ReadSampleFile( $sampleFile, $projectName, $outdir );

# #Read Title File
# my (
# 	$barcode,      $pool,      $titleSampleId, $collabId,
# 	$patientId,    $class,     $sampleType,    $inputNg,
# 	$libraryYeild, $poolInput, $baitVersion
# ) = &ReadTitleFile( $titleFile, $outdir );

#Check for Project Name
#if(!$projectName)
#{
   	#print "Project name is not given in the configuration file, pipeline will use the project name stated in the SampleSheet.csv file.\n";
    #$projectName = @$sampleProject[0];
    #if (!$projectName) {
#       print "Can not read project name in the configuration file.\n";
#        exit(1);
        
    #}
    
#}

if ( !$mvFiles ) {
	print "Folders will be created and Files will be moved.\n";
	$mvFiles = 1;
}
else {
	if ( $mvFiles == 1 ) {
		print "Folders will be created and Files will be moved.\n";
	}
	if ( $mvFiles == 2 ) {
		print "Folders will not be created and Files will not be moved.\n";
	}
}
if ( !$stdNormal ) {
	$stdNormal = "NULL";
}
else {
	print "Starndard Normal is given: $stdNormal\n";
}

if ( !$fastqSource ) {
	print "Assume fastq files came from GCL.\n";
	$fastqSource = "GCL";
}
else {
	if ( $fastqSource ne "GCL" && $fastqSource ne "DMP" ) {
		print "Please indicate fastqSource. See Usage\n";
		Usage();
		exit;
	}
}
if ($barcodeFile) {
	print "The barcode file in use is $barcodeFile.\n";
}
else {
	print "Please enter the barcode file.See Usage\n";
	Usage();
	exit;
}
if ($adaptorFile) {
	print "The barcode file in use is $adaptorFile.\n";
}
else {
	print "Please enter the adaptor file.See Usage\n";
	Usage();
	exit;
}

if ( !$TMPDIR ) {
	print "Path to temporary directory is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "TMPDIR=$TMPDIR\n";
}
if ( !$PICARD ) {
	print "Path to samtools executables is not given.See Usage\n";
	Usage();
	exit;
}
else {
	print "SAMTOOLS=$samtools\n";
}
if ( !$PICARD ) {
	print "Path to Picard's executables is not given, See Usage.\n";
	Usage();
	exit;
}
else {
	print "PICARD=$PICARD\n";
}
if ( !$bwa ) {
	print "Path to BWA MEM executables is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "BWA=$bwa\n";
}
if ( !$baitIntervalFile ) {
	print "Bait Interval file is not given. See Usage.\n";
	Usage();
	exit;
}
else {
	print "BAIT_INTERVAL=$baitIntervalFile\n";
}
if ( !$targetIntervalFile ) {
	print "Target Interval file is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "Target_INTERVAL=$targetIntervalFile\n";
}
if ( !$refFile ) {
	print "Reference Sequence file is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "Reference_Sequence=$refFile\n";
}
###if ( !$ExcludeRegions ) {
###	print "ExcludeRegions file is not given. See Usage\n";
###	Usage();
###	exit;
###}
###else {
###	print "ExcludeRegions=$ExcludeRegions\n";
###}
if ( !$HotspotFile ) {
	print "Hotspot location BED4 file is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "HotspotFile=$HotspotFile\n";
}
if ( !$JAVA ) {
	print "Path to java executables is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "JAVA=$JAVA\n";
}

#Check if path to FilterSV script is given
if ( !$FilterSV ) {
	$FilterSV = $binPath . "/strvar/FilterSV.pl";
	print
	  "FilterSV.pl script path is not given, default location will be used.\n";
}
else {
	print "Filter SV script location: $FilterSV.\n";
}
#Check if path to FilterSV script is given
if ( !$AnnotateSV ) {
	$AnnotateSV = $binPath . "/strvar/AnnotateSVs.py";
	print
	  "AnnotateSVs.py script path is not given, default location will be used.\n";
}
else {
	print "Annotate SV script location: $AnnotateSV.\n";
}
if ( !$DELLY ) {
	print "Path to Delly is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "Delly=$DELLY\n";
}

#overallSupportingReadsPE
if ( !$OverallSupportingReads ) {
	print "OverallSupportingReads not given default will be used.\n";
	$OverallSupportingReads = 5;
	print "OverallSupportingReads:$OverallSupportingReads\n";
}
else {
	print "OverallSupportingReadsPE:$OverallSupportingReads.\n";
}

#sampleTumorSupportingReadsPE
if ( !$SampleTumorSupportingReads ) {
	print "SampleTumorSupportingReadsPE not given default will be used.\n";
	$SampleTumorSupportingReads = 2;
	print "SampleTumorSupportingReadsPE:$SampleTumorSupportingReads\n";
}
else {
	print "SampleTumorSupportingReadsPE:$SampleTumorSupportingReads.\n";
}

#sampleTumorSupportingReadsPE_Hotspot
if ( !$SampleTumorSupportingReadsHotspot ) {
	print "SampleTumorSupportingReadsHotspot not given default will be used.\n";
	$SampleTumorSupportingReadsHotspot = 1;
	print
	  "SampleTumorSupportingReadsHotspot:$SampleTumorSupportingReadsHotspot\n";
}
else {
	print
	  "SampleTumorSupportingReadsHotspot:$SampleTumorSupportingReadsHotspot.\n";
}

#sampleNormalSupportingReadsPE
if ( !$SampleNormalSupportingReads ) {
	print "SampleNormalSupportingReadsPE not given default will be used.\n";
	$SampleNormalSupportingReads = 2;
	print "SampleNormalSupportingReadsPE:$SampleNormalSupportingReads\n";
}
else {
	print "SampleNormalSupportingReadsPE:$SampleNormalSupportingReads.\n";
}

#sampleNormalSupportingReadsPE_Hotspot
if ( !$SampleNormalSupportingReadsHotspot ) {
	print
	  "SampleNormalSupportingReadsHotspot not given default will be used.\n";
	$SampleNormalSupportingReadsHotspot = 3;
	print
"SampleNormalSupportingReadsHotspot:$SampleNormalSupportingReadsHotspot\n";
}
else {
	print
"SampleNormalSupportingReadsHotspot:$SampleNormalSupportingReadsHotspot.\n";
}

#overallSupportingReadsSR
if ( !$OverallSupportingSplitReads ) {
	print "OverallSupportingSplitReads not given default will be used.\n";
	$OverallSupportingSplitReads = 0;
	print "OverallSupportingReadsSR:$OverallSupportingSplitReads\n";
}
else {
	print "OverallSupportingReadsSR:$OverallSupportingSplitReads.\n";
}

#overallSupportingReadsPE_Hotspot
if ( !$OverallSupportingReadsHotspot ) {
	print "OverallSupportingReadsHotspot not given default will be used.\n";
	$OverallSupportingReadsHotspot = 3;
	print "OverallSupportingReadsPE_Hotspot:$OverallSupportingReadsHotspot\n";
}
else {
	print "OverallSupportingReadsPE_Hotspot:$OverallSupportingReadsHotspot.\n";
}

#overallSupportingReadsSR_Hotspot
if ( !$OverallSupportingSplitReadsHotspot ) {
	print
	  "OverallSupportingSplitReadsHotspot not given default will be used.\n";
	$OverallSupportingSplitReadsHotspot = 0;
	print
"OverallSupportingSplitReadsHotspott:$OverallSupportingSplitReadsHotspot\n";
}
else {
	print
"OverallSupportingSplitReadsHotspot:$OverallSupportingSplitReadsHotspot.\n";
}

#lengthOfSV
if ( !$LengthOfSV ) {
	print "LengthOfSV not given default will be used.\n";
	$LengthOfSV = 300;
	print "LengthOfSV:$LengthOfSV\n";
}
else {
	print "LengthOfSV:$LengthOfSV.\n";
}

#overallMAPQ
if ( !$OverallMapq ) {
	print "OverallMapq not given default will be used.\n";
	$OverallMapq = 10;
	print "OverallMapq:$OverallMapq\n";
}
else {
	print "OverallMapq:$OverallMapq.\n";
}

#overallMAPQ_Hotspot
if ( !$OverallMapqHotspot ) {
	print "OverallMapqHotspot not given default will be used.\n";
	$OverallMapqHotspot = 5;
	print "OverallMapqHotspot:$OverallMapqHotspot\n";
}
else {
	print "OverallMapqHotspot:$OverallMapqHotspot.\n";
}

#sampleTumorGenotypeQualityFilter
if ( !$SampleTumorGenotypeQualityFilter ) {
	print "SampleTumorGenotypeQualityFilter not given default will be used.\n";
	$SampleTumorGenotypeQualityFilter = 15;
	print
	  "SampleTumorGenotypeQualityFilter = $SampleTumorGenotypeQualityFilter\n";
}
else {
	print
	  "SampleTumorGenotypeQualityFilter = $SampleTumorGenotypeQualityFilter\n";
}

#sampleTumorGenotypeQualityFilterHotspot
if ( !$SampleTumorGenotypeQualityFilterHotspot ) {
	print
"SampleTumorGenotypeQualityFilterHotspot not given default will be used.\n";
	$SampleTumorGenotypeQualityFilterHotspot = 5;
	print
"SampleTumorGenotypeQualityFilterHotspot = $SampleTumorGenotypeQualityFilterHotspot\n";
}
else {
	print
"SampleTumorGenotypeQualityFilterHotspot = $SampleTumorGenotypeQualityFilterHotspot\n";
}

#sampleNormalGenotypeQualityFilter
if ( !$SampleNormalGenotypeQualityFilter ) {
	print "SampleNormalGenotypeQualityFilter not given default will be used.\n";
	$SampleNormalGenotypeQualityFilter = 15;
	print
"SampleNormalGenotypeQualityFilter = $SampleNormalGenotypeQualityFilter\n";
}
else {
	print
"SampleNormalGenotypeQualityFilter = $SampleNormalGenotypeQualityFilter\n";
}

#sampleGenotypeQualityFilterHotspot
if ( !$SampleNormalGenotypeQualityFilterHotspot ) {
	print
"SampleNormalGenotypeQualityFilterHotspot not given default will be used.\n";
	$SampleNormalGenotypeQualityFilterHotspot = 20;
	print
"SampleNormalGenotypeQualityFilterHotspot = $SampleNormalGenotypeQualityFilterHotspot\n";
}
else {
	print
"SampleNormalGenotypeQualityFilterHotspot = $SampleNormalGenotypeQualityFilterHotspot\n";
}

#DistanceBtwTumorNormalCTX
if ( !$DistanceBtwTumorNormalCTX ) {
	print "DistanceBtwTumorNormalCTX not given default will be used.\n";
	$DistanceBtwTumorNormalCTX = 5;
	print "DistanceBtwTumorNormalCTX = $DistanceBtwTumorNormalCTX\n";
}
else {
	print "SDistanceBtwTumorNormalCTX = $DistanceBtwTumorNormalCTX\n";
}


#tie( my %classPerBarcode, 'Tie::IxHash' );
#for ( my $i = 0 ; $i < scalar(@$barcode) ; $i++ ) {

	#print "$$barcode[$i] => $$class[$i]\n";
#	$classPerBarcode{ $$barcode[$i] . "_" . $$titleSampleId[$i] } = $$class[$i];
#}

# Make a dummy wait file to run after process
&MakeCSH($outdir);
#Split porces to know how many to run
my @allProcess = split( ",", $process );

#Check for raw data directory
# if ( $allProcess[0] == 1 && !$datadir ) {
# 	print
# "Please enter the directory that contains the data to be processed.See Usage\n";
# 	Usage();
# 	exit(1);
# }
# elsif ( $allProcess[0] == 1 && $datadir ) {
# 	print "The raw data directory given is $datadir\n";
# }
# else {
# 	print
# "Process 1 is not selected so I am not checking to see if the raw data directory is given or not.\n";
# }
my $allProcessList  = join( ",", @allProcess );
my $numberOfProcess = scalar(@allProcess);
my $processCount    = 0;
my $parseFilenames;
while ( $processCount < $numberOfProcess ) {

	#print "@allProcess\n";
	my $runProcess = shift(@allProcess);

	#print "$runProcess\n";
	($parseFilenames) = &Select( $runProcess, $parseFilenames );
	$processCount++;
}
# if ($fof) {
# 	print "Removing Softlinked files form $outdir\n";
# 	my @softlinks = ();
# 	eval {
# 		@softlinks =
# 		  `find $outdir -maxdepth 1 -type l -print0 | xargs -0 ls -d`;
# 	};
# 	if ($@) { print "Cannot find all softlinks in $outdir. Error:$@\n"; }
# 	foreach my $sfLink (@softlinks) {
# 		eval { `unlink $sfLink`; };
# 		if ($@) { print "Cannot unlink $sfLink. Error:$@\n"; }
# 	}
# }

#--Clean up intermediate files

CleanUpFiles($outdir . "/StrVarAnalysis", $outdir . "/unfiltered");


#--Calculate total runtime (current time minus start time)
$now = time - $now;

#--Print runtime
printf(
	"\n\nTotal running time: %02d:%02d:%02d\n\n",
	int( $now / 3600 ),
	int( ( $now % 3600 ) / 60 ),
	int( $now % 60 )
);
print "\n!!!!Done, Thanks for using the script!!!\n";
exit;
#####################################
#####################################
#Run only single process at a time.
sub Select {
	my ( $process, $parseFilenames ) = @_;
	my $SVoutput;
	if ( $process == 1 ) {
		#$parseFilenames = &MergeDataFromDirectory();
	}
	elsif ( $process == 2 ) {
		#$parseFilenames = &DoMapping($parseFilenames);
	}
	elsif ( $process == 3 ) {
		#$parseFilenames = &CalcHsMetrics($parseFilenames);
	}
	elsif ( $process == 4 ) {
		$parseFilenames = &CallStructuralVariants($parseFilenames);
	}
	elsif ( $process == 5 ) {
		$parseFilenames = &FilterStructuralVariants($parseFilenames);
	}
	elsif ( $process == 6 ) {
		( $parseFilenames, $SVoutput ) =
		  &AnnotateStructuralVariants($parseFilenames);
	}
	elsif (( $process eq "all" )
		or ( $process eq "ALL" )
		or ( $process eq "All" ) )
	{
		#$parseFilenames = &MergeDataFromDirectory();
		#$parseFilenames = &DoMapping($parseFilenames);
		#$parseFilenames = &CalcHsMetrics($parseFilenames);
		$parseFilenames = &CallStructuralVariants($parseFilenames);
		$parseFilenames = &FilterStructuralVariants($parseFilenames);
		$parseFilenames = &AnnotateStructuralVariants($parseFilenames);
	}
	else {
		print
"Select:The process number entered ($process) does not exist, See Usage.\n";
		exit(1);
	}
	return ($parseFilenames);
}
#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);
	print "\nUsage : StructuralVariantFinder.pl [options]
        [--config|c                        S Path to configration file(required;Template:/home/shahr2/Scripts/All/StrVarConfig_(script_version)_(bait_version.txt)]
        [--dataDirectory|d                 S Path where all the files to be processed are located (required)]
        [--outputDirectory|o               S Path where all the output files will be written (required)]
        Inside Config File:
        >Locations:
        barcodeFile                        S tab-delimited barcode file describing name of the barcode(Required)
        adaptorFile                        S tab-delimited adaptor file describing sequence of the adaptor for corresponding barcode(Required)
        BWA                                S Path to bwa mem program. (Required)
        PICARD                             S Path to picard tools (Required)
       	DELLY							   S Path to delly executables (Required)
        SAMTOOLS                           S Path to samtools source folder.(Required)
        BaitInterval                       S Path to baits interval file.(Required)
        TargetInterval                     S Path to targets interval file. (Required)
        Reference                          S Path to genome reference file. (Required)
        TMPDIR                             S Path to temporary directory. (Required)
        JAVA                               S Path to java executables. (Required)
        ExcludeRegions					   S Path to bed file indicating regions to exclude(Required)
        >Parmeters
        SampleFile                         S csv file describing details about the sample (required and submit with full path)
        TitleFile                          S tab-delimited title file for the samples (required and submit with full path)
        stdNormal                          S file to be used as standard normal #full path of bam file, *.bai file to be located in same folder)
        ListOfStandardNormals              S Name of the normal files with there path (fof:Paired files;one per line;one after another) (Required)
        projectName                        S Name of the project(required,e.g:Colons).
        poolName                           S Name of the pool(required,e.g:Colon5P1).]
        Process                            S 1 => MergingFastq 2 => Run Mapping. 3 => Run HsMetrics. 4 => Call Structural Variants. 5 => Filter Structural Varaints 6 => Annotate Structural Variants\"1,2,3,4,5,6\" for all or in valid combination.]
        Datadir                            S Path where all the files to be processed are located. (required)
        Outdir                             S Path where all the output files will be written. (required)
        moveFiles                          I 2 => Skip Making Folders & Moving Files. 1 => Make Folders and Move Files. (default:1,optional)
       	NumberOfProcessors                 I Number of processors to use for analysis (default:1,optional)
        FOF                                S Name of the files with there path. (fof:Paired files;one per line;one after another) (optional)
        >Versions
        This would contain the versions for all things used in the pipeline.
        >Parmeters
        This would contain the parameters for the Pipleine
	\n";
	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}
###################################################
###################################################
#--GET CONFIGRATION DETAIL
sub GetConfiguration {
	my ($config_file) = @_;
	my @data = ();
	tie( my %config,     'Tie::IxHash' );
	tie( my %location,   'Tie::IxHash' );
	tie( my %version,    'Tie::IxHash' );
	tie( my %parameters, 'Tie::IxHash' );
	print "Reading the configuration file\n";

	# change the default Input Record Separator.
	$/ = ">";
	open( CONFIG, "$config_file" )
	  or die "Cannot open $config_file. Error: $!\n";
	my $junk = <CONFIG>;
	while (<CONFIG>) {
		next if ( $_ =~ /^#/ );
		chomp($_);
		my ( $defLine, @configLines ) = split /\n/, $_;
		if ( $defLine =~ /Locations/ ) {
			foreach my $config (@configLines) {
				next if ( $config =~ /^#/ );
				@data = split( "=", $config, 2 );
				$data[0] =~ s/\s//g;
				$data[1] =~ s/\s//g;
				$location{ $data[0] } = $data[1];
			}
		}
		if ( $defLine =~ /Parameters/ ) {
			foreach my $config (@configLines) {
				next if ( $config =~ /^#/ );
				@data = split( "=", $config, 2 );
				$data[0] =~ s/\s//g;
				$data[1] =~ s/\s//g;
				$parameters{ $data[0] } = $data[1];
			}
		}
		if ( $defLine =~ /Versions/ ) {
			foreach my $config (@configLines) {
				next if ( $config =~ /^#/ );
				@data = split( "=", $config, 2 );
				$data[0] =~ s/\s//g;
				$data[1] =~ s/\s//g;
				$version{ $data[0] } = $data[1];
			}
		}
	}
	close(CONFIG);
	print "Completed reading the configuration file\n";

	# Change back the Input Record Separator to default.
	$/ = "\n";
	###Set Variables
	eval {
		$TMPDIR             = $location{"TMPDIR"};
		$PERL               = $location{"PERL"};
		$JAVA               = $location{"JAVA"};
		$GATK               = $location{"GATK"};
		$refFile            = $location{"Reference"};
		$PICARD             = $location{"PICARD"};
		$baitIntervalFile   = $location{"BaitInterval"};
		$targetIntervalFile = $location{"TargetInterval"};
		$bwa                = $location{"BWA"};
		$PYTHON				= $location{"PYTHON"};
		$PYTHONPATH			= $location{"PYTHONPATH"};
		$samtools           = $location{"SAMTOOLS"};
		$DELLY              = $location{"DELLY"};
		$barcodeFile        = $location{"barcodeFile"};
		$RepeatRegionFile   = $location{"RepeatRegionAnnotation"};
		$DGvFile     		= $location{"DGvAnnotations"};
		$CancerCensusFile = $location{"CosmicCensus"};
		#$RepeatRegionFile   = $location{"RepeatRegionAnnotation"};
		$RefGeneFile        = $location{"RefGeneFile"};
		$barcodeFile        = $location{"BarcodeKey"};
		$adaptorFile        = $location{"AdaptorKey"};
		$QSUB               = $location{"QSUB"};
		$TrimGalore         = $location{"TrimGalore"};
		$ZCAT               = $location{"ZCAT"};
		$GZIP               = $location{"GZIP"};
		$FilterSV           = $location{"FilterSV"};
		$AnnotateSV			= $location{"AnnotateSV"};
		$HotspotFile        = $location{"HotspotFile"};
		$dRANGER            = $location{"dRANGER"};
		$MCR                = $location{"MCR"};
		$queue              = $parameters{"SGE_QUEUE"};
		$fastqSource        = $parameters{"fastqSource"};
		$sampleFile         = $parameters{"SampleFile"};
		$titleFile          = $parameters{"TitleFile"};
		$standardNormalList = $parameters{"ListOfStandardNoramlsForGenotyping"};
		#$outdir             = $parameters{"Outdir"};
		#$datadir            = $parameters{"Datadir"};
		$stdNormal          = $parameters{"stdNormal"};
		$fof                = $parameters{"FOF"};
		$MAPQ               = $parameters{"MAPQ"};
		$BASEQ              = $parameters{"BASEQ"};
		#$poolName           = $parameters{"poolName"};
		#$projectName        = $parameters{"projectName"};
		$mvFiles            = $parameters{"moveFiles"};
		$process            = $parameters{"Process"};
		$prog               = $parameters{"Program"};
		$nprocessors        = $parameters{"NumberOfProcessors"};
		$OverallSupportingReads = $parameters{"OverallSupportingReads"};
		$OverallSupportingReadsHotspot =
		  $parameters{"OverallSupportingReadsHotspot"};
		$SampleTumorSupportingReads = $parameters{"SampleTumorSupportingReads"};
		$SampleTumorSupportingReadsHotspot =
		  $parameters{"SampleTumorSupportingReadsHotspot"};
		$SampleNormalSupportingReads =
		  $parameters{"SampleNormalSupportingReads"};
		$SampleNormalSupportingReadsHotspot =
		  $parameters{"SampleNormalSupportingReadsHotspot"};
		$OverallSupportingSplitReads =
		  $parameters{"OverallSupportingSplitReads"};
		$OverallSupportingSplitReadsHotspot =
		  $parameters{"OverallSupportingSplitReadsHotspot"};
		$LengthOfSV         = $parameters{"LengthOfSV"};
		$OverallMapq        = $parameters{"OverallMapq"};
		$OverallMapqHotspot = $parameters{"OverallMapqHotspot"};
		$SampleTumorGenotypeQualityFilter =
		  $parameters{"SampleTumorGenotypeQualityFilter"};
		$SampleTumorGenotypeQualityFilterHotspot =
		  $parameters{"SampleTumorGenotypeQualityFilterHotspot"};
		$SampleNormalGenotypeQualityFilter =
		  $parameters{"SampleNormalGenotypeQualityFilter"};
		$SampleNormalGenotypeQualityFilterHotspot =
		  $parameters{"SampleNormalGenotypeQualityFilterHotspot"};
	};
	if ($@) {
		print "Did not find a variable in configuration file.Error: $@\n";
		exit(1);
	}
	return ( \%version );
}
###################################################
###################################################
#--Make Notification file
sub MakeCSH {
	my ($outdir) = @_;
	my $filename = $outdir . "/Notify.csh";
	if ( !-e $filename ) {
		my $ntmp = new IO::File(">$filename");
		print $ntmp "#!/bin/csh\n";
		print $ntmp "#Notification File\n";
		print $ntmp "echo", " This is Done", "\n";
		$ntmp->close();
		`chmod +x $filename`;
	}
	else {
		print "$filename exists and wont be created.\n";
	}
	return;
}
###################################################
###################################################
#--Waiting for the process to finish
sub WaitToFinish {
	my ( $outdir, @waitfilenames ) = @_;
	print "Waiting for the Process to finish...\n";
	foreach my $wfile (@waitfilenames) {
		next if ( $wfile eq "NULL" );
		wait while ( !-e "$outdir/$wfile" );

		#print "$outdir/$wfile\n";
		while ( -e "$outdir/$wfile" ) {

			#print "$outdir/$wfile\n";
			open( FH, "<", "$outdir/$wfile" );
			while (<FH>) {
				if ( $_ =~ /This is Done/ig ) {

					#print "\nFinished: $wfile\n";
					last;
				}
				else {
					wait;
				}
			}
			last;
		}
		close(FH);
	}
	foreach my $wfile (@waitfilenames) {
		next if ( $wfile eq "NULL" );
		`rm $outdir/$wfile`;
	}
	return;
}
###################################################
###################################################
#--Make array of file of files list from the outdir
# sub GetNames {
# 	my ( $fof, $outdir ) = @_;
# 	my (@filenames) = ();
# 	open( FOF, "$outdir/$fof" )
# 	  || die "Cannot open ListFile: \"$outdir/$fof\"\n";
# 	while (<FOF>) {
# 		$_ =~ s/\s//g;
# 		my $filename = pop @{ [ split( "/", $_ ) ] };
# 		push( @filenames, $filename );
# 	}
# 	close(FOF);
# 	return (@filenames);
# }
# ###################################################
# ###################################################
# #--Make Pairs of the files.
# sub MAKEPAIRS {
# 	my ( $filenames, $outdir ) = @_;
# 	my @names      = @$filenames;
# 	my $count      = scalar @names;
# 	my (@newnames) = ();
# 	if ( $count % 2 != 0 ) {
# 		print STDERR "\nOdd number of files given, please check Input file.\n";
# 		exit;
# 	}
# 	else {
# 		for ( my $i = 0 ; $i < scalar(@names) ; $i += 2 ) {
# 			chomp( $names[$i] );
# 			chomp( $names[ $i + 1 ] );
# 			push( @newnames, "$names[$i],$names[$i+1]" );
# 		}
# 	}
# 	return (@newnames);
# }
#####################################
#####################################
#Read data related to samples as well as barcodes.
# sub ReadSampleFile {
# 	my ( $sampleFile, $projectName, $outdir ) = @_;
# 	my (
# 		@fcId,        @lane,    @sampleId, @sampleRef, @index,
# 		@description, @control, @recipe,   @operator,  @sampleProject
# 	) = ();
# 	my $sampleFileName = "";
# 	if ( $sampleFile =~ /\// ) {
# 		$sampleFileName = pop @{ [ split( "/", $sampleFile ) ] };
# 	}
# 	else {
# 		$sampleFileName = $sampleFile;
# 	}
# 	open( SAMPLEFILE, $sampleFile )
# 	  || die "Cannot open SAMPLEFILE:$sampleFile,$!\n";
# 	while (<SAMPLEFILE>) {
# 		next if $. == 1;
# 		my @dataCols = split( ",", $_ );
# 		if ( $dataCols[0] ) { push( @fcId,          $dataCols[0] ); }
# 		if ( $dataCols[1] ) { push( @lane,          $dataCols[1] ); }
# 		if ( $dataCols[2] ) { push( @sampleId,      $dataCols[2] ); }
# 		if ( $dataCols[3] ) { push( @sampleRef,     $dataCols[3] ); }
# 		if ( $dataCols[4] ) { push( @index,         $dataCols[4] ); }
# 		if ( $dataCols[5] ) { push( @description,   $dataCols[5] ); }
# 		if ( $dataCols[6] ) { push( @control,       $dataCols[6] ); }
# 		if ( $dataCols[7] ) { push( @recipe,        $dataCols[7] ); }
# 		if ( $dataCols[8] ) { push( @operator,      $dataCols[8] ) }
# 		if ( $dataCols[9] ) { push( @sampleProject, $dataCols[9] ); }
# 	}
# 	close(SAMPLEFILE);
# 	if ( !-e "$outdir/$sampleFileName" ) {
# 		`cp $sampleFile $outdir/$sampleFileName`;
# 	}
# 	return (
# 		\@fcId,     \@lane,        \@sampleId, \@sampleRef,
# 		\@index,    \@description, \@control,  \@recipe,
# 		\@operator, \@sampleProject
# 	);
# }
# #####################################
# #####################################
# #Read data related to samples as well as barcodes from title file.
# sub ReadTitleFile {
# 	my ( $titleFile, $outdir ) = @_;
# 	my @barcode      = ();
# 	my @pool         = ();
# 	my @sampleId     = ();
# 	my @collabId     = ();
# 	my @patientId    = ();
# 	my @class        = ();
# 	my @sampleType   = ();
# 	my @inputNg      = ();
# 	my @libraryYeild = ();
# 	my @poolInput    = ();
# 	my @baitVersion  = ();
# 	my @fof          = ();
# 	my @newfof       = ();
# 	open( TFH, $titleFile )
# 	  || die "Cannot open file TitleFile:$titleFile, $!\n";

# 	while (<TFH>) {
# 		next if ( $. == 1 );
# 		my @dataCols = split( "\t", $_ );
# 		my @newDatacols =
# 		  grep( s/\s*$//g, @dataCols );    #remove whitespace if any
# 		push( @barcode,      $newDatacols[0] );
# 		push( @pool,         $newDatacols[1] );
# 		push( @sampleId,     $newDatacols[2] );
# 		push( @collabId,     $newDatacols[3] );
# 		push( @patientId,    $newDatacols[4] );
# 		push( @class,        $newDatacols[5] );
# 		push( @sampleType,   $newDatacols[6] );
# 		push( @inputNg,      $newDatacols[7] );
# 		push( @libraryYeild, $newDatacols[8] );
# 		push( @poolInput,    $newDatacols[9] );
# 		push( @baitVersion,  $newDatacols[10] );
# 	}
# 	close(TFH);

# 	my $poolName         = $pool[0];
# 	my $newtitleFileName = $poolName . "_title.txt";
# 	if ( !-e "$outdir/$newtitleFileName" ) {
# 		`cp $titleFile $outdir/$newtitleFileName`;
# 	}
# 	return (
# 		\@barcode,      \@pool,      \@sampleId,   \@collabId,
# 		\@patientId,    \@class,     \@sampleType, \@inputNg,
# 		\@libraryYeild, \@poolInput, \@baitVersion
# 	);
# }
#####################################
#####################################
#sort by barcode name:
# sub lowestNumber {
# 	my $files     = shift;
# 	my @filenames = split( ",", $files );
# 	my ($number)  = $filenames[0] =~ m/.*_bc(\d{1,4})_.*/g;
# 	return $number;
# }
#####################################
#####################################
#Merge data from reading data from the directory
# sub MergeDataFromDirectory {
# 	my @lane          = @$lane;
# 	my @sampleId      = @$sampleId;
# 	my @index         = @$index;
# 	my @titleBarcode  = @$barcode;
# 	my @titlePool     = @$pool;
# 	my @titleSampleId = @$titleSampleId;
# 	my %barcodes      = ();
# 	my %indexHash     = ();
# 	my %titleInfo     = ();
# 	my @notifyNames   = ();
# 	my $newIndex;
# 	my $name;
# 	my $Null           = "NULL";
# 	my @parseFilenames = ();
# 	my $now            = time;
# 	open( BARCODEFILE, $barcodeFile )
# 	  || die "Cannot open BARCODEFILE:$barcodeFile,$!\n";

# 	while (<BARCODEFILE>) {
# 		next if ( $. == 1 );
# 		my @dataCols = split( "\t", $_ );
# 		$dataCols[0] =~ s/\s//g;
# 		$dataCols[1] =~ s/\s//g;
# 		$barcodes{ $dataCols[0] }  = $dataCols[1];
# 		$indexHash{ $dataCols[1] } = $dataCols[0];
# 	}
# 	close(BARCODEFILE);
# 	print "Running merge jobs on SGE at " . localtime() . "\n";
# 	if ( $fastqSource =~ /DMP|PATH/i ) {
# 		foreach my $i ( 0 .. $#titleBarcode ) {
# 			my $newIndex = $titleBarcode[$i];
# 			my $name =
# 			    $titleSampleId[$i] . "_"
# 			  . $titleBarcode[$i] . "_"
# 			  . $titlePool[$i];
# 			my $read1ListName = "";
# 			my $read2ListName = "";
# 			foreach my $j ( 0 .. $#sampleId ) {
# 				if ( $index[$j] eq $indexHash{ $titleBarcode[$i] } ) {
# 					$read1ListName .=
# 					    $datadir . "/"
# 					  . $sampleId[$j] . "_"
# 					  . $index[$j] . "_L00"
# 					  . $lane[$j]
# 					  . "_R1_001.fastq.gz ";
# 					$read2ListName .=
# 					    $datadir . "/"
# 					  . $sampleId[$j] . "_"
# 					  . $index[$j] . "_L00"
# 					  . $lane[$j]
# 					  . "_R2_001.fastq.gz ";
# 				}
# 			}
# 			my $read1Name = $outdir . "/" . $name . "_L000_R1_mrg.fastq.gz";
# 			my $read2Name = $outdir . "/" . $name . "_L000_R2_mrg.fastq.gz";
# 			if (    ( -e $read1Name )
# 				and ( ( -s $read1Name ) != 0 )
# 				and ( -e $read2Name )
# 				and ( ( -s $read2Name ) != 0 ) )
# 			{
# 				print
# "Files:\n$read1Name\n$read2Name\n they exists and process will not run to merge files.\n";
# 				push( @notifyNames,    $Null );
# 				push( @notifyNames,    $Null );
# 				push( @parseFilenames, "$read1Name,$read2Name" );
# 			}
# 			else {
# 				my $notifyMergeCMD = "$outdir/Notify.csh";

# 				#Read1
# 				my $read1MergeCMD =
# 				  "'$ZCAT $read1ListName | $GZIP > $read1Name'";

# #`qsub -q all.q -V -wd $outdir -N MergeRead1.$newIndex.$i.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead1.$newIndex.$i.$$.err -o /dev/null -b y "/bin/zcat $read1ListName | gzip > $read1Name"`;
# 				launchQsub(
# 					$read1MergeCMD,
# 					$outdir,
# 					"8G",
# 					"/dev/null",
# 					"MergeRead1.$newIndex.$i.$$.stderr",
# 					"1",
# 					"$queue",
# 					"MergeRead1.$newIndex.$i.$$",
# 					"Null"
# 				);

# #`qsub -q all.q -V -wd $outdir -hold_jid MergeRead1.$newIndex.$i.$$ -N NotifyMR.Read1.$i.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read1.$i.$$.stat -b y "$outdir/Notify.csh"`;
# 				launchQsub(
# 					$notifyMergeCMD,
# 					$outdir,
# 					"2G",
# 					"NotifyMR.Read1.$i.$$.stat",
# 					"NotifyMR.Read1.$i.$$.stderr",
# 					"1",
# 					"$queue",
# 					"NotifyMR.Read1.$i.$$",
# 					"MergeRead1.$newIndex.$i.$$"
# 				);

# 				#Read2
# 				my $read2MergeCMD =
# 				  "'$ZCAT $read2ListName | $GZIP > $read2Name'";

# #`qsub -q all.q -V -wd $outdir -N MergeRead2.$newIndex.$i.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead2.$newIndex.$i.$$.err -o /dev/null -b y "/bin/zcat $read2ListName | gzip > $read2Name"`;
# 				launchQsub(
# 					$read2MergeCMD,
# 					$outdir,
# 					"8G",
# 					"/dev/null",
# 					"MergeRead2.$newIndex.$i.$$.stderr",
# 					"1",
# 					"$queue",
# 					"MergeRead2.$newIndex.$i.$$",
# 					"Null"
# 				);

# #`qsub -q all.q -V -wd $outdir -hold_jid MergeRead2.$newIndex.$i.$$ -N NotifyMR.Read2.$i.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read2.$i.$$.stat -b y "$outdir/Notify.csh"`;
# 				launchQsub(
# 					$notifyMergeCMD,
# 					$outdir,
# 					"2G",
# 					"NotifyMR.Read2.$i.$$.stat",
# 					"NotifyMR.Read2.$i.$$.stderr",
# 					"1",
# 					"$queue",
# 					"NotifyMR.Read2.$i.$$",
# 					"MergeRead2.$newIndex.$i.$$"
# 				);
# 				push( @notifyNames,    "NotifyMR.Read1.$i.$$.stat" );
# 				push( @notifyNames,    "NotifyMR.Read2.$i.$$.stat" );
# 				push( @parseFilenames, "$read1Name,$read2Name" );
# 			}
# 		}
# 		&WaitToFinish( $outdir, @notifyNames );
# 		$now = time - $now;
# 		print "Finished running merge jobs on SGE at " . localtime() . "\n";
# 		printf(
# 			"Total running time: %02d:%02d:%02d\n\n",
# 			int( $now / 3600 ),
# 			int( ( $now % 3600 ) / 60 ),
# 			int( $now % 60 )
# 		);
# 		my (@sortedparseFilenames) =
# 		  sort { lowestNumber($a) <=> lowestNumber($b) } @parseFilenames;
# 		return ( \@sortedparseFilenames );
# 	}
# 	else {
# 		for ( my $i = 0 ; $i < scalar(@titleBarcode) ; $i++ ) {
# 			$titleInfo{ $titleBarcode[$i] } =
# 			    $titleSampleId[$i] . "_"
# 			  . $titleBarcode[$i] . "_"
# 			  . $titlePool[$i];
# 		}
# 		for ( my $sampleNum = 0 ;
# 			$sampleNum < scalar(@sampleId) ; $sampleNum++ )
# 		{
# 			my $read1ListName =
# 			    $datadir . "/"
# 			  . $sampleId[$sampleNum] . "_"
# 			  . $index[$sampleNum] . "_L00"
# 			  . $lane[$sampleNum]
# 			  . "_R1_*.fastq.gz";
# 			my $read2ListName =
# 			    $datadir . "/"
# 			  . $sampleId[$sampleNum] . "_"
# 			  . $index[$sampleNum] . "_L00"
# 			  . $lane[$sampleNum]
# 			  . "_R2_*.fastq.gz";
# 			if ( exists $barcodes{ $index[$sampleNum] } ) {
# 				$newIndex = $barcodes{ $index[$sampleNum] };
# 				if ( exists $titleInfo{$newIndex} ) {
# 					$name = $titleInfo{$newIndex};
# 				}
# 				else {
# 					print
# "The barcode $newIndex doesnot exists in the title file. Cannot move ahead. Please check and rerun.\n";
# 					exit;
# 				}
# 			}
# 			else {
# 				print
# "The barcode sequence $barcodes{$index[$sampleNum]} does not exists in barcode file. Cannot move ahead. Please check and rerun.\n";
# 				exit;
# 			}
# 			my $read1Name =
# 			    $outdir . "/"
# 			  . $name . "_L00"
# 			  . $lane[$sampleNum]
# 			  . "_R1_mrg.fastq.gz";
# 			my $read2Name =
# 			    $outdir . "/"
# 			  . $name . "_L00"
# 			  . $lane[$sampleNum]
# 			  . "_R2_mrg.fastq.gz";

# 			#Run the qsub command to merge the files.
# 			#Read1
# 			if (    ( -e $read1Name )
# 				and ( ( -s $read1Name ) != 0 )
# 				and ( -e $read2Name )
# 				and ( ( -s $read2Name ) != 0 ) )
# 			{
# 				print
# "Files:\n$read1Name\n$read2Name\n they exists and process will not run to merge files.\n";
# 				push( @notifyNames,    $Null );
# 				push( @notifyNames,    $Null );
# 				push( @parseFilenames, "$read1Name,$read2Name" );
# 				next;
# 			}
# 			else {
# 				my $notifyMergeCMD = "$outdir/Notify.csh";

# 				#Read1
# 				my $read1MergeCMD =
# 				  "'$ZCAT $read1ListName | $GZIP > $read1Name'";

# #`qsub -q all.q -V -wd $outdir -N MergeRead1.$newIndex.$sampleNum.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead1.$newIndex.$sampleNum.$$.err -o /dev/null -b y "/bin/zcat $read1ListName | gzip > $read1Name"`;
# 				launchQsub(
# 					$read1MergeCMD,
# 					$outdir,
# 					"8G",
# 					"/dev/null",
# 					"MergeRead1.$newIndex.$sampleNum.$$.stderr",
# 					"1",
# 					"$queue",
# 					"MergeRead1.$newIndex.$sampleNum.$$",
# 					"Null"
# 				);

# #`qsub -q all.q -V -wd $outdir -hold_jid MergeRead1.$newIndex.$sampleNum.$$ -N NotifyMR.Read1.$sampleNum.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read1.$sampleNum.$$.stat -b y "$outdir/Notify.csh"`;
# 				launchQsub(
# 					$notifyMergeCMD,
# 					$outdir,
# 					"2G",
# 					"NotifyMR.Read1.$sampleNum.$$.stat",
# 					"NotifyMR.Read1.$sampleNum.$$.stderr",
# 					"1",
# 					"$queue",
# 					"NotifyMR.Read1.$sampleNum.$$",
# 					"MergeRead1.$newIndex.$sampleNum.$$"
# 				);

# 				#Read2
# 				my $read2MergeCMD =
# 				  "'$ZCAT $read2ListName | $GZIP > $read2Name'";

# #`qsub -q all.q -V -wd $outdir -N MergeRead2.$newIndex.$sampleNum.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead2.$newIndex.$sampleNum.$$.err -o /dev/null -b y "/bin/zcat $read2ListName | gzip > $read2Name"`;
# 				launchQsub(
# 					$read2MergeCMD,
# 					$outdir,
# 					"8G",
# 					"/dev/null",
# 					"MergeRead2.$newIndex.$sampleNum.$$.stderr",
# 					"1",
# 					"$queue",
# 					"MergeRead2.$newIndex.$sampleNum.$$",
# 					"Null"
# 				);

# #`qsub -q all.q -V -wd $outdir -hold_jid MergeRead2.$newIndex.$sampleNum.$$ -N NotifyMR.Read2.$sampleNum.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read2.$sampleNum.$$.stat -b y "$outdir/Notify.csh"`;
# 				launchQsub(
# 					$notifyMergeCMD,
# 					$outdir,
# 					"2G",
# 					"NotifyMR.Read2.$sampleNum.$$.stat",
# 					"NotifyMR.Read2.$sampleNum.$$.stderr",
# 					"1",
# 					"$queue",
# 					"NotifyMR.Read2.$sampleNum.$$",
# 					"MergeRead2.$newIndex.$sampleNum.$$"
# 				);
# 				push( @notifyNames,    "NotifyMR.Read1.$sampleNum.$$.stat" );
# 				push( @notifyNames,    "NotifyMR.Read2.$sampleNum.$$.stat" );
# 				push( @parseFilenames, "$read1Name,$read2Name" );
# 			}
# 		}
# 		&WaitToFinish( $outdir, @notifyNames );
# 		$now = time - $now;
# 		print "Finished running merge jobs on SGE at " . localtime() . "\n";
# 		printf(
# 			"Total running time: %02d:%02d:%02d\n\n",
# 			int( $now / 3600 ),
# 			int( ( $now % 3600 ) / 60 ),
# 			int( $now % 60 )
# 		);
# 	}
# 	my (@sortedparseFilenames) =
# 	  sort { lowestNumber($a) <=> lowestNumber($b) } @parseFilenames;
# 	return ( \@sortedparseFilenames );
# }
# #####################################
# #####################################
# #Do Mapping which includes:
# #Run BWA mem
# #Sort SAM
# sub DoMapping {
# 	my ($filenames) = @_;
# 	my ( @names, @notifyNames, @SAFilenames, @SamFilenames, @sortedBamFilenames,
# 		@MarkDuplicatesBamFilenames )
# 	  = ();
# 	if ($filenames) { (@names) = @$filenames; }
# 	if ( ( scalar(@names) == 0 ) and ($fof) ) {
# 		my @fnames = &GetNames( $fof, $outdir );
# 		@names = &MAKEPAIRS( \@fnames, $outdir );
# 	}
# 	my (@sortedparseFilenames) =
# 	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
# 	@names = @sortedparseFilenames;
# 	my @clippedFilenames = ();

# 	#my @SAFilenames = ();
# 	#my @SamFilenames = ();
# 	#my @sortedBamFilenames = ();
# 	#my @notifyNames = ();
# 	my $now         = time;
# 	my %adaptorList = ();
# 	open( ADAPTORFILE, $adaptorFile )
# 	  or die "DoMapping:Cannot open $adaptorFile,Error:$!\n";
# 	while (<ADAPTORFILE>) {
# 		my @dataCols = split( "\t", $_ );
# 		$dataCols[0] =~ s/\s//g;
# 		$dataCols[1] =~ s/\s//g;
# 		$adaptorList{ $dataCols[0] } = $dataCols[1];
# 	}
# 	close(ADAPTORFILE);

# 	#Running Cutadapt through Trim Galore
# 	print "Started runing clipping jobs on SGE\n";
# 	for ( my $i = 0 ; $i < scalar @names ; $i++ ) {
# 		my ( $file1, $file2 ) = split( ",", $names[$i] );

# 		#print "$file1\n$file2\n";
# 		my ( $read1clipped, $read2clipped, $notifyname ) =
# 		  &RunTrimGalore( $file1, $file2, $outdir, \%adaptorList, $i );
# 		push( @notifyNames,      $notifyname );
# 		push( @clippedFilenames, "$read1clipped,$read2clipped" );
# 	}

# 	#waiting for adapter trimming to finish
# 	&WaitToFinish( $outdir, @notifyNames );
# 	print "Finished running clipping jobs on SGE \n";
# 	$now = time - $now;
# 	printf(
# 		"Total Adaptor Clipping run time: %02d:%02d:%02d\n\n",
# 		int( $now / 3600 ),
# 		int( ( $now % 3600 ) / 60 ),
# 		int( $now % 60 )
# 	);

# 	#Running BwaMem
# 	$now = time;
# 	print "Started runing bwa mem jobs on SGE at " . localtime() . "\n";
# 	@notifyNames = ();
# 	for ( my $i = 0 ; $i < scalar(@clippedFilenames) ; $i++ ) {
# 		my ( $file1, $file2 ) = split( ",", $clippedFilenames[$i] );

# 		#print "$file1\n$file2\n";
# 		my ( $samFilename, $notifyname ) =
# 		  &RunBwaMem( $file1, $file2, $outdir, $i );
# 		push( @notifyNames, $notifyname );
# 		$samFilename = basename($samFilename);
# 		push( @SamFilenames, "$samFilename" );
# 	}

# 	#waiting for bwa aln to finish
# 	&WaitToFinish( $outdir, @notifyNames );
# 	print "Finished running bwa mem jobs on SGE at " . localtime() . "\n";
# 	$now = time - $now;
# 	printf(
# 		"Total running time: %02d:%02d:%02d\n\n",
# 		int( $now / 3600 ),
# 		int( ( $now % 3600 ) / 60 ),
# 		int( $now % 60 )
# 	);

# 	#Run Sort Sam
# 	$now = time;
# 	print "Started running Sort Sam jobs on SGE at " . localtime() . "\n";
# 	@notifyNames = ();
# 	for ( my $i = 0 ; $i < scalar(@SamFilenames) ; $i++ ) {
# 		my ( $sortedBamFile, $notifyname ) =
# 		  &RunSortSam( $SamFilenames[$i], $outdir, $i );
# 		push( @notifyNames,        $notifyname );
# 		push( @sortedBamFilenames, $sortedBamFile );
# 	}

# 	#waiting for sort sam to finish
# 	&WaitToFinish( $outdir, @notifyNames );
# 	$now = time - $now;
# 	print "Finished running Sort Sam jobs on SGE at " . localtime() . "\n";

# 	#Run MarkDuplicates
# 	$now = time;
# 	print "Started running Mark Duplicates jobs on SGE at "
# 	  . localtime() . "\n";
# 	@notifyNames = ();
# 	for ( my $i = 0 ; $i < scalar(@sortedBamFilenames) ; $i++ ) {
# 		my ( $MarkDuplicatesBamFile, $notifyname ) =
# 		  &RunMarkDuplicates( $sortedBamFilenames[$i], $outdir, $i );
# 		push( @notifyNames,                $notifyname );
# 		push( @MarkDuplicatesBamFilenames, $MarkDuplicatesBamFile );
# 	}

# 	#waiting for Mark Duplicates to finish
# 	&WaitToFinish( $outdir, @notifyNames );
# 	$now = time - $now;
# 	print "Finished running Mark Duplicates jobs on SGE at "
# 	  . localtime() . "\n";
# 	printf(
# 		"Total running time: %02d:%02d:%02d\n\n",
# 		int( $now / 3600 ),
# 		int( ( $now % 3600 ) / 60 ),
# 		int( $now % 60 )
# 	);
# 	return ( \@MarkDuplicatesBamFilenames );
# }
# #####################################
# #####################################
# #Clip adapter sequences.
# sub RunTrimGalore {
# 	my ( $file1, $file2, $outdir, $adaptorList, $id ) = @_;
# 	my %barcodeList = %$adaptorList;
# 	my ($barcode)   = $file1 =~ /.*_(bc\d+)_.*/;
# 	my $adapter1    = $barcodeList{$barcode};
# 	my $adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
# 	my ($basename1)  = $file1 =~ /(.*)\.fastq.gz/;
# 	my $outFilename1 = "$basename1" . "_cl.fastq.gz";
# 	my ($basename2)  = $file2 =~ /(.*)\.fastq.gz/;
# 	my $outFilename2 = "$basename2" . "_cl.fastq.gz";

# 	if (    ( -e "$outFilename1" )
# 		and ( ( -s "$outFilename1" ) != 0 )
# 		and ( -e "$outFilename2" )
# 		and ( ( -s "$outFilename2" ) != 0 ) )
# 	{
# 		print
# "Files:\n$outFilename1\n$outFilename2\n they exists and process will not run to clip adapters in them.\n";
# 		return ( "$outFilename1", "$outFilename2", 'NULL' );
# 	}
# 	else {
# 		my $clipReadCMD =
# "$PERL $TrimGalore --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o $outdir -a $adapter1 -a2 $adapter2 $file1 $file2";

# #`$QSUB -q $queue -V -wd $outdir -N Clipping.$id.$$ -o Clipping.$id.$$.stdout -e Clipping.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1  -b y "$PERL $TrimGalore --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o $outdir -a $adapter1 -a2 $adapter2 $file1 $file2"`;
# 		launchQsub( $clipReadCMD, $outdir, "2G", "Clipping.$id.$$.stdout",
# 			"Clipping.$id.$$.stderr", "1", $queue, "Clipping.$id.$$", "Null" );

# #`$QSUB -q $queue -V -wd $outdir -hold_jid Clipping.$id.$$ -N NotifyCR.$id.$$ -e NotifyCR.$id.$$.stderr -o NotifyCR.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
# 		my $notifyClipCMD = "$outdir/Notify.csh";
# 		launchQsub(
# 			$notifyClipCMD,           $outdir,
# 			"2G",                     "NotifyCR.$id.$$.stat",
# 			"NotifyCR.$id.$$.stderr", "1",
# 			"$queue",                 "NotifyCR.$id.$$",
# 			"Clipping.$id.$$"
# 		);
# 	}
# 	return ( "$outFilename1", "$outFilename2", "NotifyCR.$id.$$.stat" );
# }
# #####################################
# #####################################
# #BWA MEM to align fastq.
# sub RunBwaMem {
# 	my ( $fastq1, $fastq2, $outdir, $id ) = @_;
# 	my ($basename) = $fastq1 =~ /(.*)_R1.*\.fastq.gz/;
# 	my $outFilename = "$basename" . "_mrg_cl_aln.sam";
# 	if ( $basename =~ /\// ) {
# 		$basename = basename($basename);
# 	}
# 	my @sampleDetails = split( "_bc", $basename );
# 	my $sampleId      = $sampleDetails[0];
# 	my ($barcode)     = $basename =~ /.*_(bc\d+)_.*/;
# 	my ($pool)        = $basename =~ /.*bc\d+_(.*)_L\d{1,3}_.*/;
# 	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
# 		print
# "Files:\n$outFilename\n they exists and process will not run to make \"_aln.sam\" file.\n";
# 		return ( "$outFilename", 'NULL' );
# 	}
# 	else {
# 		my $bwaCMD =
# "\"$bwa mem -t 4 -PM -R \'\@RG\\tID:$basename\\tLB:$id\\tSM:$sampleId\\tPL:Illumina\\tPU:$barcode\\tCN:MSKCC\' $refFile $fastq1 $fastq2\"";

# #`qsub -q all.q -wd $outdir -N bwaMem.$id.$$ -l h_vmem=6G,virtual_free=6G -pe smp $nprocessors -o $outFilename -e /dev/null -b y "$bwa mem -t 4 -PM -R \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' $refFile $fastq1 $fastq2"`;
# 		launchQsub( $bwaCMD, $outdir, "6G", $outFilename,
# 			"bwaMem.$id.$$.stderr", $nprocessors, $queue,
# 			"bwaMem.$id.$$", "Null" );

# #`qsub -q all.q -V -wd $outdir -hold_jid bwaMem.$id.$$ -N NotifyBwaMem.$id.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyBwaMem.$id.$$.stat -b y "$outdir/Notify.csh"`;
# 		my $notifyBwaCMD = "$outdir/Notify.csh";
# 		launchQsub(
# 			$notifyBwaCMD,                $outdir,
# 			"2G",                         "NotifyBwaMem.$id.$$.stat",
# 			"NotifyBwaMem.$id.$$.stderr", "1",
# 			"$queue",                     "NotifyBwaMem.$id.$$",
# 			"bwaMem.$id.$$"
# 		);
# 	}
# 	return ( "$outFilename", "NotifyBwaMem.$id.$$.stat" );
# }
# #####################################
# #####################################
# #Sort Sam file
# sub RunSortSam {
# 	my ( $samFile, $outdir, $id ) = @_;
# 	my $outFilename = $samFile;
# 	$outFilename =~ s/\.sam/_srt\.bam/;
# 	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
# 		print
# "Files:\n$outFilename\n they exists and process will not run to make \"_srt.bam\" file.\n";
# 		return ( "$outFilename", 'NULL' );
# 	}
# 	else {
# 		my $sortsamCMD =
# "$JAVA -Xmx4g -jar $PICARD/SortSam.jar I=$samFile O=$outFilename SO=coordinate TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT";
# 		launchQsub( $sortsamCMD, $outdir, "8G", "SortSam.$id.$$.stdout",
# 			"SortSam.$id.$$.stderr", "1", $queue, "SortSam.$id.$$", "Null" );

# #`qsub -q all.q -wd $outdir -N SortSam.$id.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -o /dev/null -e /dev/null -b y "$JAVA -Xmx4g -jar $PICARD/SortSam.jar I=$samFile O=$outFilename SO=coordinate TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
# #`qsub -q all.q -V -wd $outdir -hold_jid SortSam.$id.$$ -N NotifySortSam.$id.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifySortSam.$id.$$.stat -b y "$outdir/Notify.csh"`;
# 		my $notifySortsamCMD = "$outdir/Notify.csh";
# 		launchQsub(
# 			$notifySortsamCMD,             $outdir,
# 			"2G",                          "NotifySortSam.$id.$$.stat",
# 			"NotifySortSam.$id.$$.stderr", "1",
# 			"$queue",                      "NotifySortSam.$id.$$",
# 			"SortSam.$id.$$"
# 		);
# 	}
# 	return ( "$outFilename", "NotifySortSam.$id.$$.stat" );
# }
# #####################################
# #####################################
# #Mark Duplicates in Bam
# sub RunMarkDuplicates {
# 	my ( $bamFile, $outdir, $id ) = @_;
# 	my $outFilename     = $bamFile;
# 	my $metricsFilename = $bamFile;
# 	$outFilename =~ s/\.bam/_MD\.bam/g;
# 	$metricsFilename =~ s/\.bam/_MD\.metrics/g;
# 	if ( ( -e "$outFilename" ) and ( ( -s "$outFilename" ) != 0 ) ) {
# 		print
# "Files:\n$outFilename\n they exists and process will not run to make \"_MD.bam\" file.\n";
# 		return ( "$outFilename", 'NULL' );
# 	}
# 	else {

# #`qsub -q all.q -V -wd $outdir -N MD.$id.$$  -o /dev/null -e /dev/null -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA -Xmx4g -jar $PICARD/MarkDuplicates.jar I=$bamFile O=$outFilename ASSUME_SORTED=true METRICS_FILE=$metricsFilename TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
# 		my $mdCMD =
# "$JAVA -Xmx4g -jar $PICARD/MarkDuplicates.jar I=$bamFile O=$outFilename ASSUME_SORTED=true METRICS_FILE=$metricsFilename TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT";
# 		launchQsub( $mdCMD, $outdir, "8G", "MD.$id.$$.stdout",
# 			"MD.$id.$$.stderr", "1", $queue, "MD.$id.$$", "Null" );

# #`qsub -q all.q -V -wd $outdir -hold_jid MD.$id.$$ -N NotifyMD.$id.$$ -e NotifyMD.$id.$$.stderr -o NotifyMD.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
# 		my $notifyMdCMD = "$outdir/Notify.csh";
# 		launchQsub( $notifyMdCMD, $outdir, "2G", "NotifyMD.$id.$$.stat",
# 			"NotifyMD.$id.$$.stderr", "1", "$queue", "NotifyMD.$id.$$",
# 			"MD.$id.$$" );
# 	}
# 	return ( "$outFilename", "NotifyMD.$id.$$.stat" );
# }
# #####################################
# #####################################
# #This will calculate and compile metrics for BAM files:
# sub CalcHsMetrics {
# 	my ($filenames) = @_;
# 	my @names = ();
# 	if ($filenames) { (@names) = @$filenames; }
# 	if ( ( scalar(@names) == 0 ) and ($fof) ) {
# 		@names = &GetNames( $fof, $outdir );
# 	}
# 	my @notifyNames = ();
# 	my (@sortedparseFilenames) =
# 	  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
# 	@names = @sortedparseFilenames;
# 	##################
# 	#Calculate HsMetrics
# 	$now = time;
# 	print "Started running metrics calculation jobs on SGE at "
# 	  . localtime() . "\n";
# 	for ( my $i = 0 ; $i < scalar(@names) ; $i++ ) {
# 		my $waitFileNames = &RunHsMetrics( $names[$i], $outdir, $i );
# 		foreach my $waitName (@$waitFileNames) {
# 			push( @notifyNames, $waitName );
# 		}
# 	}

# 	#waiting for metrics calculations to finish
# 	&WaitToFinish( $outdir, @notifyNames );
# 	$now = time - $now;
# 	print "Finished running metrics calculation jobs on SGE at "
# 	  . localtime() . "\n";
# 	printf(
# 		"Total running time: %02d:%02d:%02d\n\n",
# 		int( $now / 3600 ),
# 		int( ( $now % 3600 ) / 60 ),
# 		int( $now % 60 )
# 	);
# 	return ( \@names );
# }
# #####################################
# #####################################
# #Run Picard HsMetrics
# sub RunHsMetrics {
# 	my ( $bamFile, $outdir, $id ) = @_;
# 	my ($basename)        = $bamFile =~ /(.*)\.bam/;
# 	my $HSmetricsFilename = $basename . ".HSmetrics.txt";
# 	my @notifynames       = ();

# 	#Calculate Hybrid Selection specific metrics
# 	if ( ( -e "$HSmetricsFilename" ) and ( ( -s "$HSmetricsFilename" ) != 0 ) )
# 	{
# 		print
# "Files:\n$HSmetricsFilename\n they exists and process will not run to make \".HSmetrics.txt\" file.\n";
# 		push( @notifynames, "NULL" );
# 	}
# 	else {
# 		my $hsMetricsCMD =
# "$JAVA -Xmx4g -jar $PICARD/CalculateHsMetrics.jar I=$bamFile O=$HSmetricsFilename BI=$baitIntervalFile TI=$targetIntervalFile REFERENCE_SEQUENCE=$refFile TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT";

# #`qsub -q all.q -wd $outdir -N HSmetrics.$id.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -o /dev/null -e /dev/null -b y "$JAVA -Xmx4g -jar $PICARD/CalculateHsMetrics.jar I=$bamFile O=$HSmetricsFilename BI=$baitIntervalFile TI=$targetIntervalFile REFERENCE_SEQUENCE=$refFile TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
# 		launchQsub( $hsMetricsCMD, $outdir, "8G", "HSmetrics.$id.$$.stdout",
# 			"HSmetrics.$id.$$.stderr", "1", $queue, "HSmetrics.$id.$$",
# 			"Null" );

# #`qsub -q all.q -V -wd $outdir -hold_jid HSmetrics.$id.$$ -N NotifyHSmetrics.$id.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyHSmetrics.$id.$$.stat -b y "$outdir/Notify.csh"`;
# 		my $notifyHsMetricsCMD = "$outdir/Notify.csh";
# 		launchQsub(
# 			$notifyHsMetricsCMD,             $outdir,
# 			"2G",                            "NotifyHSmetrics.$id.$$.stat",
# 			"NotifyHSmetrics.$id.$$.stderr", "1",
# 			"$queue",                        "NotifyHSmetrics.$id.$$",
# 			"HSmetrics.$id.$$"
# 		);
# 		push( @notifynames, "NotifyHSmetrics.$id.$$.stat" );
# 	}
# 	return ( \@notifynames );
# }
# #####################################
# #####################################
# #This will help to call:
# #Somatic SVs: Delly
sub CallStructuralVariants {
	#my ($filenames) = @_;
	my @names       = keys(%BamSampleHash);
	my @FilterData  = ();
	#if ($filenames) { (@names) = @$filenames; }

	#print "F:$fof\n";
	#if ( ( scalar(@names) == 0 ) and ($fof) ) {
	#	@names = &GetNames( $fof, $outdir );
	#}

	 ### Modified Code To Use Pairing File for Delly begin ###
    #my %BamPatientMap = ();
    #foreach my $file (@names)
    #{
    #        my ($cur_sample_id) = $file =~ /(.*)_bc\d{1,4}_/;
    #        $BamPatientMap{$cur_sample_id} = $file;
    #}
     ### Modified Code To Use Pairing File for Delly end ###


	my @notifyNames = ();
	#tie( my %groupedFilenames, 'Tie::IxHash' );
	#my @somaticSVfiles     = ();
	#my @CoveragePerSample  = ();
	#my @meanCoverageValues = ();
	#tie( my %coverageForNormal, 'Tie::IxHash' );
	tie( my %NormalPerFile,     'Tie::IxHash' );
	#my $standardNormal;
	#my @stdNormalsList;
	#my @stdNormalsListIds;
	my $now = time;
	#my (@sortedparseFilenames) =
	#  sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	#@names = @sortedparseFilenames;

	#my ($poolName) = $names[0] =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
	my $NormalUsed = $poolName . "_NormalUsedInSVcalling.txt";
	my @fileNames  = ();
	my ( $waitFileNames, $svOutdir );

	# if ( $names[0] =~ /\// ) {
	# 	foreach (@names) {
	# 		my $filename = pop @{ [ split( "/", $_ ) ] };
	# 		push( @fileNames, $filename );
	# 	}
	# 	@names = @fileNames;
	# }

	#Call Somatic SVs
	print "Started running Somatics Variant jobs on SGE at " . localtime() . "\n";
#	my $fCount = 0;

# 	#Group the files
# 	foreach my $file (@names) {
# 		if ( exists $groupedFilenames{ @$patientId[$fCount] } ) {

# 			#print "$file:$fCount:@$patientId[$fCount]\n";
# 			my $files = $groupedFilenames{ @$patientId[$fCount] };
# 			$files = "$files" . ",$file";
# 			$groupedFilenames{ @$patientId[$fCount] } = "$files";
# 		}
# 		else {
# 			#print "$file:$fCount:@$patientId[$fCount]\n";
# 			$groupedFilenames{ @$patientId[$fCount] } = "$file";
# 		}
# 		$fCount++;
# 	}

# 	#Get Mean Coverage from HSmetrics file.
# 	my $poolNormalMeanCov;
# 	my $poolNormal;
# 	foreach my $file (@names) {

# 		#print "$file\n";
# 		my ($fileBarcode)  = $file =~ /.*_(bc\d+)_.*/;
# 		my ($fileSampleId) = $file =~ /(.*)_bc\d+_/;

# 		#print $fileBarcode . "_" . $fileSampleId, "\n";
# 		my $fileClass = $classPerBarcode{ $fileBarcode . "_" . $fileSampleId };

# 		#print "$fileClass\n";
# 		if ( $fileClass =~ m/Normal/i ) {
# 			if ( $fileClass =~ m/PoolN/i ) {
# 				$poolNormal = $file;
# 				my $HSmetricsFile = $file;
# 				$HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;
# 				open( FH, "$outdir/$HSmetricsFile" )
# 				  or die
# "CallSomaticSV:Cannot Open $outdir/$HSmetricsFile, Error:$!\n";
# 				while (<FH>) {
# 					next until ( $_ =~ /^BAIT_SET/ );
# 					while (<FH>) {
# 						next if ( ( $_ =~ /^BAIT_SET/ ) or ( $_ =~ /^\s$/ ) );
# 						my (@values) = split( "\t", $_ );
# 						$poolNormalMeanCov = $values[21];
# 					}
# 				}
# 				close(FH);
# 				next;
# 			}
# 			else {
# 				my $HSmetricsFile = $file;
# 				$HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;

# 				#print "HS:$HSmetricsFile\n";
# 				open( FH, "$outdir/$HSmetricsFile" )
# 				  or die
# 				  "Cannot Open HSmetricsFile:$outdir/$HSmetricsFile, $!\n";
# 				my $meanCov;
# 				my $CovForFile;
# 				while (<FH>) {
# 					next until ( $_ =~ /^BAIT_SET/ );
# 					while (<FH>) {
# 						next if ( ( $_ =~ /^BAIT_SET/ ) or ( $_ =~ /^\s$/ ) );
# 						my (@values) = split( "\t", $_ );

# 						#print "MeanCOv:$values[21]\n";
# 						$CovForFile = $values[21];
# 						$meanCov    = $values[21];
# 					}
# 				}
# 				close(FH);
# 				$coverageForNormal{$file} = $CovForFile;
# 				push( @CoveragePerSample, $meanCov );
# 			}
# 		}
# 		else {
# 			next;
# 		}
# 	}

# 	#Get file that will be used as standard normal
# 	my $maxCoverage = max @CoveragePerSample;

# 	#print "MAX:$maxCoverage\n";
# 	while ( ( my $key, my $value ) = each(%coverageForNormal) ) {
# 		if ( $value == $maxCoverage ) {
# 			$standardNormal = $key;
# 		}
# 		else {
# 			next;
# 		}
# 	}
# 	if ( !$standardNormal ) {
# 		$standardNormal = $stdNormal;
# 	}

# 	#print "SN:$standardNormal:$maxCoverage\n";
# 	#print "PN:$poolNormal:$poolNormalMeanCov\n";
# 	#Running Mutect and Somatic Indel Caller
 	my $count = 0;
# 	while ( ( my $key, my $value ) = each(%groupedFilenames) ) {
# 		my @files = split( ",", $value );

# 		# Section of Normal
# 		my @normalSamples = ();
# 		tie( my %coverageForSampleNormals, 'Tie::IxHash' );
# 		my @CoverageForMultipleNormal = ();
# 		my $normal;

# 		#	my $poolNormal;
# 		foreach my $file (@files) {
# 			my ($fileBarcode)  = $file =~ /.*_(bc\d+)_.*/;
# 			my ($fileSampleId) = $file =~ /(.*)_bc\d+_/;

# 			#print $fileBarcode . "_" . $fileSampleId, "\n";
# 			my $fileClass =
# 			  $classPerBarcode{ $fileBarcode . "_" . $fileSampleId };
# 			if ( $fileClass =~ m/Normal/i ) {
# 				push( @normalSamples, $file );
# 			}
# 			else {
# 				next;
# 			}
# 		}
# 		if ( scalar @normalSamples != 0 ) {
# 			foreach my $file (@normalSamples) {
# 				my $HSmetricsFile = $file;
# 				$HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;

# 				#print "HS:$HSmetricsFile\n";
# 				open( FH, "$outdir/$HSmetricsFile" )
# 				  or die
# "CallSomaticSV:Cannot Open $outdir/$HSmetricsFile, Error:$!\n";
# 				my $meanCov;
# 				my $CovForFile;
# 				while (<FH>) {
# 					next until ( $_ =~ /^BAIT_SET/ );
# 					while (<FH>) {
# 						next if ( ( $_ =~ /^BAIT_SET/ ) or ( $_ =~ /^\s$/ ) );
# 						my (@values) = split( "\t", $_ );

# 						#print "MeanCOv:$values[21]\n";
# 						$CovForFile = $values[21];
# 						$meanCov    = $values[21];
# 					}
# 					$coverageForSampleNormals{$file} = $CovForFile;
# 					push( @CoverageForMultipleNormal, $meanCov );
# 				}
# 				close(FH);
# 			}

# 			#Get file that will be used as normal
# 			my $maxCoverage = max @CoverageForMultipleNormal;

# 			#print "MAX:$maxCoverage\n";
# 			if ( scalar @normalSamples > 1 ) {
# 				while ( ( my $key, my $value ) =
# 					each(%coverageForSampleNormals) )
# 				{
# 					if ( ( $value == $maxCoverage ) and ( $value >= 50 ) ) {
# 						$normal = $key;
# 					}
# 					else {
# 						if ( $poolNormalMeanCov <= 50 ) {
# 							$normal = $standardNormal;
# 						}
# 						else {
# 							$normal = $poolNormal;
# 						}
# 					}
# 				}
# 			}
# 			else {
# 				if ( scalar @normalSamples == 1 ) {
# 					my $coverage =
# 					  $coverageForSampleNormals{ $normalSamples[0] };
# 					if ( $coverage >= 50 ) {
# 						$normal = $normalSamples[0];
# 					}
# 					else {
# 						if ( $poolNormalMeanCov <= 50 ) {
# 							$normal = $standardNormal;
# 						}
# 						else {
# 							$normal = $poolNormal;
# 						}
# 					}
# 				}
# 				else {
# 					$normal = $poolNormal;
# 				}
# 			}
# 		}
# 		else {
# 			if ( $poolNormalMeanCov <= 50 ) {
# 				$normal = $standardNormal;
# 			}
# 			else {
# 				$normal = $poolNormal;
# 			}
# 		}

		#RUN SV Calling Jobs
#		foreach my $file (@files) {
			foreach my $file (@names) {
			#my ($fileBarcode)  = $file =~ /.*_(bc\d+)_.*/;
			my $fileSampleId = $BamSampleHash{$file};
			my $fileClass;
			#  $classPerBarcode{ $fileBarcode . "_" . $fileSampleId };
			my $normalSampleId = "";
			if(exists($NormalTumorPair{$fileSampleId}) )
			{
				$fileClass = "Tumor";
				$normalSampleId = $NormalTumorPair{$fileSampleId};
			}
			else
			{
				$fileClass = "Normal"
			}
			next if ( $fileClass =~ m/Normal/i );
            my $normal = $SampleBamHash{$normalSampleId};


            #Check if the normal file is with full path
            #if ( $normal =~ /\// ) {
            #        $normal = pop @{ [ split( "/", $normal ) ] };
            #}
            #else {
            #        $normal = $normal;
            #}
			print "Final2:Tumor->$file\nNormal->$normal\n\n";
			#my ($tFileId)   = $file =~ /(.*)_$poolName\_/;
			#my ($nPoolName) = $normal =~ /.*_bc\d+_(.*)_L\d{1,3}.*/;
			#my ($nFileId)   = $normal =~ /(.*)_$nPoolName\_/;
			
			$NormalPerFile{$fileSampleId} = $normalSampleId;
			( $waitFileNames, $svOutdir ) =
			  &RunDelly( $normal, $file, $normalSampleId, $fileSampleId, $outdir, $count );

			foreach my $waitName (@$waitFileNames) {
				push( @notifyNames, $waitName );
			}
			push( @FilterData, "$svOutdir,$normalSampleId,$fileSampleId" );
			$count++;
		}
	#}

	exit 0 if (scalar(@FilterData) == 0); #no valid pairing found, no need to run the following steps

	&WaitToFinish( $outdir, @notifyNames );
	open( NFH, ">", "$outdir/$NormalUsed" )
	  || die "Cannot open NormalUsedinSVFile:$outdir/$NormalUsed;$!\n";
	while ( my ( $key, $value ) = each(%NormalPerFile) ) {
		print NFH "$key\t$value\n";
	}
	close(NFH);
	$now = time - $now;
	print "Finished running Germline and Somatic Variant jobs on SGE at "
	  . localtime() . "\n";
	printf(
		"Total running time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@names, \@FilterData );
}
#######################################
#######################################
#Run Delly
sub RunDelly {
	my ( $normal, $tumor, $nId, $tId, $outdir, $count ) = @_;
	my @waitFilenames     = ();
	my $date              = `date "+%Y%m%d"`;
	my $dellyOutdir       = $outdir . "/StrVarAnalysis";
	my $sampleTumorOutput = $dellyOutdir . "/" . $tId;
	my ( $tFlag, $nFlag );

	#my $stdNormals = join( " ", @$stdNormalsList );
	#Make Ouput dir
	if ( !( -d "$dellyOutdir" ) ) {
		`mkdir $dellyOutdir`;
	}
	else {
		if ( $count == 0 ) {
			warn "$dellyOutdir exists !!\n";
		}
	}

	#for making link to files
	my $NormalBai = $normal;
	$NormalBai =~ s/\.bam/\.bai/;
	my ($TumorBai) = $tumor;
	$TumorBai =~ s/\.bam/\.bai/;

	my $lnNormal = pop @{ [ split( "/", $normal ) ] };
	chomp($lnNormal);
	my ($lnNormalBai) = $lnNormal . ".bai";
	my $lnTumor = pop @{ [ split( "/", $tumor ) ] };
	chomp($lnTumor);
	my ($lnTumorBai) = $lnTumor . ".bai";


	#Make Tumor Sample Output Dir
	if ( -d "$sampleTumorOutput" ) {
		warn "$sampleTumorOutput exists !!\n";
		$tFlag = 1;
	}
	else {
		`mkdir $sampleTumorOutput`;
		`ln -s $tumor $sampleTumorOutput/`;
		`ln -s $normal $sampleTumorOutput/`;
		`ln -s $TumorBai $sampleTumorOutput/$lnTumorBai`;
		`ln -s $NormalBai $sampleTumorOutput/$lnNormalBai`;
		$tFlag = 2;

		#foreach (@$stdNormalsList)
		#{
		#	chomp($_);
		#	`ln -s $_ $sampleTumorOutput/`;
		#	my ($stdNormalBai) = $_ =~ s/\.bam/\.bai/;
		#	my ($lnStdNormalBai) = $_ . ".bai";
		#	`ln -s $outdir/$stdNormalBai $sampleTumorOutput/$lnStdNormalBai`;
		#}
	}
	
	$normal = $lnNormal;
	$tumor =  $lnTumor;	

	#Assign Queue
	my $runQueue = $queue;

	#Notify CMD
	my $notify_cmd = "$outdir/Notify.csh";

	#Delly Tumor CMD
	my $dellyT_cmd =
"$DELLY -t DEL -g $refFile -x $ExcludeRegions -q $MAPQ -o $tId\_del.vcf $tumor $normal";
	my $dellyT_jname  = "delly_$tId.$$.$count";
	my $dellyT_stdout = $dellyT_jname . ".stdout";
	my $dellyT_stderr = $dellyT_jname . ".stderr";

	#Duppy Tumor CMD
	my $duppyT_cmd =
"$DELLY -t DUP -g $refFile -x $ExcludeRegions -q $MAPQ -o $tId\_dup.vcf $tumor $normal";
	my $duppyT_jname  = "duppy_$tId.$$.$count";
	my $duppyT_stdout = $duppyT_jname . ".stdout";
	my $duppyT_stderr = $duppyT_jname . ".stderr";

	#Invy Tumor CMD
	my $invyT_cmd =
"$DELLY -t INV -g $refFile -x $ExcludeRegions -q $MAPQ -o $tId\_inv.vcf  $tumor $normal";
	my $invyT_jname  = "invy_$tId.$$.$count";
	my $invyT_stdout = $invyT_jname . ".stdout";
	my $invyT_stderr = $invyT_jname . ".stderr";

	#Jumpy Tumor CMD
	my $jumpyT_cmd =
"$DELLY -t TRA -g $refFile -x $ExcludeRegions -q $MAPQ -o $tId\_jmp.vcf  $tumor $normal";
	my $jumpyT_jname  = "jumpy_$tId.$$.$count";
	my $jumpyT_stdout = $jumpyT_jname . ".stdout";
	my $jumpyT_stderr = $jumpyT_jname . ".stderr";

=begin
	#Jumpy Tumor CMD
	my $jumpyT_cmd =
"$DELLY_TransLocations -p -q $MAPQ -g $refFile -i $tId -x $ExcludeRegionsCTX -o $tId\_jmp.txt -r $tId\_jmpmerged.txt -b $tId\_jmpbrkpts.txt -k $tId\_jmpbrkptsmerged.txt $tumor";
	my $jumpyT_jname  = "jumpy_$tId.$$.$count";
	my $jumpyT_stdout = $jumpyT_jname . ".stdout";
	my $jumpyT_stderr = $jumpyT_jname . ".stderr";

	#Jumpy Normal CMD
	my $jumpyN_cmd =
"$DELLY_TransLocations -p -q $MAPQ -g $refFile -i $nId  -x $ExcludeRegionsCTX -o $nId\_jmp.txt -r $nId\_jmpmerged.txt -b $nId\_jmpbrkpts.txt -k $nId\_jmpbrkptsmerged.txt $normal";
	my $jumpyN_jname  = "jumpy_$nId.$$.$count";
	my $jumpyN_stdout = $jumpyN_jname . ".stdout";
	my $jumpyN_stderr = $jumpyN_jname . ".stderr";
=cut

	#Notify CMD values
	my $notifyT_hjname =
	  "$dellyT_jname,$duppyT_jname,$invyT_jname,$jumpyT_jname";
	my $notifyT_jname  = "NotifyDelly.$tId.$$.$count";
	my $notifyT_stdout = $notifyT_jname . ".stat";
	my $notifyT_stderr = $notifyT_jname . ".stderr";

	#Launch only if folder does not exists
	if ( $tFlag == 2 ) {
		&launchQsub(
			$dellyT_cmd,    $sampleTumorOutput,
			"8G",           $dellyT_stdout,
			$dellyT_stderr, "2",
			$runQueue,      $dellyT_jname,
			"Null"
		);
		&launchQsub(
			$duppyT_cmd,    $sampleTumorOutput,
			"8G",           $duppyT_stdout,
			$duppyT_stderr, "2",
			$runQueue,      $duppyT_jname,
			"Null"
		);
		&launchQsub(
			$invyT_cmd,    $sampleTumorOutput,
			"8G",          $invyT_stdout,
			$invyT_stderr, "2",
			$runQueue,     $invyT_jname,
			"Null"
		);
		&launchQsub(
			$jumpyT_cmd,    $sampleTumorOutput,
			"8G",           $jumpyT_stdout,
			$jumpyT_stderr, "2",
			$runQueue,      $jumpyT_jname,
			"Null"
		);

		#&launchQsub(
		#	$jumpyN_cmd,    $sampleTumorOutput,
		#	"10G",          $jumpyN_stdout,
		#	$jumpyN_stderr, "1",
		#	$runQueue,      $jumpyN_jname,
		#	"Null"
		#);
		&launchQsub(
			$notify_cmd,     $outdir,
			"2G",            $notifyT_stdout,
			$notifyT_stderr, "1",
			$runQueue,       $notifyT_jname,
			$notifyT_hjname
		);
		push( @waitFilenames, $notifyT_stdout );
	}
	else {
		print
		  "Resuts for Tumor:$tId sample exists. Thus Delly would not be ran\n";
		push( @waitFilenames, "NULL" );
	}
	return ( \@waitFilenames, "$sampleTumorOutput" );
}
#####################################
#####################################
#This will help to Call Filter
#SV module
sub FilterStructuralVariants {
	my ($filenames) = @_;
	my $now = time;
	my @names;
	my @notifyNames    = ();
	my @somaticSVfiles = ();
	my $count          = 0;
	if ($filenames) { (@names) = @$filenames; }

	#print "F:$fof\n";
	my $NormalUsed = $poolName . "_NormalUsedInSVcalling.txt";
	#if ( ( scalar(@names) == 0 ) and ($fof) ) {
	#	@names = &GetNames( $fof, $outdir );
	#}
	open( NFH, "$outdir/$NormalUsed" )
	  || die "Cannot open NormalUsedinSVFile:$outdir/$NormalUsed;$!\n";
	while (<NFH>) {
		chomp;
		my ( $tumorId, $normalId ) = split( "\t", $_ );
		my $dellyOutdir       = $outdir . "/StrVarAnalysis";
		my $sampleTumorOutput = $dellyOutdir . "/" . $tumorId;
		my $entry             = "$sampleTumorOutput,$normalId,$tumorId";
		my (
			$waitFileName, $delFilterVcf, $dupFilterVcf,
			$invFilterVcf, $jmpOutFile
		) = RunFilterStructuralVariants( $outdir, $entry, $count );
		my ( $svOutdir, $nFileId, $tFileId ) = split( ",", $entry );
		push( @somaticSVfiles,
			"$entry,$delFilterVcf,$dupFilterVcf,$invFilterVcf,$jmpOutFile" );
		push( @notifyNames, $waitFileName );
		$count++;
	}
	close(NFH);
	&WaitToFinish( $outdir, @notifyNames );
	$now = time - $now;
	print "Finished Filtering Variant jobs on SGE at " . localtime() . "\n";
	printf(
		"Total running time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@names, \@somaticSVfiles );
}
#######################################
#######################################
#Run the the cmd as qsub
sub launchQsub {
	my (
		$cmd,        $outdir, $mem,     $stdout, $stderr,
		$processors, $queue,  $jobname, $holdjobname
	) = @_;

	
	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", work_dir => "$outdir");
	my $additionalParams = Schedule::additionalParams(%addParams);


	my $qcmd = "";

	#Run Job with hold job id

	if ( $holdjobname ne "Null" ) {
	    ###$qcmd = "$QSUB -q $queue -V -v OMP_NUM_THREADS=2 -wd $outdir -hold_jid $holdjobname -N $jobname -o $stdout -e $stderr -l h_vmem=$mem,virtual_free=$mem -pe smp $processors -b y $cmd";
	    
	    my %stdParams = (scheduler => "$scheduler", job_name => "$jobname", job_hold => "$holdjobname", cpu => "$processors", mem => "$mem", cluster_out => "$stdout", cluster_error => "$stderr");
	    my $standardParams = Schedule::queuing(%stdParams);
	    $qcmd = "$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $cmd";
	    
	    eval {
		print "CMD:$qcmd\n";
		`$qcmd`;
	    };
	    if ($@) {
		print "Problem Sumitting the Qusb Command.Error: $@\n";
		exit(1);
	    }
	}
	
	#Run Jobs without hold job Id
	if ( $holdjobname eq "Null" ) {
	    ###$qcmd = "$QSUB -q $queue -V -v OMP_NUM_THREADS=2 -wd $outdir -N $jobname -o $stdout -e $stderr -l h_vmem=$mem,virtual_free=$mem -pe smp $processors -b y $cmd";

	    my %stdParams = (scheduler => "$scheduler", job_name => "$jobname", cpu => "$processors", mem => "$mem", cluster_out => "$stdout", cluster_error => "$stderr");
	    my $standardParams = Schedule::queuing(%stdParams);
	    $qcmd = "$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $cmd";

	    eval {
		print "CMD:$qcmd\n";
		`$qcmd`;
	    };
	    if ($@) {
		print "Problem Sumitting the Qusb Command.Error: $@\n";
		exit(1);
	    }
	}
	return;
}
#####################################
#####################################
#This will help to Filter:
#Somatic SVs
sub RunFilterStructuralVariants {
	my ( $outdir, $entry, $count ) = @_;
	my ( $svOutdir, $nFileId, $tFileId ) = split( ",", $entry );
	my $id           = basename($svOutdir);
	my $delVcf       = $id . "_del.vcf";
	my $delFilterVcf = $id . "_del_stdfilter.vcf";
	my $dupVcf       = $id . "_dup.vcf";
	my $dupFilterVcf = $id . "_dup_stdfilter.vcf";
	my $invVcf       = $id . "_inv.vcf";
	my $invFilterVcf = $id . "_inv_stdfilter.vcf";

	#my $jmpOutFile   = $id . "_jmp_stdfilter.txt";
	my $jmpVcf       = $id . "_jmp.vcf";
	my $jmpFilterVcf = $id . "_jmp_stdfilter.vcf";
	my $runQueue     = $queue;

	#Notify CMD
	my $notify_cmd = "$outdir/Notify.csh";

	#Variables for Deletion
	my $delcmd =
"$PERL $FilterSV -i $svOutdir/$delVcf -hsf $HotspotFile -tid $tFileId -nid $nFileId -tos DEL -outdir $svOutdir -o $delFilterVcf -ope $OverallSupportingReads -osr $OverallSupportingSplitReads -opeh $OverallSupportingReadsHotspot -osrh $OverallSupportingSplitReadsHotspot -stpe $SampleTumorSupportingReads -snpe $SampleNormalSupportingReads -stpeh $SampleTumorSupportingReadsHotspot -snpeh $SampleNormalSupportingReadsHotspot -svlen $LengthOfSV -omq $OverallMapq -omqh $OverallMapqHotspot -stgqf $SampleTumorGenotypeQualityFilter -stgqfh $SampleTumorGenotypeQualityFilterHotspot -sngqf $SampleNormalGenotypeQualityFilter -sngqfh $SampleNormalGenotypeQualityFilterHotspot";
	my $delFilter_jname  = "DelFilter_" . $id . "_" . $$;
	my $delFilter_stdout = $delFilter_jname . ".stdout";
	my $delFilter_stderr = $delFilter_jname . ".stderr";

	#Variables for Duplication
	my $dupcmd =
"$PERL $FilterSV -i $svOutdir/$dupVcf -hsf $HotspotFile -tid $tFileId -nid $nFileId -tos DUP -outdir $svOutdir -o $dupFilterVcf -ope $OverallSupportingReads -osr $OverallSupportingSplitReads -opeh $OverallSupportingReadsHotspot -osrh $OverallSupportingSplitReadsHotspot -stpe $SampleTumorSupportingReads -snpe $SampleNormalSupportingReads -stpeh $SampleTumorSupportingReadsHotspot -snpeh $SampleNormalSupportingReadsHotspot -svlen $LengthOfSV -omq $OverallMapq -omqh $OverallMapqHotspot -stgqf $SampleTumorGenotypeQualityFilter -stgqfh $SampleTumorGenotypeQualityFilterHotspot -sngqf $SampleNormalGenotypeQualityFilter -sngqfh $SampleNormalGenotypeQualityFilterHotspot";
	my $dupFilter_jname  = "DupFilter_" . $id . "_" . $$;
	my $dupFilter_stdout = $dupFilter_jname . ".stdout";
	my $dupFilter_stderr = $dupFilter_jname . ".stderr";

	#Variables for Inversion
	my $invcmd =
"$PERL $FilterSV -i $svOutdir/$invVcf -hsf $HotspotFile -tid $tFileId -nid $nFileId -tos INV -outdir $svOutdir -o $invFilterVcf -ope $OverallSupportingReads -osr $OverallSupportingSplitReads -opeh $OverallSupportingReadsHotspot -osrh $OverallSupportingSplitReadsHotspot -stpe $SampleTumorSupportingReads -snpe $SampleNormalSupportingReads -stpeh $SampleTumorSupportingReadsHotspot -snpeh $SampleNormalSupportingReadsHotspot -svlen $LengthOfSV -omq $OverallMapq -omqh $OverallMapqHotspot -stgqf $SampleTumorGenotypeQualityFilter -stgqfh $SampleTumorGenotypeQualityFilterHotspot -sngqf $SampleNormalGenotypeQualityFilter -sngqfh $SampleNormalGenotypeQualityFilterHotspot";
	my $invFilter_jname  = "InvFilter_" . $id . "_" . $$;
	my $invFilter_stdout = $invFilter_jname . ".stdout";
	my $invFilter_stderr = $invFilter_jname . ".stderr";

	#Variables for Translocations
	my $jmpcmd =
"$PERL $FilterSV -i $svOutdir/$jmpVcf -hsf $HotspotFile -tid $tFileId -nid $nFileId -tos TRA -outdir $svOutdir -o $jmpFilterVcf -ope $OverallSupportingReads -osr $OverallSupportingSplitReads -opeh $OverallSupportingReadsHotspot -osrh $OverallSupportingSplitReadsHotspot -stpe $SampleTumorSupportingReads -snpe $SampleNormalSupportingReads -stpeh $SampleTumorSupportingReadsHotspot -snpeh $SampleNormalSupportingReadsHotspot -svlen $LengthOfSV -omq $OverallMapq -omqh $OverallMapqHotspot -stgqf $SampleTumorGenotypeQualityFilter -stgqfh $SampleTumorGenotypeQualityFilterHotspot -sngqf $SampleNormalGenotypeQualityFilter -sngqfh $SampleNormalGenotypeQualityFilterHotspot";
	my $jmpFilter_jname  = "JmpFilter_" . $id . "_" . $$;
	my $jmpFilter_stdout = $jmpFilter_jname . ".stdout";
	my $jmpFilter_stderr = $jmpFilter_jname . ".stderr";

	#Notify CMD values
	my $notifyT_hjname =
	  "$delFilter_jname,$dupFilter_jname,$invFilter_jname,$jmpFilter_jname";
	my $notifyT_jname  = "NotifyFilter.$tFileId.$$.$count";
	my $notifyT_stdout = $notifyT_jname . ".stat";
	my $notifyT_stderr = $notifyT_jname . ".stderr";
	&launchQsub(
		$delcmd,           $svOutdir,
		"10G",             $delFilter_stdout,
		$delFilter_stderr, "1",
		$runQueue,         $delFilter_jname,
		"Null"
	);
	&launchQsub(
		$dupcmd,           $svOutdir,
		"10G",             $dupFilter_stdout,
		$dupFilter_stderr, "1",
		$runQueue,         $dupFilter_jname,
		"Null"
	);
	&launchQsub(
		$invcmd,           $svOutdir,
		"10G",             $invFilter_stdout,
		$invFilter_stderr, "1",
		$runQueue,         $invFilter_jname,
		"Null"
	);
	&launchQsub(
		$jmpcmd,           $svOutdir,
		"10G",             $jmpFilter_stdout,
		$jmpFilter_stderr, "1",
		$runQueue,         $jmpFilter_jname,
		"Null"
	);
	&launchQsub(
		$notify_cmd,     $outdir, "2G",      $notifyT_stdout,
		$notifyT_stderr, "1",     $runQueue, $notifyT_jname,
		$notifyT_hjname
	);
	return (
		$notifyT_stdout, $delFilterVcf, $dupFilterVcf,
		$invFilterVcf,   $jmpFilterVcf
	);
}
#####################################
#####################################
#This will help to Annotate:
#Somatic SVs
sub AnnotateStructuralVariants {
	my ( $filenames, $output ) = @_;
	my $now = time;
	my @names;
	my @notifyNames = ();
	my @AnnoSVfiles = ();
	my $count       = 0;
	if ($filenames) { (@names) = @$filenames; }
	my $NormalUsed = $poolName . "_NormalUsedInSVcalling.txt";

	#if ( ( scalar(@names) == 0 ) and ($fof) ) {
	#	@names = &GetNames( $fof, $outdir );
	#}
	open( NFH, "$outdir/$NormalUsed" )
	  || die "Cannot open NormalUsedinSVFile:$outdir/$NormalUsed;$!\n";
	while (<NFH>) {
		chomp;
		my ( $tumorId, $normalId ) = split( "\t", $_ );
		my $dellyOutdir       = $outdir . "/StrVarAnalysis";
		my $sampleTumorOutput = $dellyOutdir . "/" . $tumorId;
		my $entry             = "$sampleTumorOutput,$normalId,$tumorId";
		my ( $waitFileName, $delAnno, $dupAnno, $invAnno, $jmpAnno ) =
		  RunAnnotateStructuralVariants( $outdir, $entry, $count );
		my ( $svOutdir, $nFileId, $tFileId ) = split( ",", $entry );
		push( @AnnoSVfiles, "$entry,$delAnno,$dupAnno,$invAnno,$jmpAnno" );
		push( @notifyNames, $waitFileName );
		$count++;
	}
	close(NFH);
	&WaitToFinish( $outdir, @notifyNames );
	@notifyNames = ();
	my($dRangerAnnotatedFile) = &MergeAllFiles( \@AnnoSVfiles );
	my $finalFilePreFix = $poolName . "_AllAnnotatedSVs";
	#Annotate CMD values

	###my $annotate_cmd = "-v PYTHONPATH=$PYTHONPATH $PYTHON $AnnotateSV -r $RepeatRegionFile -d $DGvFile -c $CancerCensusFile -s $dRangerAnnotatedFile -o $outdir/$finalFilePreFix";

	my $annotate_cmd = "$PYTHON $AnnotateSV -r $RepeatRegionFile -d $DGvFile -c $CancerCensusFile -s $dRangerAnnotatedFile -o $outdir/$finalFilePreFix";

	my $annotate_jname  = "AnnotateSVs_" . $$;
	my $annotate_stdout = $annotate_jname . ".stdout";
	my $annotate_stderr = $annotate_jname . ".stderr";
	#Notify CMD values
	my $notify_cmd = "$outdir/Notify.csh";
	my $notify_hjname = "$annotate_jname";
	my $notify_jname  = "NotifyAnnotateSv.$$";
	my $notify_stdout = $notify_jname . ".stat";
	my $notify_stderr = $notify_jname . ".stderr";
	&launchQsub(
		$annotate_cmd,     $outdir, "2G",      $annotate_stdout,
		$annotate_stderr, "1",     $queue, $annotate_jname,
		"Null"
	);
	&launchQsub(
		$notify_cmd,     $outdir, "2G",      $notify_stdout,
		$notify_stderr, "1",     $queue, $notify_jname,
		$notify_hjname
	);
	push( @notifyNames, $notify_stdout );
	&WaitToFinish( $outdir, @notifyNames );
	$now = time - $now;
	print "Finished Filtering Variant jobs on SGE at " . localtime() . "\n";
	printf(
		"Total running time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);
	return ( \@names, \@AnnoSVfiles );
}
#####################################
#####################################
#MergeVCF CTX & dRangerAnnotation
sub MergeAllFiles {
	my ($annoFiles) = @_;
	my @AnnoFiles   = @$annoFiles;
	my $NormalUsed  = $poolName . "_NormalUsedInSVcalling.txt";
	my $outFile     = $poolName . "_All_dRangerAnnotatedSVs.txt";
	open( NFH, "$outdir/$NormalUsed" )
	  || die "Cannot open NormalUsedinSVFile:$outdir/$NormalUsed;$!\n";
	open( OFH, ">", "$outdir/$outFile" )
	  || die "Cannot open $outdir/$outFile;$!\n";
	print OFH
"TumorId\tNormalId\tChr1\tPos1\tChr2\tPos2\tSV_Type\tGene1\tGene2\tTranscript1\tTranscript2\tSite1Description\tSite2Description\tFusion\tConfidence\tComments\tConnection_Type\tSV_LENGTH\tMAPQ\tPairEndReadSupport\tSplitReadSupport\tBrkptType\tConsensusSequence\tTumorVariantCount\tTumorSplitVariantCount\tTumorReadCount\tTumorGenotypeQScore\tNormalVariantCount\tNormalSplitVariantCount\tNormalReadCount\tNormalGenotypeQScore\n";
	while (<NFH>) {
		chomp;
		my ( $tumorId, $normalId ) = split( "\t", $_ );
		my $dellyOutdir       = $outdir . "/StrVarAnalysis";
		my $sampleTumorOutput = $dellyOutdir . "/" . $tumorId;
		my $entry             = "$sampleTumorOutput,$normalId,$tumorId";
		my $delOut            = $tumorId . "_del_stdfilter.vcf";
		my $dupOut            = $tumorId . "_dup_stdfilter.vcf";
		my $invOut            = $tumorId . "_inv_stdfilter.vcf";
		my $jmpOut            = $tumorId . "_jmp_stdfilter.vcf";
		my $delAnnoOut        = $tumorId . "_del_stdfilter_dRangerOut.txt";
		my $dupAnnoOut        = $tumorId . "_dup_stdfilter_dRangerOut.txt";
		my $invAnnoOut        = $tumorId . "_inv_stdfilter_dRangerOut.txt";
		my $jmpAnnoOut        = $tumorId . "_jmp_stdfilter_dRangerOut.txt";

		#Merge VCF to dRanger
		my $delData =
		  MergeVCFwithdRanger( $sampleTumorOutput, $delOut, $delAnnoOut,
			$tumorId, $normalId );
		my $dupData =
		  MergeVCFwithdRanger( $sampleTumorOutput, $dupOut, $dupAnnoOut,
			$tumorId, $normalId );
		my $invData =
		  MergeVCFwithdRanger( $sampleTumorOutput, $invOut, $invAnnoOut,
			$tumorId, $normalId );
		my $jmpData =
		  MergeVCFwithdRanger( $sampleTumorOutput, $jmpOut, $jmpAnnoOut,
			$tumorId, $normalId );
		if ( $delData ne "Null" ) {
			tie( my %delHash, 'Tie::IxHash' );
			%delHash = %$delData;
			foreach my $record ( sort keys %delHash ) {
				my ( $chr1, $pos1, $chr2, $pos2 ) = split( ":", $record );
				my (
					$svType, $gene1,   $gene2, $transcript1, $transcript2, $site1,
					$site2,  $fusion,  $connectionType, $svlen,
					$mapq,   $peReads, $srReads,$bkptType,$consensusSeq,$tDV,$tRV,
					$tRC,    $tGQ,     $nDV, $nRV,            $nRC,
					$nGQ
				) = split( ";", $delHash{$record} );
				print OFH
"$tumorId\t$normalId\t$chr1\t$pos1\t$chr2\t$pos2\t$svType\t$gene1\t$gene2\t$transcript1\t$transcript2\t$site1\t$site2\t$fusion\t\t\t$connectionType\t$svlen\t$mapq\t$peReads\t$srReads\t$bkptType\t$consensusSeq\t$tDV\t$tRV\t$tRC\t$tGQ\t$nDV\t$nRV\t$nRC\t$nGQ\n";
			}
		}
		if ( $dupData ne "Null" ) {
			tie( my %dupHash, 'Tie::IxHash' );
			%dupHash = %$dupData;
			foreach my $record ( sort keys %dupHash ) {
				my ( $chr1, $pos1, $chr2, $pos2 ) = split( ":", $record );
				my (
					$svType, $gene1,   $gene2, $transcript1, $transcript2, $site1,
					$site2,  $fusion,  $connectionType, $svlen,
					$mapq,   $peReads, $srReads,$bkptType,$consensusSeq,$tDV,$tRV,
					$tRC,    $tGQ,     $nDV, $nRV,            $nRC,
					$nGQ
				) = split( ";", $dupHash{$record} );
				print OFH
"$tumorId\t$normalId\t$chr1\t$pos1\t$chr2\t$pos2\t$svType\t$gene1\t$gene2\t$transcript1\t$transcript2\t$site1\t$site2\t$fusion\t\t\t$connectionType\t$svlen\t$mapq\t$peReads\t$srReads\t$bkptType\t$consensusSeq\t$tDV\t$tRV\t$tRC\t$tGQ\t$nDV\t$nRV\t$nRC\t$nGQ\n";
			}
		}
		if ( $invData ne "Null" ) {
			tie( my %invHash, 'Tie::IxHash' );
			%invHash = %$invData;
			foreach my $record ( sort keys %invHash ) {
				my ( $chr1, $pos1, $chr2, $pos2 ) = split( ":", $record );
				my (
					$svType, $gene1,   $gene2, $transcript1, $transcript2, $site1,
					$site2,  $fusion,  $connectionType, $svlen,
					$mapq,   $peReads, $srReads,$bkptType,$consensusSeq,$tDV,$tRV,
					$tRC,    $tGQ,     $nDV, $nRV,            $nRC,
					$nGQ
				) = split( ";", $invHash{$record} );
				print OFH
"$tumorId\t$normalId\t$chr1\t$pos1\t$chr2\t$pos2\t$svType\t$gene1\t$gene2\t$transcript1\t$transcript2\t$site1\t$site2\t$fusion\t\t\t$connectionType\t$svlen\t$mapq\t$peReads\t$srReads\t$bkptType\t$consensusSeq\t$tDV\t$tRV\t$tRC\t$tGQ\t$nDV\t$nRV\t$nRC\t$nGQ\n";
			}
		}
		if ( $jmpData ne "Null" ) {
			tie( my %jmpHash, 'Tie::IxHash' );
			%jmpHash = %$jmpData;
			foreach my $record ( sort keys %jmpHash ) {
				my ( $chr1, $pos1, $chr2, $pos2 ) = split( ":", $record );
				my (
					$svType, $gene1,   $gene2, $transcript1, $transcript2, $site1,
					$site2,  $fusion,  $connectionType, $svlen,
					$mapq,   $peReads, $srReads,$bkptType,$consensusSeq,$tDV,$tRV,
					$tRC,    $tGQ,     $nDV, $nRV,            $nRC,
					$nGQ
				) = split( ";", $jmpHash{$record} );
				print OFH
"$tumorId\t$normalId\t$chr1\t$pos1\t$chr2\t$pos2\t$svType\t$gene1\t$gene2\t$transcript1\t$transcript2\t$site1\t$site2\t$fusion\t\t\t$connectionType\t$svlen\t$mapq\t$peReads\t$srReads\t$bkptType\t$consensusSeq\t$tDV\t$tRV\t$tRC\t$tGQ\t$nDV\t$nRV\t$nRC\t$nGQ\n";
			}
		}
	}
	close(NFH);
	close(OFH);
	return($outFile);
}
#####################################
#####################################
#MergeVCF & dRangerAnnotation
sub MergeVCFwithdRanger {
	my ( $svOutdir, $vcf, $file, $tumorId, $normalId ) = @_;
	tie( my %outHash, 'Tie::IxHash' );
	my $tag = CheckIfFileHasContent( $svOutdir, $file );
	if ( $tag == 1 ) {
		my %dRangerData = &Read_dRangerOut( $svOutdir, $file );
		my %vcfData = &Read_VCFout( $svOutdir, $vcf, $tumorId, $normalId );
		foreach my $key ( sort keys %vcfData ) {
			my $vcfValue    = $vcfData{$key};
			my $dRangeValue = $dRangerData{$key};

			#print "Key=>$key\tVCF=>$vcfValue\tDranger=>$dRangeValue\n";
			my ( $chr1, $pos1, $chr2, $pos2 ) = split( ":", $key );
			my ( $str1, $str2, $gene1, $gene2, $transcript1, $transcript2,$site1, $site2, $fusion ) =
			  split( ";", $dRangeValue );
			my ( $svlen, $mapq, $svType, $peReads, $srReads, $connectionType,
				$bkptType,$consensusSeq,$tGQ, $tFT, $tRC, $tDR, $tDV, $tRR, $tRV, $nGQ, $nFT, $nRC, $nDR, $nDV,$nRR, $nRV )
			  = split( ";", $vcfValue );
			$outHash{$key} =
"$svType;$gene1;$gene2;$transcript1;$transcript2;$site1;$site2;$fusion;$connectionType;$svlen;$mapq;$peReads;$srReads;$bkptType;$consensusSeq;$tDV;$tRV;$tRC;$tGQ;$nDV;$nRV;$nRC;$nGQ";
		}
	}
	else {
		return ("Null");
	}
	return ( \%outHash );
}
#####################################
#####################################
#MergeVCF & dRangerAnnotation
sub MergeCTXwithdRanger {
	my ( $svOutdir, $jmp, $file, $tumorId, $normalId ) = @_;
	tie( my %outHash, 'Tie::IxHash' );
	my $tag = CheckIfFileHasContent( $svOutdir, $file );
	if ( $tag == 1 ) {
		my %dRangerData = &Read_dRangerOut( $svOutdir, $file );
		my %ctxData = &Read_CTXout( $svOutdir, $jmp );
		foreach my $key ( sort keys %ctxData ) {
			my ( $chr1, $pos1, $chr2, $pos2 ) = split( ":", $key );
			my ( $str1, $str2, $gene1, $gene2, $transcript1, $transcript2, $site1, $site2, $fusion ) =
			  split( ":", $dRangerData{$key} );
			my (
				$Direction,    $TumorPE,    $TumorSR,    $TumorMAPQ,
				$CanBeSomatic, $NormalPos1, $NormalPos2, $NormalPE,
				$NormalSR,     $NormalMAPQ
			) = split( ":", $ctxData{$key} );
			$outHash{$key} =
"CTX;$gene1;$gene2;$transcript1;$transcript2;$site1;$site2;$fusion;$Direction;-;$TumorMAPQ;$TumorPE;$TumorSR;-;-;-;$NormalPE;-;-";
		}
	}
	else {
		return ("Null");
	}
	return ( \%outHash );
}
#####################################
#####################################
#Read dRanger Annotation & store them as hash
sub Read_dRangerOut {
	my ( $dir, $file ) = @_;
	my %outHash = ();
	open( DFH, "$dir/$file" ) || die "Cannot open file $dir/$file, Error:$!\n";
	while (<DFH>) {
		chomp;
		next if ( $. == 1 );
		my (
			$chr1,  $pos1,  $str1,  $chr2,  $pos2, $str2,
			$gene1, $site1,$transcript1, $gene2, $site2, $transcript2, $fusion
		) = split( "\t", $_ );
		$chr1 = "X" if ( $chr1 == 23 );
		$chr2 = "X" if ( $chr2 == 23 );
		$chr1 = "Y" if ( $chr1 == 24 );
		$chr2 = "Y" if ( $chr2 == 24 );
		$chr1 = "chr" . $chr1;
		$chr2 = "chr" . $chr2;
		$outHash{"$chr1:$pos1:$chr2:$pos2"} =
		  "$str1;$str2;$gene1;$gene2;$transcript1;$transcript2;$site1;$site2;$fusion";
	}
	close(DFH);
	return (%outHash);
}
#####################################
#####################################
#Read  vcf & store them as hash
sub Read_CTXout {
	my ( $dir, $file ) = @_;
	my %outHash = ();
	open( FH, "$dir/$file" ) || die "Cannot Open $dir/$file. Error:$!\n";
	while (<FH>) {
		chomp;
		next if ( $. == 1 );
		my (
			$Chr1,         $TumorPos1,  $Chr2,       $TumorPos2,
			$Direction,    $TumorPE,    $TumorSR,    $TumorMAPQ,
			$CanBeSomatic, $NormalPos1, $NormalPos2, $NormalPE,
			$NormalSR,     $NormalMAPQ
		) = split( "\t", $_ );
		$outHash{"$Chr1:$TumorPos1:$Chr2:$TumorPos2"} =
"$Direction;$TumorPE;$TumorSR;$TumorMAPQ;$CanBeSomatic;$NormalPos1;$NormalPos2;$NormalPE;$NormalSR;$NormalMAPQ";
	}
	close(FH);
	return (%outHash);
}
#####################################
#####################################
#Read  vcf & store them as hash
sub Read_VCFout {
	my ( $dir, $vcf, $tumorId, $normalId ) = @_;
	my %outHash = ();
	my @header  = ();
	my (
		$tumorIndex, $normalIndex, $infoIndex,
		$chrIndex,   $startIndex,  $filterIndex
	);
	open( VFH, "$dir/$vcf" ) || die "Cannot open file $dir/$vcf, Error:$!\n";
	while (<VFH>) {
		chomp($_);

		#Get the Header
		if ( $_ =~ m/^#/ ) {

			#Get what is tumor what is normal
			if ( $_ =~ m/^#CHROM/ ) {
				@header = split( "\t", $_ );
				for ( my $i = 0 ; $i < scalar(@header) ; $i++ ) {
					if ( $header[$i] =~ /$tumorId/ ) {
						$tumorIndex = $i;
					}
					if ( $header[$i] =~ /$normalId/ ) {
						$normalIndex = $i;
					}
					if ( $header[$i] =~ /INFO/ ) {
						$infoIndex = $i;
					}
					if ( $header[$i] =~ /CHROM/ ) {
						$chrIndex = $i;
					}
					if ( $header[$i] =~ /POS/ ) {
						$startIndex = $i;
					}
					if ( $header[$i] =~ /FILTER/ ) {
						$filterIndex = $i;
					}
				}
			}
		}
		else {

			#Get All the Lines
			my (@line) = split( "\t", $_ );

			#Get Chromosome
			my $svChr = $line[$chrIndex];

			#Get Start Position
			my $svStart = $line[$startIndex];

			#Get Filter String
			my $filterString = $line[$filterIndex];

			#Get Info column
			my @info = split( ";", $line[$infoIndex] );
			my %infoData = ();
			foreach my $infoTypes (@info) {
				if ( $infoTypes =~ m/=/ ) {
					my ( $key, $value ) = split( "=", $infoTypes );
					$infoData{$key} = $value;
				}
				else {
					$infoData{$infoTypes} = "";
				}
			}

			#get Tumor Genotype Information
			my ( $tGT, $tGL, $tGQ, $tFT, $tRC, $tDR, $tDV, $tRR, $tRV ) =
			  split( ":", $line[$tumorIndex] );

			#Get Normal Genotype Information
			my ( $nGT, $nGL, $nGQ, $nFT, $nRC, $nDR, $nDV, $nRR, $nRV ) =
			  split( ":", $line[$normalIndex] );
			my ( $svlen, $mapq, $peReads, $srReads, $svType, $svTNratio,
				$svChr2, $svEnd )
			  = 0.0;
			my ( $connectionType, $str1, $str2, $bkptType, $consensusSeq ) = "";

			#Calculate Tumor Noraml SV read ratio if its 0
			#($svTNratio) = 5 * $nDV;
			$bkptType = "IMPRECISE" if exists $infoData{"IMPRECISE"};
			$bkptType = "PRECISE" if exists $infoData{"PRECISE"};
			$svlen          = $infoData{"SVLEN"}  if exists $infoData{"SVLEN"};
			$mapq           = $infoData{"MAPQ"}   if exists $infoData{"MAPQ"};
			$svType         = $infoData{"SVTYPE"} if exists $infoData{"SVTYPE"};
			$peReads        = $infoData{"PE"}     if exists $infoData{"PE"};
			$srReads        = $infoData{"SR"}     if exists $infoData{"SR"};
			$svEnd          = $infoData{"END"}    if exists $infoData{"END"};
			$connectionType = $infoData{"CT"}     if exists $infoData{"CT"};
			$svChr2         = $infoData{"CHR2"}   if exists $infoData{"CHR2"};
			$consensusSeq   = $infoData{"CONSENSUS"} if exists $infoData{"CONSENSUS"};
			if ( !$svChr2 ) { $svChr2 = $svChr }
			$outHash{"$svChr:$svStart:$svChr2:$svEnd"} =
"$svlen;$mapq;$svType;$peReads;$srReads;$connectionType;$bkptType;$consensusSeq;$tGQ;$tFT;$tRC;$tDR;$tDV;$tRR;$tRV;$nGQ;$nFT;$nRC;$nDR;$nDV;$nRR;$nRV";
		}
	}
	close(VFH);
	return (%outHash);
}
#####################################
#####################################
#This will help to Filter:
#Somatic SVs
sub RunAnnotateStructuralVariants {
	my ( $outdir, $entry, $count ) = @_;
	my ( $svOutdir, $nFileId, $tFileId ) = split( ",", $entry );
	my $id         = basename($svOutdir);
	my $delAnnoIn  = $id . "_del_stdfilter_dRangerInput.txt";
	my $delAnnoOut = $id . "_del_stdfilter_dRangerOut.txt";
	my $dupAnnoIn  = $id . "_dup_stdfilter_dRangerInput.txt";
	my $dupAnnoOut = $id . "_dup_stdfilter_dRangerOut.txt";
	my $invAnnoIn  = $id . "_inv_stdfilter_dRangerInput.txt";
	my $invAnnoOut = $id . "_inv_stdfilter_dRangerOut.txt";
	my $jmpAnnoIn  = $id . "_jmp_stdfilter_dRangerInput.txt";
	my $jmpAnnoOut = $id . "_jmp_stdfilter_dRangerOut.txt";
	my $runQueue   = $queue;
	my $delTag     = &CheckIfFileHasContent( $svOutdir, $delAnnoIn );
	my $dupTag     = &CheckIfFileHasContent( $svOutdir, $dupAnnoIn );
	my $invTag     = &CheckIfFileHasContent( $svOutdir, $invAnnoIn );
	my $jmpTag     = &CheckIfFileHasContent( $svOutdir, $jmpAnnoIn );

	#Notify CMD
	my $notify_cmd     = "$outdir/Notify.csh";
	my $delAnno_jname  = "DelAnno_" . $id . "_" . $$;
	my $dupAnno_jname  = "DupAnno_" . $id . "_" . $$;
	my $invAnno_jname  = "InvAnno_" . $id . "_" . $$;
	my $jmpAnno_jname  = "JmpAnno_" . $id . "_" . $$;
	my $notify_jname   = "NotifyAnno.$tFileId.$$.$count";
	my (@checkProcess) = ();
	my (@annoOutFiles) = ();
	if ( $delTag == 1 ) {
		my $delAnnoCMD     = "$dRANGER $MCR $delAnnoIn $HG19MAT $delAnnoOut";
		my $delAnno_stdout = $delAnno_jname . ".stdout";
		my $delAnno_stderr = $delAnno_jname . ".stderr";
		push( @checkProcess, $delAnno_jname );
		&launchQsub(
			$delAnnoCMD,     $svOutdir,
			"5G",            $delAnno_stdout,
			$delAnno_stderr, "1",
			$runQueue,       $delAnno_jname,
			"Null"
		);
	}
	if ( $dupTag == 1 ) {
		my $dupAnnoCMD     = "$dRANGER $MCR $dupAnnoIn $HG19MAT $dupAnnoOut";
		my $dupAnno_stdout = $dupAnno_jname . ".stdout";
		my $dupAnno_stderr = $dupAnno_jname . ".stderr";
		push( @checkProcess, $dupAnno_jname );
		&launchQsub(
			$dupAnnoCMD,     $svOutdir,
			"5G",            $dupAnno_stdout,
			$dupAnno_stderr, "1",
			$runQueue,       $dupAnno_jname,
			"Null"
		);
	}
	if ( $invTag == 1 ) {
		my $invAnnoCMD     = "$dRANGER $MCR $invAnnoIn $HG19MAT $invAnnoOut";
		my $invAnno_stdout = $invAnno_jname . ".stdout";
		my $invAnno_stderr = $invAnno_jname . ".stderr";
		push( @checkProcess, $invAnno_jname );
		&launchQsub(
			$invAnnoCMD,     $svOutdir,
			"5G",            $invAnno_stdout,
			$invAnno_stderr, "1",
			$runQueue,       $invAnno_jname,
			"Null"
		);
	}
	if ( $jmpTag == 1 ) {
		my $jmpAnnoCMD     = "$dRANGER $MCR $jmpAnnoIn $HG19MAT $jmpAnnoOut";
		my $jmpAnno_stdout = $jmpAnno_jname . ".stdout";
		my $jmpAnno_stderr = $jmpAnno_jname . ".stderr";
		push( @checkProcess, $jmpAnno_jname );
		&launchQsub(
			$jmpAnnoCMD,     $svOutdir,
			"5G",            $jmpAnno_stdout,
			$jmpAnno_stderr, "1",
			$runQueue,       $jmpAnno_jname,
			"Null"
		);
	}
	my $notify_hjname;
	my $notify_stdout = $notify_jname . ".stat";
	my $notify_stderr = $notify_jname . ".stderr";

	#Notify CMD values
	if ( scalar @checkProcess > 1 ) {
		$notify_hjname = join( ",", @checkProcess );
	}
	elsif ( scalar @checkProcess == 1 ) {
		$notify_hjname = $checkProcess[0];
	}
	else {
		return ( "NULL", $delAnnoOut, $dupAnnoOut, $invAnnoOut, $jmpAnnoOut );
	}
	&launchQsub(
		$notify_cmd,    $outdir, "2G",      $notify_stdout,
		$notify_stderr, "1",     $runQueue, $notify_jname,
		$notify_hjname
	);
	return ( $notify_stdout, $delAnnoOut, $dupAnnoOut, $invAnnoOut,
		$jmpAnnoOut );
}
#####################################
#####################################
#Check if file has more then one line
sub CheckIfFileHasContent {
	my ( $svDir, $file ) = @_;
	my $tag;
	if ( ( -f "$svDir/$file" ) and ( -s "$svDir/$file" ) ) {
		open( FH, "$svDir/$file" )
		  || die "Cannot Open $svDir/$file.Error:$!\n";
		my $lines = 0;
		$lines++ while (<FH>);
		close(FH);
		if ( $lines > 1 ) {
			$tag = "1";
		}
		else {
			$tag = "2";
		}
	}
	else {
		$tag = "2";
	}
	return ($tag);
}

sub ReadNormalTumorPair
{
    open(PF, "$pair_file") or die "can't open pairing file $pair_file $!";
    print "reading pairing file: $pair_file\n";
    while(<PF>){
            chomp;
            my @data = split(/\s+/);
            if($data[0] && $data[1] && ($data[0] !~ /^NA$/i) && ($data[1] !~ /^NA$/i) )     
            {
                    $NormalTumorPair{$data[1]} = $data[0];
            }
    }
    close PF;
}

sub ReadBamList
{
    open(BL, "$bam_list") or die "can't open ban list file $bam_list $!";
    print "reading bam list file: $bam_list\n";
    while(<BL>){
            chomp;
            my @data = split(/\s+/);
            if($data[0] && $data[1])
            {
                    die "Could not find bam file $data[1]" if (!-e $data[1]);
                    $SampleBamHash{$data[0]} = $data[1];
                    $BamSampleHash{$data[1]} = $data[0];
            }
    }
    close BL;
}


sub CleanUpFiles
{
	my($src_dir, $dst_dir) = @_;
        `/bin/mkdir -p $dst_dir`;
        die "[ERROR]: Fail to create directory: $dst_dir\n" if(!-d $dst_dir);
	`find $src_dir -name \"*.vcf\" ! -name \"*stdfilter.vcf\" -exec mv {} $dst_dir \\;`;
	`rm -r $src_dir`;
}








