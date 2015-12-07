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
    $refFile_hg19,
    $refFile_b37,
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
    $DistanceBtwTumorNormalCTX,
    $dRANGER,
    $MCR,
    $pair_file,
    $bam_list,
    $genome,
    $scheduler,
    $priority_project,
    $priority_group
);

my $binPath = $Bin;
my $variant_config_file = $binPath  . "/variants_pipeline_config.txt";   #default variant config file
my $sv_config_file = $binPath . "/strvar/template_dmp_sv.conf"; #default sv config file

#--This variable holds the current time
my $now = time;

if (@ARGV < 1 or !GetOptions (
	'pre=s'                             => \$poolName,
	'variant_config=s'                => \$variant_config_file,  	
	'sv_config=s'                     => \$sv_config_file,       
	'outputDirectory=s'               => \$outdir,
	'pair=s'                            => \$pair_file,
	'bam_list=s'                        => \$bam_list,
	'genome=s'							=> \$genome,
	'scheduler=s'                       => \$scheduler,
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

if(!$genome || ($genome ne "hg19" && $genome ne "b37") )
{
    die "Please specify the genome build by using -genome. Valid values are 'hg19' and 'b37'\n";
}


if ($variant_config_file) {
	print "The variant configration file in use is $variant_config_file.\n";
}
else {
	print "Please speciy varaint configuration file by -variant_config\n";
	Usage();
	exit;
}

if ($sv_config_file) {
	print "The sv configration file in use is $sv_config_file.\n";
}
else {
	print "Please speciy sv configuration file by -sv_config\n";
	Usage();
	exit;
}

if ($outdir) {
	print "Results will be written in $outdir\n";
	`/bin/mkdir -p $outdir`;
	die "[ERROR]: Fail to create outdir: $outdir\n" if(!-d $outdir);
}
else {
	print "Please enter the directory where to write the data while/after being processed.See Usage\n";
	Usage();
	exit;
}

my %NormalTumorPair = ();
if($pair_file)
{
    ReadNormalTumorPair();
}
else
{

    die "Please specify the pairing file by using -pair PAIR_FILE\n";
}

my %SampleBamHash = ();
my %BamSampleHash = ();
if($bam_list)
{
    ReadBamList();
}
else
{
	die "Please specify the bam list file by using -bam_list BAM_LIST_FILE\n";
}

# This is variable is the path to the bin folder
my $genomeDataPath = $binPath . "/data/" . $genome . "/";
my $HG19MAT = $genomeDataPath  . "RefSeq_hg19.mat";
my $ExcludeRegions = $genomeDataPath . "human.hg19.excl.tsv";
my $HotspotFile = $genomeDataPath . "v3clin_hg19_structuralvariants_geneInterval.txt";
my $CancerCensusFile = $genomeDataPath . "cancer_gene_census.tsv";
my $DGvFile = $genomeDataPath . "DGv_Annotation.tsv";
my $RepeatRegionFile = $genomeDataPath . "repeatRegion.tsv";




#Load Variant Pipeline Main Configuration File
&verifyConfig($variant_config_file);

#Get SV Configration File details
my ($Version) = &GetConfiguration( $sv_config_file, $outdir );
my $PrintConfigFile = "RunConfigration_StrVar.txt";
open( VERSION, ">", "$outdir/$PrintConfigFile" ) || die "Cannot open $outdir/$PrintConfigFile, Reason:$!\n";
#Prin Version of tools and Files used
print VERSION "Tools|Files\tVersion\n";
while ( my ( $tools_files, $version ) = each(%$Version) ) {
	print VERSION "$tools_files\t$version\n";
}
close(VERSION);

#use default value to override config file
$TMPDIR = $outdir . "/tmp/";
`/bin/mkdir -p $TMPDIR`;
die "[ERROR]: Fail to create TMPDIR: $TMPDIR\n" if(!-d $TMPDIR);

# if ( !$mvFiles ) {
# 	print "Folders will be created and Files will be moved.\n";
# 	$mvFiles = 1;
# }
# else {
# 	if ( $mvFiles == 1 ) {
# 		print "Folders will be created and Files will be moved.\n";
# 	}
# 	if ( $mvFiles == 2 ) {
# 		print "Folders will not be created and Files will not be moved.\n";
# 	}
# }
# if ( !$stdNormal ) {
# 	$stdNormal = "NULL";
# }
# else {
# 	print "Starndard Normal is given: $stdNormal\n";
# }

# if ( !$fastqSource ) {
# 	print "Assume fastq files came from GCL.\n";
# 	$fastqSource = "GCL";
# }
# else {
# 	if ( $fastqSource ne "GCL" && $fastqSource ne "DMP" ) {
# 		print "Please indicate fastqSource. See Usage\n";
# 		Usage();
# 		exit;
# 	}
# }
# if ($barcodeFile) {
# 	print "The barcode file in use is $barcodeFile.\n";
# }
# else {
# 	print "Please enter the barcode file.See Usage\n";
# 	Usage();
# 	exit;
# }
# if ($adaptorFile) {
# 	print "The barcode file in use is $adaptorFile.\n";
# }
# else {
# 	print "Please enter the adaptor file.See Usage\n";
# 	Usage();
# 	exit;
# }

if ( !$TMPDIR ) {
	print "Path to temporary directory is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "TMPDIR=$TMPDIR\n";
}

# if ( !$baitIntervalFile ) {
# 	print "Bait Interval file is not given. See Usage.\n";
# 	Usage();
# 	exit;
# }
# else {
# 	print "BAIT_INTERVAL=$baitIntervalFile\n";
# }
# if ( !$targetIntervalFile ) {
# 	print "Target Interval file is not given. See Usage\n";
# 	Usage();
# 	exit;
# }
# else {
# 	print "Target_INTERVAL=$targetIntervalFile\n";
# }
if ( !$refFile ) {
	print "Reference Sequence file is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "Reference_Sequence=$refFile\n";
}

if ( !$HotspotFile ) {
	print "Hotspot location BED4 file is not given. See Usage\n";
	Usage();
	exit;
}
else {
	print "HotspotFile=$HotspotFile\n";
}
#if ( !$JAVA ) {
#	print "Path to java executables is not given. See Usage\n";
#	Usage();
#	exit;
#}
#else {
#	print "JAVA=$JAVA\n";
#}

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


&MakeCSH($outdir);

my @allProcess = split( ",", $process );

my $allProcessList  = join( ",", @allProcess );
my $numberOfProcess = scalar(@allProcess);
my $processCount    = 0;
my $parseFilenames;
while ( $processCount < $numberOfProcess ) {
	my $runProcess = shift(@allProcess);
	($parseFilenames) = &Select( $runProcess, $parseFilenames );
	$processCount++;
}

CleanUpFiles($outdir);


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
#--GET VARIANT CONFIGRATION DETAIL
sub verifyConfig{
    my ($paths) = @_;
    open(CONFIG, "$paths") || die "Can't open config file $paths $!";
    while(<CONFIG>){
        chomp;
        my @conf = split(/\s+/, $_);
        if($conf[0] =~ /DELLY/i){
            if(!-e "$conf[1]/delly"){
                die "CAN'T FIND delly IN $conf[1] $!";
            }
            $DELLY = "$conf[1]/delly";
        }
        elsif($conf[0] =~ /^PERL/i){
            if(!-e "$conf[1]/perl"){
                die "CAN'T FIND perl IN $conf[1] $!";
            }
            $PERL = "$conf[1]/perl";
        }
        elsif($conf[0] =~ /^HG19_FASTA/i && $genome eq "hg19"){
            if(!-e "$conf[1]"){
                die "CAN'T FIND hg19 fasta at $conf[1] $!";
            }
            $refFile = $conf[1];
        }
        elsif($conf[0] =~ /^B37_FASTA/i && $genome eq "b37"){
            if(!-e "$conf[1]"){
                die "CAN'T FIND b37 fasta at $conf[1] $!";
            }
            $refFile = $conf[1];
        }
        elsif($conf[0] =~ /^PYTHON/i){
            if(!-e "$conf[1]/python"){
                die "CAN'T FIND python IN $conf[1] $!";
            }
            $PYTHON = "$conf[1]/python";
        }
        elsif($conf[0] =~ /^dRANGER/i){
            if(!-e "$conf[1]/AnnotateSVs/run_AnnotateSVs_v2.sh"){
                die "CAN'T FIND dRANGER IN $conf[1] $!";
            }
            $dRANGER= "$conf[1]/AnnotateSVs/run_AnnotateSVs_v2.sh";
        }
        elsif($conf[0] =~ /^MCR/i){
            if(!-d "$conf[1]"){
                die "CAN'T FIND MCR at $conf[1] $!";
            }
            $MCR= $conf[1];
        }
   }
   close CONFIG;
}



###################################################
###################################################
#--GET SV CONFIGRATION DETAIL
sub GetConfiguration {
	my ($sv_config_file) = @_;
	my @data = ();
	tie( my %config,     'Tie::IxHash' );
	tie( my %location,   'Tie::IxHash' );
	tie( my %version,    'Tie::IxHash' );
	tie( my %parameters, 'Tie::IxHash' );
	print "Reading the configuration file\n";

	# change the default Input Record Separator.
	$/ = ">";
	open( CONFIG, "$sv_config_file" ) or die "Cannot open $sv_config_file. Error: $!\n";
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
		#$TMPDIR             = $location{"TMPDIR"};
		#$PERL               = $location{"PERL"};
		#$JAVA               = $location{"JAVA"};
		#$GATK               = $location{"GATK"};
		#$refFile            = $location{"Reference"};
		#$PICARD             = $location{"PICARD"};
		#$baitIntervalFile   = $location{"BaitInterval"};
		#$targetIntervalFile = $location{"TargetInterval"};
		#$bwa                = $location{"BWA"};
		#$PYTHON				= $location{"PYTHON"};
		#$PYTHONPATH			= $location{"PYTHONPATH"};
		#$samtools           = $location{"SAMTOOLS"};
		#$DELLY              = $location{"DELLY"};
		#$barcodeFile        = $location{"barcodeFile"};
		#$RepeatRegionFile   = $location{"RepeatRegionAnnotation"};
		#$DGvFile     		= $location{"DGvAnnotations"};
		#$CancerCensusFile = $location{"CosmicCensus"};
		#$RepeatRegionFile   = $location{"RepeatRegionAnnotation"};
		#$RefGeneFile        = $location{"RefGeneFile"};
		#$barcodeFile        = $location{"BarcodeKey"};
		#$adaptorFile        = $location{"AdaptorKey"};
		#$QSUB               = $location{"QSUB"};
		#$TrimGalore         = $location{"TrimGalore"};
		#$ZCAT               = $location{"ZCAT"};
		#$GZIP               = $location{"GZIP"};
		#$FilterSV           = $location{"FilterSV"};
		#$AnnotateSV			= $location{"AnnotateSV"};
		#$HotspotFile        = $location{"HotspotFile"};
		#$dRANGER            = $location{"dRANGER"};
		#$MCR                = $location{"MCR"};
		$queue              = $parameters{"SGE_QUEUE"};
		#$fastqSource        = $parameters{"fastqSource"};
		#$sampleFile         = $parameters{"SampleFile"};
		#$titleFile          = $parameters{"TitleFile"};
		#$standardNormalList = $parameters{"ListOfStandardNoramlsForGenotyping"};
		#$outdir             = $parameters{"Outdir"};
		#$datadir            = $parameters{"Datadir"};
		#$stdNormal          = $parameters{"stdNormal"};
		#$fof                = $parameters{"FOF"};
		$MAPQ               = $parameters{"MAPQ"};
		$BASEQ              = $parameters{"BASEQ"};
		#$poolName           = $parameters{"poolName"};
		#$projectName        = $parameters{"projectName"};
		$mvFiles            = $parameters{"moveFiles"};
		$process            = $parameters{"Process"};
		$prog               = $parameters{"Program"};
		$nprocessors        = $parameters{"NumberOfProcessors"};
		$OverallSupportingReads = $parameters{"OverallSupportingReads"};
		$OverallSupportingReadsHotspot = $parameters{"OverallSupportingReadsHotspot"};
		$SampleTumorSupportingReads = $parameters{"SampleTumorSupportingReads"};
		$SampleTumorSupportingReadsHotspot = $parameters{"SampleTumorSupportingReadsHotspot"};
		$SampleNormalSupportingReads = $parameters{"SampleNormalSupportingReads"};
		$SampleNormalSupportingReadsHotspot = $parameters{"SampleNormalSupportingReadsHotspot"};
		$OverallSupportingSplitReads = $parameters{"OverallSupportingSplitReads"};
		$OverallSupportingSplitReadsHotspot = $parameters{"OverallSupportingSplitReadsHotspot"};
		$LengthOfSV         = $parameters{"LengthOfSV"};
		$OverallMapq        = $parameters{"OverallMapq"};
		$OverallMapqHotspot = $parameters{"OverallMapqHotspot"};
		$SampleTumorGenotypeQualityFilter = $parameters{"SampleTumorGenotypeQualityFilter"};
		$SampleTumorGenotypeQualityFilterHotspot = $parameters{"SampleTumorGenotypeQualityFilterHotspot"};
		$SampleNormalGenotypeQualityFilter = $parameters{"SampleNormalGenotypeQualityFilter"};
		$SampleNormalGenotypeQualityFilterHotspot = $parameters{"SampleNormalGenotypeQualityFilterHotspot"};
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

# #####################################
# #####################################
# #This will help to call:
# #Somatic SVs: Delly
sub CallStructuralVariants {
	my @names       = keys(%BamSampleHash);
	my @FilterData  = ();

	my @notifyNames = ();
	tie( my %NormalPerFile,     'Tie::IxHash' );
	my $now = time;
	my $NormalUsed = $poolName . "_NormalUsedInSVcalling.txt";
	my @fileNames  = ();
	my ( $waitFileNames, $svOutdir );

	#Call Somatic SVs
	print "Started running Somatics Variant jobs on SGE at " . localtime() . "\n";

	my $count = 0;

	foreach my $file (@names) {
	my $fileSampleId = $BamSampleHash{$file};
	my $fileClass;
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

	print "Final2:Tumor->$file\nNormal->$normal\n\n";

	$NormalPerFile{$fileSampleId} = $normalSampleId;
	( $waitFileNames, $svOutdir ) = &RunDelly( $normal, $file, $normalSampleId, $fileSampleId, $outdir, $count );

	foreach my $waitName (@$waitFileNames) {
		push( @notifyNames, $waitName );
	}
	push( @FilterData, "$svOutdir,$normalSampleId,$fileSampleId" );
	$count++;
	}

	&WaitToFinish( $outdir, @notifyNames );
	open( NFH, ">", "$outdir/$NormalUsed" ) || die "Cannot open NormalUsedinSVFile:$outdir/$NormalUsed;$!\n";
	while ( my ( $key, $value ) = each(%NormalPerFile) ) {
		print NFH "$key\t$value\n";
	}
	close(NFH);
	$now = time - $now;
	print "Finished running Germline and Somatic Variant jobs on SGE at " . localtime() . "\n";
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
		if($genome eq "hg19")
		{
			$chr1 = "chr" . $chr1;
			$chr2 = "chr" . $chr2;
		}
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
            if($data[0] && $data[1])     
            {
                if( ($data[0] !~ /^NA$/i) && ($data[1] !~ /^NA$/i) )    
                {
                    $NormalTumorPair{$data[1]} = $data[0];
                }
                else
                {
                    print "Skipping invalid pair: $data[0] : $data[1]\n";
                }
            }
    }
    close PF;

    if(scalar(keys %NormalTumorPair) == 0) #no valid pairing found
    {
	print "No valid pairing found, the structural variation pipeline won't run.\n";
        exit 0;
    }
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
	my ($result_dir) = @_;
        my $src_dir = $result_dir . "/StrVarAnalysis";
	my $dst_dir = $result_dir . "/unfiltered";
	my $tmp_dir = $result_dir . "/tmp";
	`/bin/mkdir -p $dst_dir`;
	`/bin/mkdir -p $tmp_dir`;
        die "[ERROR]: Fail to create directory: $dst_dir\n" if(!-d $dst_dir);
	die "[ERROR]: Fail to create directory: $tmp_dir\n" if(!-d $tmp_dir);

	`find $src_dir -name \"*.vcf\" ! -name \"*stdfilter.vcf\" -exec mv {} $dst_dir \\;`;

	`find $result_dir -maxdepth 1 -type f ! -name \"*AllAnnotatedSVs*\" -exec mv {} $tmp_dir \\;`;

	`mv $src_dir $tmp_dir`;

}








