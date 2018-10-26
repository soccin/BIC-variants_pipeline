#!/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl
no if $] >= 5.017011, warnings => 'experimental::smartmatch'; #disable smartmatch warning for newer version of perl
use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use Cwd;
use Cwd 'abs_path';
use Tie::IxHash;
use lib "$Bin/lib";
use Schedule;
use Cluster;


my ($pre, $template_config, $result_dir, $berger, $bamlist, $patient, $title, $std_bamlist, $std_patient, $std_title, $std_covg, $genome, $scheduler, $priority_project, $priority_group);
my ($config, $custom_design_dir, $project_dir, $std_dir);
#my $lsf_queue = 0;

my $root_root;
BEGIN {($root_root) = $Bin =~ /(.*)\/.*\/.*/;}    #get the path above root


my $wrapperPath = $Bin . "/impact_dmp/wrapper/";
my $dmpCnvPath = $Bin . "/impact_dmp/bin/"; 
my $configFilePath = $Bin . "/impact_dmp/configuration/";
my $addCustomGenePath = $Bin . "/impact_dmp/ParseRefseqGene/ParseRefseqGene";
my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my $singularityParams = '';
my $singularityBind = '';
my $singularityenv_prepend_path = "";

my (
	#$sampleFile,                    #$titleFile,                     
	#$stdNormal,										 #$fof,
	#$runBQSR,                       #$mailID,
	#$list,                          ##$poolName,
	#$projectName,                   #$barcodeFile,
	#$process,                       #$standardNormalList,
	#$mvFiles,                       #$mergeDinucleotide,
	#$fastqSource,                   #$adaptorFile,
	#$TMPDIR,                        #$JAVA_1_6,
	$JAVA_1_7,                      #$ExonToGenCov,
	#$FPGenotypesScript,             #$GenerateMetricsFilesScript,
	$FP_genotypes,                  #$GATK_SomaticIndel,
	#$GATK,                          
	$Reference,
	$Reference_b37,			$Reference_hg19,
	$PICARD,						 						 #$Refseq,                         
	#$Mutect,                        #$filter_Mutect,
	#$filter_SomaticIndel,           
	$BaitInterval,
	$TargetInterval,                #$CompileMetrics,
	#$CAT,                           
	$PYTHON,
	$TrimGalore,                     $PERL,
	$BWA,                            $GeneInterval,
	$GeneIntervalAnn,                $GeneCoord,
	$TilingInterval,                 $TilingIntervalAnn,
	$FingerPrintInterval,           #$canonicalExonIntervalsFile,
	#$dbSNP,                         #$COSMIC,
	#$Mills_1000G_Indels,            #$dbSNP_bitset,
	#$dbProperties,                  #$Oncotator,
	#$Mutation_Assessor,             #$AnnotateAssessFilterVariants,
	#$LoessNormalization,            #$BestCopyNumber,
	#$NormVsNormCopyNumber,          #$StdNormalLoess,
	#$AllMetrics,                    #$SAMTOOLS,
	#$BEDTOOLS,                      #$GenotypeAllele,
	#$cosmicHotspotsVcf,             
	$GCBiasFile,
	#$HistNormDir,                   #$TNfreqRatio_MutectStdFilter,
	#$TNfreqRatio_SomIndelStdFilter, #$TNfreqRatioThreshold,
	#$ad_SomIndelStdFilter,          #$dp_SomIndelStdFilter,
	#$vf_SomIndelStdFilter,          #$ad_MutectStdFilter,
	#$dp_MutectStdFilter,            #$vf_MutectStdFilter,
	#$dp_snv,                        #$ad_snv,
	#$vf_snv,                        #$dp_snvHS,
	#$ad_snvHS,                      #$vf_snvHS,
	#$dp_indel,                      #$ad_indel,
	#$vf_indel,                      #$dp_indelHS,
	#$ad_indelHS,                    #$vf_indelHS,
	#$MAFthreshold,                  #$occurrencePercent,
	#$NormalVariantsVCF,             
	$queue,
	#$deleteIntermediateFiles,       #$MAPQ,
	#$BASQ,                          
	#$QSUB,
	$RHOME,                          $RLIBS,
	#$RSYNC,                         #$filterGenotypedVariants,
	#$createPatientVCFfile,          #$createCoverageFiles,
	#$ClinicalExons,                 #$HotspotMutations,
	#$IGVtools,                      #$translationFolder,
	$exonIntervalsFile,             #$validatedExons,
	#$coverageThreshold, 			 			 #$SVpipeline,
	#$stdNormalPath,
	$Bam2Fastq, 					 
	$REFSEQ,						 $REFSEQ_CANONICAL,
	$CYTOBAND,                       $UCSC_DBSNP,
	$FASTX_TOOLKIT,			 $maxium_readlength,
        $SINGULARITY
);

GetOptions ('pre=s' => \$pre,
		'genome=s' => \$genome,
		'template_config=s' => \$template_config,
		'result=s' => \$result_dir,
		'berger=s' => \$berger,
		'bamlist=s' => \$bamlist,
	    	'patient=s' => \$patient,
	    	'title=s' => \$title,
	    	'std_bamlist=s' => \$std_bamlist,
	    	'std_patient=s' => \$std_patient,
	    	'std_title=s' => \$std_title,
	    	'std_covg=s' => \$std_covg,
		'scheduler=s' => \$scheduler,
        	'priority_project=s' => \$priority_project,
        	'priority_group=s' => \$priority_group) or exit(1);

if(!$pre || !$genome || !$result_dir || !$berger || !$bamlist || (!$patient && !$title) ){
    print STDERR <<HELP;

    USAGE: ./exome_cnv.pl [options]

	* -pre <string>: output prefix (Project Name)
	* -genome <string>: genome version, b37 or hg19 
	* -template_config <string>: template config file, default is /configuration/template_exome_cnv.conf
	* -result <string>: result dir
	* -berger <string>: berger file
	* -bamlist <string>: bam list file for project samples
	* -patient <string>: patient file for project samples, used to create title file, mutually exclusive with -title
	* -title <string>: title file for project samples, mutually exclusive with -patient
	* -std_bamlist <string>: bam list file for standard normal samples, mutually exclusive with -std_covg
	* -std_patient <string>: patient file for standard normal samples, used to create standard normal title file, mutually exclusive with -std_title
	* -std_title <string>: title file for standard normal samples, mutually exclusive with -std_patient
	* -std_covg <string>: processed standard normal coverage file prefix, mutually exclusive with -std_bamlist and -std_title/-std_patient

HELP
exit;
}

die "[ERROR]: -scheduler/-priority_project/-priority_group unspecified\n" if(!$scheduler || !$priority_project || !$priority_group);
die "[ERROR]: Please specify the genome build by using -genome. Valid values are 'hg19' and 'b37'\n" if($genome ne "hg19" && $genome ne "b37");
die "[ERROR]: Please provide standard normal bam files with -std_bamlist/-std_patient or provide standard normal coverage file with -std_covg\n" if(!$std_bamlist && !$std_covg);
die "[ERROR]: -std_bamlist is provided, but -std_patient is missing\n" if($std_bamlist && (!$std_patient && !$std_title) );
die "[ERROR]: -patient and -title are mutually exclusive\n" if($patient && $title);
die "[ERROR]: -std_patient and -std_title are mutually exclusive\n" if($std_patient && $std_title);
die "[ERROR]: -std_covg is mutually exclusive with -std_bamlist and -std_title/-std_patient\n" if($std_covg && ($std_bamlist || $std_patient || $std_title) );


#CheckQueueSystem();

CheckFileDirectory();

GetConfiguration($template_config);

AddGeneInterval();

PrintConfig();

RunCopyNumber();



# check whether is SGE or LSF
#sub CheckQueueSystem
#{
#	if($ENV{'SGE_ROOT'})
#	{
#		print STDERR "[INFO]: Running jobs on SGE cluster\n";
#	}
#	else
#	{
#		print STDERR "[INFO]: Running jobs on LSF cluster\n";
#		$lsf_queue = 1;
#	}
#}



sub CheckFileDirectory
{
	if(!$template_config)
	{
		$template_config = "$configFilePath/template_exome_cnv.conf"; #default template config file
	}
	print STDERR "[INFO]: Config file in use: $template_config\n";


	$project_dir = "$result_dir/project_tmp/";
	`/bin/mkdir -p $project_dir`;
	die "[ERROR]: Fail to create result directory: $project_dir\n" if(!-d $project_dir);

	$std_dir = "$result_dir/standard_normal_tmp/";
	`/bin/mkdir -p $std_dir`;
	die "[ERROR]: Fail to create result directory: $std_dir\n" if(!-d $std_dir);
	die "[ERROR]: Berger file does not exist: $berger\n" if(!-e $berger);
	die "[ERROR]: Bam list file does not exist: $bamlist\n" if(!-e $bamlist); 
	die "[ERROR]: Patient file does not exist: $patient\n" if($patient && !-e $patient);
	die "[ERROR]: title file does not exist: $title\n" if($title && !-e $title); 
	die "[ERROR]: Template config file does not exist: $template_config\n" if(!-e $template_config);

	if($std_covg)
	{
		my $std_covg_title = $std_covg . "_title.txt";
		die "[ERROR]: Standard normal coverage title file does not exist: $std_covg_title\n" if(!-e $std_covg_title);
		my $std_covg_loess = $std_covg . "_ALL_intervalnomapqcoverage_loess.txt";
		die "[ERROR]: Standard normal loess coverage file does not exist: $std_covg_loess\n" if(!-e $std_covg_loess);
		$std_covg = abs_path($std_covg);
	}
	else
	{
		die "[ERROR]: Standard normal bam list file does not exist: $std_bamlist\n" if(!-e $std_bamlist); 
		die "[ERROR]: Standard normal patient file does not exist: $std_patient\n" if($std_patient && !-e $std_patient);
		die "[ERROR]: Standard normal title file does not exist: $std_title\n" if($std_title && !-e $std_title);
	}


	### get absolute path of files and directory
	$result_dir = abs_path($result_dir);
	$berger = abs_path($berger);
	$bamlist = abs_path($bamlist);
	if($patient)
	{
		$patient = abs_path($patient);
	}
	if($title)
	{
		$title = abs_path($title);
	}
	if($std_bamlist)
	{
		$std_bamlist = abs_path($std_bamlist);
	}
	if($std_patient)
	{
		$std_patient = abs_path($std_patient);
	}
	if($std_title)
	{
		$std_title = abs_path($std_title);
	}
	$template_config = abs_path($template_config);
	
	$custom_design_dir = "$result_dir/custom_design/";
	`/bin/mkdir -p $custom_design_dir`;
	die "[ERROR]: Fail to create impact plus directory: $custom_design_dir\n" if(!-d $custom_design_dir);

	$config = "$result_dir/${pre}_exome_cnv.conf";
}




#--GET CONFIGRATION DETAIL
sub GetConfiguration {
	my ($config_file) = @_;
	my @data = ();
	tie( my %config,     'Tie::IxHash' );
	tie( my %location,   'Tie::IxHash' );
	tie( my %version,    'Tie::IxHash' );
	tie( my %parameters, 'Tie::IxHash' );
	# change the default Input Record Separator.
	$/ = ">";
	open( CONFIG, "$config_file" ) or die("Cannot open $config_file. Error: $!\n");
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
				#iif($lsf_queue)
				#{
        			#	if($data[0] eq "QSUB")
        			#	{
                		#		$data[1] = "$wrapperPath/sge2lsf.pl";
        			#	}
        			#	else
        			#	{
                		#		$data[1] =~ s:/ifs/data/zeng/dmp/resources/:/ifs/work/zeng/dmp/resources/:g;
                		#		$data[1] =~ s:/ifs/data/dmp/data/:/ifs/work/zeng/dmp/data/data/:g;
       				#	}
				#}
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

	# Change back the Input Record Separator to default.
	$/ = "\n";
	eval {
		##Set Locations
		#$TMPDIR              = $location{"TMPDIR"};
		#$JAVA_1_6            = $location{"JAVA_1_6"};
		$JAVA_1_7            = $location{"JAVA_1_7"};
		$RHOME               = $location{"RHOME"};
		$RLIBS               = $location{"RLIBS"};
		#$ExonToGenCov        = $location{"ExonToGenCov"};
		#$FPGenotypesScript   = $location{"FPGenotypesScript"};
		$FP_genotypes        = $location{"FP_genotypes"};
		#$GATK_SomaticIndel   = $location{"GATK_SomaticIndel"};
		#$GATK                = $location{"GATK"};
		$Reference_b37           = $location{"Reference_b37"};	
		$Reference_hg19           = $location{"Reference_hg19"};
		#$Reference           = $location{"Reference"};
		#$Refseq              = $location{"Refseq"};
		$PICARD              = $location{"PICARD"};
		#$Mutect              = $location{"Mutect"};
		#$filter_Mutect       = $location{"filter_Mutect"};
		#$filter_SomaticIndel = $location{"filter_SomaticIndel"};
		$BaitInterval        = $location{"BaitInterval"};
		$TargetInterval      = $location{"TargetInterval"};
		#$CompileMetrics      = $location{"CompileMetrics"};
		#$SVpipeline			 = $location{"SVpipeline"};
		#$CAT                 = $location{"CAT"};
		$PYTHON              = $location{"PYTHON"};
		$TrimGalore          = $location{"TrimGalore"};
		$PERL                = $location{"PERL"};
		$BWA                 = $location{"BWA"};
		$GeneInterval        = $location{"GeneInterval"};
		$GeneIntervalAnn     = $location{"GeneIntervalAnn"};
		$GeneCoord           = $location{"GeneCoord"};
		$TilingInterval      = $location{"TilingInterval"};
		$TilingIntervalAnn   = $location{"TilingIntervalAnn"};
		$FingerPrintInterval = $location{"FingerPrintInterval"};
		#$dbSNP               = $location{"dbSNP"};
		#$COSMIC              = $location{"COSMIC"};
		#$Mills_1000G_Indels  = $location{"Mills_1000G_Indels"};
		#$dbSNP_bitset        = $location{"dbSNP_bitset"};
		#$dbProperties        = $location{"dbProperties"};
		#$Oncotator           = $location{"Oncotator"};
		#$Mutation_Assessor   = $location{"Mutation_Assessor"};
		#$AnnotateAssessFilterVariants = $location{"AnnotateAssessFilterVariants"};
		#$LoessNormalization   = $location{"LoessNormalization"};
		#$BestCopyNumber       = $location{"BestCopyNumber"};
		#$NormVsNormCopyNumber = $location{"NormVsNormCopyNumber"};
		#$StdNormalLoess       = $location{"StdNormalLoess"};
		#$NormalVariantsVCF    = $location{"NormalVariantsVCF"};
		#$AllMetrics           = $location{"AllMetrics"};
		#$SAMTOOLS             = $location{"SAMTOOLS"};
		#$BEDTOOLS             = $location{"BEDTOOLS"};
		#$GenotypeAllele       = $location{"GenotypeAllele"};
		#$cosmicHotspotsVcf    = $location{"CosmicHotspotVcf"};
		$GCBiasFile           = $location{"GCBiasFile"};
		#$HistNormDir          = $location{"HistNormDir"};
		#$barcodeFile          = $location{"BarcodeKey"};
		#$adaptorFile          = $location{"AdaptorKey"};
		#$QSUB                 = $location{"QSUB"};
		#$RSYNC                = $location{"RSYNC"};
		#$ClinicalExons        = $location{"clinicalExons"};
		#$HotspotMutations     = $location{"HotSpot_mutations"};
		#$IGVtools             = $location{"IGVtools"};
		#$translationFolder    = $location{"TranslationFolder"};
		$exonIntervalsFile = $location{"Canonical_Exon_Interval_table_with_aa"};
		#$canonicalExonIntervalsFile = $location{"Canonical_Exon_Interval_list"};    # this is for DoC
		#$validatedExons = $location{"Validated_Exons"};
		#$stdNormalPath = $location{"StandardNormalsDirectory"};
		$Bam2Fastq = $location{"Bam2Fastq"};
		$REFSEQ = $location{"REFSEQ"};
		$REFSEQ_CANONICAL = $location{"REFSEQ_CANONICAL"};
		$CYTOBAND = $location{"CYTOBAND"};
		$UCSC_DBSNP = $location{"UCSC_DBSNP"};
		$FASTX_TOOLKIT = $location{"FASTX_TOOLKIT"};
                $SINGULARITY = $location{"SINGULARITY"};

		##Set Parameters
		#$sampleFile  = $parameters{"SampleFile"};
		#$titleFile   = $parameters{"TitleFile"};
		#$fastqSource = $parameters{"FastqSource"};
		#$outdir = $parameters{"OutputDir"};
		#$datadir = $parameters{"RawDataDir"};
		#$MAPQ               = $parameters{"MAPQ"};
		#$BASQ               = $parameters{"BASQ"};
		#$stdNormal          = $parameters{"StdNormalForMutationCalling"};
		#$fof                = $parameters{"ListOfFiles"};
		#$poolName           = $parameters{"PoolName"};
		#$projectName        = $parameters{"ProjectName"};
		#$runBQSR            = $parameters{"RunBQSR"};
		#$process            = $parameters{"Process"};
		#$standardNormalList = $parameters{"ListOfStandardNoramlsForGenotyping"};
		#$mvFiles            = $parameters{"MoveFiles"};
		#$mergeDinucleotide  = $parameters{"MergeDinucleotide"};
		#$MAFthreshold       = $parameters{"MAFthreshold_AnnotationFilter"};
		#$TNfreqRatio_MutectStdFilter = $parameters{"TNfreqRatio_MutectStdFilter"};
		#$TNfreqRatio_SomIndelStdFilter = $parameters{"TNfreqRatio_SomIndelStdFilter"};
		#$TNfreqRatioThreshold = $parameters{"TNfreqRatio_AnnotationFilter"};
		#$ad_SomIndelStdFilter = $parameters{"AD_SomIndelSTDfilter"};
		#$dp_SomIndelStdFilter = $parameters{"DP_SomIndelSTDfilter"};
		#$vf_SomIndelStdFilter = $parameters{"VF_SomIndelSTDfilter"};
		#$ad_MutectStdFilter   = $parameters{"AD_MutectSTDfilter"};
		#$dp_MutectStdFilter   = $parameters{"DP_MutectSTDfilter"};
		#$vf_MutectStdFilter   = $parameters{"VF_MutectSTDfilter"};
		#$dp_snv               = $parameters{"minimumDPforSNV"};
		#$ad_snv               = $parameters{"minimumADforSNV"};
		#$vf_snv               = $parameters{"minimumVFforSNV"};
		#$dp_snvHS             = $parameters{"minimumDPforSNVhs"};
		#$ad_snvHS             = $parameters{"minimumADforSNVhs"};
		#$vf_snvHS             = $parameters{"minimumVFforSNVhs"};
		#$dp_indel             = $parameters{"minimumDPforINDEL"};
		#$ad_indel             = $parameters{"minimumADforINDEL"};
		#$vf_indel             = $parameters{"minimumVFforINDEL"};
		#$dp_indelHS           = $parameters{"minimumDPforINDELhs"};
		#$ad_indelHS           = $parameters{"minimumADforINDELhs"};
		#$vf_indelHS           = $parameters{"minimumVFforINDELhs"};
		#$occurrencePercent    = $parameters{"occurrencePercent"};
		#$coverageThreshold    = $parameters{"Coverage_threshold_darwin_report"};
		#$mailID               = $parameters{"Email"};
		$queue                = $parameters{"SGE_QUEUE"};
		#$deleteIntermediateFiles = $parameters{"DeleteIntermediateFiles"};
		$maxium_readlength    = $parameters{"MAXIMUM_READ_LENGTH"};
	};
	if($genome eq "b37")
	{
		$Reference = $Reference_b37;
	}
	else
	{
		$Reference = $Reference_hg19;
	}

        my %sinParams = (singularity_exec => "$SINGULARITY/singularity", singularity_image => "$root_root/variants_pipeline_singularity_prod.simg");
        $singularityParams = Schedule::singularityParams(%sinParams);
        $singularityBind = Schedule::singularityBind();

        $ENV{'SINGULARITYENV_PREPEND_PATH'} = $singularityenv_prepend_path;
        $ENV{'SINGULARITY_BINDPATH'} = $singularityBind;

	die ("[ERROR]: Did not find a variable in configuration file.Error: $@\n") if ($@);
	return ( \%version );
}



### add custom gene intervals to the pipeline
sub AddGeneInterval
{
	if(-e "$custom_design_dir/done.txt")
	{
		print "[INFO]: $custom_design_dir/done.txt exists, skip generating gene intervals\n" 
	}
	else
	{
		open(AL, ">$custom_design_dir/log.txt") or die "[ERROR]: Can't write to add custom gene log file $!\n";
		print STDERR "[INFO]: $addCustomGenePath --reference $Reference --refseq $REFSEQ --refseq_canonical $REFSEQ_CANONICAL --cytoband $CYTOBAND --output $custom_design_dir --gene_interval $GeneInterval --gene_interval_annotated $GeneIntervalAnn --gene_coord $GeneCoord --gc_bias $GCBiasFile --target_ilist $TargetInterval --bait_ilist $BaitInterval --canonical_aa $exonIntervalsFile --tiling_interval $TilingInterval --tiling_interval_annotated $TilingIntervalAnn --fp_genotype $FP_genotypes --fp_interval $FingerPrintInterval --ucsc_dnsnp $UCSC_DBSNP --custom_bed $berger --custom_only\n\n";
		print AL "[INFO]: $addCustomGenePath --reference $Reference --refseq $REFSEQ --refseq_canonical $REFSEQ_CANONICAL --cytoband $CYTOBAND --output $custom_design_dir --gene_interval $GeneInterval --gene_interval_annotated $GeneIntervalAnn --gene_coord $GeneCoord --gc_bias $GCBiasFile --target_ilist $TargetInterval --bait_ilist $BaitInterval --canonical_aa $exonIntervalsFile --tiling_interval $TilingInterval --tiling_interval_annotated $TilingIntervalAnn --fp_genotype $FP_genotypes --fp_interval $FingerPrintInterval --ucsc_dnsnp $UCSC_DBSNP --custom_bed $berger --custom_only\n\n";
		close AL;
		
		`$addCustomGenePath --reference $Reference --refseq $REFSEQ --refseq_canonical $REFSEQ_CANONICAL --cytoband $CYTOBAND --output $custom_design_dir --gene_interval $GeneInterval --gene_interval_annotated $GeneIntervalAnn --gene_coord $GeneCoord --gc_bias $GCBiasFile --target_ilist $TargetInterval --bait_ilist $BaitInterval --canonical_aa $exonIntervalsFile --tiling_interval $TilingInterval --tiling_interval_annotated $TilingIntervalAnn --fp_genotype $FP_genotypes --fp_interval $FingerPrintInterval --ucsc_dnsnp $UCSC_DBSNP --custom_bed $berger --custom_only 2>&1 | tee -a $custom_design_dir/log.txt /dev/stderr`;
		die "[ERROR]: Failed to generate custom interval files for IMPACT+\n" if(!-e "$custom_design_dir/done.txt");
	}
}




sub PrintConfig
{
### print config file
	open(TC, "$template_config") or die "[ERROR]: Can't open template config file $!\n";
	open(CF, ">$config") or die "[ERROR]: Can't write config file $!\n";
	while(<TC>){
		#if($lsf_queue)
 		#{
 		#	if($_ =~ /^QSUB =.*/)
 		#	{
 		#		$_ = "QSUB = $wrapperPath/sge2lsf.pl\n";
 		#	}
 		#	else
 		#	{
 		#		$_ =~ s:/ifs/data/zeng/dmp/resources/:/ifs/work/zeng/dmp/resources/:g;
 		#		$_ =~ s:/ifs/data/dmp/data/:/ifs/work/zeng/dmp/data/data/:g;
 		#	}
 		#}
		if($_ =~ /^Reference_b37 =.*/)
                {       
                        print CF "Reference = $Reference\n";
                }
		elsif($_ =~ /^Reference_hg19 =.*/)
                {
                	# do not print
		}
		elsif($_ =~ /^LoessNormalization =.*/)
		{
			print CF "LoessNormalization = $dmpCnvPath/loessnormalize_nomapq.R\n";
		}
		elsif($_ =~ /^BestCopyNumber =.*/)
		{
			print CF "BestCopyNumber = $dmpCnvPath/copynumber_tm.batchdiff.bic.R\n";
		}
		elsif($_ =~ /^NormVsNormCopyNumber =.*/)
		{
			print CF "NormVsNormCopyNumber = $dmpCnvPath/copynumber_testclass.batchdiff.bic.R\n";
		}
		elsif($_ =~ /^TMPDIR =.*/)
		{
			print CF "TMPDIR = /scratch/$uID\n";
		}
		elsif($_ =~ /^ProjectName =.*/)
		{	
			print CF "ProjectName = $pre\n";
		}
		elsif($_ =~ /^PoolName =.*/)
		{	
			print CF "PoolName = $pre\n";
		}
		elsif($_ =~ /^OutputDir =.*/)
		{	
			print CF "OutputDir = $result_dir\n";
		}
		elsif($_ =~ /^GeneInterval =.*/)
		{
			print CF "GeneInterval = $custom_design_dir/combined_gene_intervals.list\n";
		}
		elsif($_ =~ /^GeneIntervalAnn =.*/)
		{
			print CF "GeneIntervalAnn = $custom_design_dir/combined_gene_intervals.list.annotated\n";			
		}	
		elsif($_ =~ /^GeneCoord =.*/)
		{
			print CF "GeneCoord = $custom_design_dir/combined_gene_coords.txt\n";			
		}	
		elsif($_ =~ /^GCBiasFile =.*/)
		{
			print CF "GCBiasFile = $custom_design_dir/combined_gc_bias.txt\n";		
		}
		elsif($_ =~ /^TargetInterval =.*/)
		{
			print CF "TargetInterval = $custom_design_dir/combined_target.ilist\n";		
		}
		elsif($_ =~ /^BaitInterval =.*/)
		{
			print CF "BaitInterval = $custom_design_dir/combined_bait.ilist\n";		
		}
		elsif($_ =~ /^Canonical_Exon_Interval_table_with_aa =.*/)
		{
			print CF "Canonical_Exon_Interval_table_with_aa = $custom_design_dir/combined_canonical_exon_with_aa.list\n";		
		}
		elsif($_ =~ /^Canonical_Exon_Interval_list =.*/)
		{
			print CF "Canonical_Exon_Interval_list = $custom_design_dir/combined_canonical_exon.list\n";		
		}
		elsif($_ =~ /^FP_genotypes =.*/)
		{
			print CF "FP_genotypes = $custom_design_dir/combined_fp_tiling_genotypes.txt\n";		
		}
		elsif($_ =~ /^TilingInterval =.*/)
		{
			print CF "TilingInterval = $custom_design_dir/combined_tiling_interval.list\n";		
		}
		elsif($_ =~ /^TilingIntervalAnn =.*/)
		{
			print CF "TilingIntervalAnn = $custom_design_dir/combined_tiling_interval.list.annotated\n";		
		}
		elsif($_ =~ /^FingerPrintInterval =.*/)
		{
			print CF "FingerPrintInterval = $custom_design_dir/combined_fp_tiling_interval.list\n";		
		}
		else
		{
 			print CF $_;
		}
	}
	close TC;
	close CF;
}



sub RunCopyNumber
{
	my $hold_job_id = "$pre\_$uID\_NormalizedCoverage";

	my $extra_para;
	if($title)
	{
	 	$extra_para = "-title $title";
	}
	else
	{
	 	$extra_para = "-patient $patient";
	}
	if($genome eq "hg19")
	{
		$extra_para .= " -hg19"
	}

	my %addParams = (scheduler => "$scheduler", runtime => "60", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", work_dir => "$project_dir", iounits => "1");
	my $additionalParams = Schedule::additionalParams(%addParams);

	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_NormalizedCoverage", cpu => "1", mem => "8", cluster_out => "$pre\_$uID\_NormalizedCoverage.stdout", cluster_error => "$pre\_$uID\_NormalizedCoverage.stderr");
	my $standardParams = Schedule::queuing(%stdParams);

	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $PERL $wrapperPath/RunNormalizedCoverage.pl -pre $pre -config $config -bamlist $bamlist -result $project_dir $extra_para -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;

	if(!$std_covg)
	{
 		my $std_extra_para;
	 	if($std_title)
	 	{
	 		$std_extra_para = "-title $std_title";
	 	}
	 	else
	 	{
	 		$std_extra_para = "-patient $std_patient";
	 	}
        	if($genome eq "hg19")
        	{
        	        $std_extra_para .= " -hg19"
	        }

		%addParams = (scheduler => "$scheduler", runtime => "60", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", work_dir => "$std_dir", iounits => "1");
		$additionalParams = Schedule::additionalParams(%addParams);
		%stdParams = (scheduler => "$scheduler", job_name => "$pre\_std\_$uID\_NormalizedCoverage", cpu => "1", mem => "8", cluster_out => "pre\_std\_$uID\_NormalizedCoverage.stdout", cluster_error => "$pre\_std\_$uID\_NormalizedCoverage.stderr");
		$standardParams = Schedule::queuing(%stdParams);

 		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $PERL $wrapperPath/RunNormalizedCoverage.pl -pre $pre\_std -config $config -bamlist $std_bamlist -result $std_dir $std_extra_para -add_barcode_prefix -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;
	
 		$std_covg = $std_dir . "$pre\_std_fixed";
 		$hold_job_id = $hold_job_id . ",$pre\_std\_$uID\_NormalizedCoverage";
 	}

	$ENV{"R_LIBS_USER"}=$RLIBS;

	%addParams = (scheduler => "$scheduler", runtime => "60", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", work_dir => "$project_dir", iounits => "1");
	$additionalParams = Schedule::additionalParams(%addParams);

	%stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_TumorCopyNumber", job_hold => "$hold_job_id", cpu => "1", mem => "8", cluster_out => "$pre\_$uID\_TumorCopyNumber.stdout", cluster_error => "$pre\_$uID\_TumorCopyNumber.stderr");
	$standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams "$singularityParams $RHOME/R --slave --vanilla --args $pre $std_covg $custom_design_dir/combined_gene_intervals.list.annotated $custom_design_dir/combined_tiling_interval.list.annotated MIN < $dmpCnvPath/copynumber_tm.batchdiff.bic.R"`;


        %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_NormalCopyNumber", job_hold => "$hold_job_id", cpu => "1", mem => "8", cluster_out => "$pre\_$uID\_NormalCopyNumber.stdout", cluster_error => "$pre\_$uID\_NormalCopyNumber.stderr");
        $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams "$singularityParams $RHOME/R --slave --vanilla --args $pre $std_covg $custom_design_dir/combined_gene_intervals.list.annotated $custom_design_dir/combined_tiling_interval.list.annotated Normal MIN < $dmpCnvPath/copynumber_testclass.batchdiff.bic.R"`;


	my $post_processing_extra_para;
	if($genome eq "hg19")
	{
		$post_processing_extra_para = "-hg19"
	}

	%stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CopyNumberPostProcessing", job_hold => "$pre\_$uID\_TumorCopyNumber,$pre\_$uID\_NormalCopyNumber", cpu => "1", mem => "8", cluster_out => "$pre\_$uID\_CopyNumberPostProcessing.stdout", cluster_error => "$pre\_$uID\_CopyNumberPostProcessing.stderr");
        $standardParams = Schedule::queuing(%stdParams);
	
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $PERL $singularityParams $wrapperPath/exome_cnv_post_processing.pl -pre $pre -input $project_dir -output $result_dir -perl $PERL $post_processing_extra_para`;

	print "[INFO]: Submitted copynumber analysis job to cluster\n";

}



