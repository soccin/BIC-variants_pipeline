#!/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl
no if $] >= 5.017011, warnings => 'experimental::smartmatch'; #disable smartmatch warning for newer version of perl
use strict;
use Getopt::Long;
use IO::File;
use List::Util qw(max sum);
use Tie::IxHash;
use Vcf;
use File::Basename;
use MSKCC_DMP_Logger;
use Data::Dumper;
use FindBin qw($Bin);
use version;
my $logger = MSKCC_DMP_Logger->get_logger('IMPACT_Pipeline_Logger');
$logger->start_local();
undef $!;
undef $@;
my $now = time;
my ($root) = $Bin =~ /(.*\/).*/;    #get the path above bin
my $dmpCnvPath = $root . "/bin/"; 
my $root_root;
BEGIN {($root_root) = $Bin =~ /(.*)\/.*\/.*/;}    #get the path above root
use lib "$root_root/lib";
use Schedule;
use Cluster;


my @projectBamFiles = ();
my %sampleToBarcode = ();
my %barcodeToFixedBarcode = ();
my $add_barcode_prefix = 0;
my @blood_type = ("Blood","Frozen","Fresh","OCT");
my @ffpe_type = ("FFPE", "Cytology");

my $singularityParams = '';
my $singularityBind = '';
my $singularityenv_prepend_path = "";


my ($pre, $config, $bamlist, $patient, $title, $result_dir, $hg19, $scheduler, $priority_project, $priority_group);
if (
	 @ARGV < 1
	 or !GetOptions(
					 'pre=s' => \$pre,
					 'hg19' =>\$hg19,
					 'config=s'  => \$config,
					 'bamlist=s' => \$bamlist,
					 'patient=s' => \$patient,
					 'title=s' => \$title,
					 'result=s'  => \$result_dir,
					 'add_barcode_prefix' => \$add_barcode_prefix,
					 'scheduler=s' => \$scheduler,
                			 'priority_project=s' => \$priority_project,
                			 'priority_group=s' => \$priority_group)
  )
{
	Usage();
}

die "[ERROR]: argument error" if(!$pre || !$config || !$bamlist || (!$patient && !$title) || !$result_dir);
die "[ERROR]: -patient and -title are mutually exclusive" if($patient && $title);


my (
	#$sampleFile,                    #$titleFile,                     
	#$stdNormal,										 #$fof,
	#$runBQSR,                       #$mailID,
	#$list,                          ##$poolName,
	#$projectName,                   #$barcodeFile,
	#$process,                       #$standardNormalList,
	#$mvFiles,                       #$mergeDinucleotide,
	#$fastqSource,                   #$adaptorFile,
	$TMPDIR,                        #$JAVA_1_6,
	$JAVA_1_7,                      #$ExonToGenCov,
	#$FPGenotypesScript,             #$GenerateMetricsFilesScript,
	#$FP_genotypes,                  #$GATK_SomaticIndel,
	$GATK,                          
	$Reference,
	#$PICARD,						 						 #$Refseq,                         
	#$Mutect,                        #$filter_Mutect,
	#$filter_SomaticIndel,           
	#$BaitInterval,
	#$TargetInterval,                #$CompileMetrics,
	#$CAT,                           
	#$PYTHON,
	#$TrimGalore,                     
	$PERL,
	#$BWA,                            
	$GeneInterval,
	$GeneIntervalAnn,                $GeneCoord,
	$TilingInterval,                 $TilingIntervalAnn,
	$FingerPrintInterval,           #$canonicalExonIntervalsFile,
	#$dbSNP,                         #$COSMIC,
	#$Mills_1000G_Indels,            #$dbSNP_bitset,
	#$dbProperties,                  #$Oncotator,
	#$Mutation_Assessor,             #$AnnotateAssessFilterVariants,
	$LoessNormalization,            #$BestCopyNumber,
	#$NormVsNormCopyNumber,          #$StdNormalLoess,
	#$AllMetrics,                    
	$SAMTOOLS,
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
	$BASQ,                          
	#$QSUB,
	$RHOME,                         $RLIBS,
	#$RSYNC,                         #$filterGenotypedVariants,
	#$createPatientVCFfile,          #$createCoverageFiles,
	#$ClinicalExons,                 #$HotspotMutations,
	#$IGVtools,                      #$translationFolder,
	#$exonIntervalsFile,             #$validatedExons,
	#$coverageThreshold, 			 			 #$SVpipeline,
	#$stdNormalPath,
	#$Bam2Fastq, 					 
	#$REFSEQ,						 $REFSEQ_CANONICAL,
	#$CYTOBAND,                       $UCSC_DBSNP,
	#$FASTX_TOOLKIT,			 		$maxium_readlength,
	#$ReferenceHG19
        $SINGULARITY
);


GetConfiguration($config);

MakeCSH($result_dir);

if($patient)
{
	ProcessPatientFile($patient, $result_dir);
}
else
{
	ProcessTitleFile($title, $result_dir);
}

ReadBamList($bamlist, $result_dir);

RunNormalizedCoverage(\@projectBamFiles, $result_dir);




#How to use the script.
sub Usage
{
	print "\nUsage : RunNormalizedCoverage.pl [options]
             * -pre <string>: project prefix
             * -config <string>: config file
             * -bamlist <string>: bam list file contains Sample_ID and Bam_File_Path
             * -patient <string>: project patient file
             * -result <string>: path to write the output
	\n";
	exit();
}
###################################################

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
		$TMPDIR              = $location{"TMPDIR"};
		#$JAVA_1_6            = $location{"JAVA_1_6"};
		$JAVA_1_7            = $location{"JAVA_1_7"};
		$RHOME               = $location{"RHOME"};
		$RLIBS               = $location{"RLIBS"};
		#$ExonToGenCov        = $location{"ExonToGenCov"};
		#$FPGenotypesScript   = $location{"FPGenotypesScript"};
		#$FP_genotypes        = $location{"FP_genotypes"};
		#$GATK_SomaticIndel   = $location{"GATK_SomaticIndel"};
		$GATK                = $location{"GATK"};
                $SINGULARITY               = $location{"SINGULARITY"};
		$Reference           = $location{"Reference"};
		#$ReferenceHG19       = $location{"ReferenceHG19"};
		#$Refseq              = $location{"Refseq"};
		#$PICARD              = $location{"PICARD"};
		#$Mutect              = $location{"Mutect"};
		#$filter_Mutect       = $location{"filter_Mutect"};
		#$filter_SomaticIndel = $location{"filter_SomaticIndel"};
		#$BaitInterval        = $location{"BaitInterval"};
		#$TargetInterval      = $location{"TargetInterval"};
		#$CompileMetrics      = $location{"CompileMetrics"};
		#$SVpipeline			 = $location{"SVpipeline"};
		#$CAT                 = $location{"CAT"};
		#$PYTHON              = $location{"PYTHON"};
		#$TrimGalore          = $location{"TrimGalore"};
		$PERL                = $location{"PERL"};
		#$BWA                 = $location{"BWA"};
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
		$LoessNormalization   = $location{"LoessNormalization"};
		#$BestCopyNumber       = $location{"BestCopyNumber"};
		#$NormVsNormCopyNumber = $location{"NormVsNormCopyNumber"};
		#$StdNormalLoess       = $location{"StdNormalLoess"};
		#$NormalVariantsVCF    = $location{"NormalVariantsVCF"};
		#$AllMetrics           = $location{"AllMetrics"};
		$SAMTOOLS             = $location{"SAMTOOLS"};
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
		#$exonIntervalsFile = $location{"Canonical_Exon_Interval_table_with_aa"};
		#$canonicalExonIntervalsFile = $location{"Canonical_Exon_Interval_list"};    # this is for DoC
		#$validatedExons = $location{"Validated_Exons"};
		#$stdNormalPath = $location{"StandardNormalsDirectory"};
		#$Bam2Fastq = $location{"Bam2Fastq"};
		#$REFSEQ = $location{"REFSEQ"};
		#$REFSEQ_CANONICAL = $location{"REFSEQ_CANONICAL"};
		#$CYTOBAND = $location{"CYTOBAND"};
		#$UCSC_DBSNP = $location{"UCSC_DBSNP"};
		#$FASTX_TOOLKIT = $location{"FASTX_TOOLKIT"};

		##Set Parameters
		#$sampleFile  = $parameters{"SampleFile"};
		#$titleFile   = $parameters{"TitleFile"};
		#$fastqSource = $parameters{"FastqSource"};
		#$outdir = $parameters{"OutputDir"};
		#$datadir = $parameters{"RawDataDir"};
		#$MAPQ               = $parameters{"MAPQ"};
		$BASQ               = $parameters{"BASQ"};
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
		#$maxium_readlength    = $parameters{"MAXIMUM_READ_LENGTH"};
	};
        my %sinParams = (singularity_exec => "$SINGULARITY/singularity", singularity_image => "$root_root/variants_pipeline_singularity_prod.simg");
        $singularityParams = Schedule::singularityParams(%sinParams);
        $singularityBind = Schedule::singularityBind();

        $ENV{'SINGULARITYENV_PREPEND_PATH'} = $singularityenv_prepend_path;
        $ENV{'SINGULARITY_BINDPATH'} = $singularityBind;

	die ("[ERROR]: Did not find a variable in configuration file.Error: $@\n") if ($@);
	return ( \%version );
}


#--Make Notification file
sub MakeCSH
 {
        my($outdir) = @_;
        my $filename = $outdir . "/Notify.csh";
        my $ntmp = new IO::File(">$filename");
        print $ntmp "#!/bin/csh\n";
        print $ntmp "#Notification File\n";
        print $ntmp "echo"," This is Done","\n";
        $ntmp->close();
        eval{`chmod +x $filename`;};
		if($@){$logger->fatal("Cannot change permissions for $filename. Error:$@");exit(1);}
        return;
}



###################################################
#--Waiting for the process to finish
sub WaitToFinish
{
	my ( $outdir, @waitfilenames ) = @_;
	$logger->info("Waiting for the Process to finish...");
	foreach my $wfile (@waitfilenames)
	{
		next if ( $wfile eq "NULL" );
		sleep 10 while ( !( -e "$outdir/$wfile" ) );

		#print "$outdir/$wfile\n";
		while ( -e "$outdir/$wfile" )
		{

			#print "$outdir/$wfile\n";
			open( FH, "<", "$outdir/$wfile" );
			while (<FH>)
			{
				if ( $_ =~ /This is Done/ig )
				{

					#print "\nFinished: $wfile\n";
					last;
				} else
				{
					wait;
				}
			}
			last;
		}
		close(FH);
	}
	foreach my $wfile (@waitfilenames)
	{
		next if ( $wfile eq "NULL" );
		eval { `rm $outdir/$wfile`; };
		if ($@) { $logger->warn("Cannot remove $outdir/$wfile. Error:$@"); }
	}
	return;
}
###################################################



###################################################
#--Check Output Files
sub CheckOutputFiles
{
	my ( $outdir, $prog, @Outputfilenames ) = @_;
	$logger->info("Checking Output Files");
	$prog = "NULL" if ( !$prog );
	if ( ( $prog eq "STDfilter" ) or ( $prog eq "RawMutationFiles" ) )
	{
		foreach my $ofile (@Outputfilenames)
		{
			if ( ( $ofile =~ /\// ) and ( !( $ofile =~ /,/ ) ) )
			{
				if ( -e "$ofile" )
				{
					next;
				} else
				{
					$logger->fatal(    "CheckFile:Please check, $ofile file was not created;Something went wrong here."
					);
					exit(1);
				}
			}
			if ( -e "$outdir/$ofile" )
			{
				next;
			} else
			{
				$logger->fatal( "CheckFile:Please check, $outdir/$ofile file was not created;Something went wrong here."
				);
				exit(1);
			}
		}
	} else
	{
		foreach my $ofile (@Outputfilenames)
		{
			if ( ( $ofile =~ /\// ) and ( !( $ofile =~ /,/ ) ) )
			{
				if ( ( -e "$ofile" ) and ( ( -s "$ofile" ) != 0 ) )
				{
					next;
				} else
				{
					$logger->fatal(    "CheckFile:Please check, $ofile file was not created or the file size is zero; Something went wrong here."
					);
					exit(1);
				}
			} elsif ( ( $ofile =~ /,/ ) and ( !( $ofile =~ /\// ) ) )
			{
				my ( $file1, $file2 ) = split( ",", $ofile );
				if (     ( -e "$outdir/$file1" )
					 and ( ( -s "$outdir/$file1" ) != 0 )
					 and ( -e "$outdir/$file2" )
					 and ( ( -s "$outdir/$file2" ) != 0 ) )
				{
					next;
				} else
				{
					$logger->fatal(    "CheckFile:Please check,$file1 & $file2 files was not created or the file size is zero; Something went wrong here."
					);
					exit(1);
				}
			} elsif ( ( $ofile =~ /,/ ) and ( $ofile =~ /\// ) )
			{
				my ( $file1, $file2 ) = split( ",", $ofile );
				if (     ( -e "$file1" )
					 and ( ( -s "$file1" ) != 0 )
					 and ( -e "$file2" )
					 and ( ( -s "$file2" ) != 0 ) )
				{
					next;
				} else
				{
					$logger->fatal(    "CheckFile:Please check,$file1 & $file2 files was not created or the file size is zero; Something went wrong here."
					);
					exit(1);
				}
			} else
			{
				if (     ( -e "$outdir/$ofile" )
					 and ( ( -s "$outdir/$ofile" ) != 0 ) )
				{
					next;
				} else
				{
					$logger->fatal(    "CheckFile:Please check, $outdir/$ofile file was not created or the file size is zero;Something went wrong here."
					);
					exit(1);
				}
			}
		}
	}
	$logger->info("Finished Checking Output Files");
	return;
}
###################################################


#sort by barcode name:
sub lowestNumber
{
	my $files = shift;
	my @filenames = split( ",", $files );
	my ($number) = $filenames[0] =~ m/.*_bc(\d{1,4})_.*/g;
	return $number;
}


#sub CheckBamHeaderChr
#{
#	my ($filename) = @_;
#	my $first_seq_name = `$SAMTOOLS view -H $filename | grep \@SQ | head -n 1`;
#	chomp $first_seq_name;
#	if ($first_seq_name =~ /chr/)
#	{
#		print "Bam file contains hg19 chromomesome header: $first_seq_name\n";
#		$hg19 = 1;
#	}
#}




##### read bam list file
sub ReadBamList
{
	my ($filename, $outdir) = @_;
	open(BL, "$filename") or die "can't open bam list file $filename $!";
	while(<BL>){
		chomp;
		my @data = split(/\s+/);
		my $sample_id = $data[0];
		die "Could not find $sample_id in patient file: $patient\n" if(!exists $sampleToBarcode{$sample_id});
		my $barcode_id = $sampleToBarcode{$sample_id};
		my $bam_file = $data[1];
		my $bai_file1 = $bam_file . ".bai";
		my $bai_file2 = $bam_file;
		$bai_file2 =~ s/^(.*)\.bam$/$1\.bai/; 
		
		my $bam_file_ln = $sample_id . "_" . $barcode_id . "_" . $pre . "_L000_mrg_cl_aln_srt_MD_IR_BR.bam";
		my $bai_file_ln = $bam_file_ln . ".bai";

		die "BAM file doest not exist: $bam_file\n" if(!-e $bam_file);
		`ln -s $bam_file $outdir/$bam_file_ln`;
		if(-e $bai_file1)
		{
			`ln -s $bai_file1 $outdir/$bai_file_ln`;
		}
		elsif(-e $bai_file2)
		{
			`ln -s $bai_file2 $outdir/$bai_file_ln`;
		}
		else
		{
			die "Neither BAI file exists: $bai_file1, $bai_file2\n";
		}
		push( @projectBamFiles, $bam_file_ln);
	}
	close BL;
}



### read information from patient file
sub ProcessPatientFile
{
	my ($patient_file, $outdir) = @_;
	my $titlefile = $pre . "_title.txt";
	my $barcode_count = 0;

	open(PF, "$patient_file") or die "[ERROR]: Can't open patient file $patient_file $!\n";
	my $pf_header = <PF>;
	my @header_data = split(/\s+/, $pf_header);
	die "[ERROR]: Sex information not found in patient file.\n" if(scalar(@header_data) < 11 || $header_data[10] ne "Sex");
	
	open(TF, ">$outdir/$titlefile") or die "[ERROR]: Can't write to $outdir/$titlefile $!\n";
	my @titleFileHeader = ("Barcode", "Pool", "Sample_ID", "Collab_ID", "Patient_ID","Class", "Sample_type", "Input_ng", "Library_yield", "Pool_input", "Bait_version", "Gender");
	print TF join("\t", @titleFileHeader), "\n";

	while(<PF>){
		chomp;
		my @data = split(/\s+/);
		if($data[4] eq "Normal" || $data[4] eq "Tumor" || $data[4] eq "PoolNormal" )
		{
			my $sample_id = $data[1];
			my $collab_id = $data[2];
			my $patient_id = $data[3];
			my $tumor_type = $data[4];
			my $sample_type = $data[5];
			my $input_ng = $data[6];
			my $library_yield = $data[7];
			my $pool_input = $data[8];
			my $bait_version = $data[9];
			my $gender;
			if($data[10] eq "M")
			{
				$gender = "Male";
			}
			elsif($data[10] eq "F")
			{
				$gender = "Female";
			}
			else
			{
				$gender = "-";
			}
			
			$barcode_count ++;
			my $barcode_str = "bc" . sprintf( "%04d", $barcode_count);
			$sampleToBarcode{$sample_id} = $barcode_str;

			my $fixed_barcode_str;
			if($sample_type ~~ @blood_type)
			{
				$fixed_barcode_str = "R0." . $barcode_str;
			}
			elsif($sample_type ~~ @ffpe_type)
			{
				$fixed_barcode_str = "FFPE." . $barcode_str;
			}
			else
			{
				$fixed_barcode_str = $barcode_str;
			}
			$barcodeToFixedBarcode{$barcode_str} = $fixed_barcode_str;
			print TF "$barcode_str\t$pre\t$sample_id\t$collab_id\t$patient_id\t$tumor_type\t$sample_type\t$input_ng\t$library_yield\t$pool_input\t$bait_version\t$gender\n";
		}
	}
	close PF;
	close TF;
}



### read information from patient file
sub ProcessTitleFile
{
	my ($src_title_file, $outdir) = @_;
	my $titlefile = $pre . "_title.txt";
	my $barcode_count = 0;

	open(STF, "$src_title_file") or die "[ERROR]: Can't open source title file $src_title_file $!\n";
	my $stf_header = <STF>;
	chomp $stf_header;
	
	open(TF, ">$outdir/$titlefile") or die "[ERROR]: Can't write to $outdir/$titlefile $!\n";
	print TF "$stf_header\n";

	while(<STF>){
		chomp;
		my @data = split(/\s+/);
			my $sample_id = $data[2];
			my $original_barcode = $data[0];
			$barcode_count ++;
			my $barcode_str = "bc" . sprintf( "%04d", $barcode_count);
			$data[0] = $barcode_str;
			$sampleToBarcode{$sample_id} = $barcode_str;
			print TF join("\t", @data), "\n";

			my $fixed_barcode_str;
			if(index($original_barcode, ".") != -1)
			{
				my @prefix_data = split(/\./, $original_barcode);
				$fixed_barcode_str = $prefix_data[0] . "." . $barcode_str;
			}
			else
			{
				$fixed_barcode_str = $barcode_str;
			}
			$barcodeToFixedBarcode{$barcode_str} = $fixed_barcode_str;
	}
	close STF;
	close TF;
}


#Metrics Calculations for Bam
sub RunDepthOfCoverage
{
	my ( $bamFile, $outdir, $id ) = @_;
	my ($basename) = $bamFile =~ /(.*)\.bam/;
	my $gene_nomapqCoverage         = $basename . ".gene_nomapq.covg";
	my $tiling_nomapqCoverage = $basename . ".tiling_nomapq.covg";

	my @notifynames   = ();
	my @metricsOutput = ();

	my $cur_gene_interval = $GeneInterval;
	my $cur_tiling_interval = $TilingInterval;
	if($hg19)
	{
		$cur_gene_interval = $cur_gene_interval . ".chr.list";
		$cur_tiling_interval = $cur_tiling_interval . ".chr.list";
	}

	#Gene_Nomapq Coverage
	if (
		     ( -e "${outdir}/${gene_nomapqCoverage}.sample_interval_summary" ) and
		 ( ( -s "${outdir}/${gene_nomapqCoverage}.sample_interval_summary" ) != 0 ))
	{
		$logger->info("File:\n${gene_nomapqCoverage}.sample_interval_summary\n they exists and process will not run to make \".gene_nomapq.covg\.sample_interval_summary\" file.");
		push( @notifynames, "NULL" );
		push( @metricsOutput, "${gene_nomapqCoverage}.sample_interval_summary" );
	} else
	{
		eval {
	        	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", work_dir => "$outdir", iounits => "1");
        		my $additionalParams = Schedule::additionalParams(%addParams);

        		my %stdParams = (scheduler => "$scheduler", job_name => "GeneNomapqCoverage.$id.$$", cpu => "1", mem => "8", cluster_out => "GeneNomapqCoverage.$id.$$.stdout", cluster_error => "GeneNomapqCoverage.$id.$$.stderr");
        		my $standardParams = Schedule::queuing(%stdParams);		

			my $cmd = "$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $JAVA_1_7 -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $gene_nomapqCoverage -L $cur_gene_interval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput --includeRefNSites";
			$logger->debug("COMMAND: $cmd");
			`$cmd`;

                        %stdParams = (scheduler => "$scheduler", job_name => "NotifyGeneNomapqCoverage.$id.$$", job_hold => "GeneNomapqCoverage.$id.$$", cpu => "1", mem => "2", cluster_out => "NotifyGeneNomapqCoverage.$id.$$.stat", cluster_error => "NotifyGeneNomapqCoverage.$id.$$.stderr");
                        $standardParams = Schedule::queuing(%stdParams); 
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $outdir/Notify.csh`;
		};
		if ($@)
		{
			$logger->fatal("GeneNomapqCoverage:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyGeneNomapqCoverage.$id.$$.stat" );
		push( @metricsOutput, "${gene_nomapqCoverage}.sample_interval_summary" );
	}

	#Tiling_Nomapq Coverage
	if (
		     ( -e "${outdir}/${tiling_nomapqCoverage}.sample_interval_summary" ) and
		 ( ( -s "${outdir}/${tiling_nomapqCoverage}.sample_interval_summary" ) != 0 ))
	{
		$logger->info("File:\n${tiling_nomapqCoverage}.sample_interval_summary\n they exists and process will not run to make \".tiling_nomapq.covg\.sample_interval_summary\" file.");
		push( @notifynames, "NULL" );
		push( @metricsOutput, "${tiling_nomapqCoverage}.sample_interval_summary" );
	} else
	{
		eval {
			my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", work_dir => "$outdir", iounits => "1");
                        my $additionalParams = Schedule::additionalParams(%addParams);

                        my %stdParams = (scheduler => "$scheduler", job_name => "TilingNomapqCoverage.$id.$$", cpu => "1", mem => "8", cluster_out => "TilingNomapqCoverage.$id.$$.stdout", cluster_error => "TilingNomapqCoverage.$id.$$.stderr");
                        my $standardParams = Schedule::queuing(%stdParams);

			my $cmd = "$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $JAVA_1_7 -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $GATK -T DepthOfCoverage -R $Reference -I $bamFile -o $tiling_nomapqCoverage -L $cur_tiling_interval  -rf BadCigar -mmq 0 -mbq $BASQ -omitLocusTable -omitSampleSummary -omitBaseOutput --includeRefNSites";
			$logger->debug("COMMAND: $cmd");
			`$cmd`;


                        %stdParams = (scheduler => "$scheduler", job_name => "NotifyTilingNomapqCoverage.$id.$$", job_hold => "TilingNomapqCoverage.$id.$$", cpu => "1", mem => "2", cluster_out => "NotifyTilingNomapqCoverage.$id.$$.stat", cluster_error => "NotifyTilingNomapqCoverage.$id.$$.stderr");
                        $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $outdir/Notify.csh`;
		};
		if ($@)
		{
			$logger->fatal("TilingNomapqCoverage:Job Submission Failed, Error:$@");
			exit(1);
		}
		push( @notifynames, "NotifyTilingNomapqCoverage.$id.$$.stat" );
		push( @metricsOutput, "${tiling_nomapqCoverage}.sample_interval_summary" );
	}
	return ( \@metricsOutput, \@notifynames );
}


#####################################
#####################################
#Compile Nomapq Exon Coverage Metrics
sub RunExonNomapqCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @ECmetricsFiles = ();
    my $covgheader;
    my @covgvalues;
    my @intervals; 
    my($ECout) =  $bams[0] =~ /.*_bc\d{1,4}_(.*)_L\d{1,3}.*/;
    $ECout =  $ECout . "_ALL_exonnomapqcoverage.txt"; 
    foreach my $file (@bams)
    {
	my $ecFile = $file;
	$ecFile =~ s/\.bam/\.gene_nomapq\.covg\.sample_interval_summary/;
	#print "ec:$ecFile\n";
	push(@ECmetricsFiles, "$outdir/$ecFile");
    }
    open (OUT, ">$outdir/$ECout") or die $logger->fatal("can't open output file $ECout, $!\n");
    for (my $j=0; $j<=$#ECmetricsFiles; $j++)
    {
	@intervals = ();
	open (FH, "<$ECmetricsFiles[$j]") or die $logger->info("can't open $ECmetricsFiles[$j], Skip $!\n");
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    if($hg19)
	    {
	    	$cov[0] =~ s/chr//;
	    }
	    push @intervals, $cov[0];
	    push @{$covgvalues[$j]}, $cov[2];
	}
	close (FH);
    }
    print OUT "0\tTarget";
    for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++)
    {
	my($barcode) = $ECmetricsFiles[$bc] =~ m/.*_(bc\d{1,4}).*/g;
	print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $i=0; $i<=$#{$covgvalues[0]}; $i++)
    {
	my $tally = $i+1;
	print OUT "$tally\t$intervals[$i]";
	for (my $bc=0; $bc<=$#ECmetricsFiles; $bc++) 
	{
	    print OUT "\t$covgvalues[$bc][$i]";
	}
	print OUT "\n";
    }
    close(OUT);
    return($ECout);
}


#####################################
#Compile Finger Print Tiling Coverage Metrics
sub RunFPtilingNomapqCoverageMetrics
{
    my ($bamList,$outdir) = @_;
    my @bams = @$bamList;
    my @FPTmetricsFiles = ();
    my $tilingheader;
    my @tilingvalues;
    my @tilingintervals; 
    my($FPTout) =  $bams[0] =~ /.*_bc\d{1,4}_(.*)_L\d{1,3}.*/;
    $FPTout =  $FPTout . "_ALL_tilingnomapqcoverage.txt";
    foreach my $file (@bams)
    {
	my $fptFile = $file;
	$fptFile =~ s/\.bam/\.tiling_nomapq\.covg\.sample_interval_summary/;
	#print "fpt:$fptFile\n";
	push(@FPTmetricsFiles, "$outdir/$fptFile");
    }
    open (OUT, ">$outdir/$FPTout") or die $logger->fatal("Can't open output file $FPTout. Error: $!\n");
    for (my $j=0; $j<=$#FPTmetricsFiles; $j++)
    {
	@tilingintervals = ();
	open (FH, "<$FPTmetricsFiles[$j]") or die $logger->info("Can't open $FPTmetricsFiles[$j].Skip:  $!\n");
	my $topline = <FH>;
	while (my $text = <FH>)
	{
	    chomp $text;
	    my @cov = split "\t", $text;
	    if($hg19)
	    {
	    	$cov[0] =~ s/chr//;
	    }
	    push @tilingintervals, $cov[0];
	    push @{$tilingvalues[$j]}, $cov[2];
	}
	close (FH);
    }
    print OUT "0\tTarget";
    for (my $bc=0; $bc<=$#FPTmetricsFiles; $bc++)
    {
	my($barcode) = $FPTmetricsFiles[$bc] =~ m/.*_(bc\d{1,4})_.*/g;
	print OUT "\t$barcode";
    }
    print OUT "\n";
    for (my $i=0; $i<=$#{$tilingvalues[0]}; $i++)
    {
	my $tally = $i+1;
	print OUT "$tally\t$tilingintervals[$i]";
	for (my $bc=0; $bc<=$#FPTmetricsFiles; $bc++)
	{
	    print OUT "\t$tilingvalues[$bc][$i]";
	}
	print OUT "\n";
    }
    close(OUT);
    return($FPTout);
}


sub RunLoessNormalization
{
	my $outdir = $result_dir;
	`cp $dmpCnvPath/textplot.R $outdir/`;
	my @notifynames   = ();
	
	$ENV{"R_LIBS_USER"} = $RLIBS;
	
        my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", work_dir => "$outdir", iounits => "1");
        my $additionalParams = Schedule::additionalParams(%addParams);
        my %stdParams = (scheduler => "$scheduler", job_name => "LOESS.$$", cpu => "1", mem => "8", cluster_out => "LOESS.$$.stdout", cluster_error => "LOESS.$$.stderr");
        my $standardParams = Schedule::queuing(%stdParams);
	my $cmd = "$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $RHOME/R --slave --vanilla --args $pre $outdir $GCBiasFile < $LoessNormalization";
    $logger->debug( "COMMAND: $cmd");
    
	eval{
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams "$singularityParams $RHOME/R --slave --vanilla --args $pre $outdir $GCBiasFile < $LoessNormalization"`;
		
		my %stdParams = (scheduler => "$scheduler", job_name => "NotifyLOESS.$$", job_hold => "LOESS.$$", cpu => "1", mem => "2", cluster_out => "NotifyLOESS.$$.stat", cluster_error => "/dev/null");
        	my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{cluster_error} $additionalParams $singularityParams $outdir/Notify.csh`;
    };
    if ($@) {
        $logger->fatal(
            "PostCompileWrapper:Cannot run PostCompileWrapper Error:$@"
        );
        exit(1);
    }
    push(@notifynames,"NotifyLOESS.$$.stat");
    WaitToFinish($outdir,@notifynames);
    
    return;
}



sub AddBarcodePrefix
{
    my ($bamList,$outdir) = @_;
	my $loess_file = $outdir . "/". $pre . "_ALL_intervalnomapqcoverage_loess.txt";
	my $title_file = $outdir . "/". $pre . "_title.txt";
	my $fixed_loess_file = $outdir . "/". $pre . "_fixed_ALL_intervalnomapqcoverage_loess.txt";
	my $fixed_title_file = $outdir . "/". $pre . "_fixed_title.txt";

	open(LF, "$loess_file") or die "[ERROR]: Can't open loess file $loess_file $!\n";
	open(FLF, ">$fixed_loess_file") or die "[ERROR]: Can't write to $fixed_loess_file $!\n";
	my $lf_header = <LF>;
	chomp $lf_header;
	my @lf_header_data = split(/\s+/, $lf_header);
	for(my $i = 0; $i < scalar(@lf_header_data); $i++)
	{
		$lf_header_data[$i] =~ tr/"//d;
		if($i > 0)
		{
			print FLF "\t";
		}
		if(exists $barcodeToFixedBarcode{$lf_header_data[$i]})
		{
			print FLF "\"$barcodeToFixedBarcode{$lf_header_data[$i]}\"";
		}
		else
		{
			print FLF "\"$lf_header_data[$i]\"";
		}
	}
	print FLF "\n";
	while(<LF>){
		chomp;
		print FLF $_, "\n";
	}
	close LF;
	close FLF;



	open(TF, "$title_file") or die "[ERROR]: Can't open title file $title_file $!\n";
	open(FTF, ">$fixed_title_file") or die "[ERROR]: Can't write to $fixed_title_file $!\n";
	my $tf_header = <TF>;
	chomp $tf_header;
	print FTF $tf_header, "\n";
	while(<TF>){
		chomp;
		my @data = split(/\s+/);
		if(exists( $barcodeToFixedBarcode{$data[0]} ) )
		{
			$data[0] = $barcodeToFixedBarcode{$data[0]};
		}
		print FTF join("\t", @data), "\n";
	}
	close TF;
	close FTF;
}



sub RunNormalizedCoverage
{
	my ($filenames, $outdir) = @_;
	my @names = ();
	if ($filenames) { (@names) = @$filenames; }
	my ( @notifyNames, @metricsOutNames ) = ();
	my (@sortedparseFilenames) = sort { lowestNumber($a) <=> lowestNumber($b) } @names;
	@names = @sortedparseFilenames;
	###################Calculate DepthOfCoverage
	$now = time;
	$logger->info("Started running GATK DepthOfCoverage jobs on SGE");
	for ( my $i = 0 ; $i < scalar(@names) ; $i++ )
	{
		my ( $metricsOutFiles, $waitFileNames ) = &RunDepthOfCoverage($names[$i], $outdir, $i);
		foreach my $waitName (@$waitFileNames)
		{
			push( @notifyNames, $waitName );
		}
		foreach my $outName (@$metricsOutFiles)
		{
			push( @metricsOutNames, $outName );
		}
	}

	###################waiting for metrics calculations to finish
	&WaitToFinish( $outdir, @notifyNames );
	$now = time - $now;
	$logger->info("Finished running metrics calculation jobs on SGE");
	printf( "Total Metrics Calculation run time: %02d:%02d:%02d\n\n",
			int( $now / 3600 ),
			int( ( $now % 3600 ) / 60 ),
			int( $now % 60 ) );
	##################Check metrics out files
	my $prog;
	&CheckOutputFiles( $outdir, $prog, @metricsOutNames );
	##################

	RunExonNomapqCoverageMetrics(\@names, $outdir);

	RunFPtilingNomapqCoverageMetrics(\@names, $outdir);

	RunLoessNormalization();

	if($add_barcode_prefix)
	{
		AddBarcodePrefix(\@names, $outdir);
	}
}





































