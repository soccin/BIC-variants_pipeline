#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

my ($pair, $bamgroup, $config, $nosnps, $target_bed, $ug, $hybrid);
my $pre = 'TEMP';
my $output = "results";
GetOptions ('pre=s' => \$pre,
	    'pair=s' => \$pair,
	    'config=s' => \$config,
	    'target_bed=s' => \$target_bed,
	    'nosnps' => \$nosnps,
	    'ug|unifiedgenotyper' => \$ug,
	    'hybrid' => \$hybrid,
	    'bamgroup=s' => \$bamgroup,
 	    'output|out|o=s' => \$output) or exit(1);


### pre: output prefix

### pairing: sample pairing information for mutect/maf conversion
###          NOTE: if no pairing file provided, will assume no somatic analysis

### config: file listing paths to all used programs; otherwise expects them in /opt/bin

### nosnp: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed

### bamgroup - LIST OF SAMPLE PAIRS THAT YOU WANT TO HAVE REALIGNED/RECALIBRATED TOGETHER
###          CAN BE NORMAL/TUMOR PAIRING, TRIPLETS, OR ANY OF PAIRING
###          CAN BE A GROUPING OF ANY LENGTH >= 1
###         ***  WARNING ***
###              MAKE SURE THAT THE SAME SAMPLE DOESN'T SHOW UP IN MULTIPLE LINES I.E. MULTIPLE COMPARISONS
###              BECAUSE THAT SAMPLE WILL SHOW UP MULTIPLE TIMES IN FINAL BAM
###              COLLAPSE ALL SUCH SAMPLES INTO A SINGLE COMPARISON
###          ***         ***
###          HOWEVER, GATK SAYS OTHER THAN NORMAL/TUMOR PAIRING, EVERYTHING ELSE SHOULD BE DONE AT SAMPLE LEVEL, NOT COHORT LEVEL
###          ONLY AFFECTS REALIGNMENT/RECALIBRATION PORTION
###          ALL SAMPLES IN FILE WILL THEN BE USED TO CALL SNPS/INDELS TOGETHER 
###          NOTE: MAKE SURE THAT ALL INPUTS HAVE BEEN MARKDUPed AT THE SAMPLE LEVEL (USE mergeBamsSameLib.pl)
###          IN THIS VERSION, REDUCED READS BAMS ISN'T FULLY OPTIMIZED BECAUSE IT REDUCES READ AT THE 
###          COMPARISON LEVEL AND NOT FULL COHORT LEVEL

### NOTE: TARGETS MUST HAVE BED AND ILIST FILES
### NOTE: BAIT MUST BE AN ILIST FILE

### NOTE: LINUX HAS A ARG_MAX OF 131072 CHARS IN THE COMMAND LINE

### GROUP: sample grouping information

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

my $GATK = '';
my $PICARD = '';
my $MUTECT = '';
my $SAMTOOLS = '';
my $SOMATIC_SNIPER = '';
my $VARSCAN = '';
my $STRELKA = '';
my $SCALPEL = '';
my $VIRMID = '';

my $REF_SEQ = '/ifs/data/bio/assemblies/H.sapiens/hg19/hg19.fasta';
if($hybrid){
    $REF_SEQ = '/ifs/data/bio/assemblies/hybrid_H.sapiens_M.musculus/hybrid_hg19_mm10/hybrid_hg19_mm10.fasta';
}

my $curDir = `pwd`;
chomp $curDir;
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;
    
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /gatk/i){
	$GATK = $conf[1];
    }
    elsif($conf[0] =~ /mutect/i){
	$MUTECT = $conf[1];
    }
    elsif($conf[0] =~ /picard/i){
	$PICARD = $conf[1];
    }
    elsif($conf[0] =~ /samtools/i){
	$SAMTOOLS = $conf[1];
    }
    elsif($conf[0] =~ /somaticsniper/i){
	$SOMATIC_SNIPER = $conf[1];
    }
    elsif($conf[0] =~ /varscan/i){
	$VARSCAN = $conf[1];
    }
    elsif($conf[0] =~ /strelka/i){
	$STRELKA = $conf[1];
    }
    elsif($conf[0] =~ /scalpel/i){
	$SCALPEL = $conf[1];
    }
    elsif($conf[0] =~ /virmid/i){
	$VIRMID = $conf[1];
    }
}
close CONFIG;

### make sure all markdup bam files are there before proceeding
open(BGR, "$bamgroup") || die "CAN'T OPEN GROUPING FILE OF MARKDUP BAMS $bamgroup $!";
while(<BGR>){
    chomp;

    my @bgr = split(/\s+/, $_);
    foreach my $bg (@bgr){
	if(!-e $bg){
	    die "file $bg does not exist";
	}
    }
}
close BGR;

if(!$target_bed){
    $target_bed = "$Bin/targets/hg19__MegaGene__v131104.bed";
    if($hybrid){
	$target_bed = "$Bin/targets/hg19_mm10_hybrid_MegaGene__v131104.bed";
    }
}
else{
    if(!-e $target_bed){
	die "CAN'T LOCATE $target_bed $!";
    }
}

my $multipleTargets = '';
### NOTE: THERE SEEMS TO BE SOMETHING WRONG WITH DIPS IN COVERAGE AT TARGET BOUNDARY
###       SO NOT USING THE TARGETS FOR ANALYSI
###if($target_bed){
###   if(!-e $target_bed){
###	die "target file $target_bed cannot be found $!";
###    }
###    $multipleTargets = "-L $target_bed --interval_set_rule INTERSECTION";
###}

my $count = 0;
my %inputFiles = ();
my %processedBams = ();
my @finalBams = ();
my %ran_pr_glob = 0;
my $ran_ssf = 0;

`/bin/mkdir -m 775 -p $output`;
`/bin/mkdir -m 775 -p $output/intFiles`;
`/bin/mkdir -m 775 -p $output/alignments`;
`/bin/mkdir -m 775 -p $output/progress`;

open(IN, "$bamgroup") || die "CAN'T OPEN GROUPING FILE OF MARKDUP BAMS $bamgroup $!";
while(<IN>){
    chomp;
    
    my @pair = split(/\s+/, $_);
    my @pins = ();
    foreach my $pai (@pair){
	my @sn = split(/\//, $pai);
	if($inputFiles{$pai}){
	    ### IN CASE A FILE SHOWS UP MULTIPLE TIMES DUE TO BEING IN MULTIPLE COMPARISONS AND WASN'T COLLAPSED e.g. MET ANALYSIS
	    ### THIS MAKES SURE THAT A FILE ISN'T INCLUDED MULTIPLE TIMES IN PROCESSING
	    next;
	}
	push @pins, "-I $pai";
	$inputFiles{$pai} = 1;

        ## store final bams for post-recalibration stats
        my $samp = $sn[2];
        push @finalBams, "$output/alignments/$pre\_indelRealigned_recal_$samp.bam"; 
    }

    if(scalar(@pins) == 0){
	### IN CASE ALL ENTRIES IN LINE CONTAINED BAMS ALREADY SEEN IN OTHER COMPARISONS
	next;
    }

    my $bgroup = join(" ", @pins);
    $count++;
    
    my @indelBams = ();
    my $ran_ir == 0;
    foreach my $c (1..22, 'X', 'Y', 'M'){
	my $ran_ir_chr = 0;
	if(!-e "$output/progress/$pre\_$uID\_group_$count\_CHR$c\_RTC.done"){
	    ### mad.q,nce.q have timeout issues with this step
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_CHR$c\_RTC -pe alloc 10 -l virtual_free=1G -q lau.q,lcg.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ -L chr$c $multipleTargets --known $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf --known $Bin/data/dbsnp_135.hg19__ReTag.vcf -S LENIENT -nt 10 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigner.intervals $bgroup`;
	    `/bin/touch $output/progress/$pre\_$uID\_group_$count\_CHR$c\_RTC.done`;
	    $ran_ir = 1;
	    $ran_ir_chr = 1;
	}
	
	### seems to randomly have issues with file locking timeout on the m nodes so not submitting to it anymore
	if(!-e "$output/progress/$pre\_$uID\_group_$count\_CHR$c\_IR.done" || $ran_ir_chr){
	    sleep(3);
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_IR -hold_jid $pre\_$uID\_group_$count\_CHR$c\_RTC -pe alloc 1 -l virtual_free=15G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ -L chr$c $multipleTargets --knownAlleles $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf --knownAlleles $Bin/data/dbsnp_135.hg19__ReTag.vcf -S LENIENT --targetIntervals $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar --out $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned.bam $bgroup`;
	    `/bin/touch $output/progress/$pre\_$uID\_group_$count\_CHR$c\_IR.done`;
	    $ran_ir = 1;
	}
	
    	push @indelBams, "-I $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned.bam";
    }

    my $iBams = join(" ", @indelBams);
    my $ran_br = 0;
    if(!-e "$output/progress/$pre\_$uID\_group_$count\_BR.done" || $ran_ir){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_BR -hold_jid $pre\_$uID\_group_$count\_IR -pe alloc 6 -l virtual_free=5G $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO -R $REF_SEQ -S LENIENT --knownSites $Bin/data/dbsnp_135.hg19__ReTag.vcf --knownSites $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf --knownSites $Bin/data/hapmap_3.3.hg19.vcf --knownSites $Bin/data/1000G_omni2.5.hg19.vcf --knownSites $Bin/data/1000G_phase1.snps.high_confidence.hg19.vcf --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 6 --out $output/intFiles/$pre\_group_$count\_recal_data.grp $iBams`;
	`/bin/touch $output/progress/$pre\_$uID\_group_$count\_BR.done`;
	$ran_br = 1;
    }
    
    my @indelRecalBams1 = ();
    my $ran_pr = 0;
    foreach my $c (1..22, 'X', 'Y', 'M'){
	if(!-e "$output/progress/$pre\_$uID\_group_$count\_CHR$c\_PR.done" || $ran_br){
	    sleep(3);
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_PR -hold_jid $pre\_$uID\_group_$count\_BR -pe alloc 6 -l virtual_free=5G $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R $REF_SEQ -L chr$c $multipleTargets --emit_original_quals -BQSR $output/intFiles/$pre\_group_$count\_recal_data.grp --num_cpu_threads_per_data_thread 6 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned_recal.bam -I $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned.bam`;
	    `/bin/touch $output/progress/$pre\_$uID\_group_$count\_CHR$c\_PR.done`;
	    $ran_pr = 1;
	    $ran_pr_glob{"chr$c"} = 1;
	}
    
	push @indelRecalBams1, "I=$output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned_recal.bam";
	push @{$processedBams{"chr$c"}}, "-I $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned_recal.bam";
    }

    my $irBams1 = join(" ", @indelRecalBams1);
    my $ran_m = 0;
    if(!-e "$output/progress/$pre\_$uID\_group_$count\_MERGE.done" || $ran_pr){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE -hold_jid $pre\_$uID\_PR -pe alloc 8 -l virtual_free=4G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $irBams1 O=$output/intFiles/$pre\_group_$count\_indelRealigned_recal.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
	`/bin/touch $output/progress/$pre\_$uID\_group_$count\_MERGE.done`;
	$ran_m = 1;
    }

    if(!-e "$output/progress/$pre\_$uID\_group_$count\_SSF.done" || $ran_m){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SSF -hold_jid $pre\_$uID\_MERGE -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T SplitSamFile -R $REF_SEQ -I $output/intFiles/$pre\_group_$count\_indelRealigned_recal.bam --outputRoot $output/alignments/$pre\_indelRealigned_recal_`;
	`/bin/touch $output/progress/$pre\_$uID\_group_$count\_SSF.done`;
	$ran_ssf = 1;
    }

}

if(!-e "$output/progress/$pre\_$uID\_MQ.done" || $ran_ssf){
    my $ran_mq = 0;
    foreach my $finalBam (@finalBams){
        my @sn = split(/\//, $finalBam);
        my $samp = $sn[-1];
        $samp =~ s/\.bam//g;
        $samp =~ s/$pre\_indelRealigned_recal_//g;
        `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MQ_METRICS -hold_jid $pre\_$uID\_SSF -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MeanQualityByCycle INPUT=$finalBam OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.txt CHART_OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.pdf REFERENCE_SEQUENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true TMP_DIR=/scratch/$uID`;
	`/bin/touch $output/progress/$pre\_$uID\_MQ.done`;
	$ran_mq = 1;
    }

    
    if(!-e "$output/progress/$pre\_$uID\_MERGE_MQ.done" || $ran_mq){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_MQ -hold_jid '$pre\_$uID\_MQ_METRICS' -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergeMeanQualityHistograms.py . '*_MeanQualityByCycle_*.txt' $output/metrics/$pre\_post_recal_MeanQualityByCycle.txt $output/metrics/$pre\_pre_recal_MeanQualityByCycle.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_MQ.done`;
    }
}


if($nosnps){
    exit;
}

`/bin/mkdir -m 775 -p $output/variants`;
`/bin/mkdir -m 775 -p $output/variants/haplotypecaller`;
my @ugVariants = ();
my @hcVariants = ();
my $ran_ug = 0;
my $ran_hc = 0;
foreach my $c (1..22, 'X', 'Y', 'M'){
    ### NOTE1: GET THIS RANDOM ERROR NOW FOR SNP CALLS
    ###       Unexpectedly couldn't find valid codec for temporary output file
    ###       RERUNNING THE SAME COMMAND SEEMS TO FIX IT
    ###
    ### NOTE2: GET SYSTEM MEMORY ISSUE WHEN RUNNING ON THE CLUSTER AND LOTS
    ###        OF INPUT FILES TO UNIFIEDGENOTYPE WITH NUM THREADS > 1 -nt
    ###        NOT A PROBLEM WHEN RUNNING IT WITH LOTS OF INPUT FILES ON RAY
    
    my $irBams2 = join(" ", @{$processedBams{"chr$c"}});

    if($ug){
	if(!-e "$output/progress/$pre\_$uID\_CHR$c\_UG.done" || $ran_pr_glob{"chr$c"}){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG -hold_jid $pre\_$uID\_PR -pe alloc 8 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name hg19 -L chr$c $multipleTargets --dbsnp $Bin/data/dbsnp_135.hg19__ReTag.vcf --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model BOTH --read_filter BadCigar --num_cpu_threads_per_data_thread 8 --out $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper.vcf $irBams2`;
	`/bin/touch $output/progress/$pre\_$uID\_CHR$c\_UG.done`;
	$ran_ug = 1;
	}
    }

    ### NOTE: ANNOTATIONS THAT DON'T WORK:AlleleBalance, HardyWeinberg,IndelType
    if(!-e "$output/progress/$pre\_$uID\_CHR$c\_HC.done" || $ran_pr_glob{"chr$c"}){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HC -hold_jid $pre\_$uID\_PR -pe alloc 24 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx75g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L chr$c $multipleTargets --dbsnp $Bin/data/dbsnp_135.hg19__ReTag.vcf --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 24 --out $output/intFiles/$pre\_CHR$c\_HaplotypeCaller.vcf $irBams2`;
	`/bin/touch $output/progress/$pre\_$uID\_CHR$c\_HC.done`;
	$ran_hc = 1;
    }
    
    push @ugVariants, "--variant $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper.vcf";
    push @hcVariants, "--variant $output/intFiles/$pre\_CHR$c\_HaplotypeCaller.vcf";
    
}

my $hcVars = join(" " , @hcVariants);
my $ran_cv_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_CV_HC.done" || $ran_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_HC -hold_jid $pre\_$uID\_HC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf --assumeIdenticalSamples $hcVars`;
    `/bin/touch $output/progress/$pre\_$uID\_CV_HC.done`;
    $ran_cv_hc = 1;
}

### NOTE: InbreedingCoeff requires >= 10 samples
###       doing a stupid count of just number of unique MD bams that are inputs
###       this makes the assumption that there is only 1 MD bam per sample
###       which isn't always true e.g. when a sample has multiple libraries
###       we treat it as 2 "samples" but gatk will not

my $ran_vr_snp_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_VR_SNP_HC.done" || $ran_cv_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VR_SNP_HC -hold_jid $pre\_$uID\_CV_HC -pe alloc 4 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $Bin/data/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 $Bin/data/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 $Bin/data/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $Bin/data/dbsnp_135.hg19__ReTag.vcf -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_SNP.plots.R -nt 8`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_SNP_HC.done`;
    $ran_vr_snp_hc = 1;
}

### NOTE: sometimes throws an error when running with multiple threads
###       about not being able to find a tmp file; running with -nt 1 to avoid errors
my $ran_ar_snp_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_AR_SNP_HC.done" || $ran_vr_snp_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_AR_SNP_HC -hold_jid $pre\_$uID\_VR_SNP_HC -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -o $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_SNP_HC.done`;
    $ran_ar_snp_hc = 1;
}

my $ran_vr_indel_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_VR_INDEL_HC.done" || $ran_ar_snp_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VR_INDEL_HC -hold_jid $pre\_$uID\_AR_SNP_HC -pe alloc 4 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.plots.R -nt 8`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_INDEL_HC.done`;
    $ran_vr_indel_hc = 1;   
}

my $ran_ar_indel_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_AR_INDEL_HC.done" || $ran_vr_indel_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_AR_INDEL_HC -hold_jid $pre\_$uID\_VR_INDEL_HC -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -o $output/intFiles/$pre\_HaplotypeCaller_vqsr.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_INDEL_HC.done`;
    $ran_ar_indel_hc = 1;   
}

### Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP) (URL: http://evs.gs.washington.edu/EVS/)
my $ran_vf_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_VF_HC.done" || $ran_ar_indel_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VF_HC -hold_jid $pre\_$uID\_AR_INDEL_HC -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $Bin/data/ESP6500SI-V2-SSA137.updatedRsIds.snps_indels.vcf --maskName GO_ESP --variant $output/intFiles/$pre\_HaplotypeCaller_vqsr.vcf -o $output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf`;
    `/bin/touch $output/progress/$pre\_$uID\_VF_HC.done`;
    $ran_vf_hc = 1;   
}

if(!-e "$output/progress/$pre\_$uID\_HC_MAF.done" || $ran_vf_hc){
    sleep(3);
    &generateMaf("$output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'haplotypecaller', "$pre\_$uID\_VF_HC");
    `/bin/touch $output/progress/$pre\_$uID\_HC_MAF.done`;
}

if($ug){
    `/bin/mkdir -m 775 -p $output/variants/unifiedgenotyper`;
    my $ugVars = join(" " , @ugVariants);
    my $ran_cv_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_CV_UG.done" || $ran_ug){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_UG -hold_jid $pre\_$uID\_UG -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper_RAW.vcf --assumeIdenticalSamples $ugVars`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG.done`;
	$ran_cv_ug = 1;
    }

    my $ran_vr_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_VR_UG.done" || $ran_cv_ug){  
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VR_SNP_UG -hold_jid $pre\_$uID\_CV_UG -pe alloc 4 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $Bin/data/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 $Bin/data/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 $Bin/data/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $Bin/data/dbsnp_135.hg19__ReTag.vcf -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.plots.R -nt 8`;
	`/bin/touch $output/progress/$pre\_$uID\_VR_UG.done`;
	$ran_vr_ug = 1;
    }

    my $ran_ar_snp_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_AR_SNP_UG.done" || $ran_vr_ug){  
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_AR_SNP_UG -hold_jid $pre\_$uID\_VR_SNP_UG -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -o $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -nt 1`;
	`/bin/touch $output/progress/$pre\_$uID\_AR_SNP_UG.done`;
	$ran_ar_snp_ug = 1;
    }

    my $ran_vr_indel_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_VR_INDEL_UG.done" || $ran_ar_snp_ug){  
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VR_INDEL_UG -hold_jid $pre\_$uID\_AR_SNP_UG -pe alloc 4 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.plots.R -nt 8`;
	`/bin/touch $output/progress/$pre\_$uID\_VR_INDEL_UG.done`;
	$ran_vr_indel_ug = 1;
    }

    my $ran_ar_indel_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_AR_INDEL_UG.done" || $ran_vr_indel_ug){  
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_AR_INDEL_UG -hold_jid $pre\_$uID\_VR_INDEL_UG -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -o $output/intFiles/$pre\_UnifiedGenotyper_vqsr.vcf -nt 1`;
	`/bin/touch $output/progress/$pre\_$uID\_AR_INDEL_UG.done`;
	$ran_ar_indel_ug = 1;
    }

    my $ran_vf_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_VF_UG.done" || $ran_ar_indel_ug){  
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VF_UG -hold_jid $pre\_$uID\_AR_INDEL_UG -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $Bin/data/ESP6500SI-V2-SSA137.updatedRsIds.snps_indels.vcf --maskName GO_ESP --variant $output/intFiles/$pre\_UnifiedGenotyper_vqsr.vcf -o $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_VF_UG.done`;
	$ran_vf_ug = 1;
    }

    if(!-e "$output/progress/$pre\_$uID\_UG_MAF.done" || $ran_vf_ug){  
	sleep(3);
	&generateMaf("$output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf", 'unifiedgenotyper', "$pre\_$uID\_VF_UG");
	`/bin/touch $output/progress/$pre\_$uID\_UG_MAF.done`;
    }
}

if($pair){
    ### QC pairing
    if(!-e "$output/progress/$pre\_$uID\_FP_WRAP.done" || $ran_vf_hc){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_FP_WRAP -hold_jid $pre\_$uID\_VF_HC $Bin/qCMD $Bin/qc/pairing_qc.pl -pair $pair -vcf $output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf -pre $pre -bin $Bin -outdir $output/metrics/fingerprint`;
	`/bin/touch $output/progress/$pre\_$uID\_FP_WRAP.done`;
    }
    
    `/bin/mkdir -m 775 -p $output/variants/mutect`;
    `/bin/mkdir -m 775 -p $output/variants/somaticsniper`;
    `/bin/mkdir -m 775 -p $output/variants/virmid`;
    `/bin/mkdir -m 775 -p $output/variants/scalpel`;
    `/bin/mkdir -m 775 -p $output/variants/strelka`;
    ###`/bin/mkdir -m 775 -p $output/variants/varscan`;

    open(PAIR, "$pair") or die "Can't open $pair file";
    ### NOTE: THE SAMPLE NAMES IN THE PAIRING FILE MUST MATCH EXACTLY THE SAMPLE NAMES IN THE REALIGNED/RECALIBRATED BAM FILE
    while(<PAIR>){
	chomp;
	
	my @data = split(/\s+/, $_);

	### pairing file also contains unpaired samples with NA for the other pai
	### so we can keep track of what sample is tumor or normal
	if($data[0] =~ /^NA$/i || $data[1] =~ /^NA$/i){
	    next;
	}
	
	### NOTE: NOT AUTOMATICALLY RUNNING EACH SOMATIC CALLER WHEN $pre\_$uID\_SSF RUNS
	###       BCAUSE DON'T WANT TO KICK IT OFF FOR ALL SAMPLE PAIRS
	###       IN CASE ONLY ONE SAMPLE GROUP GOT MODIFIED ABOVE

	###       WILL RUN SOMATIC ANALYSIS FOR ALL SAMPLE PAIRS
	###       IF JUST ONE SAMPLE HAD ITS BAM MODIFIED
	my $ran_mutect = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done" || $ran_ssf){  
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MUTECT -hold_jid $pre\_$uID\_SSF -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xmx4g -Djava.io.tmpdir=/scratch/$uID -jar $MUTECT/muTect.jar --analysis_type MuTect --reference_sequence $REF_SEQ --dbsnp $Bin/data/dbsnp_135.hg19__ReTag.vcf --cosmic $Bin/data/CosmicCodingMuts_v67_20131024.vcf --input_file:normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --input_file:tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --vcf $output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf --out $output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.txt -rf BadCigar --enable_extended_output --downsampling_type NONE`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done`;
	    $ran_mutect = 1;
	}

	my $ran_somatic_sniper = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.done" || $ran_ssf){  
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER -hold_jid $pre\_$uID\_SSF -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SOMATIC_SNIPER/bam-somaticsniper -F vcf -f $REF_SEQ -q 1 $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam $output/variants/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.done`;
	    $ran_somatic_sniper = 1;
	}	    

	###my $ran_mpileup = 0;
	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MPILEUP.done" || $ran_ssf){  
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $pre\_$uID\_SSF -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -f $REF_SEQ -d 500000 $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam ">$output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup"`;

	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $pre\_$uID\_SSF -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -f $REF_SEQ -d 500000 $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam ">$output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup"`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MPILEUP.done`;
	    ###$ran_mpileup = 1;
	###}	    

	### BECAUSE THE OUTPUT PREFIX FOR VIRMID IS THE DISEASE BAM,
	### NEED TO CREATE A DIRECTORY FOR EACH TUMOR/NORMAL PAIR
	### SINCE A DISEASE BAM CAN BE USED IN MULTIPLE COMPARISONS
	### virmid fails if directory already exists
	my $ran_virmid = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done" || $ran_ssf){  
	    if(-d "$output/variants/virmid/$data[0]\_$data[1]\_virmid"){
		`/bin/rm -rf $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    }
	    `/bin/mkdir -m 775 -p $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VIRMID -hold_jid $pre\_$uID\_SSF -pe alloc 4 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.6.0_24/bin/java -Xms256m -Xmx12g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $VIRMID/Virmid.jar -R $REF_SEQ -D $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam -N $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam -t 4 -o $pre\_$data[0]\_$data[1]\_virmid -w $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done`;
	    $ran_virmid = 1;
	}

	my $ran_strelka_config = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done" || $ran_ssf){  
	    if(-d "$output/variants/strelka/$data[0]\_$data[1]\_strelka"){
		### strelka DIES IF DIR ALREADY EXISTS
		`/bin/rm -rf $output/variants/strelka/$data[0]\_$data[1]\_strelka`;
	    }

	    ### NOTE: Strelka only recognizes X.bam.bai as the index for X.bam, not X.bai
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_LNS -hold_jid $pre\_$uID\_SSF -pe alloc 1 -l virtual_free=1G $Bin/qCMD /bin/ln -s $pre\_indelRealigned_recal\_$data[0]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam.bai`;
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_LNS -hold_jid $pre\_$uID\_SSF -pe alloc 1 -l virtual_free=1G $Bin/qCMD /bin/ln -s $pre\_indelRealigned_recal\_$data[1]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam.bai`;
	    
	    ### NOTE: strelka ONLY HAS CONFIG FOR BWA ALN, NOT SURE HOW IT WILL WORK WITH BWA MEM
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG -hold_jid $pre\_$uID\_LNS -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $STRELKA/bin/configureStrelkaWorkflow.pl --normal=$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor=$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --ref=$REF_SEQ --config=$STRELKA/etc/strelka_config_bwa_default.ini --output-dir=$output/variants/strelka/$data[0]\_$data[1]\_strelka`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done`;
	    $ran_strelka_config = 1;
	}

	my $ran_scalpel = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done" || $ran_ssf){
	    `/bin/mkdir -m 775 -p $output/variants/scalpel/$data[0]\_$data[1]\_scalpel`;
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_SCALPEL -hold_jid $pre\_$uID\_SSF -pe alloc 24 -l virtual_free=3G $Bin/qCMD $SCALPEL/scalpel --somatic --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --bed $target_bed --ref $REF_SEQ --dir $output/variants/scalpel/$data[0]\_$data[1]\_scalpel --numprocs 24`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done`;
	    $ran_scalpel = 1;
	}

	my $ran_strelka_run = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done" || $ran_strelka_config){
	    sleep(3);
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG -pe alloc 8 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /usr/bin/make -C $output/variants/strelka/$data[0]\_$data[1]\_strelka -j 8`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done`;
	    $ran_strelka_run = 1;
	}

	###my $ran_varscan_somatic = 0;
	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC.done" || $ran_mpileup){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -pe alloc 2 -l virtual_free=5G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.6.0_24/bin/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $VARSCAN/VarScan.jar somatic $output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup $output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic --strand-filter 1 --output-vcf 1`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC.done`;
	    ###$ran_varscan_somatic = 1;
	###}

	sleep(3);

	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done" || $ran_mutect){
	    &generateMaf("$output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf", 'mutect', "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT", $data[0], $data[1]);
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done`;
	}

	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.done" || $ran_somatic_sniper){
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf -c somaticsniper -o $output/variants/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper_MAF.txt -n $data[0] -t $data[1]`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.done`;
	}

	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done" || $ran_varscan_somatic){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp_MAF.txt -n $data[0] -t $data[1]`;

	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.indel.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.indel_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done`;
	###}

	if($ran_scalpel){
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_SCALPEL_CLEANUP -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_SCALPEL -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q,mad.q $Bin/qCMD /bin/rm -rf $output/variants/scalpel/$data[0]\_$data[1]\_scalpel/main/ $output/variants/scalpel/$data[0]\_$data[1]\_scalpel/validation/`;
	}

	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.done" || $ran_strelka_run){
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.snvs.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
       
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.snvs.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;

	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.indels.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;

	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.indels.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.done`;
	    
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q,mad.q $Bin/qCMD /bin/rm -rf $output/variants/strelka/$data[0]\_$data[1]\_strelka/config $output/variants/strelka/$data[0]\_$data[1]\_strelka/chromosomes $output/variants/strelka/$data[0]\_$data[1]\_strelka/Makefile $output/variants/strelka/$data[0]\_$data[1]\_strelka/task.complete`;
	}
    }
    close PAIR;
}

sub generateMaf{
    my ($vcf, $type, $hold, $normal_sample, $tumor_sample) = @_;

    my $n_sample = '';
    my $t_sample = '';
    if($normal_sample && $tumor_sample){
	$n_sample = "-normal_sample $normal_sample";
	$t_sample = "-tumor_sample $tumor_sample";
    }

    if($pair){
	### NOTE: ASKING FOR PE ALLOC 4 TO THROTTLE NUMBER OF JOBS ACCESING ONCOTATOR
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$type\_MAF_PAIRED -hold_jid $hold -pe alloc 4 -l virtual_free=2G,internet=1 -q ito.q $Bin/qCMD $Bin/generateMAF.pl -vcf $vcf -pairing $pair -species hg19 -config $config -caller $type $n_sample $t_sample -delete_temp`;

	if($type =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$type\_MAF_UNPAIRED -hold_jid $hold -pe alloc 4 -l virtual_free=2G,internet=1 -q ito.q $Bin/qCMD $Bin/generateMAF.pl -vcf $vcf -species hg19 -config $config -caller $type -delete_temp`;
	}
    }
    else{
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$type\_MAF_UNPAIRED -hold_jid $hold -pe alloc 4 -l virtual_free=2G,internet=1 -q ito.q $Bin/qCMD $Bin/generateMAF.pl -vcf $vcf -species hg19 -config $config -caller $type -delete_temp`;
    }
}
