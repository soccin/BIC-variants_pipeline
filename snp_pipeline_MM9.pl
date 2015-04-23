#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

my ($pair, $bamgroup, $config, $nosnps, $target_bed, $ug);
my $pre = 'TEMP';
my $output = "results";
GetOptions ('pre=s' => \$pre,
	    'pair=s' => \$pair,
	    'config=s' => \$config,
	    'target_bed=s' => \$target_bed,
	    'nosnps' => \$nosnps,
	    'ug|unifiedgenotyper' => \$ug,
	    'bamgroup=s' => \$bamgroup,
 	    'output|out|o=s' => \$output) or exit(1);


### pre: output prefix
### group - LIST OF SAMPLE PAIRS THAT YOU WANT TO HAVE REALIGNED/RECALIBRATED TOGETHER
### -nosnps: if no snps are to be called; e.g. when only indelrealigned/recal bam needed

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

my $GATK = '';
my $PICARD = '';
my $MUTECT = '';
my $SNPEFF = '';
my $SAMTOOLS = '';
my $SOMATIC_SNIPER = '';
my $VARSCAN = '';
my $STRELKA = '';
my $SCALPEL = '';
my $VIRMID = '';

my $REF_SEQ = '/ifs/data/bio/assemblies/M.musculus/mm9/mm9.fasta';

my $curDir = `pwd`;
chomp $curDir;
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config";
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
    elsif($conf[0] =~ /snpeff/i){
	$SNPEFF = $conf[1];
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
    $target_bed = "$Bin/targets/mm9__MegaGene__v120228_MERGE.bed";
}
else{
    if(!-e $target_bed){
	die "CAN'T LOCATE $target_bed $!";
    }
}

my $multipleTargets = '';

###if($target){
###   if(!-e $target){
###	die "target file $target cannot be found $!";
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
	    next;
	}
	push @pins, "-I $pai";
	$inputFiles{$pai} = 1;
        
        my $samp = $sn[2];
        push @finalBams, "$output/alignments/$pre\_indelRealigned_recal_$samp.bam";
    }

    if(scalar(@pins) == 0){
	next;
    }

    my $bgroup = join(" ", @pins);
    $count++;

    my @indelBams = ();
    my $ran_ir == 0;
    foreach my $c (1..19, 'X', 'Y', 'M'){
	my $ran_ir_chr = 0;
	if(!-e "$output/progress/$pre\_$uID\_group_$count\_CHR$c\_RTC.done"){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_CHR$c\_RTC -pe alloc 10 -l virtual_free=1G -q lau.q,lcg.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ -L chr$c $multipleTargets -S LENIENT -nt 10 -rf BadCigar --out $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigner.intervals $bgroup`;
	    `/bin/touch $output/progress/$pre\_$uID\_group_$count\_CHR$c\_RTC.done`;
	    $ran_ir = 1;
	    $ran_ir_chr = 1;
	}
		
	if(!-e "$output/progress/$pre\_$uID\_group_$count\_CHR$c\_IR.done" || $ran_ir_chr){
	    sleep(3);
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_IR -hold_jid $pre\_$uID\_group_$count\_CHR$c\_RTC -pe alloc 1 -l virtual_free=15G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ -L chr$c $multipleTargets -S LENIENT --targetIntervals $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar --out $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned.bam $bgroup`;
	    `/bin/touch $output/progress/$pre\_$uID\_group_$count\_CHR$c\_IR.done`;
	    $ran_ir = 1;
	}
	
    	push @indelBams, "-I $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned.bam";
    }
    
    my $iBams = join(" ", @indelBams);
    my $ran_br = 0;
    if(!-e "$output/progress/$pre\_$uID\_group_$count\_BR.done" || $ran_ir){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_BR -hold_jid $pre\_$uID\_group_$count\_IR -pe alloc 6 -l virtual_free=5G $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO -R $REF_SEQ -S LENIENT --knownSites $Bin/data/UCSC_dbSNP128_MM9.bed --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 6 --out $output/intFiles/$pre\_group_$count\_recal_data.grp $iBams`;
    	`/bin/touch $output/progress/$pre\_$uID\_group_$count\_BR.done`;
	$ran_br = 1;
    }

    sleep(3);
    
    my @indelRecalBams1 = ();
    my $ran_pr = 0;
    foreach my $c (1..19, 'X', 'Y', 'M'){
	if(!-e "$output/progress/$pre\_$uID\_group_$count\_CHR$c\_PR.done" || $ran_br){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_PR -hold_jid $pre\_$uID\_group_$count\_BR -pe alloc 6 -l virtual_free=5G $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R $REF_SEQ -L chr$c $multipleTargets --emit_original_quals -BQSR $output/intFiles/$pre\_group_$count\_recal_data.grp --num_cpu_threads_per_data_thread 6 -rf BadCigar --out $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned_recal.bam -I $output/intFiles/$pre\_group_$count\_CHR$c\_indelRealigned.bam`;
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
my @iVariants = ();
my @sVariants = ();
my @hcVariants = ();
my $ran_hc = 0;
my $ran_ug_snp = 0;
my $ran_ug_indel = 0;

foreach my $c (1..19, 'X', 'Y', 'M'){
    my $irBams2 = join(" ", @{$processedBams{"chr$c"}});
    
    if($ug){
	if(!-e "$output/progress/$pre\_$uID\_CHR$c\_UG_SNP.done" || $ran_pr_glob{"chr$c"}){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_SNP -hold_jid $pre\_$uID\_PR -pe alloc 8 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name mm9 -L chr$c $multipleTargets --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model SNP --read_filter BadCigar --num_cpu_threads_per_data_thread 8 --out $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper_SNP.vcf $irBams2`;
	    `/bin/touch $output/progress/$pre\_$uID\_CHR$c\_UG_SNP.done`;
	    $ran_ug_snp = 1;
	}
	
	if(!-e "$output/progress/$pre\_$uID\_CHR$c\_UG_INDEL.done" || $ran_pr_glob{"chr$c"}){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_INDEL -hold_jid $pre\_$uID\_PR -pe alloc 8 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name mm9 -L chr$c $multipleTargets --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model INDEL --read_filter BadCigar --num_cpu_threads_per_data_thread 8 --out $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper_INDEL.vcf $irBams2`;
	    `/bin/touch $output/progress/$pre\_$uID\_CHR$c\_UG_INDEL.done`;
	    $ran_ug_indel = 1;
	}
	
	push @sVariants, "--variant $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper_SNP.vcf";
	push @iVariants, "--variant $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper_INDEL.vcf";
    }
    
    if(!-e "$output/progress/$pre\_$uID\_CHR$c\_HC.done" || $ran_pr_glob{"chr$c"}){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HC -hold_jid $pre\_$uID\_PR -pe alloc 12 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L chr$c $multipleTargets --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_CHR$c\_HaplotypeCaller.vcf $irBams2`;
	`/bin/touch $output/progress/$pre\_$uID\_CHR$c\_HC.done`;
	$ran_hc = 1;
    }
    
    push @hcVariants, "--variant $output/intFiles/$pre\_CHR$c\_HaplotypeCaller.vcf";
}

my $hcVars = join(" " , @hcVariants);
my $ran_cv_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_CV_HC.done" || $ran_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_HC -hold_jid $pre\_$uID\_HC -pe alloc 1 -l virtual_free=2G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf --assumeIdenticalSamples $hcVars`;
    `/bin/touch $output/progress/$pre\_$uID\_CV_HC.done`;
    $ran_cv_hc = 1;
}

my $ran_snpeff_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_SNPEFF_HC.done" || $ran_cv_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SNPEFF_HC -hold_jid  $pre\_$uID\_CV_HC -pe alloc 1 -l virtual_free=8G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx8g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -v -i vcf -o gatk -no-downstream -no-intergenic -no-intron -no-upstream -no-utr GRCm37.67 $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf ">$output/intFiles/$pre\_HaplotypeCaller_RAW_snpEff_output.vcf"`;
    `/bin/touch $output/progress/$pre\_$uID\_SNPEFF_HC.done`;
    $ran_snpeff_hc = 1;
}

my $ran_va_hc = 0;
if(!-e "$output/progress/$pre\_$uID\_VA_HC.done" || $ran_snpeff_hc){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VA_HC -hold_jid $pre\_$uID\_SNPEFF_HC -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q  $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantAnnotator -R $REF_SEQ -A SnpEff --variant $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf --snpEffFile $output/intFiles/$pre\_HaplotypeCaller_RAW_snpEff_output.vcf -o $output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf`;
    `/bin/touch $output/progress/$pre\_$uID\_VA_HC.done`;
    $ran_va_hc = 1;   
}

if(!-e "$output/progress/$pre\_$uID\_MAF_HC.done" || $ran_va_hc){
    sleep(3);
    &generateMaf("$output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'haplotypecaller', "$pre\_$uID\_VA_HC");
    `/bin/touch $output/progress/$pre\_$uID\_MAF_HC.done`;
}

if($ug){
    `/bin/mkdir -m 775 -p $output/variants/unifiedgenotyper`;
    my $sVars = join(" " , @sVariants);
    my $iVars = join(" " , @iVariants);

    my $ran_cv_ug_snp = 0;
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_SNP.done" || $ran_ug_snp){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_UG_SNP -hold_jid $pre\_$uID\_UG_SNP -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf --assumeIdenticalSamples $sVars`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_SNP.done`;
	$ran_cv_ug_snp = 1;
    }

    my $ran_cv_ug_indel = 0;
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_INDEL.done" || $ran_ug_indel){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_UG_INDEL -hold_jid $pre\_$uID\_UG_INDEL -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf --assumeIdenticalSamples $iVars`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_INDEL.done`;
	$ran_cv_ug_indel = 1;
    }

    if(!-e "$output/progress/$pre\_$uID\_CV_UG_RAW.done" || $ran_cv_ug_snp || $ran_cv_ug_indel){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_UG_RAW -hold_jid $pre\_$uID\_CV_UG_SNP,$pre\_$uID\_CV_UG_INDEL -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper_RAW.vcf --assumeIdenticalSamples --variant $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf --variant $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_RAW.done`;
    }

     my $ran_vf_ug_snp = 0;
    if(!-e "$output/progress/$pre\_$uID\_VF_UG_SNP.done" || $ran_cv_ug_snp){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_VF -hold_jid $pre\_$uID\_CV_UG_SNP,$pre\_$uID\_CV_UG_INDEL -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf --maskName nearIndel --variant $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf -o $output/intFiles/$pre\_UnifiedGenotyper_SNP_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"MQ \\< 40.0\\" --filterExpression \\"FS \\> 60.0\\" --filterExpression \\"HaplotypeScore \\> 13.0\\" --filterExpression \\"MQRankSum \\< -12.5\\" --filterExpression \\"ReadPosRankSum \\< -8.0\\" --filterName QDFilter --filterName MQFilter --filterName FSFilter --filterName HSFilter --filterName MQRSFilter --filterName ReadPosFilter`;
	`/bin/touch $output/progress/$pre\_$uID\_VF_UG_SNP.done`;
	$ran_vf_ug_snp = 1;
    }

     my $ran_vf_ug_indel = 0;
    if(!-e "$output/progress/$pre\_$uID\_VF_UG_INDEL.done" || $ran_cv_ug_indel){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_VF -hold_jid $pre\_$uID\_CV_UG_INDEL -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --variant $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf -o $output/intFiles/$pre\_UnifiedGenotyper_INDEL_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"ReadPosRankSum \\< -20.0\\" --filterExpression \\"InbreedingCoeff \\< -0.8\\" --filterExpression \\"FS \\> 200.0\\" --filterName QDFilter --filterName ReadPosFilter --filterName InbreedingFilter --filterName FSFilter`;
	`/bin/touch $output/progress/$pre\_$uID\_VF_UG_INDEL.done`;
	$ran_vf_ug_indel = 1;
    }

    my $ran_cv_ug_si = 0;
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_SI.done" || $ran_cv_ug_snp || $ran_cv_ug_indel){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_UG_SI -hold_jid $pre\_$uID\_UG_VF -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf --assumeIdenticalSamples --variant $output/progress/$pre\_UnifiedGenotyper_vf.vcf --variant $output/progress/$pre\_UnifiedGenotyper_INDEL_vf.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_SI.done`;
	$ran_cv_ug_si = 1;
    }

    my $ran_snpeff_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_SNPEFF_UG.done" || $ran_cv_ug_si){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SNPEFF_UG -hold_jid $pre\_$uID\_CV_UG_SI -pe alloc 1 -l virtual_free=8G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx8g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -v -i vcf -o gatk -no-downstream -no-intergenic -no-intron -no-upstream -no-utr GRCm37.67 $output/progress/$pre\_UnifiedGenotyper_vf.vcf ">$output/progress/$pre\_UnifiedGenotyper_vf_snpEff_output.vcf"`;
	`/bin/touch $output/progress/$pre\_$uID\_SNPEFF_UG.done`;
	$ran_snpeff_ug = 1;
    }

    my $ran_va_ug = 0;
    if(!-e "$output/progress/$pre\_$uID\_VA_UG.done" || $ran_snpeff_ug){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_VA_UG -hold_jid $pre\_$uID\_SNPEFF_UG -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantAnnotator -R $REF_SEQ -A SnpEff --variant  $output/progress/$pre\_UnifiedGenotyper_vf.vcf--snpEffFile $output/progress/$pre\_UnifiedGenotyper_vf_snpEff_output.vcf -o $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_VA_UG.done`;
	$ran_va_ug = 1;
    }
        
    if(!-e "$output/progress/$pre\_$uID\_MAF_UG.done" || $ran_va_ug){  
	sleep(3);
	&generateMaf("$pre\_UnifiedGenotyper.vcf", 'unifiedgenotyper', "$pre\_$uID\_VA_UG");
 	`/bin/touch $output/progress/$pre\_$uID\_MAF_UG.done`;
    }
}

if($pair){
    open(PAIR, "$pair") or die "Can't open $pair file";
    my %seen = ();
    
    ### NOTE: THE SAMPLE NAMES IN THE PAIRING FILE MUST MATCH EXACTLY THE SAMPLE NAMES IN THE REALIGNED/RECALIBRATED BAM FILE
    `/bin/mkdir -m 775 -p $output/variants/mutect`;
    `/bin/mkdir -m 775 -p $output/variants/somaticsniper`;
    `/bin/mkdir -m 775 -p $output/variants/virmid`;
    `/bin/mkdir -m 775 -p $output/variants/scalpel`;
    `/bin/mkdir -m 775 -p $output/variants/strelka`;
    ###`/bin/mkdir -m 775 -p $output/variants/varscan`;
    while(<PAIR>){
	chomp;
	
	my @data = split(/\s+/, $_);

	if($data[0] =~ /^NA$/i || $data[1] =~ /^NA$/i){
	    next;
	}

	my $ran_mutect = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done" || $ran_ssf){  
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MUTECT -hold_jid $pre\_$uID\_SSF -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xmx4g -Djava.io.tmpdir=/scratch/$uID -jar $MUTECT/muTect.jar --analysis_type MuTect --reference_sequence $REF_SEQ --input_file:normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --input_file:tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --vcf $output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf --out $output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.txt -rf BadCigar --enable_extended_output --downsampling_type NONE`;
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
	###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $pre\_$uID\_SSFI -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -f $REF_SEQ $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam ">$output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup"`;

	###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $pre\_$uID\_SSF -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -f $REF_SEQ $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam ">$output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup"`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MPILEUP.done`;
	    ###$ran_mpileup = 1;
	###}	    

	my $ran_virmid = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done" || $ran_ssf){  
	    if(-d "$output/variants/virmid/$data[0]\_$data[1]\_virmid"){
		`/bin/rm -rf $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    }
	    `/bin/mkdir -m 775 -p $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VIRMID -hold_jid $pre\_$uID\_SSF -pe alloc 4 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.6.0_24/bin/java -Xmx12g -Djava.io.tmpdir=/scratch/$uID -jar $VIRMID/Virmid.jar -R $REF_SEQ -D $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam -N $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam -t 4 -o $pre\_$data[0]\_$data[1]\_virmid -w $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done`;
	    $ran_virmid = 1;
	}
	
	my $ran_strelka_config = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done" || $ran_ssf){  
	    if(-d "$output/variants/strelka/$data[0]\_$data[1]\_strelka"){
		### STRELKA DIES IF DIR ALREADY EXISTS
		`/bin/rm -rf $output/variants/strelka//$data[0]\_$data[1]\_strelka`;
	    }
	    
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_LNS -hold_jid $pre\_$uID\_SSF -pe alloc 1 -l virtual_free=1G $Bin/qCMD /bin/ln -s $pre\_indelRealigned_recal\_$data[0]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam.bai`;
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_LNS -hold_jid $pre\_$uID\_SSF -pe alloc 1 -l virtual_free=1G $Bin/qCMD /bin/ln -s $pre\_indelRealigned_recal\_$data[1]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam.bai`;
	    
	    ### NOTE: STRELKA ONLY HAS CONFIG FOR BWA ALN, NOT SURE HOW IT WILL WORK WITH BWA MEM
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG -hold_jid $pre\_$uID\_LNS -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $STRELKA/bin/configureStrelkaWorkflow.pl --normal=$output/alignments//$pre\_indelRealigned_recal\_$data[0]\.bam --tumor=$output/alignments//$pre\_indelRealigned_recal\_$data[1]\.bam --ref=$REF_SEQ --config=$STRELKA/etc/strelka_config_bwa_default.ini --output-dir=$output/variants/strelka//$data[0]\_$data[1]\_strelka`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done`;
	    $ran_strelka_config = 1;
	}

	my $ran_scalpel = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done" || $ran_ssf){
	    `/bin/mkdir -m 775 -p $output/variants/scalpel/$data[0]\_$data[1]\_scalpel`;
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_SCALPEL -hold_jid $pre\_$uID\_SSF -pe alloc 4 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SCALPEL/scalpel --somatic --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --bed $target_bed --ref $REF_SEQ --dir $output/variants/scalpel/$data[0]\_$data[1]\_scalpel --numprocs 4`;
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
	###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/java/jdk1.6.0_24/bin/java -Xmx2g -Djava.io.tmpdir=/scratch/$uID -jar $VARSCAN/VarScan.jar somatic $output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup $output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic --strand-filter 1 --output-vcf 1`;
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
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.snvs.vcf -c strelka -o $output/variants/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
	    
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.snvs.vcf -c strelka -o $output/variants/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
	    
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.indels.vcf -c strelka -o $output/variants/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
	    
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.indels.vcf -c strelka -o $output/variants/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
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
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$type\_MAF_PAIRED -hold_jid $hold -pe alloc 4 -l virtual_free=2G,internet=1 -q ito.q $Bin/qCMD $Bin/generateMAF_MM9.pl -vcf $vcf -pairing $pair -species mm9 -config $config -caller $type -target $target_bed $n_sample $t_sample`;

	if($type =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
	    `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$type\_MAF_UNPAIRED -hold_jid $hold -pe alloc 4 -l virtual_free=2G,internet=1 -q ito.q $Bin/qCMD $Bin/generateMAF_MM9.pl -vcf $vcf -species mm9 -config $config -caller $type -target $target_bed`;
	}
    }
    else{
	`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$type\_MAF_UNPAIRED -hold_jid $hold -pe alloc 4 -l virtual_free=2G,internet=1 -q ito.q $Bin/qCMD $Bin/generateMAF_MM9.pl -vcf $vcf -species mm9 -config $config -caller $type -target $target_bed`;
    }
}
