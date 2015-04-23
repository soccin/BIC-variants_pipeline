#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

my $pre = 'TEMP';
my ($group, $config, $nosnps, $ug);
GetOptions ('pre=s' => \$pre,
	    'config=s' => \$config,
	    'nosnps' => \$nosnps,
	    'ug|unifiedgenotyper' => \$ug,
	    'group=s' => \$group);


### pre: output prefix
### group - LIST OF SAMPLE PAIRS THAT YOU WANT TO HAVE REALIGNED/RECALIBRATED TOGETHER
### -nosnps: if no snps are to be called; e.g. when only indelrealigned/recal bam needed

### NOTE: CAN'T RECALIBRATE BECAUSE OF LACK OF KNOWN SITES 

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

my $GATK = '';
my $PICARD = '';
my $SNPEFF = '';

my $REF_SEQ = '/ifs/data/bio/assemblies/D.melanogaster/dm3/dm3.fasta';

open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /gatk/i){
	$GATK = $conf[1];
    }
    elsif($conf[0] =~ /picard/i){
	$PICARD = $conf[1];
    }
    elsif($conf[0] =~ /snpeff/i){
	$SNPEFF = $conf[1];
    }
}
close CONFIG;

my @indelBams = ();
my @indelBams2 = ();  
my $count = 0;
my %inputFiles = ();

open(IN, "$group") || die "CAN'T OPEN GROUPING FILE OF MARKDUP BAMS $group $!";
while(<IN>){
    chomp;
    
    my @pair = split(/\s+/, $_);
    my @pins = ();
    foreach my $pai (@pair){
	if($inputFiles{$pai}){
	    next;
	}
	push @pins, "-I $pai";
	$inputFiles{$pai} = 1;
    }

    if(scalar(@pins) == 0){
	next;
    }

    my $group = join(" ", @pins);
    $count++;
    
    my $cancerMode = '';
    if(scalar(@pins) > 1){
	$cancerMode = "--cancer_mode";
    }

    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_RTC -pe alloc 10 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ -S LENIENT -nt 10 -rf BadCigar --out $pre\_group_$count\_indelRealigner.intervals $group`;
    
    sleep(5);
    
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_group_$count\_IR -hold_jid $pre\_$uID\_group_$count\_RTC -pe alloc 2 -l virtual_free=26G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx48g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ -S LENIENT --targetIntervals $pre\_group_$count\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 5000 -rf BadCigar --out $pre\_group_$count\_indelRealigned.bam $group`;
    
    push @indelBams, "I=$pre\_group_$count\_indelRealigned.bam";
    push @indelBams2, "-I $pre\_group_$count\_indelRealigned.bam";
}

sleep(5);

my $iBams = join(" ", @indelBams);
my $iBams2 = join(" ", @indelBams2);

`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE -hold_jid $pre\_$uID\_IR -pe alloc 24 -l virtual_free=3G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $iBams O=$pre\_indelRealigned.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=true`;

if($nosnps){
    exit;
}
 
sleep(5);
       
`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HC -hold_jid $pre\_$uID\_IR -pe alloc 12 -l virtual_free=3G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $pre\_HaplotypeCaller.vcf $iBams2`;
    
sleep(5);

`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HC_SNPEFF -hold_jid  $pre\_$uID\_HC -pe alloc 1 -l virtual_free=8G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx8g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -v -i vcf -o gatk -no-downstream -no-intergenic -no-intron -no-upstream -no-utr dm5.48 $pre\_HaplotypeCaller.vcf ">$pre\_HaplotypeCaller_snpEff_output.vcf"`;

sleep(5);

`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HC_VA -hold_jid $pre\_$uID\_HC_SNPEFF -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantAnnotator -R $REF_SEQ -A SnpEff --variant $pre\_HaplotypeCaller.vcf --snpEffFile $pre\_HaplotypeCaller_snpEff_output.vcf -o $pre\_HaplotypeCaller_va_snpEff.vcf`;

if($ug){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_SNP -hold_jid $pre\_$uID\_IR -pe alloc 8 -l virtual_free=2G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name dm3 --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model SNP --read_filter BadCigar --num_cpu_threads_per_data_thread 8 --out $pre\_UnifiedGenotyper_SNP.vcf $iBams2`;

    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_INDEL -hold_jid $pre\_$uID\_IR -pe alloc 8 -l virtual_free=2G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name dm3 --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model INDEL --read_filter BadCigar --num_cpu_threads_per_data_thread 8 --out $pre\_UnifiedGenotyper_INDEL.vcf $iBams2`;

    sleep(5);

    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_VF -hold_jid $pre\_$uID\_UG_SNP,$pre\_$uID\_UG_INDEL -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $pre\_UnifiedGenotyper_INDEL.vcf --maskName nearIndel --variant $pre\_UnifiedGenotyper_SNP.vcf -o $pre\_UnifiedGenotyper_SNP_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"MQ \\< 40.0\\" --filterExpression \\"FS \\> 60.0\\" --filterExpression \\"HaplotypeScore \\> 13.0\\" --filterExpression \\"MQRankSum \\< -12.5\\" --filterExpression \\"ReadPosRankSum \\< -8.0\\" --filterName QDFilter --filterName MQFilter --filterName FSFilter --filterName HSFilter --filterName MQRSFilter --filterName ReadPosFilter`;

    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_VF -hold_jid $pre\_$uID\_UG_INDEL -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --variant $pre\_UnifiedGenotyper_INDEL.vcf -o $pre\_UnifiedGenotyper_INDEL_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"ReadPosRankSum \\< -20.0\\" --filterExpression \\"InbreedingCoeff \\< -0.8\\" --filterExpression \\"FS \\> 200.0\\" --filterName QDFilter --filterName ReadPosFilter --filterName InbreedingFilter --filterName FSFilter`;

    sleep(5);

    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CV_UG_SI -hold_jid $pre\_$uID\_UG_VF -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $pre\_UnifiedGenotyper.vcf --assumeIdenticalSamples --variant $pre\_UnifiedGenotyper_SNP_vf.vcf --variant $pre\_UnifiedGenotyper_INDEL_vf.vcf`;

sleep(5);

    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_SNPEFF -hold_jid $pre\_$uID\_CV_UG_SI -pe alloc 1 -l virtual_free=8G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Xms256m -Xmx8g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -v -i vcf -o gatk -no-downstream -no-intergenic -no-intron -no-upstream -no-utr dm5.48 $pre\_UnifiedGenotyper.vcf ">$pre\_UnifiedGenotyper_snpEff_output.vcf"`;

    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_UG_VA -hold_jid $pre\_$uID\_UG_SNPEFF -pe alloc 1 -l virtual_free=1G -q lau.q $Bin/qCMD /opt/java/jdk1.7.0_25/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantAnnotator -R $REF_SEQ -A SnpEff --variant $pre\_UnifiedGenotyper.vcf --snpEffFile $pre\_UnifiedGenotyper_snpEff_output.vcf -o $pre\_UnifiedGenotyper_va_snpEff.vcf`;
}
