#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Schedule;
use Cluster;

my ($pair, $group, $bamgroup, $config, $nosnps, $targets, $ug, $hybrid, $scheduler, $priority_project, $priority_group, $abra, $help, $step1);

my $pre = 'TEMP';
my $output = "results";

GetOptions ('pre=s' => \$pre,
	    'pair=s' => \$pair,
	    'group=s' => \$group,
	    'config=s' => \$config,
	    'targets=s' => \$targets,
	    'nosnps' => \$nosnps,
	    'ug|unifiedgenotyper' => \$ug,
	    'hybrid' => \$hybrid,
	    'abra' => \$abra,
	    'step1' => \$step1,
	    'bamgroup=s' => \$bamgroup,
 	    'scheduler=s' => \$scheduler,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
	    'help' => \$help,
 	    'output|out|o=s' => \$output) or exit(1);


if(!$group || !$config || !$scheduler || !$targets || !$bamgroup || $help){
    print <<HELP;

    USAGE: process_alignments_hg19.pl -group GROUP -pair PAIR -pre PRE -config CONFIG -species SPECIES -scheduler SCHEDULER -targets TARGETS
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED)
	* BAMGROUP: files listing bams to be processed together; every bam for each group on 1 line, comma-separated (required)
	* TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
	* PRE: output prefix (default: TEMP)
	* OUTPUT: output results directory (default: results)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* -nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed
	* -abra: run abra instead of GATK indelrealigner
	* -step1: forece the pipeline to start from the first step in pipeline
	* -hybrid: data aligned to hybrid assembly
	* haplotypecaller is default; -ug || -unifiedgenotyper to also make unifiedgenotyper variant calls	
HELP
exit;
}

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

my $ABRA = '';
my $GATK = '';
my $PICARD = '';
my $MUTECT = '';
my $SAMTOOLS = '';
my $SOMATIC_SNIPER = '';
my $VARSCAN = '';
my $STRELKA = '';
my $SCALPEL = '';
my $VIRMID = '';
my $JAVA = '';
my $PYTHON = '';
my $PERL = '';
my $HG19_FASTA = '';
my $HG19_MM10_HYBRID_FASTA = '';
my $HG19_BWA_INDEX = '';
my $HG19_MM10_HYBRID_BWA_INDEX = '';

my $curDir = `pwd`;
chomp $curDir;
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;
    
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /abra/i){
	if(!-e "$conf[1]/abra.jar"){
	    die "CAN'T FIND GenomeAnalysisTK.jar IN $conf[1] $!";
	}
	$ABRA = $conf[1];
    }
    elsif($conf[0] =~ /gatk/i){
	if(!-e "$conf[1]/GenomeAnalysisTK.jar"){
	    die "CAN'T FIND GenomeAnalysisTK.jar IN $conf[1] $!";
	}
	$GATK = $conf[1];
    }
    elsif($conf[0] =~ /mutect/i){
	if(!-e "$conf[1]/muTect.jar"){
	    die "CAN'T FIND muTect.jar IN $conf[1] $!";
	}
	$MUTECT = $conf[1];
    }
    elsif($conf[0] =~ /picard/i){
	if(!-e "$conf[1]/picard.jar"){
	    die "CAN'T FIND picard.jar IN $conf[1] $!";
	}
	$PICARD = $conf[1];
    }
    elsif($conf[0] =~ /samtools/i){
	if(!-e "$conf[1]/samtools"){
	    die "CAN'T FIND samtools IN $conf[1] $!";
	}
	$SAMTOOLS = $conf[1];
    }
    elsif($conf[0] =~ /somaticsniper/i){
	if(!-e "$conf[1]/bam-somaticsniper"){
	    die "CAN'T FIND bam-somaticsniper IN $conf[1] $!";
	}
	$SOMATIC_SNIPER = $conf[1];
    }
    elsif($conf[0] =~ /varscan/i){
	if(!-e "$conf[1]/VarScan.jar"){
	    die "CAN'T FIND VarScan.jar IN $conf[1] $!";
	}
	$VARSCAN = $conf[1];
    }
    elsif($conf[0] =~ /strelka/i){
	if(!-e "$conf[1]/bin/configureStrelkaWorkflow.pl"){
	    die "CAN'T FIND bin/configureStrelkaWorkflow.pl IN $conf[1] $!";
	}
	$STRELKA = $conf[1];
    }
    elsif($conf[0] =~ /scalpel/i){
	if(!-e "$conf[1]/scalpel"){
	    die "CAN'T FIND scalpel IN $conf[1] $!";
	}
	$SCALPEL = $conf[1];
    }
    elsif($conf[0] =~ /virmid/i){
	if(!-e "$conf[1]/Virmid.jar"){
	    die "CAN'T FIND Virmid.jar IN $conf[1] $!";
	}
	$VIRMID = $conf[1];
    }
    elsif($conf[0] =~ /java/i){
	if(!-e "$conf[1]/java"){
	    die "CAN'T FIND java IN $conf[1] $!";
	}
	$JAVA = $conf[1];
    }
    elsif($conf[0] =~ /perl/i){
	if(!-e "$conf[1]/perl"){
	    die "CAN'T FIND perl IN $conf[1] $!";
	}
	$PERL = $conf[1];
    }
    elsif($conf[0] =~ /python/i){
	if(!-e "$conf[1]/python"){
	    die "CAN'T FIND python IN $conf[1] $!";
	}
	$PYTHON = $conf[1];
    }
    elsif($conf[0] =~ /hg19_fasta/i){
	if(!-e "$conf[1]"){
	    die "CAN'T FIND $conf[1] $!";
	}
	$HG19_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm10_hybrid_fasta/i){
	if(!-e "$conf[1]"){
	    die "CAN'T FIND $conf[1] $!";
	}
	$HG19_MM10_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /hg19_bwa_index/i){
	if(!-e "$conf[1]"){
	    die "CAN'T FIND hg19 bwa index with prefix $conf[1] $!";
	}
	$HG19_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm10_hybrid_bwa_index/i){
	if(!-e "$conf[1]"){
	    ###die "CAN'T FIND hg19-mm10 hybrid bwa index with prefix $conf[1] $!";
	}
	$HG19_MM10_HYBRID_BWA_INDEX = $conf[1];
    }
}
close CONFIG;

my $REF_SEQ = "$HG19_FASTA";
my $BWA_INDEX = "$HG19_BWA_INDEX";
if($hybrid){
    $REF_SEQ = "$HG19_MM10_HYBRID_FASTA";
    $BWA_INDEX = "$HG19_MM10_HYBRID_BWA_INDEX";
}

### make sure all markdup bam files are there before proceeding
open(BGR, "$bamgroup") || die "CAN'T OPEN GROUPING FILE OF MARKDUP BAMS $bamgroup $!";
while(<BGR>){
    chomp;

    my @bgro = split(/\s+/, $_);
    my @bgr = split(/,/, $bgro[1]);
    foreach my $bg (@bgr){
	if(!-e $bg){
	    die "file $bg does not exist";
	}
    }
}
close BGR;

my $targets_bed_padded = "$Bin/targets/$targets/$targets\_targets_plus5bp.bed";
if(!-e "$targets_bed_padded"){
    die "CAN'T LOCATE $targets_bed_padded FOR $targets; REQUIRED FOR SCALPEL $!";
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
my @prg_jids = ();
my $ran_ssf = 0;
my @ssf_jids = ();

`/bin/mkdir -m 775 -p $output`;
`/bin/mkdir -m 775 -p $output/intFiles`;
`/bin/mkdir -m 775 -p $output/alignments`;
`/bin/mkdir -m 775 -p $output/progress`;

my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);

open(IN, "$bamgroup") || die "CAN'T OPEN GROUPING FILE OF MARKDUP BAMS $bamgroup $!";
while(<IN>){
    chomp;
    

    my @gpair = split(/\s+/, $_);
    my @pair = split(/,/, $gpair[1]);
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
    my @indelBams = ();
    my $ran_ir = 0;
    my @ir_jids = ();

    if($abra){
	my @inBams = ();
	my @outBams = ();
	foreach my $pin (@pins){
	    my @inB = split(/\s+/, $pin);
	    push @inBams, $inB[1];
	    push @outBams, "$inB[1]\_ABRA.bam";
	    push @indelBams, "-I $inB[1]\_ABRA.bam\_FM.bam";
	}
	
	my $aiBams = join(",", @inBams);
	my $aoBams = join(",", @outBams);
	my $ran_abra = 0;
	my $abra_jid = '';
	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_ABRA.done" || $step1){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_ABRA", cpu => "12", mem => "500", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_ABRA.log");
	    my $standardParams = Schedule::queuing(%stdParams);	    
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx100g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $ABRA/abra.jar --in $aiBams --out $aoBams --ref $REF_SEQ --bwa-ref $BWA_INDEX --targets $Bin/targets/abra_hg19.bed --working $output/intFiles/abra_$gpair[0] --threads 12`;
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/abra_wrapper.pl -inBams $aiBams -outBams $aoBams -refSeq $REF_SEQ -bwaRef $BWA_INDEX -targets $Bin/targets/abra_hg19.bed -working $output/intFiles/abra_$gpair[0] -config $config -log $output/progress/$pre\_$uID\_$gpair[0]\_ABRA_WRAPPER.log`;

	    $abra_jid = "$pre\_$uID\_$gpair[0]\_ABRA";
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_ABRA.done`;
	    $ran_abra = 1;
	}

	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_FIXMATE.done" || $ran_abra){
	    my $bcount = 0;
	    foreach my $outBam (@outBams){
		$bcount++;
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$bcount\_FIXMATE", job_hold => "$abra_jid", cpu => "1", mem => "50", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_FIXMATE.log");
		my $standardParams = Schedule::queuing(%stdParams);
		my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "5");
		my $additionalParams = Schedule::additionalParams(%addParams);

		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar FixMateInformation I=$outBam O=$outBam\_FM.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true`;
		push @ir_jids, "$pre\_$uID\_$gpair[0]\_$bcount\_FIXMATE";
	    }
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_FIXMATE.done`;
	    $ran_ir = 1;
	}
    }
    else{
	foreach my $c (1..22, 'X', 'Y', 'M'){
	    my $ran_rtc = 0;
	    my $rtc_jid = '';
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_RTC.done" || $step1){
		### mad.q,nce.q have timeout issues with this step
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_CHR$c\_RTC", cpu => "10", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_RTC.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ -L chr$c $multipleTargets --known $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf --known $Bin/data/dbsnp_135.hg19__ReTag.vcf -S LENIENT -nt 10 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigner.intervals $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_RTC.done`;
		$rtc_jid = "$pre\_$uID\_$gpair[0]\_CHR$c\_RTC";
		$ran_rtc = 1;
	    }
	    
	    ### seems to randomly have issues with file locking timeout on the m nodes so not submitting to it anymore
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_IR.done" || $ran_rtc){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_CHR$c\_IR", job_hold => "$rtc_jid", cpu => "1", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_IR.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ -L chr$c $multipleTargets --knownAlleles $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf --knownAlleles $Bin/data/dbsnp_135.hg19__ReTag.vcf -S LENIENT --targetIntervals $output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigned.bam $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_IR.done`;
		push @ir_jids, "$pre\_$uID\_$gpair[0]\_CHR$c\_IR";
		$ran_ir = 1;
	    }

	    push @indelBams, "-I $output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigned.bam";
	}
    }

    my $irBams = join(" ", @indelBams);
    my $ran_br = 0;
    my $irj = join(",", @ir_jids);
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_BR.done" || $ran_ir){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_BR", job_hold => "$irj", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_BR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO -R $REF_SEQ -S LENIENT --knownSites $Bin/data/dbsnp_135.hg19__ReTag.vcf --knownSites $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf --knownSites $Bin/data/hapmap_3.3.hg19.vcf --knownSites $Bin/data/1000G_omni2.5.hg19.vcf --knownSites $Bin/data/1000G_phase1.snps.high_confidence.hg19.vcf --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 6 --out $output/intFiles/$pre\_$gpair[0]\_recal_data.grp $irBams`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_BR.done`;
	$ran_br = 1;
    }
    
    my @indelRecalBams1 = ();
    my $ran_pr = 0;
    my @pr_jids = ();
    foreach my $c (1..22, 'X', 'Y', 'M'){
	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_PR.done" || $ran_br){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_CHR$c\_PR", job_hold => "$pre\_$uID\_$gpair[0]\_BR", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_PR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R $REF_SEQ -L chr$c $multipleTargets --emit_original_quals -BQSR $output/intFiles/$pre\_$gpair[0]\_recal_data.grp --num_cpu_threads_per_data_thread 6 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigned_recal.bam -I $output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigned.bam`;
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_CHR$c\_PR.done`;
	    push @pr_jids, "$pre\_$uID\_$gpair[0]\_CHR$c\_PR";
	    push @prg_jids, "$pre\_$uID\_$gpair[0]\_CHR$c\_PR";
	    $ran_pr = 1;
	    $ran_pr_glob{"chr$c"} = 1;
	}
    
	push @indelRecalBams1, "I=$output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigned_recal.bam";
	push @{$processedBams{"chr$c"}}, "-I $output/intFiles/$pre\_$gpair[0]\_CHR$c\_indelRealigned_recal.bam";
    }

    my $irBams1 = join(" ", @indelRecalBams1);
    my $ran_m = 0;
    my $prj = join(",", @pr_jids);
    my @merge_jids = ();
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE.done" || $ran_pr){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_MERGE", job_hold => "$prj", cpu => "8", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE.log");
	my $standardParams = Schedule::queuing(%stdParams);
	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "10");
	my $additionalParams = Schedule::additionalParams(%addParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $irBams1 O=$output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_MERGE.done`;
	push @merge_jids, "$pre\_$uID\_$gpair[0]\_MERGE";
	$ran_m = 1;
    }

    my $mj = join(",", @merge_jids);
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.done" || $ran_m){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_SSF", job_hold => "$mj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T SplitSamFile -R $REF_SEQ -I $output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam --outputRoot $output/alignments/$pre\_indelRealigned_recal_`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_SSF.done`;
	push @ssf_jids, "$pre\_$uID\_$gpair[0]\_SSF";
	$ran_ssf = 1;
    }

}

my $ssfj = join(",", @ssf_jids);
if(!-e "$output/progress/$pre\_$uID\_MQ.done" || $ran_ssf){
    my $ran_mq = 0;
    my @mq_metrics_jid = ();
    foreach my $finalBam (@finalBams){
        my @sn = split(/\//, $finalBam);
        my $samp = $sn[-1];
        $samp =~ s/\.bam//g;
        $samp =~ s/$pre\_indelRealigned_recal_//g;

	if(!-e "$output/progress/$pre\_$uID\_MQ_METRICS_$samp\.done" || $ran_ssf){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MQ_METRICS_$samp", job_hold => "$ssfj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_MQ_METRICS_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MeanQualityByCycle INPUT=$finalBam OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.txt CHART_OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.pdf REFERENCE_SEQUENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true TMP_DIR=/scratch/$uID`;
	    push @mq_metrics_jid, "$pre\_$uID\_MQ_METRICS_$samp";
	    `/bin/touch $output/progress/$pre\_$uID\_MQ_METRICS_$samp\.done`;
	    $ran_mq = 1;
	}
    }

    my $mqmj = join(",", @mq_metrics_jid);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_MQ.done" || $ran_mq){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_MQ", job_hold => "$mqmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_MQ.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeMeanQualityHistograms.py . '*_MeanQualityByCycle_*.txt' $output/metrics/$pre\_post_recal_MeanQualityByCycle.txt $output/metrics/$pre\_pre_recal_MeanQualityByCycle.txt`;
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
my @ug_jids = ();
my @hc_jids = ();
my $prgj = join(",", @prg_jids);
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
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CHR$c\_UG", job_hold => "$prgj", cpu => "8", mem => "24", cluster_out => "$output/progress/$pre\_$uID\_CHR$c\_UG.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name hg19 -L chr$c $multipleTargets --dbsnp $Bin/data/dbsnp_135.hg19__ReTag.vcf --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model BOTH --read_filter BadCigar --num_cpu_threads_per_data_thread 8 --out $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper.vcf $irBams2`;
	    `/bin/touch $output/progress/$pre\_$uID\_CHR$c\_UG.done`;
	    push @ug_jids, "$pre\_$uID\_CHR$c\_UG";
	    $ran_ug = 1;
	}
    }

    ### NOTE: ANNOTATIONS THAT DON'T WORK:AlleleBalance, HardyWeinberg,IndelType
    if(!-e "$output/progress/$pre\_$uID\_CHR$c\_HC.done" || $ran_pr_glob{"chr$c"}){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CHR$c\_HC", job_hold => "$prgj", cpu => "12", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_CHR$c\_HC.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L chr$c $multipleTargets --dbsnp $Bin/data/dbsnp_135.hg19__ReTag.vcf --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_CHR$c\_HaplotypeCaller.vcf $irBams2`;
	`/bin/touch $output/progress/$pre\_$uID\_CHR$c\_HC.done`;
	push @hc_jids, "$pre\_$uID\_CHR$c\_HC";
	$ran_hc = 1;
    }
    
    push @ugVariants, "--variant $output/intFiles/$pre\_CHR$c\_UnifiedGenotyper.vcf";
    push @hcVariants, "--variant $output/intFiles/$pre\_CHR$c\_HaplotypeCaller.vcf";
    
}

my $hcVars = join(" " , @hcVariants);
my $ran_cv_hc = 0;
my $hcj = join(",", @hc_jids);
my $cvhcj = '';
if(!-e "$output/progress/$pre\_$uID\_CV_HC.done" || $ran_hc){
    sleep(3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_HC", job_hold => "$hcj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf --assumeIdenticalSamples $hcVars`;
    `/bin/touch $output/progress/$pre\_$uID\_CV_HC.done`;
    $cvhcj = "$pre\_$uID\_CV_HC";
    $ran_cv_hc = 1;
}

### NOTE: InbreedingCoeff requires >= 10 samples
###       doing a stupid count of just number of unique MD bams that are inputs
###       this makes the assumption that there is only 1 MD bam per sample
###       which isn't always true e.g. when a sample has multiple libraries
###       we treat it as 2 "samples" but gatk will not

my $ran_vr_snp_hc = 0;
my $vrshcj = '';
if(!-e "$output/progress/$pre\_$uID\_VR_SNP_HC.done" || $ran_cv_hc){
    sleep(3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_SNP_HC", job_hold => "$cvhcj", cpu => "4", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_VR_SNP_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $Bin/data/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 $Bin/data/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 $Bin/data/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $Bin/data/dbsnp_135.hg19__ReTag.vcf -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_SNP.plots.R -nt 4`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_SNP_HC.done`;
    $vrshcj = "$pre\_$uID\_VR_SNP_HC";
    $ran_vr_snp_hc = 1;
}

### NOTE: sometimes throws an error when running with multiple threads
###       about not being able to find a tmp file; running with -nt 1 to avoid errors
my $ran_ar_snp_hc = 0;
my $arshcj = '';
if(!-e "$output/progress/$pre\_$uID\_AR_SNP_HC.done" || $ran_vr_snp_hc){
    sleep(3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_SNP_HC", job_hold => "$vrshcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_SNP_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/variants/haplotypecaller/$pre\_HaplotypeCaller_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -o $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_SNP_HC.done`;
    $arshcj = "$pre\_$uID\_AR_SNP_HC";
    $ran_ar_snp_hc = 1;
}

my $ran_vr_indel_hc = 0;
my $vrihcj = '';
if(!-e "$output/progress/$pre\_$uID\_VR_INDEL_HC.done" || $ran_ar_snp_hc){
    sleep(3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_INDEL_HC", job_hold => "$arshcj", cpu => "4", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_VR_INDEL_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.plots.R -nt 4`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_INDEL_HC.done`;
    $vrihcj = "$pre\_$uID\_VR_INDEL_HC";
    $ran_vr_indel_hc = 1;   
}

my $ran_ar_indel_hc = 0;
my $arihcj = '';
if(!-e "$output/progress/$pre\_$uID\_AR_INDEL_HC.done" || $ran_vr_indel_hc){
    sleep(3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_INDEL_HC", job_hold => "$vrihcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_INDEL_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -o $output/intFiles/$pre\_HaplotypeCaller_vqsr.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_INDEL_HC.done`;
    $arihcj = "$pre\_$uID\_AR_INDEL_HC";
    $ran_ar_indel_hc = 1;   
}

### Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP) (URL: http://evs.gs.washington.edu/EVS/)
my $ran_vf_hc = 0;
my $vfhcj = '';
if(!-e "$output/progress/$pre\_$uID\_VF_HC.done" || $ran_ar_indel_hc){
    sleep(3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VF_HC", job_hold => "$arihcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_VF_HC.log");
   my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $Bin/data/ESP6500SI-V2-SSA137.updatedRsIds.snps_indels.vcf --maskName GO_ESP --variant $output/intFiles/$pre\_HaplotypeCaller_vqsr.vcf -o $output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf`;
    `/bin/touch $output/progress/$pre\_$uID\_VF_HC.done`;
    $vfhcj = "$pre\_$uID\_VF_HC";
    $ran_vf_hc = 1;   
}

if(!-e "$output/progress/$pre\_$uID\_HC_MAF.done" || $ran_vf_hc){
    sleep(3);
    &generateMaf("$output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'haplotypecaller', "$vfhcj");
    `/bin/touch $output/progress/$pre\_$uID\_HC_MAF.done`;
}

if($ug){
    `/bin/mkdir -m 775 -p $output/variants/unifiedgenotyper`;
    my $ugVars = join(" " , @ugVariants);
    my $ran_cv_ug = 0;
    my $cvugj = '';
    my $ugj = join(",", @ug_jids);
    if(!-e "$output/progress/$pre\_$uID\_CV_UG.done" || $ran_ug){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_UG", job_hold => "$ugj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper_RAW.vcf --assumeIdenticalSamples $ugVars`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG.done`;
	$cvugj = "$pre\_$uID\_CV_UG";
	$ran_cv_ug = 1;
    }

    my $ran_vr_snp_ug = 0;
    my $vrsugj = '';
    if(!-e "$output/progress/$pre\_$uID\_VR_SNP_UG.done" || $ran_cv_ug){  
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_SNP_UG", job_hold => "$cvugj", cpu => "4", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_VR_SNP_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $Bin/data/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 $Bin/data/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 $Bin/data/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $Bin/data/dbsnp_135.hg19__ReTag.vcf -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.plots.R -nt 4`;
	`/bin/touch $output/progress/$pre\_$uID\_VR_SNP_UG.done`;
	$vrsugj = "$pre\_$uID\_VR_SNP_UG";
	$ran_vr_snp_ug = 1;
    }

    my $ran_ar_snp_ug = 0;
    my $arsugj = '';
    if(!-e "$output/progress/$pre\_$uID\_AR_SNP_UG.done" || $ran_vr_snp_ug){  
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_SNP_UG", job_hold => "$vrsugj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_SNP_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -o $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -nt 1`;
	`/bin/touch $output/progress/$pre\_$uID\_AR_SNP_UG.done`;
	$arsugj = "$pre\_$uID\_AR_SNP_UG";
	$ran_ar_snp_ug = 1;
    }

    my $ran_vr_indel_ug = 0;
    my $vriugj = '';
    if(!-e "$output/progress/$pre\_$uID\_VR_INDEL_UG.done" || $ran_ar_snp_ug){  
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_INDEL_UG", job_hold => "$arsugj", cpu => "4", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_VR_INDEL_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $Bin/data/Mills_and_1000G_gold_standard.indels.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.plots.R -nt 4`;
	`/bin/touch $output/progress/$pre\_$uID\_VR_INDEL_UG.done`;
	$vriugj = "$pre\_$uID\_VR_INDEL_UG";
	$ran_vr_indel_ug = 1;
    }

    my $ran_ar_indel_ug = 0;
    my $ariugj = '';
    if(!-e "$output/progress/$pre\_$uID\_AR_INDEL_UG.done" || $ran_vr_indel_ug){  
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_INDEL_UG", job_hold => "$vriugj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_INDEL_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -o $output/intFiles/$pre\_UnifiedGenotyper_vqsr.vcf -nt 1`;
	`/bin/touch $output/progress/$pre\_$uID\_AR_INDEL_UG.done`;
	$ariugj = "$pre\_$uID\_AR_INDEL_UG";
	$ran_ar_indel_ug = 1;
    }

    my $ran_vf_ug = 0;
    my $vfugj = '';
    if(!-e "$output/progress/$pre\_$uID\_VF_UG.done" || $ran_ar_indel_ug){  
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VF_UG", job_hold => "$ariugj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_VF_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $Bin/data/ESP6500SI-V2-SSA137.updatedRsIds.snps_indels.vcf --maskName GO_ESP --variant $output/intFiles/$pre\_UnifiedGenotyper_vqsr.vcf -o $output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_VF_UG.done`;
	$vfugj = "$pre\_$uID\_VF_UG";
	$ran_vf_ug = 1;
    }

    ###if(!-e "$output/progress/$pre\_$uID\_UG_MAF.done" || $ran_vf_ug){  
	###sleep(3);
	###&generateMaf("$output/variants/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf", 'unifiedgenotyper', "$vfugj");
	###`/bin/touch $output/progress/$pre\_$uID\_UG_MAF.done`;
    ###}
}

if($pair){
    ### QC pairing
    if(!-e "$output/progress/$pre\_$uID\_FP_WRAP.done" || $ran_vf_hc){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_FP_WRAP", job_hold => "$vfhcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_FP_WRAP.log");
	my $standardParams = Schedule::queuing(%stdParams);

	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/pairing_qc.pl -pair $pair -vcf $output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf -pre $pre -outdir $output/metrics/fingerprint -scheduler $scheduler -config $config -progress $output/progress $vfhcj`;
	`/bin/touch $output/progress/$pre\_$uID\_FP_WRAP.done`;
    }
    
    `/bin/mkdir -m 775 -p $output/variants/mutect`;
    `/bin/mkdir -m 775 -p $output/variants/somaticsniper`;
    `/bin/mkdir -m 775 -p $output/variants/virmid`;
    `/bin/mkdir -m 775 -p $output/variants/scalpel`;
    `/bin/mkdir -m 775 -p $output/variants/strelka`;
    ###`/bin/mkdir -m 775 -p $output/variants/varscan`;

    open(PAIR, "$pair") or die "Can't open $pair file";
    my %submitted_lns = ();
    my @mu_jids = ();
    my $ran_mutect_glob = 0;
    ### NOTE: THE SAMPLE NAMES IN THE PAIRING FILE MUST MATCH EXACTLY THE SAMPLE NAMES IN THE REALIGNED/RECALIBRATED BAM FILE
    while(<PAIR>){
	chomp;
	
	my @data = split(/\s+/, $_);
	### pairing file also contains unpaired samples with NA for the other pai
	### so we can keep track of what sample is tumor or normal
	if($data[0] =~ /^NA$/i || $data[1] =~ /^NA$/i){
	    next;
	}
	
	### NOTE: NOT AUTOMATICALLY RUNNING EACH SOMATIC CALLER WHEN SSF RUNS
	###       BCAUSE DON'T WANT TO KICK IT OFF FOR ALL SAMPLE PAIRS
	###       IN CASE ONLY ONE SAMPLE GROUP GOT MODIFIED ABOVE

	###       WILL RUN SOMATIC ANALYSIS FOR ALL SAMPLE PAIRS
	###       IF JUST ONE SAMPLE HAD ITS BAM MODIFIED
	my $ran_mutect = 0;
	my $mutectj = '';
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done" || $ran_ssf){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT", job_hold => "$ssfj", cpu => "2", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xmx4g -Djava.io.tmpdir=/scratch/$uID -jar $MUTECT/muTect.jar --analysis_type MuTect --reference_sequence $REF_SEQ --dbsnp $Bin/data/dbsnp_135.hg19__ReTag.vcf --cosmic $Bin/data/CosmicCodingMuts_v67_20131024.vcf --input_file:normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --input_file:tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --vcf $output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf --out $output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.txt -rf BadCigar --enable_extended_output --downsampling_type NONE`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done`;
	    $mutectj = "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT";
	    push @mu_jids, "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT";
	    $ran_mutect = 1;
	    $ran_mutect_glob = 1;
	}

	my $ran_somatic_sniper = 0;
	my $ssj = '';
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.done" || $ran_ssf){  
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER", job_hold => "$ssfj", cpu => "2", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $SOMATIC_SNIPER/bam-somaticsniper -F vcf -f $REF_SEQ -q 1 $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam $output/variants/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.done`;
	    $ssj = "$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER";
	    $ran_somatic_sniper = 1;
	}	    

	###my $ran_mpileup = 0;
	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MPILEUP.done" || $ran_ssf){  
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $ssfj -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -f $REF_SEQ -d 500000 $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam ">$output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup"`;

	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $ssfj -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -f $REF_SEQ -d 500000 $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam ">$output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup"`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MPILEUP.done`;
	    ###$ran_mpileup = 1;
	###}	    

	### BECAUSE THE OUTPUT PREFIX FOR VIRMID IS THE DISEASE BAM,
	### NEED TO CREATE A DIRECTORY FOR EACH TUMOR/NORMAL PAIR
	### SINCE A DISEASE BAM CAN BE USED IN MULTIPLE COMPARISONS
	### virmid fails if directory already exists
	my $ran_virmid = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done" || $ran_ssf){  
	    sleep(3);
	    if(-d "$output/variants/virmid/$data[0]\_$data[1]\_virmid"){
		`/bin/rm -rf $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    }
	    `/bin/mkdir -m 775 -p $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_VIRMID", job_hold => "$ssfj", cpu => "4", mem => "12", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx12g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $VIRMID/Virmid.jar -R $REF_SEQ -D $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam -N $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam -t 4 -o $pre\_$data[0]\_$data[1]\_virmid -w $output/variants/virmid/$data[0]\_$data[1]\_virmid`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done`;
	    $ran_virmid = 1;
	}

	my $ran_strelka_config = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done" || $ran_ssf){  
	    if(-d "$output/variants/strelka/$data[0]\_$data[1]\_strelka"){
		### strelka DIES IF DIR ALREADY EXISTS
		`/bin/rm -rf $output/variants/strelka/$data[0]\_$data[1]\_strelka`;
	    }

	    my @lns_jids = ();
	    ### NOTE: Strelka only recognizes X.bam.bai as the index for X.bam, not X.bai
	    if((!-e "$output/progress/$pre\_indelRealigned_recal\_$data[0]\_LNS.done" || $ran_ssf)  && !-e "$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam.bai" && !$submitted_lns{$data[0]}){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_indelRealigned_recal\_$data[0]\_LNS", job_hold => "$ssfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_indelRealigned_recal\_$data[0]\_LNS.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/ln -s $pre\_indelRealigned_recal\_$data[0]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam.bai`;
		#push @lns_jids, "$pre\_indelRealigned_recal\_$data[0]\_LNS";
		$submitted_lns{$data[0]} = 1;
		`/bin/touch $output/progress/$pre\_indelRealigned_recal\_$data[0]\_LNS.done`;
	    }

	    if((!-e "$output/progress/$pre\_indelRealigned_recal\_$data[1]\_LNS.done" || $ran_ssf) && !-e "$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam.bai" && !$submitted_lns{$data[1]}){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_indelRealigned_recal\_$data[1]\_LNS", job_hold => "$ssfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_indelRealigned_recal\_$data[1]\_LNS.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/ln -s $pre\_indelRealigned_recal\_$data[1]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam.bai`;
		#push @lns_jids, "$pre\_indelRealigned_recal\_$data[1]\_LNS";
		$submitted_lns{$data[1]} = 1;
		`/bin/touch $output/progress/$pre\_indelRealigned_recal\_$data[1]\_LNS.done`;
	    }

	    if($submitted_lns{$data[0]}){
		push @lns_jids, "$pre\_$uID\_indelRealigned_recal\_$data[0]\_LNS";
	    }
	    if($submitted_lns{$data[1]}){
		push @lns_jids, "$pre\_$uID\_indelRealigned_recal\_$data[1]\_LNS";
	    }
	    my $lnsj = join(",", @lns_jids, $ssfj);
	    ### NOTE: strelka ONLY HAS CONFIG FOR BWA ALN, NOT SURE HOW IT WILL WORK WITH BWA MEM
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", job_hold => "$lnsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STRELKA/bin/configureStrelkaWorkflow.pl --normal=$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor=$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --ref=$REF_SEQ --config=$STRELKA/etc/strelka_config_bwa_default.ini --output-dir=$output/variants/strelka/$data[0]\_$data[1]\_strelka`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done`;
	    $ran_strelka_config = 1;
	}

	my $ran_scalpel = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done" || $ran_ssf){
	    sleep(3);
	    `/bin/mkdir -m 775 -p $output/variants/scalpel/$data[0]\_$data[1]\_scalpel`;
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL", job_hold => "$ssfj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $SCALPEL/scalpel --somatic --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --bed $targets_bed_padded --ref $REF_SEQ --dir $output/variants/scalpel/$data[0]\_$data[1]\_scalpel --numprocs 24`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done`;
	    $ran_scalpel = 1;
	}

	my $ran_strelka_run = 0;
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done" || $ran_strelka_config){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", cpu => "8", mem => "16", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /usr/bin/make -C $output/variants/strelka/$data[0]\_$data[1]\_strelka -j 8`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done`;
	    $ran_strelka_run = 1;
	}

	###my $ran_varscan_somatic = 0;
	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC.done" || $ran_mpileup){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -pe alloc 2 -l virtual_free=5G -q lau.q,lcg.q,nce.q $Bin/qCMD $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $VARSCAN/VarScan.jar somatic $output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup $output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic --strand-filter 1 --output-vcf 1`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC.done`;
	    ###$ran_varscan_somatic = 1;
	###}

	###sleep(3);

	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done" || $ran_mutect){
	   ### &generateMaf("$output/variants/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf", 'mutect', "$mutectj", $data[0], $data[1]);
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done`;
	###}

	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.done" || $ran_somatic_sniper){
	   ### my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF", job_hold => "$ssj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.log");
	    ###my $standardParams = Schedule::queuing(%stdParams);
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf -c somaticsniper -o $output/variants/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.done`;
	###}

	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done" || $ran_varscan_somatic){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp_MAF.txt -n $data[0] -t $data[1]`;

	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.indel.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.indel_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done`;
	###}

	if($ran_scalpel){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL_CLEANUP", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL", cpu => "1", mem => "1");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/rm -rf $output/variants/scalpel/$data[0]\_$data[1]\_scalpel/main/ $output/variants/scalpel/$data[0]\_$data[1]\_scalpel/validation/`;
	}

	if($ran_strelka_run){
	    sleep(3);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/rm -rf $output/variants/strelka/$data[0]\_$data[1]\_strelka/config $output/variants/strelka/$data[0]\_$data[1]\_strelka/chromosomes $output/variants/strelka/$data[0]\_$data[1]\_strelka/Makefile $output/variants/strelka/$data[0]\_$data[1]\_strelka/task.complete`;
	}

	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.done" || $ran_strelka_run){
	   ### my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.log");
	    ###my $standardParams = Schedule::queuing(%stdParams);
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.snvs.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
       
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.snvs.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;

	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.indels.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;

	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.indels.vcf -c strelka -o $output/variants/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.done`;
	    ###	}
    }
    close PAIR;

    if(!-e "$output/progress/$pre\_HAPLOTECT.done" || $ran_mutect_glob || $ran_vf_hc){
	sleep(3);
	my $muj = join(",", @mu_jids);
	my %addParams_I = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "ito.q");
	my $additionalParams_I = Schedule::additionalParams(%addParams_I);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HAPLOTECT", job_hold => "$vfhcj,$muj", cpu => "4", mem => "8", internet => "1", cluster_out => "$output/progress/$pre\_HAPLOTECT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{internet} $additionalParams_I $PERL/perl $Bin/haploTect_merge.pl -pair $pair -hc_vcf $output/variants/haplotypecaller/$pre\_HaplotypeCaller.vcf -pre $pre -output $output/variants/haplotect -mutect_dir $output/variants/mutect -config $config -delete_temp`;
	`/bin/touch $output/progress/$pre\_HAPLOTECT.done`;
    }
}

open(GROUP, "$group") || die "CAN'T OPEN SAMPLE GROUPING FILE $group $!";
open(BLIST, ">$output/intFiles/$pre\_bam_list.txt") || die "CAN'T WRITE TO SAMPLE BAM LIST FILE $output/intFiles/$pre\_bam_list.txt $!";
while(<GROUP>){
    chomp;

    my @data = split(/\s+/, $_);
    print BLIST "$data[0]\t$curDir/$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam\n";
}
close GROUP;
close BLIST;

if(!-e "$output/progress/$pre\_$uID\_STRVAR.done" || $ran_ssf){
    sleep(3);
    if(-d "$curDir/$output/strvar"){
	`/bin/rm -rf $curDir/$output/strvar`;
    }

    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_STRVAR", job_hold => "$ssfj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_STRVAR.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/RunStructuralVariantPipeline_Delly.pl -pre $pre -pair $pair -bam_list $curDir/$output/intFiles/$pre\_bam_list.txt -out $curDir/$output/strvar -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;
    `/bin/touch $output/progress/$pre\_$uID\_STRVAR.done`;
}

sub generateMaf{
    my ($vcf, $type, $hold, $normal_sample, $tumor_sample) = @_;

    my $n_sample = '';
    my $t_sample = '';
    if($normal_sample && $tumor_sample){
	$n_sample = "-normal_sample $normal_sample";
	$t_sample = "-tumor_sample $tumor_sample";
    }

    my $jna = $vcf;
    $jna =~ s/\//_/g;

    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "ito.q");
    my $additionalParams = Schedule::additionalParams(%addParams);

    if($pair){
	### NOTE: ASKING FOR PE ALLOC 4 TO THROTTLE NUMBER OF JOBS ACCESING ONCOTATOR
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_MAF_PAIRED", job_hold => "$hold", cpu => "4", mem => "8", internet => "1", cluster_out => "$output/progress/$pre\_$uID\_$jna\_MAF_PAIRED.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{internet} $additionalParams $Bin/generateMAF.pl -vcf $vcf -pairing $pair -species hg19 -config $config -caller $type $n_sample $t_sample -delete_temp`;

	if($type =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_MAF_UNPAIRED", job_hold => "$hold", cpu => "4", mem => "8", internet => "1", cluster_out => "$output/progress/$pre\_$uID\_$jna\_MAF_UNPAIRED.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{internet} $additionalParams $Bin/generateMAF.pl -vcf $vcf -species hg19 -config $config -caller $type -delete_temp`;
	}
    }
    else{
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_MAF_UNPAIRED", job_hold => "$hold", cpu => "4", mem => "8", internet => "1", cluster_out => "$output/progress/$pre\_$uID\_$jna\_MAF_UNPAIRED.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $standardParams->{internet} $additionalParams $Bin/generateMAF.pl -vcf $vcf -species hg19 -config $config -caller $type -delete_temp`;
    }
}