#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

use lib "$Bin/lib";
use Schedule;
use Cluster;
use POSIX qw(strftime);
my $cur_date = strftime "%Y%m%d", localtime;

my ($pair, $svnRev, $email, $patient, $group, $bamgroup, $config, $nosnps, $targets, $ug, $scheduler, $priority_project, $priority_group, $abra, $indelrealigner, $help, $step1, $DB_SNP, $allSomatic, $scalpel, $somaticsniper, $strelka, $varscan, $virmid, $lancet, $vardict, $pindel, $abra_target);

my $pre = 'TEMP';
my $output = "results";
my $species = 'species_custom';

my $uID = `/usr/bin/id -u -n`;
chomp $uID;
my $rsync = "/ifs/res/$uID";
my $tempdir = "/scratch/$uID";

GetOptions ('email=s' => \$email,
	    'pre=s' => \$pre,
	    'pair=s' => \$pair,
            'patient=s'=> \$patient,
	    'group=s' => \$group,
	    'config=s' => \$config,
	    'targets=s' => \$targets,
	    'nosnps' => \$nosnps,
	    'step1' => \$step1,
	    'ug|unifiedgenotyper' => \$ug,
	    'abra' => \$abra,
	    'indelrealigner' => \$indelrealigner,
	    'bamgroup=s' => \$bamgroup,
 	    'scheduler=s' => \$scheduler,
            'svnRev=s' => \$svnRev,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
 	    'db_snp|dbsnp=s' => \$DB_SNP,
	    'help' => \$help,
	    'rsync=s' => \$rsync,
	    'allsomatic|allSomatic|all_somatic' => \$allSomatic,
	    'scalpel' => \$scalpel,
	    'somaticsniper' => \$somaticsniper,
	    'strelka' => \$strelka,
	    'varscan' => \$varscan,
	    'virmid' => \$virmid,
	    'lancet' => \$lancet,
 	    'vardict' => \$vardict,
            'pindel' => \$pindel,
 	    'tempdir=s' => \$tempdir,
 	    'output|out|o=s' => \$output,
            'abratarget|abra_target=s' => \$abra_target) or exit(1);


if(!$group || !$config || !$scheduler || !$targets || !$bamgroup || $help){
    print <<HELP;

    USAGE: process_alignments_hg19.pl -group GROUP -pair PAIR -pre PRE -config CONFIG -species SPECIES -scheduler SCHEDULER -targets TARGETS
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED)
	* BAMGROUP: files listing bams to be processed together; every bam for each group on 1 line, comma-separated (required)
	* TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS unless given full path to targets directory (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE, LUNA, and JUNO (REQUIRED)
	* DB_SNP: VCF file containing known sites of genetic variation (REQUIRED; either here or in config)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
	* PRE: output prefix (default: TEMP)
	* OUTPUT: output results directory (default: results)
	* RSYNC:  path to rsync data for archive (default: /ifs/res/USER_ID)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* -nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed
	* -abra: run abra to realign indels (default)
	* -indelrealigner: run GATK indelrealigner (abra runs as default)
	* -step1: forece the pipeline to start from the first step in pipeline
	* haplotypecaller is default; -ug || -unifiedgenotyper to also make unifiedgenotyper variant calls	
	* TEMPDIR:  temp directory (default: /scratch/$uID)
	* ALLSOMATIC: run all somatic callers; mutect/haplotypecaller always run; otherwise -scalpel, -somaticsniper, -strelka, -varscan, -virmid, -lancet, -vardict, -pindel to run them individually	
HELP
exit;
}

my $curDir = `pwd`;
chomp $curDir;
if($output !~ /^\//){
    $output = "$curDir/$output";
}

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

if($abra && $indelrealigner){
    die "Cannot run both abra and gatk indelrealigner $!";
}
elsif(!$abra && !$indelrealigner){
    $abra = 1;
}

if($allSomatic){
    $scalpel = 1;
    $somaticsniper = 1;
    $strelka = 1;
    $varscan = 1;
    $virmid = 1;
    $lancet = 1;
    $vardict = 1;
    $pindel = 1;
}

my $ABRA = '';
my $GATK = '';
my $LANCET = '';
my $MUTECT = '';
my $PICARD = '';
my $PINDEL = '';
my $SAMTOOLS = '';
my $SCALPEL = '';
my $SOMATIC_SNIPER = '';
my $STRELKA = '';
my $VARDICT_JAVA = '';
my $VARDICT_PERL = '';
my $VARSCAN = '';
my $VIRMID = '';

my $JAVA = '';
my $JAVA7_MUTECT = '';
my $PYTHON = '';
my $PERL = '';

my $SINGULARITY = '';
my $singularityParams = '';
my $singularityBind = '';
my $singularityenv_prepend_path = "";

my $REF_SEQ = "";
my $REF_FAI = '';
my $BWA_INDEX = "";

open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /singularity/i){
        if(!-e "$conf[1]/singularity"){
            die "CAN'T FIND singularity IN $conf[1] $!";
        }
        $SINGULARITY = $conf[1];
    }
    elsif($conf[0] =~ /abra/i){
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
    elsif($conf[0] =~ /lancet/i){
	if(!-e "$conf[1]/lancet"){
	    die "CAN'T FIND lancet IN $conf[1] $!";
	}
	$LANCET = $conf[1];
    }
    elsif($conf[0] =~ /pindel/i){
        if(!-e "$conf[1]/pindel" or !-e "$conf[1]/pindel2vcf"){
            die "CAN'T FIND pindel or pindel2vcf IN $conf[1] $!";
        }
        $PINDEL = $conf[1];
    }
    elsif($conf[0] =~ /^mutect$/i){
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
    elsif($conf[0] =~ /scalpel/i){
	if(!-e "$conf[1]/scalpel"){
	    die "CAN'T FIND scalpel IN $conf[1] $!";
	}
	$SCALPEL = $conf[1];
    }
    elsif($conf[0] =~ /somaticsniper/i){
	if(!-e "$conf[1]/bam-somaticsniper"){
	    die "CAN'T FIND bam-somaticsniper IN $conf[1] $!";
	}
	$SOMATIC_SNIPER = $conf[1];
    }
    elsif($conf[0] =~ /strelka/i){
	if(!-e "$conf[1]/bin/configureStrelkaWorkflow.pl"){
	    die "CAN'T FIND bin/configureStrelkaWorkflow.pl IN $conf[1] $!";
	}
	$STRELKA = $conf[1];
    }
    elsif($conf[0] =~ /vardict_java/i){
	if(!-e "$conf[1]/VarDict"){
	    die "CAN'T FIND VarDict IN $conf[1] $!";
	}
	$VARDICT_JAVA = $conf[1];
    }
    elsif($conf[0] =~ /vardict_perl/i){
	if(!-e "$conf[1]/testsomatic.R" || !-e "$conf[1]/var2vcf_paired.pl"){
	    die "CAN'T FIND testsomatic.R OR var2vcf_paired.pl IN $conf[1] $!";
	}
	$VARDICT_PERL = $conf[1];
    }
    elsif($conf[0] =~ /varscan/i){
	if(!-e "$conf[1]/VarScan.jar"){
	    die "CAN'T FIND VarScan.jar IN $conf[1] $!";
	}
	$VARSCAN = $conf[1];
    }
    elsif($conf[0] =~ /virmid/i){
	if(!-e "$conf[1]/Virmid.jar"){
	    die "CAN'T FIND Virmid.jar IN $conf[1] $!";
	}
	$VIRMID = $conf[1];
    }
    elsif($conf[0] =~ /^java$/i){
	if(!-e "$conf[1]/java"){
	    die "CAN'T FIND java IN $conf[1] $!";
	}
	$JAVA = $conf[1];
    }
    elsif($conf[0] =~ /^java7_mutect$/i){ 
        if(!-e "$conf[1]/java"){
            die "CAN'T FIND java IN $conf[1] $!";
        }    
        $JAVA7_MUTECT = $conf[1];
    }
    elsif($conf[0] =~ /^perl$/i){
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
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
        $singularityenv_prepend_path .= ":$conf[1]";
    }
    elsif($conf[0] =~ /^r$/i){
	if(!-e "$conf[1]/R"){
	    die "CAN'T FIND R IN $conf[1] $!";
	}
	my $path_tmp = $ENV{'PATH'};
	$ENV{'PATH'} = "$conf[1]:$path_tmp";
        $singularityenv_prepend_path .= ":$conf[1]";
    }
    elsif($conf[0] =~ /species_custom_fasta/i){
	if(!-e "$conf[1]"){
	    die "CAN'T FIND $conf[1] $!";
	}
	$REF_SEQ = $conf[1];
    }
    elsif($conf[0] =~ /species_custom_fai/i){
	if(!-e "$conf[1]"){
	    die "CAN'T FIND $conf[1] $!";
	}
	$REF_FAI = $conf[1];
    }
    elsif($conf[0] =~ /species_custom_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR SPECIES CUSTOM WITH PREFIX $conf[1] $!";
	}
	$BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /species_custom_db_snp/i){
	if(!-e "$conf[1]"){
	    if($DB_SNP){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	if($DB_SNP && ($DB_SNP ne $conf[1])){
	    die "INPUT PARAM FOR DB_SNP $DB_SNP IS NOT THE SAME AS THAT IN THE CONFIG FILE $conf[1] $!";
	}
	$DB_SNP = $conf[1];
    }
}
my %sinParams = (singularity_exec => "$SINGULARITY/singularity", singularity_image => "$Bin/variants_pipeline_singularity_prod.simg");
$singularityParams = Schedule::singularityParams(%sinParams);
$singularityBind = Schedule::singularityBind();

$ENV{'SINGULARITYENV_PREPEND_PATH'} = $singularityenv_prepend_path;
$ENV{'SINGULARITY_BINDPATH'} = $singularityBind; 
close CONFIG;


my $ABRA_TARGETS = '';

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
if(-d $targets){
    my @path = split(/\//, $targets);
    my $assay = pop @path;
    $targets_bed_padded = "$targets/$assay\_targets_plus5bp.bed";
}
if(!-e "$targets_bed_padded"){
    die "CAN'T LOCATE $targets_bed_padded FOR $targets; REQUIRED FOR SCALPEL $!";
}
if(!$abra_target){
    $abra_target = $targets_bed_padded;
}
if(!-e "$abra_target"){
    die "CAN'T LOCATE $abra_target; REQUIRED FOR ABRA $!";
}

my $multipleTargets = '';

###if($target){
###   if(!-e $target){
###	die "target file $target cannot be found $!";
###    }
###    $multipleTargets = "-L $target_bed --interval_set_rule INTERSECTION";
###}

my @ref_chr = ();
open(FAI, "$REF_FAI") or die "CAN'T OPEN FAI FILE $REF_FAI $!";
while(<FAI>){
    chomp;

    my @line = split(/\s+/, $_);
    push @ref_chr, $line[0];
}
close FAI;

my $count = 0;
my %inputFiles = ();
my @processedBams = ();
my @finalBams = ();
my $ran_pr_glob = 0;
my @prg_jids = ();
my $ran_ssf = 0;
my @ssf_jids = ();
my @all_jids = ();

if(!-d $output){
    mkdir("$output", 0775) or die "Can't make $output";
}
if(!-d "$output/intFiles"){
    mkdir("$output/intFiles", 0775) or die "Can't make $output/intFiles";
}
if(!-d "$output/alignments"){
    mkdir("$output/alignments", 0775) or die "Can't make $output/alignments";
}
if(!-d "$output/progress"){
    mkdir("$output/progress", 0775) or die "Can't make $output/progress";
}
if(!-d $tempdir){
    mkdir("$tempdir", 0775) or die "Can't make $tempdir";
}

my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);

open(IN, "$bamgroup") || die "CAN'T OPEN GROUPING FILE OF MARKDUP BAMS $bamgroup $!";
while(<IN>){
    chomp;
    
    my @gpair = split(/\s+/, $_);
    my @pair = split(/,/, $gpair[1]);
    my @pins = ();
    foreach my $pai (@pair){
        my @sp = split(/intFiles\//, $pai);
        my @sn = split(/\//, $sp[1]);
	if($inputFiles{$pai}){
	    next;
	}
	push @pins, "-I $pai";
	$inputFiles{$pai} = 1;
        
        my $samp = $sn[0];
        push @finalBams, "$output/alignments/$pre\_indelRealigned_recal_$samp.bam";
    }

    if(scalar(@pins) == 0){
	next;
    }

    my $bgroup = join(" ", @pins);
    my $indelBam = '';
    my $ran_ir == 0;
    my $irj = '';

    if($abra){
	my @inBams = ();
	my @outBams = ();
	foreach my $pin (@pins){
	    my @inB = split(/\s+/, $pin);
	    push @inBams, $inB[1];
	    push @outBams, "$inB[1]\_ABRA.bam";
	    $indelBam = "$inB[1]\_ABRA.bam\_FM.bam";
	}
	
	my $aiBams = join(",", @inBams);
	my $aoBams = join(",", @outBams);
	my $ran_abra = 0;
	my $abraj = '';
	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_ABRA.done" || $step1){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_ABRA", cpu => "12", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_ABRA.log");
	    my $standardParams = Schedule::queuing(%stdParams);	    
	    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "4");
	    my $additionalParams = Schedule::additionalParams(%addParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/abra_wrapper.pl -inBams $aiBams -outBams $aoBams -refSeq $REF_SEQ -bwaRef $BWA_INDEX -targets $abra_target -working $output/intFiles/abra_$gpair[0] -config $config -tempdir $tempdir -log $output/progress/$pre\_$uID\_$gpair[0]\_ABRA_WRAPPER.log`;

	    $abraj = "$pre\_$uID\_$gpair[0]\_ABRA";
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_ABRA.done`;
	    $ran_abra = 1;
	}

	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_FIXMATE.done" || $ran_abra){
	    sleep(2);
	    my $bcount = 0;
	    foreach my $outBam (@outBams){
		$bcount++;
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$bcount\_FIXMATE", job_hold => "$abraj", cpu => "1", mem => "50", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_FIXMATE.log");
		my $standardParams = Schedule::queuing(%stdParams);
		my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "5");
		my $additionalParams = Schedule::additionalParams(%addParams);

		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar FixMateInformation I=$outBam O=$outBam\_FM.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true`;
		$irj = "$pre\_$uID\_$gpair[0]\_$bcount\_FIXMATE";
	    }
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_FIXMATE.done`;
	    $ran_ir = 1;
	}
    }
    else{
	my $ran_rtc = 0;
	my $rtc_jid = '';
	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_RTC.done" || $step1){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_RTC", cpu => "10", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_RTC.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ $multipleTargets -S LENIENT --known $DB_SNP -nt 10 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_indelRealigner.intervals $bgroup`;
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_RTC.done`;
	    $rtc_jid = "$pre\_$uID\_$gpair[0]\_RTC";
	    $ran_rtc = 1;
	}
	
	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_IR.done" || $ran_rtc){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_IR", job_hold => "$rtc_jid", cpu => "1", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_IR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ $multipleTargets -S LENIENT --knownAlleles $DB_SNP --targetIntervals $output/intFiles/$pre\_$gpair[0]\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_indelRealigned.bam $bgroup`;
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_IR.done`;
	    $irj = "$pre\_$uID\_$gpair[0]\_IR";
	    $ran_ir = 1;
	}
	
	$indelBam = "$output/intFiles/$pre\_$gpair[0]\_indelRealigned.bam";
    }
	
    my $ran_br = 0;
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_BR.done" || $ran_ir){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_BR", job_hold => "$irj", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_BR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO -R $REF_SEQ -S LENIENT --knownSites $DB_SNP --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 6 --out $output/intFiles/$pre\_$gpair[0]\_recal_data.grp -I $indelBam`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_BR.done`;
	$ran_br = 1;
    }
        
    my $ran_pr = 0;
    my $prj = '';
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_PR.done" || $ran_br){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_PR", job_hold => "$pre\_$uID\_$gpair[0]\_BR", cpu => "12", mem => "40", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_PR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R $REF_SEQ $multipleTargets --emit_original_quals -BQSR $output/intFiles/$pre\_$gpair[0]\_recal_data.grp --num_cpu_threads_per_data_thread 12 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam -I $output/intFiles/$pre\_$gpair[0]\_indelRealigned.bam`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_PR.done`;
	$prj = "$pre\_$uID\_$gpair[0]\_PR";
	push @prg_jids, "$pre\_$uID\_$gpair[0]\_PR";
	$ran_pr = 1;
	$ran_pr_glob = 1;
    }
    
    push @processedBams, "-I $output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam";
        
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.done" || $ran_pr){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_SSF", job_hold => "$prj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T SplitSamFile -R $REF_SEQ -I $output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam --outputRoot $output/alignments/$pre\_indelRealigned_recal_`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_SSF.done`;
	push @ssf_jids, "$pre\_$uID\_$gpair[0]\_SSF";
	$ran_ssf = 1;
    }
}

my $ssfj = join(",", @ssf_jids);
push @all_jids, @ssf_jids;
my @mq_metrics_jid = ();
my $ran_mqm = 0;
foreach my $finalBam (@finalBams){
    my @sn = split(/\//, $finalBam);
    my $samp = $sn[-1];
    $samp =~ s/\.bam//g;
    $samp =~ s/$pre\_indelRealigned_recal_//g;
    
    if(!-e "$output/progress/$pre\_$uID\_MQ_METRICS_$samp\.done" || $ran_ssf){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MQ_METRICS_$samp", job_hold => "$ssfj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_MQ_METRICS_$samp\.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MeanQualityByCycle INPUT=$finalBam OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.txt CHART_OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.pdf REFERENCE_SEQUENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true TMP_DIR=$tempdir`;
	push @mq_metrics_jid, "$pre\_$uID\_MQ_METRICS_$samp";
	`/bin/touch $output/progress/$pre\_$uID\_MQ_METRICS_$samp\.done`;
	$ran_mqm = 1; 
    }
}

my $mqmj = join(",", @mq_metrics_jid);
my $ran_mmqm = 0;
if(!-e "$output/progress/$pre\_$uID\_MERGE_MQ.done" || $ran_mqm){
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_MQ", job_hold => "$mqmj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_MQ.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/mergeMeanQualityHistograms.py $output '*_MeanQualityByCycle_*.txt' $output/metrics/$pre\_post_recal_MeanQualityByCycle.txt $output/metrics/$pre\_pre_recal_MeanQualityByCycle.txt`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_MQ.done`;
    push @all_jids, "$pre\_$uID\_MERGE_MQ";
    $ran_mmqm = 1;
}

my $allj = join(",", @all_jids);
if(!-e "$output/progress/$pre\_$uID\_RSYNC_1.done" || $ran_ssf || $ran_mmqm){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC_1", job_hold => "$allj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_1.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' --exclude 'variants' --exclude 'metrics' $curDir $rsync`;
    push @all_jids, "$pre\_$uID\_RSYNC_1";
    `/bin/touch $output/progress/$pre\_$uID\_RSYNC_1.done`;
}

if($nosnps){
    exit;
}

if(!-d "$output/variants"){
    mkdir("$output/variants", 0775) or die "Can't make $output/variants";
}
if(!-d "$output/variants/snpsIndels"){
    mkdir("$output/variants/snpsIndels", 0775) or die "Can't make $output/variants/snpsIndels";
}
if(!-d "$output/variants/snpsIndels/haplotypecaller"){
    mkdir("$output/variants/snpsIndels/haplotypecaller", 0775) or die "Can't make $output/variants/snpsIndels/haplotypecaller";
}
my $ran_hc = 0;
my $ran_ug_snp = 0;
my $ran_ug_indel = 0;
my $hcj = '';
my $ugsj = '';
my $ugij = '';
my $prgj = join(",", @prg_jids);

my $irBams2 = join(" ", @processedBams);    
if($ug){
    if(!-e "$output/progress/$pre\_$uID\_UG_SNP.done" || $ran_pr_glob){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_UG_SNP", job_hold => "$prgj", cpu => "12", mem => "48", cluster_out => "$output/progress/$pre\_$uID\_UG_SNP.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name SPECIES_CUSTOM $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model SNP --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf $irBams2`;
	`/bin/touch $output/progress/$pre\_$uID\_UG_SNP.done`;
	$ugsj = "$pre\_$uID\_UG_SNP";
	$ran_ug_snp = 1;
    }
    
    if(!-e "$output/progress/$pre\_$uID\_UG_INDEL.done" || $ran_pr_glob){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_UG_INDEL", job_hold => "$prgj", cpu => "12", mem => "48", cluster_out => "$output/progress/$pre\_$uID\_UG_INDEL.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name SPECIES_CUSTOM $multipleTargets --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model INDEL --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf $irBams2`;
	`/bin/touch $output/progress/$pre\_$uID\_UG_INDEL.done`;
	$ugij = "$pre\_$uID\_UG_INDEL";
	$ran_ug_indel = 1;
    }    
}

if(!-e "$output/progress/$pre\_$uID\_HC.done" || $ran_pr_glob){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HC", job_hold => "$prgj", cpu => "30", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "7");
	my $additionalParams = Schedule::additionalParams(%addParams);
    my $response = "";
    my $HC_submission_num = 0;
    while($response !~ /is submitted to default queue/ && $HC_submission_num < 10){
        my $time = 60 * $HC_submission_num;
        `sleep $time `;
        $HC_submission_num++;
        print "HC attempt number $HC_submission_num\n"; 
        print "Command: $standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 30 --out $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf $irBams2";
        $response = `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 30 --out $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf $irBams2       2>&1`;
        print "Response: $response\n";
        if($response !~ /is submitted to default queue/ && $HC_submission_num == 10){
            `mail -s "HaplotypeCaller Job Not Working" $email <<< "This is an e-mail letting you know that your variants pipeline run for $pre is not submitting HC jobs to the cluster. Please remove all .done files for haplotype caller and rerun the pipeline."  `;
        }
    }
    `/bin/touch $output/progress/$pre\_$uID\_HC.done`;
    $hcj = "$pre\_$uID\_HC";
    $ran_hc = 1;
}

if(!-e "$output/progress/$pre\_$uID\_MAF_HC.done" || $ran_hc){
    sleep(2);
    &generateMaf("$output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'haplotypecaller', "$hcj");
    `/bin/touch $output/progress/$pre\_$uID\_MAF_HC.done`;
}

if($ug){
    if(!-d "$output/variants/snpsIndels/unifiedgenotyper"){
	mkdir("$output/variants/snpsIndels/unifiedgenotyper", 0775) or die "Can't make $output/variants/snpsIndels/unifiedgenotyper";
    }
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_RAW.done" || $ran_ug_snp || $ran_ug_indel){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_UG_RAW", job_hold => "$ugsj,$ugij", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG_RAW.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf --assumeIdenticalSamples --variant $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf --variant $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_UG_RAW.done`;
    }

    my $ran_vf_ug_snp = 0;
    my @vfug_jids = ();
    if(!-e "$output/progress/$pre\_$uID\_VF_UG_SNP.done" || $ran_ug_snp){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VF_UG_SNP", job_hold => "$ugsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_VF_UG_SNP.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf --maskName nearIndel --variant $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf -o $output/intFiles/$pre\_UnifiedGenotyper_SNP_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"MQ \\< 40.0\\" --filterExpression \\"FS \\> 60.0\\" --filterExpression \\"HaplotypeScore \\> 13.0\\" --filterExpression \\"MQRankSum \\< -12.5\\" --filterExpression \\"ReadPosRankSum \\< -8.0\\" --filterName QDFilter --filterName MQFilter --filterName FSFilter --filterName HSFilter --filterName MQRSFilter --filterName ReadPosFilter`;
	`/bin/touch $output/progress/$pre\_$uID\_VF_UG_SNP.done`;
	push @vfug_jids, "$pre\_$uID\_VF_UG_SNP";
	$ran_vf_ug_snp = 1;
    }

    my $ran_vf_ug_indel = 0;
    if(!-e "$output/progress/$pre\_$uID\_VF_UG_INDEL.done" || $ran_ug_indel){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VF_UG_INDEL", job_hold => "$ugij", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_VF_UG_INDEL.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --variant $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf -o $output/intFiles/$pre\_UnifiedGenotyper_INDEL_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"ReadPosRankSum \\< -20.0\\" --filterExpression \\"InbreedingCoeff \\< -0.8\\" --filterExpression \\"FS \\> 200.0\\" --filterName QDFilter --filterName ReadPosFilter --filterName InbreedingFilter --filterName FSFilter`;
	`/bin/touch $output/progress/$pre\_$uID\_VF_UG_INDEL.done`;
	push @vfug_jids, "$pre\_$uID\_VF_UG_INDEL";
	$ran_vf_ug_indel = 1;
    }

    my $ran_cv_ug_si = 0;
    my $cvugsij = '';
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_SI.done" || $ran_ug_snp || $ran_ug_indel){
	sleep(2);
	my $vfugj = join(",", @vfug_jids);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_UG_SI", job_hold => "$vfugj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG_SI.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/snpsIndels/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf --assumeIdenticalSamples --variant $output/progress/$pre\_UnifiedGenotyper_SNP_vf.vcf --variant $output/progress/$pre\_UnifiedGenotyper_INDEL_vf.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_SI.done`;
        $cvugsij = "$pre\_$uID\_CV_UG_SI";
	$ran_cv_ug_si = 1;
    }
        
    if(!-e "$output/progress/$pre\_$uID\_MAF_UG.done" || $ran_cv_ug_si){  
	sleep(2);
	&generateMaf("$pre\_UnifiedGenotyper.vcf", 'unifiedgenotyper', "$cvugsij");
 	`/bin/touch $output/progress/$pre\_$uID\_MAF_UG.done`;
    }
}

if($pair){
    if(!-d "$output/variants/snpsIndels/haplotect"){
	mkdir("$output/variants/snpsIndels/haplotect", 0775) or die "Can't make $output/variants/snpsIndels/haplotect";
    }

    open(PAIR, "$pair") or die "Can't open $pair file";
    my %submitted_lns = ();
    my @mu_jids = ();
    my $ran_mutect_glob = 0;
    my $hasPair = 0;
    my $haplotect_run=0;
    while(<PAIR>){
	chomp;
	
	my @data = split(/\s+/, $_);
	
	if($data[0] =~ /^NA$/i || $data[1] =~ /^NA$/i){
	    next;
	}
	
        $hasPair = 1;
	
	if(!-d "$output/variants/snpsIndels/mutect"){
	    mkdir("$output/variants/snpsIndels/mutect", 0775) or die "Can't make $output/variants/snpsIndels/mutect";
	}
	my $ran_mutect = 0;
	my $mutectj = '';
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done" || $ran_ssf){  
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT", job_hold => "$ssfj", cpu => "2", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA7_MUTECT/java -Xmx16g -Djava.io.tmpdir=$tempdir -jar $MUTECT/muTect.jar --analysis_type MuTect --reference_sequence $REF_SEQ --dbsnp $DB_SNP --input_file:normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --input_file:tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --vcf $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf --out $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.txt -rf BadCigar --enable_extended_output --downsampling_type NONE`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done`;
	    $mutectj = "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT";
	    push @mu_jids, "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT";
	    push @all_jids, $mutectj;
	    $ran_mutect = 1;
	    $ran_mutect_glob = 1;
	}

	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done" || $ran_mutect){
	### sleep(2);
	###&generateMaf("$output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf", 'mutect', "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT", $data[0], $data[1]);
	###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done`;
	###}
	
	if($somaticsniper){
	    if(!-d "$output/variants/snpsIndels/somaticsniper"){
		mkdir("$output/variants/snpsIndels/somaticsniper", 0775) or die "Can't make $output/variants/snpsIndels/somaticsniper";
	    }
	    my $ran_somatic_sniper = 0;
	    my $ssj = '';
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.done" || $ran_ssf){  
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER", job_hold => "$ssfj", cpu => "2", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $SOMATIC_SNIPER/bam-somaticsniper -F vcf -f $REF_SEQ -q 1 $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER.done`;
		$ssj = "$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER";
		push @all_jids, $ssj;
		$ran_somatic_sniper = 1;
	    }

	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.done" || $ran_somatic_sniper){
	    ### `/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf -c somaticsniper -o $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.done`;
	    ###}
	}

	if($virmid){
	    my $ran_virmid = 0;
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done" || $ran_ssf){  
		sleep(2);
		if(-d "$output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid"){
		    `/bin/rm -rf $output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid`;
		}
		###mkdir("$output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid", 0775) or die "Can't make $output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid";
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_VIRMID", job_hold => "$ssfj", cpu => "4", mem => "12", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx12g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $VIRMID/Virmid.jar -R $REF_SEQ -D $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam -N $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam -t 4 -o $pre\_$data[0]\_$data[1]\_virmid -w $output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done`;
		push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_VIRMID";
		$ran_virmid = 1;
	    }
	}
	
	if($strelka){
	    my $ran_strelka_config = 0;
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done" || $ran_ssf){  
		if(-d "$output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka"){
		    ### STRELKA DIES IF DIR ALREADY EXISTS
		    `/bin/rm -rf $output/variants/snpsIndels/strelka//$data[0]\_$data[1]\_strelka`;
		}
		
		my @lns_jids = ();
		### NOTE: Strelka only recognizes X.bam.bai as the index for X.bam, not X.bai
		if((!-e "$output/progress/$pre\_indelRealigned_recal\_$data[0]\_LNS.done" || $ran_ssf)  && !-e "$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam.bai" && !$submitted_lns{$data[0]}){
		    sleep(2);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_indelRealigned_recal\_$data[0]\_LNS", job_hold => "$ssfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_indelRealigned_recal\_$data[0]\_LNS.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/ln -s $pre\_indelRealigned_recal\_$data[0]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam.bai`;
		    #push @lns_jids, "$pre\_indelRealigned_recal\_$data[0]\_LNS";
		    $submitted_lns{$data[0]} = 1;
		    `/bin/touch $output/progress/$pre\_indelRealigned_recal\_$data[0]\_LNS.done`;
		}
		
		if((!-e "$output/progress/$pre\_indelRealigned_recal\_$data[1]\_LNS.done" || $ran_ssf) && !-e "$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam.bai" && !$submitted_lns{$data[1]}){
		    sleep(2);
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_indelRealigned_recal\_$data[1]\_LNS", job_hold => "$ssfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_indelRealigned_recal\_$data[1]\_LNS.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/ln -s $pre\_indelRealigned_recal\_$data[1]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam.bai`;
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
		sleep(2);
		### NOTE: strelka ONLY HAS CONFIG FOR BWA ALN, NOT SURE HOW IT WILL WORK WITH BWA MEM
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", job_hold => "$lnsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.log");
		my $standardParams = Schedule::queuing(%stdParams);	    
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $STRELKA/bin/configureStrelkaWorkflow.pl --normal=$output/alignments//$pre\_indelRealigned_recal\_$data[0]\.bam --tumor=$output/alignments//$pre\_indelRealigned_recal\_$data[1]\.bam --ref=$REF_SEQ --config=$STRELKA/etc/strelka_config_bwa_default.ini --output-dir=$output/variants/snpsIndels/strelka//$data[0]\_$data[1]\_strelka`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done`;
		$ran_strelka_config = 1;
	    }

	    my $ran_strelka_run = 0;
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done" || $ran_strelka_config){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", cpu => "8", mem => "16", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /usr/bin/make -C $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka -j 8`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done`;
		$ran_strelka_run = 1;
	    }

	    if($ran_strelka_run){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/rm -rf $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/config $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/chromosomes $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/Makefile $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/task.complete`;
		push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP";
	    }
	    
	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.done" || $ran_strelka_run){
	    ### my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.log");
	    ###my $standardParams = Schedule::queuing(%stdParams);
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.snvs.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
	    
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.snvs.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
	    
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.indels.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
	    
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.indels.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.done`;
	    ###	}
	}
	
	if($scalpel){
	    my $ran_scalpel = 0;
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done" || $ran_ssf){
		sleep(2);
                if(!-d "$output/variants/snpsIndels/scalpel/"){
                    mkdir("$output/variants/snpsIndels/scalpel/", 0775) or die "Can't make $output/variants/snpsIndels/scalpel/";
                }
		if(!-d "$output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel"){
		    mkdir("$output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel", 0775) or die "Can't make $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel";
		}
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL", job_hold => "$ssfj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $SCALPEL/scalpel --somatic --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --bed $targets_bed_padded --ref $REF_SEQ --dir $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel --numprocs 24`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done`;
		$ran_scalpel = 1;
	    }
	
	    if($ran_scalpel){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL_CLEANUP", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL_CLEANUP.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /bin/rm -rf $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel/main/ $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel/validation/`;
		push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL_CLEANUP";
	    }
	}
	
	if($lancet){
	    if(!-d "$output/variants/snpsIndels/lancet"){
		mkdir("$output/variants/snpsIndels/lancet", 0775) or die "Can't make $output/variants/snpsIndels/lancet";
	    }

	    my @filterLancetVariants = ();
	    my @filter_lancet_jids = ();
	    my $ran_fl = 0;
	    foreach my $c (@ref_chr){
		my $ran_lancet = 0;
		my $lancet_jid = '';
		if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_LANCET.done" || $ran_ssf){
		    sleep(2);
		    
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_$c\_LANCET", job_hold => "$ssfj", cpu => "12", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_LANCET.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $LANCET/lancet --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --ref $REF_SEQ --reg $c --num-threads 12 ">$output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls.vcf"`;
		    
		    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_LANCET.done`;
		    $lancet_jid = "$pre\_$uID\_$data[0]\_$data[1]\_$c\_LANCET";
		    $ran_lancet = 1;
		}
	 
		if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET.done" || $ran_ssf){
		    sleep(2);
		    
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET", job_hold => "$lancet_jid", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $Bin/vcfFilterLancet.pl -input $output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls.vcf -output $output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls_FILTERED.vcf -log $output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls_DISCARDED.txt`;
		    
		    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET.done`;
		    push @filter_lancet_jids, "$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET";
		    $ran_fl = 1
		}
		
		push @filterLancetVariants, "--variant $output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls_FILTERED.vcf";
	    }
	    
	    
	    my $filterLancetVars = join(" " , @filterLancetVariants);
	    my $flj = join(",", @filter_lancet_jids);
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_CV_LANCET.done" || $ran_fl){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_CV_LANCET", job_hold => "$flj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_CV_LANCET.log");
		my $standardParams = Schedule::queuing(%stdParams);

		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/snpsIndels/lancet/$pre\_$data[0]\_$data[1]\_lancet_calls.vcf --assumeIdenticalSamples $filterLancetVars`;

		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_CV_LANCET.done`;
		push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_CV_LANCET";
	    }
	}
	
	if($vardict){

	    if(!-d "$output/variants/snpsIndels/vardict"){
		mkdir("$output/variants/snpsIndels/vardict", 0775) or die "Can't make $output/variants/snpsIndels/vardict";
	    }
	    
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT.done" || $ran_ssf){
		
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_VARDICT", job_hold => "$ssfj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT.log");
		my $standardParams = Schedule::queuing(%stdParams);
		my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "3");
		my $additionalParams = Schedule::additionalParams(%addParams);
		
		open(FILE, ">$output/intFiles/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT_command.sh") or die "Can't write to file $output/intFiles/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT_command.sh $!";
		
		print FILE "#! /bin/bash\n\n";
		print FILE "/opt/common/CentOS_6/VarDictJava/VarDict-1.5.1/bin/VarDict -G $REF_SEQ -f 0.01 -N $data[1] -b \"$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam|$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam\" -z -F 0 -c 1 -S 2 -E 3 -g 4 $targets_bed_padded | /opt/common/CentOS_6/VarDict/VarDict_85cc3f6/testsomatic.R | /opt/common/CentOS_6/VarDict/VarDict_85cc3f6/var2vcf_paired.pl -N \"$data[1]|$data[0]\" -f 0.01";
		close FILE;
		chmod 0750, "$output/intFiles/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT_command.sh";
		
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $output/intFiles/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT_command.sh ">$output/variants/snpsIndels/vardict/$pre\_$data[0]\_$data[1]\_vardict_calls.vcf"`;
		
		
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT.done`;
		push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_VARDICT";
	    }
	}


        if($pindel){
            my $ran_pindel = 0;
            if(!-d "$output/variants/snpsIndels/pindel"){
                mkdir("$output/variants/snpsIndels/pindel", 0775) or die "Can't make $output/variants/snpsIndels/pindel";
            }
            if(!-d "$output/intFiles/pindel"){
                mkdir("$output/intFiles/pindel", 0775) or die "Can't make $output/intFiles/pindel";
            }
            if(!-d "$output/intFiles/pindel/$data[0]\_$data[1]"){
                mkdir("$output/intFiles/pindel/$data[0]\_$data[1]", 0775) or die "Can't make $output/intFiles/pindel/$data[0]\_$data[1]";
            }
            if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL.done" || $ran_ssf){
                sleep(2);

                open(PC, ">$output/intFiles/$data[0]\_$data[1]\_pindel_config.txt") or die "Can't write to file $output/intFiles/$data[0]\_$data[1]\_pindel_config.txt $!";
                print PC "$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam\t200\t$data[1]\n";
                print PC "$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam\t200\t$data[0]\n";
                close PC;

                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_PINDEL", job_hold => "$ssfj", cpu => "12", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL.log");
                my $standardParams = Schedule::queuing(%stdParams);
                my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "3");
                my $additionalParams = Schedule::additionalParams(%addParams);

                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PINDEL/pindel -i $output/intFiles/$data[0]\_$data[1]\_pindel_config.txt -f $REF_SEQ -c ALL -o $output/intFiles/pindel/$data[0]\_$data[1]/$data[0]\_$data[1] -r false -t false -I false -T 12`;

                `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL.done`;
                $ran_pindel = 1;

            }
            if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF.done" || $ran_pindel){
                sleep(2);
                    
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_PINDEL", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF.log");
                my $standardParams = Schedule::queuing(%stdParams);

                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PINDEL/pindel2vcf --pindel_output_root $output/intFiles/pindel/$data[0]\_$data[1]/$data[0]\_$data[1] --reference $REF_SEQ --reference_name $species --reference_date $cur_date --vcf $output/variants/snpsIndels/pindel/$pre\_$data[0]\_$data[1]\_pindel_calls.vcf -b true --gatk_compatible`;

                `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF.done`;
                push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF";
            }
        }


	if($varscan){
	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done" || $ran_varscan_somatic){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp_MAF.txt -n $data[0] -t $data[1]`;
	}
    }
    close PAIR;
    
    
    if($hasPair && (!-e "$output/progress/$pre\_$uID\_HAPLOTECT.done" || $ran_mutect_glob || $ran_hc)){
	sleep(2);
	my $patientFile = "";
	if($patient){
	    $patientFile = "-patient $patient";
	}
	my $muj = join(",", @mu_jids);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HAPLOTECT", job_hold => "$hcj,$muj", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_HAPLOTECT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/haploTect_merge.pl -pair $pair -hc_vcf $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf -species $species -pre $pre -output $output/variants/snpsIndels/haplotect -mutect_dir $output/variants/snpsIndels/mutect -config $config $patientFile -align_dir $output/alignments/ -tempdir $tempdir -delete_temp`;
	
	push @all_jids, "$pre\_$uID\_HAPLOTECT";
	$haplotect_run = 1;
	`/bin/touch $output/progress/$pre\_$uID\_HAPLOTECT.done`;
	#$facets_haplotect_jid .= ",$pre\_$uID\_HAPLOTECT";
    }
}

my $allj2 = join(",", @all_jids);
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC_2", job_hold => "$allj2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_2.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' $curDir $rsync`;
`/bin/touch $output/progress/$pre\_$uID\_RSYNC_2.done`;

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

    if($pair){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_MAF_PAIRED", job_hold => "$hold", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_$jna\_MAF_PAIRED.log");
	my $standardParams = Schedule::queuing(%stdParams);
	###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $Bin/generateMAF.pl -vcf $vcf -pairing $pair -species $species -config $config -caller $type $n_sample $t_sample -delete_temp`;

	if($type =~ /unifiedgenotyper|ug|haplotypecaller|hc/i){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_MAF_UNPAIRED", job_hold => "$hold", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_$jna\_MAF_UNPAIRED.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $Bin/generateMAF.pl -vcf $vcf -species $species -config $config -caller $type -delete_temp`;
	}
    }
    else{
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_MAF_UNPAIRED", job_hold => "$hold", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_$jna\_MAF_UNPAIRED.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $Bin/generateMAF.pl -vcf $vcf -species $species -config $config -caller $type -delete_temp`;
    }

    push @all_jids, "$pre\_$uID\_$jna\_MAF_UNPAIRED";
}
