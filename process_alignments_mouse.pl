#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

use lib "$Bin/lib";
use Schedule;
use Cluster;
use File::Basename;
use POSIX qw(strftime);
my $cur_date = strftime "%Y%m%d", localtime;


my ($pair, $svnRev, $email, $patient, $group, $bamgroup, $config, $nosnps, $targets, $ug, $scheduler, $priority_project, $priority_group, $abra, $indelrealigner, $help, $step1, $allSomatic, $scalpel, $somaticsniper, $strelka, $varscan, $virmid, $lancet, $vardict, $pindel);

my $pre = 'TEMP';
my $output = "results";
my $species = 'mm10';

my $uID = `/usr/bin/id -u -n`;
chomp $uID;
my $rsync = "/ifs/solres/$uID";
my $tempdir = "/scratch/$uID";

GetOptions ('email=s' => \$email,
        'pre=s' => \$pre,
	    'pair=s' => \$pair,
	    'patient=s' => \$patient,
            'group=s' => \$group,
	    'config=s' => \$config,
	    'targets=s' => \$targets,
	    'species=s' => \$species,
	    'nosnps' => \$nosnps,
	    'step1' => \$step1,
	    'ug|unifiedgenotyper' => \$ug,
	    'abra' => \$abra,
	    'indelrealigner' => \$indelrealigner,
	    'bamgroup=s' => \$bamgroup,
            'svnRev=s' => \$svnRev,
 	    'scheduler=s' => \$scheduler,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
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
 	    'output|out|o=s' => \$output) or exit(1);


if(!$group || !$config || !$scheduler || !$targets || !$bamgroup || $help){
    print <<HELP;

    USAGE: process_alignments_hg19.pl -group GROUP -pair PAIR -pre PRE -config CONFIG -species SPECIES -scheduler SCHEDULER -targets TARGETS
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED)
	* BAMGROUP: files listing bams to be processed together; every bam for each group on 1 line, comma-separated (required)
	* TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS unless given full path to targets directory (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
	* PRE: output prefix (default: TEMP)
	* SPECIES: mm10 (default), mm10_custom, and mm9
	* OUTPUT: output results directory (default: results)
	* RSYNC:  path to rsync data for archive (default: /ifs/solres/USER_ID)
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
my $BCFTOOLS = '';
my $GATK = '';
my $FACETS_LIB = '';
my $FACETS_SUITE = '';
my $LANCET = '';
my $MUTECT = '';
my $PICARD = '';
my $PINDEL = '';
my $SAMTOOLS = '';
my $SCALPEL = '';
my $SOMATIC_SNIPER = '';
my $STRELKA = '';
my $TABIX = '';
my $VARDICT_JAVA = '';
my $VARDICT_PERL = '';
my $VARSCAN = '';
my $VIRMID = '';

my $JAVA = '';
my $JAVA7_MUTECT = '';
my $PYTHON = '';
my $PERL = '';

my $MM9_FASTA = '';
my $MM9_FAI = '';
my $MM9_BWA_INDEX = '';
my $MM10_FASTA = '';
my $MM10_FAI = '';
my $MM10_BWA_INDEX = '';
my $MM10_CUSTOM_FASTA = '';
my $MM10_CUSTOM_FAI = '';
my $MM10_CUSTOM_BWA_INDEX = '';

open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /abra/i){
	if(!-e "$conf[1]/abra.jar"){
	    die "CAN'T FIND GenomeAnalysisTK.jar IN $conf[1] $!";
	}
	$ABRA = $conf[1];
    }
    elsif($conf[0] =~ /bcftools/i){
        if(!-e "$conf[1]/bcftools"){
            die "CAN'T FIND bcftools IN $conf[1] $!";
        }
        $BCFTOOLS = $conf[1];
    }
    elsif($conf[0] =~ /facets_lib/i){
        if(!-e "$conf[1]/facets"){
            die "CAN'T FIND facets_lib IN $conf[1] $!";
        }
        $FACETS_LIB = $conf[1];
    }
    elsif($conf[0] =~ /facets_suite/i){
        if(!-e "$conf[1]/facets"){
            die "CAN'T FIND facets_suite IN $conf[1] $!";
        }
        $FACETS_SUITE = $conf[1];
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
    elsif($conf[0] =~ /tabix/i){
        if(!-e "$conf[1]/bgzip"){
            die "CAN'T FIND tabix IN $conf[1] $!";
        }
        $TABIX = $conf[1];
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
    }
     elsif($conf[0] =~ /^r$/i){
	if(!-e "$conf[1]/R"){
	    die "CAN'T FIND R IN $conf[1] $!";
	}
	my $path_tmp = $ENV{'PATH'};
	$ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
   elsif($conf[0] =~ /mm9_fasta/i){
	if(!-e "$conf[1]"){
	    if($species =~ /mm9/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$MM9_FASTA = $conf[1];
    }
   elsif($conf[0] =~ /mm9_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /mm9/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$MM9_FAI = $conf[1];
    }
    elsif($conf[0] =~ /mm9_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /mm9/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM9 WITH PREFIX $conf[1] $!";
	    }
	}
	$MM9_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /mm10_fasta/i){
	if(!-e "$conf[1]"){
	    if($species =~ /^mm10$/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$MM10_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /mm10_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /^mm10$/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$MM10_FAI = $conf[1];
    }
    elsif($conf[0] =~ /mm10_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /^mm10$/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM10 WITH PREFIX $conf[1] $!";
	    }
	}
	$MM10_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /mm10_custom_fasta/i){
	if(!-e "$conf[1]"){
	    if($species =~ /mm10_custom/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$MM10_CUSTOM_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /mm10_custom_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /mm10_custom/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$MM10_CUSTOM_FAI = $conf[1];
    }
    elsif($conf[0] =~ /mm10_custom_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /mm10_custom/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM10_CUSTOM WITH PREFIX $conf[1] $!";
	    }
	}
	$MM10_CUSTOM_BWA_INDEX = $conf[1];
    }
}
close CONFIG;

my $REF_SEQ = "$MM10_FASTA";
my $REF_FAI = "$MM10_FAI";
my $BWA_INDEX = "$MM10_BWA_INDEX";
my $ExAC_VCF;
my $DB_SNP = "$Bin/data/mm10/mm10_snp142.vcf";
my $CHR_M = 'M'; #as of right now default to M. Once we add NCBI assemblies, we can chang it below.
my $CHR_PREFIX = "chr"; #as of right now default to 'chr'. Once we add NCBI assemblies, we can change it below.

if($species =~ /mm10_custom/i){
    $REF_SEQ = "$MM10_CUSTOM_FASTA";
    $REF_FAI = "$MM10_CUSTOM_FAI";
    $BWA_INDEX = "$MM10_CUSTOM_BWA_INDEX";
}
elsif($species =~ /mm9/i){
    $REF_SEQ = "$MM9_FASTA";
    $REF_FAI = "$MM9_FAI";
    $BWA_INDEX = "$MM9_BWA_INDEX";
    $DB_SNP = "";
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
if(-d $targets){
    my @path = split(/\//, $targets);
    my $assay = pop @path;
    $targets_bed_padded = "$targets/$assay\_targets_plus5bp.bed";
}

if(!-e "$targets_bed_padded"){
    die "CAN'T LOCATE $targets_bed_padded FOR $targets; REQUIRED FOR SCALPEL $!";
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
my %processedBams = ();
my @finalBams = ();
my %ran_pr_glob = 0;
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
    my @indelBams = ();
    my $ran_ir == 0;
    my @ir_jids = ();

    foreach my $c (@ref_chr){
	if($abra){
	    my @inBams = ();
	    my @outBams = ();
	    my $ran_chr_split_index = 0;
	    my @chri_jids = ();
	    foreach my $pin (@pins){
		my @inB = split(/\s+/, $pin);
		my @bpath = split(/\//, $inB[1]);
		my $boi = pop @bpath;
		
		my $ran_chr_split = 0;
		my $chrsj = '';
		if(!-e "$output/progress/$pre\_$uID\_$boi\_$c\.done" || $step1){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$boi\_$c", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$boi\_$c\.log");
		    my $standardParams = Schedule::queuing(%stdParams);	    
		    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "2");
		    my $additionalParams = Schedule::additionalParams(%addParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $SAMTOOLS/samtools view -b -o $inB[1]\_$c\.bam $inB[1] $c`;
		
		    $chrsj = "$pre\_$uID\_$boi\_$c";
		    `/bin/touch $output/progress/$pre\_$uID\_$boi\_$c\.done`;
		    $ran_chr_split = 1;
		}

		if(!-e "$output/progress/$pre\_$uID\_$boi\_$c\_INDEX.done" || $ran_chr_split){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$boi\_$c\_INDEX", job_hold => "$chrsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$boi\_$c\_INDEX.log");
		    my $standardParams = Schedule::queuing(%stdParams);	    
		    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "2");
		    my $additionalParams = Schedule::additionalParams(%addParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $SAMTOOLS/samtools index $inB[1]\_$c\.bam`;
		    
		    push @chri_jids, "$pre\_$uID\_$boi\_$c\_INDEX";
		    `/bin/touch $output/progress/$pre\_$uID\_$boi\_$c\_INDEX.done`;
		    $ran_chr_split_index = 1;
		}

		push @inBams, "$inB[1]\_$c\.bam";
		push @outBams, "$inB[1]\_$c\_ABRA.bam";
	    }
	
	    my $aiBams = join(",", @inBams);
	    my $aoBams = join(",", @outBams);
	    my $cij = join(",", @chri_jids);
	    my $ran_abra = 0;
	    my $abra_jid = '';
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA.done" || $ran_chr_split_index){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_ABRA", job_hold => "$cij", cpu => "12", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA.log");
		my $standardParams = Schedule::queuing(%stdParams);	    
		my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "4");
		my $additionalParams = Schedule::additionalParams(%addParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/abra_wrapper.pl -inBams $aiBams -outBams $aoBams -refSeq $REF_SEQ -bwaRef $BWA_INDEX -targets $targets_bed_padded -working $output/intFiles/abra_$gpair[0]\_$c -config $config -log $output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA_WRAPPER.log`;
		
		$abra_jid = "$pre\_$uID\_$gpair[0]\_$c\_ABRA";
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA.done`;
		$ran_abra = 1;
	    }
	    
	    my $ran_fm = 0;
	    my @fm_jids = ();
	    my @fm_bams = ();
	    my $bcount = 0;
	    foreach my $outBam (@outBams){
		$bcount++;
		if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE.done" || $ran_abra){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE", job_hold => "$abra_jid", cpu => "1", mem => "50", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "5");
		    my $additionalParams = Schedule::additionalParams(%addParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar FixMateInformation I=$outBam O=$outBam\_FM.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true`;
		    push @fm_jids, "$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE";
		    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE.done`;
		    $ran_fm = 1;
		}
		push @fm_bams, "I=$outBam\_FM.bam";
	    }
	    	    
	    my $fmj = join(",", @fm_jids);
	    my $fmb = join(" ", @fm_bams);
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE_$c\_FM.done" || $ran_fm){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_MERGE_$c\_FM", job_hold => "$fmj", cpu => "8", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE_$c\_FM.log");
		my $standardParams = Schedule::queuing(%stdParams);
		my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "10");
		my $additionalParams = Schedule::additionalParams(%addParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $fmb O=$output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_MERGE_$c\_FM.done`;
		push @ir_jids, "$pre\_$uID\_$gpair[0]\_MERGE_$c\_FM";
		$ran_ir = 1;
	    }
	}
	elsif($indelrealigner){
	    my $ran_rtc = 0;
	    my $rtc_jid = '';
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_RTC.done" || $step1){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_RTC", cpu => "10", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_RTC.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ -L $c $multipleTargets -S LENIENT --known $DB_SNP -nt 10 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigner.intervals $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_RTC.done`;
		$rtc_jid = "$pre\_$uID\_$gpair[0]\_$c\_RTC";
		$ran_rtc = 1;
	    }
	    
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.done" || $ran_rtc){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_IR", job_hold => "$rtc_jid", cpu => "1", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ -L $c $multipleTargets -S LENIENT --knownAlleles $DB_SNP --targetIntervals $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.done`;
		push @ir_jids, "$pre\_$uID\_$gpair[0]\_$c\_IR";
		$ran_ir = 1;
	    }	    
	}
	else{
	    die "no indel realinger chosen #!";
	}
	push @indelBams, "-I $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam";
    }
	
    my $irBams = join(" ", @indelBams);
    my $ran_br = 0;
    my $brj = '';
    my $irj = join(",", @ir_jids);
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_BR.done" || $ran_ir){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_BR", job_hold => "$irj", cpu => "12", mem => "40", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_BR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO -R $REF_SEQ -S LENIENT --knownSites $DB_SNP --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_$gpair[0]\_recal_data.grp $irBams`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_BR.done`;
	$brj = "$pre\_$uID\_$gpair[0]\_BR";
	$ran_br = 1;
    }
        
    my @indelRecalBams1 = ();
    my $ran_pr = 0;
    my @pr_jids = ();
    foreach my $c (@ref_chr){
	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_PR.done" || $ran_br){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_PR", job_hold => "$brj", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_PR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R $REF_SEQ -L $c $multipleTargets --emit_original_quals -BQSR $output/intFiles/$pre\_$gpair[0]\_recal_data.grp --num_cpu_threads_per_data_thread 6 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned_recal.bam -I $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam`;
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_PR.done`;
	    push @pr_jids, "$pre\_$uID\_$gpair[0]\_$c\_PR";
	    push @prg_jids, "$pre\_$uID\_$gpair[0]\_$c\_PR";
	    $ran_pr = 1;
	    $ran_pr_glob{"$c"} = 1;
	}
    
	push @indelRecalBams1, "I=$output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned_recal.bam";
	push @{$processedBams{"$c"}}, "-I $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned_recal.bam";
    }

    my $irBams1 = join(" ", @indelRecalBams1);    
    my $ran_m = 0;
    my $prj = join(",", @pr_jids);
    my @merge_jids = ();
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE_PR.done" || $ran_pr){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_MERGE_PR", job_hold => "$prj", cpu => "8", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE_PR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "10");
	my $additionalParams = Schedule::additionalParams(%addParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $irBams1 O=$output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_MERGE_PR.done`;
	push @merge_jids, "$pre\_$uID\_$gpair[0]\_MERGE_PR";
	$ran_m = 1;
    }
        
    my $mj = join(",", @merge_jids);
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.done" || $ran_m){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_SSF", job_hold => "$mj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T SplitSamFile -R $REF_SEQ -I $output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam --outputRoot $output/alignments/$pre\_indelRealigned_recal_`;
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
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MQ_METRICS_$samp", job_hold => "$ssfj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_MQ_METRICS_$samp\.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MeanQualityByCycle INPUT=$finalBam OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.txt CHART_OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.pdf REFERENCE_SEQUENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true TMP_DIR=$tempdir`;
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
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeMeanQualityHistograms.py $output '*_MeanQualityByCycle_*.txt' $output/metrics/$pre\_post_recal_MeanQualityByCycle.txt $output/metrics/$pre\_pre_recal_MeanQualityByCycle.txt`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_MQ.done`;
    push @all_jids, "$pre\_$uID\_MERGE_MQ";
    $ran_mmqm = 1;
}

my $allj = join(",", @all_jids);
if(!-e "$output/progress/$pre\_$uID\_RSYNC_1.done" || $ran_ssf || $ran_mmqm){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC_1", job_hold => "$allj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_1.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' --exclude 'variants' --exclude 'metrics' $curDir $rsync`;
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
my @iVariants = ();
my @sVariants = ();
my @hcVariants = ();
my $ran_hc = 0;
my $ran_ug_snp = 0;
my $ran_ug_indel = 0;
my @hc_jids = ();
my @ugs_jids = ();
my @ugi_jids = ();
my $prgj = join(",", @prg_jids);

foreach my $c (@ref_chr){
    my $irBams2 = join(" ", @{$processedBams{"$c"}});

    if($ug){
	if(!-e "$output/progress/$pre\_$uID\_$c\_UG_SNP.done" || $ran_pr_glob{"$c"}){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$c\_UG_SNP", job_hold => "$prgj", cpu => "12", mem => "48", cluster_out => "$output/progress/$pre\_$uID\_$c\_UG_SNP.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name $species -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model SNP --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_$c\_UnifiedGenotyper_SNP.vcf $irBams2`;
	    `/bin/touch $output/progress/$pre\_$uID\_$c\_UG_SNP.done`;
	    push @ugs_jids, "$pre\_$uID\_$c\_UG_SNP";
	    $ran_ug_snp = 1;
	}
	
	if(!-e "$output/progress/$pre\_$uID\_$c\_UG_INDEL.done" || $ran_pr_glob{"$c"}){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$c\_UG_INDEL", job_hold => "$prgj", cpu => "12", mem => "48", cluster_out => "$output/progress/$pre\_$uID\_$c\_UG_INDEL.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name $species -L $c $multipleTargets --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model INDEL --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_$c\_UnifiedGenotyper_INDEL.vcf $irBams2`;
	    `/bin/touch $output/progress/$pre\_$uID\_$c\_UG_INDEL.done`;
	    push @ugi_jids, "$pre\_$uID\_$c\_UG_INDEL";
	    $ran_ug_indel = 1;
	}

	push @sVariants, "--variant $output/intFiles/$pre\_$c\_UnifiedGenotyper_SNP.vcf";
	push @iVariants, "--variant $output/intFiles/$pre\_$c\_UnifiedGenotyper_INDEL.vcf";
    }
    
    if(!-e "$output/progress/$pre\_$uID\_$c\_HC.done" || $ran_pr_glob{"$c"}){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$c\_HC", job_hold => "$prgj", cpu => "30", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_$c\_HC.log");
	my $standardParams = Schedule::queuing(%stdParams);
	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "7");
	my $additionalParams = Schedule::additionalParams(%addParams);
	my $response = "";
    my $HC_submission_num = 0;
    while($response !~ /is submitted to default queue/ && $HC_submission_num < 10){
        my $time = 60 * $HC_submission_num;
        `sleep $time `;
        $HC_submission_num++;
        print "HC Chromosome $c attempt number $HC_submission_num\n"; 
        print "Command: $standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 30 --out $output/intFiles/$pre\_$c\_HaplotypeCaller.vcf $irBams2";
        $response = `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 30 --out $output/intFiles/$pre\_$c\_HaplotypeCaller.vcf $irBams2     2>&1`;
        print "Response: $response\n";
        if($response !~ /is submitted to default queue/ && $HC_submission_num == 10 && $c eq '1'){
            `mail -s "HaplotypeCaller Job Not Working" $email <<< "This is an e-mail letting you know that your variants pipeline run for $pre is not submitting HC jobs to the cluster. Please remove all .done files for haplotype caller and rerun the pipeline."  `;
        }
    }
    `/bin/touch $output/progress/$pre\_$uID\_$c\_HC.done`;
	push @hc_jids, "$pre\_$uID\_$c\_HC";
	$ran_hc = 1;
    }

    push @hcVariants, "--variant $output/intFiles/$pre\_$c\_HaplotypeCaller.vcf";
}

my $hcVars = join(" " , @hcVariants);
my $ran_cv_hc = 0;
my $hcj = join(",", @hc_jids);
my $cvhcj = '';
if(!-e "$output/progress/$pre\_$uID\_CV_HC.done" || $ran_hc){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_HC", job_hold => "$hcj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf --assumeIdenticalSamples $hcVars`;
    `/bin/touch $output/progress/$pre\_$uID\_CV_HC.done`;
    $cvhcj = "$pre\_$uID\_CV_HC";
    $ran_cv_hc = 1;
}

sleep(2);
# Run this anyway, it will go through the function and make sure everything was ran
&generateMaf("$output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'haplotypecaller', "$cvhcj,$ssfj", $ran_cv_hc);

if($ug){
    if(!-d "$output/variants/snpsIndels/unifiedgenotyper"){
	mkdir("$output/variants/snpsIndels/unifiedgenotyper", 0775) or die "Can't make $output/variants/snpsIndels/unifiedgenotyper";
    }
    my $sVars = join(" " , @sVariants);
    my $iVars = join(" " , @iVariants);

    my $ran_cv_ug_snp = 0;
    my @cvug_jids = ();
    my $ugsj = join(",", @ugs_jids);
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_SNP.done" || $ran_ug_snp){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_UG_SNP", job_hold => "$ugsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG_SNP.log");
	my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf --assumeIdenticalSamples $sVars`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_SNP.done`;
	push @cvug_jids, "$pre\_$uID\_CV_UG_SNP";
	$ran_cv_ug_snp = 1;
    }

    my $ran_cv_ug_indel = 0;
    my $ugij = join(",", @ugi_jids);
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_INDEL.done" || $ran_ug_indel){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_UG_INDEL", job_hold => "$ugij", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG_INDEL.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf --assumeIdenticalSamples $iVars`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_INDEL.done`;
	push @cvug_jids, "$pre\_$uID\_CV_UG_INDEL";
	$ran_cv_ug_indel = 1;
    }

    my $cvugj = join(",", @cvug_jids);
    if(!-e "$output/progress/$pre\_$uID\_CV_UG_RAW.done" || $ran_cv_ug_snp || $ran_cv_ug_indel){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_UG_RAW", job_hold => "$cvugj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG_RAW.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf --assumeIdenticalSamples --variant $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf --variant $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_UG_RAW.done`;
    }

    my $ran_vf_ug_snp = 0;
    my @vfug_jids = ();
    if(!-e "$output/progress/$pre\_$uID\_VF_UG_SNP.done" || $ran_ug_snp){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VF_UG_SNP", job_hold => "$ugsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_VF_UG_SNP.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --mask $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf --maskName nearIndel --variant $output/intFiles/$pre\_UnifiedGenotyper_SNP.vcf -o $output/intFiles/$pre\_UnifiedGenotyper_SNP_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"MQ \\< 40.0\\" --filterExpression \\"FS \\> 60.0\\" --filterExpression \\"HaplotypeScore \\> 13.0\\" --filterExpression \\"MQRankSum \\< -12.5\\" --filterExpression \\"ReadPosRankSum \\< -8.0\\" --filterName QDFilter --filterName MQFilter --filterName FSFilter --filterName HSFilter --filterName MQRSFilter --filterName ReadPosFilter`;
	`/bin/touch $output/progress/$pre\_$uID\_VF_UG_SNP.done`;
	push @vfug_jids, "$pre\_$uID\_VF_UG_SNP";
	$ran_vf_ug_snp = 1;
    }

    my $ran_vf_ug_indel = 0;
    if(!-e "$output/progress/$pre\_$uID\_VF_UG_INDEL.done" || $ran_ug_indel){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VF_UG_INDEL", job_hold => "$ugij", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_VF_UG_INDEL.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ --variant $output/intFiles/$pre\_UnifiedGenotyper_INDEL.vcf -o $output/intFiles/$pre\_UnifiedGenotyper_INDEL_vf.vcf --clusterWindowSize 10 --filterExpression \\"QD \\< 2.0\\" --filterExpression \\"ReadPosRankSum \\< -20.0\\" --filterExpression \\"InbreedingCoeff \\< -0.8\\" --filterExpression \\"FS \\> 200.0\\" --filterName QDFilter --filterName ReadPosFilter --filterName InbreedingFilter --filterName FSFilter`;
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx2g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/snpsIndels/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf --assumeIdenticalSamples --variant $output/intFiles/$pre\_UnifiedGenotyper_SNP_vf.vcf --variant $output/intFiles/$pre\_UnifiedGenotyper_INDEL_vf.vcf`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG_SI.done`;
        $cvugsij = "$pre\_$uID\_CV_UG_SI";
	$ran_cv_ug_si = 1;
    }
        
    if(!-e "$output/progress/$pre\_$uID\_MAF_UG.done" || $ran_cv_ug_si){  
	sleep(2);
	&generateMaf("$pre\_UnifiedGenotyper.vcf", 'unifiedgenotyper', "$cvugsij", $ran_cv_ug_si);
 	`/bin/touch $output/progress/$pre\_$uID\_MAF_UG.done`;
    }
}

my $hasPair = 0;

if($pair){
    if(!-d "$output/variants/snpsIndels/haplotect"){
	mkdir("$output/variants/snpsIndels/haplotect", 0775) or die "Can't make $output/variants/snpsIndels/haplotect";
    }
    if(!-d "$output/variants/copyNumber"){
	#mkdir("$output/variants/copyNumber", 0775) or die "Can't make $output/variants/copyNumber";
    }
    if(!-d "$output/variants/copyNumber/facets"){
	#mkdir("$output/variants/copyNumber/facets", 0775) or die "Can't make $output/variants/copyNumber/facets";
    }

    open(PAIR, "$pair") or die "Can't open $pair file";
    my %submitted_lns = ();
    my @mu_jids = ();
    my $ran_mutect_glob = 0;
    my $haplotect_run = 0;
    my $facets_run = 0;
    my @facets_jid = ();

       while(<PAIR>){
	chomp;
	
	my @data = split(/\s+/, $_);

	if($data[0] =~ /^NA$/i || $data[1] =~ /^NA$/i){
	    next;
	}
        ## This means there really is a sample pair, haplotect should run.
        $hasPair=1;

	if(!-d "$output/variants/snpsIndels/mutect"){
	    mkdir("$output/variants/snpsIndels/mutect", 0775) or die "Can't make $output/variants/snpsIndels/mutect";
	}
	my $ran_mutect = 0;
	my $mutectj = '';
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done" || $ran_ssf){  
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT", job_hold => "$ssfj", cpu => "2", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA7_MUTECT/java -Xmx16g -Djava.io.tmpdir=$tempdir -jar $MUTECT/muTect.jar --analysis_type MuTect --reference_sequence $REF_SEQ --dbsnp $DB_SNP --input_file:normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --input_file:tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --vcf $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf --out $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.txt -rf BadCigar --enable_extended_output --downsampling_type NONE`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $SOMATIC_SNIPER/bam-somaticsniper -F vcf -f $REF_SEQ -q 1 $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx12g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $VIRMID/Virmid.jar -R $REF_SEQ -D $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam -N $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam -t 4 -o $pre\_$data[0]\_$data[1]\_virmid -w $output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid`;
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
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/ln -s $pre\_indelRealigned_recal\_$data[0]\.bai $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam.bai`;
		    #push @lns_jids, "$pre\_indelRealigned_recal\_$data[0]\_LNS";
		    $submitted_lns{$data[0]} = 1;
		    `/bin/touch $output/progress/$pre\_indelRealigned_recal\_$data[0]\_LNS.done`;
		}
		
		if((!-e "$output/progress/$pre\_indelRealigned_recal\_$data[1]\_LNS.done" || $ran_ssf) && !-e "$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam.bai" && !$submitted_lns{$data[1]}){
		    sleep(2);
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
		sleep(2);
		### NOTE: strelka ONLY HAS CONFIG FOR BWA ALN, NOT SURE HOW IT WILL WORK WITH BWA MEM
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", job_hold => "$lnsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.log");
		my $standardParams = Schedule::queuing(%stdParams);	    
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STRELKA/bin/configureStrelkaWorkflow.pl --normal=$output/alignments//$pre\_indelRealigned_recal\_$data[0]\.bam --tumor=$output/alignments//$pre\_indelRealigned_recal\_$data[1]\.bam --ref=$REF_SEQ --config=$STRELKA/etc/strelka_config_bwa_default.ini --output-dir=$output/variants/snpsIndels/strelka//$data[0]\_$data[1]\_strelka`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done`;
		$ran_strelka_config = 1;
	    }

	    my $ran_strelka_run = 0;
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done" || $ran_strelka_config){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", cpu => "8", mem => "16", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /usr/bin/make -C $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka -j 8`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN.done`;
		$ran_strelka_run = 1;
	    }

	    if($ran_strelka_run){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/rm -rf $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/config $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/chromosomes $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/Makefile $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/task.complete`;
		push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CLEANUP";
	    }
	    
	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.done" || $ran_strelka_run){
	    ### my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_RUN_MAF.log");
	    ###my $standardParams = Schedule::queuing(%stdParams);
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.snvs.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
	    
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.snvs.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.snvs_MAF.txt -n $data[0] -t $data[1]`;
	    
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/all.somatic.indels.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_all.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
	    
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/passed.somatic.indels.vcf -c strelka -o $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka/results/$pre\_$data[0]\_$data[1]\_STRELKA_passed.somatic.indels_MAF.txt -n $data[0] -t $data[1]`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $SCALPEL/scalpel --somatic --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --bed $targets_bed_padded --ref $REF_SEQ --dir $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel --numprocs 24`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL.done`;
		$ran_scalpel = 1;
	    }
	    
	    if($ran_scalpel){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL_CLEANUP", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SCALPEL_CLEANUP.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /bin/rm -rf $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel/main/ $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel/validation/`;
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
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $LANCET/lancet --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --ref $REF_SEQ --reg $c --num-threads 12 ">$output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls.vcf"`;
		    
		    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_LANCET.done`;
		    $lancet_jid = "$pre\_$uID\_$data[0]\_$data[1]\_$c\_LANCET";
		    $ran_lancet = 1;
		}
	 
		if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET.done" || $ran_ssf){
		    sleep(2);
		    
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET", job_hold => "$lancet_jid", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_$c\_FILTER_LANCET.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/vcfFilterLancet.pl -input $output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls.vcf -output $output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls_FILTERED.vcf -log $output/intFiles/$pre\_$data[0]\_$data[1]\_$c\_lancet_calls_DISCARDED.txt`;
		    
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

		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/variants/snpsIndels/lancet/$pre\_$data[0]\_$data[1]\_lancet_calls.vcf --assumeIdenticalSamples $filterLancetVars`;

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
		
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $output/intFiles/$pre\_$uID\_$data[0]\_$data[1]\_VARDICT_command.sh ">$output/variants/snpsIndels/vardict/$pre\_$data[0]\_$data[1]\_vardict_calls.vcf"`;
				
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
                my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "3");+                my $additionalParams = Schedule::additionalParams(%addParams);
                my $additionalParams = Schedule::additionalParams(%addParams);

                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PINDEL/pindel -i $output/intFiles/$data[0]\_$data[1]\_pindel_config.txt -f $REF_SEQ -c ALL -o $output/intFiles/pindel/$data[0]\_$data[1]/$data[0]\_$data[1] -r false -t false -I false -T 12`;

                `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL.done`;
                $ran_pindel = 1;

            }
            if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF.done" || $ran_pindel){
                sleep(2);
                    
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_PINDEL", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF.log");
                my $standardParams = Schedule::queuing(%stdParams);

                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PINDEL/pindel2vcf --pindel_output_root $output/intFiles/pindel/$data[0]\_$data[1]/$data[0]\_$data[1] --reference $REF_SEQ --reference_name $species --reference_date $cur_date --vcf $output/variants/snpsIndels/pindel/$pre\_$data[0]\_$data[1]\_pindel_calls.vcf -b true --gatk_compatible`;

                `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF.done`;
                push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_PINDEL2VCF";
            }
        }


	if($varscan){
	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done" || $ran_varscan_somatic){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp_MAF.txt -n $data[0] -t $data[1]`;
	}


=begin FOR_FACETS
        ## Here we will add the facets scripts
        ## These are the #'s Nick uses
        my $MINCOV=0;
        my $BASEQ=20;
        my $MAPQ=15;

        ## Set up tumor and normal counts
	if(!-d "$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets"){
	mkdir("$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets", 0775) or die "Can't make $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets";
}
	if(!-d "$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp"){
	mkdir("$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp", 0775) or die "Can't make $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp";
}
        my $facetsSETUP_jid = '';
        my $facets_setup = 0;
        if($hasPair && (! -e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.done" || $ssfj )) {
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_facets_SETUP",  cpu => "4", mem => "5", job_hold => "$ssfj", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_facets_SETUP.log");
            my $standardParams = Schedule::queuing(%stdParams);
            my %addParams = (runtime => "30");
            my $additionalParams = Schedule::additionalParams(%addParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $Bin/facets/bin/GetBaseCounts --filter_improper_pair --sort_output --fasta $REF_SEQ --vcf $DB_SNP --maq $MAPQ --baq $BASEQ --cov $MINCOV --bam $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --out $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_indelRealigned_recal\_$data[0].dat`;

            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[1]\_facets_SETUP",  cpu => "4", mem => "5", job_hold => "$ssfj", cluster_out => "$output/progress/$pre\_$uID\_$data[1]\_facets_SETUP.log");
            my $standardParams2 = Schedule::queuing(%stdParams);
            %addParams = (runtime => "30");
            my $additionalParams2 = Schedule::additionalParams(%addParams);
            `$standardParams2->{submit} $standardParams2->{job_name} $standardParams2->{cpu} $standardParams2->{mem} $standardParams2->{job_hold} $standardParams2->{cluster_out} $additionalParams2 $Bin/facets/bin/GetBaseCounts --filter_improper_pair --sort_output --fasta $REF_SEQ --vcf $DB_SNP --maq $MAPQ --baq $BASEQ --cov $MINCOV --bam $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --out $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_indelRealigned_recal\_$data[1].dat`;


            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_merge_counts_facets_SETUP",  cpu => "4", mem => "18", job_hold => "$pre\_$uID\_$data[0]\_facets_SETUP,$pre\_$uID\_$data[1]\_facets_SETUP", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.log");
            my $standardParams3 = Schedule::queuing(%stdParams);
            %addParams = (runtime => "30");
            my $additionalParams3 = Schedule::additionalParams(%addParams);
            `$standardParams3->{submit} $standardParams3->{job_name} $standardParams3->{cpu} $standardParams3->{mem} $standardParams3->{job_hold} $standardParams3->{cluster_out} $additionalParams3 $Bin/facets/mergeTN.R  $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_indelRealigned_recal\_$data[1].dat  $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_indelRealigned_recal\_$data[0].dat $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_countsMerged_$data[0]\_$data[1].dat.gz`;

            $facets_setup = 1;
            `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.done`;
            $facetsSETUP_jid = "$pre\_$uID\_$data[0]\_$data[1]\_merge_counts_facets_SETUP";
        }
        ## now facets
        if($hasPair && (! -e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.done" || $facets_setup)){

            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN",  cpu => "3", mem => "2", job_hold => "$facetsSETUP_jid", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.log");
            my $standardParams = Schedule::queuing(%stdParams);
            my %addParams = (runtime => "10");
            my $additionalParams = Schedule::additionalParams(%addParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $Bin/facets/facets_RUN.sh $FACETS_SUITE $FACETS_LIB $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets $data[0]\_$data[1] $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_countsMerged_$data[0]\_$data[1].dat $species 300 100`;
            push @facets_jid, "$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN" ;
            $facets_run = 1;
            `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.done`;
        }
        `/bin/echo "$data[1]\t$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/$data[0]\_$data[1]\_hisens.Rdata" >> $output/variants/copyNumber/facets/facets_mapping.txt`;
=end FOR_FACETS
=cut
    
    }
    close PAIR;
    my $facets_haplotect_jid = '';
=begin FOR_FACETS
    if($hasPair && (! -e "$output/progress/$pre\_$uID\_merge_facets_seg.done" || $facets_run)){
        my $seg_outfile = "$output/variants/copyNumber/facets/$pre\_facets_merge_hisens.seg";
        if( -f "$seg_outfile"){
            unlink("$seg_outfile") or die "Cannot delete? $!";
        }
        my $facets_js = join(",", @facets_jid);

        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_merge_facets_seg",  cpu => "1", mem => "1", job_hold => "$facets_js", cluster_out => "$output/progress/$pre\_$uID\_merge_facets_seg.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my %addParams = (scheduler => "$scheduler", runtime => "1", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "0");
        my $additionalParams = Schedule::additionalParams(%addParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/facets/merge_facets_seg.pl -facets_dir $output/variants/copyNumber/facets -outfile $seg_outfile`;

        `/bin/touch $output/progress/$pre\_$uID\_merge_facets_seg.done`;
        $facets_haplotect_jid = "$pre\_$uID\_merge_facets_seg";
    }
=end FOR_FACETS
=cut

    if($hasPair && (!-e "$output/progress/$pre\_$uID\_HAPLOTECT.done" || $ran_mutect_glob || $ran_hc)){
        sleep(2);
        my $patientFile = "";
        my $addOptions = "";
        if($ExAC_VCF){
            $addOptions = "-exac_vcf $ExAC_VCF";
        }
        if($patient){
            $patientFile = "-patient $patient";
        }
        my $muj = join(",", @mu_jids);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HAPLOTECT", job_hold => "$hcj,$muj", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_HAPLOTECT.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/haploTect_merge.pl -pair $pair -hc_vcf $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf -species $species -pre $pre -output $output/variants/snpsIndels/haplotect -mutect_dir $output/variants/snpsIndels/mutect -config $config $patientFile -align_dir $output/alignments/ -svnRev $svnRev $addOptions -delete_temp`;

        $haplotect_run = 1;
        `/bin/touch $output/progress/$pre\_$uID\_HAPLOTECT.done`;
        $facets_haplotect_jid .= ",$pre\_$uID\_HAPLOTECT";
         push @all_jids, "$pre\_$uID\_HAPLOTECT";
    }

=begin FOR_FACETS

    if($hasPair && (!-e "$output/progress/$pre\_$uID\_join_maf.done" || $haplotect_run || $facets_run)){
        sleep(2);

        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_join_maf", job_hold => "$facets_haplotect_jid", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_join_maf.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $FACETS_SUITE/facets mafAnno -m $output/variants/snpsIndels/haplotect/$pre\_haplotect_VEP_MAF.txt -f $output/variants/copyNumber/facets/facets_mapping.txt -o $output/variants/$pre\_CMO_MAF.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_join_maf.done`;
        push @all_jids, "$pre\_$uID\_join_maf";
    }

=end FOR_FACETS
=cut

}


my $allj2 = join(",", @all_jids);
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC_2", job_hold => "$allj2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_2.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' $curDir $rsync`;
`/bin/touch $output/progress/$pre\_$uID\_RSYNC_2.done`;


sub generateMaf{
    my ($vcf, $type, $hold, $ran_hc, $normal_sample, $tumor_sample) = @_;
    
    my $n_sample = '';
    my $t_sample = '';
    if($normal_sample && $tumor_sample){
        $n_sample = "-normal_sample $normal_sample";
        $t_sample = "-tumor_sample $tumor_sample";
    }
    
    my $vcf_dir = dirname($vcf);
    my $jna = basename($vcf);
    $jna =~ s/\//_/g;
    
    my $patientFile = "";
    if($patient){
        $patientFile = "-patient $patient";
    }
    
    my $bgz_jid = '';
    my $bgzipped = 0;
    
    my @vcf_files;
    my @chr_maf_jids;
    # split and send each split thing to generate maf separately
    foreach my $c (@ref_chr){
        if( !-e "$output/progress/$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED.done" || $ran_hc ){
	    
            if((! -e "$vcf.gz" || $ran_hc) && !$bgzipped){
                $bgzipped = 1;
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_bgzip", job_hold => "$hold", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_bgzip.log");
                my $standardParams = Schedule::queuing(%stdParams);
                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$TABIX/bgzip -cf $vcf > $vcf.gz"`;
		
		%stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_bgzip_index", job_hold => "$pre\_$uID\_$jna\_bgzip", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_bgzip_index.log");
		$standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $BCFTOOLS/bcftools index $vcf.gz`;
		
		$bgz_jid = "$pre\_$uID\_$jna\_bgzip,$pre\_$uID\_$jna\_bgzip_index";
            }
	    if(!-d "$vcf_dir/chrom_$c"){
                mkdir("$vcf_dir/chrom_$c", 0775) or die "Can't make $vcf_dir/chrom_$c";
            }
            my $addOptions = "";
            if($ExAC_VCF){
                $addOptions = "-exac_vcf $ExAC_VCF";
            }

            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_split__$c", job_hold => $bgz_jid, cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_split_$c.log");
            my $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $BCFTOOLS/bcftools filter -r $c $vcf.gz -O v -o $vcf_dir/chrom_$c/$jna\_$c.vcf`;
	    
            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED", job_hold => "$pre\_$uID\_$jna\_split__$c", cpu => "4", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED.log");
            $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/generateMAF.pl -vcf $vcf_dir/chrom_$c/$jna\_$c.vcf -species $species -config $config -caller $type $patientFile -align_dir $output/alignments $addOptions -delete_temp`;
          `/bin/touch $output/progress/$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED.done`;
            push @all_jids, "$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED";
            push @chr_maf_jids, "$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED";
            push @vcf_files, "$vcf_dir/chrom_$c/$jna\_$c.vcf";
        }
        #push @vcf_files, "$vcf_dir/chrom_$c/$jna\_$c.vcf";
        #push @chr_maf_jids, "$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED";
    }
    
    my $jid_holds = join(",", @chr_maf_jids);
    my @merge_files = @vcf_files;
    s/vcf$/vcf_UNPAIRED_TCGA_MAF.txt/g for @merge_files;
    my $merge_files = join(" -i ", @merge_files);
    my @merge_jids;
    
    if(@chr_maf_jids){
        # MAKE A JOB HERE THAT WILL MERGE ALL THE MAF FILES FOR TCGA_MAF, TCGA_PORTAL_MAF, TCGA_PORTAL_MAF_fillout, and VEP_MAF
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_TCGA_MAF", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_TCGA_MAF.log");
        my $standardParams = Schedule::queuing(%stdParams);
	###print "$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_MAF.txt\n";
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_MAF";
	
        $merge_files =~ s/UNPAIRED_TCGA_MAF.txt/UNPAIRED_TCGA_PORTAL_MAF.txt/g;
        %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF.log");
	$standardParams = Schedule::queuing(%stdParams);
	###print "$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF.txt\n";
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF";

        $merge_files =~ s/UNPAIRED_TCGA_PORTAL_MAF.txt/UNPAIRED_VEP_MAF.txt/g;
        %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_VEP_MAF", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_VEP_MAF.log");
        $standardParams = Schedule::queuing(%stdParams);
	###print "$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_VEP_MAF.txt\n";
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_VEP_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_VEP_MAF";
	
        if($patient){
            s/vcf_UNPAIRED_TCGA_MAF.txt$/vcf_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt/g for @merge_files;
            $merge_files = join(" -i ", @merge_files);
            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF_fillout", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF_fillout.log");
            $standardParams = Schedule::queuing(%stdParams);
            ###print "$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt\n";
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt`;
            push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF_fillout";
        }
	
        #Now to clean up!
        my $merge_holds = join(",", @merge_jids);
        push @all_jids, "$merge_holds";
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_cleanup", job_hold => "$merge_holds", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$jna\_cleanup.log");
        $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams rm -rf $vcf_dir/chrom_* $vcf.gz $vcf.gz.csi`;
    }
}
