#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Schedule;
use Cluster;
use File::Basename;

my ($patient, $svnRev, $pair, $group, $bamgroup, $config, $nosnps, $targets, $ug, $scheduler, $priority_project, $priority_group, $abra, $help, $step1, $allSomatic, $scalpel, $somaticsniper, $strelka, $varscan, $virmid);

my $pre = 'TEMP';
my $output = "results";
my $species = 'b37';

my $uID = `/usr/bin/id -u -n`;
chomp $uID;
my $rsync = "/ifs/solres/$uID";

GetOptions ('pre=s' => \$pre,
	    'pair=s' => \$pair,
	    'patient=s' => \$patient,
            'group=s' => \$group,
	    'config=s' => \$config,
	    'targets=s' => \$targets,
	    'species=s' => \$species,
	    'nosnps' => \$nosnps,
	    'ug|unifiedgenotyper' => \$ug,
	    'abra' => \$abra,
	    'step1' => \$step1,
	    'bamgroup=s' => \$bamgroup,
 	    'scheduler=s' => \$scheduler,
 	    'svnRev=s' => \$svnRev,
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
 	    'output|out|o=s' => \$output) or exit(1);


if(!$group || !$config || !$scheduler || !$targets || !$bamgroup || $help){
    print <<HELP;

    USAGE: process_alignments_human.pl -group GROUP -bampgroup BAMBROUP -config CONFIG -scheduler SCHEDULER -targets TARGETS
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED)
	* BAMGROUP: files listing bams to be processed together; every bam for each group on 1 line, comma-separated (required)
	* SPECIES: b37 (default), hg19, hybrid (b37_mm10)
	* TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS unless given full path to targets directory (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
	* PRE: output prefix (default: TEMP)
	* OUTPUT: output results directory (default: results)
	* RSYNC:  path to rsync data for archive (default: /ifs/solres/USER_ID)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* -nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed
	* -abra: run abra instead of GATK indelrealigner
	* -step1: forece the pipeline to start from the first step in pipeline
	* haplotypecaller is default; -ug || -unifiedgenotyper to also make unifiedgenotyper variant calls	
	* ALLSOMATIC: run all somatic callers; mutect/haplotypecaller always run; otherwise -scalpel, -somaticsniper, -strelka, -varscan, -virmid to run them individually	
HELP
exit;
}

if(!-d "$rsync"){
    die "Can't rsync to $rsync. Please pick a differnt location using the rsync parameter $!";
}

my $curDir = `pwd`;
chomp $curDir;
if($output !~ /^\//){
    $output = "$curDir/$output";
}

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
my $FACETS_LIB = '';
my $FACETS_SUITE = '';
my $BCFTOOLS= '';
my $TABIX = '';
my $PYTHON = '';
my $PERL = '';
my $B37_FASTA = '';
my $B37_FAI = '';
my $B37_MM10_HYBRID_FASTA = '';
my $B37_MM10_HYBRID_FAI = '';
my $B37_BWA_INDEX = '';
my $B37_MM10_HYBRID_BWA_INDEX = '';
my $HG19_FASTA = '';
my $HG19_FAI = '';
my $HG19_MM10_HYBRID_FASTA = '';
my $HG19_MM10_HYBRID_FAI = '';
my $HG19_BWA_INDEX = '';
my $HG19_MM10_HYBRID_BWA_INDEX = '';

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
    elsif($conf[0] =~ /tabix/i){
        if(!-e "$conf[1]/bgzip"){
            die "CAN'T FIND tabix IN $conf[1] $!";
        }
        $TABIX = $conf[1];
    }
    elsif($conf[0] =~ /bcftools/i){
        if(!-e "$conf[1]/bcftools"){
            die "CAN'T FIND bcftools IN $conf[1] $!";
        }
        $BCFTOOLS = $conf[1];
    }
    elsif($conf[0] =~ /facets_suite/i){
        if(!-e "$conf[1]/facets"){
            die "CAN'T FIND facets_suite IN $conf[1] $!";
        }
        $FACETS_SUITE = $conf[1];
    }
    elsif($conf[0] =~ /facets_lib/i){
        if(!-e "$conf[1]/facets"){
            die "CAN'T FIND facets_lib IN $conf[1] $!";
        }
        $FACETS_LIB = $conf[1];
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
	my $path_tmp = $ENV{'PATH'};
	$ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /perl/i){
	if(!-e "$conf[1]/perl"){
	    die "CAN'T FIND perl IN $conf[1] $!";
	}
	$PERL = $conf[1];
	my $path_tmp = $ENV{'PATH'};
	$ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /python/i){
	if(!-e "$conf[1]/python"){
	    die "CAN'T FIND python IN $conf[1] $!";
	}
	$PYTHON = $conf[1];
	my $path_tmp = $ENV{'PATH'};
	$ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /^r$/i){
	if(!-e "$conf[1]/R"){
	    die "CAN'T FIND R IN $conf[1] $!";
	}
	my $path_tmp = $ENV{'PATH'};
	$ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /b37_fasta/i){
	if(!-e "$conf[1]"){
	    if($species =~ /human|^b37$/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$B37_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /b37_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /human|^b37$/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$B37_FAI = $conf[1];
    }
    elsif($conf[0] =~ /b37_mm10_hybrid_fasta/i){
	if(!-e "$conf[1]"){
	    if($species =~ /hybrid|b37_mm10/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$B37_MM10_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /b37_mm10_hybrid_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /hybrid|b37_mm10/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$B37_MM10_HYBRID_FAI = $conf[1];
    }
    elsif($conf[0] =~ /b37_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /^b37$/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR B37 WITH PREFIX $conf[1] $!";
	    }
	}
	$B37_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /b37_mm10_hybrid_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /hybrid|b37_mm10/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR b37-MM10 HYBRID WITH PREFIX $conf[1] $!";
	    }
	}
	$B37_MM10_HYBRID_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /hg19_fasta/i){
	if(!-e "$conf[1]"){
	    if($species =~ /^hg19$/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$HG19_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /hg19_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /^hg19$/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$HG19_FAI = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm10_hybrid_fasta/i){
	if(!-e "$conf[1]"){
	    if($species =~ /hybrid|hg19_mm10/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$HG19_MM10_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm10_hybrid_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /hybrid|hg19_mm10/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$HG19_MM10_HYBRID_FAI = $conf[1];
    }
    elsif($conf[0] =~ /hg19_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /^hg19$/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19 WITH PREFIX $conf[1] $!";
	    }
	}
	$HG19_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm10_hybrid_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /hybrid|hg19_mm10/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19-MM10 HYBRID WITH PREFIX $conf[1] $!";
	    }
	}
	$HG19_MM10_HYBRID_BWA_INDEX = $conf[1];
    }
}
close CONFIG;

my $REF_SEQ = '';
my $REF_FAI = '';
my $BWA_INDEX = '';
my $DB_SNP = "";
my $MILLS_1000G = '';
my $HAPMAP = '';
my $OMNI_1000G = '';
my $PHASE1_SNPS_1000G = '';
my $COSMIC = '';
my $ABRA_TARGETS = '';
if($species =~ /^b37$|human/i){
    $REF_SEQ = "$B37_FASTA";
    $REF_FAI = "$B37_FAI";
    $BWA_INDEX = "$B37_BWA_INDEX";
    $DB_SNP = "$Bin/data/b37/dbsnp_138.b37.excluding_sites_after_129.vcf";
    $MILLS_1000G = "$Bin/data/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
    $HAPMAP = "$Bin/data/b37/hapmap_3.3.b37.vcf";
    $OMNI_1000G = "$Bin/data/b37/1000G_omni2.5.b37.vcf";
    $PHASE1_SNPS_1000G = "$Bin/data/b37/1000G_phase1.snps.high_confidence.b37.vcf";
    $COSMIC = "$Bin/data/b37/CosmicCodingMuts_v67_b37_20131024__NDS.vcf";
    $ABRA_TARGETS = "$Bin/targets/abra/abra_target_regions_b37.bed";
}
elsif($species =~ /hybrid|b37_mm10/i){
    $REF_SEQ = "$B37_MM10_HYBRID_FASTA";
    $REF_FAI = "$B37_MM10_HYBRID_FAI";
    $BWA_INDEX = "$B37_MM10_HYBRID_BWA_INDEX";
    $DB_SNP = "$Bin/data/b37/dbsnp_138.b37.excluding_sites_after_129.vcf";
    $MILLS_1000G = "$Bin/data/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
    $HAPMAP = "$Bin/data/b37/hapmap_3.3.b37.vcf";
    $OMNI_1000G = "$Bin/data/b37/1000G_omni2.5.b37.vcf";
    $PHASE1_SNPS_1000G = "$Bin/data/b37/1000G_phase1.snps.high_confidence.b37.vcf";
    $COSMIC = "$Bin/data/b37/CosmicCodingMuts_v67_b37_20131024__NDS.vcf";
    $ABRA_TARGETS = "$Bin/targets/abra/abra_target_regions_b37.bed";
}
elsif($species =~ /^hg19$/i){
    $REF_SEQ = "$HG19_FASTA";
    $REF_FAI = "$HG19_FAI";
    $BWA_INDEX = "$HG19_BWA_INDEX";
    $DB_SNP = "$Bin/data/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf";
    $MILLS_1000G = "$Bin/data/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf";
    $HAPMAP = "$Bin/data/hg19/hapmap_3.3.hg19.vcf";
    $OMNI_1000G = "$Bin/data/hg19/1000G_omni2.5.hg19.vcf";
    $PHASE1_SNPS_1000G = "$Bin/data/hg19/1000G_phase1.snps.high_confidence.hg19.vcf";
    $COSMIC = "$Bin/data/hg19/CosmicCodingMuts_v67_20131024.vcf";
    $ABRA_TARGETS = "$Bin/targets/abra/abra_target_regions_hg19.bed";
}
elsif($species =~ /hg19_mm10/i){
    $REF_SEQ = "$HG19_MM10_HYBRID_FASTA";
    $REF_FAI = "$HG19_MM10_HYBRID_FAI";
    $BWA_INDEX = "$HG19_MM10_HYBRID_BWA_INDEX";
    $DB_SNP = "$Bin/data/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf";
    $MILLS_1000G = "$Bin/data/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf";
    $HAPMAP = "$Bin/data/hg19/hapmap_3.3.hg19.vcf";
    $OMNI_1000G = "$Bin/data/hg19/1000G_omni2.5.hg19.vcf";
    $PHASE1_SNPS_1000G = "$Bin/data/hg19/1000G_phase1.snps.high_confidence.hg19.vcf";
    $COSMIC = "$Bin/data/hg19/CosmicCodingMuts_v67_20131024.vcf";
    $ABRA_TARGETS = "$Bin/targets/abra/abra_target_regions_hg19.bed";
}
elsif($species !~ /b37|hg19|hybrid|b37_mm10|hg19_mm10/){
    die "ONLY SUPPORT FOR b37, hg19 or hybrd assemblies $!";
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
my $target_design = "$Bin/targets/$targets/$targets\__DESIGN.berger";
### std normals prefix; $targets\_norms_title.txt and $targets\_norms\_ALL_intervalnomapqcoverage_loess.txt required
my $target_std_normals = "$Bin/targets/$targets/StdNormals/$targets\_norms";
if(-d $targets){
    my @path = split(/\//, $targets);
    my $assay = pop @path;
    $targets_bed_padded = "$targets/$assay\_targets_plus5bp.bed";
    $target_design = "$targets/$assay\__DESIGN.berger";
    $target_std_normals = "$targets/StdNormals/$assay\_norms";
}

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

my @ref_chr = ();
open(FAI, "$REF_FAI") or die "CAN'T OPEN FAI FILE $REF_FAI $!";
while(<FAI>){
    chomp;

    my @line = split(/\s+/, $_);
    push @ref_chr, $line[0];
}
close FAI;

if($allSomatic){
    $scalpel = 1;
    $somaticsniper = 1;
    $strelka = 1;
    $varscan = 1;
    $virmid = 1;
}

my $count = 0;
my %inputFiles = ();
my %processedBams = ();
my @finalBams = ();
my %ran_pr_glob = 0;
my @prg_jids = ();
my $ran_ssf = 0;
my @ssf_jids = ();
my @all_jids = ();

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
        my $samp = $sn[-3];
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
		###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx100g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $ABRA/abra.jar --in $aiBams --out $aoBams --ref $REF_SEQ --bwa-ref $BWA_INDEX --targets $ABRA_TARGETS --working $output/intFiles/abra_$gpair[0] --threads 12`;
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/abra_wrapper.pl -inBams $aiBams -outBams $aoBams -refSeq $REF_SEQ -bwaRef $BWA_INDEX -targets $ABRA_TARGETS -working $output/intFiles/abra_$gpair[0]\_$c -config $config -log $output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA_WRAPPER.log`;
		
		$abra_jid = "$pre\_$uID\_$gpair[0]\_$c\_ABRA";
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA.done`;
		$ran_abra = 1;
	    }
	    
	    my $ran_fm = 0;
	    my @fm_jids = ();
	    my @fm_bams = ();
	    my $bcount = 0;
	    my $bcount = 0;
	    foreach my $outBam (@outBams){
		$bcount++;
		if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE.done" || $ran_abra){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE", job_hold => "$abra_jid", cpu => "1", mem => "50", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "5");
		    my $additionalParams = Schedule::additionalParams(%addParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar FixMateInformation I=$outBam O=$outBam\_FM.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $fmb O=$output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_MERGE_$c\_FM.done`;
		push @ir_jids, "$pre\_$uID\_$gpair[0]\_MERGE_$c\_FM";
		$ran_ir = 1;
	    }
	}
	else{
	    my $ran_rtc = 0;
	    my $rtc_jid = '';
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_RTC.done" || $step1){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_RTC", cpu => "10", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_RTC.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ -L $c $multipleTargets --known $MILLS_1000G --known $DB_SNP -S LENIENT -nt 10 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigner.intervals $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_RTC.done`;
		$rtc_jid = "$pre\_$uID\_$gpair[0]\_$c\_RTC";
		$ran_rtc = 1;
	    }
	    
	    ### seems to randomly have issues with file locking timeout on the m nodes so not submitting to it anymore
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.done" || $ran_rtc){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_IR", job_hold => "$rtc_jid", cpu => "1", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ -L $c $multipleTargets --knownAlleles $MILLS_1000G --knownAlleles $DB_SNP -S LENIENT --targetIntervals $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.done`;
		push @ir_jids, "$pre\_$uID\_$gpair[0]\_$c\_IR";
		$ran_ir = 1;
	    }
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO -R $REF_SEQ -S LENIENT --knownSites $DB_SNP --knownSites $MILLS_1000G --knownSites $HAPMAP --knownSites $OMNI_1000G --knownSites $PHASE1_SNPS_1000G --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_$gpair[0]\_recal_data.grp $irBams`;
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
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R $REF_SEQ -L $c $multipleTargets --emit_original_quals -BQSR $output/intFiles/$pre\_$gpair[0]\_recal_data.grp --num_cpu_threads_per_data_thread 6 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned_recal.bam -I $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam`;
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $irBams1 O=$output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_MERGE_PR.done`;
	push @merge_jids, "$pre\_$uID\_$gpair[0]\_MERGE_PR";
	$ran_m = 1;
    }

    my $mj = join(",", @merge_jids);
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.done" || $ran_m){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_SSF", job_hold => "$mj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T SplitSamFile -R $REF_SEQ -I $output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam --outputRoot $output/alignments/$pre\_indelRealigned_recal_`;
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MeanQualityByCycle INPUT=$finalBam OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.txt CHART_OUTPUT=$output/intFiles/$pre\_MeanQualityByCycle_$samp.pdf REFERENCE_SEQUENCE=$REF_SEQ VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true TMP_DIR=/scratch/$uID`;
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
    exit(0);
}

`/bin/mkdir -m 775 -p $output/variants`;
`/bin/mkdir -m 775 -p $output/variants/snpsIndels`;    
`/bin/mkdir -m 775 -p $output/variants/snpsIndels/haplotypecaller`;
my @ugVariants = ();
my @hcVariants = ();
my $ran_ug = 0;
my $ran_hc = 0;
my @ug_jids = ();
my @hc_jids = ();
my $prgj = join(",", @prg_jids);
foreach my $c (@ref_chr){
    ### NOTE1: GET THIS RANDOM ERROR NOW FOR SNP CALLS
    ###       Unexpectedly couldn't find valid codec for temporary output file
    ###       RERUNNING THE SAME COMMAND SEEMS TO FIX IT
    ###
    ### NOTE2: GET SYSTEM MEMORY ISSUE WHEN RUNNING ON THE CLUSTER AND LOTS
    ###        OF INPUT FILES TO UNIFIEDGENOTYPE WITH NUM THREADS > 1 -nt
    ###        NOT A PROBLEM WHEN RUNNING IT WITH LOTS OF INPUT FILES ON RAY
    
    my $irBams2 = join(" ", @{$processedBams{"$c"}});

    if($ug){
	if(!-e "$output/progress/$pre\_$uID\_$c\_UG.done" || $ran_pr_glob{"$c"}){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$c\_UG", job_hold => "$prgj", cpu => "12", mem => "24", cluster_out => "$output/progress/$pre\_$uID\_$c\_UG.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name $species -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model BOTH --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_$c\_UnifiedGenotyper.vcf $irBams2`;
	    `/bin/touch $output/progress/$pre\_$uID\_$c\_UG.done`;
	    push @ug_jids, "$pre\_$uID\_$c\_UG";
	    $ran_ug = 1;
	}
    }

    ### NOTE: ANNOTATIONS THAT DON'T WORK:AlleleBalance, HardyWeinberg,IndelType
    if(!-e "$output/progress/$pre\_$uID\_$c\_HC.done" || $ran_pr_glob{"$c"}){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$c\_HC", job_hold => "$prgj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_$c\_HC.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx90g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 24 --out $output/intFiles/$pre\_$c\_HaplotypeCaller.vcf $irBams2`;
	`/bin/touch $output/progress/$pre\_$uID\_$c\_HC.done`;
	push @hc_jids, "$pre\_$uID\_$c\_HC";
	$ran_hc = 1;
    }
    
    push @ugVariants, "--variant $output/intFiles/$pre\_$c\_UnifiedGenotyper.vcf";
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
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_HaplotypeCaller_RAW.vcf --assumeIdenticalSamples $hcVars`;
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
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_SNP_HC", job_hold => "$cvhcj", cpu => "4", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_VR_SNP_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1000G -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE1_SNPS_1000G -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DB_SNP -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_SNP.plots.R -nt 4`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_SNP_HC.done`;
    $vrshcj = "$pre\_$uID\_VR_SNP_HC";
    $ran_vr_snp_hc = 1;
}

### NOTE: sometimes throws an error when running with multiple threads
###       about not being able to find a tmp file; running with -nt 1 to avoid errors
my $ran_ar_snp_hc = 0;
my $arshcj = '';
if(!-e "$output/progress/$pre\_$uID\_AR_SNP_HC.done" || $ran_vr_snp_hc){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_SNP_HC", job_hold => "$vrshcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_SNP_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -o $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_SNP_HC.done`;
    $arshcj = "$pre\_$uID\_AR_SNP_HC";
    $ran_ar_snp_hc = 1;
}

my $ran_vr_indel_hc = 0;
my $vrihcj = '';
if(!-e "$output/progress/$pre\_$uID\_VR_INDEL_HC.done" || $ran_ar_snp_hc){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_INDEL_HC", job_hold => "$arshcj", cpu => "4", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_VR_INDEL_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS_1000G -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.plots.R -nt 4`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_INDEL_HC.done`;
    $vrihcj = "$pre\_$uID\_VR_INDEL_HC";
    $ran_vr_indel_hc = 1;   
}

my $ran_ar_indel_hc = 0;
my $arihcj = '';
if(!-e "$output/progress/$pre\_$uID\_AR_INDEL_HC.done" || $ran_vr_indel_hc){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_INDEL_HC", job_hold => "$vrihcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_INDEL_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -o $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_INDEL_HC.done`;
    $arihcj = "$pre\_$uID\_AR_INDEL_HC";
    $ran_ar_indel_hc = 1;   
}

# This runs regardless because it has a .done check in the function
sleep(2);
&generateMaf("$output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'haplotypecaller', "$arihcj,$ssfj", $ran_ar_indel_hc);

if($ug){
    `/bin/mkdir -m 775 -p $output/variants/snpsIndels/unifiedgenotyper`;
    my $ugVars = join(" " , @ugVariants);
    my $ran_cv_ug = 0;
    my $cvugj = '';
    my $ugj = join(",", @ug_jids);
    if(!-e "$output/progress/$pre\_$uID\_CV_UG.done" || $ran_ug){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_UG", job_hold => "$ugj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf --assumeIdenticalSamples $ugVars`;
	`/bin/touch $output/progress/$pre\_$uID\_CV_UG.done`;
	$cvugj = "$pre\_$uID\_CV_UG";
	$ran_cv_ug = 1;
    }

    my $ran_vr_snp_ug = 0;
    my $vrsugj = '';
    if(!-e "$output/progress/$pre\_$uID\_VR_SNP_UG.done" || $ran_cv_ug){  
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_SNP_UG", job_hold => "$cvugj", cpu => "4", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_VR_SNP_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1000G -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE1_SNPS_1000G -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DB_SNP -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.plots.R -nt 4`;
	`/bin/touch $output/progress/$pre\_$uID\_VR_SNP_UG.done`;
	$vrsugj = "$pre\_$uID\_VR_SNP_UG";
	$ran_vr_snp_ug = 1;
    }

    my $ran_ar_snp_ug = 0;
    my $arsugj = '';
    if(!-e "$output/progress/$pre\_$uID\_AR_SNP_UG.done" || $ran_vr_snp_ug){  
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_SNP_UG", job_hold => "$vrsugj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_SNP_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -o $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -nt 1`;
	`/bin/touch $output/progress/$pre\_$uID\_AR_SNP_UG.done`;
	$arsugj = "$pre\_$uID\_AR_SNP_UG";
	$ran_ar_snp_ug = 1;
    }

    my $ran_vr_indel_ug = 0;
    my $vriugj = '';
    if(!-e "$output/progress/$pre\_$uID\_VR_INDEL_UG.done" || $ran_ar_snp_ug){  
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_INDEL_UG", job_hold => "$arsugj", cpu => "4", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_VR_INDEL_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS_1000G -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.plots.R -nt 4`;
	`/bin/touch $output/progress/$pre\_$uID\_VR_INDEL_UG.done`;
	$vriugj = "$pre\_$uID\_VR_INDEL_UG";
	$ran_vr_indel_ug = 1;
    }

    my $ran_ar_indel_ug = 0;
    my $ariugj = '';
    if(!-e "$output/progress/$pre\_$uID\_AR_INDEL_UG.done" || $ran_vr_indel_ug){  
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_INDEL_UG", job_hold => "$vriugj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_AR_INDEL_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -o $output/variants/snpsIndels/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf -nt 1`;
	`/bin/touch $output/progress/$pre\_$uID\_AR_INDEL_UG.done`;
	$ariugj = "$pre\_$uID\_AR_INDEL_UG";
	$ran_ar_indel_ug = 1;
    }

    ###if(!-e "$output/progress/$pre\_$uID\_UG_MAF.done" || $ran_ar_indel_ug){  
	###sleep(2);
	###&generateMaf("$output/variants/snpsIndels/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf", 'unifiedgenotyper', "$ariugj");
	###`/bin/touch $output/progress/$pre\_$uID\_UG_MAF.done`;
    ###}
}

my $hasPair = 0;

if($pair){
    `/bin/mkdir -m 775 -p $output/variants/snpsIndels/haplotect`;
    `/bin/mkdir -m 775 -p $output/variants/copyNumber`;
    `/bin/mkdir -m 775 -p $output/variants/copyNumber/facets`;
    `/bin/mkdir -m 775 -p $output/variants/structVar`;
    `/bin/echo "Tumor_Sample_Barcode\tRdata_filename" > $output/variants/copyNumber/facets/facets_mapping.txt `;
    ###`/bin/mkdir -m 775 -p $output/variants/varscan`;

    open(PAIR, "$pair") or die "Can't open $pair file";
    my %submitted_lns = ();
    my @mu_jids = ();
    my $ran_mutect_glob = 0;
    my $haplotect_run = 0;
    my $facets_run = 0;
    my @facets_jid = ();
    ### NOTE: THE SAMPLE NAMES IN THE PAIRING FILE MUST MATCH EXACTLY THE SAMPLE NAMES IN THE REALIGNED/RECALIBRATED BAM FILE
    while(<PAIR>){
	chomp;
	
	my @data = split(/\s+/, $_);
	### pairing file also contains unpaired samples with NA for the other pai
	### so we can keep track of what sample is tumor or normal
	if($data[0] =~ /^NA$/i || $data[1] =~ /^NA$/i){
	    next;
	}
        ## this means that there really is a paired sample set
        $hasPair=1;
	
	### NOTE: NOT AUTOMATICALLY RUNNING EACH SOMATIC CALLER WHEN SSF RUNS
	###       BCAUSE DON'T WANT TO KICK IT OFF FOR ALL SAMPLE PAIRS
	###       IN CASE ONLY ONE SAMPLE GROUP GOT MODIFIED ABOVE

	###       WILL RUN SOMATIC ANALYSIS FOR ALL SAMPLE PAIRS
	###       IF JUST ONE SAMPLE HAD ITS BAM MODIFIED
	`/bin/mkdir -m 775 -p $output/variants/snpsIndels/mutect`;
	my $ran_mutect = 0;
	my $mutectj = '';
	### MUTECT will fail if order of ref seq contigs and vcf don't match
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done" || $ran_ssf){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT", job_hold => "$ssfj", cpu => "2", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xmx4g -Djava.io.tmpdir=/scratch/$uID -jar $MUTECT/muTect.jar --analysis_type MuTect --reference_sequence $REF_SEQ --dbsnp $DB_SNP --cosmic $COSMIC --input_file:normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --input_file:tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --vcf $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf --out $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.txt -rf BadCigar --enable_extended_output --downsampling_type NONE`;
	    `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done`;
	    $mutectj = "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT";
	    push @mu_jids, "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT";
	    push @all_jids, $mutectj;
	    $ran_mutect = 1;
	    $ran_mutect_glob = 1;
	}

	###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done" || $ran_mutect){
	### &generateMaf("$output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf", 'mutect', "$mutectj", $data[0], $data[1]);
	###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT_MAF.done`;
	###}

	if($somaticsniper){
	    `/bin/mkdir -m 775 -p $output/variants/snpsIndels/somaticsniper`;
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
	    ### my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF", job_hold => "$ssj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.log");
	    ###my $standardParams = Schedule::queuing(%stdParams);
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf -c somaticsniper -o $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.done`;
	    ###}
	}

	if($varscan){
	    ###my $ran_mpileup = 0;
	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MPILEUP.done" || $ran_ssf){  
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $ssfj -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -A -f $REF_SEQ -d 500000 -o $output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam`;
	    
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -hold_jid $ssfj -pe alloc 2 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $SAMTOOLS/samtools mpileup -A -f $REF_SEQ -d 500000 -o $output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MPILEUP.done`;
	    ###$ran_mpileup = 1;
	    ###}	    

	    ###my $ran_varscan_somatic = 0;
	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC.done" || $ran_mpileup){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -pe alloc 2 -l virtual_free=5G -q lau.q,lcg.q,nce.q $Bin/qCMD $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $VARSCAN/VarScan.jar somatic $output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup $output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic --strand-filter 1 --output-vcf 1`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC.done`;
	    ###$ran_varscan_somatic = 1;
	    ###}

	    ###if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done" || $ran_varscan_somatic){
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.snp_MAF.txt -n $data[0] -t $data[1]`;
	    
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -pe alloc 1 -l virtual_free=2G -q lau.q,lcg.q,nce.q $Bin/qCMD $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.indel.vcf -c varscan -o $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic\.indel_MAF.txt -n $data[0] -t $data[1]`;
	    ###`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC_MAF.done`;
	    ###}

	}
	### BECAUSE THE OUTPUT PREFIX FOR VIRMID IS THE DISEASE BAM,
	### NEED TO CREATE A DIRECTORY FOR EACH TUMOR/NORMAL PAIR
	### SINCE A DISEASE BAM CAN BE USED IN MULTIPLE COMPARISONS
	### virmid fails if directory already exists
	### NOTE: strelka doesn't work on hybrid because asembly is too large
	if($virmid){
	    my $ran_virmid = 0;
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done" || $ran_ssf){  
		sleep(2);
		if(-d "$output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid"){
		    `/bin/rm -rf $output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid`;
		}
		`/bin/mkdir -m 775 -p $output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid`;
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_VIRMID", job_hold => "$ssfj", cpu => "4", mem => "12", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx12g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $VIRMID/Virmid.jar -R $REF_SEQ -D $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam -N $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam -t 4 -o $pre\_$data[0]\_$data[1]\_virmid -w $output/variants/snpsIndels/virmid/$data[0]\_$data[1]\_virmid`;
		`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_VIRMID.done`;
		push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_VIRMID";
		$ran_virmid = 1;
	    }
	}

	if($strelka){
	    my $ran_strelka_config = 0;
	    if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.done" || $ran_ssf){  
		if(-d "$output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka"){
		    ### strelka DIES IF DIR ALREADY EXISTS
		    `/bin/rm -rf $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka`;
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
		### NOTE: strelka ONLY HAS CONFIG FOR BWA ALN, NOT SURE HOW IT WILL WORK WITH BWA MEM
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", job_hold => "$lnsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $STRELKA/bin/configureStrelkaWorkflow.pl --normal=$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor=$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --ref=$REF_SEQ --config=$STRELKA/etc/strelka_config_bwa_default.ini --output-dir=$output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka`;
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
		`/bin/mkdir -m 775 -p $output/variants/snpsIndels/scalpel/$data[0]\_$data[1]\_scalpel`;
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


        ## Here we will add the facets scripts
        ## These are the #'s Nick uses
        my $MINCOV=0;
        my $BASEQ=20;
        my $MAPQ=15;
	
        ## Set up tumor and normal counts
        `/bin/mkdir -m 775 -p $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets`;
        `/bin/mkdir -m 775 -p $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp`;
        my $facetsSETUP_jid = '';
        my $facets_setup = 0;
        if($hasPair && (! -e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.done" || $ssfj )) {
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP_T",  cpu => "4", mem => "5", job_hold => "$ssfj", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP_T.log");
            my $standardParams = Schedule::queuing(%stdParams);
            my %addParams = (runtime => "30");
            my $additionalParams = Schedule::additionalParams(%addParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $Bin/facets/bin/GetBaseCounts --thread 4 --filter_improper_pair 0 --sort_output --fasta $REF_SEQ --vcf $DB_SNP --maq $MAPQ --baq $BASEQ --cov $MINCOV --bam $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --out $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_indelRealigned_recal\_$data[1].dat`;
	    
            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP_N",  cpu => "4", mem => "5", job_hold => "$ssfj", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP_N.log");
            my $standardParams2 = Schedule::queuing(%stdParams);
            %addParams = (runtime => "30");
            my $additionalParams2 = Schedule::additionalParams(%addParams);
            `$standardParams2->{submit} $standardParams2->{job_name} $standardParams2->{cpu} $standardParams2->{mem} $standardParams2->{job_hold} $standardParams2->{cluster_out} $additionalParams2 $Bin/facets/bin/GetBaseCounts --thread 4 --filter_improper_pair 0 --sort_output --fasta $REF_SEQ --vcf $DB_SNP --maq $MAPQ --baq $BASEQ --cov $MINCOV --bam $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --out $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_indelRealigned_recal\_$data[0].dat`;

            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_merge_counts_facets_SETUP",  cpu => "4", mem => "18", job_hold => "$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP_N,$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP_T", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.log");
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
	    
            my $sub_species = $species eq 'b37' ? 'hg19' : $species;

            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN",  cpu => "3", mem => "2", job_hold => "$facetsSETUP_jid", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.log");
            my $standardParams = Schedule::queuing(%stdParams);
            my %addParams = (runtime => "10");
            my $additionalParams = Schedule::additionalParams(%addParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $Bin/facets/facets_RUN.sh $FACETS_SUITE $FACETS_LIB $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets $data[0]\_$data[1] $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_countsMerged_$data[0]\_$data[1].dat $sub_species 300 100`;
	    push @facets_jid, "$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN" ;
            $facets_run = 1;
            `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.done`; 
        }
    `/bin/echo "$data[1]\t$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/$data[0]\_$data[1]\_hisens.Rdata" >> $output/variants/copyNumber/facets/facets_mapping.txt`;
    }
    close PAIR;

    ## run merge seg script
    my $facets_haplotect_jid = ''; 
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

    if($hasPair && (!-e "$output/progress/$pre\_$uID\_HAPLOTECT.done" || $ran_mutect_glob || $ran_ar_indel_hc)){
	sleep(2);
        my $patientFile = "";
        if($patient){
            $patientFile = "-patient $patient";
        }
	my $muj = join(",", @mu_jids);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HAPLOTECT", job_hold => "$arihcj,$muj", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_HAPLOTECT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/haploTect_merge.pl -pair $pair -hc_vcf $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf -species $species -pre $pre -output $output/variants/snpsIndels/haplotect -mutect_dir $output/variants/snpsIndels/mutect -config $config $patientFile -align_dir $output/alignments/ -svnRev $svnRev -delete_temp`;

        $haplotect_run = 1;
	`/bin/touch $output/progress/$pre\_$uID\_HAPLOTECT.done`;
        $facets_haplotect_jid .= ",$pre\_$uID\_HAPLOTECT";
    }

    ## Now join maf
    my $mafAnnoRun=0;
    if($hasPair && (!-e "$output/progress/$pre\_$uID\_join_maf.done" || $haplotect_run || $facets_run)){
        sleep(2);

       my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_join_maf", job_hold => "$facets_haplotect_jid", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_join_maf.log");
        my $standardParams = Schedule::queuing(%stdParams); 
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $FACETS_SUITE/facets mafAnno -m $output/variants/snpsIndels/haplotect/$pre\_haplotect_VEP_MAF.txt -f $output/variants/copyNumber/facets/facets_mapping.txt -o $output/intFiles/$pre\_CMO_MAF_intermediate.txt`; 
        `/bin/touch $output/progress/$pre\_$uID\_join_maf.done`;
	push @all_jids, "$pre\_$uID\_join_maf";
        $mafAnnoRun = 1;
    }

    if($hasPair && (! -e "$output/progress/$pre\_$uID\_wes_filters.done" || $mafAnnoRun )) {
        my $holdVar = '';
        if( $mafAnnoRun) { 
            $holdVar = "$pre\_$uID\_join_maf";
        }
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_wes_filters", job_hold => "$holdVar", cpu => "2", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_wes_filters.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my %addParams = (scheduler => "$scheduler", runtime => "7", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "0");
        my $additionalParams = Schedule::additionalParams(%addParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/maf/post_filters.pl -in_maf $output/intFiles/$pre\_CMO_MAF_intermediate.txt -out_maf $output/variants/$pre\_CMO_MAF.txt -config $config -species $species -filter_ffpe -blacklist -low_conf`;
        `/bin/touch $output/progress/$pre\_$uID\_wes_filters.done`;
        push @all_jids, "$pre\_$uID\_wes_filters";
    }

    open(GROUP, "$group") || die "CAN'T OPEN SAMPLE GROUPING FILE $group $!";
    open(BLIST, ">$output/intFiles/$pre\_sv_bam_list.txt") || die "CAN'T WRITE TO SAMPLE BAM LIST FILE $output/intFiles/$pre\_sv_bam_list.txt $!";
    while(<GROUP>){
	chomp;
	
	my @data = split(/\s+/, $_);
	print BLIST "$data[0]\t$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam\n";
    }
    close GROUP;
    close BLIST;

    if(!-e "$output/progress/$pre\_$uID\_DMP_CNV.done" || $ran_ssf){
	sleep(2);
  
	### NOTE: DMP CNV only supported for impact410/hemepactv3 because they are the only ones that we have standard normals bams for
	if($targets =~ /IMPACT410|HemePACT_v3/i && -e $patient && -e "$target_std_normals\_title.txt" && -e "$target_std_normals\_ALL_intervalnomapqcoverage_loess.txt"){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DMP_CNV", job_hold => "$ssfj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_DMP_CNV.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/exome_cnv.pl -pre $pre -result $output/variants/copyNumber/dmp_cnv -berger $target_design -bamlist $output/intFiles/$pre\_sv_bam_list.txt -patient $patient -std_covg $target_std_normals -genome $species -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;
	    `/bin/touch $output/progress/$pre\_$uID\_DMP_CNV.done`;
	    push @all_jids, "$pre\_$uID\_DMP_CNV";
	}
    }
    
    my $ran_strvar = 0;
    if($hasPair && (!-e "$output/progress/$pre\_$uID\_DELLY.done" || $ran_ssf)){
	sleep(2);
	if(-d "$output/variants/structVar/delly"){
	    `/bin/rm -rf $output/variants/structVar/delly`;
	}
	
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DELLY", job_hold => "$ssfj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_DELLY.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/RunStructuralVariantPipeline_Delly.pl -pre $pre -out $output/variants/structVar/delly -pair $pair -bam_list $output/intFiles/$pre\_sv_bam_list.txt -genome $species -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;
	`/bin/touch $output/progress/$pre\_$uID\_DELLY.done`;
	$ran_strvar = 1;
    }
    
    if($hasPair && (!-e "$output/progress/$pre\_$uID\_CDNA_CONTAM.done" || $ran_strvar)){
	sleep(2);

        my $hold_value = $ran_strvar ? "$pre\_$uID\_DELLY" : "";
	
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CDNA_CONTAM", job_hold => "$hold_value", cpu => "1", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_CDNA_CONTAM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/check_cDNA_contamination.py -s $output/variants/structVar/delly/$pre\_AllAnnotatedSVs.txt -o $output/metrics/$pre\_cDNA_contamination.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_CDNA_CONTAM.done`;
	push @all_jids, "$pre\_$uID\_CDNA_CONTAM";
    }
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
        if(!-e "$output/progress/$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED.done" || $ran_hc){
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

            `/bin/mkdir -m 775 -p $vcf_dir/chrom_$c`;
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_split_CHR_$c", job_hold => $bgz_jid, cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_split_$c.log");
            my $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $BCFTOOLS/bcftools filter -r $c $vcf.gz -O v -o $vcf_dir/chrom_$c/$jna\_$c.vcf`;

            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED", job_hold => "$pre\_$uID\_$jna\_split_CHR_$c", cpu => "4", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED.log");
            $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/generateMAF.pl -vcf $vcf_dir/chrom_$c/$jna\_$c.vcf -species $species -config $config -caller $type $patientFile -align_dir $output/alignments -delete_temp`;
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
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_MAF";

        $merge_files =~ s/UNPAIRED_TCGA_MAF.txt/UNPAIRED_TCGA_PORTAL_MAF.txt/g;
        %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF.log");
        $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF";

        $merge_files =~ s/UNPAIRED_TCGA_PORTAL_MAF.txt/UNPAIRED_VEP_MAF.txt/g;
        %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_VEP_MAF", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_VEP_MAF.log");
        $standardParams = Schedule::queuing(%stdParams);
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
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_cleanup", job_hold => "$merge_holds", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$jna\_cleanup.log");
        $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams rm -rf $vcf_dir/chrom_* $vcf.gz $vcf.gz.csi`;
        push @all_jids, "$pre\_$uID\_$jna\_cleanup";
    }
}
