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

my ($patient, $email, $impact, $wes, $svnRev, $pair, $group, $bamgroup, $config, $nosnps, $targets, $ug, $scheduler, $priority_project, $priority_group, $abra, $indelrealigner, $noindelrealign, $help, $step1, $allSomatic, $scalpel, $somaticsniper, $strelka, $varscan, $virmid, $lancet, $vardict, $pindel, $abra_target, $rna, $mutect2, $nofacets);

my $pre = 'TEMP';
my $output = "results";
my $species = 'b37';

my $uID = `/usr/bin/id -u -n`;
chomp $uID;
my $rsync = "/juno/res/bic/$uID";
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
	    'ug|unifiedgenotyper' => \$ug,
	    'abra' => \$abra,
	    'impact' => \$impact,
            'wes' => \$wes,
	    'indelrealigner' => \$indelrealigner,
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
	    'lancet' => \$lancet,
 	    'vardict' => \$vardict,
            'pindel' => \$pindel,
	    'tempdir=s' => \$tempdir,
 	    'output|out|o=s' => \$output,
            'abratarget|abra_target=s' => \$abra_target,
            'rna=s' => \$rna,
            'mutect2' => \$mutect2,
            'noindelrealign|noir' => \$noindelrealign,
            'nofacets' => \$nofacets) or exit(1);


if(!$group || !$config || !$scheduler || !$targets || !$bamgroup || $help){
    print <<HELP;

    USAGE: process_alignments_human.pl -group GROUP -bampgroup BAMBROUP -config CONFIG -scheduler SCHEDULER -targets TARGETS
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED)
	* BAMGROUP: files listing bams to be processed together; every bam for each group on 1 line, comma-separated (required)
	* SPECIES: b37, human_hybrid|xenograft (b37_mm10)
	* TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS unless given full path to targets directory (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE, LUNA, and JUNO (REQUIRED)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
	* PRE: output prefix (default: TEMP)
	* OUTPUT: output results directory (default: results)
	* RSYNC:  path to rsync data for archive (default: /juno/res/bic/USER_ID)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* -nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed
	* -abra: run abra to realign indels (default)
	* -indelrealigner: run GATK indelrealigner (abra runs as default)
        * -noindelrealign: skip indel realignment
        * -nofacets: skip facets
	* -step1: forece the pipeline to start from the first step in pipeline
	* haplotypecaller is default; -ug || -unifiedgenotyper to also make unifiedgenotyper variant calls	
	* TEMPDIR:  temp directory (default: /scratch/$uID)
	* ALLSOMATIC: run all somatic callers; mutect/haplotypecaller always run; otherwise -scalpel, -somaticsniper, -strelka, -varscan, -virmid, -lancet, -vardict, -pindel, -mutect2 to run them individually	
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

if($wes && $impact){
    die "Cannot run as both a wes and an impact. Please select -wes OR -impact";
}elsif(!$wes && !$impact){
    die "Must specify to run the project as either -wes or -impact";
}

if($abra && $indelrealigner){
    die "Cannot run both abra and gatk indelrealigner $!";
}
if(($abra || $indelrealigner) && $noindelrealign)
{
    die "Cannot skip indel realignment when -abra or -indelrealigner is specified $!";
}
if(!$abra && !$indelrealigner && !$noindelrealign){
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
    $mutect2 = 1;
}

my $ABRA = '';
my $BCFTOOLS= '';
my $BEDTOOLS = '';
my $GATK = '';
my $GATK4 = '';
my $HTSTOOLS = '';
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
my $VEP = '';
my $VIRMID = '';

my $JAVA = '';
my $JAVA7_MUTECT = '';
my $PYTHON = '';
my $PERL = '';

my $SINGULARITY = '';
my $singularityParams = '';
my $singularityBind = '';
my $singularityenv_prepend_path = "";

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
    elsif($conf[0] =~ /bcftools/i){
        if(!-e "$conf[1]/bcftools"){
            die "CAN'T FIND bcftools IN $conf[1] $!";
        }
        $BCFTOOLS = $conf[1];
    }
    elsif($conf[0] =~ /bedtools/i){
        if(!-e "$conf[1]/bedtools"){
            die "CAN'T FIND bedtools IN $conf[1] $!";
        }
        $BEDTOOLS = $conf[1];
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
    elsif($conf[0] =~ /^gatk$/i){
	if(!-e "$conf[1]/GenomeAnalysisTK.jar"){
	    die "CAN'T FIND GenomeAnalysisTK.jar IN $conf[1] $!";
	}
	$GATK = $conf[1];
    }
    elsif($conf[0] =~ /^gatk4$/i){
        if(!-e "$conf[1]/gatk"){
        #    die "CAN'T gatk IN $conf[1] $!";
        }
        $GATK4 = $conf[1];
    }
    elsif($conf[0] =~ /htstools/i){
        if(!-e "$conf[1]/snp-pileup"){
            die "CAN'T FIND htstools IN $conf[1] $!";
        }
        $HTSTOOLS = $conf[1];
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
    elsif($conf[0] =~ /^vep/i){
        if(!-e "$conf[1]/variant_effect_predictor.pl"){
            die "CAN'T FIND VEP IN $conf[1] $!";
        }
        $VEP = $conf[1];
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
	my $path_tmp = $ENV{'PATH'};
	$ENV{'JAVA_HOME'} = "";
        $ENV{'JAVA_JRE'} = "";
        $ENV{'JAVA_LD_LIBRARY_PATH'} = "";
        $ENV{'JAVA_LIBS'} = "";
        $singularityenv_prepend_path .= ":$conf[1]";
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
	my $path_tmp = $ENV{'PATH'};
	$ENV{'PATH'} = "$conf[1]:$path_tmp";
        $singularityenv_prepend_path .= ":$conf[1]";
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
	    if($species =~ /hybrid|xenograft|b37_mm10/i){
		die "CAN'T FIND $conf[1] $!";
	    }
	}
	$B37_MM10_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /b37_mm10_hybrid_fai/i){
	if(!-e "$conf[1]"){
	    if($species =~ /hybrid|xenograft|b37_mm10/i){
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
    elsif($conf[0] =~ /b37_mm10_hybrid_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /human_hybrid|b37_mm10/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $B37_MM10_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /b37_mm10_hybrid_bwa_index/i){
        if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
            if($species =~ /human_hybrid|b37_mm10/i){
                die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR b37-MM10 HYBRID WITH PREFIX $conf[1] $!";
            }
        }
        $B37_MM10_HYBRID_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /hg19_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /hg19/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $HG19_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /hg19_bwa_index/i){
        if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
            if($species =~ /hg19/i){
                die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19 WITH PREFIX $conf[1] $!";
            }
        }
        $HG19_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm10_hybrid_fasta/i){
        if(!-e "$conf[1]"){
            if($species =~ /hg19_hybrid/i){
                die "CAN'T FIND $conf[1] $!";
            }
        }
        $HG19_MM10_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm10_hybrid_bwa_index/i){
        if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
            if($species =~ /hg19_hybrid/i){
                die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19-MM9 HYBRID WITH PREFIX $conf[1] $!";
            }
        }
        $HG19_MM10_HYBRID_BWA_INDEX = $conf[1];
    }

}
my %sinParams = (singularity_exec => "$SINGULARITY/singularity", singularity_image => "$Bin/variants_pipeline_singularity_prod.simg");
$singularityParams = Schedule::singularityParams(%sinParams);
$singularityBind = Schedule::singularityBind($scheduler);

$ENV{'SINGULARITYENV_PREPEND_PATH'} = $singularityenv_prepend_path;
$ENV{'SINGULARITY_BINDPATH'} = $singularityBind;
close CONFIG;

my $REF_SEQ = '';
my $REF_FAI = '';
my $BWA_INDEX = '';
my $DB_SNP = "";
my $ExAC_VCF = '';
my $FACETS_DB_SNP = "";
my $MILLS_1000G = '';
my $HAPMAP = '';
my $OMNI_1000G = '';
my $PHASE1_SNPS_1000G = '';
my $COSMIC = '';
my $COSMIC_HOTSPOTS = '';
if($species =~ /^b37$|human/i){
    $DB_SNP = "$Bin/data/b37/dbsnp_138.b37.excluding_sites_after_129.vcf";
    $ExAC_VCF = "$VEP/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
    $FACETS_DB_SNP = "$Bin/data/b37/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf";
    $MILLS_1000G = "$Bin/data/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
    $HAPMAP = "$Bin/data/b37/hapmap_3.3.b37.vcf";
    $OMNI_1000G = "$Bin/data/b37/1000G_omni2.5.b37.vcf";
    $PHASE1_SNPS_1000G = "$Bin/data/b37/1000G_phase1.snps.high_confidence.b37.vcf";
    $COSMIC = "$Bin/data/b37/CosmicCodingMuts_v67_b37_20131024__NDS.vcf";
    $COSMIC_HOTSPOTS = "$Bin/data/b37/dmp_cosmic_for_hotspots.vcf";

    if (! -s $ExAC_VCF ){
        die "$ExAC_VCF not present! $!";
    }

    if($species =~ /human_hybrid|xenograft|b37_mm10/i){
        $species = "human_hybrid";
        $REF_SEQ = "$B37_MM10_HYBRID_FASTA";
        $REF_FAI = "$B37_MM10_HYBRID_FAI";
        $BWA_INDEX = "$B37_MM10_HYBRID_BWA_INDEX";
    } else {
        $species = "b37";
        $REF_SEQ = "$B37_FASTA";
        $REF_FAI = "$B37_FAI";
        $BWA_INDEX = "$B37_BWA_INDEX";
    }
    
}
elsif($species =~ /^hg19$/i){

    #die "hg19 is no longer supported in the variants pipeline";

    $REF_SEQ = "$HG19_FASTA";
    $REF_FAI = "$HG19_FAI";
    $BWA_INDEX = "$HG19_BWA_INDEX";
    $DB_SNP = "$Bin/data/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf";
    $ExAC_VCF = "$VEP/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
    $FACETS_DB_SNP = "$Bin/data/hg19/dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT.vcf";
    $MILLS_1000G = "$Bin/data/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf";
    $HAPMAP = "$Bin/data/hg19/hapmap_3.3.hg19.vcf";
    $OMNI_1000G = "$Bin/data/hg19/1000G_omni2.5.hg19.vcf";
    $PHASE1_SNPS_1000G = "$Bin/data/hg19/1000G_phase1.snps.high_confidence.hg19.vcf";
    $COSMIC = "$Bin/data/hg19/CosmicCodingMuts_v67_20131024.vcf";
    $COSMIC_HOTSPOTS = "$Bin/data/b37/dmp_cosmic_for_hotspots.vcf";
}
elsif($species =~ /hg19_mm10/i){

    die "hg19 is no longer supported in the variants pipeline";

    $REF_SEQ = "$HG19_MM10_HYBRID_FASTA";
    $REF_FAI = "$HG19_MM10_HYBRID_FAI";
    $BWA_INDEX = "$HG19_MM10_HYBRID_BWA_INDEX";
    $DB_SNP = "$Bin/data/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf";
    $FACETS_DB_SNP = "$Bin/data/hg19/dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT.vcf";
    $MILLS_1000G = "$Bin/data/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf";
    $HAPMAP = "$Bin/data/hg19/hapmap_3.3.hg19.vcf";
    $OMNI_1000G = "$Bin/data/hg19/1000G_omni2.5.hg19.vcf";
    $PHASE1_SNPS_1000G = "$Bin/data/hg19/1000G_phase1.snps.high_confidence.hg19.vcf";
    $COSMIC = "$Bin/data/hg19/CosmicCodingMuts_v67_20131024.vcf";
    $COSMIC_HOTSPOTS = "$Bin/data/b37/dmp_cosmic_for_hotspots.vcf";
}
elsif($species !~ /b37|hybrid|xenograft|b37_mm10/){
    die "ONLY SUPPORT FOR b37 or hybrd|xenograft assemblies $!";
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
my $targets_facet = "$Bin/targets/$targets/$targets\_targets_FACETS.ilist";

if(-d $targets){
    my @path = split(/\//, $targets);
    my $assay = pop @path;
    $targets_bed_padded = "$targets/$assay\_targets_plus5bp.bed";
    $target_design = "$targets/$assay\__DESIGN.berger";
    $targets_facet = "$targets/$assay\_targets_FACETS.ilist";
    if($abra){
        $target_std_normals = "$targets/StdNormals/$assay\_abra_norms";
    }
    else{
        $target_std_normals = "$targets/StdNormals/$assay\_gatk_norms";
    }
}
if(!-e "$targets_bed_padded"){
    die "CAN'T LOCATE $targets_bed_padded FOR $targets; REQUIRED FOR SCALPEL or respectively $!";
}
if(!$abra_target){
    $abra_target = $targets_bed_padded;
}
if(!-e "$abra_target"){
    die "CAN'T LOCATE $abra_target; REQUIRED FOR ABRA $!";
}

print "\n\n";
print "-abra_target $abra_target\n";
print "-targets $targets\n";
print "\n\n";

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
	    ### IN CASE A FILE SHOWS UP MULTIPLE TIMES DUE TO BEING IN MULTIPLE COMPARISONS AND WASN'T COLLAPSED e.g. MET ANALYSIS
	    ### THIS MAKES SURE THAT A FILE ISN'T INCLUDED MULTIPLE TIMES IN PROCESSING
	    next;
	}
	push @pins, "-I $pai";
	$inputFiles{$pai} = 1;

        ## store final bams for post-recalibration stats
        my $samp = $sn[0];
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
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $SAMTOOLS/samtools view -b -o $inB[1]\_$c\.bam $inB[1] $c`;
		
		    $chrsj = "$pre\_$uID\_$boi\_$c";		    
		    `/bin/touch $output/progress/$pre\_$uID\_$boi\_$c\.done`;
		    $ran_chr_split = 1;
		}

		if(!-e "$output/progress/$pre\_$uID\_$boi\_$c\_INDEX.done" || $ran_chr_split){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$boi\_$c\_INDEX", job_hold => "$chrsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$boi\_$c\_INDEX.log");
		    my $standardParams = Schedule::queuing(%stdParams);	    
		    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "2");
		    my $additionalParams = Schedule::additionalParams(%addParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $SAMTOOLS/samtools index $inB[1]\_$c\.bam`;
		    
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
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_ABRA", job_hold => "$cij", cpu => "12", mem => "120", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA.log");
		my $standardParams = Schedule::queuing(%stdParams);	    
		my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "4");
		my $additionalParams = Schedule::additionalParams(%addParams);
		###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx100g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $ABRA/abra.jar --in $aiBams --out $aoBams --ref $REF_SEQ --bwa-ref $BWA_INDEX --targets $ABRA_TARGETS --working $output/intFiles/abra_$gpair[0] --threads 12`;
		print "\$abra_target being passed to abra_wrapper.pl = $abra_target\n\n";
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/abra_wrapper.pl -inBams $aiBams -outBams $aoBams -refSeq $REF_SEQ -bwaRef $BWA_INDEX -targets $abra_target -working $output/intFiles/abra_$gpair[0]\_$c -config $config -tempdir $tempdir -log $output/progress/$pre\_$uID\_$gpair[0]\_$c\_ABRA_WRAPPER.log`;
		
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
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE", job_hold => "$abra_jid", cpu => "4", mem => "50", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$bcount\_$c\_FIXMATE.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "5");
		    my $additionalParams = Schedule::additionalParams(%addParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar FixMateInformation I=$outBam O=$outBam\_FM.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $fmb O=$output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx5g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF_SEQ -L $c $multipleTargets --known $MILLS_1000G --known $DB_SNP -S LENIENT -nt 10 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigner.intervals $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_RTC.done`;
		$rtc_jid = "$pre\_$uID\_$gpair[0]\_$c\_RTC";
		$ran_rtc = 1;
	    }
	    
	    ### seems to randomly have issues with file locking timeout on the m nodes so not submitting to it anymore
	    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.done" || $ran_rtc){
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_IR", job_hold => "$rtc_jid", cpu => "2", mem => "15", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx15g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF_SEQ -L $c $multipleTargets --knownAlleles $MILLS_1000G --knownAlleles $DB_SNP -S LENIENT --targetIntervals $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigner.intervals --maxReadsForRealignment 500000 --maxReadsInMemory 3000000 --maxReadsForConsensuses 500000 -rf BadCigar --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam $bgroup`;
		`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_IR.done`;
		push @ir_jids, "$pre\_$uID\_$gpair[0]\_$c\_IR";
		$ran_ir = 1;
	    }
	}
	else{ #skip indel realignment

	}
	push @indelBams, "-I $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam";
    }

    my $irBams ='';
    my $irj ='';
    if($abra || $indelrealigner){
        $irBams = join(" ", @indelBams);
        $irj = join(",", @ir_jids);
    }
    else{  #skipped indel realignment previously
        $irBams = $bgroup;
    }

    my $ran_br = 0;
    my $brj = '';
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_BR.done" || $ran_ir){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_BR", job_hold => "$irj", cpu => "12", mem => "60", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_BR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx50g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO -R $REF_SEQ -S LENIENT --knownSites $DB_SNP --knownSites $MILLS_1000G --knownSites $HAPMAP --knownSites $OMNI_1000G --knownSites $PHASE1_SNPS_1000G --covariate ContextCovariate --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -rf BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_$gpair[0]\_recal_data.grp $irBams`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_BR.done`;
	$brj = "$pre\_$uID\_$gpair[0]\_BR";
	$ran_br = 1;
    }
    
    my @indelRecalBams1 = ();
    my $ran_pr = 0;
    my @pr_jids = ();
    foreach my $c (@ref_chr){
	if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_PR.done" || $ran_br){
            my $prBams = '';
            if($abra || $indelrealigner){
                $prBams = "-I $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned.bam";
            }
            else{  #skipped indel realignment previously
                $prBams = $bgroup;
            }          
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_$c\_PR", job_hold => "$brj", cpu => "6", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_$c\_PR.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -R $REF_SEQ -L $c $multipleTargets --emit_original_quals -BQSR $output/intFiles/$pre\_$gpair[0]\_recal_data.grp --num_cpu_threads_per_data_thread 6 -rf BadCigar --downsampling_type NONE --out $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned_recal.bam $prBams`;
	    `/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_$c\_PR.done`;
	    push @pr_jids, "$pre\_$uID\_$gpair[0]\_$c\_PR";
	    #push @prg_jids, "$pre\_$uID\_$gpair[0]\_$c\_PR";    ### HC jobs sometimes failed to submit due to holding on too many PR jobs, so we now make a dummy holding job per group, and let HC holdiing on group dummy jobs
	    $ran_pr = 1;
	    $ran_pr_glob{"$c"} = 1;
	}
    
	push @indelRecalBams1, "I=$output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned_recal.bam";
	push @{$processedBams{"$c"}}, "-I $output/intFiles/$pre\_$gpair[0]\_$c\_indelRealigned_recal.bam";
    }


### HC jobs sometimes failed to submit due to holding on too many PR jobs, so we now create a dummy holding job per group, and let HC holdiing on group dummy jobs
    my $prj = join(",", @pr_jids);
    if($ran_pr)
    {
        sleep(2);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_HOLD_PR", job_hold => "$prj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_HOLD_PR.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my %addParams = (scheduler => "$scheduler", runtime => "10", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "0", iounits => "0");
        my $additionalParams = Schedule::additionalParams(%addParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams sleep 1`;
        push @prg_jids, "$pre\_$uID\_$gpair[0]\_HOLD_PR";
    }


    my $irBams1 = join(" ", @indelRecalBams1);
    my $ran_m = 0;
    my @merge_jids = ();
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE_PR.done" || $ran_pr){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_MERGE_PR", job_hold => "$prj", cpu => "8", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_MERGE_PR.log");
	my $standardParams = Schedule::queuing(%stdParams);
	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "10");
	my $additionalParams = Schedule::additionalParams(%addParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $irBams1 O=$output/intFiles/$pre\_$gpair[0]\_indelRealigned_recal.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
	`/bin/touch $output/progress/$pre\_$uID\_$gpair[0]\_MERGE_PR.done`;
	push @merge_jids, "$pre\_$uID\_$gpair[0]\_MERGE_PR";
	$ran_m = 1;
    }

    my $mj = join(",", @merge_jids);
    if(!-e "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.done" || $ran_m){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$gpair[0]\_SSF", job_hold => "$mj", cpu => "2", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_$gpair[0]\_SSF.log");
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
my @ad_jids = ();
my $ran_mqm = 0;
my $ran_ad = 0;
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

    ## generate *alleledepth file for 'hotspots in normals' qc
    if(!-e "$output/progress/$pre\_$uID\_MUT_AD_$samp\.done" || $ran_ssf){
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MUT_AD_$samp", job_hold => "$ssfj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_MUT_AD_$samp\.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/qc/dmp_genotype_allele.pl -fmv $COSMIC_HOTSPOTS -bam $finalBam -rf $REF_SEQ -s $SAMTOOLS -b $BEDTOOLS -o $output/intFiles/$samp -of $pre\_indelRealigned_recal_$samp\.mpileup.alleledepth -mof $pre\_indelRealigned_recal_$samp\.mpileup`;
        push @ad_jids, "$pre\_$uID\_MUT_AD_$samp";
        `/bin/touch $output/progress/$pre\_$uID\_MUT_AD_$samp\.done`;
        $ran_ad = 1;
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

my $adj = join(",", @ad_jids);
my $ran_hn = 0;
if(!-e "$output/progress/$pre\_$uID\_HOTSPOTS_NORMALS.done" || $ran_ad){
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HOTSPOTS_NORMALS", job_hold => "$adj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_HOTSPOTS_NORMALS.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/hotspots_in_normals.py $output/intFiles '*.alleledepth' $COSMIC_HOTSPOTS $pair $output/metrics/$pre\_HotspotsInNormals.txt`;
    `/bin/touch $output/progress/$pre\_$uID\_HOTSPOTS_NORMALS.done`;
    push @all_jids, "$pre\_$uID\_HOTSPOTS_NORMALS";
    $ran_hn = 1;
}

my $allj = join(",", @all_jids);
#if(!-e "$output/progress/$pre\_$uID\_RSYNC_1.done" || $ran_ssf || $ran_mmqm || $ran_hn){
sleep(2);
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC_1", job_hold => "$allj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_1.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' --exclude 'variants' --exclude 'metrics' --exclude 'variants_pipeline' --exclude 'rna' $curDir $rsync`;
push @all_jids, "$pre\_$uID\_RSYNC_1";
`/bin/touch $output/progress/$pre\_$uID\_RSYNC_1.done`;
#}

if($nosnps){
    exit(0);
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
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx24g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_SEQ --reference_sample_name $species -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotateNDA --annotation AlleleBalance --annotation AlleleBalanceBySample --annotation HardyWeinberg --genotype_likelihoods_model BOTH --read_filter BadCigar --num_cpu_threads_per_data_thread 12 --out $output/intFiles/$pre\_$c\_UnifiedGenotyper.vcf $irBams2`;
	    `/bin/touch $output/progress/$pre\_$uID\_$c\_UG.done`;
	    push @ug_jids, "$pre\_$uID\_$c\_UG";
	    $ran_ug = 1;
	}
    }

    ### NOTE: ANNOTATIONS THAT DON'T WORK:AlleleBalance, HardyWeinberg,IndelType
    if(!-e "$output/progress/$pre\_$uID\_$c\_HC.done" || $ran_pr_glob{"$c"}){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$c\_HC", job_hold => "$prgj", cpu => "36", mem => "180", cluster_out => "$output/progress/$pre\_$uID\_$c\_HC.log");
	my $standardParams = Schedule::queuing(%stdParams);
	my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "7");
	my $additionalParams = Schedule::additionalParams(%addParams);
        my $rna_param = '';
        if($rna){
            $rna_param = "-dontUseSoftClippedBases -stand_call_conf 20.0";
        }
	my $response = "";
    my $HC_submission_num = 0;
    while($response !~ /is submitted to default queue/ && $HC_submission_num < 10){
        my $time = 60 * $HC_submission_num;
        `sleep $time `;
        $HC_submission_num++;
        my $stamp = localtime(time);
        print "$stamp : HC Chromosome $c attempt number $HC_submission_num.\n"; 
        print "Command: $standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx180g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 30 --out $output/intFiles/$pre\_$c\_HaplotypeCaller.vcf $irBams2\n";
        $response = `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx180g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF_SEQ -L $c $multipleTargets --dbsnp $DB_SNP --downsampling_type NONE --annotation AlleleBalanceBySample --annotation ClippingRankSumTest --read_filter BadCigar --num_cpu_threads_per_data_thread 30 $rna_param --out $output/intFiles/$pre\_$c\_HaplotypeCaller.vcf $irBams2        2>&1`;
        print "Response: $response\n";
        if($response !~ /is submitted to default queue/ && $HC_submission_num == 10 && $c eq '1'){
            `mail -s "HaplotypeCaller Job Not Working" $email <<< "This is an e-mail letting you know that your variants pipeline run for $pre is not submitting HC jobs to the cluster. Please remove all .done files for haplotype caller and rerun the pipeline."  `;
        }
    }
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
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_HC", job_hold => "$hcj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_CV_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_HaplotypeCaller_RAW.vcf --assumeIdenticalSamples $hcVars`;
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
if((!-e "$output/progress/$pre\_$uID\_VR_SNP_HC.done" || $ran_cv_hc) && !$rna){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_SNP_HC", job_hold => "$cvhcj", cpu => "4", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_VR_SNP_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1000G -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE1_SNPS_1000G -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DB_SNP -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_SNP.plots.R -nt 4`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_SNP_HC.done`;
    $vrshcj = "$pre\_$uID\_VR_SNP_HC";
    $ran_vr_snp_hc = 1;
}

### NOTE: sometimes throws an error when running with multiple threads
###       about not being able to find a tmp file; running with -nt 1 to avoid errors
my $ran_ar_snp_hc = 0;
my $arshcj = '';
if((!-e "$output/progress/$pre\_$uID\_AR_SNP_HC.done" || $ran_vr_snp_hc) && !$rna){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_SNP_HC", job_hold => "$vrshcj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_AR_SNP_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_SNP.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_SNP.recal -o $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_SNP_HC.done`;
    $arshcj = "$pre\_$uID\_AR_SNP_HC";
    $ran_ar_snp_hc = 1;
}

my $ran_vr_indel_hc = 0;
my $vrihcj = '';
if((!-e "$output/progress/$pre\_$uID\_VR_INDEL_HC.done" || $ran_ar_snp_hc) && !$rna){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_VR_INDEL_HC", job_hold => "$arshcj", cpu => "4", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_VR_INDEL_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS_1000G -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -rscriptFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.plots.R -nt 4`;
    `/bin/touch $output/progress/$pre\_$uID\_VR_INDEL_HC.done`;
    $vrihcj = "$pre\_$uID\_VR_INDEL_HC";
    $ran_vr_indel_hc = 1;   
}

my $ran_ar_indel_hc = 0;
my $arihcj = '';
if((!-e "$output/progress/$pre\_$uID\_AR_INDEL_HC.done" || $ran_vr_indel_hc) && !$rna){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AR_INDEL_HC", job_hold => "$vrihcj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_AR_INDEL_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_HaplotypeCaller_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.tranches -recalFile $output/intFiles/$pre\_HaplotypeCaller_INDEL.recal -o $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf -nt 1`;
    `/bin/touch $output/progress/$pre\_$uID\_AR_INDEL_HC.done`;
    $arihcj = "$pre\_$uID\_AR_INDEL_HC";
    $ran_ar_indel_hc = 1;   
}


#use hard variants filtering for rnaseq
my $run_rna_vf_hc = 0;
my $rvfhcj = '';
if((!-e "$output/progress/$pre\_$uID\_RNA_VF_HC.done" || $ran_cv_hc) && $rna){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RNA_VF_HC", job_hold => "$cvhcj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_RNA_VF_HC.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_SEQ -V $output/intFiles/$pre\_HaplotypeCaller_RAW.vcf -window 35 -cluster 3 -filterName FS -filter \"FS \> 30.0\" -filterName QD -filter \"QD \< 2.0\" -o $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf`;
    `/bin/touch $output/progress/$pre\_$uID\_RNA_VF_HC.done`;
    $rvfhcj = "$pre\_$uID\_RNA_VF_HC";
    $run_rna_vf_hc = 1;   
} 

# This runs regardless because it has a .done check in the function
sleep(2);
if(!$rna){
    &generateMaf("$output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'HaplotypeCaller', "$arihcj,$ssfj", $ran_ar_indel_hc);
}
else{
    &generateMaf("$output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf", 'HaplotypeCaller', "$rvfhcj,$ssfj", $run_rna_vf_hc);
}
if($ug){
    if(!-d "$output/variants/snpsIndels/unifiedgenotyper"){
	mkdir("$output/variants/snpsIndels/unifiedgenotyper", 0775) or die "Can't make $output/variants/snpsIndels/unifiedgenotyper";
    }
    my $ugVars = join(" " , @ugVariants);
    my $ran_cv_ug = 0;
    my $cvugj = '';
    my $ugj = join(",", @ug_jids);
    if(!-e "$output/progress/$pre\_$uID\_CV_UG.done" || $ran_ug){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CV_UG", job_hold => "$ugj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_CV_UG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R $REF_SEQ -o $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf --assumeIdenticalSamples $ugVars`;
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1000G -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE1_SNPS_1000G -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DB_SNP -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.plots.R -nt 4`;
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_RAW.vcf --ts_filter_level 99.0 -mode SNP -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_SNP.recal -o $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -nt 1`;
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS_1000G -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -rscriptFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.plots.R -nt 4`;
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
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_SEQ -input $output/intFiles/$pre\_UnifiedGenotyper_SNP_vqsr.vcf --ts_filter_level 99.0 -mode INDEL -tranchesFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.tranches -recalFile $output/intFiles/$pre\_UnifiedGenotyper_INDEL.recal -o $output/variants/snpsIndels/unifiedgenotyper/$pre\_UnifiedGenotyper.vcf -nt 1`;
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

    if(!-d "$output/variants/snpsIndels/haplotect"){
	mkdir("$output/variants/snpsIndels/haplotect", 0775) or die "Can't make $output/variants/snpsIndels/haplotect";
    }
    if(!-d "$output/variants/copyNumber"){
	mkdir("$output/variants/copyNumber", 0775) or die "Can't make $output/variants/copyNumber";
    }
    if(!-d "$output/variants/copyNumber/facets"){
	mkdir("$output/variants/copyNumber/facets", 0775) or die "Can't make $output/variants/copyNumber/facets";
    }
    if(!-d "$output/variants/structVar"){
	mkdir("$output/variants/structVar", 0775) or die "Can't make $output/variants/structVar";
    }

    `/bin/echo "Tumor_Sample_Barcode\tRdata_filename" > $output/variants/copyNumber/facets/facets_mapping.txt `;
    ###mkdir("$output/variants/varscan", 0775) or die "Can't make $output/variants/varscan";

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
	if(!-d "$output/variants/snpsIndels/mutect"){
	    mkdir("$output/variants/snpsIndels/mutect", 0775) or die "Can't make $output/variants/snpsIndels/mutect";
	}
	my $ran_mutect = 0;
	my $mutectj = '';
	### MUTECT will fail if order of ref seq contigs and vcf don't match
	if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.done" || $ran_ssf){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT", job_hold => "$ssfj", cpu => "2", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA7_MUTECT/java -Xmx16g -Djava.io.tmpdir=$tempdir -jar $MUTECT/muTect.jar --analysis_type MuTect --reference_sequence $REF_SEQ --dbsnp $DB_SNP --cosmic $COSMIC --input_file:normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --input_file:tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --vcf $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.vcf --out $output/variants/snpsIndels/mutect/$pre\_$data[0]\_$data[1]\_mutect_calls.txt -rf BadCigar --enable_extended_output --downsampling_type NONE`;
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
	    ### my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF", job_hold => "$ssj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_SOMATIC_SNIPER_MAF.log");
	    ###my $standardParams = Schedule::queuing(%stdParams);
	    ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/maf/vcf2maf0.py -i $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper.vcf -c somaticsniper -o $output/variants/snpsIndels/somaticsniper/$pre\_indelRealigned_recal\_$data[0]\_$data[1]\_somatic_sniper_MAF.txt -n $data[0] -t $data[1]`;
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
	    ###`/common/sge/bin/lx24-amd64/qsub -N $pre\_$uID\_$data[0]\_$data[1]\_VARSCAN_SOMATIC -hold_jid $pre\_$uID\_$data[0]\_$data[1]\_MPILEUP -pe alloc 2 -l virtual_free=5G -q lau.q,lcg.q,nce.q $Bin/qCMD $JAVA/java -Xms256m -Xmx10g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $VARSCAN/VarScan.jar somatic $output/intFiles/$pre\_indelRealigned_recal\_$data[0]\.bam.mpileup $output/intFiles/$pre\_indelRealigned_recal\_$data[1]\.bam.mpileup $output/variants/varscan/$pre\_$data[0]\_$data[1]\_varscan_somatic --strand-filter 1 --output-vcf 1`;
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
		    ### strelka DIES IF DIR ALREADY EXISTS
		    `/bin/rm -rf $output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka`;
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
		### NOTE: strelka ONLY HAS CONFIG FOR BWA ALN, NOT SURE HOW IT WILL WORK WITH BWA MEM
		sleep(2);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG", job_hold => "$lnsj", cpu => "1", mem => "2", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_STRELKA_CONFIG.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $STRELKA/bin/configureStrelkaWorkflow.pl --normal=$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --tumor=$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --ref=$REF_SEQ --config=$STRELKA/etc/strelka_config_bwa_default.ini --output-dir=$output/variants/snpsIndels/strelka/$data[0]\_$data[1]\_strelka`;
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

	    ### NOTE: THERE ARE REGIONS OF THE ASSEMBLY THAT HAVE NON ATGC BASES, SUCH AS M ON CHR3
	    ###       OTHER CALLERS SEEM TO SKIP THESE BASES; LANCET DOES NTO AND MAKES CALLS
	    ###       THESE VCF FILES FAIL TO MERGE IN GATK BECAUSE IT IS CONSIDERED MALFORMED WITH AN UNPARSABLE VCF WITH THAT ALLELE

	    if(!-d "$output/variants/snpsIndels/lancet"){
		mkdir("$output/variants/snpsIndels/lancet", 0775) or die "Can't make $output/variants/snpsIndels/lancet";
	    }

	    #### LANCET ON TARGETED REGION BED
	    #my $ran_lancet = 0;
	    #if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_LANCET.done" || $ran_ssf){
	#	sleep(2);
		
	#	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_LANCET", job_hold => "$ssfj", cpu => "12", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_LANCET.log");
	#	my $standardParams = Schedule::queuing(%stdParams);
	#	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $LANCET/lancet --tumor $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam --normal $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam --ref $REF_SEQ --bed $targets_bed_padded --num-threads 12 ">$output/variants/snpsIndels/lancet/$pre\_$data[0]\_$data[1]\_lancet_calls.vcf"`;
	#	`/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_LANCET.done`;
	#	push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_LANCET";
	#	$ran_lancet = 1;
	 #   }
	    
	    
	    #### END LANCET ON TARGETED REGION
	    
	    
	    ### LANCET WHOLE GENOME SPLIT BY CHR

	    ### NOTE: human assemblies have non ATGC bases
	    ###       lancet makes calls at those positions; haven't found a case where other callers do
	    ###       these alleles cause gatk CombineVariants to fail because it thinks it's a malformed vcf
	    ###       added step to remove those positions from the vcf

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
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_CV_LANCET", job_hold => "$flj", cpu => "1", mem => "20", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_CV_LANCET.log");
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
		
		
		#`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /opt/common/CentOS_6/VarDictJava/VarDict-1.5.1/bin/VarDict -G $REF_SEQ -f 0.01 -N $data[1] -b "$output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam|$output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam" -z -F 0 -c 1 -S 2 -E 3 -g 4 $targets_bed_padded | /opt/common/CentOS_6/VarDict/VarDict_85cc3f6/testsomatic.R | /opt/common/CentOS_6/VarDict/VarDict_85cc3f6/var2vcf_paired.pl -N "$data[1]|$data[0]" -f 0.01 ">$output/variants/snpsIndels/vardict/$pre\_$data[0]\_$data[1]\_vardict_calls.vcf"`;
		
		
		### NOTE: need to write out command to file and submit command
		###       beause i can't figure out how to get double quotes 
		###       to be included in bsub command; 
		###       it always gets removed, no matter how i try to escape it with backslash
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

                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PINDEL/pindel -i $output/intFiles/$data[0]\_$data[1]\_pindel_config.txt -f $REF_SEQ -c ALL -o $output/intFiles/pindel/$data[0]\_$data[1]/$data[0]\_$data[1] -r true -t true -I true -T 12`;

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

        if($mutect2){
            if(!-d "$output/variants/snpsIndels/mutect2"){
                mkdir("$output/variants/snpsIndels/mutect2", 0775) or die "Can't make $output/variants/snpsIndels/mutect2";
            }
            if(!-e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT2.done" || $ran_ssf){
                sleep(2);
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT2", job_hold => "$ssfj", cpu => "8", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT2.log"); 
                my $standardParams = Schedule::queuing(%stdParams);

                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $JAVA/java -Xmx25g -Djava.io.tmpdir=$tempdir -jar $GATK4/gatk-package-4.1.1.0-local.jar Mutect2 --native-pair-hmm-threads 4 -R $REF_SEQ -I $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam -I $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam -O $output/variants/snpsIndels/mutect2/$pre\_$data[0]\_$data[1]\_mutect2_calls.vcf -normal $data[0] -tumor $data[1] --max-reads-per-alignment-start 0 --read-filter GoodCigarReadFilter --read-filter ReadLengthEqualsCigarLengthReadFilter`;

                `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_MUTECT2.done`;
                push @all_jids, "$pre\_$uID\_$data[0]\_$data[1]\_MUTECT2";
            }
        }

        ## Here we will add the facets scripts
        ## Set up tumor and normal counts
        my $facetsSETUP_jid = '';
        my $facets_setup = 0;
        if(!$nofacets && $hasPair && (! -e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.done" || $ssfj )) {
            if(-d "$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets"){
                `/bin/rm -rf -d "$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets"`;
            }
            mkdir("$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets", 0775) or die "Can't make $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets";
            mkdir("$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp", 0775) or die "Can't make $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp";

            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP",  cpu => "4", mem => "4", job_hold => "$ssfj", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.log");
            my $standardParams = Schedule::queuing(%stdParams);
            my %addParams_local = (scheduler => "$scheduler", runtime => "30", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
            my $additionalParams_local = Schedule::additionalParams(%addParams_local);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams_local $singularityParams $HTSTOOLS/snp-pileup -A -g -P 50 -r15,0 $FACETS_DB_SNP $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_countsMerged_$data[0]\_$data[1].dat.gz $output/alignments/$pre\_indelRealigned_recal\_$data[0]\.bam $output/alignments/$pre\_indelRealigned_recal\_$data[1]\.bam`;
	    
            $facets_setup = 1;
            `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP.done`;
            $facetsSETUP_jid = "$pre\_$uID\_$data[0]\_$data[1]\_facets_SETUP";
        }
        ## now facets
        if(!$nofacets && $hasPair && (! -e "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.done" || $facets_setup)){ 
            #parameters for hisens
            my $cval = $impact ? 100 : 150; 
            my $snp_nbhd = $impact ? 250 : 150;
            my $min_nhet = $impact ? 10 : 15;
            my $ndepth = 15;          
 
            # parameters for purity
            my $p_cval = $impact ? 150 : 250;
            my $p_snp_nbhd = $impact ? 250 : 150;
            my $p_min_nhet = $impact ? 15 : 25;
            my $p_ndepth = 15;


            my $sub_species = $species eq 'b37' ? 'hg19' : $species;

            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN",  cpu => "3", mem => "60", job_hold => "$facetsSETUP_jid", cluster_out => "$output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.log");
            my $standardParams = Schedule::queuing(%stdParams);
            my %addParams = (scheduler => "$scheduler", runtime => "10", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1");
            my $additionalParams = Schedule::additionalParams(%addParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $FACETS_SUITE/facets doFacets -c $cval -s $snp_nbhd -n $ndepth -m $min_nhet -pc $p_cval -ps $p_snp_nbhd -pn $p_ndepth -pm $p_min_nhet -f $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/tmp/$pre\_countsMerged_$data[0]\_$data[1].dat.gz -t $data[0]\_$data[1] -D $output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets -r $FACETS_LIB`;
	    push @facets_jid, "$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN" ;
            $facets_run = 1;
            `/bin/touch $output/progress/$pre\_$uID\_$data[0]\_$data[1]\_facets_RUN.done`; 
        }
    `/bin/echo "$data[1]\t$output/variants/copyNumber/facets/$data[0]\_$data[1]\_facets/$data[0]\_$data[1]\_hisens.Rdata" >> $output/variants/copyNumber/facets/facets_mapping.txt`;
    }
    close PAIR;

    ## run merge seg script
    my $facets_haplotect_jid = ''; 
    if(!$nofacets && $hasPair && (! -e "$output/progress/$pre\_$uID\_merge_facets_seg.done" || $facets_run)){
        my $seg_outfile = "$output/variants/copyNumber/facets/$pre\_facets_merge_hisens.seg";
        if( -f "$seg_outfile"){
            unlink("$seg_outfile") or die "Cannot delete? $!";
        }
        my $facets_js = join(",", @facets_jid);
	
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_merge_facets_seg", cpu => "1", mem => "1", job_hold => "$facets_js", cluster_out => "$output/progress/$pre\_$uID\_merge_facets_seg.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my %addParams = (scheduler => "$scheduler", runtime => "1", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
        my $additionalParams = Schedule::additionalParams(%addParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/facets/merge_facets_seg.pl -facets_dir $output/variants/copyNumber/facets -outfile $seg_outfile`;
	
        `/bin/touch $output/progress/$pre\_$uID\_merge_facets_seg.done`;
        $facets_haplotect_jid = "$pre\_$uID\_merge_facets_seg";
        $facets_run = 1;
    }

    ### create cna file for portal upload
    ### NOTE: the ___HISENS_GeneCalls_v2.txt chr field doesn't have chr even for hg19
    ###       so the target has to be b37 chr convention even if input is hg19 chr convention with chr
    my $facets_geneLevel_jid = '';
    my $ran_fgl = 0;
    if(!$nofacets && $hasPair && (! -e "$output/progress/$pre\_$uID\_facets_genelevel.done" || $facets_run)){
        my $facets_js = join(",", @facets_jid);
	
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_facets_genelevel", cpu => "1", mem => "10", job_hold => "$facets_js", cluster_out => "$output/progress/$pre\_$uID\_facets_genelevel.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my %addParams = (scheduler => "$scheduler", runtime => "1", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
        my $additionalParams = Schedule::additionalParams(%addParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $FACETS_SUITE/facets geneLevel -f $output/variants/copyNumber/facets/*/*_hisens.cncf.txt -t $targets_facet -o $output/intFiles/$pre\___HISENS_GeneCalls_v2.txt`;
	
        `/bin/touch $output/progress/$pre\_$uID\_facets_genelevel.done`;
        $facets_geneLevel_jid = "$pre\_$uID\_facets_genelevel";
        $ran_fgl = 1;
    }

    if(!$nofacets && $hasPair && (! -e "$output/progress/$pre\_$uID\_facets_split_gene.done" || $ran_fgl)){	
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_facets_split_gene", cpu => "1", mem => "10", job_hold => "$facets_geneLevel_jid", cluster_out => "$output/progress/$pre\_$uID\_facets_split_gene.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my %addParams = (scheduler => "$scheduler", runtime => "1", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
        my $additionalParams = Schedule::additionalParams(%addParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{job_hold} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/facets/split_gene_file.py -g $output/intFiles/$pre\___HISENS_GeneCalls_v2.txt -p $pair -d $output/variants/copyNumber/facets/`;
	
        `/bin/touch $output/progress/$pre\_$uID\_facets_split_gene.done`;
        push @all_jids, "$pre\_$uID\_facets_split_gene";
    }

    if($hasPair && (!-e "$output/progress/$pre\_$uID\_HAPLOTECT.done" || $ran_mutect_glob || $ran_ar_indel_hc)){
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
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HAPLOTECT", job_hold => "$arihcj,$muj", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_HAPLOTECT.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/haploTect_merge.pl -pair $pair -hc_vcf $output/variants/snpsIndels/haplotypecaller/$pre\_HaplotypeCaller.vcf -species $species -pre $pre -output $output/variants/snpsIndels/haplotect -mutect_dir $output/variants/snpsIndels/mutect -config $config $patientFile -align_dir $output/alignments/ -svnRev $svnRev $addOptions -tempdir $tempdir -delete_temp`;

        $haplotect_run = 1;
	`/bin/touch $output/progress/$pre\_$uID\_HAPLOTECT.done`;
        $facets_haplotect_jid .= ",$pre\_$uID\_HAPLOTECT";
    }

    ## Now join maf
    my $mafAnnoRun=0;
    if(!$nofacets && $hasPair && (!-e "$output/progress/$pre\_$uID\_join_maf.done" || $haplotect_run || $facets_run)){
        sleep(2);

       my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_join_maf", job_hold => "$facets_haplotect_jid", cpu => "4", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_join_maf.log");
        my $standardParams = Schedule::queuing(%stdParams); 
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $FACETS_SUITE/facets mafAnno -m $output/variants/snpsIndels/haplotect/$pre\_haplotect_VEP_MAF.txt -f $output/variants/copyNumber/facets/facets_mapping.txt -o $output/intFiles/$pre\_CMO_MAF_intermediate.txt`; 
        `/bin/touch $output/progress/$pre\_$uID\_join_maf.done`;
	push @all_jids, "$pre\_$uID\_join_maf";
        $mafAnnoRun = 1;
    }

    if(!$nofacets && $hasPair && (! -e "$output/progress/$pre\_$uID\_wes_filters.done" || $mafAnnoRun )) {
        my $holdVar = '';
        if( $mafAnnoRun) { 
            $holdVar = "$pre\_$uID\_join_maf";
        }
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_wes_filters", job_hold => "$holdVar", cpu => "2", mem => "8", cluster_out => "$output/progress/$pre\_$uID\_wes_filters.log");
        my $standardParams = Schedule::queuing(%stdParams);
        my %addParams = (scheduler => "$scheduler", runtime => "7", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "0");
        my $additionalParams = Schedule::additionalParams(%addParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $Bin/maf/post_filters.pl -in_maf $output/intFiles/$pre\_CMO_MAF_intermediate.txt -out_maf $output/variants/$pre\_CMO_MAF.txt -config $config -species $species -filter_ffpe -blacklist -low_conf`;
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
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/exome_cnv.pl -pre $pre -result $output/variants/copyNumber/dmp_cnv -berger $target_design -bamlist $output/intFiles/$pre\_sv_bam_list.txt -patient $patient -std_covg $target_std_normals -genome $species -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;
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
        #my $tempGenome = $species == "hybrid" ? "b37" : $species;
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/RunStructuralVariantPipeline_Delly.pl -pre $pre -out $output/variants/structVar/delly -pair $pair -bam_list $output/intFiles/$pre\_sv_bam_list.txt -genome $species -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group`;
	`/bin/touch $output/progress/$pre\_$uID\_DELLY.done`;
	$ran_strvar = 1;
    }
    
    if($hasPair && (!-e "$output/progress/$pre\_$uID\_CDNA_CONTAM.done" || $ran_strvar)){
	sleep(2);

        my $hold_value = $ran_strvar ? "$pre\_$uID\_DELLY" : "";
	
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CDNA_CONTAM", job_hold => "$hold_value", cpu => "1", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_CDNA_CONTAM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PYTHON/python $Bin/qc/check_cDNA_contamination.py -s $output/variants/structVar/delly/$pre\_AllAnnotatedSVs.txt -o $output/metrics/$pre\_cDNA_contamination.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_CDNA_CONTAM.done`;
	push @all_jids, "$pre\_$uID\_CDNA_CONTAM";
    }
}

# FOR TESTING
#die;

my $allj2 = join(",", @all_jids);
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC_2", job_hold => "$allj2", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_2.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' --exclude 'variants_pipeline' --exclude 'rna' $curDir $rsync`;
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
        next if($c =~ /^MM/);
        if(!-e "$output/progress/$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED.done" || $ran_hc){
            if((! -e "$vcf.gz" || $ran_hc) && !$bgzipped){
                $bgzipped = 1;
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_bgzip", job_hold => "$hold", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_bgzip.log");
                my $standardParams = Schedule::queuing(%stdParams);
                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$singularityParams $TABIX/bgzip -cf $vcf > $vcf.gz"`;
                %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_bgzip_index", job_hold => "$pre\_$uID\_$jna\_bgzip", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_bgzip_index.log");
                $standardParams = Schedule::queuing(%stdParams);
                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $BCFTOOLS/bcftools index $vcf.gz`;
                $bgz_jid = "$pre\_$uID\_$jna\_bgzip,$pre\_$uID\_$jna\_bgzip_index";
            }

            my $addOptions = "";
            if($ExAC_VCF){
                $addOptions = "-exac_vcf $ExAC_VCF";
            }
	    
	    if(!-d "$vcf_dir/chrom_$c"){
		mkdir("$vcf_dir/chrom_$c", 0775) or die "Can't make $vcf_dir/chrom_$c";
	    }
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_split_CHR_$c", job_hold => $bgz_jid, cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_split_$c.log");
            my $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $BCFTOOLS/bcftools filter -r $c $vcf.gz -O v -o $vcf_dir/chrom_$c/$jna\_$c.vcf`;

            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED", job_hold => "$pre\_$uID\_$jna\_split_CHR_$c", cpu => "4", mem => "60", cluster_out => "$output/progress/$pre\_$uID\_$jna\_$c\_MAF_UNPAIRED.log");
            $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $Bin/generateMAF.pl -vcf $vcf_dir/chrom_$c/$jna\_$c.vcf -species $species -config $config -caller $type $patientFile -align_dir $output/alignments $addOptions -delete_temp`;
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
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_MAF";

        $merge_files =~ s/UNPAIRED_TCGA_MAF.txt/UNPAIRED_TCGA_PORTAL_MAF.txt/g;
        %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF.log");
        $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF";

        $merge_files =~ s/UNPAIRED_TCGA_PORTAL_MAF.txt/UNPAIRED_VEP_MAF.txt/g;
        %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_VEP_MAF", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_VEP_MAF.log");
        $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_VEP_MAF.txt`;
        push @merge_jids, "$pre\_$uID\_$jna\_merge_VEP_MAF";

        if($patient){
            s/vcf_UNPAIRED_TCGA_MAF.txt$/vcf_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt/g for @merge_files;
            $merge_files = join(" -i ", @merge_files);
            %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF_fillout", job_hold => "$jid_holds", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF_fillout.log");
            $standardParams = Schedule::queuing(%stdParams);
            ###print "$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt\n";
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams $PERL/perl $Bin/maf/mergeMaf.pl -i $merge_files -o $vcf_dir/$jna\_UNPAIRED_TCGA_PORTAL_MAF_fillout.txt`;
            push @merge_jids, "$pre\_$uID\_$jna\_merge_TCGA_PORTAL_MAF_fillout";
        }

        #Now to clean up!
        my $merge_holds = join(",", @merge_jids);
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$jna\_cleanup", job_hold => "$merge_holds", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_$jna\_cleanup.log");
        $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $singularityParams rm -rf $vcf_dir/chrom_* $vcf.gz $vcf.gz.csi`;
        push @all_jids, "$pre\_$uID\_$jna\_cleanup";
    }
}
