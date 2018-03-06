#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Schedule;
use Cluster;

### takes a list of directories containing fastq.gz files
### 1) makes directories for each sample
### 2) soft lins to all fastq.gz files in directory
### 3) runs solexa_HG19_PE.pl

### ex. input list $file
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14028 PE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14278 PE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14314 PE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14348 SE
### LIB_ID SAMPLE_ID RUN_id /ifs/archive/GCL/hiseq/FASTQ/LIZ_0096_AB07L5ABXX/Project_2144/Sample_14512 SE
### etc.

### ex. input list $group
###     how to group samples for realign/recal steps
###     s_WD1345  Group_1
###     s_WD1346  Group_1
###     s_WD1349  Group_1
###     s_NF1341  Group_2
###     s_WD1340  Group_2
###     s_WD1342  Group_2
###     s_NF1311  Group_3
###     s_WD1306  Group_3

### ex. input list $pair
###     tumor/normal pairing for mutect/maf conversion
###     NOTE: if no pairing file provided, will assume no somatic analysis
###     s_NF1236        s_WD1216
###     s_NF1236        s_DD1278
###     s_WD1216        s_DD1278
###     s_NF1236        s_WD1216
###     s_NF1236        s_WD1210

### nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed

### config file: list of programs
### PROGRAM_NAME PATH_TO_PROGRAM
### ex.
### PICARD /ifs/data/mpirun/bin/picard
### GATK /ifs/data/mpirun/bin/gatk/GenomeAnalysisTK.jar


###  CONFIG FILE MUST CONTAIN PATHS TO DIRECTORIES THE FOLLOWING PROGRAMS; DEFAULTS TO WHAT IS IN /opt/bin IF NOT SPECIFIED
### BWA, PICARD, GATK, MUTECT, BEDTOOLS, ONCOTATOR, CUTADAPT
### NOTE: MUST GIVE FULL PATH TO CONFIG FILE


### NOTE: SAMPLE NAMES MUST BE CONSISTENT
###       THIS MAY BE A PROBLEM WHEN SAMPLE NAME STARTS WITH DIGIT BECAUSE 
###       IN SCRIPT IT WILL GET PREPEPENDED WITH s_; i.e. sample 34fsh will change to s_34fsh

### NOTE: ASSUMES MARK DUP BAMS ARE IN $sample/$lib/$sample\_$lib\_MD.bam
###       OTHERWISE, 2ND HALF OF PIPELINE WILL FAIL
###       CAN MANUALLY CREATE LIST OF BAMS TO GROUP AND RERUN THAT PORTION IF NEEDED

### NOTE: DO NOT SUBMIT THIS WRAPPER SCRIPT TO THE CLUSTER BECAUSE ONCOTATOR STEP WILL FAIL
###       DUE TO NODES NOT HAVING NETWORK ACCESS

my ($map, $group, $pair, $patient, $impact, $wes, $config, $help, $nosnps, $removedups, $species, $ug, $scheduler, $abra, $indelrealigner, $targets, $mdOnly, $noMD, $DB_SNP, $noClip, $request, $allSomatic, $scalpel, $somaticsniper, $strelka, $varscan, $virmid, $chip, $lancet, $vardict, $pindel);

my $pre = 'TEMP';
my $output = "results";
my $priority_project = "ngs";
my $priority_group = "Pipeline";
my $r1adaptor = 'AGATCGGAAGAGCACACGTCT';
my $r2adaptor = 'AGATCGGAAGAGCGTCGTGTA';
my $bqTrim = '3';

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my $email = "$uID\@cbio.mskcc.org";
my $rsync = "/ifs/solres/$uID";
my $tempdir = "/scratch/$uID";

GetOptions ('map=s' => \$map,
	    'group=s' => \$group,
	    'pair=s' => \$pair,
	    'patient=s' => \$patient,
	    'pre=s' => \$pre,
	    'config=s' => \$config,
	    'targets=s' => \$targets,
	    'request=s' => \$request,
	    'species=s' => \$species,
 	    'scheduler=s' => \$scheduler,
	    'help' => \$help,
	    'nosnps' => \$nosnps,
	    'removedups' => \$removedups,
	    'abra' => \$abra,
	    'impact' => \$impact,
            'wes' => \$wes,
            'indelrealigner' => \$indelrealigner,
	    'mdOnly|mdonly' => \$mdOnly,
	    'chip|chipseq|chip-seq' => \$chip,
	    'noMD|nomd|nomarkdups|noMarkdups' => \$noMD,
	    'ug|unifiedgenotyper' => \$ug,
 	    'output|out|o=s' => \$output,
	    'rsync=s' => \$rsync,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
	    'r1adaptor=s' => \$r1adaptor,
	    'r2adaptor=s' => \$r2adaptor,
	    'noClip' => \$noClip,
	    'bqTrim=s' => \$bqTrim,
 	    'db_snp|dbsnp=s' => \$DB_SNP,
	    'allsomatic|allSomatic|all_somatic' => \$allSomatic,
	    'scalpel' => \$scalpel,
	    'somaticsniper' => \$somaticsniper,
	    'strelka' => \$strelka,
	    'varscan' => \$varscan,
	    'virmid' => \$virmid,
	    'lancet' => \$lancet,
	    'vardict' => \$vardict,
            'pindel' => \$pindel,
	    'email' => \$email,
	    'tempdir=s' => \$tempdir) or exit(1);

if(!$map || !$species || !$config || !$scheduler || !$request || $help){
    print <<HELP;

    USAGE: variants_pipeline.pl -wes -config CONFIG -species SPECIES -scheduler SCHEDULER
	* MAP: file listing sample information for processing (REQUIRED)
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED, unless using -mdOnly flag)
	* SPECIES: b37 (default: human), mm9, mm10 (default: mouse), hybrid (b37+mm10), mm10_custom, species_custom and dm3 currently supported (REQUIRED)
	* TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS unless given full path to targets directory; required for non-chipseq projects
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
        * REQUEST: file containing request information such as PI, investigator, etc. (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* EMAIL: email to send notication of finished final job of pipeline (default: $uID\@cbio.mskcc.org)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
        * PATIENT: if a patient file is given, patient wide fillout will be added to maf file
	* PRE: output prefix (default: TEMP)
	* OUTPUT: output results directory (default: results)
	* RSYNC:  path to rsync data for archive (default: /ifs/solres/$uID)
	* TEMPDIR:  temp directory (default: /scratch/$uID)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* -wes: run pipeline as whole exome (must use either -wes or -impact)
	* -impact: run pipeline as an impact project (must use either -wes or -impact) 
	* -nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed
	* -removedups: remove duplicate reads instead of just marking them
 	* -chip: will stop after markdups step
	* -abra: run abra to realign indels (default)
	* -indelrealigner: run GATK indelrealigner (abra runs as default)
	* -mdOnly: will stop after markdups step (e.g. chip seq analysis)
	* -chip: will stop after markdups step; will also run bwa without -PM options
	* -noMD: will not run MarkDups
	* R1ADAPTOR to specify R1 adaptor sequence (default: AGATCGGAAGAGCACACGTCT)
	* R2ADAPTOR to specify R1 adaptor sequence (default: AGATCGGAAGAGCGTCGTGTA)
	* NOCLIP: no clipping of adaptor sequences
	* BASEQTRIM: base quality to trim in reads (default: 3)
	* DB_SNP: VCF file containing known sites of genetic variation (REQUIRED for species_custom; either here or in config)
	* haplotypecaller is default; -ug || -unifiedgenotyper to also make unifiedgenotyper variant calls	
	* ALLSOMATIC: run all somatic callers; mutect/haplotypecaller always run; otherwise -scalpel, -somaticsniper, -strelka, -varscan, -virmid, -lancet, -vardict, -pindel to run them individually	
HELP
exit;
}

my $curDir = `pwd`;
chomp $curDir;
my $cd = $curDir;
$cd =~ s/\//_/g;

if($output !~ /^\//){
    $output = "$curDir/$output";
}

my $commandLine = reconstructCL();

my $REF_SEQ = '';
my $FP_INT = '';
my $FP_TG = '';

my $B37_FASTA = '';
my $B37_MM10_HYBRID_FASTA = '';
my $HG19_FASTA = '';
my $HG19_MM10_HYBRID_FASTA = '';
my $MM9_FASTA = '';
my $MM10_FASTA = '';
my $MM10_CUSTOM_FASTA = '';
my $SPECIES_CUSTOM_FASTA = '';
my $SPECIES_CUSTOM_DB_SNP = '';
my $DM3_FASTA = '';

my $targets_ilist = "$Bin/targets/$targets/$targets\_targets.ilist";
my $baits_ilist = "$Bin/targets/$targets/$targets\_baits.ilist";
my $targets_bed_padded = "$Bin/targets/$targets/$targets\_targets_plus5bp.bed";
my $targets_facet = "$Bin/targets/$targets/$targets\_targets_FACETS.ilist";
my $assay = $targets;
if(-d $targets){
    my @path = split(/\//, $targets);
    $assay = pop @path;

    $targets_ilist = "$targets/$assay\_targets.ilist";
    $baits_ilist = "$targets/$assay\_baits.ilist";
    $targets_bed_padded = "$targets/$assay\_targets_plus5bp.bed";
    $targets_facet = "$targets/$assay\_targets_FACETS.ilist";
}

my %grouping = ();
my %grouping_samples = ();
my @run_merge_jids = ();
my $PICARD = '';
my $JAVA = '';
my $PERL = '';
my $PYTHON = '';
my $GATK = '';

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

if($wes && ($impact || $chip)){
    die "Cannot run as wes and impact/chip. Please select -chip OR -wes OR -impact";
} elsif($impact && $chip){
    die "Cannot run as impact and chip. Please select -chip OR -wes OR -impact";
}elsif(!$wes && !$impact && !$chip){
    die "Must specify to run the project as either -chip, -wes or -impact";
}

if($abra && $indelrealigner){
    die "Cannot run both abra and gatk indelrealigner $!";
}
elsif(!$abra && !$indelrealigner){
    $abra = 1;
}

if(!$targets){
    if(!$chip){
	die "NON-CHIPSEQ runs require a target file $!"
    }
}

if(!-d $output){
    mkdir("$output", 0775) or die "Can't make $output";
}

if(!-d $tempdir){
    mkdir("$tempdir", 0775) or die "Can't make $tempdir";
}

&verifyRequest($request);
&verifyConfig($config);
processInputs();

my $svnRev = `svn info $Bin | grep Revision | cut -d " " -f 2`;
chomp $svnRev;

if(!$svnRev){
    my @currentTime = &getTime();
    $svnRev = "$currentTime[5]$currentTime[4]$currentTime[3]";
}

my %samp_libs_run = ();
my $slr_count = 0;
my %ran_solexa = ();
my $projectHasParedEndSamp = 0;

open(LOG, ">$cd\_variants_pipeline.log") or die "can't write to output log";
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSTARTING VARIANTS PIPELINE FOR $pre\n";
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tCOMMAND LINE: $commandLine\n";

if(!-d "$output/intFiles"){
    mkdir("$output/intFiles", 0775) or die "Can't make $output/intFiles";
}
if(!-d "$output/progress"){
    mkdir("$output/progress", 0775) or die "Can't make $output/progress";
}
if(!-d "$output/metrics"){
    mkdir("$output/metrics", 0775) or die "Can't make $output/metrics";
}
if(!-d "$output/metrics/fingerprint"){
    mkdir("$output/metrics/fingerprint", 0775) or die "Can't make $output/metrics/fingerprint";
}
if(!-d "$output/alignments"){
    mkdir("$output/alignments", 0775) or die "Can't make $output/alignments";
}

my $ran_sol = 0;
### NOTE: if a sample consists of both paired and single end runs, it will take the last one
###       this only affects InsertSizeMetrics
###       InsertSizeMetrics doesn't output a file for single end runs; not even a zero-byte file
my %seq_type = ();
alignReads();


### write sample list for qc db
open(SLIST, ">$output/metrics/$pre\_sample_list.txt");
foreach my $samp (keys %samp_libs_run){
    print SLIST "$samp\n";
}
close SLIST;

my $rmdups = 'false';
if($removedups){
    $rmdups = 'true';
}

foreach my $run_merge_jid (@run_merge_jids){
    `$Bin/jobSync $scheduler $run_merge_jid`;
}

my $rmj = join(",", @run_merge_jids);

my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\t$pre\tREADS PROCESSING/ALIGNMENTS FOR ALL SAMPLES COMPLETED\n";
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\t$pre\_$uID\_MARKDUPS RUNNING\n";

my @hsm = ();
my @ism = ();
my @asm = ();
my @mdm = ();
my @cogm = ();
my @docm = ();
my @hfm = ();
my $ran_hs = 0;
my $ran_is = 0;
my $ran_as = 0;
my $ran_md_glob = 0;
my $ran_cog = 0;
my $ran_doc = 0;
my $ran_gcm = 0;
my $ran_hf = 0;
my @md_jids = ();
my @hs_jids = ();
my @is_jids = ();
my @as_jids = ();
my @gc_jids = ();
my @cog_jids = ();
my @doc_jids = ();
my @hf_jids = ();
my %bamsggf = ();
my @r3 = ();

processBams();
foreach my $md_jid (@md_jids){
    `$Bin/jobSync $scheduler $md_jid`;
}
my $mdj = join(",", @md_jids);

if($mdOnly || $chip){
    mergeStats();
    finalSync();
    exit(0);
}

generateGroupFile();
callSNPS();
if (! $nosnps){
    push @r3, "$pre\_$uID\_RSYNC_2";
}

mergeStats();
finalSync();
close LOG;

sub finalSync {
    my $r3j = join(",", @r3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSYNC_3", job_hold => "$r3j", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_RSYNC_3.log");
    my $standardParams = Schedule::queuing(%stdParams);
    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1", mail => "$email");
    my $additionalParams = Schedule::additionalParams(%addParams);
    $ENV{'LSB_JOB_REPORT_MAIL'} = 'Y';
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams /usr/bin/rsync -azvP --exclude 'intFiles' --exclude 'progress' $curDir $rsync`;
    `/bin/touch $output/progress/$pre\_$uID\_RSYNC_3.done`;
}

sub reconstructCL {
### dumb reconstruction of command line
    my $rCL = "$Bin/variants_pipeline.pl";

    $rCL .= " -pre $pre";
    
    if($map){
	$rCL .= " -map $map";
    }
    if($group){
	$rCL .= " -group $group";
    }
    if($pair){
	$rCL .= " -pair $pair";
    }
    if($patient){
        $rCL .= " -patient $patient";
    }
    if($config){
	$rCL .= " -config $config";
    }
    if($request){
        $rCL .= " -request $request";
    }
    if($species){
	$rCL .= " -species $species";
    }
    if($scheduler){
	$rCL .= " -scheduler $scheduler";
    }
    if($targets){
	$rCL .= " -targets $targets";
    }
    if($output){
	$rCL .= " -output $output";
    }
    if($rsync){
	$rCL .= " -rsync $rsync";
    }
    if($email){
	$rCL .= " -email $email";
    }
    if($nosnps){
	$rCL .= " -nosnps";
    }
    if($removedups){
	$rCL .= " -removedups";
    }
    if($wes){
        $rCL .= " -wes";
    }
    if($impact){
        $rCL .= " -impact";
    }
    if($abra){
	$rCL .= " -abra";
    }
    if($mdOnly){
	$rCL .= " -mdOnly";
    }
    if($chip){
	$rCL .= " -chip";
    }
    if($noMD){
	$rCL .= " -noMD";
    }
    if($ug){
	$rCL .= " -unifiedgenotyper";
    }
    
    $rCL .= " -priority_project $priority_project";
    
    $rCL .= " -priority_group $priority_group";
    
    if($species){
	$rCL .= " -species $species";
    }
    
    if($r1adaptor){
	$rCL .= " -r1adaptor $r1adaptor";
    }

    if($r2adaptor){
	$rCL .= " -r2adaptor $r2adaptor";
    }

    if($noClip){
	$rCL .= " -noClip";
    }

    if($bqTrim){
	$rCL .= " -bqTrim $bqTrim";
    }

    if($DB_SNP){
	$rCL .= " -db_snp $DB_SNP";
    }

    if($allSomatic){
	$rCL .= " -allSomatic";
    }

    if($scalpel){
	$rCL .= " -scalpel";
    }

    if($somaticsniper){
	$rCL .= " -somaticsniper";
    }

    if($strelka){
	$rCL .= " -strelka";
    }

    if($varscan){
	$rCL .= " -varscan";
    }

    if($virmid){
	$rCL .= " -virmid";
    }

    if($lancet){
	$rCL .= " -lancet";
    }

    if($vardict){
	$rCL .= " -vardict";
    }

    if($pindel){
        $rCL .= " -pindel";
    }

    my $numArgs = $#ARGV + 1;
    foreach my $argnum (0 .. $#ARGV) {
	$rCL .= " $ARGV[$argnum]";
    }
    
    return $rCL;
}

sub verifyConfig{
    my $paths = shift;

    open(CONFIG, "$paths") || die "Can't open config file $paths $!";
    while(<CONFIG>){
	chomp;
	
	my @conf = split(/\s+/, $_);	
	if($conf[0] =~ /cutadapt/i){
	    if(!-e "$conf[1]/cutadapt"){
		die "CAN'T FIND cutadapt IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /^bwa/i){
	    if(!-e "$conf[1]/bwa"){
		die "CAN'T FIND bwa IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /^abra/i){
	    if(!-e "$conf[1]/abra.jar"){
		die "CAN'T FIND abra.jar IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /picard/i){
	    if(!-e "$conf[1]/picard.jar"){
		die "CAN'T FIND picard.jar IN $conf[1] $!";
	    }
	    $PICARD = $conf[1];
	}
	elsif($conf[0] =~ /gatk/i){
	    if(!-e "$conf[1]/GenomeAnalysisTK.jar"){
		die "CAN'T FIND GenomeAnalysisTK.jar IN $conf[1] $!";
	    }
            $GATK = $conf[1];
	}
	elsif($conf[0] =~ /^mutect$/i){
	    if(!-e "$conf[1]/muTect.jar"){
		die "CAN'T FIND muTect.jar IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /bedtools/i){
	    if(!-e "$conf[1]/bedtools"){
		die "CAN'T FIND bedtools IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /samtools/i){
	    if(!-e "$conf[1]/samtools"){
		die "CAN'T FIND samtools IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /somaticsniper/i){
	    if(!-e "$conf[1]/bam-somaticsniper"){
		if($somaticsniper){
		    die "CAN'T FIND bam-somaticsniper IN $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /varscan/i){
	    if(!-e "$conf[1]/VarScan.jar"){
		if($varscan){
		    die "CAN'T FIND VarScan.jar IN $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /strelka/i){
	    if(!-e "$conf[1]/bin/configureStrelkaWorkflow.pl"){
		if($strelka){
		    die "CAN'T FIND bin/configureStrelkaWorkflow.pl IN $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /scalpel/i){
	    if(!-e "$conf[1]/scalpel"){
		if($scalpel){
		    die "CAN'T FIND scalpel IN $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /virmid/i){
	    if(!-e "$conf[1]/Virmid.jar"){
		if($virmid){
		    die "CAN'T FIND Virmid.jar IN $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /lancet/i){
	    if(!-e "$conf[1]/lancet"){
		die "CAN'T FIND lancet IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /vardict_java/i){
	    if(!-e "$conf[1]/VarDict"){
		die "CAN'T FIND VarDict IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /vardict_perl/i){
	    if(!-e "$conf[1]/testsomatic.R" || !-e "$conf[1]/var2vcf_paired.pl"){
		die "CAN'T FIND testsomatic.R OR var2vcf_paired.pl IN $conf[1] $!";
	    }
	}
        elsif($conf[0] =~ /pindel/i){
            if(!-e "$conf[1]/pindel" or !-e "$conf[1]/pindel2vcf"){
                die "CAN'T FIND pindel or pindel2vcf IN $conf[1] $!";
            }
        }
	elsif($conf[0] =~ /^perl$/i){
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
	elsif($conf[0] =~ /^java$/i){
	    if(!-e "$conf[1]/java"){
		die "CAN'T FIND java IN $conf[1] $!";
	    }
	    $JAVA = $conf[1];
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
	elsif($conf[0] =~ /b37_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /^b37$|human/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR B37 WITH PREFIX $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /b37_mm10_hybrid_fasta/i){
	    if(!-e "$conf[1]"){
		if($species =~ /hybrid|b37_mm10/i){
		    die "CAN'T FIND $conf[1] $!";
		}
	    }
	    $B37_MM10_HYBRID_FASTA = $conf[1];
	}
	elsif($conf[0] =~ /b37_mm10_hybrid_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /hybrid|b37_mm10/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR b37-MM10 HYBRID WITH PREFIX $conf[1] $!";
		}
	    }
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
	}
	elsif($conf[0] =~ /hg19_mm10_hybrid_fasta/i){
	    if(!-e "$conf[1]"){
		if($species =~ /hybrid/i){
		    die "CAN'T FIND $conf[1] $!";
		}
	    }
	    $HG19_MM10_HYBRID_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /hg19_mm10_hybrid_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /hybrid/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19-MM9 HYBRID WITH PREFIX $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /mm9_fasta/i){
	    if(!-e "$conf[1]"){
		if($species =~ /mm9/i){
		    die "CAN'T FIND $conf[1] $!";
		}
	    }
	    $MM9_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /mm9_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /mm9/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM9 WITH PREFIX $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /mm10_fasta/i){
	    if(!-e "$conf[1]"){
		if($species =~ /^mm10$/i){
		    die "CAN'T FIND $conf[1] $!";
		}
	    }
	    $MM10_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /mm10_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /^mm10$/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM10 WITH PREFIX $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /mm10_custom_fasta/i){
	    if(!-e "$conf[1]"){
		if($species =~ /mm10_custom/i){
		    die "CAN'T FIND $conf[1] $!";
		}
	    }
	    $MM10_CUSTOM_FASTA = $conf[1];
	}
	elsif($conf[0] =~ /mm10_custom_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /mm10_custom/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM10_CUSTOM WITH PREFIX $conf[1] $!";
		}
	    }
	}
	elsif($conf[0] =~ /species_custom_fasta/i){
	    if(!-e "$conf[1]"){
		if($species =~ /species_custom/i){
		    die "CAN'T FIND $conf[1] $!";
		}
	    }
	    $SPECIES_CUSTOM_FASTA = $conf[1];
	}
	elsif($conf[0] =~ /species_custom_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /species_custom/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR SPECIES_CUSTOM WITH PREFIX $conf[1] $!";
		}
	    }
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
	    $SPECIES_CUSTOM_DB_SNP = $conf[1];
	}
	elsif($conf[0] =~ /dm3_fasta/i){
	    if(!-e "$conf[1]"){
		if($species =~ /dm3/i){
		    die "CAN'T FIND $conf[1] $!";
		}
	    }
	    $DM3_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /dm3_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		if($species =~ /dm3/i){
		    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR DM3 WITH PREFIX $conf[1] $!";
		}
	    }
	}
   }
    close CONFIG;
}

sub verifyRequest{
    my $reqfile = shift;

    my($pi,$piname,$inv,$invname,$pid,$runnum,$assay,$pipelines,$rerunReason,$origPi,$origInv);  

    open(REQUEST, "$reqfile") || die "Can't open request file $reqfile $!";
    while(<REQUEST>){
        chomp;

        my @req = split(/:\s/, $_);
        if($req[0] =~ /^PI$/i){
            $pi = $req[1];
            $pi =~ s/\@mskcc\.org//;
        }
        elsif($req[0] =~ /^ORIG_PI$/i){
            $origPi = $req[1];
            $origPi =~ s/\@mskcc\.org//;
        }
        elsif($req[0] =~ /^PI_Name$/i){
            $piname = $req[1];
        }
        elsif($req[0] =~ /^Investigator$/i){
            $inv = $req[1];
            $inv =~ s/\@mskcc\.org//;
        } 
        elsif($req[0] =~ /^ORIG_Investigator$/i){
            $origInv = $req[1];
            $origInv =~ s/\@mskcc\.org//;
        }
        elsif($req[0] =~ /^Investigator_Name$/i){
            $invname = $req[1];
        }
        elsif($req[0] =~ /^ProjectID$/i){
            $pid = $req[1];
        }
        elsif($req[0] =~ /^RunNumber$/i){
            $runnum = $req[1];
        }
        elsif($req[0] =~ /^Assay$/i){
            $assay = $req[1];
        }
        elsif($req[0] =~ /^Pipelines$/i){
            $pipelines = $req[1];
        }
        elsif($req[0] =~ /^Reason_for_rerun$/i){
            $rerunReason = $req[1];
        }
    }
    close REQUEST;

    if(!$pi || !$piname || !$inv || !$invname || !$pid || !$runnum || !$assay || !$pipelines) {
        die "\nERROR: Info missing from request file.\n\nREQUIRED fields are:\n  PI\n  PI_Name\n  Investigator\n  Investigator_Name\n  ProjectID\n  RunNumber\n  Assay\n  Pipelines\n\n";
    }

    my $delDir = "/ifs/solres/seq/$pi/$inv/$pid/";
    if($runnum > 1) {
        if(!$rerunReason) {
            die "\nERROR: This is a rerun, but there was no rerun reason given in the request file.\n\n";
        }
        if( ! -d $delDir ){
            if( $origPi && $origInv){
                $delDir = "/ifs/solres/seq/$origPi/$origInv/$pid/";
                if( ! -d $delDir ){
                    die "\nERROR: This looks like a rerun (RunNum: $runnum) but the delivered directory does not exists: $delDir \n\n";
                }
            } else {
                die "\nERROR: This looks like a rerun (RunNum: $runnum) but the delivered directory does not exists: $delDir \n\n";
            }
        }
    } else {
        my $x = sprintf("%03d", $runnum);
        $delDir .= "r_$x";
        if( -d  "$delDir"){ 
            if( $origPi && $origInv){
                my $x = sprintf("%03d", $runnum);
                $delDir = "/ifs/solres/seq/$origPi/$origInv/$pid/r_$x";
                if( -d  "$delDir"){
                #    die "\nERROR: This project has a delivered directory with this rerun folder already present: $delDir \n\n";
                }
            }
            #die "\nERROR: This project has a delivered directory with this rerun folder already present: $delDir \n\n";
        } 

    }
    
    if( $origPi && $origInv){
        $delDir = "/ifs/solres/seq/$origPi/$origInv/$pid/";
	my $x = sprintf("%03d", $runnum);
	$delDir .= "r_$x";
	if( -d  "$delDir"){
	#    die "\nERROR: This project has a delivered directory with this rerun folder already present: $delDir \n\n";
	}
    }
 
    `/bin/cp $request $output/`;
}


sub processInputs {
    if($pre =~ /^\d+/){
	$pre = "s_$pre";
    }
    
    if($species !~ /human|b37|mouse|mm9|mm10|mm10_custom|species_custom|drosophila|dm3|hybrid|xenograft/i){
	die "Species must be human (b37), mouse (mm9|mm10|mm10_custom), species_custom or drosophila (dm3) $!";
    }
    
    if($scheduler =~ /lsf/i){
	$scheduler = 'lsf';
    }
    elsif($scheduler =~ /sge/i){
	$scheduler = 'sge';
    }
    
    if($species =~ /human|^b37$/i){
	$species = 'b37';
	$REF_SEQ = "$B37_FASTA";
	$DB_SNP = "$Bin/data/b37/dbsnp_138.b37.vcf";
	$FP_INT = "$Bin/data/b37/Agilent51MBExome__b37__FP_intervals.list";
	$FP_TG = "$Bin/data/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt";
    }
    elsif($species =~ /hg19/i){

	die "hg19 is no longer supported in the variants pipeline";

	$species = 'hg19';
	$REF_SEQ = "$HG19_FASTA";
	$DB_SNP = "$Bin/data/hg19/dbsnp_138.hg19.vcf";
	$FP_INT = "$Bin/data/hg19/Agilent51MBExome__hg19__FP_intervals.list";
	$FP_TG = "$Bin/data/hg19/Agilent51MBExome__hg19__FP_tiling_genotypes.txt";
    }
    elsif($species =~ /mm9/i){
	$species = 'mm9';
	$REF_SEQ = "$MM9_FASTA";
	$DB_SNP = "";
    }
    elsif($species =~ /mouse|^mm10$/i){
	$species = 'mm10';
	$REF_SEQ = "$MM10_FASTA";
	$DB_SNP = "$Bin/data/mm10/mm10_snp142.vcf";
    }
    elsif($species =~ /mm10_custom/i){
	$species = 'mm10_custom';
	$REF_SEQ = "$MM10_CUSTOM_FASTA";
	$DB_SNP = "$Bin/data/mm10/mm10_snp142.vcf";
    }
    elsif($species =~ /hybrid|xenograft/i){
	$species = 'hybrid';
	$REF_SEQ = "$B37_MM10_HYBRID_FASTA";
	$DB_SNP = "$Bin/data/b37/dbsnp_138.b37.vcf";
	$FP_INT = "$Bin/data/b37/Agilent51MBExome__b37__FP_intervals.list";
	$FP_TG = "$Bin/data/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt";
    }
    elsif($species =~ /species_custom/i){
	$species = 'species_custom';
	$REF_SEQ = "$SPECIES_CUSTOM_FASTA";
	$DB_SNP = "$SPECIES_CUSTOM_DB_SNP";
    }
    elsif($species =~ /drosophila|dm3/i){
	$species = 'dm3';
	$REF_SEQ = "$DM3_FASTA";
	$DB_SNP = "";
    }
      
    if(!-e "$DB_SNP"){
	die "CAN'T FIND dbSNP file $DB_SNP $!";
    }

    if(!-e "$REF_SEQ"){
	die "CAN'T FIND REF SEQ file $REF_SEQ $!";
    }


    if($group){
	open(GROUP, "$group") or die "Can't open grouping file $group $!";
	while(<GROUP>){
	    chomp;
	    
	    my @data = split(/\s+/, $_);
	    $grouping{$data[1]}{$data[0]} = 1;
	    $grouping_samples{$data[0]} = 1;
	}
	close GROUP;
	`/bin/cp $group $output/`;
    }
    else{
	if(!$mdOnly || $chip){
	    die "grouping file required to run pipeline $!";
	}
    }
    
    my %pairing_samples = ();
    if($pair){
	open(PAIR, "$pair") or die "Can't open pairing file $pair $!";
	while(<PAIR>){
	    chomp;
	    
	    my @data = split(/\s+/, $_);
	    foreach my $dat (@data){
		if($dat =~ /^NA$/i){
		    next;
		}
		$pairing_samples{$dat} = 1;
	    }
	}
	close PAIR;
	`/bin/cp $pair $output/`;
    }
    
    my %mapping_samples = ();
    my %mapaths = ();
    open(MA, "$map") or die "Can't open mapping file $map $!";
    while(<MA>){
	chomp;
	
	my @data = split(/\s+/, $_);
	
	$mapping_samples{$data[1]} = 1;
	if(!-d $data[3]){
	    die "$data[3] does not exist\n";
	}

	if(!$mapaths{$data[3]}){
	    $mapaths{$data[3]} = 1;
	}
	else{
	    die "fastq path $data[3] is in the mapping file $map multiple times $!";
	}

	if($group){
	    if(!$grouping_samples{$data[1]}){
		die "grouping file $group missing sample $data[1] found in mapping file $map $!";
	    }
	}
    }
    close MA;
    `/bin/cp $map $output/`;
        
    # Check patient sample has fields needed
    # Patient file needs only "Sample_ID", "Patient_ID", "Class", and "Bait_version"
    my %patient_samples = ();
    if($patient){
        open(PATIENT, "$patient") or die "Can't open patient file $patient $!";
        my $h = <PATIENT>;
        my @h = split(/\s+/, $h);
        my @mandatory_headers = ("Sample_ID", "Patient_ID", "Class", "Bait_version");
        foreach my $loopVal (@mandatory_headers){
            if(! grep(/^$loopVal$/, @h)) {
                die "Missing header value $loopVal in patient file $patient $!";
            }
        }
        my ($sID_index) = grep {$h[$_] =~ /Sample_ID/} 0..$#h;
        my ($class_index) = grep {$h[$_] =~ /Class/} 0..$#h;

        while(<PATIENT>){
            my @line = split(/\s+/,$_);
            my $sid = $line[$sID_index];

            $patient_samples{$sid}=1;

            if(!$mapping_samples{$sid}) {
                 die "mapping file $map missing sample $sid found in patient file $patient $!";
            }

            if(!$grouping_samples{$sid}) {
                 die "grouping file $group missing sample $sid found in patient file $patient $!";
            }

            if(!$pairing_samples{$sid} &&  $line[$class_index] ne "PoolNormal") {
                 die "pairing file $pair missing sample $sid found in patient file $patient $!";
            }
        }
        close PATIENT;
	`/bin/cp $patient $output/`;

    }
 
    if($group){
	foreach my $gro (keys %grouping_samples){
	    if(!$mapping_samples{$gro}){
		die "grouping file $group contains sample $gro that isn't in mapping file $map $!";
	    }
	    
	    if($pair){
                # add this extra check because pooled normals do not have to be paired if all samples have
                # matched normals
		if(!$pairing_samples{$gro}){
                    if($patient){
                        if( !$patient_samples{$gro} ){
                            die "grouping file $group contains sample $gro that isn't in pairing file $pair or patient file $patient $!";
                        }
                    } elsif ($gro !~ /^s_Normal_Pooled/) {
                        # If there is no patient file, just exit because we can't tell if the sample
                        # is a pooled normal or not
		        die "grouping file $group contains sample $gro that isn't in pairing file $pair $!";
                   }  
	       }
           }
	}
    }

    if($pair){
	foreach my $pai (keys %pairing_samples){
	    if(!$mapping_samples{$pai}){
		die "pairing file $pair contains sample $pai that isn't in mapping file $map $!";
	    }
	    
	    if(!$grouping_samples{$pai}){
		die "pairing file $pair contains sample $pai that isn't in grouping file $group $!";
	    }
            if($patient){
                if(!$patient_samples{$pai}){
                    die "pairing file $pair contains sample $pai that is not found in patient file $patient $!";
                }
            }
        }
    }
    
    if(!-e $targets_ilist || !-e $baits_ilist || !-e $targets_bed_padded){
	if(!$chip){
	    die "directory $targets OR $Bin/targets/$targets MUST CONTAIN THE FOLLOWING FILES: $targets\_targets.ilist; $targets\_baits.ilist; $targets\_targets_plus5bp.bed $!";
	}	
    }

    if($species =~ /b37|hg19|hybrid|xenograft/i){
	if(!-e "$targets_facet"){
	    if(!$chip){
		die "directory $targets OR $Bin/targets/$targets MUST CONTAIN THE FOLLOWING FILES FOR FACETS TO RUN: $targets\_targets_FACETS.ilist $!";
	    }	
	}
    }

}


sub getTime(){

    my @timeData = localtime();

    if($timeData[0] < 10){
	$timeData[0] = "0" . $timeData[0];
    }
    if($timeData[1] < 10){
	$timeData[1] = "0" . $timeData[1];
    }
    if($timeData[2] < 10){
	$timeData[2] = "0" . $timeData[2];
    }
    if($timeData[3] < 10){
	$timeData[3] = "0" . $timeData[3];
    }
    
    my $month = $timeData[4] + 1;
    $timeData[4] += 1;
    if($timeData[4] < 10){
	$timeData[4] = "0" . $timeData[4];
    }
    
    $timeData[5] += 1900;

    return(@timeData);
}

sub alignReads {

    my $no_clip = '';
    my $aSeq = '';
    if($noClip){
	$no_clip = '-noClip'
    }
    else{
	$aSeq = "-r1adaptor $r1adaptor -r2adaptor $r2adaptor";
    }

    my $defaultBWA = '';
    if($chip){
	$defaultBWA = '-defaultBWA';
    }

    open(IN, "$map") or die "Can't open $map $!";
    while(<IN>){
	chomp;
	
	my @data = split(/\s+/, $_);
	if($data[1] =~ /^\d+/){
	    $data[1] = "s_" . "$data[1]";
	}
	
	### sometimes in the mapping file, 
	### the lib, sample, and run id isn't a unique identifier
	if($samp_libs_run{$data[1]}{$data[0]}{$data[2]}){
	    $slr_count++;
	    my @currentTime = &getTime();
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tWARNING: $data[1]\t$data[0]\t$data[2] ISN'T UNIQUE; ";
	    $data[2] = "$data[2]\_$slr_count";
	    print LOG "WRITING INSTEAD TO $data[1]\/$data[0]\/$data[2]\n";
	}
	
	if(!-d "$output/intFiles/$data[1]"){
	    mkdir("$output/intFiles/$data[1]", 0775) or die "Can't make $output/intFiles/$data[1]";
	}

	if(!-d "$output/intFiles/$data[1]/$data[0]"){
	    mkdir("$output/intFiles/$data[1]/$data[0]", 0775) or die "Can't make $output/intFiles/$data[1]/$data[0]";
	}

	if(!-d "$output/intFiles/$data[1]/$data[0]/$data[2]"){
	    mkdir("$output/intFiles/$data[1]/$data[0]/$data[2]", 0775) or die "Can't make $output/intFiles/$data[1]/$data[0]/$data[2]";
	}

	$samp_libs_run{$data[1]}{$data[0]}{$data[2]} = 1;	
	
	`ln -s $data[3]/* $output/intFiles/$data[1]/$data[0]/$data[2]/`;
	chdir "$output/intFiles/$data[1]/$data[0]/$data[2]";
	
	opendir(workDir, "./");
	my @unsorted = readdir workDir;
	closedir workDir;
	my @files = sort @unsorted;
	
	open(OUT, ">files_$data[1]\_$data[0]\_$data[2]");
	foreach my $file (@files){
	    if($file =~ /fastq\.gz$/ && $file =~ m/^(.*)(R\d+)(.*)$/){
		if($2 eq 'R1'){
		    my $file_R2 = $file;
		    $file_R2 =~ s/^(.*)R1(.*)$/$1R2$2/;

		    ### if PE or SE designation in mapping file, will assume that's what's wanted
		    ### if it's not there, will determine it based on what's in reads dir
		    if($data[4]){
			if($data[4] =~ /pe/i){
			    $seq_type{$data[1]} = "PE";
			    print OUT "$file\t$file_R2\n";
			}
			elsif($data[4] =~ /se/i){
			    $seq_type{$data[1]} = "SE";
			    print OUT "$file\n";
			}
                        else{ 
                            die "CAN'T DETERMINE WHETHER RUN IS PAIRED OR SINGLE ENDED for sample $data[1]\n"; 
                        }	
		    }
		    else{
			if(-e "$file_R2"){
			    $seq_type{$data[1]} = "PE";
			    print OUT "$file\t$file_R2\n";
			}
			else{
			    $seq_type{$data[1]} = "SE";
			    print OUT "$file\n";
			}
		    }	
		}
	    }
	}
	close OUT;
	
	my @currentTime = &getTime();
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSTARTING READS PROCESSING/ALIGNMENT FOR $data[1]\_$data[0]\_$data[2]\n";

	if(!-e "$output/progress/reads_files_$data[1]\_$data[0]\_$data[2]\.done"){	    
	    `$Bin/process_reads.pl -file files_$data[1]\_$data[0]\_$data[2] -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\tID:$data[1]\_$data[0]\_$data[2]\_$seq_type{$data[1]}\\tPL:Illumina\\tPU:$data[1]\_$data[0]\_$data[2]\\tLB:$data[1]\_$data[0]\\tSM:$data[1]" -species $species -config $config -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $no_clip $aSeq -bqTrim $bqTrim $defaultBWA > files_$data[1]\_$data[0]\_$data[2]\_process_reads_$seq_type{$data[1]}\.log 2>&1`;
	    ###`/common/sge/bin/lx24-amd64/qsub /home/mpirun/tools/qCMD $Bin/solexa_PE.pl -file files -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\\tID:$data[1]\_$data[0]\_$data[2]\_PE\\\tPL:Illumina\\\tPU:$data[1]\_$data[0]\_$data[2]\\\tLB:$data[1]\_$data[0]\\\tSM:$data[1]" -species $species -config $config $targeted -scheduler $scheduler`;
	    $ran_sol = 1;
	    $ran_solexa{$data[1]} = 1;
	    push @run_merge_jids, "$pre\_$uID\_$data[1]\_$data[0]\_$data[2]\_MERGE";
	    `/bin/touch $output/progress/reads_files_$data[1]\_$data[0]\_$data[2]\.done`;
	}
	else{
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSKIPPING READS PROCESSING/ALIGNMENT FOR $data[1]\_$data[0]\_$data[2]; PREVIOUSLY RAN TO COMPLETION\n";
	}
	chdir $curDir;
    }
    close IN;
}

sub processBams {
    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
    my $additionalParams = Schedule::additionalParams(%addParams);

    foreach my $samp (keys %samp_libs_run){
	my @sBams = ();
	my @lib_merge_samp_jids = ();
	my @mdBams = ();
	foreach my $lib (keys %{$samp_libs_run{$samp}}){
	    my @lBams = ();
	    my $ran_md = 0;
	    foreach my $run (keys %{$samp_libs_run{$samp}{$lib}}){
		if(!-e "$output/intFiles/$samp/$lib/$run/$samp\_$lib\_$run\.bam"){
		    my @currentTime = &getTime();
		    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tCan't locate $output/intFiles/$samp/$lib/$run/$samp\_$lib\_$run\.bam ...EXITING VARIANTS PIPELINE FOR PROJECT $pre";
		    die;
		}
		push @lBams, "I=$output/intFiles/$samp/$lib/$run/$samp\_$lib\_$run\.bam";
	    }
	    
	    my $fin = join(" ", @lBams);
	    my $ran_lb_merge = 0;
	    ### NOTE: DOESN'T CHECK TO SEE IF SOLEXA_PE OR SOLEXA_SE WAS RUN
	    ###       WILL AUTOMATICALLY RUN IF IT DOESN'T FIND THE .done FILE FOR THE SAMPLE LIB
	    ###       THIS IS TO AVOID RUNNING THIS CPU COSTLY STEP IF JUST 1 OR SO SAMPLE
	    ###       IN A COHORT NEEDED TO HAVE THEIR READS REPROCESSED
	    if(scalar(@lBams) == 1){
		push @sBams, "$fin";
		if($noMD){
		    my @lbt = split(/=/, $lBams[0]);
		    $bamsggf{$samp}{$lib} = "$lbt[1]";
		}
		else{
		    if(!-e "$output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done" || $ran_solexa{$samp}){
			my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$samp\_$lib\_MARKDUPS", job_hold => "$rmj", cpu => "3", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.log");
			my $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MarkDuplicates $fin OUTPUT=$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam METRICS_FILE=$output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt TMP_DIR=$tempdir VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$rmdups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;
			
			`/bin/touch $output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done`;
			push @md_jids, "$pre\_$uID\_$samp\_$lib\_MARKDUPS";
			$ran_md = 1;
			$ran_md_glob = 1;
		    }
		    $bamsggf{$samp}{$lib} = "$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam"; 
		}
	    }
	    else{
		push @sBams, "I=$output/intFiles/$samp/$lib/$samp\_$lib\.bam";
		if(!-e "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done" || $ran_solexa{$samp}){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_LIB_MERGE_$samp\_$lib", job_hold => "$rmj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $fin O=$output/intFiles/$samp/$lib/$samp\_$lib\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;

		    `/bin/touch $output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done`;
		    $ran_lb_merge = 1;
		    push @lib_merge_samp_jids, "$pre\_$uID\_LIB_MERGE_$samp\_$lib";
		    push @md_jids, "$pre\_$uID\_LIB_MERGE_$samp\_$lib";		    
		    sleep(2);
		}
		
		if($noMD){
		    $bamsggf{$samp}{$lib} = "$output/intFiles/$samp/$lib/$samp\_$lib\.bam";
		}
		else{
		    if(!-e "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done" || $ran_lb_merge){
			my $lb_merge_hold = "";
			if($ran_lb_merge){
			    $lb_merge_hold = "$pre\_$uID\_LIB_MERGE_$samp\_$lib";
			}		
			
			my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$samp\_$lib\_MARKDUPS", job_hold => "$lb_merge_hold", cpu => "3", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.log");
			my $standardParams = Schedule::queuing(%stdParams);
			`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MarkDuplicates INPUT=$output/intFiles/$samp/$lib/$samp\_$lib\.bam OUTPUT=$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam METRICS_FILE=$output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt TMP_DIR=$tempdir VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$rmdups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;
			
			`/bin/touch $output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done`;
			push @md_jids, "$pre\_$uID\_$samp\_$lib\_MARKDUPS";
			$ran_md = 1;
			$ran_md_glob = 1;
		    }
		    $bamsggf{$samp}{$lib} = "$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam";  
		}
	    }
	    
	    push @mdm, "-metrics $output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt";
	    push @mdBams, "I=$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam";
	}
	
	if($mdOnly || $chip){
	    my $mdb = join(" ", @mdBams);
	    if(!-e "$output/progress/$pre\_$uID\_SAMP_MD_MERGE_$samp\.done" || $ran_md_glob){
		my $mdj = join(",", @md_jids);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SAMP_MD_MERGE_$samp", job_hold => "$mdj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_SAMP_MD_MERGE_$samp\.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $mdb O=$output/alignments/$pre\_$samp\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_SAMP_MD_MERGE_$samp\.done`;
                push @r3, "$pre\_$uID\_SAMP_MD_MERGE_$samp";
	    }
	}
	
	my $rin = join(" ", @sBams);
	my $bamForStats = '';
	my $bamForHybridFilter = '';
	my $ran_samp_merge = 0;
	my $ran_samp_query_sort = 0;	
	my @samp_merge_samp_jids = ();
	my @samp_query_sort_jids = ();
	my $lmsj = join(",", @lib_merge_samp_jids);
	if(scalar(@sBams) == 1){
	    my @bname = split(/=/, $sBams[0]);
	    $bamForStats = "$bname[1]";
	}
	else{
	    if(!-e "$output/progress/$pre\_$uID\_SAMP_MERGE_$samp\.done" || $ran_solexa{$samp}){
		my $lmsj = join(",", @lib_merge_samp_jids);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SAMP_MERGE_$samp", job_hold => "$rmj,$lmsj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_SAMP_MERGE_$samp\.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $rin O=$output/intFiles/$samp/$samp\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_SAMP_MERGE_$samp\.done`;
		push @samp_merge_samp_jids, "$pre\_$uID\_SAMP_MERGE_$samp";
		$ran_samp_merge = 1;
	    }
	    $bamForStats = "$output/intFiles/$samp/$samp\.bam";
	}
	
	if($species =~ /hybrid/)
	{
		if(!-e "$output/progress/$pre\_$uID\_SAMP_QUERY_SORT_$samp\.done" || $ran_solexa{$samp}){
                	my $lmsj = join(",", @lib_merge_samp_jids);
                	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SAMP_QUERY_SORT_$samp", job_hold => "$rmj,$lmsj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_SAMP_QUERY_SORT_$samp\.log");
                	my $standardParams = Schedule::queuing(%stdParams);
                	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $rin O=$output/intFiles/$samp/$samp\_query_sort.bam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
                	`/bin/touch $output/progress/$pre\_$uID\_SAMP_QUERY_SORT_$samp\.done`;
                	push @samp_query_sort_jids, "$pre\_$uID\_SAMP_QUERY_SORT_$samp";
                	$ran_samp_query_sort = 1;
		}
		$bamForHybridFilter = "$output/intFiles/$samp/$samp\_query_sort.bam";	


	        my $sqsj = join(",", @samp_query_sort_jids);
	        if((!-e "$output/progress/$pre\_$uID\_HYBRID_FILTER_$samp\.done" || $ran_solexa{$samp} || $ran_samp_query_sort) && !$chip){
	            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HYBRID_FILTER\_$samp", job_hold => "$rmj,$sqsj,$lmsj", cpu => "1", mem => "5", cluster_out => "$output/progress/$pre\_$uID\_HYBRID_FILTER_$samp\.log");
	            my $standardParams = Schedule::queuing(%stdParams);
	            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/Bam2FastqFilterContamination --input $bamForHybridFilter ">$output/intFiles/$pre\_HybridFilter_$samp\.txt"`;
	            `/bin/touch $output/progress/$pre\_$uID\_HYBRID_FILTER_$samp\.done`;
	            push @hf_jids, "$pre\_$uID\_HYBRID_FILTER\_$samp";
	            $ran_hf = 1;
	        }
	        push @hfm, "-metrics $output/intFiles/$pre\_HybridFilter_$samp\.txt";
	}

	
	my $smsj = join(",", @samp_merge_samp_jids);
	### HS will fail if seq dict of ref seq and bait/target intervals don't match
        if((!-e "$output/progress/$pre\_$uID\_HS_METRICS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge) && !$chip){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HS_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_HS_METRICS_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar CalculateHsMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_HsMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE BAIT_INTERVALS=$baits_ilist BAIT_SET_NAME=$assay TARGET_INTERVALS=$targets_ilist VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_HS_METRICS_$samp\.done`;
	    push @hs_jids, "$pre\_$uID\_HS_METRICS\_$samp";
	    $ran_hs = 1;
	}
	push @hsm, "-metrics $output/intFiles/$pre\_HsMetrics_$samp\.txt";
	
	### NOTE: will only run if the sample is paired end
	if($seq_type{$samp} eq "PE"){
            $projectHasParedEndSamp=1;
	    if(!-e "$output/progress/$pre\_$uID\_IS_METRICS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_IS_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_IS_METRICS_$samp\.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar CollectInsertSizeMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_InsertSizeMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE HISTOGRAM_FILE=$output/intFiles/$pre\_InsertSizeMetrics_Histogram_$samp\.txt VALIDATION_STRINGENCY=LENIENT`;
		`/bin/touch $output/progress/$pre\_$uID\_IS_METRICS_$samp\.done`;
		push @is_jids, "$pre\_$uID\_IS_METRICS\_$samp";
		$ran_is = 1;
	    }	
	    push @ism, "-metrics $output/intFiles/$pre\_InsertSizeMetrics_$samp\.txt";
	}
	
	if(!-e "$output/progress/$pre\_$uID\_AS_METRICS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AS_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_AS_METRICS_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_AlignmentSummaryMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_AS_METRICS_$samp\.done`;
	    push @as_jids, "$pre\_$uID\_AS_METRICS\_$samp";
	    $ran_as = 1;
	}
	push @asm, "-metrics $output/intFiles/$pre\_AlignmentSummaryMetrics_$samp\.txt";

	if(!-e "$output/progress/$pre\_$uID\_COG_METRICS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_COG_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_COG_METRICS_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar CollectOxoGMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_CollectOxoGMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ DB_SNP=$DB_SNP VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_COG_METRICS_$samp\.done`;
	    push @cog_jids, "$pre\_$uID\_COG_METRICS\_$samp";
	    $ran_cog = 1;
	}
	push @cogm, "-metrics $output/intFiles/$pre\_CollectOxoGMetrics_$samp\.txt";

	if($species =~ /b37|hg19|hybrid/){
            if(!-e "$output/progress/$pre\_$uID\_DOC_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
                my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_DOC\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "4", cluster_out => "$output/progress/$pre\_$uID\_DOC_$samp\.log");
                my $standardParams = Schedule::queuing(%stdParams);
                `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -I $bamForStats -R $REF_SEQ -o $bamForStats\_FP_base_counts\.txt -L $FP_INT -rf BadCigar -mmq 20 -mbq 0 -omitLocusTable -omitSampleSummary -baseCounts --includeRefNSites -omitIntervals`;
                `/bin/touch $output/progress/$pre\_$uID\_DOC_$samp\.done`;
                push @doc_jids, "$pre\_$uID\_DOC\_$samp";
                $ran_doc = 1;
            }
            push @docm, "$bamForStats\_FP_base_counts\.txt";
        }

        if(!-e "$output/progress/$pre\_$uID\_GC_METRICS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_GC_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_GC_METRICS_$samp\.log");
            my $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar CollectGcBiasMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_GcBiasMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ SUMMARY_OUTPUT=$output/intFiles/$pre\_gcbias_$samp\_summary_tmp.txt CHART_OUTPUT=$output/intFiles/$pre\_gcbias_$samp\_tmp.pdf`;
            `/bin/touch $output/progress/$pre\_$uID\_GC_METRICS_$samp\.done`;
            push @gc_jids, "$pre\_$uID\_GC_METRICS\_$samp";
            $ran_gcm = 1;
        }

    }
}

sub mergeStats {
    my @qcpdf_jids = ();

    if($nosnps){
        push @qcpdf_jids, "$pre\_$uID\_RSYNC_1"; 
    }
    elsif(!$mdOnly && !$chip){
	push @qcpdf_jids, "$pre\_$uID\_RSYNC_2";
    }
    

    my $ran_merge = 0;
    my $ran_merge_ism = 0;
    my %addParams = (scheduler => "$scheduler", runtime => "50", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "1");
    my $additionalParams = Schedule::additionalParams(%addParams);

    if(!$noMD || $chip){
	my $mdfiles = join(" ", @mdm);
	if(!-e "$output/progress/$pre\_$uID\_MERGE_MDM.done" || $ran_md_glob){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_MDM", job_hold => "$mdj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_MDM.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $mdfiles ">$output/metrics/$pre\_markDuplicatesMetrics.txt"`;
	    `/bin/touch $output/progress/$pre\_$uID\_MERGE_MDM.done`;
	    push @qcpdf_jids, "$pre\_$uID\_MERGE_MDM";
	    $ran_merge = 1;
	}
    }
   
    my $hffiles = join(" ", @hfm);
    my $hfj = join(",", @hf_jids);
    if((!-e "$output/progress/$pre\_$uID\_MERGE_HFM.done" || $ran_hf) && !$chip && $species =~ /hybrid/)
    {
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_HFM", job_hold => "$hfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_HFM.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergeHybridMetrics.pl $hffiles ">$output/metrics/$pre\_HybridMetrics.txt"`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_HFM.done`;
        push @qcpdf_jids, "$pre\_$uID\_MERGE_HFM";
        $ran_merge = 1;
    }

    my $hsfiles = join(" ", @hsm);
    my $hsj = join(",", @hs_jids);
    if((!-e "$output/progress/$pre\_$uID\_MERGE_HSM.done" || $ran_hs) && !$chip){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_HSM", job_hold => "$hsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_HSM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $hsfiles ">$output/metrics/$pre\_HsMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_HSM.done`;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_HSM";
	$ran_merge = 1;
    }
    
    my $isfiles = join(" ", @ism);
    my $isj = join(",", @is_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_ISM.done" && $projectHasParedEndSamp || $ran_is) {
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ISM", job_hold => "$isj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ISM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $isfiles ">$output/metrics/$pre\_InsertSizeMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_ISM.done`;
	$ran_merge_ism = 1;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_ISM";
	$ran_merge = 1;
    }

    if((!-e "$output/progress/$pre\_$uID\_ISM_MATRIX.done" || $ran_merge_ism) && $projectHasParedEndSamp){    
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ISM_MATRIX", job_hold => "$pre\_$uID\_MERGE_ISM", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ISM_MATRIX.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeInsertSizeHistograms.py $output/intFiles InsertSizeMetrics_*.txt $output/metrics/$pre\_InsertSizeMetrics_Histograms.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_ISM_MATRIX.done`;
	push @qcpdf_jids, "$pre\_$uID\_ISM_MATRIX";
	$ran_merge = 1;
    }
    
    my $asfiles = join(" ", @asm);
    my $asj = join(",", @as_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_ASM.done" || $ran_as){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ASM", job_hold => "$asj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ASM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $asfiles ">$output/metrics/$pre\_AlignmentSummaryMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_ASM.done`;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_ASM";
	$ran_merge = 1;
    }
    
    my $cogfiles = join(" ", @cogm);
    my $cogj = join(",", @cog_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_COGM.done" || $ran_cog){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_COGM", job_hold => "$cogj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_COG.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $cogfiles ">$output/metrics/$pre\_CollectOxoGMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_COGM.done`;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_COGM";
	$ran_merge = 1;
    }

    my $gcj = join(",", @gc_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_GCM.done" || $ran_gcm){
       ### NOTE: all of the gcbiasmetrics jobs should be done because of the job sync on merge jobs
       ###       so don't need this to hold on gcbiasmetrics
       ###       this hold causes issues on lsf because of the *
       my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_GCM", job_hold => "$gcj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_GCM.log");
       my $standardParams = Schedule::queuing(%stdParams);
       `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeGcBiasMetrics.py $output/intFiles GcBiasMetrics* $output/metrics/$pre\_GcBiasMetrics.txt`;
       `/bin/touch $output/progress/$pre\_$uID\_MERGE_GCM.done`;
       push @qcpdf_jids, "$pre\_$uID\_MERGE_GCM";
       $ran_merge = 1;
    }

    if($species =~ /b37|hg19|hybrid/){
        my $docFiles = join(" ", @docm);
        my $docj = join(",", @doc_jids);
        if(!-e "$output/progress/$pre\_$uID\_FINGERPRINTING.done" || $ran_doc){
            my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_FINGERPRINTING", job_hold => "$docj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_FINGERPRINTING.log");
            my $standardParams = Schedule::queuing(%stdParams);
            `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/analyzeFingerprint.py -pattern '*_FP_base_counts.txt' -pre $pre -fp $FP_TG -group $group -pair $pair -outdir $output/metrics/fingerprint -dir $output/intFiles`; 
            `/bin/touch $output/progress/$pre\_$uID\_FINGERPRINTING.done`;
            push @qcpdf_jids, "$pre\_$uID\_FINGERPRINTING";
        }
    }

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CAS.done" || $ran_sol){
	###my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CAS", job_hold => "$pre\_$uID\_CUTADAPT*", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CAS.log");
	### NOTE: all of the cutadapt jobs should be done because of the job sync on merge jobs
	###       so don't need this to hold on cutadapt
	###       this hold causes issues on lsf because of the *
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CAS", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CAS.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeCutAdaptStats.py $output CUTADAPT_STATS.txt $output/metrics/$pre\_CutAdaptStats.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CAS.done`;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_CAS";
    }

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CQS.done" || $ran_sol){
       ### NOTE: all of the cqs jobs should be done because of the job sync on merge jobs
       ###       so don't need this to hold on cutadapt
       ###       this hold causes issues on lsf because of the *
       my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CQS", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CQS.log");
       my $standardParams = Schedule::queuing(%stdParams);
       `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeConvertQualityScoreMetrics.py $output _cqs_metrics $output/metrics/$pre\_ConvertQualityScoreMetrics.txt`;
       `/bin/touch $output/progress/$pre\_$uID\_MERGE_CQS.done`;
        #push @qcpdf_jids, "$pre\_$uID\_MERGE_CQS";
    }

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CQS_LOGS.done" || $ran_sol){
       ### NOTE: all of the cqs jobs should be done because of the job sync on merge jobs
       ###       so don't need this to hold on cutadapt
       ###       this hold causes issues on lsf because of the *
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CQS_LOGS", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CQS_LOGS.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeConvertQualityScoreLogs.py $output _cqs_log $output/metrics/$pre\_ConvertQualityScoreLog.txt`;
        `/bin/touch $output/progress/$pre\_$uID\_MERGE_CQS_LOGS.done`;
     }

    my $ran_qcpdf = 0; 
    my $qcpdfj = join(",", @qcpdf_jids);
    #if(!-e "$output/progress/$pre\_$uID\_QCPDF.done" || $ran_merge){
    if(!$noMD){
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QCPDF", job_hold => "$qcpdfj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_QCPDF.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/qc/qcPDF.pl -path $output/metrics -pre $pre -config $config -request $request -log $output/progress/$pre\_$uID\_QCPDF_ERRORS.log -v $svnRev`;
        `/bin/touch $output/progress/$pre\_$uID\_QCPDF.done`;
        push @r3, "$pre\_$uID\_QCPDF";
        $ran_qcpdf = 1;
    }
    #}

    my $qcdbj = join(",", @qcpdf_jids);
    my $chp = '';
    my $se = '';
    if($chip){
        $chp = '--chip';
    }
    if( grep( /^SE$/, (values %seq_type))){
        $se = "--se";
    }

    if(!-e "$output/progress/$pre\_$uID\_QCDB.done" || $ran_merge || $ran_qcpdf){
        my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QCDB", job_hold => "$qcdbj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_QCDB.log");
        my $standardParams = Schedule::queuing(%stdParams);
        `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/db/load_exome_project.py -o $output/progress/$pre\_$uID\_QCDB_LOG.log $chp $se $request $svnRev $output\n`;
        `/bin/touch $output/progress/$pre\_$uID\_QCDB.done`;
        push @r3, "$pre\_$uID\_QCDB";
    }
            
}

sub generateGroupFile {
    my @currentTime = &getTime();
    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\t$pre\_$uID\_MARKDUPS COMPLETED\n";
    
    open(GR, ">$output/intFiles/$pre\_MDbams_groupings.txt") or die "Can't write $output/intFiles/$pre\_MDbams_groupings.txt $!";
    foreach my $grou (keys %grouping){
	my @groupings = ();
	foreach my $sample (keys %{$grouping{$grou}}){
	    foreach my $lib (keys %{$samp_libs_run{$sample}}){
		if(!-e "$bamsggf{$sample}{$lib}"){
		    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tCAN'T LOCATE $bamsggf{$sample}{$lib}...EXITING VARIANTS PIPELINE FOR PROJECT $pre";
		    die;
		}
		###push @groupings, "$output/intFiles/$sample/$lib/$sample\_$lib\_MD.bam";
		push @groupings, "$bamsggf{$sample}{$lib}";
	    }
	}
	
	if(scalar(@groupings) < 1){
	    my @currentTime = &getTime();
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tCAN'T LOCATE ANY MARKDUP BAMS FOR GROUP $grou ...EXITING VARIANTS PIPELINE FOR PROJECT $pre";
	    die;
	}
	
	my $gro = join(",", @groupings);
	print GR "$grou\t$gro\n";
    }
    close GR;
}

sub callSNPS {
    my $callSnps = '';
    if($nosnps){
	$callSnps = '-nosnps';
    }
    
    
    my $paired = '';
    if($pair){
	$paired = "-pair $pair";
    }

    my $patientFile = '';
    if($patient){
        $patientFile = "-patient $patient";
    }
    
    my $run_ug = '';
    if($ug){
	$run_ug = "-ug";
    }

    ###my $run_abra = '';
    ###if($abra){
###	$run_abra = "-abra";
   ### }

    my $wesImpact = '';
    if($wes){
        $wesImpact = '-wes';
    } else {
        $wesImpact = '-impact';
    }

    my $run_ir = '';
    if($indelrealigner){
	$run_ir = '-indelrealigner';
    }
    
    my $run_step1 = '';
    if($ran_md_glob){
	$run_step1 = '-step1';
    }

    my $somaticCallers = '';
    if($allSomatic){
	$somaticCallers .= " -allSomatic";
    }
    if($scalpel){
	$somaticCallers .= " -scalpel";
    }
    if($somaticsniper){
	$somaticCallers .= " -somaticsniper";
    }
    if( $strelka){
	$somaticCallers .= " -strelka";
    }
    if( $varscan){
	$somaticCallers .= " -varscan";
    }
    if($virmid){
	$somaticCallers .= " -virmid";
    }
    if($lancet){
	$somaticCallers .= " -lancet";
    }
    if($pindel){
        $somaticCallers .= " -pindel";
    }

    my @currentTime = &getTime();
    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSNP CALLING RUNNING\n";
    ###my %addParams = (scheduler => "$scheduler", runtime => "50", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "1");
    ###my $additionalParams = Schedule::additionalParams(%addParams);
    ###my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SNP_PIPE", job_hold => "$mdj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_SNP_PIPE.log");
    ###my $standardParams = Schedule::queuing(%stdParams);
    
    if($species =~ /b37|hg19|hybrid|xenograft/i){
        ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_human.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_abra $run_step1 $patientFile -species $species`;

	print LOG "$Bin/process_alignments_human.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps $wesImpact -targets $targets $run_ug -email $email -output $output -svnRev $svnRev -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_ir $run_step1 $patientFile -species $species -rsync $rsync $somaticCallers\n";
        `$Bin/process_alignments_human.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps $wesImpact -targets $targets $run_ug -email $email -output $output -svnRev $svnRev -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_ir $run_step1 $patientFile -species $species -rsync $rsync $somaticCallers > $output/progress/$pre\_process_alignments_human.log 2>&1`;
    }
    elsif($species =~ /mm9|^mm10$|mm10_custom/i){
        ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_mouse.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group -species $species $run_abra $run_step1`;

	print LOG "$Bin/process_alignments_mouse.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -email $email -output $output -svnRev $svnRev -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_ir $run_step1 -species $species -rsync $rsync $somaticCallers\n";
        `$Bin/process_alignments_mouse.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -email $email -output $output -svnRev $svnRev -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_ir $run_step1 -species $species -rsync $rsync $somaticCallers > $output/progress/$pre\_process_alignments_mouse.log 2>&1`;
    }
    elsif($species =~ /species_custom/i){
        ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_species_custom.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_abra $run_step1`;

	print LOG "$Bin/process_alignments_species_custom.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -email $email -output $output -svnRef $svnRev -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_step1 -rsync $rsync $somaticCallers\n";
        `$Bin/process_alignments_species_custom.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -email $email -output $output -svnRef $svnRev -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_step1 -rsync $rsync $somaticCallers > $output/progress/$pre\_process_alignments_species_custom.log 2>&1`;
    }
    elsif($species =~ /dm3/i){
        ###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_dm3.pl -pre $pre -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps $run_ug -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_abra $run_step1 -species $species -rsync $rsync $somaticCallers > $output/progress/process_alignments_drosophila.log 2>&1`;
    }
} 
