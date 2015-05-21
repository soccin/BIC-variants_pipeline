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

my ($map, $group, $pair, $config, $help, $nosnps, $removedups, $species, $ug, $scheduler, $abra, $targets);

my $pre = 'TEMP';
my $output = "results";
my $priority_project = "ngs";
my $priority_group = "Pipeline";

GetOptions ('map=s' => \$map,
	    'group=s' => \$group,
	    'pair=s' => \$pair,
	    'pre=s' => \$pre,
	    'config=s' => \$config,
	    'targets=s' => \$targets,
	    'help' => \$help,
	    'nosnps' => \$nosnps,
	    'removedups' => \$removedups,
	    'abra' => \$abra,
	    'ug|unifiedgenotyper' => \$ug,
 	    'output|out|o=s' => \$output,
 	    'scheduler=s' => \$scheduler,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
	    'species=s' => \$species) or exit(1);

if(!$map || !$group || !$species || !$config || !$scheduler || !$targets || $help){
    print <<HELP;

    USAGE: variants_pipeline.pl -map MAP -group GROUP -pair PAIR -pre PRE -config CONFIG -species SPECIES -scheduler SCHEDULER -targets TARGETS
	* MAP: file listing sample information for processing (REQUIRED)
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED)
	* SPECIES: only hg19 (default: human), mm9, mm10 (default: mouse), hybrid (hg19+mm10) and dm3 currently supported (REQUIRED)
	* TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
	* PRE: output prefix (default: TEMP)
	* OUTPUT: output results directory (default: results)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* -nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed
	* -removedups: remove duplicate reads instead of just marking them
	* -abra: run abra instead of GATK indelrealigner
	* haplotypecaller is default; -ug || -unifiedgenotyper to also make unifiedgenotyper variant calls	
HELP
exit;
}


my $commandLine = reconstructCL();

my $REF_SEQ = '';
my $DB_SNP = '';

my $HG19_FASTA = '';
my $HG19_MM10_HYBRID_FASTA = '';
my $MM9_FASTA = '';
my $MM10_FASTA = '';
my $DM3_FASTA = '';

my $targets_ilist = "$Bin/targets/$targets/$targets\_targets.ilist";
my $baits_ilist = "$Bin/targets/$targets/$targets\_baits.ilist";
my $targets_bed_padded = "$Bin/targets/$targets/$targets\_targets_plus5bp.bed";

my %grouping = ();
my %grouping_samples = ();
my @run_merge_jids = ();
my $PICARD = '';
my $JAVA = '';
my $PERL = '';
my $PYTHON = '';

&verifyConfig($config);
processInputs();

my $curDir = `pwd`;
chomp $curDir;
my $cd = $curDir;
$cd =~ s/\//_/g;

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my %samp_libs_run = ();
my $slr_count = 0;
my %ran_solexa = ();


open(LOG, ">$cd\_variants_pipeline.log") or die "can't write to output log";
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSTARTING VARIANTS PIPELINE FOR $pre\n";
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tCOMMAND LINE: $commandLine\n";


### mkdir -p does not apply the -m permissions to the parent dir; only subdir
`/bin/mkdir -m 775 -p $output`; 
`/bin/mkdir -m 775 -p $output/intFiles`; 
`/bin/mkdir -m 775 -p $output/progress`;
`/bin/mkdir -m 775 -p $output/metrics`;
`/bin/mkdir -m 775 -p $output/metrics/fingerprint`;

my $ran_sol = 0;
alignReads();

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
my $ran_hs = 0;
my $ran_is = 0;
my $ran_as = 0;
my $ran_md = 0;
my $ran_cog = 0;
my @md_jids = ();
my @hs_jids = ();
my @is_jids = ();
my @as_jids = ();
my @cog_jids = ();

processBams();
foreach my $md_jid (@md_jids){
    `$Bin/jobSync $scheduler $md_jid`;
}
my $mdj = join(",", @md_jids);

mergeStats();
generateGroupFile();
callSNPS();

close LOG;

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
    if($config){
	$rCL .= " -config $config";
    }
    if($species){
	$rCL .= " -species $species";
    }
    if($scheduler){
	$rCL .= " -scheduler $scheduler";
    }
    if($output){
	$rCL .= " -output $output";
    }
    if($targets){
	$rCL .= " -targets $targets";
    }
    if($nosnps){
	$rCL .= " -nosnps";
    }
    if($removedups){
	$rCL .= " -removedups";
    }
    if($abra){
	$rCL .= " -abra";
    }
    if($ug){
	$rCL .= " -unifiedgenotyper";
    }
    
    $rCL .= " -priority_project $priority_project";
    
    $rCL .= " -priority_group $priority_group";
    
    if($species){
	$rCL .= " -species $species";
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
	}
	elsif($conf[0] =~ /mutect/i){
	    if(!-e "$conf[1]/muTect.jar"){
		die "CAN'T FIND muTect.jar IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /bedtools/i){
	    if(!-e "$conf[1]/bedtools"){
		die "CAN'T FIND bedtools IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /oncotator/i){
	    if(!-e "$conf[1]/oncotateMaf.sh"){
		die "CAN'T FIND oncotateMaf.sh IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /samtools/i){
	    if(!-e "$conf[1]/samtools"){
		die "CAN'T FIND samtools IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /snpeff/i){
	    if(!-e "$conf[1]/snpEff.jar"){
		die "CAN'T FIND snpEff.jar IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /somaticsniper/i){
	    if(!-e "$conf[1]/bam-somaticsniper"){
		die "CAN'T FIND bam-somaticsniper IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /varscan/i){
	    if(!-e "$conf[1]/VarScan.jar"){
		die "CAN'T FIND VarScan.jar IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /strelka/i){
	    if(!-e "$conf[1]/bin/configureStrelkaWorkflow.pl"){
		die "CAN'T FIND bin/configureStrelkaWorkflow.pl IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /scalpel/i){
	    if(!-e "$conf[1]/scalpel"){
		die "CAN'T FIND scalpel IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /virmid/i){
	    if(!-e "$conf[1]/Virmid.jar"){
		die "CAN'T FIND Virmid.jar IN $conf[1] $!";
	    }
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
	elsif($conf[0] =~ /java/i){
	    if(!-e "$conf[1]/java"){
		die "CAN'T FIND java IN $conf[1] $!";
	    }
	    $JAVA = $conf[1];
	    my $path_tmp = $ENV{'PATH'};
	    $ENV{'PATH'} = "$conf[1]:$path_tmp";
	}
	elsif($conf[0] =~ /hg19_fasta/i){
	    if(!-e "$conf[1]"){
		die "CAN'T FIND $conf[1] $!";
	    }
	    $HG19_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /hg19_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19 WITH PREFIX $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /hg19_mm10_hybrid_fasta/i){
	    if(!-e "$conf[1]"){
		die "CAN'T FIND $conf[1] $!";
	    }
	    $HG19_MM10_HYBRID_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /hg19_mm10_hybrid_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19-MM9 HYBRID WITH PREFIX $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /mm9_fasta/i){
	    if(!-e "$conf[1]"){
		die "CAN'T FIND $conf[1] $!";
	    }
	    $MM9_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /mm9_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM9 WITH PREFIX $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /mm10_fasta/i){
	    if(!-e "$conf[1]"){
		die "CAN'T FIND $conf[1] $!";
	    }
	    $MM10_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /mm10_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM10 WITH PREFIX $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /dm3_fasta/i){
	    if(!-e "$conf[1]"){
		###die "CAN'T FIND $conf[1] $!";
	    }
	    $DM3_FASTA = $conf[1];
	}
 	elsif($conf[0] =~ /dm3_bwa_index/i){
	    if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
		###die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR DM3 WITH PREFIX $conf[1] $!";
	    }
	}
   }
    close CONFIG;
}

sub processInputs {
    if($pre =~ /^\d+/){
	$pre = "s_$pre";
    }
    
    if($species !~ /human|hg19|mouse|mm9|mm10|drosophila|dm3|hybrid/i){
	die "Species must be human (hg19), mouse (mm9/mm10) or drosophila (dm3)";
    }
    
    if($scheduler =~ /lsf/i){
	$scheduler = 'lsf';
    }
    elsif($scheduler =~ /sge/i){
	$scheduler = 'sge';
    }
    
    if($species =~ /human|hg19/i){
	$species = 'hg19';
	$REF_SEQ = "$HG19_FASTA";
	$DB_SNP = "$Bin/data/dbsnp_135.hg19__ReTag.vcf";
    }
    elsif($species =~ /mm9/i){
	$species = 'mm9';
	$REF_SEQ = "$MM9_FASTA";
	$DB_SNP = "";
    }
    elsif($species =~ /mouse|mm10/i){
	$species = 'mm10';
	$REF_SEQ = "$MM10_FASTA";
	$DB_SNP = "";
    }
    elsif($species =~ /hybrid/i){
	$species = 'hybrid';
	$REF_SEQ = "$HG19_MM10_HYBRID_FASTA";
	$DB_SNP = "$Bin/data/dbsnp_135.hg19__ReTag.vcf";
    }
    elsif($species =~ /drosophila|dm3/i){
	$species = 'dm3';
	$REF_SEQ = "$DM3_FASTA";
	$DB_SNP = "";
    }
        
    open(GROUP, "$group") or die "Can't open grouping file $group $!";
    while(<GROUP>){
	chomp;
	
	my @data = split(/\s+/, $_);
	$grouping{$data[1]}{$data[0]} = 1;
	$grouping_samples{$data[0]} = 1;
    }
    close GROUP;
    
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
    }
    
    my %mapping_samples = ();
    open(MA, "$map") or die "Can't open mapping file $map $!";
    while(<MA>){
	chomp;
	
	my @data = split(/\s+/, $_);
	
	$mapping_samples{$data[1]} = 1;
	if(!-d $data[3]){
	    die "$data[3] does not exist\n";
	}
	
	if(!$grouping_samples{$data[1]}){
	    die "grouping file $group missing sample $data[1] found in mapping file $map $!";
	}
	
	if($pair){
	    if(!$pairing_samples{$data[1]}){
		die "pairing file $pair missing sample $data[1] found in mapping file $map $!";
	    }
	}
    }
    close MA;
    
    foreach my $gro (keys %grouping_samples){
	if(!$mapping_samples{$gro}){
	    die "grouping file $group contains sample $gro that isn't in mapping file $map $!";
	}
	
	if($pair){
	    if(!$pairing_samples{$gro}){
		die "grouping file $group contains sample $gro that isn't in pairing file $pair $!";
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
	}
    }
    
    if(!-e $targets_ilist || !-e $baits_ilist || !-e $targets_bed_padded){
	die "$Bin/targets/$targets MUST CONTAIN THE FOLLOWING FILES: $targets\_targets.ilist; $targets\_baits.ilist; $targets\_targets_plus5bp.bed $!";	
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
	
	`/bin/mkdir -m 775 -p $output/intFiles/$data[1]`;
	`/bin/mkdir -m 775 -p $output/intFiles/$data[1]/$data[0]`;
	`/bin/mkdir -m 775 -p $output/intFiles/$data[1]/$data[0]/$data[2]`;
	$samp_libs_run{$data[1]}{$data[0]}{$data[2]} = 1;	
	
	if(!-e "$output/progress/reads_files_$data[1]\_$data[0]\_$data[2]\.done"){
	    $ran_sol = 1;
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
			
			if($data[4] =~ /pe/i){
			    print OUT "$file\t$file_R2\n";
			}
			elsif($data[4] =~ /se/i){
			    print OUT "$file\n";
			}
		    }
		}
	    }
	    close OUT;
	
	    my @currentTime = &getTime();
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSTARTING READS PROCESSING/ALIGNMENT FOR $data[1]\_$data[0]\_$data[2]\n";
	    
	    if($data[4] =~ /pe/i){
		`$Bin/process_reads_pe.pl -file files_$data[1]\_$data[0]\_$data[2] -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\tID:$data[1]\_$data[0]\_$data[2]\_PE\\tPL:Illumina\\tPU:$data[1]\_$data[0]\_$data[2]\\tLB:$data[1]\_$data[0]\\tSM:$data[1]" -species $species -config $config -scheduler $scheduler > files_$data[1]\_$data[0]\_$data[2]\_process_reads_pe.log 2>&1`;
		
		###`/common/sge/bin/lx24-amd64/qsub /home/mpirun/tools/qCMD $Bin/solexa_PE.pl -file files -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\\tID:$data[1]\_$data[0]\_$data[2]\_PE\\\tPL:Illumina\\\tPU:$data[1]\_$data[0]\_$data[2]\\\tLB:$data[1]\_$data[0]\\\tSM:$data[1]" -species $species -config $config $targeted -scheduler $scheduler`;
	    }
	    elsif($data[4] =~ /se/i){
		`$Bin/process_reads_se.pl -file files_$data[1]\_$data[0]\_$data[2] -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\tID:$data[1]\_$data[0]\_$data[2]\_SE\\tPL:Illumina\\tPU:$data[1]\_$data[0]\_$data[2]\\tLB:$data[1]\_$data[0]\\tSM:$data[1]" -species $species -config $config -scheduler $scheduler > files_$data[1]\_$data[0]\_$data[2]\_process_reads.log 2>&1`;
	    }
	    $ran_solexa{$data[1]} = 1;
	    chdir $curDir;
	    push @run_merge_jids, "$pre\_$uID\_$data[1]\_$data[0]\_$data[2]\_MERGE";
	    `/bin/touch $output/progress/reads_files_$data[1]\_$data[0]\_$data[2]\.done`;
	}
	else{
	    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSKIPPING READS PROCESSING/ALIGNMENT FOR $data[1]\_$data[0]\_$data[2]; PREVIOUSLY RAN TO COMPLETION\n";
	}
    }
}

sub processBams {
    my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "1");
    my $additionalParams = Schedule::additionalParams(%addParams);

    foreach my $samp (keys %samp_libs_run){
	my @sBams = ();
	my @lib_merge_samp_jids = ();
	foreach my $lib (keys %{$samp_libs_run{$samp}}){
	    my @lBams = ();
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
		if(!-e "$output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done" || $ran_solexa{$samp}){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$samp\_$lib\_MARKDUPS", job_hold => "$rmj", cpu => "3", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MarkDuplicates $fin OUTPUT=$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam METRICS_FILE=$output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$rmdups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;

		    `/bin/touch $output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done`;
		    push @md_jids, "$pre\_$uID\_$samp\_$lib\_MARKDUPS";
		    $ran_md = 1;
		}
	    }
	    else{
		push @sBams, "I=$output/intFiles/$samp/$lib/$samp\_$lib\.bam";
		if(!-e "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done" || $ran_solexa{$samp}){
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_LIB_MERGE_$samp\_$lib", job_hold => "$rmj", cpu => "24", mem => "90", cluster_out => "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $fin O=$output/intFiles/$samp/$lib/$samp\_$lib\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;

		    `/bin/touch $output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done`;
		    $ran_lb_merge = 1;
		    push @lib_merge_samp_jids, "$pre\_$uID\_LIB_MERGE_$samp\_$lib";
		    sleep(3);
		}
		
		if(!-e "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done" || $ran_lb_merge){
		    my $lb_merge_hold = "";
		    if($ran_lb_merge){
			$lb_merge_hold = "$pre\_$uID\_LIB_MERGE_$samp\_$lib";
		    }		
			
		    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$samp\_$lib\_MARKDUPS", job_hold => "$lb_merge_hold", cpu => "3", mem => "30", cluster_out => "$output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.log");
		    my $standardParams = Schedule::queuing(%stdParams);
		    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MarkDuplicates INPUT=$output/intFiles/$samp/$lib/$samp\_$lib\.bam OUTPUT=$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam METRICS_FILE=$output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$rmdups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;

		    `/bin/touch $output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done`;
		    push @md_jids, "$pre\_$uID\_$samp\_$lib\_MARKDUPS";
		    $ran_md = 1;
		}
	    }
	    
	    push @mdm, "-metrics $output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt";
	}
	
	my $rin = join(" ", @sBams);
	my $bamForStats = '';
	my $ran_samp_merge = 0;	
	my @samp_merge_samp_jids = ();
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
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $rin O=$output/intFiles/$samp/$samp\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_SAMP_MERGE_$samp\.done`;
		push @samp_merge_samp_jids, "$pre\_$uID\_SAMP_MERGE_$samp";
		$ran_samp_merge = 1;
	    }
	    $bamForStats = "$output/intFiles/$samp/$samp\.bam";
	}
	
	my $smsj = join(",", @samp_merge_samp_jids);
	if(!-e "$output/progress/$pre\_$uID\_HS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_HS_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_HS_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CalculateHsMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_HsMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE BAIT_INTERVALS=$baits_ilist BAIT_SET_NAME=$targets TARGET_INTERVALS=$targets_ilist VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_HS_$samp\.done`;
	    push @hs_jids, "$pre\_$uID\_HS_METRICS\_$samp";
	    $ran_hs = 1;
	}
	push @hsm, "-metrics $output/intFiles/$pre\_HsMetrics_$samp\.txt";
	
	if(!-e "$output/progress/$pre\_$uID\_IS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_IS_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_IS_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectInsertSizeMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_InsertSizeMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE HISTOGRAM_FILE=$output/intFiles/$pre\_InsertSizeMetrics_Histogram_$samp\.txt VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_IS_$samp\.done`;
	    push @is_jids, "$pre\_$uID\_IS_METRICS\_$samp";
	    $ran_is = 1;
	}
	push @ism, "-metrics $output/intFiles/$pre\_InsertSizeMetrics_$samp\.txt";
	
	if(!-e "$output/progress/$pre\_$uID\_AS_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_AS_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_AS_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_AlignmentSummaryMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_AS_$samp\.done`;
	    push @as_jids, "$pre\_$uID\_AS_METRICS\_$samp";
	    $ran_as = 1;
	}
	push @asm, "-metrics $output/intFiles/$pre\_AlignmentSummaryMetrics_$samp\.txt";

	if(!-e "$output/progress/$pre\_$uID\_COG_$samp\.done" || $ran_solexa{$samp} || $ran_samp_merge){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_COG_METRICS\_$samp", job_hold => "$rmj,$smsj,$lmsj", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_COG_$samp\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectOxoGMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_CollectOxoGMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ DB_SNP=$DB_SNP VALIDATION_STRINGENCY=LENIENT`;
	    `/bin/touch $output/progress/$pre\_$uID\_COG_$samp\.done`;
	    push @cog_jids, "$pre\_$uID\_COG_METRICS\_$samp";
	    $ran_cog = 1;
	}
	push @cogm, "-metrics $output/intFiles/$pre\_CollectOxoGMetrics_$samp\.txt";
    }
}

sub mergeStats {
    my @qcpdf_jids = ();
    my $ran_merge = 0;
    my $ran_merge_ism = 0;
    my %addParams = (scheduler => "$scheduler", runtime => "50", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "1");
    my $additionalParams = Schedule::additionalParams(%addParams);

    my $mdfiles = join(" ", @mdm);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_MDM.done" || $ran_md){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_MDM", job_hold => "$mdj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_MDM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $mdfiles ">$output/metrics/$pre\_markDuplicatesMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_MDM.done`;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_MDM";
	$ran_merge = 1;
    }
    
    my $hsfiles = join(" ", @hsm);
    my $hsj = join(",", @hs_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_HSM.done" || $ran_hs){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_HSM", job_hold => "$hsj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_HSM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $hsfiles ">$output/metrics/$pre\_HsMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_HSM.done`;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_HSM";
	$ran_merge = 1;
    }
    
    my $isfiles = join(" ", @ism);
    my $isj = join(",", @is_jids);
    if(!-e "$output/progress/$pre\_$uID\_MERGE_ISM.done" || $ran_is){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_ISM", job_hold => "$isj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_ISM.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/mergePicardMetrics.pl $isfiles ">$output/metrics/$pre\_InsertSizeMetrics.txt"`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_ISM.done`;
	$ran_merge_ism = 1;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_ISM";
	$ran_merge = 1;
    }

    if(!-e "$output/progress/$pre\_$uID\_ISM_MATRIX.done" || $ran_merge_ism){    
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ISM_MATRIX", job_hold => "$pre\_$uID\_MERGE_ISM", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_ISM_MATRIX.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeInsertSizeHistograms.py $output/intFiles '*InsertSizeMetrics_*.txt' $output/metrics/$pre\_InsertSizeMetrics_Histograms.txt`;
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

    if(!-e "$output/progress/$pre\_$uID\_MERGE_CAS.done" || $ran_sol){
	###my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CAS", job_hold => "$pre\_$uID\_CUTADAPT*", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CAS.log");
	### NOTE: all of the cutadapt jobs should be done because of the job sync on merge jobs
	###       so don't need this to hold on cutadapt
	###       this hold causes issues on lsf because of the *
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_MERGE_CAS", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_MERGE_CAS.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/mergeCutAdaptStats.py . '*CUTADAPT_STATS.txt' $output/metrics/$pre\_CutAdaptStats.txt`;
	`/bin/touch $output/progress/$pre\_$uID\_MERGE_CAS.done`;
	push @qcpdf_jids, "$pre\_$uID\_MERGE_CAS";
    }

    if(!-e "$output/progress/$pre\_$uID\_QCPDF.done" || $ran_merge){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_QCPDF", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_QCPDF.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PERL/perl $Bin/qc/qcPDF.pl -path $output/metrics -pre $pre -config $config`;
	`/bin/touch $output/progress/$pre\_$uID\_QCPDF.done`;
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
		if(!-e "$output/intFiles/$sample/$lib/$sample\_$lib\_MD.bam"){
		    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tCAN'T LOCATE $output/intFiles/$sample/$lib/$sample\_$lib\_MD.bam...EXITING VARIANTS PIPELINE FOR PROJECT $pre";
		    die;
		}
		push @groupings, "$output/intFiles/$sample/$lib/$sample\_$lib\_MD.bam";
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
    
    my $run_ug = '';
    if($ug){
	$run_ug = "-ug";
    }

    my $run_abra = '';
    if($abra){
	$run_abra = "-abra";
    }
    
    my $run_step1 = '';
    if($ran_md){
	$run_step1 = '-step1';
    }

    my @currentTime = &getTime();
    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[3]\/$currentTime[4]\/$currentTime[5]\tSNP CALLING RUNNING\n";
    my %addParams = (scheduler => "$scheduler", runtime => "50", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "1");
    my $additionalParams = Schedule::additionalParams(%addParams);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_SNP_PIPE", job_hold => "$mdj", cpu => "1", mem => "1", cluster_out => "$output/progress/$pre\_$uID\_SNP_PIPE.log");
    my $standardParams = Schedule::queuing(%stdParams);
    
    if($species =~ /hg19/i){
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_hg19.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_abra $run_step1`;
    }
    elsif($species =~ /hybrid/i){
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_hg19.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -hybrid -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_abra $run_step1`;
    }
    elsif($species =~ /mm9|mm10/i){
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_mouse.pl -pre $pre $paired -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -targets $targets $run_ug -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group -species $species $run_abra $run_step1`;
    }
    elsif($species =~ /dm3/i){
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/process_alignments_dm3.pl -pre $pre -group $group -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps $run_ug -output $output -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group $run_abra $run_step1`;
    }
}
