#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

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

my ($map, $group, $pair, $pre, $config, $help, $nosnps, $removedups, $species, $target_bed, $target_ilist, $bait, $ug, $output);

$pre = 'TEMP';
$output = "results";
GetOptions ('map=s' => \$map,
	    'group=s' => \$group,
	    'pair=s' => \$pair,
	    'pre=s' => \$pre,
	    'config=s' => \$config,
	    'target_bed=s' => \$target_bed,
	    'target_ilist=s' => \$target_ilist,
	    'bait=s' => \$bait,
	    'help' => \$help,
	    'nosnps' => \$nosnps,
	    'removedups' => \$removedups,
	    'ug|unifiedgenotyper' => \$ug,
 	    'output|out|o=s' => \$output,
	    'species=s' => \$species) or exit(1);

if(!$map || !$group || !$species || !$config || $help){
    print <<HELP;

    USAGE: exome_pipeline.pl -map MAP -group GROUP -pair PAIR -pre PRE -config CONFIG -species SPECIES
	* MAP: file listing sample information for processing (REQUIRED)
	* GROUP: file listing grouping of samples for realign/recal steps (REQUIRED)
	* PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
	* PRE: output prefix (default: TEMP)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SPECIES: only hg19, mm9, hybrid (hg19+mm10) and dm3 currently supported (REQUIRED)
	* OUTPUT: output results directory (default: results)
	* -nosnps: if no snps to be called; e.g. when only indelrealigned/recalibrated bams needed
	* -removedups: remove duplicate reads instead of just marking them
	* TARGET_BED: target bed file (default: hg19__MegaGene__v131104.bed)
	* TARGET_ILIST: target ilist file (default:: hg19__MegaGene__v131104.ilist)
	* BAIT: bait ilist file (default: hg19__MegaGene__v131104)
	* haplotypecaller is default; -ug || -unifiedgenotyper to also make unifiedgenotyper variant calls	
HELP
exit;
}

### dumb reconstruction of command line
my $commandLine = "$Bin/exome_pipeline.pl";
if($pre){
    $commandLine .= " -pre $pre";
}
if($map){
    $commandLine .= " -map $map";
}
if($group){
    $commandLine .= " -group $group";
}
if($pair){
    $commandLine .= " -pair $pair";
}
if($config){
    $commandLine .= " -config $config";
}
if($species){
    $commandLine .= " -species $species";
}
if($output){
    $commandLine .= " -output $output";
}
if($target_bed){
    $commandLine .= " -target_bed $target_bed";
}
if($target_ilist){
    $commandLine .= " -target_ilist $target_ilist";
}
if($bait){
    $commandLine .= " -bait $bait";
}
if($nosnps){
    $commandLine .= " -nosnps";
}
if($removedups){
    $commandLine .= " -removedups";
}
if($ug){
    $commandLine .= " -unifiedgenotyper";
}
if($output){
    $commandLine .= " -output $output";
}
if($species){
    $commandLine .= " -species $species";
}

my $numArgs = $#ARGV + 1;
foreach my $argnum (0 .. $#ARGV) {
    $commandLine .= " $ARGV[$argnum]";
}

if($pre =~ /^\d+/){
    $pre = "s_$pre";
}

if($species !~ /human|hg19|mouse|mm9|drosophila|dm3|hybrid/i){
    die "Species must be human (hg19), mouse (mm9) or drosophila (dm3)";
}

my $REF_SEQ = '';
my $DB_SNP = '';
if($species =~ /human|hg19/i){
    $species = 'hg19';
    $REF_SEQ = '/ifs/data/bio/assemblies/H.sapiens/hg19/hg19.fasta';
    $DB_SNP = "$Bin/data/dbsnp_135.hg19__ReTag.vcf";
}
elsif($species =~ /mouse|mm9/i){
    $species = 'mm9';
    $REF_SEQ = '/ifs/data/bio/assemblies/M.musculus/mm9/mm9.fasta';
}
elsif($species =~ /hybrid/i){
    $species = 'hybrid';
    $REF_SEQ = '/ifs/data/bio/assemblies/hybrid_H.sapiens_M.musculus/hybrid_hg19_mm10/hybrid_hg19_mm10.fasta';
}
elsif($species =~ /drosophila|dm3/i){
    $species = 'dm3';
    $REF_SEQ = '/ifs/data/bio/assemblies/D.melanogaster/dm3/dm3.fasta';
}

my $PICARD = '';

### Check that all programs are available
&verifyConfig($config);

my %grouping = ();
my %grouping_samples = ();
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
	die "grouping file $group contains samplee $gro that isn't in mapping file $map $!";
    }

    if($pair){
	if(!$pairing_samples{$gro}){
	    die "grouping file $group contains samplee $gro that isn't in pairing file $pair $!";
	}
    }
}

if($pair){
    foreach my $pai (keys %pairing_samples){
	if(!$mapping_samples{$pai}){
	    die "pairing file $pair contains samplee $pai that isn't in mapping file $map $!";
	}

	if(!$grouping_samples{$pai}){
	    die "pairing file $pair contains samplee $pai that isn't in grouping file $group $!";
	}
    }
}

if(($target_bed && !$target_ilist) || ($target_ilist && !$target_bed)){
    die "IF SPECIFYING TARGETS MUST SPECIFY target_bed AND target_list FILES $!";
}
elsif(!$target_bed && !$target_ilist){
    if($species =~ /human|hg19/i){
	$target_bed = "$Bin/targets/hg19__MegaGene__v131104.bed";
	$target_ilist = "$Bin/targets/hg19__MegaGene__v131104.ilist";
    }
    elsif($species =~ /hybrid/i){
	$target_bed = "$Bin/targets/hg19_mm10_hybrid_MegaGene__v131104.bed";
	$target_ilist = "$Bin/targets/hg19_mm10_hybrid_MegaGene__v131104.ilist";
    }
    elsif($species =~ /mouse|mm9/i){
	$target_bed = "$Bin/targets/mm9__MegaGene__v120228_MERGE.bed";
	$target_ilist = "$Bin/targets/mm9__MegaGene__v120228_MERGE.ilist";
    }
}
else{
    if(!-e $target_bed || !-e $target_ilist){
	die "CAN'T LOCATE $target_bed || $target_ilist. MUST HAVE BOTH $!";
    }
    else{
	if($target_bed !~ /\.bed/i || $target_ilist !~ /\.ilist|\.interval_list/i){
	    die "TARGET BED MUST END WITH .bed AND TARGET INTERVAL MUST END IN .ilist OR .interval_list EXTENSION $!";
	}
    }
}

if($bait){
    if($bait !~ /\.ilist|\.interval_list/i || !-e $bait){
	die "BAIT FILE MUST END WITH .ilist OR .interval_list EXTENSION OR bait file $bait COULD NOT BE FOUND $!";
    }
}
else{
    $bait = $target_ilist;
}

my $curDir = `pwd`;
chomp $curDir;
my $cd = $curDir;
$cd =~ s/\//_/g;

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my %samp_libs_run = ();
my $slr_count = 0;
my %ran_solexa = ();


open(LOG, ">$cd\_exome_pipeline.log") or die "can't write to output log";
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSTARTING EXOME PIPELINE FOR $pre\n";
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tCOMMAND LINE: $commandLine\n";

### mkdir -p does not apply the -m permissions to the parent dir; only subdir
`/bin/mkdir -m 775 -p $output`; 
`/bin/mkdir -m 775 -p $output/intFiles`; 
`/bin/mkdir -m 775 -p $output/progress`;
`/bin/mkdir -m 775 -p $output/metrics`;
`/bin/mkdir -m 775 -p $output/metrics/fingerprint`;

my $ran_sol = 0;
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
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tWARNING: $data[1]\t$data[0]\t$data[2] ISN'T UNIQUE; ";
	$data[2] = "$data[2]\_$slr_count";
	print LOG "WRITING INSTEAD TO $data[1]\/$data[0]\/$data[2]\n";
    }

    `/bin/mkdir -m 775 -p $output/intFiles/$data[1]`;
    `/bin/mkdir -m 775 -p $output/intFiles/$data[1]/$data[0]`;
    `/bin/mkdir -m 775 -p $output/intFiles/$data[1]/$data[0]/$data[2]`;
    $samp_libs_run{$data[1]}{$data[0]}{$data[2]} = 1;


    if(!-e "$output/progress/solexa_files_$data[1]\_$data[0]\_$data[2]\.done"){
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
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSTARTING READS PROCESSING/ALIGNMENT FOR $data[1]\_$data[0]\_$data[2]\n";
	
	if($data[4] =~ /pe/i){
	    `$Bin/solexa_PE.pl -file files_$data[1]\_$data[0]\_$data[2] -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\tID:$data[1]\_$data[0]\_$data[2]\_PE\\tPL:Illumina\\tPU:$data[1]\_$data[0]\_$data[2]\\tLB:$data[1]\_$data[0]\\tSM:$data[1]" -species $species -config $config > files_$data[1]\_$data[0]\_$data[2]\_solexa_PE.log 2>&1`;
	    
	    ###`/common/sge/bin/lx24-amd64/qsub /home/mpirun/tools/qCMD $Bin/solexa_PE.pl -file files -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\\tID:$data[1]\_$data[0]\_$data[2]\_PE\\\tPL:Illumina\\\tPU:$data[1]\_$data[0]\_$data[2]\\\tLB:$data[1]\_$data[0]\\\tSM:$data[1]" -species $species -config $config $targeted`;
	}
	elsif($data[4] =~ /se/i){
	    `$Bin/solexa_SE.pl -file files_$data[1]\_$data[0]\_$data[2] -pre $pre -run $data[1]\_$data[0]\_$data[2] -readgroup "\@RG\\tID:$data[1]\_$data[0]\_$data[2]\_SE\\tPL:Illumina\\tPU:$data[1]\_$data[0]\_$data[2]\\tLB:$data[1]\_$data[0]\\tSM:$data[1]" -species $species -config $config > files_$data[1]\_$data[0]\_$data[2]\_solexa_SE.log 2>&1`;
	}
	$ran_solexa{$data[1]} = 1;
	chdir $curDir;
	`/bin/touch $output/progress/solexa_files_$data[1]\_$data[0]\_$data[2]\.done`;
    }
    else{
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSKIPPING READS PROCESSING/ALIGNMENT FOR $data[1]\_$data[0]\_$data[2]; PREVIOUSLY RAN TO COMPLETION\n";
    }
}

my $rmdups = 'false';
if($removedups){
    $rmdups = 'true';
}

`$Bin/qSYNC $pre\_$uID\_RUN_MERGE`;
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\t$pre\tREADS PROCESSING/ALIGNMENTS FOR ALL SAMPLES COMPLETED\n";
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\t$pre\_$uID\_MARKDUPS RUNNING\n";

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

foreach my $samp (keys %samp_libs_run){
    my @sBams = ();
    foreach my $lib (keys %{$samp_libs_run{$samp}}){
	my @lBams = ();
	foreach my $run (keys %{$samp_libs_run{$samp}{$lib}}){
	    if(!-e "$output/intFiles/$samp/$lib/$run/$samp\_$lib\_$run\.bam"){
		my @currentTime = &getTime();
		print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tCan't locate $output/intFiles/$samp/$lib/$run/$samp\_$lib\_$run\.bam ...exiting exome pipeline for $pre";
		die;
	    }
	    push @lBams, "I=$output/intFiles/$samp/$lib/$run/$samp\_$lib\_$run\.bam";
	    push @sBams, "I=$output/intFiles/$samp/$lib/$run/$samp\_$lib\_$run\.bam";
	}

	my $fin = join(" ", @lBams);
	my $ran_lb_merge = 0;
	### NOTE: DOESN'T CHECK TO SEE IF SOLEXA_PE OR SOLEXA_SE WAS RUN
	###       WILL AUTOMATICALLY RUN IF IT DOESN'T FIND THE .done FILE FOR THE SAMPLE LIB
	###       THIS IS TO AVOID RUNNING THIS CPU COSTLY STEP IF JUST 1 OR SO SAMPLE
	###       IN A COHORT NEEDED TO HAVE THEIR READS REPROCESSED
	if(scalar(@lBams) == 1){
	    if(!-e "$output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done" || $ran_solexa{$samp}){
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MARKDUPS -hold_jid $pre\_$uID\_RUN_MERGE -pe alloc 3 -l virtual_free=10G $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MarkDuplicates $fin OUTPUT=$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam METRICS_FILE=$output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$rmdups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done`;
		$ran_md = 1;
	    }
	}
	else{
	    if(!-e "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done" || $ran_solexa{$samp}){
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_LIB_MERGE_$samp\_$lib -hold_jid $pre\_$uID\_RUN_MERGE -pe alloc 24 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $fin O=$output/intFiles/$samp/$lib/$samp\_$lib\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done`;
		$ran_lb_merge = 1;
                sleep(3);
	    }

	    if(!-e "$output/progress/$pre\_$uID\_LIB_MERGE_$samp\_$lib\.done" || $ran_lb_merge){	
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MARKDUPS -hold_jid $pre\_$uID\_LIB_MERGE\_$samp\_$lib -pe alloc 3 -l virtual_free=10G $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MarkDuplicates INPUT=$output/intFiles/$samp/$lib/$samp\_$lib\.bam OUTPUT=$output/intFiles/$samp/$lib/$samp\_$lib\_MD.bam METRICS_FILE=$output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$rmdups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;
		`/bin/touch $output/progress/$pre\_$uID\_$samp\_$lib\_MARKDUPS.done`;
		$ran_md = 1;
	    }
	}
	push @mdm, "-metrics $output/intFiles/$samp/$lib/$samp\_$lib\_markDuplicatesMetrics.txt";
    }

    my $rin = join(" ", @sBams);
    my $bamForStats = '';
    if(scalar(@sBams) == 1){
	my @bname = split(/=/, $sBams[0]);
	$bamForStats = "$bname[1]";
    }
    else{
	if(!-e "$output/progress/$pre\_$uID\_SAMP_MERGE_$samp\.done"){
	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SAMP_MERGE_$samp -hold_jid $pre\_$uID\_RUN_MERGE -pe alloc 24 -l virtual_free=3G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $rin O=$output/intFiles/$samp/$samp\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
	    `/bin/touch $output/progress/$pre\_$uID\_SAMP_MERGE_$samp\.done`;
	}
	$bamForStats = "$output/intFiles/$samp/$samp\.bam";
    }

    sleep(3);

    if(!-e "$output/progress/$pre\_$uID\_HS_$samp\.done"){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_HS_METRICS -hold_jid $pre\_$uID\_RUN_MERGE,$pre\_$uID\_SAMP_MERGE_$samp -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CalculateHsMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_HsMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE BAIT_INTERVALS=$bait TARGET_INTERVALS=$target_ilist VALIDATION_STRINGENCY=LENIENT`;
	`/bin/touch $output/progress/$pre\_$uID\_HS_$samp\.done`;
	$ran_hs = 1;
    }
    push @hsm, "-metrics $output/intFiles/$pre\_HsMetrics_$samp\.txt";

    if(!-e "$output/progress/$pre\_$uID\_IS_$samp\.done"){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_IS_METRICS -hold_jid $pre\_$uID\_RUN_MERGE,$pre\_$uID\_SAMP_MERGE_$samp -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectInsertSizeMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_InsertSizeMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE HISTOGRAM_FILE=$output/intFiles/$pre\_InsertSizeMetrics_Histogram_$samp\.txt VALIDATION_STRINGENCY=LENIENT`;
	`/bin/touch $output/progress/$pre\_$uID\_IS_$samp\.done`;
	$ran_is = 1;
    }
    push @ism, "-metrics $output/intFiles/$pre\_InsertSizeMetrics_$samp\.txt";

    if(!-e "$output/progress/$pre\_$uID\_AS_$samp\.done"){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_AS_METRICS -hold_jid $pre\_$uID\_RUN_MERGE,$pre\_$uID\_SAMP_MERGE_$samp -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_AlignmentSummaryMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT`;
	`/bin/touch $output/progress/$pre\_$uID\_AS_$samp\.done`;
	$ran_as = 1;
    }
    push @asm, "-metrics $output/intFiles/$pre\_AlignmentSummaryMetrics_$samp\.txt";

    if(!-e "$output/progress/$pre\_$uID\_COG_$samp\.done"){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_COG_METRICS -hold_jid $pre\_$uID\_RUN_MERGE,$pre\_$uID\_SAMP_MERGE_$samp -pe alloc 1 -l virtual_free=10G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar CollectOxoGMetrics INPUT=$bamForStats OUTPUT=$output/intFiles/$pre\_CollectOxoGMetrics_$samp\.txt REFERENCE_SEQUENCE=$REF_SEQ DB_SNP=$DB_SNP VALIDATION_STRINGENCY=LENIENT`;
	`/bin/touch $output/progress/$pre\_$uID\_COG_$samp\.done`;
	$ran_cog = 1;
    }
    push @cogm, "-metrics $output/intFiles/$pre\_CollectOxoGMetrics_$samp\.txt";
}
sleep(3);

my $mdfiles = join(" ", @mdm);
if(!-e "$output/progress/$pre\_$uID\_MERGE_MDM.done" || $ran_md){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_MDM -hold_jid $pre\_$uID\_MARKDUPS -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergePicardMetrics.pl $mdfiles ">$output/metrics/$pre\_markDuplicatesMetrics.txt"`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_MDM.done`;
}

my $hsfiles = join(" ", @hsm);
if(!-e "$output/progress/$pre\_$uID\_MERGE_HSM.done" || $ran_hs){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_HSM -hold_jid $pre\_$uID\_HS_METRICS -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergePicardMetrics.pl $hsfiles ">$output/metrics/$pre\_HsMetrics.txt"`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_HSM.done`;
}

my $isfiles = join(" ", @ism);
if(!-e "$output/progress/$pre\_$uID\_MERGE_ISM.done" || $ran_is){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_ISM -hold_jid $pre\_$uID\_IS_METRICS -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergePicardMetrics.pl $isfiles ">$output/metrics/$pre\_InsertSizeMetrics.txt"`;
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_ISM_MATRIX -hold_jid $pre\_$uID\_MERGE_ISM -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergeInsertSizeHistograms.py $output/intFiles '*InsertSizeMetrics_*.txt' $output/metrics/$pre\_InsertSizeMetrics_Histograms.txt`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_ISM.done`;
}

my $asfiles = join(" ", @asm);
if(!-e "$output/progress/$pre\_$uID\_MERGE_ASM.done" || $ran_as){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_ASM -hold_jid $pre\_$uID\_AS_METRICS -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergePicardMetrics.pl $asfiles ">$output/metrics/$pre\_AlignmentSummaryMetrics.txt"`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_ASM.done`;
}

my $cogfiles = join(" ", @cogm);
if(!-e "$output/progress/$pre\_$uID\_MERGE_COGM.done" || $ran_cog){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_COGM -hold_jid $pre\_$uID\_COG_METRICS -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergePicardMetrics.pl $cogfiles ">$output/metrics/$pre\_CollectOxoGMetrics.txt"`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_COGM.done`;
}

if(!-e "$output/progress/$pre\_$uID\_MERGE_CAS.done" || $ran_sol){
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_CAS -hold_jid '$pre\_$uID\_CUTADAPT*' -pe alloc 1 -l virtual_free=1G -q lau.q,lcg.q,nce.q $Bin/qCMD $Bin/mergeCutAdaptStats.py . '*CUTADAPT_STATS.txt' $output/metrics/$pre\_CutAdaptStats.txt`;
    `/bin/touch $output/progress/$pre\_$uID\_MERGE_CAS.done`;
}

`$Bin/qSYNC $pre\_$uID\_MARKDUPS`;
my @currentTime = &getTime();
print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\t$pre\_$uID\_MARKDUPS COMPLETED\n";

open(GR, ">$output/intFiles/$pre\_MDbams_groupings.txt") or die "Can't write $output/intFiles/$pre\_MDbams_groupings.txt $!";
foreach my $grou (keys %grouping){
    my @groupings = ();
    foreach my $sample (keys %{$grouping{$grou}}){
	foreach my $lib (keys %{$samp_libs_run{$sample}}){
	    if(!-e "$output/intFiles/$sample/$lib/$sample\_$lib\_MD.bam"){
		print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tCAN'T LOCATE $output/intFiles/$sample/$lib/$sample\_$lib\_MD.bam...EXITING EXOME PIPELINE FOR $pre";
		die;
	    }
	    push @groupings, "$output/intFiles/$sample/$lib/$sample\_$lib\_MD.bam";
	}
    }
    
    if(scalar(@groupings) < 1){
	my @currentTime = &getTime();
	print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tCAN'T LOCATE ANY MARKDUP BAMS FOR GROUP $grou ...EXITING EXOME PIPELINE FOR $pre";
	die;
    }

    my $gro = join("\t", @groupings);
    print GR "$gro\n";
}
close GR;

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

my @currentTime = &getTime();
if(!-e "$output/progress/$pre\_$uID\_SNP_PIPE.done" || $ran_md){
    print LOG "$currentTime[2]:$currentTime[1]:$currentTime[0], $currentTime[5]\/$currentTime[4]\/$currentTime[3]\tSNP CALLING RUNNING\n";

    if($species =~ /hg19/i){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SN_HG19 -hold_jid $pre\_$uID\_MARKDUPS -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/snp_pipeline_HG19.pl -pre $pre $paired -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -target_bed $target_bed $run_ug -output $output`;
    }
    elsif($species =~ /hybrid/i){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SN_HG19 -hold_jid $pre\_$uID\_MARKDUPS -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/snp_pipeline_HG19.pl -pre $pre $paired -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -target_bed $target_bed $run_ug -hybrid -output $output`;
    }
    elsif($species =~ /mm9/i){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SN_MM9 -hold_jid $pre\_$uID\_MARKDUPS -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/snp_pipeline_MM9.pl -pre $pre $paired -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps -target_bed $target_bed $run_ug -output $output`;
    }
    elsif($species =~ /dm3/i){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_SN_DM3 -hold_jid $pre\_$uID\_MARKDUPS -pe alloc 1 -l virtual_free=1G $Bin/qCMD $Bin/snp_pipeline_DM3.pl -pre $pre -bamgroup $output/intFiles/$pre\_MDbams_groupings.txt -config $config $callSnps $run_ug -output $output`;
    }
    `/bin/touch $output/progress/$pre\_$uID\_SNP_PIPE.done`;
}

close LOG;

sub verifyConfig{
    my $paths = shift;

    open(CONFIG, "$paths") || die "Can't open config file $paths $!";
    while(<CONFIG>){
	chomp;
	
	my @conf = split(/\s+/, $_);	
	if($conf[0] =~ /cutadapt/i){
	    if(!-e "$conf[1]/bin/cutadapt"){
		die "CAN'T FIND bin/cutadapt IN $conf[1] $!";
	    }
	}
	elsif($conf[0] =~ /bwa/i){
	    if(!-e "$conf[1]/bwa"){
		die "CAN'T FIND bwa IN $conf[1] $!";
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
	    if(!-e "$conf[1]/bin/bedtools"){
		die "CAN'T FIND bin/bedtools IN $conf[1] $!";
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
    }
    close CONFIG;
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
