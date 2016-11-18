#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Schedule;
use Cluster;

#### new casava splits fastqs into batches of 4m reads
### this will process each fastq
### and create a single bam
### that can go trhough the second half of pipeline

### takes in a file that list of gzipped read pairs
## e.g.LID46438_NoIndex_L001_R1_001.fastq.gz	LID46438_NoIndex_L001_R2_001.fastq.gz
##     LID46438_NoIndex_L001_R1_002.fastq.gz	LID46438_NoIndex_L001_R2_002.fastq.gz
##     LID46438_NoIndex_L001_R1_003.fastq.gz	LID46438_NoIndex_L001_R2_003.fastq.gz
##     LID46438_NoIndex_L001_R1_004.fastq.gz	LID46438_NoIndex_L001_R2_004.fastq.gz

### '@RG\tID:TEST_SFT20_PE\tPL:Illumina\tPU:TEST_SFT20\tLB:TEST_SFT20\tSM:20'

### ASSUMES THAT YOU HAVE /scratch/$uID ON ALL NODES SET UP

my ($file, $species, $config, $scheduler, $help, $priority_project, $priority_group, $noClip, $noTrim, $defaultBWA);
my $pre = 'TEMP';
my $run = 'TEMP_RUN';
my $readgroup = "\@RG\tID:TEMP_ID\tPL:Illumina\tPU:TEMP_PU\tLB:TEMP_LB\tSM:TEMP_SM";
my $r1adaptor = 'AGATCGGAAGAGCACACGTCT';
my $r2adaptor = 'AGATCGGAAGAGCGTCGTGTA';
my $bqTrim = '3';

my $uID = `/usr/bin/id -u -n`;
chomp $uID;
my $tempdir = "/scratch/$uID";

GetOptions ('file=s' => \$file,
	    'pre=s' => \$pre,
	    'species=s' => \$species,
	    'run=s' => \$run,
	    'config=s' => \$config,
 	    'scheduler=s' => \$scheduler,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
	    'r1adaptor=s' => \$r1adaptor,
	    'r2adaptor=s' => \$r2adaptor,
	    'noClip' => \$noClip,
	    'bqTrim=s' => \$bqTrim,
	    'noTrim' => \$noTrim,
	    'defaultbwa|defaultBWA' => \$defaultBWA,
 	    'tempdir=s' => \$tempdir,
	    'readgroup=s' => \$readgroup) or exit(1);

if(!$file || !$config || !$species || !$scheduler || $help){
    print <<HELP;

    USAGE: process_reads_pe.pl -file FILE -pre PRE -species SPECIES -config CONFIG -run RUN -readgroup READGROUP -scheduler SCHEDULER
	* FILE: file listing read pairs to process (REQUIRED)
	* PRE: output prefix (default: TEMP)
	* SPECIES: b37, mm9, mm10, mm10_custom, species_custom, and dm3 currently supported (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* RUN: RUN IDENTIFIER (default: TEMP_RUN)
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* R1ADAPTOR to specify R1 adaptor sequence (default: AGATCGGAAGAGCACACGTCT)
	* R2ADAPTOR to specify R1 adaptor sequence (default: AGATCGGAAGAGCGTCGTGTA)
	* NOCLIP: no clipping of adaptor sequences
	* BASEQTRIM: base quality to trim in reads (default: 3)
	* NOTRIM: no quality trimming of reads
	* DEFAULTBWA: runs bwa with default parameters; useful for chipseq when -PM causes issues for peakcallers like MACS; -PM is default otherwise
	* TEMPDIR:  temp directory (default: /scratch/$uID)
	* READGROUP: string containing information for ID, PL, PU, LB, and SM e.g. '\@RG\\tID:s_CH_20__1_LOLA_1018_PE\\tPL:Illumina\\tPU:s_CH_20__1_LOLA_1018\\tLB:s_CH_20__1\\t SM:s_CH_20'; WARNING: INCOMPLETE INFORMATION WILL CAUSE GATK ISSSUES (default: '\@RG\\tID:TEMP_ID\\tPL:Illumina\\tPU:TEMP_PU\\tLB:TEMP_LB\\tSM:TEMP_SM')
HELP
exit;
}

my @raw = ();
my $count = 0;
$readgroup =~ s/\s+/\\t/g;

my $bwaDB = '';
my $B37_BWA_INDEX = '';
my $HG19_BWA_INDEX = '';
my $B37_MM10_HYBRID_BWA_INDEX = '';
my $MM9_BWA_INDEX = '';
my $MM10_BWA_INDEX = '';
my $MM10_CUSTOM_BWA_INDEX = '';
my $SPECIES_CUSTOM_BWA_INDEX = '';
my $DM3_BWA_INDEX = '';
my $CUTADAPT = '';
my $BWA = '';
my $PICARD = '';
my $JAVA = '';
my $PYTHON = '';
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config SO USING DEFAULT SETTINGS $!";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /cutadapt/i){
	if(!-e "$conf[1]/cutadapt"){
	    die "CAN'T FIND cutadapt IN $conf[1] $!";
	}
	$CUTADAPT = $conf[1];
    }
    elsif($conf[0] =~ /^bwa/i){
	if(!-e "$conf[1]/bwa"){
	    die "CAN'T FIND bwa IN $conf[1] $!";
	}
	$BWA = $conf[1];
    }
    elsif($conf[0] =~ /picard/i){
	if(!-e "$conf[1]/picard.jar"){
	    die "CAN'T FIND picard.jar IN $conf[1] $!";
	}
	$PICARD = $conf[1];
    }
    elsif($conf[0] =~ /java/i){
	if(!-e "$conf[1]/java"){
	    die "CAN'T FIND java IN $conf[1] $!";
	}
	$JAVA = $conf[1];
    }
    elsif($conf[0] =~ /python/i){
	if(!-e "$conf[1]/python"){
	    die "CAN'T FIND python IN $conf[1] $!";
	}
	$PYTHON = $conf[1];
    }
    elsif($conf[0] =~ /b37_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /^b37$/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR B37 WITH PREFIX $conf[1] $!";
	    }
	}
	$B37_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /hg19_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /^hg19$/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19 WITH PREFIX $conf[1] $!";
	    }
	}
	$HG19_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /b37_mm10_hybrid_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /hybrid/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19-MM10 HYBRID WITH PREFIX $conf[1] $!";
	    }
	}
	$B37_MM10_HYBRID_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /mm9_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /mm9/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM9 WITH PREFIX $conf[1] $!";
	    }
	}
	$MM9_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /mm10_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /^mm10$/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM10 WITH PREFIX $conf[1] $!";
	    }
	}
	$MM10_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /mm10_custom_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /mm10_custom/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR MM10 WITH PREFIX $conf[1] $!";
	    }
	}
	$MM10_CUSTOM_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /species_custom_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /species_custom/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR SPECIES CUSTOM WITH PREFIX $conf[1] $!";
	    }
	}
	$SPECIES_CUSTOM_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /dm3_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    if($species =~ /dm3/i){
		die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR DM3 WITH PREFIX $conf[1] $!";
	    }
	}
 	$DM3_BWA_INDEX = $conf[1];
    }
}
close CONFIG;

if($species =~/^b37$|human/i){
    $bwaDB = "$B37_BWA_INDEX";
}
elsif($species =~/^hg19$/i){

    die "hg19 is no longer supported in the variants pipeline";
    
    $bwaDB = "$HG19_BWA_INDEX";
}
elsif($species =~ /^mm10$|mouse/i){
    $bwaDB = "$MM10_BWA_INDEX";
}
elsif($species =~ /mm10_custom/i){
    $bwaDB = "$MM10_CUSTOM_BWA_INDEX";
}
elsif($species =~ /mm9/i){
    $bwaDB = "$MM9_BWA_INDEX";
}
elsif($species =~ /hybrid/i){
    $bwaDB = "$B37_MM10_HYBRID_BWA_INDEX";
}
elsif($species =~ /species_custom/i){
    $bwaDB = "$SPECIES_CUSTOM_BWA_INDEX";
}
elsif($species =~ /dm3/i){
    $bwaDB = "$DM3_BWA_INDEX";
}
else{
    die "SPECIES $species ISN'T CURRENTLY SUPPORTED; ONLY SUPPORT FOR b37, mm9|mm10|mm10_custom, and species_custom\n";
}

if(!-d "progress"){
    mkdir("progress", 0775) or die "Can't make progress";
}

open(IN, "$file") or die "Can't open $file file listing fastqs $!";
open(MISS, ">$file\_MISSING_READS.txt") or die "Can't write $file\_MISSING_READS.txt file listing fastqs $!";

my %addParams = (scheduler => "$scheduler", runtime => "50", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);

my $clipR1 = "-a $r1adaptor";
my $clipR2 = "-A $r2adaptor";
if($noClip){
    $clipR1 = "";
    $clipR2 = "";
}

my $bwaOptions = "-MP";
if($defaultBWA){
    $bwaOptions = "";
}

my $ran_bwa = 0;
my @bwa_jids = ();
while (<IN>){
    chomp;

    my @data = split(/\s+/, $_);
    my $pe = 0;
    ### if both ends don't exist, then skip the pair and output missing read name
    if(!-e $data[0]){
	if(!-e $data[0]){
	    print MISS "$data[0]\n";
	}
	next;
    }
    if(scalar(@data) == 2){
	if(!-e $data[1]){
	    if(!-e $data[1]){
		print MISS "$data[1]\n";
	    }
	    next;
	}
	$pe = 1;
    }

    my @nameR1 = split(/\.gz/, $data[0]);
    my @nameR2 = split(/\.gz/, $data[1]);

    $count++;

    if(!-d "$count"){
	mkdir("$count", 0775) or die "Can't make $count";
    }
    my $ran_zcatR1 = 0;
    my $zcatR1_jid = '';
    if(!-e "progress/$pre\_$uID\_ZCAT_$nameR1[0]\.done"){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT_$nameR1[0]", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_ZCAT_$nameR1[0]\.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "/bin/zcat $data[0] >$count/$nameR1[0]"`;
	`/bin/touch progress/$pre\_$uID\_ZCAT_$nameR1[0]\.done`;
	$zcatR1_jid = "$pre\_$uID\_ZCAT_$nameR1[0]";
	$ran_zcatR1 = 1;
    }

    my $ran_zcatR2 = 0;
    my $zcatR2_jid = '';
    if($pe){
	if(!-e "progress/$pre\_$uID\_ZCAT_$nameR2[0]\.done"){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT_$nameR2[0]", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_ZCAT_$nameR2[0]\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "/bin/zcat $data[1] >$count/$nameR2[0]"`;
	    `/bin/touch progress/$pre\_$uID\_ZCAT_$nameR2[0]\.done`;
	    $zcatR2_jid = "$pre\_$uID\_ZCAT_$nameR2[0]";
	    $ran_zcatR2 = 1;
	}
    }

    ###`$Bin/jobSync $scheduler $pre\_$uID\_ZCAT_$nameR1[0],$pre\_$uID\_ZCAT_$nameR2[0]`;
    if($ran_zcatR1){
	`$Bin/jobSync $scheduler $pre\_$uID\_ZCAT_$nameR1[0]`;
    }    
    if($ran_zcatR2){
    `$Bin/jobSync $scheduler $pre\_$uID\_ZCAT_$nameR2[0]`;
    }

    ### skipping if either fastq is empty
    if(!-s "$count/$nameR1[0]"){
	print "$count/$nameR1[0] is empty\n";
	next;
    }
    if($pe){
	if(!-s "$count/$nameR2[0]"){
	    print "$count/$nameR2[0] is empty\n";
	    next;
	}
    }

    ### discard reads less than half of original read length by cutadapt
    my $readO = `/usr/bin/head -2 $count/$nameR1[0]`;
    chomp $readO;
    my @dataO = split(/\n/, $readO);
    my $readLength = length($dataO[1]);
    my $minReadLength = int(0.5*$readLength);

    my $ran_cqs1 = 0;
    my $cqs1_jid = '';
    if(!-e "progress/$pre\_$uID\_CQS_$nameR1[0]\.done" || $ran_zcatR1){
	sleep(2);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CQS_$nameR1[0]", job_hold => "$zcatR1_jid", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_CQS_$nameR1[0]\.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/ConvertQualityScore --input $count/$nameR1[0] --output $count/$nameR1[0]\_CQS`;
	`/bin/touch progress/$pre\_$uID\_CQS_$nameR1[0]\.done`;
	$cqs1_jid = "$pre\_$uID\_CQS_$nameR1[0]";
	$ran_cqs1 = 1;	
    }

    my $ran_cqs2 = 0;
    my $cqs2_jid = '';
    if($pe){
	if(!-e "progress/$pre\_$uID\_CQS_$nameR2[0]\.done" || $ran_zcatR2){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CQS_$nameR2[0]", job_hold => "$zcatR2_jid", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_CQS_$nameR2[0]\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/ConvertQualityScore --input $count/$nameR2[0] --output $count/$nameR2[0]\_CQS`;
	    `/bin/touch progress/$pre\_$uID\_CQS_$nameR2[0]\.done`;
	    $cqs2_jid = "$pre\_$uID\_CQS_$nameR2[0]";
	    $ran_cqs2 = 1;	
	}
    }
    

    if($pe){
	my $ran_cutadapt = 0;
	my $ca_jid = '';
	if(!-e "progress/$pre\_$uID\_CUTADAPT_$nameR1[0]\_$nameR2[0]\_CQS.done" || $ran_cqs1 || $ran_cqs2){
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUTADAPT_$nameR1[0]\_$nameR2[0]\_CQS", job_hold => "$cqs1_jid,$cqs2_jid", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_CUTADAPT_$nameR1[0]\_$nameR2[0]\_CQS.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$PYTHON/python $CUTADAPT/cutadapt -f fastq $clipR1 $clipR2 --overlap 10 --minimum-length $minReadLength --quality-cutoff $bqTrim -o $count/$nameR1[0]\_CQS_CT_PE.fastq --paired-output $count/$nameR2[0]\_CQS_CT_PE.fastq $count/$nameR1[0]\_CQS $count/$nameR2[0]\_CQS >$count/$nameR1[0]\_$nameR2[0]\_CQS_CUTADAPT\_STATS.txt"`;
	    `/bin/touch progress/$pre\_$uID\_CUTADAPT_$nameR1[0]\_$nameR2[0]\_CQS.done`;
	    $ca_jid = "$pre\_$uID\_CUTADAPT_$nameR1[0]\_$nameR2[0]\_CQS";
	    $ran_cutadapt = 1;	
	}
	
	if(!-e "progress/$pre\_$uID\_BWA_$nameR1[0]\_$nameR2[0]\.done" || $ran_cutadapt){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_BWA_$nameR1[0]\_$nameR2[0]", job_hold => "$ca_jid", cpu => "12", mem => "12", cluster_out => "progress/$pre\_$uID\_BWA_$nameR1[0]\_$nameR2[0]\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$BWA/bwa mem $bwaOptions -R \'$readgroup\' -t 12 $bwaDB  $count/$nameR1[0]\_CQS_CT_PE.fastq $count/$nameR2[0]\_CQS_CT_PE.fastq >$count/$nameR1[0]\_$nameR2[0]\_CQS_CT_PE.fastq_$species\.bwa.sam"`;
	    `/bin/touch progress/$pre\_$uID\_BWA_$nameR1[0]\_$nameR2[0]\.done`;
	    push @bwa_jids, "$pre\_$uID\_BWA_$nameR1[0]\_$nameR2[0]";
	    $ran_bwa = 1;
	}
	push @raw, "I=$count/$nameR1[0]\_$nameR2[0]\_CQS_CT_PE.fastq_$species\.bwa.sam";
    }
    else{
	my $ran_cutadapt = 0;
	my $ca_jid = '';
	if(!-e "progress/$pre\_$uID\_CUTADAPT_$nameR1[0]\_CQS.done" || $ran_cqs1){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUTADAPT_$nameR1[0]\_CQS", job_hold => "$cqs1_jid", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_CUTADAPT_$nameR1[0]\_CQS.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$PYTHON/python $CUTADAPT/cutadapt -f fastq $clipR1 --overlap 10 --minimum-length $minReadLength --quality-cutoff $bqTrim -o $count/$nameR1[0]\_CQS_CT_SE.fastq $count/$nameR1[0]\_CQS >$count/$nameR1[0]\_CQS_CUTADAPT\_STATS.txt"`;
	    `/bin/touch progress/$pre\_$uID\_CUTADAPT_$nameR1[0]\_CQS.done`;
	    $ca_jid = "$pre\_$uID\_CUTADAPT_$nameR1[0]\_CQS";
	    $ran_cutadapt = 1;
	}
	
	if(!-e "progress/$pre\_$uID\_BWA_$nameR1[0]\.done" || $ran_cutadapt){
	    sleep(2);
	    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_BWA_$nameR1[0]", job_hold => "$ca_jid", cpu => "12", mem => "12", cluster_out => "progress/$pre\_$uID\_BWA_$nameR1[0]\.log");
	    my $standardParams = Schedule::queuing(%stdParams);
	    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "$BWA/bwa mem $bwaOptions -R \'$readgroup\' -t 12 $bwaDB $count/$nameR1[0]\_CQS_CT_SE.fastq >$count/$nameR1[0]\_CQS_CT_SE.fastq_$species\.bwa.sam"`;
	    `/bin/touch progress/$pre\_$uID\_BWA_$nameR1[0]\.done`;
	    push @bwa_jids, "$pre\_$uID\_BWA_$nameR1[0]";
	    $ran_bwa = 1;
	}
	push @raw, "I=$count/$nameR1[0]\_CQS_CT_SE.fastq_$species\.bwa.sam";
    }
}
close IN;
close MISS;

my $inputRaw = join(" ", @raw);
### NOTE: MERGE DOESN'T TAKE UP MUCH MEMORY
###       BUT I/O IS SUCH THAT IT DOESN'T SEEM TO LIKE
###       HAVING OTHER JOBS RUNNING AT THE SAME TIME;
###       OTHERWISE IT RISKS HANGING OR DYING WITH NO ERROR MESSAGE

my %addParams = (scheduler => "$scheduler", runtime => "300", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "6");
my $additionalParams = Schedule::additionalParams(%addParams);
my $bwa_holds = join(",", @bwa_jids);
if(!-e "progress/$pre\_$uID\_$run\_MERGE.done" || $ran_bwa){
    sleep(2);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$run\_MERGE", job_hold => "$bwa_holds", cpu => "8", mem => "30", cluster_out => "progress/$pre\_$uID\_$run\_MERGE.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tempdir -jar $PICARD/picard.jar MergeSamFiles $inputRaw O=$run\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tempdir CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
    `/bin/touch progress/$pre\_$uID\_$run\_MERGE.done`;
}
