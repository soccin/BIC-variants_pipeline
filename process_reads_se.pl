#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

#### new casava splits fastqs into batches of 4m reads
### this will process each fastq
### and create a single bam
### that can go trhough the second half of pipeline

### takes in a file that list of gzipped read pairs
## e.g.LID46438_NoIndex_L001_R1_001.fastq.gz
##     LID46438_NoIndex_L001_R1_002.fastq.gz
##     LID46438_NoIndex_L001_R1_003.fastq.gz
##     LID46438_NoIndex_L001_R1_004.fastq.gz

### '@RG\tID:TEST_SFT20_SE\tPL:Illumina\tPU:TEST_SFT20\tLB:TEST_SFT20\tSM:20'

my ($file, $species, $config, $scheduler, $help, $priority_project, $priority_group);
my $pre = 'TEMP';
my $run = 'TEMP_RUN';
my $readgroup = '@RG\tID:TEMP_ID\tPL:Illumina\tPU:TEMP_PU\tLB:TEMP_LBSM:TEMP_SM';
GetOptions ('file=s' => \$file,
	    'pre=s' => \$pre,
	    'species=s' => \$species,
	    'run=s' => \$run,
	    'config=s' => \$config,
 	    'scheduler=s' => \$scheduler,
	    'readgroup=s' => \$readgroup) or exit(1);

if(!$file || !$config || !$species || !$scheduler || $help){
    print <<HELP;

    USAGE: process_reads_se.pl -file FILE -pre PRE -species SPECIES -config CONFIG -run RUN -readgroup READGROUP -scheduler SCHEDULER
	* FILE: file listing read pairs to process (REQUIRED)
	* PRE: output prefix (default: TEMP)
	* SPECIES: only hg19, mm9, mm10, hybrid (hg19+mm10) and dm3 currently supported (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
	* RUN: RUN IDENTIFIER (default: TEMP_RUN)
	* READGROUP: string containing information for ID, PL, PU, LB, and SM e.g. '\@RG\\tID:s_CH_20__1_LOLA_1018_PE\\tPL:Illumina\\tPU:s_CH_20__1_LOLA_1018\\tLB:s_CH_20__1\\t SM:s_CH_20'; WARNING: INCOMPLETE INFORMATION WILL CAUSE GATK ISSSUES (default: '\@RG\\tID:TEMP_ID\\tPL:Illumina\\tPU:TEMP_PU\\tLB:TEMP_LB\\tSM:TEMP_SM')
HELP
exit;
}

my @raw = ();
my $count = 0;
$readgroup =~ s/\s+/\\t/g;

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my $bwaDB = '';
my $knownSites = '';
my $refSeq = '';
my $HG19_FASTA = '';
my $HG19_BWA_INDEX = '';
my $HG19_MM9_HYBRID_FASTA = '';
my $HG19_MM9_HYBRID_BWA_INDEX = '';
my $MM9_FASTA = '';
my $MM9_BWA_INDEX = '';
my $MM10_FASTA = '';
my $MM10_BWA_INDEX = '';
my $DM3_FASTA = '';
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
	$HG19_BWA_INDEX = $conf[1];
    }
    elsif($conf[0] =~ /hg9_mm9_hybrid_fasta/i){
	if(!-e "$conf[1]"){
	    die "CAN'T FIND $conf[1] $!";
	}
	$HG19_MM9_HYBRID_FASTA = $conf[1];
    }
    elsif($conf[0] =~ /hg19_mm9_hybrid_bwa_index/i){
	if(!-e "$conf[1]\.bwt" || !-e "$conf[1]\.pac" || !-e "$conf[1]\.ann" || !-e "$conf[1]\.amb" || !-e "$conf[1]\.sa"){
	    die "CAN'T FIND ALL NECESSARY BWA INDEX FILES FOR HG19-MM9 HYBRID WITH PREFIX $conf[1] $!";
	}
	$HG19_MM9_HYBRID_BWA_INDEX = $conf[1];
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
	$MM9_BWA_INDEX = $conf[1];
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
	$MM10_BWA_INDEX = $conf[1];
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
 	$DM3_BWA_INDEX = $conf[1];
    }
}
close CONFIG;

if($species =~/hg19|human/i){
    $bwaDB = "$HG19_BWA_INDEX";
    $refSeq = "$HG19_FASTA";
}
elsif($species =~ /mm10|mouse/i){
    $bwaDB = "$MM10_BWA_INDEX";
    $refSeq = "$MM10_FASTA";
}
elsif($species =~ /mm9/i){
    $bwaDB = "$MM9_BWA_INDEX";
    $refSeq = "$MM9_FASTA";
}
elsif($species =~ /hybrid/i){
    $bwaDB = "$HG19_MM9_HYBRID_BWA_INDEX";
    $refSeq = "$HG19_MM9_HYBRID_FASTA";
}
elsif($species =~ /dm3/i){
    $bwaDB = "$DM3_BWA_INDEX";
    $refSeq = "$DM3_FASTA";
}
else{
    die "SPECIES $species ISN'T CURRENTLY SUPPORTED; ONLY SUPPORT FOR hg19 and mm9|mm10\n";
}

`/bin/mkdir -m 775 -p progress`;
sleep(3);
my $ran_bwa = 0;
open(IN, "$file") or die "Can't open $file file listing fastqs $!";

my %addParams = (scheduler => "$scheduler", runtime => "50", priority_project=> "$priority_project", priority_group=> "$priority_group", rerun => "1", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);

my @bwa_jids = ();
while (<IN>){
    chomp;

    my @data = split(/\s+/, $_);
    my @nameR1 = split(/\.gz/, $data[0]);
    $count++;
    
    `/bin/mkdir -m 775 -p $count`;
    sleep(3);
        
    my $ran_zcatR1 = 0;
    my $zcatR1_jid = '';
    if(!-e "progress/$pre\_$uID\_ZCAT_$nameR1[0]\.done"){
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ZCAT_$nameR1[0]", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_ZCAT_$nameR1[0]\.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $additionalParams "/bin/zcat $data[0] >$count/$nameR1[0]"`;
	`/bin/touch progress/$pre\_$uID\_ZCAT_$nameR1[0]\.done`;
	$zcatR1_jid = "$pre\_$uID\_ZCAT_$nameR1[0]";
	$ran_zcatR1 = 1;
    }

    my $ran_rsplit = 0;
    my $rsplit_jid = '';
    if(!-e "progress/$pre\_$uID\_RSPLIT_$nameR1[0]\.done" || $ran_zcatR1){
	sleep(3);
	my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RSPLIT_$nameR1[0]", job_hold => "$zcatR1_jid", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_RSPLIT_$nameR1[0]\.log");
	my $standardParams = Schedule::queuing(%stdParams);
	`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $additionalParams /usr/bin/split -a 3 -l 16000000 -d $count/$nameR1[0] $count/$nameR1[0]\__`;
	`/bin/touch progress/$pre\_$uID\_RSPLIT_$nameR1[0]\.done`;
	$rsplit_jid = "$pre\_$uID\_RSPLIT_$nameR1[0]";
	$ran_rsplit = 1;	
    }

    if($ran_rsplit){
	`$Bin/jobSync $scheduler $pre\_$uID\_RSPLIT_$nameR1[0]`;
    }

    ### skipping empty fastq files
    if(!-s "$count/$nameR1[0]"){
	print "$count/$nameR1[0] is empty\n";
	next;
    }

    ### discard reads less than half of original read length by cutadapt
    my $readO = `/usr/bin/head -2 $count/$nameR1[0]`;
    chomp $readO;
    my @dataO = split(/\n/, $readO);
    my $readLength = length($dataO[1]);
    my $minReadLength = int(0.5*$readLength);

    opendir(RS, "$count") or die "Can't open read split directory $count";
    my @rsplits = readdir RS;
    closedir RS;

    my $rcount = 0;
    foreach my $rfile (@rsplits){
	if($rfile =~ /$nameR1[0]\_\_/){
	    $rcount++;

	    `/bin/mkdir -m 775 -p $count/$rcount`;
	    sleep(3);
	    my $ran_cutadapt = 0;
	    my $ca_jid = '';
	    if(!-e "progress/$pre\_$uID\_CUTADAPT_$rfile\.done" || $ran_rsplit){
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CUTADAPT_$rfile", job_hold => "$rsplit_jid", cpu => "1", mem => "1", cluster_out => "progress/$pre\_$uID\_CUTADAPT_$rfile\.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $additionalParams $PYTHON/python "$CUTADAPT/cutadapt -f fastq -a AGATCGGAAGAGCACACGTCT -O 10 -m $minReadLength -o $count/$rcount/$rfile\_CT.fastq -q 3 $count/$rfile >$count/$rcount/$rfile\_CUTADAPT\_STATS.txt"`;
		`/bin/touch progress/$pre\_$uID\_CUTADAPT_$rfile\.done`;
		$ca_jid = "$pre\_$uID\_CUTADAPT_$rfile";
		$ran_cutadapt = 1;
	    }
	    	    
	    if(!-e "progress/$pre\_$uID\_BWA_$rfile\.done" || $ran_cutadapt){
		sleep(3);
		my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_BWA_$rfile", job_hold => "$ca_jid", cpu => "12", mem => "12", cluster_out => "progress/$pre\_$uID\_BWA_$rfile\.log");
		my $standardParams = Schedule::queuing(%stdParams);
		`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $additionalParams "$BWA/bwa mem -PM -R \'$readgroup\' -t 12 $bwaDB $count/$rcount/$rfile\_CT.fastq ">$count/$rcount/$rfile\_CT.fastq_$species\.bwa.sam"`;
		`/bin/touch progress/$pre\_$uID\_BWA_$rfile\.done`;
		push @bwa_jids, "$pre\_$uID\_BWA_$rfile";
		$ran_bwa = 1;
	    }
	    	    
	    push @raw, "I=$count/$rcount/$rfile\_CT.fastq_$species\.bwa.sam";
	}
    }
}
close IN;

my $inputRaw = join(" ", @raw);
my %addParams = (scheduler => "$scheduler", runtime => "300", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "1", iounits => "6");
my $additionalParams = Schedule::additionalParams(%addParams);
my $bwa_holds = join(",", @bwa_jids);
if(!-e "progress/$pre\_$uID\_$run\_MERGE.done" || $ran_bwa){
    sleep(3);
    my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_$run\_MERGE", job_hold => "$bwa_holds", cpu => "8", mem => "30", cluster_out => "progress/$pre\_$uID\_$run\_MERGE.log");
    my $standardParams = Schedule::queuing(%stdParams);
    `$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $additionalParams $JAVA/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $inputRaw O=$run\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
    `/bin/touch progress/$pre\_$uID\_$run\_MERGE.done`;
}
