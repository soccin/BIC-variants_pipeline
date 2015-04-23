#!/opt/bin/perl

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

my ($file, $readgroup, $pre, $species, $run, $config, $help);
$pre = 'TEMP';
$run = 'TEMP_RUN';
$readgroup = '@RG\tID:TEMP_ID\tPL:Illumina\tPU:TEMP_PU\tLB:TEMP_LBSM:TEMP_SM';
GetOptions ('file=s' => \$file,
	    'pre=s' => \$pre,
	    'species=s' => \$species,
	    'run=s' => \$run,
	    'config=s' => \$config,
	    'readgroup=s' => \$readgroup);

if(!$file || !$config || !$species || $help){
    print <<HELP;

    USAGE: ./exome_pipeline/solexa_PE.pl -file FILE -pre PRE -species SPECIES -config CONFIG -run RUN -readgroup READGROUP
	* FILE: file listing read pairs to process (REQUIRED)
	* PRE: output prefix (default: TEMP)
	* SPECIES: only hg19, mm9, hybrid (hg19+mm10) and dm3 currently supported (REQUIRED)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
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

if($species =~/hg19/){
    $bwaDB = '/ifs/data/bio/assemblies/H.sapiens/hg19/bwa/hg19.fasta';
    $refSeq = '/ifs/data/bio/assemblies/H.sapiens/hg19/hg19.fasta';
}
elsif($species =~ /mm9/){
    $bwaDB = '/ifs/data/bio/assemblies/M.musculus/mm9/bwa/mm9.fasta';
    $refSeq = '/ifs/data/bio/assemblies/M.musculus/mm9/mm9.fasta';
}
elsif($species =~ /hybrid/i){
    $bwaDB = '/ifs/data/bio/assemblies/hybrid_H.sapiens_M.musculus/hybrid_hg19_mm10/bwa/hybrid_hg19_mm10.fasta';
    $refSeq = '/ifs/data/bio/assemblies/hybrid_H.sapiens_M.musculus/hybrid_hg19_mm10/hybrid_hg19_mm10.fasta';
}
elsif($species =~ /dm3/){
    $bwaDB = '/ifs/data/bio/assemblies/D.melanogaster/dm3/bwa/dm3.fasta';
    $refSeq = '/ifs/data/bio/assemblies/D.melanogaster/dm3/dm3.fasta';
}
else{
    die "SPECIES $species ISN'T CURRENTLY SUPPORTED; ONLY SUPPORT FOR hg19 and mm9\n";
}

my $CUTADAPT = '/opt/bin/cutadapt';
my $BWA = '/opt/bin/bwa';
my $PICARD = '/opt/bin/picard';
open(CONFIG, "$config") or warn "CAN'T OPEN CONFIG FILE $config SO USING DEFAULT SETTINGS $!";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);

    if($conf[0] =~ /cutadapt/i){
	$CUTADAPT = $conf[1];
    }
    elsif($conf[0] =~ /bwa/i){
	$BWA = $conf[1];
    }
    elsif($conf[0] =~ /picard/i){
	$PICARD = $conf[1];
    }
}
close CONFIG;

`/bin/mkdir -m 775 -p progress`;
my $ran_bwa = 0;
open(IN, "$file") or die "Can't open $file file listing fastqs $!";
while (<IN>){
    chomp;

    my @data = split(/\s+/, $_);
    my @nameR1 = split(/\.gz/, $data[0]);
    $count++;
    
    `/bin/mkdir -m 775 -p $count`;
        
    my $ran_zcatR1 = 0;
    if(!-e "progress/$pre\_$uID\_ZCAT_$nameR1[0]\.done"){
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_ZCAT_$nameR1[0] -pe alloc 1 -l virtual_free=1G $Bin/qCMD /bin/zcat $data[0] ">$count/$nameR1[0]"`;
	`/bin/touch progress/$pre\_$uID\_ZCAT_$nameR1[0]\.done`;
	$ran_zcatR1 = 1;
    }

    my $ran_rsplit = 0;
    if(!-e "progress/$pre\_$uID\_RSPLIT_$nameR1[0]\.done" || $ran_zcatR1){
	sleep(3);
	`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_RSPLIT_$nameR1[0] -hold_jid $pre\_$uID\_ZCAT_$nameR1[0] -pe alloc 1 -l virtual_free=1G $Bin/qCMD /usr/bin/split -a 3 -l 16000000 -d $count/$nameR1[0] $count/$nameR1[0]\__`;
	$ran_rsplit = 1;	
    }

    `$Bin/qSYNC $pre\_$uID\_RSPLIT_$nameR1[0]`;

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
	    
	    my $ran_cutadapt = 0;
	    if(!-e "progress/$pre\_$uID\_CUTADAPT_$rfile\.done" || $ran_rsplit){
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_CUTADAPT_$rfile -hold_jid $pre\_$uID\_RSPLIT_$nameR1[0] -pe alloc 1 -l virtual_free=1G $Bin/qCMD $CUTADAPT/bin/cutadapt -a AGATCGGAAGAGCACACGTCT -O 10 -m $minReadLength -o $count/$rcount/$rfile\_CT.fastq -q 3 $count/$rfile > $count/$rfile\_CUTADAPT\_STATS.txt"`;
		`/bin/touch progress/$pre\_$uID\_CUTADAPT_$rfile\.done`;
		$ran_cutadapt = 1;
	    }
	    	    
	    if(!-e "progress/$pre\_$uID\_BWA_$rfile\.done" || $ran_cutadapt){
		sleep(3);
		`/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_BWA_$run -hold_jid $pre\_$uID\_CUTADAPT_$rfile -pe alloc 12 -l virtual_free=1G $Bin/qCMD $BWA/bwa mem -PM -R "'$readgroup'" -t 12 $bwaDB $count/$rcount/$rfile\_CT.fastq ">$count/$rcount/$rfile\_CT.fastq_$species\.bwa.sam"`;
		`/bin/touch progress/$pre\_$uID\_BWA_$rfile\.done`;
		$ran_bwa = 1;
	    }
	    	    
	    push @raw, "I=$count/$rcount/$rfile\_CT.fastq_$species\.bwa.sam";
	}
    }
}
close IN;

my $inputRaw = join(" ", @raw);
if(!-e "progress/$pre\_$uID\_RUN_MERGE.done" || $ran_bwa){
    sleep(3);
    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_RUN_MERGE -hold_jid $pre\_$uID\_BWA_$run -pe alloc 8 -l virtual_free=4G -q lau.q,lcg.q,nce.q $Bin/qCMD /opt/bin/java -Xms256m -Xmx30g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/picard.jar MergeSamFiles $inputRaw O=$run\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=false MAX_RECORDS_IN_RAM=5000000`;
    `/bin/touch progress/$pre\_$uID\_RUN_MERGE.done`;
}
