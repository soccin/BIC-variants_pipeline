#!/opt/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

### takes in a project dir (e.g. /ifs/data/mpirun/analysis/solid/offit/20120117)
### assumes project directory structure is as such - sample/lib(s)/run(s)
### will merge all lib/run/*_resolved.bam
###  -keepdups: marks duplicate reads, but does not remove them



my ($pre, $projectDir, $config, $keepdups);

GetOptions ('pre=s' => \$pre,
	    'config=s' => \$config,
	    'keepdups' => \$keepdups,
            'dir=s' => \$projectDir);

### DEFAULT PICARD LOCATION
my $PICARD = '/opt/bin/picard';
open(CONFIG, "$config") or die "Can't open config file $config";
while(<CONFIG>){
    chomp;

    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /picard/i){
	$PICARD = $conf[1];
    }
}
close CONFIG;

my $removedups = "true";
if($keepdups){
    $removedups = "false";
}


opendir(pDir, "$projectDir") or die "Can't open $projectDir directory";
my @pDirs = readdir pDir;
closedir pDir;

my $workDir = `pwd`;
chomp $workDir;

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

foreach my $pd (@pDirs){
    chomp $pd;

    if($pd eq "." || $pd eq ".." || !-d $pd){
	next;
    }

    opendir(lDir, "$pd") or die "Can't open $pd directory";
    my @lDirs = readdir lDir;
    closedir lDir;

    foreach my $ld (@lDirs){
	if($ld eq "." || $ld eq ".." || !-d "$pd/$ld"){
	    next;
	}

	opendir(rDir, "$pd/$ld") or die "Can't open $pd/$ld directory";
	my @rDirs = readdir rDir;
	closedir rDir;
	my $dirFound = 0;

	my @mergeBams = ();
	foreach my $rd (@rDirs){
	    if($rd eq "." || $rd eq ".." || !-d "$pd/$ld/$rd"){
		next;
	    }

	    $dirFound++;
	    opendir(hDir, "$pd/$ld/$rd") or die "Can't open $pd/$ld/$rd directory";
	    my @hDirs = readdir hDir;
	    closedir hDir;

	    foreach my $hd (@hDirs){
		if($hd eq "." || $rd eq ".."){
		    next;
		}

		if($hd =~ /_FILTERED\.bam$/){
		    push @mergeBams, "I=$pd/$ld/$rd/$hd";
		}
	    }
	}

	if($dirFound !~ scalar(@mergeBams)){
	    print "SKIPPING $pd/$ld BECAUSE IT'S MISSING FILTERED BAM FILES\n";
	    next;
	}
	

	if(scalar(@mergeBams) == 1){

	    my @name = split(/=/, $mergeBams[0]);
	    my @path = split(/\//, $name[1]);

	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MARKDUPS -pe alloc 3 -l virtual_free=10G $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MarkDuplicates.jar INPUT=$name[1] OUTPUT=$pd/$ld/$pd\_$ld\_MD.bam METRICS_FILE=$pd/$ld/$pd\_$ld\_removedDupsMetrics TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$removedups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;
	}
	elsif(scalar(@mergeBams) > 1){
	    my $merin = join(" ", @mergeBams);

	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MERGE_$pd\_$ld -pe alloc 24 -l virtual_free=3G -q lau.q $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MergeSamFiles.jar $merin O=$pd/$ld/$pd\_$ld\.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/$uID CREATE_INDEX=true USE_THREADING=true MAX_RECORDS_IN_RAM=5000000`;

	    sleep(5);

	    `/common/sge/bin/lx24-amd64/qsub -P ngs -N $pre\_$uID\_MARKDUPS -hold_jid $pre\_$uID\_MERGE_$pd\_$ld -pe alloc 3 -l virtual_free=10G $Bin/qCMD /opt/bin/java -Djava.io.tmpdir=/scratch/$uID -jar $PICARD/MarkDuplicates.jar INPUT=$pd/$ld/$pd\_$ld\.bam OUTPUT=$pd/$ld/$pd\_$ld\_MD.bam METRICS_FILE=$pd/$ld/$pd\_$ld\_removedDupsMetrics TMP_DIR=/scratch/$uID VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$removedups CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000`;

	}
    }
}
