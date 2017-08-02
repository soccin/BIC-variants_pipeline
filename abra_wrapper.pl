#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use lib "$Bin/lib";

my ($inBams, $outBams, $refSeq, $bwaRef, $targets, $working, $config, $help, $log);
my $uID = `/usr/bin/id -u -n`;
chomp $uID;

GetOptions ('inBams=s' => \$inBams,
            'outBams=s' => \$outBams,
 	    'refSeq=s' => \$refSeq,
            'bwaRef=s' => \$bwaRef,
 	    'targets=s' => \$targets,
            'working=s' => \$working,
	    'config=s' => \$config,
	    'log=s' => \$log,
            'help' => \$help) or exit(1);



my $ABRA = '';
my $JAVA = '';
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;
    
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /abra/i){
	if(!-e "$conf[1]/abra.jar"){
	    die "CAN'T FIND abra IN $conf[1] $!";
	}
	$ABRA = $conf[1];
    }
    #elsif($conf[0] =~ /^bwa/i){
    #	if(!-e "$conf[1]/bwa"){
    #	    die "CAN'T FIND bwa IN $conf[1] $!";
    #	}
    #	my $path_tmp = $ENV{'PATH'};
    #	$ENV{'PATH'} = "$conf[1]:$path_tmp";
    #}
    elsif($conf[0] =~ /^java$/i){
	if(!-e "$conf[1]/java"){
	    die "CAN'T FIND java IN $conf[1] $!";
	}
	$JAVA = $conf[1];
    }
}
close CONFIG;

#if(-d "$working"){
    ### abra DIES IF DIR ALREADY EXISTS
#    `/bin/rm -rf $working`;
#}

#abra2 DIES IF DIR DOES NOT EXIST
`mkdir -p $working`;
    


###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $JAVA/java -Xms256m -Xmx60g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $ABRA/abra.jar --in $inBams --out $outBams --ref $refSeq --bwa-ref $bwaRef --targets $targets --working $working --threads 24`;

`$JAVA/java -Xms256m -Xmx60g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=/scratch/$uID -jar $ABRA/abra.jar --in $inBams --out $outBams --ref $refSeq --targets $targets --tmpdir $working --threads 12 --mnf 5 --mbq 150 > $log 2>&1`;

my $eadj = $? >> 8;
exit($eadj);
