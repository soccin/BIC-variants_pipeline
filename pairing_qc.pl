#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Schedule;
use Cluster;

my ($pair, $vcf, $help, $outdir, $scheduler, $config, $scheduler, $priority_project, $priority_group, $holdjid, $progress);
my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my $pre = 'TEMP';
GetOptions ('pair=s' => \$pair,
            'vcf=s' => \$vcf,
            'pre=s' => \$pre,
            'outdir=s' => \$outdir,
 	    'scheduler=s' => \$scheduler,
 	    'priority_project=s' => \$priority_project,
 	    'priority_group=s' => \$priority_group,
	    'config=s' => \$config,
	    'holdjid=s' => \$holdjid,
	    'progress=s' => \$progress,
            'help' => \$help) or exit(1);

if(!$pair || !$config || !$scheduler || $help){
    print <<HELP;

    USAGE: ./pairing_qc.pl -vcf VCF -pair PAIR -pre PRE 
        * PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
        * VCF: population VCF file from which SNPs will be selected for pair evaluation
        * PRE: output prefix (default: TEMP)
        * BIN: path to root of exome pipeline
        * OUTDIR: output directory for all fingerprint files
	* PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
	* PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
	* CONFIG: file listing paths to programs needed for pipeline; full path to config file needed (REQUIRED)
	* SCHEDULER: currently support for SGE and LSF (REQUIRED)
HELP
exit;
}

my $alerts_dir = "/ifs/data/alerts/$pre";
my $PYTHON = '';
my $R = '';

open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;
    
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /^r$/i){
	if(!-e "$conf[1]/R"){
	    die "CAN'T FIND R IN $conf[1] $!";
	}
	$R = $conf[1];
    }
    elsif($conf[0] =~ /python/i){
	if(!-e "$conf[1]/python"){
	    die "CAN'T FIND python IN $conf[1] $!";
	}
	$PYTHON = $conf[1];
    }
}
close CONFIG;

## move to output dir
`cd $outdir`;
if(!$progress){
    `mkdir -m 775 $outdir/progress`;
    $progress = "$outdir/progress";
}

## make a copy of pairing file in out directory for web app
`cp $pair $outdir`;

my %addParams = (scheduler => "$scheduler", runtime => "50", priority_project=> "$priority_project", priority_group=> "$priority_group");
my $additionalParams = Schedule::additionalParams(%addParams);

## filter VCF
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_FP", job_hold => "$holdjid", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_FP.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/fingerprint.py $vcf $outdir/$pre\_FINGERPRINT.vcf`;

## flatten & further filter VCF
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_FFP", job_hold => "$pre\_$uID\_FP", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_FFP.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams "cat $outdir/$pre\_FINGERPRINT.vcf | $PYTHON/python $Bin/qc/flattenFingerprint.py > $outdir/$pre\_FINGERPRINT_FLAT.txt"`;

## plot flattened fingerprint
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_PLOT_FFP", job_hold => "$pre\_$uID\_FFP", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_PLOT_FFP.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $R/Rscript $Bin/qc/plotFingerprint.R $outdir/$pre\_FINGERPRINT_FLAT.txt`;

## pairing by pval
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_FINDPAIRS_PVAL", job_hold => "$pre\_$uID\_FP", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_FINDPAIRS_PVAL.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/pval_findPairs.py -vcf $outdir/$pre\_FINGERPRINT.vcf -out $outdir/$pre\_FINGERPRINT.vcf_pval_pairs.txt`;

## clustering by pval
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CLUSTER_PVAL", job_hold => "$pre\_$uID\_FINDPAIRS_PVAL", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_CLUSTER_PVAL.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/pval_cluster.py $outdir/$pre\_FINGERPRINT.vcf_pval_pairs.txt $outdir/$pre\_FINGERPRINT.vcf_pval_pairs_cluster.txt 0.05`;

## pairing by IBS
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_FINDPAIRS_IBS", job_hold => "$pre\_$uID\_FP", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_FINDPAIRS_IBS.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/ibs_findPairs.py $outdir/$pre\_FINGERPRINT.vcf $outdir/$pre\_FINGERPRINT.vcf_ibs_pairs.txt`;

## clustering by IBS
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_CLUSTER_IBS", job_hold => "$pre\_$uID\_FINDPAIRS_IBS", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_CLUSTER_IBS.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $PYTHON/python $Bin/qc/ibs_cluster.py $outdir/$pre\_FINGERPRINT.vcf_ibs_pairs.txt $outdir/$pre\_FINGERPRINT.vcf_ibs_pairs_cluster.txt 1.8 0.3`;

## IBS plot
my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_PLOT_IBS", job_hold => "$pre\_$uID\_FINDPAIRS_IBS", cpu => "1", mem => "2", cluster_out => "$progress/$pre\_$uID\_PLOT_IBS.log");
my $standardParams = Schedule::queuing(%stdParams);
###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $R/Rscript $Bin/qc/ibs_plot.R \"\\\"projID='$pre'\\\"\" \"\\\"manifestPairing='$pair'\\\"\" \"\\\"IBSFile='$outdir/$pre\_FINGERPRINT.vcf_ibs_pairs.txt'\\\"\"`;
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/qc/ibs_plot.pl -pre $pre -pair $pair -config $config -ibs_file $outdir/$pre\_FINGERPRINT.vcf_ibs_pairs.txt`;

## copy txt files to alerts dir to be displayed in app
###my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_ALERT", job_hold => "$pre\_$uID\_CLUSTER_IBS,$pre\_$uID\_CLUSTER_PVAL,$pre\_$uID\_FFP", cpu => "1", mem => "2");
###my $standardParams = Schedule::queuing(%stdParams);
###`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $additionalParams $Bin/qc/alert.sh $pair $outdir $alerts_dir`;
