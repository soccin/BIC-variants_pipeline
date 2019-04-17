#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
use Schedule;
use Cluster;


my ($map, $group, $pair, $config, $help, $species, $scheduler, $targets, $request, $strand, $rnaseq_pipeline);

my $pre = 'TEMP';
my $output = "results";
my $priority_project = "ngs";
my $priority_group = "Pipeline";

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

my $email = "$uID\@cbio.mskcc.org";
my $rsync = "/ifs/res/$uID";
my $tempdir = "/scratch/$uID";


GetOptions ('map=s' => \$map,
            'group=s' => \$group,
            'pair=s' => \$pair,
            'pre=s' => \$pre,
            'strand=s' => \$strand,
            'targets=s' => \$targets,
            'request=s' => \$request,
            'species=s' => \$species,
            'scheduler=s' => \$scheduler,
            'help' => \$help,
            'output|out|o=s' => \$output,
            'rsync=s' => \$rsync,
            'priority_project=s' => \$priority_project,
            'priority_group=s' => \$priority_group,
            'email' => \$email,
            'tempdir=s' => \$tempdir,
            'rnaseq_pipeline=s' => \$rnaseq_pipeline) or exit(1);

if(!$map || !$group || !$pair || !$species  || !$scheduler || !$request || !$targets || !$strand || !$rnaseq_pipeline || $help){
    print <<HELP;

    USAGE: variants_pipeline.pl -wes -config CONFIG -species SPECIES -scheduler SCHEDULER
        * MAP: file listing sample information for processing (REQUIRED)
        * GROUP: file listing grouping of samples for realign/recal steps (REQUIRED, unless using -mdOnly flag)
        * SPECIES: b37 (default: human), mm9, mm10 (default: mouse), hybrid (b37+mm10), mm10_custom, species_custom and dm3 currently supported (REQUIRED)
        * TARGETS: name of targets assay; will search for targets/baits ilists and targets padded file in $Bin/targets/TARGETS unless given full path to targets directory; required for non-chipseq projects
        * REQUEST: file containing request information such as PI, investigator, etc. (REQUIRED)
        * SCHEDULER: currently support for SGE, LUNA, and JUNO (REQUIRED)
        * EMAIL: email to send notication of finished final job of pipeline (default: $uID\@cbio.mskcc.org)
        * PAIR: file listing tumor/normal pairing of samples for mutect/maf conversion; if not specified, considered unpaired
        * PRE: output prefix (default: TEMP)
        * OUTPUT: output results directory (default: results)
        * RSYNC:  path to rsync data for archive (default: /ifs/res/$uID)
        * TEMPDIR:  temp directory (default: /scratch/$uID)
        * PRIORITY_PROJECT: sge notion of priority assigned to projects (default: ngs)
        * PRIORITY_GROUP: lsf notion of priority assigned to groups (default: Pipeline)
        * STANDARD_GENE: standard analysis - star alignment, htseq gene count
        * RNASEQ_PIPELINE: path to the rnaseq pipeline repository
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

my $rna_output = "$output/rna";


die "Can't find mapping file $map\n" if(!-e $map);
die "Can't find pairing file $pair\n" if(!-e $pair);
die "Can't find grouping file $group\n" if(!-e $group);
die "Can't find request file $request\n" if(!-e $request);
die "Can't find rsync directory $rsync\n" if(!-d $rsync);
die "Can't find rnaseq pipeline script $rnaseq_pipeline/rnaseq_pipeline.pl\n" if(!-e "$rnaseq_pipeline/rnaseq_pipeline.pl");

if(!-d $output){
    mkdir("$output", 0775) or die "Can't make $output";
    mkdir("$output/progress", 0775) or die "Can't make $output/progress";
}


if(!-d $rna_output){
    mkdir("$rna_output", 0775) or die "Can't make $rna_output";
}

if(!-d $tempdir){
    mkdir("$tempdir", 0775) or die "Can't make $tempdir";
}


my %addParams = (scheduler => "$scheduler", runtime => "500", priority_project=> "$priority_project", priority_group=> "$priority_group", queues => "lau.q,lcg.q,nce.q", rerun => "0", iounits => "1");
my $additionalParams = Schedule::additionalParams(%addParams);





my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RNA_VARIANTS_PREPROCESSING", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_RNA_VARIANTS_PREPROCESSING.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $rnaseq_pipeline/rnaseq_pipeline.pl -species $species -scheduler lsf -priority_project $priority_project -priority_group $priority_group -email $email -request $request -pre $pre -map $map -config $rnaseq_pipeline/rnaseq_pipeline_config.txt -strand $strand -star -o $rna_output -rsync null -alignment_only`;



`$Bin/jobSync $scheduler $pre\_$uID\_RNA_VARIANTS_PREPROCESSING`;



my %stdParams = (scheduler => "$scheduler", job_name => "$pre\_$uID\_RNA_VARIANTS_CALLING", job_hold => "$pre\_$uID\_RSYNC", cpu => "1", mem => "10", cluster_out => "$output/progress/$pre\_$uID\_RNA_VARIANTS_CALLING.log");
my $standardParams = Schedule::queuing(%stdParams);
`$standardParams->{submit} $standardParams->{job_name} $standardParams->{job_hold} $standardParams->{cpu} $standardParams->{mem} $standardParams->{cluster_out} $additionalParams $Bin/variants_pipeline.pl -species $species -scheduler $scheduler -priority_project $priority_project -priority_group $priority_group -tempdir $tempdir -email $email -request $request -pre $pre -map $map -config $Bin/variants_pipeline_config.txt -group $group -pair $pair -targets $targets -wes -rna $rna_output/gene/alignments/ -indelrealigner -o $output -rsync $rsync`;

























