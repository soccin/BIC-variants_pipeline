#!/opt/bin/perl

package Schedule;

use strict;
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Cluster;


my $cluster = new Cluster;

sub queuing {
    my %qSystem = @_;

    my $addParams = '';
    if($qSystem{'scheduler'} =~ /SGE/i){
	sge(%qSystem);
	#$addParams = additionalSGE(%qSystem);
    }
    elsif($qSystem{'scheduler'} =~ /lsf/i){
	lsf(%qSystem);
	#$addParams = additionalLSF(%qSystem);
    }
    else{
	die "SCHEDULER $qSystem{'scheduler'} IS NOT CURRENTLY SUPOORTED";
    }
    
    #my $submit = "$cluster->{submit} $cluster->{job_name} $cluster->{job_hold} $cluster->{cpu} $cluster->{mem} $cluster->{internet} $addParams";
    #return $submit;
    
    return $cluster;
}

sub additionalParams{
    my %aparams = @_;

    my $addParams = '';
    if($aparams{'scheduler'} =~ /SGE/i){
        $addParams = additionalSGE(%aparams);
    }
    elsif($aparams{'scheduler'} =~ /LSF/i){
        $addParams = additionalLSF(%aparams);
    }

    return $addParams;

}

sub sge {
    my %sgeParams = @_;

    my $tcpus = 0;
    $cluster->submit("qsub");

    if($sgeParams{'job_name'} =~ /^[a-zA-Z]/){
	$cluster->job_name("-N $sgeParams{'job_name'}");
    }
    else{
	$cluster->job_name("");
    }

    if($sgeParams{'job_hold'} =~ /^[a-zA-Z]/){
	$cluster->job_hold("-hold_jid $sgeParams{'job_hold'}");
    }
    else{
	$cluster->job_hold("");
    }

    if($sgeParams{'cpu'} =~ /^\d+$/){
	$cluster->cpu("-pe alloc $sgeParams{'cpu'}");
	$tcpus = $sgeParams{'cpu'};
    }
    else{
	$cluster->cpu("-pe alloc $cluster->{cpu}");
	$tcpus = $cluster->{cpu};
    }

    if($sgeParams{'mem'} =~ /^\d+$/){
	###$cluster->mem("-l virtual_free=$sgeParams{'mem'}\G");
	### NOTE: sge mem request is per cpu
	my $memPerCPU = $sgeParams{'mem'}/$tcpus;
	$cluster->mem("-l virtual_free=$memPerCPU\G");
    }
    else{
	### $cluster->mem("-l virtual_free=$cluster->{mem}\G");
	my $memPerCPU = $cluster->{mem}/$tcpus;
	$cluster->mem("-l virtual_free=$memPerCPU\G");
    }

    if($sgeParams{'internet'}){
	$cluster->internet("-l internet=1");
    }

    if($sgeParams{'cluster_out'}){
	$cluster->cluster_out("-o $sgeParams{'cluster_out'}");
    }
}

sub additionalSGE {
    my %addSGE = @_;

    my $aSGE = "-j y -b y -V";

    if($addSGE{'work_dir'}){
	$aSGE .= " -wd $addSGE{'work_dir'}"
    }
    else{
	$aSGE .= " -cwd"
    }

    if($addSGE{'queues'}){
	$aSGE .= " -q $addSGE{'queues'}"
    }

    if($addSGE{'priority_project'}){
	$aSGE .= " -P $addSGE{'priority_project'}"
    }
    else{
	$aSGE .= " -P ngs"	
    }
    
    return $aSGE;
}

sub lsf {
    my %lsfParams = @_;

    $cluster->submit("bsub");

    if($lsfParams{'job_name'} =~ /^[a-zA-Z]/){
	$cluster->job_name("-J $lsfParams{'job_name'}");
    }
    else{
	$cluster->job_name("");
    }

    if($lsfParams{'job_hold'} =~ /^[a-zA-Z|,]/){
      ###$cluster->job_hold("-w \"post_done($lsfParams{'job_hold'})\"");
	my @jobs = split(/,/, $lsfParams{'job_hold'});
	my @holds = ();
	foreach my $job (@jobs){	    
	    if(!$job){
		next;
	    }
	    push @holds, "post_done($job)";
	    ###if($job =~ /\*$/){
	###	push @holds, "(post_done($job) || post_error($job))";
	   ### }
	    ###else{
	###	push @holds, "post_done($job)";
	   ### }
	}

	if(scalar(@holds) > 0){
	    my $hold = join(" && ", @holds);
	    $cluster->job_hold("-w \"$hold\"");
	}
	else{
	    $cluster->job_hold("");
	}
    }
    else{
	$cluster->job_hold("");
    }

    if($lsfParams{'cpu'} =~ /^\d+$/){
	$cluster->cpu("-n $lsfParams{'cpu'}");
    }
    else{
	$cluster->cpu("-n 24");
    }

    if($lsfParams{'mem'} =~ /^(\d+)[Gg]?$/){
	$cluster->mem("-R \"rusage[mem=$1]\"");
    }
    else{
	###$cluster->mem("-R \"rusage[mem=$cluster->{mem}]\"");
	$cluster->mem("-R \"rusage[mem=250]\"");
    }

    if($lsfParams{'internet'}){
	$cluster ->internet("-R \"select[internet]\"");
    }

    if($lsfParams{'cluster_out'}){
	$cluster ->cluster_out("-o $lsfParams{'cluster_out'}");
    }

    if($lsfParams{'cluster_error'}){
	$cluster ->cluster_error("-e $lsfParams{'cluster_error'}");
    }
}

sub additionalLSF {
    my %addLSF = @_;

    my $aLSF = "";

    if($addLSF{'rerun'}){
	$aLSF .= " -r -Q \"all ~0\"";
    }
    
    if($addLSF{'priority_group'}){
	$aLSF .= " -sla $addLSF{'priority_group'}";
    }
###    else{
###	$aLSF .= " -sla Pipeline";
###    }

    if($addLSF{'runtime'} =~ /^\d+$/){
	$aLSF .= " -We $addLSF{'runtime'}";
    }
    else{
	### defaults to a long job if runtime isn't provided
	$aLSF .= " -We 300";
    }

    if($addLSF{'work_dir'}){
	$aLSF .= " -cwd \"$addLSF{'work_dir'}\""
    }

    if($addLSF{'iounits'} =~ /^\d+$/){
	$aLSF .= " -R \"rusage[iounits=$addLSF{'iounits'}]\"";
    }
    else{
	$aLSF .=" -R \"rusage[iounits=10]\"";
    }

    if($addLSF{'mail'}){
	$aLSF .= " -u \"$addLSF{'mail'}\" -N"
    }

    return $aLSF;
}

1;
