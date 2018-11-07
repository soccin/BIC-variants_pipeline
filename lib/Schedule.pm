#!/opt/bin/perl

package Schedule;

use strict;
use FindBin qw($Bin); 
use lib "$Bin/lib";
use Cluster;
use POSIX;
use List::Util qw(min max);

my $cluster = new Cluster;

sub queuing {
    my %qSystem = @_;

    my $addParams = '';
    if($qSystem{'scheduler'} =~ /SGE/i){
	sge(%qSystem);
	#$addParams = additionalSGE(%qSystem);
    }
    elsif($qSystem{'scheduler'} =~ /LUNA/i){
	luna(%qSystem);
	#$addParams = additionalLUNA(%qSystem);
    }
    elsif($qSystem{'scheduler'} =~ /JUNO/i){
    juno(%qSystem);
    #$addParams = additionalJUNO(%qSystem);
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
    elsif($aparams{'scheduler'} =~ /LUNA/i){
        $addParams = additionalLUNA(%aparams);
    }
    elsif($aparams{'scheduler'} =~ /JUNO/i){
        $addParams = additionalJUNO(%aparams);
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
	$cluster->mem("-l virtual_free=${memPerCPU}G");
    }
    else{
	### $cluster->mem("-l virtual_free=$cluster->{mem}\G");
	my $memPerCPU = $cluster->{mem}/$tcpus;
	$cluster->mem("-l virtual_free=${memPerCPU}G");
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


sub luna {
    my %lunaParams = @_;

    $cluster->submit("bsub");

    if($lunaParams{'job_name'} =~ /^[a-zA-Z]/){
    $cluster->job_name("-J $lunaParams{'job_name'}");
    }
    else{
    $cluster->job_name("");
    }

    if($lunaParams{'job_hold'} =~ /^[a-zA-Z|,]/){
      ###$cluster->job_hold("-w \"post_done($lunaParams{'job_hold'})\"");
    my @jobs = split(/,/, $lunaParams{'job_hold'});
    my @holds = ();
    foreach my $job (@jobs){        
        if(!$job){
        next;
        }
        push @holds, "post_done($job)";
        ###if($job =~ /\*$/){
    ### push @holds, "(post_done($job) || post_error($job))";
       ### }
        ###else{
    ### push @holds, "post_done($job)";
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

    if($lunaParams{'cpu'} =~ /^\d+$/){
    $cluster->cpu("-n $lunaParams{'cpu'}");
    }
    else{
    $cluster->cpu("-n 24");
    }

    if($lunaParams{'mem'} =~ /^(\d+)[Gg]?$/){
    $cluster->mem("-R \"rusage[mem=$1]\"");
    }
    else{
    ###$cluster->mem("-R \"rusage[mem=$cluster->{mem}]\"");
    $cluster->mem("-R \"rusage[mem=250]\"");
    }

    if($lunaParams{'internet'}){
    $cluster ->internet("-R \"select[internet]\"");
    }

    if($lunaParams{'cluster_out'}){
    $cluster ->cluster_out("-o $lunaParams{'cluster_out'}");
    }

    if($lunaParams{'cluster_error'}){
    $cluster ->cluster_error("-e $lunaParams{'cluster_error'}");
    }
}

sub additionalLUNA {
    my %addLUNA = @_;

    my $aLUNA = "";

    if($addLUNA{'rerun'}){
    $aLUNA .= " -r -Q \"all ~0\"";
    }
    
    if($addLUNA{'priority_group'}){
    $aLUNA .= " -sla $addLUNA{'priority_group'}";
    }
###    else{
### $aLUNA .= " -sla Pipeline";
###    }

    if($addLUNA{'runtime'} =~ /^\d+$/){
    $aLUNA .= " -We $addLUNA{'runtime'}";
    }
    else{
    ### defaults to a long job if runtime isn't provided
    $aLUNA .= " -We 300";
    }

    if($addLUNA{'work_dir'}){
    $aLUNA .= " -cwd \"$addLUNA{'work_dir'}\""
    }

    if($addLUNA{'iounits'} =~ /^\d+$/){
    $aLUNA .= " -R \"rusage[iounits=$addLUNA{'iounits'}]\"";
    }
    else{
    $aLUNA .=" -R \"rusage[iounits=10]\"";
    }

    if($addLUNA{'mail'}){
    $aLUNA .= " -u \"$addLUNA{'mail'}\" -N"
    }

    return $aLUNA;
}




sub juno {
    my %junoParams = @_;

    $cluster->submit("bsub");

    if($junoParams{'job_name'} =~ /^[a-zA-Z]/){
    $cluster->job_name("-J $junoParams{'job_name'}");
    }
    else{
    $cluster->job_name("");
    }

    if($junoParams{'job_hold'} =~ /^[a-zA-Z|,]/){
      ###$cluster->job_hold("-w \"post_done($junoParams{'job_hold'})\"");
    my @jobs = split(/,/, $junoParams{'job_hold'});
    my @holds = ();
    foreach my $job (@jobs){        
        if(!$job){
        next;
        }
        push @holds, "post_done($job)";
        ###if($job =~ /\*$/){
    ### push @holds, "(post_done($job) || post_error($job))";
       ### }
        ###else{
    ### push @holds, "post_done($job)";
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

    if($junoParams{'cpu'} =~ /^\d+$/){
    $cluster->cpu("-n $junoParams{'cpu'}");
    }
    else{
    $junoParams{'cpu'} = 24;
    $cluster->cpu("-n 24");
    }

    my $jobMem;
    if($junoParams{'mem'} =~ /^(\d+)[Gg]?$/){
    $jobMem = min(ceil($1 / $junoParams{'cpu'}), floor(240 / $junoParams{'cpu'}) );
    }
    else{
    $jobMem = floor(240 / $junoParams{'cpu'} );
    }
    $cluster->mem("-R \"rusage[mem=$jobMem]\"");

    if($junoParams{'internet'}){
    $cluster ->internet("-R \"select[internet]\"");
    }

    if($junoParams{'cluster_out'}){
    $cluster ->cluster_out("-o $junoParams{'cluster_out'}");
    }

    if($junoParams{'cluster_error'}){
    $cluster ->cluster_error("-e $junoParams{'cluster_error'}");
    }
}

sub additionalJUNO {
    my %addJUNO = @_;

    my $aJUNO = "";

    if($addJUNO{'rerun'}){
    $aJUNO .= " -r -Q \"all ~0\"";
    }
    
    #if($addJUNO{'priority_group'}){
    #$aJUNO .= " -sla DEVEL";
    #}
###    else{
### $aJUNO .= " -sla Pipeline";
###    }

    #if($addJUNO{'runtime'} =~ /^\d+$/){
    #$aJUNO .= " -We $addJUNO{'runtime'}";
    #}
    #else{
    ### defaults to a long job if runtime isn't provided
    $aJUNO .= " -W 300:00";
    #}

    $aJUNO .= " -app anyOS -R \"select[type==CentOS7]\"";

    if($addJUNO{'work_dir'}){
    $aJUNO .= " -cwd \"$addJUNO{'work_dir'}\""
    }

    #if($addJUNO{'iounits'} =~ /^\d+$/){
    #$aJUNO .= " -R \"rusage[iounits=$addJUNO{'iounits'}]\"";
    #}
    #else{
    #$aJUNO .=" -R \"rusage[iounits=10]\"";
    #}

    if($addJUNO{'mail'}){
    $aJUNO .= " -u \"$addJUNO{'mail'}\" -N"
    }

    return $aJUNO;
}



sub singularityBind{
    my $sBind = "/scratch/,/tmp/,/ifs/,/common/lsf/,/common/juno/";

    return $sBind;
}



sub singularityParams{
    my %sparams = @_;
 
    my $sinParams = "$sparams{'singularity_exec'} exec $sparams{'singularity_image'}";

    return $sinParams;
}
1;
