#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

###
### LSF AND SGE HANDLE QUOTES DEREFERENCING DIFFERENTLY
### SO STICKING IT IN A WRAPPER FOR SIMPLICITY
###

my ($pre, $path, $config, $log, $request, $version);
GetOptions ('pre=s' => \$pre,
	    'path=s' => \$path,
	    'config=s' => \$config,
            'log=s' => \$log,
            'request=s' => \$request,
            'version=s' => \$version
	   ) or exit(1);

if(!$pre || !$path || !$config || !$request || !$log || !$version){
    die "MUST PROVIDE PRE, PATH, CONFIG, REQUEST FILE, SVN REVISION NUMBER AND OUTPUT\n";
}

my $R = '';
my $JAVA = '';
open(CONFIG, "$config") or die "CAN'T OPEN CONFIG FILE $config $!";
while(<CONFIG>){
    chomp;
    
    my @conf = split(/\s+/, $_);

    if($conf[0] =~ /^r$/i){
        if(!-e "$conf[1]/R"){
            die "CAN'T FIND R IN $conf[1] $!";
        }
        $R = $conf[1];
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /^java$/i){
        if(!-e "$conf[1]/java"){
            die "CAN'T FIND java IN $conf[1] $!";
        }
        $JAVA = $conf[1];
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
}
close CONFIG;

## generate a PDF file for each plot, a project summary text file and a sample summary text file
#$R = '/opt/common/CentOS_7/R/R-3.2.0/bin'; ######## TEMPORARY
print "$R/R CMD BATCH \"--args path='$path' pre='$pre' bin='$Bin' logfile='$log'\" $Bin/qc_summary.R\n";
`$R/R CMD BATCH \"--args path='$path' pre='$pre' bin='$Bin' logfile='$log'\" $Bin/qc_summary.R`;

my $ec = $? >> 8;

## generate the complete, formal PDF report
print "$JAVA/java -jar $Bin/QCPDF.jar -rf $request -v $version -d $path -o $path -pl Variants\n"; 
`$JAVA/java -jar $Bin/QCPDF.jar -rf $request -v $version -d $path -o $path -pl Variants`;

my $ec2 = $? >> 8;

if($ec != 0 || $ec2 != 0){ die; }


