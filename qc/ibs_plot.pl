#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

###
### LSF AND SGE HANGLE QUOTES DEREFERENCING DIFFERENTLY
### SO STICKING IT IN A WRAPPER FOR SIMPLICITY
###

my ($pre, $pair, $ibs_file, $config);
GetOptions ('pre=s' => \$pre,
	    'pair=s' => \$pair,
	    'config=s' => \$config,
	    'ibs_file=s' => \$ibs_file) or exit(1);

if(!$pre || !$pair || !$ibs_file || !$config){
    die "MUST PROVIDE PRE, PAIR, CONFIG AND OUTPUT\n";
}

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
}
close CONFIG;

`$R/Rscript $Bin/ibs_plot.R "projID='$pre'" "manifestPairing='$pair'" "IBSFile='$ibs_file'"`;
