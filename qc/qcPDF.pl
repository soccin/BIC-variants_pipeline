#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 

###
### LSF AND SGE HANGLE QUOTES DEREFERENCING DIFFERENTLY
### SO STICKING IT IN A WRAPPER FOR SIMPLICITY
###

my ($pre, $path, $config);
GetOptions ('pre=s' => \$pre,
	    'path=s' => \$path,
	    'config=s' => \$config
	   ) or exit(1);

if(!$pre || !$path || !$config){
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

print $R;
print $Bin;

`$R/R CMD BATCH "--args path='$path' pre='$pre' bin='$Bin'" $Bin/qcPDF.R`;
