#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Path qw(make_path remove_tree);
use File::Basename;

my ($inMaf, $outMaf, $config, $species, $blist, $lowConf, $ffpe );
my @headers = ();

GetOptions('in_maf=s' => \$inMaf,
'species=s' => \$species,
'out_maf=s' => \$outMaf,
'config=s' => \$config,
'blacklist' => \$blist,
'low_conf' => \$lowConf,
'filter_ffpe' => \$ffpe) or exit(1);

my $uID = `/usr/bin/id -u -n`;
chomp $uID;

#Mandatory input
if( !$inMaf || !$outMaf || !$species  || !$config){
    die "You are missing necessary input";
}

#make sure they have at least one filter
if( !$blist && !$lowConf && !$ffpe ){
    die "You haven't picked a filter to do.";
}

#check species
#Hybrid is okay to run because we should only be looking at the human chrs
if($species !~ /b37|hg19|human|hybrid/i){
    die "You cannot run any of these filters with this species";
}

# Config
my $R = '';
my $WES_FILTER = '';

open(CONFIG, "$config") or warn "CAN'T OPEN CONFIG FILE $config SO USING DEFAULT SETTINGS";
while(<CONFIG>){
    chomp;
    my @conf = split(/\s+/, $_);
    if($conf[0] =~ /^r$/i){
        if(!-e "$conf[1]/R"){
            die "CAN'T FIND R IN $conf[1] $!";
        }
        $R=$conf[1];
        my $path_tmp = $ENV{'PATH'};
        $ENV{'PATH'} = "$conf[1]:$path_tmp";
    }
    elsif($conf[0] =~ /wes_filter/i){
        if(!-e "$conf[1]/filter_blacklist_regions.R"){
            die "CAN'T FIND wes-filter IN $conf[1] $!";
        }
        $WES_FILTER = $conf[1];
    }
}
close CONFIG;

## die if these are empty
if ( $WES_FILTER eq "" || $R eq "" ){
    die "config file is missing values";
}

my $tempOut='';
my $headers= `grep ^# $inMaf`;
foreach my $head (split(/\n/,$headers)){
    push @headers, $head;
}

my $tempIn=$inMaf;
my @delete_me_later = ();



#### START OF WES FILTERS

if($lowConf){
    $tempOut="$tempIn\_lowConf";
    print "\n$R/Rscript $WES_FILTER/filter_low_conf.R -m $tempIn -o $tempOut\n";
    `$R/Rscript $WES_FILTER/filter_low_conf.R -m $tempIn -o $tempOut`;
    $tempIn=$tempOut;
    push @headers, "#WES-FILTER: filter_low_conf.R \tLOCATION: $WES_FILTER";
    push @delete_me_later, $tempOut;
}

if($blist){
    $tempOut="$tempIn\_blacklist";
    print "\n$R/Rscript $WES_FILTER/filter_blacklist_regions.R -m $tempIn -o $tempOut\n";
    `$R/Rscript $WES_FILTER/filter_blacklist_regions.R -m $tempIn -o $tempOut`;
    $tempIn=$tempOut;
    push @headers, "#WES-FILTER: filter_blacklist_regions.R \tLOCATION: $WES_FILTER";
    push @delete_me_later, $tempOut;
}

if($ffpe){
    $tempOut="$tempIn\_filter_ffpe";
    print "\n$R/Rscript $WES_FILTER/filter_ffpe.R -m $tempIn -o $tempOut\n";
    `$R/Rscript $WES_FILTER/filter_ffpe.R -m $tempIn -o $tempOut`;
    $tempIn=$tempOut;
    push @headers, "#WES-FILTER: filter_ffpe.R \tLOCATION: $WES_FILTER";
    push @delete_me_later, $tempOut;
}

#### END OF WES FILTERS

##add Header
if( $headers[0] !~ /^#version/){
    unshift(@headers, "#version 2.4");
}

## Print headers
open ( my $fh, '>', "$outMaf\_header.txt");
while(@headers){
    print $fh shift @headers;
    print $fh "\n"; 
}
close $fh;
push @delete_me_later, "$outMaf\_header.txt";


# only save the last one
`cat $outMaf\_header.txt $tempOut > $outMaf`;

## Delete the delete_me_later stuff
unlink @delete_me_later;

