#!/usr/bin/perl 

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);

my ($facets_dir, $seg_output);

my $commandline = join " ", "\n", $0, @ARGV;

print "$commandline\n\n";

GetOptions ('facets_dir=s' => \$facets_dir,
            'outfile=s' => \$seg_output ) or exit(1);

if(!$seg_output || !$facets_dir){
    print <<HELP;
    
    USAGE: <name of this script> -facets_dir <path> -outfile <merged_seg_filename>
    
    This script will merge all facets seg files.
    You must specify:
    
    -facets_dir		The relative/absolute path to the facets path from the CWD
    -outfile		The output filename (and the relative path to it) It does not to have to be made
    
HELP
    exit;
}

if( -f "$seg_output"){
    unlink("$seg_output");
}

my @sampleDirs = `ls -d $facets_dir/*`;

foreach my $sample (@sampleDirs){
    chomp $sample;
    #skip everything except directories
    next if(! -d "$sample");

    opendir(my $dh, $sample);
    my @files = grep {/hisens.seg$/ } readdir $dh;
    closedir $dh;

    if (! @files) {
        die "Cannot find seg file in $sample directory.";
    }
    if(scalar @files > 1){
        die "More than one 'hisens.seg' file found in directory $sample. ";
    }
    my $segFile = $files[0];
    
    if(! -f "$seg_output") {
        `head -1 $sample/$segFile > $seg_output`;
    }
    `tail -n +2 $sample/$segFile >> $seg_output`;
    
}
