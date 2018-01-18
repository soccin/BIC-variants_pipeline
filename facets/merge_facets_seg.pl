#!/usr/bin/perl 

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Basename;

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

my $facetsFilename = "$facets_dir/facets_mapping.txt";
if(! -e $facetsFilename){
    die "Cannot find the facets mapping file. $!";
}

open(my $fh, '<', $facetsFilename) or die "Cannot up file $facetsFilename: $!";

while (my $line = <$fh>){
    chomp $line;

    my $sDir = (split("\t", $line))[1];
    next if( $sDir eq "Rdata_filename");

    my $sample = dirname($sDir);

    if(! -d "$sample"){
        die "Directory $sample is not a real directory: $!";
    }

    opendir(my $dh, $sample);
    my @files = grep {/hisens\.seg$/ } readdir $dh;
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
