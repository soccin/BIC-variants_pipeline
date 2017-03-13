#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use File::Basename;
use List::MoreUtils qw(first_index);

# This script should add a sample called "REF.<species>" to the end of the vcf file, with empty INFO values.
my ($in_vcf, $out_vcf, $species);

GetOptions( 'in_vcf=s' => \$in_vcf,
            'species=s' => \$species,
            'out=s' => \$out_vcf
) or die;

if(! $in_vcf){
    die "You need to include a in_vcf file\n";
}
if(!-e $in_vcf){
    die "$in_vcf not found!\n";
}
if(! $out_vcf){
    die "You need to include an -out file\n";
}
if(-e $out_vcf){
    print "Out file $out_vcf already exists. This will be overwritten.\n";
    unlink($out_vcf);
}
if(! $species){
    die "You need to include a -species value\n";
}

if($species != /b37|hg19|mm10|mouse/i){
    print "This will not be done, because we cannot do VEP for this species: $species\n";
    exit 0;
}

# Okay, now this is pretty simple. Open file
# Print "^##" without changing anything
# find "^#CHROM" add sample to the end of the header
# each line after, grab the FORMAT field, and figure out what the added info will be for the new sample

open(INVCF, "$in_vcf");
open(OUT, ">$out_vcf");
my $formatIndex;
while(my $line = <INVCF>){
    chomp $line;
    if($line =~ /^##/){
        print OUT $line . "\n";
    }elsif($line =~ /^#CHROM/i){
        # Find index for FORMAT
        my @header = split(/\s+/,$line);
        $formatIndex = first_index { $_ eq 'FORMAT' } @header;
        if (! $formatIndex){
            die "Could not find FORMAT index in header!";
        }
        print OUT $line . "\tREF." . $species . "\n";
    }else{
        # Assuming this is an actual variant record
        my @outputString;
        my @rec=split(/\s+/, $line);
        my @format=split(":",$rec[$formatIndex]);
        for my $item (@format){
            ## Once I get more formats that are different than '.', I will add them here
            if($item =~ 'GT'){
                push @outputString, "./.";
            }else{
                push @outputString, ".";
            }
        }
        print OUT $line . "\t" . join(":", @outputString) . "\n";
    }
}
close(INVCF);
close(out);

