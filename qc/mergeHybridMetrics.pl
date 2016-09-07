#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin); 
use File::Basename;


### takes in picard metrics (tested on hs/is/as/md
### through -metrics
### and/or can give it a file of files (-files)

my $pre = 'TEMP';
my @metrics = ();
my $files = '';
GetOptions ('pre=s' => \$pre,
	    'files=s' => \$files,
	    'metrics|m=s' => \@metrics) or exit(1);

if($files){
    open(FI, "$files") or die "Can't open file $files $!";
    while(<FI>){
	chomp;

	push @metrics, $_;
    }
    close FI;
}

foreach my $met (@metrics){
    if(!-s $met){
	die "METRICS FILE $met DOES NOT EXIST OR IS EMPTY";
    }
}



print "SAMPLE_NAME\tTOTAL_READ_PAIRS\tREAD_PAIRS_AFTER_FILTERING\tPCT_HUMAN_READS\n";
foreach my $met(@metrics)
{
        my $statfilename = basename($met);
        my ($samplename) = $statfilename =~ /.*_HybridFilter_(.*)\.txt/;
        
        my $prefilter_count;
        my $afterfilter_count;
        my $kept_percentage;
        
        open(MET, "$met") or die "[ERROR]: Can't open hybrid stats file $!\n";
        while(<MET>)
        {
                chomp($_);
                my @data = split(/:/);
                if($data[0] =~ /^Total read pairs before filtering$/)
                {
                        $prefilter_count = $data[1];
                }
                elsif($data[0] =~ /^Total read pairs kept after filtering$/)
                {
                        $afterfilter_count = $data[1];
                }
                elsif($data[0] =~ /^Percentage of read pairs kept$/)
                {
                        $kept_percentage = $data[1];
                }
        }
        close MET;
        print "$samplename\t$prefilter_count\t$afterfilter_count\t$kept_percentage\n";
}







