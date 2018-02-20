#!/usr/bin/perl

use strict;
use Getopt::Long qw(GetOptions);


my ($inputVcf, $outputVcf, $log, $help);
GetOptions ('input|in=s' => \$inputVcf,
	    'output|out=s' => \$outputVcf,
            'log=s' => \$log) or exit(1);

if(!$inputVcf || !$outputVcf || $help){
    print <<HELP;

    *** B37 has non ATGC nucleotides in assembly (e.g. chr 3 has an M base)
    *** ONLY known variant caller to include these positions so far is lancet
    *** gatk fails to merge chromosome split lancet calls because the vcf is considered malformed due to these positions
    *** this script removes those positions from the vcf file

    USAGE: vcfFilterLancet.pl -input INPUT -output -OUTPUT -log LOG
	* INPUT: vcf for processing (REQUIRED)
	* OUTPUT: output file name (REQUIRED)
	* LOG: discarded snp calls
HELP
exit;
}

open(IN, "$inputVcf") or die "Can't open input file $inputVcf $!";
open(OUT, ">$outputVcf") or die "Can't write output file $outputVcf $!";
open(LOG, ">$log") or die "Can't write log file $log $!";
while(my $line=<IN>){
    chomp $line;

    if($line =~ /^\#/){
	print OUT "$line\n";
	next;
    }

    my @data = split(/\t/, $line);

    if($data[3] =~ /[^atgc]/i){
	print LOG "NON_STANDARD_CHARS\t$line\n";
    }
    else{
	if ($data[5] eq 'inf'){
	    ### LANCET SOMETIMES OUTPUT "inf" in qual field, which gatk doe not like
	    ### gatk error: java.lang.NumberFormatException: For input string: "inf"

	    print LOG "INF_QUAL\t$line\n";
	}
	else{
	    print OUT "$line\n";
	}
    }

}
close IN;
close OUT;
close LOG;
