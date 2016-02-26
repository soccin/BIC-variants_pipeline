#!/usr/bin/perl
# This script will merge multiple maf files

use strict;
use warnings;
use Getopt::Long;
my @files;
my $output = "";

GetOptions(
               'i=s' => \@files,
               'o=s' => \$output
           );
       
       
       
# Different things that could be in the comments section
my $version;
my $extra_comments = "";
my $header;
           
open(my $TEMP, ">$output\_TEMP.txt") || die "Can't open temp output file $output\_TEMP.txt $!";
foreach (@files){
    my ($file) = $_;
    #print "FILE: $file\n";
    if (-z $file) {
        die "[ERROR] File $file is empty. If there were truely no variants, the file should still contain the header.";
    }
    my $header_found = 0;
    open(my $FH, "$file") || die "Can't open file $file $!";
    while(<$FH>){
        chomp;
        my $line = $_;
        if($line =~ /^#/){
            if($line =~ /^#\s*version\s*\d*.\d*$/){
                if(! $version){
                    $version = $line;
                }else{
                    my ($temp_v) = $version =~ tr/ +//;
                    my ($temp_l) = $line =~ tr/ +//;
                    if($temp_v != $temp_l){
                        die "[ERROR] Versions of maf don't match: $version $line";
                    }
                }
            }else{
                # This means it's a random comment. 
                print "EXTRA COMMENT $line";
                # TODO: Make sure it isn't in @extra_comments, if not, then add it.
                if( index($extra_comments, $line) == -1 ){
                    $extra_comments .= "$line\n";
                }
            }
        } elsif( $line =~ /^Hugo_Symbol/){
            $header_found=1;
            if(! $header){
                $header = $line;
            }elsif($header ne $line){
                die "[ERROR] Headers differ: \nHeader So Far:\n$header\nHeader in file $file:\n$line";
            }
        } else {
            # This means we are in the actual variants. Make sure we have a header. Then just print it to a temp file.
            
            if(! $header_found){
                die "[ERROR] Header was not found in current file: $file. This is a malformed maf. ";
            }
            ## PRINT OUT
            print $TEMP "$line\n";
        }
    }
}
close($TEMP);
## Here Print the COMMENTS
#Then add thme all together

open(my $OUT, ">$output") || die "Can't open output file $!";
print $OUT "$version\n" if ($version);
print $OUT $extra_comments if ($extra_comments);
print $OUT "$header\n";
close($OUT);

`cat $output\_TEMP.txt >> $output`;
`rm $output\_TEMP.txt`;



