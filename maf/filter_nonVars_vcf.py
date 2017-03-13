#!/usr/bin/env python2.7

## This script will take in vcf file that was split from a larger file (haplotect)
## If this is a somatic variant calling, I have to skip if either GT is "./.", and if GT is the same in tumor and normal
##    If Tumor is 0/0, skip (we are not keeping LOH) 
## If this is NON - Somatic, the script will filter variants where the tumor is 0/0 or ./. GT.


## Continue if
##    1: no GT or AD in INFO field for record
##    2: GT Constraints above
##    3: If for some reason AD is . or 0,0


## NOTE: I am assuming that GT is the first item in the FORMAT field. It would be nice to figure that out programmatically just incase
##       something crazy happens.

import csv
import sys
import argparse
from itertools import ifilter,dropwhile
        
parser = argparse.ArgumentParser(description="MAF filtering script")
parser.add_argument('-v', '--vcf', action='store', help="vcf file to be filtered")
parser.add_argument('-s', '--somatic', action='store_true', help="Flag for using somatic filtering")
parser.add_argument('-t', '--tumor', action='store', help='tumor id')
parser.add_argument('-n', '--normal', action='store', help='normal id')
parser.add_argument('-o', '--output', action='store', help="Output maf filename")
args = parser.parse_args()

# somewhere to store the comments
#comments=list()

with open(args.vcf) as fin:
    ofh = open(args.output, "wt")
    # first take the comments and print them to the output file
    comments = ifilter(lambda L: L.startswith("##"), fin)
    for line in comments:
        ofh.write(line)
        
    fin.seek(0)
    
    # Now grab formatted rest of file, change ofh to a csv.dictwriter, output header
    start = dropwhile(lambda L: L.lower().lstrip().startswith('##'), fin)
    cin = csv.DictReader(start, delimiter='\t')
    cout = csv.DictWriter(ofh,fieldnames=cin.fieldnames,delimiter="\t")
    cout.writeheader()
    for row in cin:
        
        # If row does not have GT or AD continue
        formatRow = row["FORMAT"].split(":")
        if "GT" not in formatRow or "AD" not in formatRow:
            continue
        
        gtIndex = formatRow.index("GT")
        adIndex = formatRow.index("AD")
            
        tumorGT = row[args.tumor].split(":")[gtIndex]
        tumorAD = row[args.tumor].split(":")[adIndex]
        normalGT = row[args.normal].split(":")[gtIndex]
        normalAD = row[args.normal].split(":")[adIndex]
        
        # filter for non-paired and paired analysis
        if tumorGT == "0/0" or tumorGT == "./.":
            continue
        if tumorAD == ".":
            continue
        
        # only for paired analysis
        if args.somatic : 
            if not normalGT == "0/0" :
                continue
            if normalAD == ".":
                continue
            
        
        
        ## If it has made it this far, print the line to output.
        cout.writerow(row)
    
    

