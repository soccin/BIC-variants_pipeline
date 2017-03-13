#!/usr/bin/env python2.7
import os
import csv
import sys
import argparse
from itertools import dropwhile


def checkFile(afile):
    if not os.path.isfile(afile):
        print("".join(["[ERROR] File ", afile, " does not exist. Please check your options and try again."]));
        sys.exit(-1)

## Take in OPtions:
parser=argparse.ArgumentParser()
parser.add_argument("-i","--inputFile",help="Input filename",required=True)
parser.add_argument("-o","--outFile",help="Output filename", required=True)
parser.add_argument("-f","--fields",help="Fields that you want to keep", required=True)
args=parser.parse_args()

checkFile(args.inputFile)
checkFile(args.fields)

flds=open(args.fields).read().strip().split()

headerCount=0
outf = open(args.outFile, 'w')
with open(args.inputFile, 'r') as f:
    for line in f:
        if line.startswith("#"):
            outf.write(line)
    f.seek(0)
    
    start = dropwhile(lambda L: L.lower().lstrip().startswith('#'),f)
    cin=csv.DictReader(start,delimiter="\t")

    outf.write("\t".join(flds) + "\n")

    print "FIELDNAMES: " + "\t".join(cin.fieldnames)

    for rec in cin:
        out=[str(rec[x]) for x in flds]
        outf.write("\t".join(out) + "\n")

    outf.close()






