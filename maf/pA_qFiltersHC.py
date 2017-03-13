#!/usr/bin/env python2.7

import csv
import sys
import argparse
from itertools import ifilter,dropwhile


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)
"""
FILTER SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT
"""

class Globals:
  pass

parser = argparse.ArgumentParser(description="MAF filtering script")
parser.add_argument('-m', '--maf', action='store', help="Maf file to be filtered", required=True)
parser.add_argument('-s', '--somatic', action='store_true', help="Flag for using somatic filtering")
parser.add_argument("-c","--caller",help="Which snp calling program used", required=True)
parser.add_argument('-o', '--output', action='store', help="Output maf filename", required=True)
args = parser.parse_args() 

## VALUES TO FILTER WITH
## SOMATIC GETS MORE VALUES
glbs=Globals()
# unpaired/ not somatic global values
# make extra values negative 1 so they are always true
glbs.MinAD=4
glbs.MinTumorDepth=-1
glbs.MinNormalDepth=-1
glbs.MultipleTumorNormalFreq=-1
glbs.MinAltFreq=-1
if args.somatic:
    glbs.MinTumorDepth=20
    glbs.MinNormalDepth=8
    glbs.MultipleTumorNormalFreq=5
    glbs.MinAD=5
    glbs.MinAltFreq=0.01

fin=open(args.maf)
ofh=open(args.output,"wt")

comments=ifilter(lambda L: L.startswith("#"), fin)
for line in comments:
    ofh.write(line)
    
fin.seek(0)

start=dropwhile(lambda L: L.lower().lstrip().startswith('#'), fin)
cin=csv.DictReader(start,delimiter="\t")
fnames = cin.fieldnames
fnames.append("Caller")
cout=csv.DictWriter(ofh,fieldnames=fnames,delimiter="\t")
cout.writeheader()

for recDict in cin:
  try:
    recDict["Caller"]=args.caller
    rec=Struct(**recDict)
    
    ## For both somatic and not somatic
    if rec.FILTER and rec.FILTER.find("LowQual")>-1:
        continue
    
    normalDepth=0
    normalFreq=0
    tumorDepth=int(rec.t_alt_count) +int(rec.t_ref_count)
    if tumorDepth <= 0:
        continue
    tumorFreq=float(rec.t_alt_count)/tumorDepth
    
    # Sometimes they may not be filled out
    if rec.n_alt_count and rec.n_ref_count:
        normalDepth=int(rec.n_alt_count)+int(rec.n_ref_count)
        if normalDepth > 0:
            normalFreq =float(rec.n_alt_count)/tumorDepth
        
    if tumorDepth>=glbs.MinTumorDepth and normalDepth>=glbs.MinNormalDepth \
       and tumorFreq>glbs.MultipleTumorNormalFreq*normalFreq \
       and int(rec.t_alt_count) > glbs.MinAD \
       and tumorFreq > glbs.MinAltFreq:
        cout.writerow(recDict)

  except:
    print recDict
    print
    print
    raise



