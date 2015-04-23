#!/usr/bin/env python2.7

import csv
import sys

class Struct:
    def __init__(self, **entries): 
        self.__dict__.update(entries)
"""
FILTER SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT
"""

DELTA_NRAF=0.10

cin=csv.DictReader(sys.stdin,delimiter="\t")
cout=csv.DictWriter(sys.stdout,cin.fieldnames,delimiter="\t")
cout.writeheader()
for recDict in cin:
  try:
    rec=Struct(**recDict)
    if rec.FILTER and rec.FILTER.find("LowQual")>-1:
        continue
    if not rec.GT in ["0/0","./."] and rec.NORM_GT=="0/0":
        if int(rec.AD_ALT)>2 and int(rec.NORM_AD_REF)>2 \
           and float(rec.ALT_FREQ)>float(rec.NORM_ALT_FREQ)+DELTA_NRAF and float(rec.ALT_FREQ)>.15:        
            cout.writerow(recDict)

  except:
    print recDict
    print
    print
    raise
