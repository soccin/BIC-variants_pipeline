#!/usr/bin/env python2.7

import csv
import sys

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)
"""
FILTER SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT
"""

class Globals:
  pass

glbs=Globals()
glbs.MinTumorDepth=20
glbs.MinNormalDepth=8
glbs.MultipleTumorNormalFreq=5
glbs.MinAD=5
glbs.MinAltFreq=0.01

cin=csv.DictReader(sys.stdin,delimiter="\t")
cout=csv.DictWriter(sys.stdout,cin.fieldnames,delimiter="\t")
cout.writeheader()
for recDict in cin:
  try:
    rec=Struct(**recDict)
    if rec.FILTER and rec.FILTER.find("LowQual")>-1:
        continue
    if not rec.GT in ["0/0","./."] and rec.NORM_GT=="0/0":
        tumorDepth=int(rec.AD_ALT)+int(rec.AD_REF)
        normalDepth=int(rec.NORM_AD_REF)+int(rec.NORM_AD_ALT)
        if tumorDepth>=glbs.MinTumorDepth and normalDepth>=glbs.MinNormalDepth \
           and float(rec.ALT_FREQ)>glbs.MultipleTumorNormalFreq*float(rec.NORM_ALT_FREQ) \
           and int(rec.AD_ALT) > glbs.MinAD \
           and float(rec.ALT_FREQ) > glbs.MinAltFreq:
            cout.writerow(recDict)

  except:
    print recDict
    print
    print
    raise
