#!/usr/bin/env python2.7

import csv
import sys

class Struct:
    def __init__(self, **entries): 
        self.__dict__.update(entries)
"""
FILTER SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT
"""

cin=csv.DictReader(sys.stdin,delimiter="\t")
cout=csv.DictWriter(sys.stdout,cin.fieldnames,delimiter="\t")
cout.writeheader()
for recDict in cin:
  try:
    rec=Struct(**recDict)
    if rec.SNPEFF_FUNCTIONAL_CLASS=="." or rec.SNPEFF_FUNCTIONAL_CLASS=="SILENT" or rec.SNPEFF_IMPACT=="MODIFIER":
        continue
    cout.writerow(recDict)
  except:
    print recDict
    print
    print
    raise
