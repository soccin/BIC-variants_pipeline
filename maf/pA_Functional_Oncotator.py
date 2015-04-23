#!/usr/bin/env python2.7

import csv
import sys

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)
"""
FILTER
ONCOTATOR_VARIANT_CLASSIFICATION
ONCOTATOR_PROTEIN_CHANGE
ONCOTATOR_COSMIC_OVERLAPPING
ONCOTATOR_DBSNP_RS
ONCOTATOR_GENE_SYMBOL
"""

cin=csv.DictReader(sys.stdin,delimiter="\t")
cout=csv.DictWriter(sys.stdout,cin.fieldnames,delimiter="\t")
cout.writeheader()
for recDict in cin:
  try:
    rec=Struct(**recDict)
    ###if rec.ONCOTATOR_VARIANT_CLASSIFICATION in ["Intron","Silent","3'UTR","5'UTR","5'Flank","3'Flank"]:
       ### continue
    if rec.ONCOTATOR_GENE_SYMBOL == "":
        continue
    cout.writerow(recDict)
  except:
    print recDict
    print
    print
    raise
