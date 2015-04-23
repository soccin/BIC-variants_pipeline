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
    ## If the value in ONCOTATOR VARIANT CLASSIFICATION is not in the list over there, don't print anything
    if not rec.ONCOTATOR_VARIANT_CLASSIFICATION in ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site"]:
       continue
    cout.writerow(recDict)
  except:
    print recDict
    print
    print
    raise
