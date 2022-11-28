#!/usr/bin/env python2.7

import csv
import sys

class Struct:
    def __init__(self, entries):
        if None in entries:
            del entries[None] 
        self.__dict__.update(**entries)

"""
FILTER SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT
"""

cin=csv.DictReader(sys.stdin,delimiter="\t")
cout=csv.DictWriter(sys.stdout,cin.fieldnames,delimiter="\t",lineterminator="\n")
cout.writeheader()
for recDict in cin:
    try:
        rec=Struct(recDict)
        if rec.FILTER.find("LowQual")>-1:
            continue
        if not rec.GT in ["0/0","./.",None] and rec.AD_ALT != "." and int(rec.AD_ALT)>4:
            cout.writerow(recDict)
    except:
        print recDict
        print
        print recDict.keys()
        raise
