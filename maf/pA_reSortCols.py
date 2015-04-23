#!/usr/bin/env python2.7

import csv
import sys

flds=open(sys.argv[1]).read().strip().split()

cin=csv.DictReader(sys.stdin,delimiter="\t")
cout=csv.DictWriter(sys.stdout,cin.fieldnames,delimiter="\t")
print "\t".join(flds)
for rec in cin:
    out=[rec[x] for x in flds]
    print "\t".join(out)