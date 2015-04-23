#!/usr/bin/env python2.7

from PAWK import PAWK
import sys

for rec in PAWK(sys.stdin):
    rec.Hugo_Symbol=rec.ONCOTATOR_GENE_SYMBOL
    rec.write()
