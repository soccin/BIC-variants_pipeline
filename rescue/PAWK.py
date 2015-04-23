#!/usr/bin/env python2.7

import csv
import sys

class CSVStruct(object):
    def __init__(self, cout, rec):
        self._cout=cout
        self._rec=rec
    def __getitem__(self,key):
        return self._rec[key]
    def __setitem__(self,key,value):
        self._rec[key]=value
    def __getattr__(self,key):
        return self._rec[key]
    def __setattr__(self,key,value):
        if key in ["_cout","_rec"]:
            object.__setattr__(self,key,value)
        else:
            self._rec[key]=value
    def write(self):
        self._cout.writerow(self._rec)
    def __repr__(self):
        out=[]
        for ki in self._rec:
            out.append("%s:%s" % (ki,self._rec[ki]))
        return "<Struct| "+", ".join(out)+">"

class PAWK(object):
    def __init__(self, fin=sys.stdin, fout=sys.stdout, header=True, delim="\t"):
        self.cin=csv.DictReader(fin,delimiter=delim)
        self.cout=csv.DictWriter(fout,self.cin.fieldnames,delimiter="\t")
        if header:
            self.cout.writeheader()
    def __iter__(self):
        return self
    def next(self):
        recDict=self.cin.next()
        return CSVStruct(self.cout, recDict)
