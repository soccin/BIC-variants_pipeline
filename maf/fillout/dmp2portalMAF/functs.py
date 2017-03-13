import TCGA_MAF
from collections import defaultdict
from lib import *
import csv
import os.path
import sys
from itertools import ifilter,dropwhile
## ifilter is not being used yet. If you want to save the comments, you need to pull them and store them with ifilter

validCallers=["hc", "haplotypecaller", "mt", "mutect", "haplotect"] # haplotect is special.  INDELs == haplotype caller, SNPS == muTect (untill I fix)
validSpOptions=["hg19", "b37", "mm10"]

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)
"""
FILTER SNPEFF_EFFECT SNPEFF_FUNCTIONAL_CLASS SNPEFF_IMPACT
"""


def validFile(filename):
    if not os.path.exists(filename):
        print>>sys.stderr,"".join(["Error, ",filename," does not exists. Please check file and rerun script."])
        sys.exit(1)

def validCaller(c):
    if c not in validCallers:
        print >> sys.stderr, "".join(["[ERROR] ", c, " is not a valid caller. Valid caller are: ", ", ".join(validCallers), "\n"])
        sys.exit(1) 

def validSpecies(s):
    if s not in validSpOptions:
        print >> sys.stderr, "".join(["[ERROR] ", s, " is not a valid species. Valid species are: ", ", ".join(validSpOptions), "\n"])
        sys.exit(1)

# Get tumor normal pairs  
def populatePairingInfo(PAIRINGFILE):
    pairs=set()
    with open(PAIRINGFILE) as fp:
        for line in fp:
            (normal,tumor)=line.strip().split("\t")
            pairs.add(";".join([tumor,normal]))
    return pairs

def populatePatientInfo(PATIENTFILE):
    sampleDb=dict()
    baits=set()
    with open(PATIENTFILE) as fp:
        cin=csv.DictReader(fp,delimiter="\t")
        for rec in bunchStream(cin):
            sampleDb[rec.Sample_ID]=rec
            baits.add(rec.Bait_version)
            
    if len(baits)>2:
        print >>sys.stderr, "Multiple Baits", baits
        sys.exit(1)
        
    return sampleDb

# Grab a bunch of information from original maf so I can match the maf record
# to what the tumor seq allele1 was in the original maf. If the records are
# fillouts, the tumor seq allele1 will ALWAYS be the reference sequence.
def populateOriginalMafInfo(mafFile):
    origMafInfo=dict()
    fp = open(mafFile)
    start=dropwhile(lambda L: L.lower().lstrip().startswith('#'), fp)
    cin=csv.DictReader(start,delimiter="\t")
    for rec in cin:
        rec=Struct(**rec)
        tum=rec.Tumor_Sample_Barcode
        norm=rec.Matched_Norm_Sample_Barcode
        chr=rec.Chromosome
        start=rec.Start_Position
        end=rec.End_Position
        varType=rec.Variant_Type
        ref=rec.Reference_Allele
        tumor1=rec.Tumor_Seq_Allele1
        tumor2=rec.Tumor_Seq_Allele2
        
        origMafInfo["-".join([tum,norm,chr,start,end,varType,ref,tumor2])]=tumor1
            
    fp.close()
    return origMafInfo
    
def generatePatientGroups(sampleDb):
    patientGroups=defaultdict(set)
    for sa in sampleDb.keys():
        patientGroups[sampleDb[sa].Patient_ID].add(sa)
    return patientGroups
