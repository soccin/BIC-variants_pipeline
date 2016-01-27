#!/usr/bin/env python2.7

import sys
import csv
import copy
import argparse
import os.path
from collections import defaultdict
from lib import *

validCallers=["hc", "haplotypecaller", "mt", "mutect", "haplotect"] # haplotect is special.  INDELs == haplotype caller, SNPS == muTect (untill I fix)
validSpOptions=["hg19", "b37", "mm10"]

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

parser=argparse.ArgumentParser()
parser.add_argument("-m","--maf",help="Original Maf File", required=True)
parser.add_argument("-p","--pairing",help="Pairing File")
parser.add_argument("-P","--patient",help="Patient File", required=True)
parser.add_argument("-b","--baseCounts",help="Base Counts File", required=True)
parser.add_argument("-o","--output",help="Output Maf Filename", required=True)
parser.add_argument("-c","--caller",help="Caller that generated variants", required=True)
parser.add_argument("-s","--species", help="Species used", required=True)
args=parser.parse_args()

# verify maf, pairing, and patient are all files.
validCaller(args.caller)
validFile(args.maf)
validFile(args.patient)
validSpecies(args.species)

### For human, if species is hg19 or b37, NCBI value is GRCh37
### For mouse it is GRCm38
NCBI="GRCh37"
if args.species == "mm10":
    NCBI="GRCm38"

pairs=dict()

if args.pairing:
    validFile(args.pairing)
    PAIRINGFILE=args.pairing

    # Get tumor normal pairs
    with open(PAIRINGFILE) as fp:
        for line in fp:
            (normal,tumor)=line.strip().split("\t")
            pairs[tumor]=normal


PATIENTFILE=args.patient

origMafInfo=dict()
# Grab a bunch of information from original maf so I can match the maf record
# to what the tumor seq allele1 was in the original maf. If the records are
# fillouts, the tumor seq allele1 will ALWAYS be the reference sequence.
with open(args.maf) as fp:
    cin=csv.DictReader(fp,delimiter="\t")
    for rec in bunchStream(cin):
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

# Take in patient file
sampleDb=dict()
baits=set()
with open(PATIENTFILE) as fp:
    cin=csv.DictReader(fp,delimiter="\t")
    for rec in bunchStream(cin):
        baits.add(rec.Bait_version)
        sampleDb[rec.Sample_ID]=rec

# If there are multiple baitsets, error!
if len(baits)>2:
    print >>sys.stderr, "Multiple Baits", baits
    sys.exit(1)

CenterTag="MSK-"+baits.pop()

# Collection of all samples from each patient
patientGroups=defaultdict(set)

#Starting output file
output=open(args.output, 'w')
output.write(TCGA_MAF_Ext.header()+"\n")

bc=open(args.baseCounts, 'r')
cin=csv.DictReader(bc,delimiter="\t")
if cin.fieldnames[33]=="Occurence_in_Normals":
    samples=set()
    for si in cin.fieldnames[34:]:
        if si in sampleDb:
            patientGroups[sampleDb[si].Patient_ID].add(si)
            samples.add(si)
else:
    print >>sys.stderr, "unexpected format for mutation file"
    sys.exit(1)

eventDb=dict()

# This is still grabbing information from base counts file
for rec in bunchStream(cin):
    caller=args.caller
    varType=getVarType(rec)
    startPos=rec.Start
    if varType=="INS":
        endPos=str(int(rec.Start)+1)
        refAllele="-"
        altAllele=rec.Alt[1:]
        if caller == "haplotect":
            caller="haplotypecaller"
    elif varType=="DEL":
        startPos=str(int(rec.Start)+1)
        refAllele=rec.Ref[1:]
        endPos=str(int(startPos)+len(refAllele)-1)
        altAllele="-"
        if caller == "haplotect":
            caller="haplotypecaller"
    else:
        endPos=str(int(rec.Start)+len(rec.Ref)-1)
        refAllele=rec.Ref
        altAllele=rec.Alt
        if caller == "haplotect":
            caller="mutect"

    maf=TCGA_MAF_Ext(Hugo_Symbol=rec.Gene,
        Center=CenterTag,
        NCBI_Build=NCBI,
        Chromosome=rec.Chrom,
        Start_Position=startPos,
        End_Position=endPos,
        Strand="+",
        Variant_Classification=rec.VariantClass,
        Variant_Type=getVarType(rec),
        Reference_Allele=refAllele,
        Tumor_Seq_Allele2=altAllele,
        Sequencer=CenterTag,
        )

    event=getEventSig(maf)
    key=(event,sampleDb[rec.Sample].Patient_ID)
    ## If the event is not recorded for that patient ID,
    if key not in eventDb:
        ## save that event (and patient ID)
        eventDb[key]=dict(MAF=maf,mutSamples=set())

        #grab set of samples that belong to patient group
        sampleGroup=set(patientGroups[sampleDb[rec.Sample].Patient_ID])


        for si in patientGroups[sampleDb[rec.Sample].Patient_ID]:
            if si in pairs:
                sampleGroup.add(pairs[si])

        for si in sampleGroup:
            eventDb[key][si]=dict([x.split("=") for x in getattr(rec,si).split(";")])

    eventDb[key]["mutSamples"].add(rec.Sample)


for ei in sorted(eventDb):
    #
    # Get non-normal samples that did not have events for MAF Fill
    # for each mutated sample at that position
    for si in eventDb[ei]["mutSamples"]:
        maf1=copy.copy(eventDb[ei])["MAF"]
        maf1=fillSampleMAFFields(maf1,si,ei,eventDb,pairs,caller,args.species)

        # changes mutaiton status based on if this is unpaired, paired with matched normal, or 
        # paired with pooled nomral
        if si not in pairs:
            maf1.Mutation_Status = "UNPAIRED"
        elif sampleDb[pairs[si]].Class=="Normal":
            maf1.Mutation_Status = "SOMATIC"
        else:
            maf1.Mutation_Status = "SOMATIC_VS_POOL"

        # Grab the TumorSeqAllele1 OR exit with a horrible error
        tum=maf1.Tumor_Sample_Barcode
        norm=maf1.Matched_Norm_Sample_Barcode
        chr=maf1.Chromosome
        varType=maf1.Variant_Type
        start=maf1.Start_Position
        end=maf1.End_Position
        ref=maf1.Reference_Allele
        tumor2=maf1.Tumor_Seq_Allele2

        key="-".join([tum,norm,chr,start,end,varType,ref,tumor2])
        if key not in origMafInfo:
            print>>sys.stderr,"".join(["ERROR: key ",key,"not found in origMafInfo."])
            sys.exit(1)

        maf1.Tumor_Seq_Allele1= origMafInfo["-".join([tum,norm,chr,start,end,varType,ref,tumor2])]

        output.write(str(maf1) + "\n")

    for si in (patientGroups[ei[1]] - eventDb[ei]["mutSamples"]):
        #
        # Only Output Non-normal samples
        #
        if sampleDb[si].Class.upper().find("NORMAL")==-1:
            maf1=copy.copy(eventDb[ei])["MAF"]
            maf1=fillSampleMAFFields(maf1,si,ei,eventDb,pairs,caller, args.species)
            maf1.Mutation_Status = "NONE"
            maf1.Tumor_Seq_Allele1=maf1.Reference_Allele
            output.write(str(maf1)+"\n")


