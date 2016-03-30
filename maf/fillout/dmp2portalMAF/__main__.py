#!/usr/bin/env python2.7

import sys
import csv
import copy
import argparse
import os.path
from collections import defaultdict
from lib import *
from functs import *

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

sampleDb=populatePatientInfo(args.patient)
CenterTag="MSK-"+sampleDb[sampleDb.keys()[0]].Bait_version

pairs=set()
if args.pairing:
    validFile(args.pairing)
    pairs = populatePairingInfo(args.pairing)

origMafInfo=populateOriginalMafInfo(args.maf)

# Collection of all samples from each patient
patientGroups=generatePatientGroups(sampleDb)

# Populate eventDb to get counts for each sample
eventDb=dict()
bc=open(args.baseCounts, 'r')
cin=csv.DictReader(bc,delimiter="\t")
if cin.fieldnames[33]!="Occurence_in_Normals":
    print >>sys.stderr, "unexpected format for mutation file"
    sys.exit(1)
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
        eventDb[key]=dict(MAF=maf,mutSamples=set(),mutTumors=set())

        #grab set of samples that belong to patient group
        sampleGroup=set(patientGroups[sampleDb[rec.Sample].Patient_ID])

        for si in sampleGroup:
            eventDb[key][si]=dict([x.split("=") for x in getattr(rec,si).split(";")])
    # add t-n pair to mutSamples, so you know which one has the variant called
    eventDb[key]["mutSamples"].add(";".join([rec.Sample, rec.NormalUsed]))
    eventDb[key]["mutTumors"].add(rec.Sample)
    if not args.pairing:
        pairs.add(";".join([rec.Sample, rec.NormalUsed]))

#Starting output file
output=open(args.output, 'w')
output.write(TCGA_MAF_Ext.header()+"\n")

for ei in sorted(eventDb):
    # for each mutated sample at that position
    tnFilledout = set()
    for sini in eventDb[ei]["mutSamples"]:
        tnFilledout.add(sini)
    	# for each mutated sample
    	(si,ni) = sini.split(";")

        if ni in sampleDb and sampleDb[ni].Class != "Normal":
    	    ni="REF"
        maf1=copy.copy(eventDb[ei])["MAF"]
        maf1=fillSampleMAFFields(maf1,si,ei,eventDb,ni,caller,args.species)

        # changes mutaiton status based on if this is unpaired, paired with matched normal, or 
        # paired with pooled nomral
        if sini not in pairs or ni not in patientGroups[sampleDb[si].Patient_ID]:
            maf1.Mutation_Status = "UNPAIRED"
        elif sampleDb[ni].Class=="Normal":
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

    for si in (patientGroups[ei[1]] - eventDb[ei]["mutTumors"]):
        #
        # Only Output Non-normal samples
        #
        print>>sys.stderr, "".join(["This is a fillout for this record: ", str(ei), " sample: ",si])
        #if sampleDb[si].Class.upper().find("NORMAL")==-1:
        maf1=copy.copy(eventDb[ei])["MAF"]
        tnPairMutations = eventDb[ei]["mutSamples"]
        print>>sys.stderr, "".join(["Pairs with mutations: ", str(tnPairMutations)]) 
        tumorPairs = [e for e in pairs if e.startswith(si + ";") or e.endswith(";" + si)]
        print>>sys.stderr, "".join(["tumorPairs: " , str(tumorPairs)])
        possiblePairs = set(tumorPairs).difference(tnPairMutations).difference(tnFilledout)
        print>>sys.stderr, "".join(["Possible pairs: ", str(possiblePairs)])
        for pp in possiblePairs:
            tnFilledout.add(pp)
            si = pp.split(";")[0]
            ni = pp.split(";")[1]
            if ni not in sampleDb or sampleDb[ni].Class != "Normal":
    	        ni="REF"
            maf1=fillSampleMAFFields(maf1,si,ei,eventDb,ni,caller, args.species)
            maf1.Mutation_Status = "NONE"
            maf1.Tumor_Seq_Allele1=maf1.Reference_Allele
            output.write(str(maf1)+"\n")


