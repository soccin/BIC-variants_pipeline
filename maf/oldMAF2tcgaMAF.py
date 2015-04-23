#!/usr/bin/env python2.7

import TCGA_MAF
import sys
import csv

class Bunch:
    def __init__(self,dictRec):
        self.__dict__.update(dictRec)
    def __str__(self):
        return str(self.__dict__)

def bunchStream(cin):
    for rec in cin:
        yield Bunch(rec)

def getALT(s):
    if(args.verbose and s.GT == "0/0"):
      return s.ALT
    alts=s.ALT.split(";")
    altGT=list(set([int(x) for x in s.GT.split("/") if x !="0"]))
    return alts[altGT[0]-1]

def getVarType(s):
    alt=getALT(s)
    if len(s.REF)==len(alt):
        if len(s.REF)==1:
            return "SNP"
        elif len(s.REF)==2:
            return "DNP"
        elif len(s.REF)==3:
            return "TNP"
        else:
            return "ONP"
    elif len(s.REF)>len(alt):
        return "DEL"
    else:
        return "INS"

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("GENOME",help="Genome build, must be specified")
parser.add_argument("maf0", help="Old maf file")
parser.add_argument("maf1", help="New maf file")
parser.add_argument('-v','--verbose',action='store_true',help='If specified, create a verbose maf that includes vcf entries with no reads')
args=parser.parse_args()

NEWFLDS="FILTER QUAL GT GQ ALT_FREQ NORM_GT NORM_GQ NORM_ALT_FREQ".split()
with open(args.maf0) as input:
  header = input.readline().strip().split()
  if "t_ref_count" in header:
    NEWFLDS = "FILTER QUAL GT GQ ALT_FREQ t_ref_count t_alt_count NORM_GT NORM_GQ NORM_ALT_FREQ n_ref_count n_alt_count".split()
  elif("HC_SNPEFF_EFFECT" in header or "HC_SNPEFF_FUNCTIONAL_CLASS" in header or "HC_SNPEFF_GENE_NAME"):
    NEWFLDS="FILTER QUAL GT GQ ALT_FREQ NORM_GT NORM_GQ NORM_ALT_FREQ HC_SNPEFF_AMINO_ACID_CHANGE HC_SNPEFF_CODON_CHANGE HC_SNPEFF_EFFECT HC_SNPEFF_EXON_ID HC_SNPEFF_FUNCTIONAL_CLASS HC_SNPEFF_GENE_BIOTYPE HC_SNPEFF_GENE_NAME HC_SNPEFF_IMPACT HC_SNPEFF_TRANSCRIPT_ID".split()


class TCGA_MAF_Ext(TCGA_MAF.TCGA_MAF):
    pass
TCGA_MAF_Ext.addFields(NEWFLDS)
with open(args.maf1, 'w') as output:
  with open(args.maf0, 'rb') as input:
    output.write(TCGA_MAF_Ext.header() + "\n")
    dreader = csv.DictReader(input, delimiter="\t")
    for rec in bunchStream(dreader):
      if not args.verbose and ("0/0" in rec.GT or "./." in rec.GT):
        sys.stderr.write("Skipping:" +  str(rec) +  ". There is 0 coverage\n")
        continue
      matchedNormSampleBarcode=rec.NORM_SAMPLE if hasattr(rec, "NORM_SAMPLE") else "REF."+args.GENOME
      maf=TCGA_MAF_Ext(
        Chromosome=rec.CHROM,
        Start_Position=rec.POS,
        End_Position=int(rec.POS)+len(rec.REF)-1,
        Reference_Allele=rec.REF,
        Tumor_Seq_Allele1=getALT(rec),
        dbSNP_RS=rec.ID,
        Tumor_Sample_Barcode=rec.SAMPLE,
        Matched_Norm_Sample_Barcode=matchedNormSampleBarcode,
        Variant_Type=getVarType(rec),
        Hugo_Symbol=rec.GENE,
        t_ref_count=rec.AD_REF,
        t_alt_count=rec.AD_ALT,
        Caller=rec.CALLER
      )
      if hasattr(rec, "NORM_AD_REF"):
        maf.n_ref_count=rec.NORM_AD_REF
        maf.n_alt_count=rec.NORM_AD_ALT
      maf.NCBI_Build=args.GENOME
      for fld in NEWFLDS:
         attr=getattr(rec,fld) if hasattr(rec,fld) else ""
         setattr(maf,fld,attr)
      output.write(str(maf) + "\n")
