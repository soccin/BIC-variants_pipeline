#!/usr/bin/env python2.7

import csv
import sys
import argparse
from itertools import ifilter,dropwhile
from numpy.dual import norm

parser = argparse.ArgumentParser(description="MutectRescue script")
parser.add_argument('-v', '--vcf', action='store', help="VCF file to be rescued", required=True)
parser.add_argument('--txt', action='store', help="txt file associated with vcf file", required=True)
parser.add_argument('-t', '--tumor_id', action='store', help='Tumor id in the vcf', required=True)
parser.add_argument('-n', '--normal_id', action='store', help='Normal id in the vcf', required=True)
parser.add_argument('-o', '--output', action='store', help="Output vcf filename", required=True)
args = parser.parse_args() 

# Nicks trick for making the dictionary a struct?
class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)
       
# Nice clean place to put global values 
class Globals:
  pass
  
muTectOkFailures=set(["alt_allele_in_normal","nearby_gap_events",
    "possible_contamination","triallelic_site","DBSNP Site"])

#
# In DMP code these are settable from the command line
# These are the defaults in the code
# which match dmp_impact.v5.conf also
glbs=Globals()

glbs.TNRatio=float(5)
glbs.TotalDepth=float(0)
glbs.AlleleDepth=float(3)
glbs.VariantFreq=float(0.01)

logfp=sys.stderr

fin=open(args.vcf)
txtin=open(args.txt)
ofh=open(args.output, "wt")

# First just skip the comments in the txt file:
startTxt = dropwhile(lambda L: L.lower().lstrip().startswith('##'), txtin)
tin=csv.DictReader(startTxt,delimiter="\t")

# grab all commented things and print them out first
comments = ifilter(lambda L: L.startswith("##"), fin)
for line in comments:
    ofh.write(line)

# Go back to beginning
fin.seek(0)
start=dropwhile(lambda L: L.lower().lstrip().startswith('##'), fin)
cin=csv.DictReader(start,delimiter="\t")
#This is incase I need to add anything. I don't think I do.
fnames = cin.fieldnames
### NOTE: lineterminator needed so a carriage return isn't added to end of each line in output csv
cout=csv.DictWriter(ofh,fieldnames=fnames,delimiter="\t",lineterminator='\n')
cout.writeheader()

for rec in cin:
    try:
        txtRec=tin.next() 
        #txtRec=Struct(**txtDict)
        #rec=Struct(**recDict)
      
        #1 : verify we are still on the same page here
        if rec["#CHROM"] != txtRec["contig"] or rec["POS"] != txtRec["position"]:
            print >>logfp, "ERROR: records don't match: vcf: " + rec["#CHROM"] + ":" + rec["POS"] + "\ttxt: " +  txtRec["contig"] + ":" + txtRec["position"] + "\n"
            sys.exit(-1)
            
        if "," in rec["ALT"]:
            print >> logfp, "ERROR: This script is not ready for multiple alleles at this point!!\n"
            sys.exit(-1)
            
        # If this is rejected but covered, look more into it
        if rec["FILTER"] == "REJECT" and txtRec["covered"] == "COVERED":
            if set(txtRec["failure_reasons"].split(",")) <=  muTectOkFailures:
                # Grab info field, if AD, good
                formatRow = rec["FORMAT"].split(":")
                if "AD" not in formatRow:
                   # Probably should error or something
                   print >>logfp, "ERROR: There is no AD in the format row: " + str(rec) + "\n"
                   sys.exit(-1)
                adIndex = formatRow.index("AD")
                
                # Now get the tumor and normal fields
                t_field = rec[args.tumor_id].split(":")
                n_field = rec[args.normal_id].split(":")
                tum_ad=t_field[adIndex]
                norm_ad=n_field[adIndex]
                
                tum_alt_ad=0
                tum_ref_ad=0
                norm_alt_ad=0
                norm_ref_ad=0
                
                # If there is no comma in AD, assume it is a ., leave the ad at zero!
                if "," in tum_ad:
                    tum_alt_ad=tum_ad.split(",")[1]
                    tum_ref_ad=tum_ad.split(",")[0]
                if "," in norm_ad:
                    norm_alt_ad=norm_ad.split(",")[1]
                    norm_ref_ad=norm_ad.split(",")[0]
                    
                tum_coverage  = int(tum_alt_ad) + int(tum_ref_ad)
                tum_alt_freq = float(tum_alt_ad)/ tum_coverage
                norm_alt_freq = float(norm_alt_ad)/ (int(norm_alt_ad) + int(norm_ref_ad))
                
                #Log values as Nick did
                print >>logfp, "OK FAILURES==>", args.tumor_id, rec["#CHROM"], rec["POS"], txtRec["failure_reasons"],
                print >>logfp, "COV,DP_AD,TNRAF,TNRatio*NRAF",
                print >>logfp, tum_coverage, tum_alt_ad, tum_alt_freq, (glbs.TNRatio * norm_alt_freq),
                
                # Now see if this passes all the thresholds!
                if tum_coverage >= glbs.TotalDepth \
                    and tum_alt_ad >= glbs.AlleleDepth \
                    and tum_alt_freq >= glbs.VariantFreq \
                    and tum_alt_freq >= glbs.TNRatio * norm_alt_freq: 
                    
                    # Change FILTER to RESCUE
                    rec["FILTER"]="RESCUE"
                    print >>logfp, "RESCUED!"
                else:
                    print >>logfp, "Not rescued."   
        cout.writerow(rec)
                
    except:
        print rec
        print
        print
        raise                
    
    #
    # Re-evaluate rejected but covered events
    # if
    #   failure reason is a subset of the muTectOkFailures
    #   ALT_FREQ >= TNRatio * NORM_ALT_FREQ
    #   T_COV >= TotalDepth
    #   AD_ALT >= AlleleDepth
    #   ALT_FREQ >= VariantFreq
    # then
    #   change REJECT to KEEP

    #if rec.MUT_KEEP=="REJECT" and rec.MUT_COVERED=="COVERED":
    #    failureReasons=set(rec.FILTER.split(","))
    #    if failureReasons <= muTectOkFailures:
    #        print >>logfp, "OK FAILURES==>", rec.SAMPLE, rec.CHROM, rec.POS, rec.FILTER,
    #        T_COV=rec.AD_ALT+rec.AD_REF
    #        print >>logfp, "COV,DP_AD,TNRAF,TNRatio*NRAF",
    #        print >>logfp, T_COV, rec.AD_ALT, rec.ALT_FREQ, rec.NORM_ALT_FREQ*TNRatio,
    ###        if T_COV>=TotalDepth and rec.AD_ALT>=AlleleDepth \
    #                and rec.ALT_FREQ>=TNRatio*rec.NORM_ALT_FREQ \
    #                and rec.ALT_FREQ>=VariantFreq:

    #            rec.MUT_VT="SNP"
    #            rec.MUT_KEEP="KEEP"
    #            rec.FILTER="PASS"
    #            rec.CALLER="".join([rec.CALLER,".Rescue"])#

    #        print >>logfp, "Final.KEEP=", rec.MUT_KEEP




