#!/usr/bin/env python2.7
import sys
import PAWK

def floatNA(x):
    if x=="NA":
        return "NA"
    else:
        return float(x)

muTectTypes={"POS":int,"ALT_FREQ":floatNA,"AD_REF":int,
    "AD_ALT":int,"NORM_ALT_FREQ":floatNA,"NORM_AD_REF":int,"NORM_AD_ALT":int}

def castFields(rec,ftypes):
    for fi in ftypes:
        rec._rec[fi]=ftypes[fi](rec._rec[fi])



muTectOkFailures=set(["alt_allele_in_normal","nearby_gap_events",
    "possible_contamination","triallelic_site","DBSNP Site"])


#
# In DMP code these are settable from the command line
# These are the defaults in the code
# which match dmp_impact.v5.conf also

#
# AD_MutectSTDFilter = 5
# DP_MutectSTDFilter = 0
# TNfreqRatio_MutectStdFilter = 5
# VF_MutectSTDFilter = 0.01
#

TNRatio=5
TotalDepth=0
AlleleDepth=3
VariantFreq=0.01


logfp=sys.stderr

for rec in PAWK.PAWK():
    castFields(rec,muTectTypes)

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

    if rec.MUT_KEEP=="REJECT" and rec.MUT_COVERED=="COVERED":
        failureReasons=set(rec.FILTER.split(","))
        if failureReasons <= muTectOkFailures:
            print >>logfp, "OK FAILURES==>", rec.SAMPLE, rec.CHROM, rec.POS, rec.FILTER,
            T_COV=rec.AD_ALT+rec.AD_REF
            print >>logfp, "COV,DP_AD,TNRAF,TNRatio*NRAF",
            print >>logfp, T_COV, rec.AD_ALT, rec.ALT_FREQ, rec.NORM_ALT_FREQ*TNRatio,
            if T_COV>=TotalDepth and rec.AD_ALT>=AlleleDepth \
                    and rec.ALT_FREQ>=TNRatio*rec.NORM_ALT_FREQ \
                    and rec.ALT_FREQ>=VariantFreq:

                rec.MUT_VT="SNP"
                rec.MUT_KEEP="KEEP"
                rec.FILTER="PASS"
                rec.CALLER="".join([rec.CALLER,".Rescue"])

            print >>logfp, "Final.KEEP=", rec.MUT_KEEP

    rec.write()

