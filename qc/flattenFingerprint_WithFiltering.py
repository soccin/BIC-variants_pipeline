#!/usr/bin/env python2.7
"""FingerPrint a VCF 

**everything commented out was done so by Caitlin (Nick wrote the script).
The input to this script is now a vcf that has already been filtered, so 
the purpose of this script at this point is strictly to flatten the vcf
format for input to Nick's plotting scripts. Otherwise, with the extra 
filtering included, everything ends up being filtered out. Will readdress
this if it becomes a problem downstream**

"""

import sys
import csv
import lib.vcf as vcf

out = csv.writer(sys.stdout, delimiter='\t')
reader = vcf.VCFReader(sys.stdin)

header = ['CHROM', 'POS', 'REF', 'ALT', 'ID', "FILTER","QUAL","INFO.DP"]+reader.samples 
numSamples=len(reader.samples)
out.writerow(header)    

def flatten(x):
    if type(x) == type([]): 
        x = ','.join(map(str, x))
    return x

for record in reader:
    depth=record.INFO["DP"]
    #print depth
    #print record.ID
    if record.FILTER!=None and "LowQual" in record.FILTER \
       or depth/numSamples<5 \
       or record.ID=="." or record.ID==None:
        #print "doesn't pass filter"
        continue
    if len(record.REF)>1 or len(record.ALT[0])>1 or len(record.ALT)>1:
        #print "indel"
        continue
    genoTypes=[] 
    for sample in record.samples:
        genoTypes.append(sample.data["GT"])
    if "./." in genoTypes or None in genoTypes:
    #    print "not called"
        continue
    if record.FILTER==None:
        filter="PASS"
    else:
        filter=flatten(record.FILTER)
    row = [record.CHROM, record.POS, record.REF, flatten(record.ALT), record.ID, filter]
    row.append(record.QUAL)
    row.append(record.INFO["DP"])
    row.extend(genoTypes)
    out.writerow(row)

