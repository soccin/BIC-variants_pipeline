#!/opt/bin/python

import sys
from Filter import SNPFilter
import time

"""
Filter a VCF file to get a subset of calls to be used to identify sample pairs.

Filtering is based on the following criteria: 
    (1) minimum percentage of samples with a call at a given position (default: 95%)
    (2) minimum and maximum alt allele frequency (default: 49%-51%)
    (3) maximum type percentage (max het/hom calls) (default: 60%)
    (4) minimum qual (default: 500)

If number of records >Y (default 250) after filtering on these criteria, every Nth
record is taken in order to reach a number between X and Y (default: 40 and 250).

If number of records <X (default 40) after filtering on these criteria, restrictions
are loosened until number of records is between X and Y (default: 40 and 250). 

Usage: %prog input_vcf filtered_vcf_name
"""

def usage():
    print>>sys.stderr,"%prog input_vcf filtered_vcf [min_recs_passing] [max_recs_passing]"
    sys.exit(-1)

## defaults
min_recs_passing = 40
max_recs_passing = 250
## start with most stringent conditions
min_pc = 0.95
min_af = 0.49
max_af = 0.51
max_tp = 0.60
min_qual = 500
pick_every = 1

## for clean output
suff = {1:'st',2:'nd',3:'rd',4:'th',5:'th',6:'th',7:'th',8:'th',9:'th',0:'th'}


if len(sys.argv) == 5:
    min_recs_passing = sys.argv[3]
    max_recs_passing = sys.argv[4]

elif not len(sys.argv) == 3:
    usage()

vcf_file = sys.argv[1]
filtered_vcf = sys.argv[2]   

filter = SNPFilter(vcf_file,filtered_vcf)
start = time.time()

## bin records by percent called, alt allele freq, type percent and qual
new_start = time.time()
print>>sys.stderr, ""
print>>sys.stderr, "Binning records in",vcf_file,"based on PC,AF,TP and QUAL..."
filter.bin_records()
print>>sys.stderr, "Done."
print>>sys.stderr, "Indexing took %d seconds" %(time.time() - new_start)
print>>sys.stderr, ""

## get list of record IDs that pass filtering criteria 
num_recs_passing = 0
ids = filter.passing_records(min_pc=min_pc,min_af=min_af,max_af=max_af,max_tp=max_tp,min_qual=min_qual)
num_recs_passing = len(ids)/pick_every

## make sure not to get caught in infinite loop, as
## method of adjusting parameters is not exactly stable
max_iter = 20
num_iter = 0

## adjust until length of list falls in desired range
while not num_recs_passing in range(min_recs_passing,max_recs_passing+1) and num_iter <= max_iter: 
    if num_recs_passing < min_recs_passing:
        print>>sys.stderr,"Loosening restrictions..."
        min_af -= 0.05
        max_af += 0.05
        max_tp += 0.05
    elif num_recs_passing > max_recs_passing:
        pick_every = (num_recs_passing/max_recs_passing)+1 ## the '+1' is to make sure pick_every !=1 at this point
        print>>sys.stderr,"Picking every %d%s record" %(pick_every,suff[pick_every%10])
    ids = filter.passing_records(min_pc=min_pc,min_af=min_af,max_af=max_af,max_tp=max_tp,min_qual=min_qual)
    num_recs_passing = len(ids)/pick_every
    print>>sys.stderr, "Num recs:",num_recs_passing
    num_iter += 1

## exit with error if unable to get within range after 20 tries
if not num_recs_passing in range(min_recs_passing,max_recs_passing+1):
    print>>sys.stderr, "ERROR: After %d tries, unable to get number of passing records within range: %d to %d.\nPlease review manually." %(max_iter,min_recs_passing,max_recs_passing)
    sys.exit(1)

## retrieve records in final list and write new VCF
print>>sys.stderr, "Retrieving %d records and writing filtered VCF..." %num_recs_passing

new_start = time.time()
filter.write_vcf_records(ids=ids,pick_every=pick_every)

## log parameters used in final filtred VCF
print>>sys.stderr, "Done."
print>>sys.stderr, "Writing filtered VCF took %d seconds" %(time.time()-new_start)
print>>sys.stderr, ""
print>>sys.stderr, "Total time: %d seconds" %(time.time()-start)
print>>sys.stderr, ""
print>>sys.stderr, ""
print>>sys.stderr, "Kept every %d%s record passing the following criteria:" %(pick_every,suff[pick_every%10])
print>>sys.stderr, "Minimum percentage of samples with call:", min_pc
print>>sys.stderr, "Minimum alternate allele frequency:", min_af
print>>sys.stderr, "Maximum alternate allele frequency:", max_af
print>>sys.stderr, "Maximum percentage of samples with same type of call (hom/het):", max_tp
print>>sys.stderr, "Minimum quality:", min_qual
print>>sys.stderr, ""


