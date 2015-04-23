#!/opt/bin/python 

import sys
import scipy.stats
import vcf
import operator
import argparse

def get_args(sysargs):
    parser = argparse.ArgumentParser()
    parser.add_argument("-vcf",dest='infile',nargs=1,required=True,help="VCF to use to evaluate sample pairs")
    parser.add_argument("-out",dest='outfile',nargs=1,required=True,help="output file name")
    parser.add_argument("-pairs",dest='pairs',nargs=1,required=False,help="tab-delimited file containing supposed sample pairs")
    args = parser.parse_args()

    return args

def all_pairs(samples):
    '''
    Creates a list containing all possible pairs of samples. 
    Each element in the list is a 2-element list where each
    element is one sample. 
    '''
    all_pairs = []
    for s in range(len(samples)):
        for x in range(s+1,len(samples)):
            all_pairs.append([samples[s],samples[x]])

    return all_pairs  

def evaluate_pairs(infile):
    """
    For each possible pair of samples in a VCF file, counts the total number
    of records for which there is a call for at least one sample in the pair, 
    and the number of assays in which the GT is the same for both samples in the pair

    Returns a dictionary where keys are all possible pairs of samples in the form
    "Sample1+Sample2" and values are a two element list where the first element is 
    the total number of records counted and the second is the number of matches counted.
    """

    pair_info = {}

    with open(infile, 'rb') as vcf_file:
        vcf_reader = vcf.VCFReader(vcf_file)
        samples = vcf_reader.samples
        # get all possible pairs of samples
        app = all_pairs(samples)

        for pair in app:
            pair_info["+".join(pair)] = [0,0] #[total_num_assays,total_num_matches]        

        for record in vcf_reader:
            for pair in pair_info:
                s1,s2 = pair.split("+")
                s1gt = record.genotype(s1)['GT']
                s2gt = record.genotype(s2)['GT']

                ## a record is added to total number of assays if there is a call for
                ## at least one sample in the pair
                if not (s1gt == None and s2gt == None): ## a "./." string is converted to None by vcfreader
                    pair_info[pair][0] += 1 ## add one to total number of assays

                    ## a record is added to number of matches if there are calls for both 
                    ## samples they have the same genotype 
                    if s1gt == s2gt:
                        pair_info[pair][1] += 1 ## add one to total number of matches

    return pair_info

def write_pairing_stats(pair_info,outfile,supposed_pairs=None):
    """
    For each pair of samples in a dictionary, calculates a p-value and Bonferroni 
    corrected p-value based on the total number of assays counted and number of
    matches counted. 

    Prints a tab-delimited outfile containing all counts and calculations.
    If a pairing file was given, print only those records whose p-values fall
    within the range of those of the supposed pairs. 
    """
    sp = []
    min_bonp_sp = None
    max_bonp_sp = None
    if supposed_pairs:
        with open(supposed_pairs,'r') as s:
            for line in s:
                s1,s2 = line.strip().split()
                sp.append(s1+"+"+s2)
                sp.append(s2+"+"+s1) 

    res = []

    num_pairs = len(pair_info.keys())
    for pair, vals in pair_info.iteritems():
        s1,s2 = pair.split("+")
        p = scipy.stats.binom_test(vals[1],vals[0],p=0.375)
        bon_p = num_pairs*p
        res.append([s1,s2,vals[1],vals[0],p,bon_p])

        if len(sp)>0 and pair in sp:
            if not min_bonp_sp and not max_bonp_sp:
                min_bonp_sp = max_bonp_sp = bon_p
            elif bon_p < min_bonp_sp:
                min_bonp_sp = bon_p
            elif bon_p > max_bonp_sp:
                max_bonp_sp = bon_p

    if max_bonp_sp:
        print>>sys.stdout,max_bonp_sp
  
    with open(outfile, 'w') as out:
        print>>out, "\t".join(['Sample','Sample','Matches','Total Assays','p-value','Bonferroni corrected p-value'])
        for p in sorted(res, key=operator.itemgetter(5)):
            ## if pairing file was given, skip records whose bonp falls
            ## outside of the range of the bonp values of the supposed pairs
            if len(sp) > 0 and (p[5]<min_bonp_sp or p[5]>max_bonp_sp):
                continue
            print>>out, "\t".join([ str(x) for x in p])
            #outstr = "\t".join([str(x) for x in p[0:4] ])
            #outstr += "\t" + format(p[4], '.50f')
            #outstr += "\t" + format(p[5], '.50f')
            #print>>out,outstr

    return

def find_pairs(args):
    """
    Calculates and outputs the likelihood of each possible pair of samples in a VCF file
    being a real tumor/normal pair. 
    """
    pair_info = evaluate_pairs(args.infile[0])
    if args.pairs:
        write_pairing_stats(pair_info,args.outfile[0],args.pairs[0])
    else:
        write_pairing_stats(pair_info,args.outfile[0])


    return


if __name__ == '__main__':
    args = get_args(sys.argv)
    find_pairs(args)
    sys.exit(0)


