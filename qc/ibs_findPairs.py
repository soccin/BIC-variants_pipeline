#!/opt/bin/python 

import sys
import scipy.stats
import numpy
import vcf
import operator

infile = sys.argv[1]
outfile = sys.argv[2]

def all_pairs(samples):
    '''
    Create a list containing all possible pairs of samples. 
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
            pair_info["+".join(pair)] = [0,[]] #[total_num_assays,score]        

        for record in vcf_reader:
            for pair in pair_info:
                s1,s2 = pair.split("+")
            
                ## a record is added to total number of assays if there is a call for
                ## at least one sample in the pair
                if not record.genotype(s1)['GT'] == "./." or not record.genotype(s2)['GT'] == "./.":
                    pair_info[pair][0] += 1 ## add one to total number of assays

                ## give record a score of 2 if genotypes are exactly equal,
                ## give record a score of 1 if only one allele is different
                ## give record a score of 0 if genotypes are AA and BB
                if record.genotype(s1).called and record.genotype(s2).called:
                    s1_1,s1_2=record.genotype(s1)['GT'].split("/")
                    s2_1,s2_2=record.genotype(s2)['GT'].split("/")

                    if len(set([s1_1,s1_2,s2_1,s2_2])) == 1: ## all are equal
                        pair_info[pair][1].append(2)
                    elif (s1_1==s1_2) and (s2_1==s2_2): ## AA and BB
                        pair_info[pair][1].append(0)
                    elif ((s1_1==s1_2) and not (s2_1==s2_2)) or ((s2_1==s2_2) and not (s1_1==s1_2)): ## AA AB or BB AB
                        pair_info[pair][1].append(1) 

    return pair_info

def write_pairing_stats(pair_info,outfile):
    """
    For each pair of samples in a dictionary, calculates a p-value and Bonferroni 
    corrected p-value based on the total number of assays counted and number of
    matches counted. 

    Prints a tab-delimited outfile containing all counts and calculations.
    """
    
    res = []
    for pair, vals in pair_info.iteritems():
        res.append(pair.split("+")+[vals[0],numpy.mean(vals[1]),numpy.var(vals[1]),numpy.std(vals[1])])

    with open(outfile, 'w') as out:
        print>>out, "\t".join(['Sample','Sample','Total Assays','Mean IBS score','IBS variance','Stdev'])
        num_pairs = len(pair_info.keys())
        for p in sorted(res, key=operator.itemgetter(4)):
            print>>out, "\t".join([ str(x) for x in p])

    return

def find_pairs(args):
    """
    Calculates and outputs the likelihood of each possible pair of samples in a VCF file
    being a real tumor/normal pair. 
    """
    if not len(args) == 2:
        print >> sys.stderr, "Usage: %prog input_vcf out_txt"
        sys.exit(-1)

    infile,outfile = args
    pair_info = evaluate_pairs(infile)
    write_pairing_stats(pair_info,outfile)

    return


if __name__ == '__main__':
    find_pairs(sys.argv[1:])
    sys.exit(0)


