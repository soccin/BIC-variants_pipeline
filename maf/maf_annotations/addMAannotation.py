#! /opt/common/CentOS_6/python/python-2.7.6/bin/python

import urllib2
import time
import re
import datetime
import argparse
import sys
import csv
from collections import namedtuple
from pymongo import MongoClient

def getArgs():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in", dest="input", nargs=1, required=True, type=str, help="input MAF file")
    parser.add_argument("-o", "--out", dest="out", nargs=1, required=True, type=str, help="output MAF file (annotated)")
    #parser.add_argument("-v", "--version", dest="version", nargs=1, required=True, type=str, help="Mutation Assessor version with which to annotate MAF")
    parser.add_argument("-d", "--detailed", dest="detailed", default=False, action="store_true", help="include ALL Mutation Assessor fields in output; by default only 'standard' fields are included in output")

    args = parser.parse_args()

    return args

def resultsComplete(results):
    if '[sent]' in results.values() or '[sent]' in results.keys():
        return 0
    else:
        return 1

def createMongoRecord(results,maVersion):
    results = results.split('\n')
    header = results[0].split('\\t')
    val_header = validateColNames(header)
    info = results[1].split('\\t')
    results = dict(zip(val_header, info))
    build,chr,start,ref,alt = results['mutation'].split(',')
    results['date'] = datetime.datetime.utcnow()
    results['source'] = 'mutation_assessor'
    results['version'] = maVersion
    results['build'] = build
    results['chr'] = "chr"+chr
    results['start'] = int(start)
    results['end'] = int(start)
    results['ref'] = ref
    results['alt'] = alt

    return results
                          
def fillInBlanks(row, total_len):
    '''
    Append an "N/A" to any blank fields at the end of a row
    '''
    for x in range(0,total_len-len(row)):
        row.append("N/A")
    return row

def validateColNames(allColNames):
    '''
    Names in namedtuples can only contain alphanumeric characters and underscores;
    Removes any invalid characters, and replace whitespace with underscore
    '''
    validHeader = []
    for colName in allColNames:
        colName = re.sub('[^0-9a-zA-Z]+','_',colName)
        colName = colName.lower()
        validHeader.append(colName)
    return validHeader

def separateHeader(iterable):
    '''
    Separate the header from the rest of the table; return each as separate entity
    '''

    it = iter(iterable)
    #return header, remainder of table
    header = '#'
    while header[0][0] == '#':
        header = it.next()
    return header, it

def get_table(header, body, ntuple=None):
    '''
    Takes as an argument a csv reader object of all data in
    flattened run log. Create a namedtuple based on header row;
    '''
    ntuple = ntuple or namedtuple('MutationRecord', header)
    for row in body:
        if len(row)<len(header):
            row = fillInBlanks(row, len(header))
        elif len(row)>len(header):
            for x in range(len(header),len(row)):
                if row[x] == '':
                    del row[x]
                else:
                    print "ERROR: Row contains extra fields."
        yield ntuple(*row)

def annotate():
    '''
    Query Mutation Assessor for every somatic mutation in a TCGA MAF file and
    append the results to the record.
    '''

    ## set up connection to db cache
    ## if connection fails, cache will just 
    ## not be used. instead, all queries will
    ## go to mutationassessor.org
    mongoserver = "eos"
    db = None
    try:
        client = MongoClient(mongoserver)
        db = client.genomic_events
        annots = db.annotations
    except ConnectionFailure, e:
        print>>sys.stderr,"Could not connect to MongoDB!"

    standard_fields = ['Func. Impact', 'FI score', 'VC score', 'VS score', 'cancer types in COSMIC@position', "gene's known role in cancer", 'regions@position', 'domain@position', 'domains']
    all_fields = ['Mutation', 'RG variant', 'RG var.type', 'Mapping issue', 'User data', 'AA variant', 'Gene', 'MSA', 'PDB', 'Func. Impact', 'FI score', 'VC score', 'VS score', 'Uniprot', 'Refseq', 'gaps in MSA', 'MSA height', 'Codon start position', 'Uniprot position', 'Uniprot residue', 'Refseq position', 'Refseq residue', 'Func. region', 'Protein bind.site', 'DNA/RNA bind.site', 'small.mol bind.site', 'N.Cosmic', 'N.SNPs', 'mutations in COSMIC@position', 'cancer types in COSMIC@position', 'SNPs@position', "gene's known role in cancer", 'regions@position', 'domain@position', 'domains']
 
    query_prefix = 'http://mutationassessor.org/?cm=var&var=hg19'
    query_suffix = '&fts=all&frm=txt'

    args = getArgs()

    ### ultimately, mutation assessor version will be a required
    ### argument to this script. Until then, set it to '2', since
    ### that is the only version available right now
    #maVersion = args.version[0]
    maVersion = '2'

    build = 'hg19'

    if not args.detailed: 
        fields = standard_fields
    else:
        fields = all_fields

    ## be sure to "validate" ONLY the field names that are to 
    ## be added to existing MAF. This code assumes headers
    ## in existing MAF are valid and should NOT be altered.
    ## Specifically, this is to avoid altering TCGA-specific headers
    fields = validateColNames(fields) 
    col_headers = ["MA_"+f for f in fields]

    #print fields

    with open(args.out[0], 'w') as out:
        with open(args.input[0], 'r') as f:

            mafData = csv.reader(f, delimiter="\t")
            writer = csv.writer(out, delimiter="\t")

            header,body = separateHeader(mafData)
            writer.writerow(list(header)+col_headers)
            numCols = len(list(header)+fields)

            query_counter = 0

            for rec in get_table(validateColNames(header),body):
                chr = rec.chromosome.replace("chr","")
                alt = ''
    
                if not rec.tumor_seq_allele1 == rec.reference_allele:
                    alt = rec.tumor_seq_allele1
                elif not rec.tumor_seq_allele2 == rec.reference_allele:
                    alt = rec.tumor_seq_allele2

                ## mutation assessor does not support indels, as far as I can tell, 
                ## so only generate a query for single nucleotide substitutions
                if len(alt)==1 and len(rec.reference_allele)==1: 
                    results = None 
                    if db:       
                        ## check cache before attemptying to run web query
                        results = annots.find_one({'build':build,\
                                               'source':'mutation_assessor',\
                                               'chr':'chr'+chr,\
                                               'start':int(rec.start_position),\
                                               'ref':rec.reference_allele,\
                                               'alt':alt})#,\
                                               #'version':maVersion})
                    if not results:
                        ## run query & parse results if not already cached
                        query = ",".join([query_prefix,chr,rec.start_position,rec.reference_allele,alt]) + query_suffix
                        results = createMongoRecord(urllib2.urlopen(query, None, 60).read(),maVersion)
                        if resultsComplete(results): ## check that query ran to completion before adding results to cache
                            annots.insert(results)
                            query_counter += 1

                    for field in fields:
                        rec += tuple([results[field]])
                    if not len(rec) == numCols:
                        print>>sys.stderr,"WARNING: Number of elements in results from query %s != number of column headers" %query
                else:
                    ## add blank fields to account for snps that are not queried
                    rec += tuple(['']*len(fields))
                
                writer.writerow(rec)

                if query_counter == 5:
                    time.sleep(0.25)
                    query_counter = 0
 

if __name__ == '__main__':
    annotate()

