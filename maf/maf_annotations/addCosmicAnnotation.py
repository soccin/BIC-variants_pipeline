#!/opt/bin/python

from SourceDataUpdater import *
from collections import namedtuple
import re
import os, fnmatch
import csv
import sys
import argparse
import datetime

def fillInBlanks(row, total_len):
    '''Append an "N/A" to any blank fields at the end of a row'''
    for x in range(0,total_len-len(row)):
        row.append("N/A")
    return row  

def validateColNames(allColNames):
    '''Names in namedtuples can only contain alphanumeric characters and underscores;
    Removes any invalid characters, and replace whitespace with underscore'''
    validHeader = []
    for colName in allColNames:
        colName = re.sub('[^0-9a-zA-Z]+','_',colName)
        colName = colName.lower()
        validHeader.append(colName)
    return validHeader        

def separateHeader(iterable):
    "Separates the header from the rest of the table; returns each as separate entity"
    it = iter(iterable)
    header = '#'
    while header[0][0] == '#':
        header = it.next()
    #return header, remainder of table
    return header, it

def get_table(header,body, ntuple=None):
    '''Takes as an argument a csv reader object of all data in
    file. Create a namedtuple based on header row;
    '''
    #header, body = separateHeader(header_with_body)
    #header = validateColNames(header)
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

def locateCosmicFile(pattern, root=os.curdir):
    '''
    Find cosmic file by pattern under the specified directory; If more than on 
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        if len(fnmatch.filter(files,pattern)) == 1:
            return os.path.join(path, fnmatch.filter(files,pattern)[0])
        else:
            print>>sys.stderr, "ERROR: more than one cosmic file found; please either rerun using the '-f' option to specify which file to use or delete all old COSMIC files."
            sys.exit(-1)
    return

def getRefAndAlt(change):
    '''
    Parse the mutation, which is in standard syntax (e.g., c.194A>G), to isolate the reference and alternate 
    '''
    sub = []
    sub = re.split(">",change)
    ref = alt = None
    if len(sub) != 2:
        sub = re.split("del",change)
        if len(sub) != 2:
            sub = re.split("ins",change)
            if len(sub) != 2:
                return ref, alt
            else: ## mut is an insertion
                ref = "?"
                alt = sub[1]
        else: ## mut is a deletion
            ref = sub[1]
            alt = "?"
    else: ## mut is a substitution
        ref = sub[0]
        alt = sub[1]

    return ref,alt

def indexCosmicRecords(cosmicFile):
    '''
    Parse COSMIC text file and index according to hg19 genomic position, reference and alternate
    '''
    file = open(cosmicFile, 'rU')
    cosmicData = csv.reader(file,delimiter='\t')

    indexedCosmic = {}
    header, body = separateHeader(cosmicData)
    header = validateColNames(header)
 
    for rec in get_table(header,body):
        ## ignore records with unknown ref/alt
        if not rec.mutation_cds == "c.?" and 'N/A' not in rec.mutation_cds:   
            ## isolate the mutation from its CDS coordinates
            change = re.sub("c\.\d+((_\d+)*((\+|-)\d+)*)*","",rec.mutation_cds)
            ## parse mutation
            ref, alt = getRefAndAlt(change)
    
            if ref and alt and rec.mutation_grch37_genome_position and len(rec.mutation_grch37_genome_position)>0:
                rec_key = "_".join([rec.mutation_grch37_genome_position,ref,alt])
                indexedCosmic[rec_key] = rec

    return indexedCosmic

def annotateMAF(inputMAF, annotatedMAF, mafHeader, cosmicHeader, detailed, indexedCosmic):
    """
    read input MAF file, match each record to associated COSMIC record
    """
    matched_records = total_records = 0

    with open(inputMAF, 'rU') as maf:
        
        with open(annotatedMAF,'w') as out:

            mafData = csv.reader(maf,delimiter="\t")
            writer = csv.writer(out,delimiter="\t")
            if detailed:
                header = mafHeader.replace("\n","").split("\t") + cosmicHeader
            else:
                header = mafHeader.replace("\n","").split("\t") + ["COS_primary_site","COS_primary_histology","COS_mutation_id","COS_mutation_satic_status"]
            writer.writerow(header)

            header2, body = separateHeader(mafData)
            #header2 = validateColNames(header)

            for rec in get_table(header2,body):
                #print rec
                total_records += 1
                chr = rec.Chromosome.replace("chr","")
                genomic_position = chr + ":" + rec.Start_Position + "-" + rec.End_Position
                ref = rec.Reference_Allele
                alt1 = rec.Tumor_Seq_Allele1
                alt2 = rec.Tumor_Seq_Allele2

                if rec.Variant_Type == "DEL":
                    alt1 = alt2 = "?"
                elif rec.Variant_Type == "INS":
                    ref = "?"

                rec_key1 = "_".join([genomic_position, ref, alt1])
                rec_key2 = "_".join([genomic_position, ref, alt2])

                matched_record = None

                if indexedCosmic.has_key(rec_key1):
                    matched_record = indexedCosmic[rec_key1]
                    print "Found COSMIC annotations for SNP: ",rec_key1
                elif indexedCosmic.has_key(rec_key2):
                    matched_record = indexedCosmic[rec_key2]
                    print "Found COSMIC annotations for SNP: ",rec_key2
                annotated_record = rec

                if matched_record:
                    print matched_record
                    matched_records += 1
                    if detailed:
                        ## add ALL COSMIC fields to MAF
                        annotated_record += matched_record 
                    else:
                        ## add only 'standard' fields to MAF
                        annotated_record += tuple([matched_record.primary_site,matched_record.primary_histology,matched_record.mutation_id,matched_record.mutation_somatic_status])
                else:
                    ## add blank fields so that if other annotations are to be added, they
                    ## will be printed in the correct columns
                    if detailed:
                        annotated_record += tuple([""]*len(cosmicHeader))
                    else:
                        annotated_record += tuple([""]*4)
                writer.writerow(annotated_record)

    print "Found COSMIC annotations for %d out of %d records in MAF" % (matched_records,total_records)
    return

def getParsedArgs():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Add COSMIC annotations to MAF")

    ## required arguments
    req = parser.add_argument_group('Input')
    req.add_argument('-i', dest='inputMAF', nargs=1, required=True, help='MAF to be annotated')
    req.add_argument('-o', dest='outMAF', nargs=1, required=True, help='output MAF')
    req.add_argument('-l', '--local', dest='local_dir', nargs=1, help='Directory to which COSMIC should be downloaded; must use this option if not specifying a local COSMIC file')
    req.add_argument('-f', '--cosmic', dest='cosmicFile', nargs=1, help='Local COSMIC file; must use this option if user does NOT wish to download the latest version of COSMIC')
    req.add_argument('-d', '--detailed', dest='detailed', default=False, action='store_true', help='Include ALL COSMIC fields in MAF; by default, only "standard" fields are included.')


    return parser

def downloadCosmic(local_dir):

    daysAgo = 1000 ## default; the COSMIC ftp dir only contains the most current version, so 1000 should always be sufficient
    alert_dir = os.getcwd() ## default; this isn't actually used in this case

    ## download latest COSMIC file 
    updater = SourceDataUpdater(daysAgo, local_dir, alert_dir)
    updater.update_cosmic()

    return

def main():
    print "\n"
    start = datetime.datetime.now()
    args = getParsedArgs().parse_args()
    inputMAF = args.inputMAF[0]
    outMAF = args.outMAF[0]
    detailed = args.detailed                 
    
    if not args.cosmicFile[0]:
        if not args.local_dir[0]:
            print>>sys.stderr,"Please specify either a COSMIC file to use or a local directory to which the latest COSMIC file should be downloaded."
        else:
            local_dir = args.local_dir[0]
            print "Downloading COMIC file..."
            downloadCosmic(local_dir)
        
            ## we don't necessarily know the exact name of the downloaded file, so look for it in the specified download dir
            ## an error will be reported if more than one COSMIC file is found
            print "Locating local COSMIC file..."
            cosmicFile = locateCosmicFile("*CosmicMutant*",local_dir)
        
    else:
        cosmicFile = args.cosmicFile[0]
    
    print "Indexing COSMIC file by genomic location and mutation..."
    indexedCosmic = indexCosmicRecords(cosmicFile)

    ## temporary hack to get headers; got errors when trying to do it using getTable generator
    mafHeader = cosmicHeader = '#'
    with open(inputMAF,'rU') as maf:
        while mafHeader[0][0] == '#':
            mafHeader = maf.next()

    with open(cosmicFile, 'rU') as cosmic:
        while cosmicHeader[0][0] == '#':
            cosmicHeader = cosmic.next()
        cosmicHeader = validateColNames(cosmicHeader.strip().split("\t"))
        cosmicHeader = ["COS_"+h for h in cosmicHeader]
        print cosmicHeader

    print "Adding COSMIC annotations to MAF..."
    ## add COSMIC annoations to MAF
    annotateMAF(inputMAF, outMAF, mafHeader, cosmicHeader, detailed, indexedCosmic)
    end = datetime.datetime.now()

    print "Script began:",start.strftime("%Y-%m-%d %H:%M")
    print "Script finished:",end.strftime("%Y-%m-%d %H:%M")
    print "\n"   
 
if __name__ == '__main__':
    sys.exit(main())
