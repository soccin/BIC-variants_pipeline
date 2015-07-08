#!/usr/bin/python

import sys
import re
import os
import fnmatch
from collections import OrderedDict

#####################################
## Find all count files matching pattern entered by user (there should be one count file
## per sample) and create a matrix with one row for each gene and one column for each sample
#####################################

#####################################
## Usage: /usr/bin/python/ makeCountMatrix.py rootDir patternToSearch outputFileName
## Example: /usr/bin/python makeCountMatrix.py /ifs/res/liang/RNASeq/Proj2983_MassagueJ .htseq_count Proj2983_MassagueJ_htseq.count_allSamples.txt
#####################################

def usage():
    print "/usr/bin/python rnaseq_count_matrix.py rootDir patternToSearch outputFileName"
    return


## create a list of full paths to files that match pattern entered
def findFiles(rootDir,pattern):
    """
    create and return a list of full paths to files that match pattern entered
    """
    filepaths = []
    for path, dirs, files in os.walk(os.path.abspath(rootDir)):
        if fnmatch.filter(files, pattern):
            for file in fnmatch.filter(files, pattern):
                filepaths.append(os.path.join(path,file))
        else:
            if "Proj" in path.split("/")[-1] and "-" in path.split("/")[-1]:
                print>>sys.stderr, "WARNING: No files matching pattern %s found in %s" %(pattern, path)

    return filepaths


def printMatrix(matrix,outFile):
    """
    """

    header = "\t".join(["Sample","R1_ProcessedReads","R1_TrimmedReads","R1_PercentTrimmedReads",\
                                 "R2_ProcessedReads","R2_TrimmedReads","R2_PercentTrimmedReads"])
    with open(outFile,'w') as out:
        print>>out,header
        for id,info in matrix.items():
            r1info = info["R1"]
            r2info = info["R2"]
            r1pct = format((r1info['trimmed']/r1info['processed'])*100,'.2f')
            try:
                r2pct = format((r2info['trimmed']/r2info['processed'])*100,'.2f')
            except ZeroDivisionError:
                r2pct = "0.0"
            r = [id]+[r1info['processed'],r1info['trimmed']]+[r1pct]+[r2info['processed'],r2info['trimmed']]+[r2pct]
            print>>out,"\t".join([str(x) for x in r])
    return

def makeMatrix(args):
    """
    Find files to parse, create one matrix of all counts and print 
    matrix to file
    """

    if len(args) == 3:
        rootDir,filePattern,outFile = args
    else:
        usage()
        sys.exit(1)

    ## store all values in an ordered dict, keyed by sample
    matrix = OrderedDict()

    filePattern = '*'+filePattern
    ## find all cutadapt stats files using pattern 
    files = findFiles(rootDir,filePattern)

    if files:
        print>>sys.stderr, "\nCombining the following files:\n"
        for file in files:
            print>>sys.stderr, file
            if os.stat(file)[6]==0: ## file is empty
                print>>sys.stderr, "WARNING: This file is empty!"
            else:
                samp = file[file.find('intFiles')+9:].split("/")[0]
                fName = file.split("/")[-1]
                rNum = fName[:fName.find(".fastq")].split("_")[-2]
                if not samp in matrix:
                    matrix[samp] = {'R1':{'processed':0.0,'trimmed':0.0}, \
                                    'R2':{'processed':0.0,'trimmed':0.0} \
                                   }
                with open(file,'r') as fl:
                    for line in fl:
                        if 'Processed reads:' in line:
                            matrix[samp][rNum]['processed'] += float(line.strip().split()[-1])
                        elif 'Trimmed reads:' in line: 
                            matrix[samp][rNum]['trimmed'] += float(line.strip().split()[-2])
                            break                            

        printMatrix(matrix,outFile)            

    else:
        print>>sys.stderr, "\nNo files found matching pattern entered. Exiting.\n"
        sys.exit(-1)

if __name__ == '__main__':
    makeMatrix(sys.argv[1:])
