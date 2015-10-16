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

    return filepaths


def printLog(allOut,outFile):
    """
    """
    with open(outFile,'w') as out:
        if allOut:
            for line in allOut:
                print>>out,line    

    return

def catLogs(args):
    """
    Find files to parse, create one matrix of all counts and print 
    matrix to file
    """

    if len(args) == 3:
        rootDir,filePattern,outFile = args
    else:
        usage()
        sys.exit(1)

    filePattern = '*'+filePattern
    ## find all cutadapt stats files using pattern 
    files = findFiles(rootDir,filePattern)

    allOut = []

    if files:
        print "\nCombining the following files:\n"
        for file in files:
            print file
            if not os.stat(file)[6]==0: ## file is empty
                with open(file,'r') as fl:
                    for line in fl:
                        allOut.append(line.strip())
     
        printLog(allOut,outFile)            

    else:
        print>>sys.stderr, "\nNo files found matching pattern entered. Exiting.\n"
        sys.exit(-1)

if __name__ == '__main__':
    catLogs(sys.argv[1:])
