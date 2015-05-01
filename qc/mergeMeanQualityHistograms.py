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
## Example: /usr/bin/python makeCountMatrix.py /ifs/res/liang/RNASeq/Proj2983_MassagueJ "*htseq.count*" Proj2983_MassagueJ_htseq.count_allSamples.txt
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


def printMatrix(matrix,allSamps,rqualOutFile,oqualOutFile):
    """
    """

    header = "\t".join(["Cycle"]+allSamps)

    with open(rqualOutFile,'w') as rout:
        with open(oqualOutFile,'w') as oout:
            print>>rout,header
            print>>oout,header
            for cycle in matrix.keys():
                for samp in allSamps:
                    if not samp in matrix[cycle]['rqual'] and not samp in matrix[cycle]['oqual']:
                        matrix[cycle]['rqual'][samp] = 0 
                        matrix[cycle]['oqual'][samp] = 0
                print>>rout,"\t".join([str(x) for x in [cycle]+[matrix[cycle]['rqual'][samp] for samp in allSamps]])
                print>>oout,"\t".join([str(x) for x in [cycle]+[matrix[cycle]['oqual'][samp] for samp in allSamps]])
    return

def makeMatrix(args):
    """
    Find files to parse, create one matrix of all counts and print 
    matrix to file
    """

    if len(args) == 4:
        rootDir,filePattern,rqualOutFile,oqualOutFile = args
    else:
        usage()
        sys.exit(1)

    ## store all values in an ordered dict, keyed by sample
    matrix = OrderedDict()
    allSamps = []

    ## find all cutadapt stats files using pattern 
    files = findFiles(rootDir,filePattern)

    if files:
        print>>sys.stderr, "\nCombining the following files:\n"
        for file in files:
            print>>sys.stderr, file
            if os.stat(file)[6]==0: ## file is empty
                print>>sys.stderr, "WARNING: This file is empty!"
            else:

                samp = file[file.find("MeanQualityByCycle_")+19:file.find(".txt")] 
                if samp in allSamps:
                    print "ERROR: sample %s found multiple times!!" %samp
                    continue
                allSamps.append(samp)

                with open(file,'r') as fl:
                    while not "## HISTOGRAM" in fl.readline():
                       next
                    fl.readline()
                    for line in fl:
                        if len(line.strip())>0:
                            cycle,rqual,oqual = line.strip().split("\t")
                            if not cycle in matrix:
                                matrix[cycle] = {'rqual':{},'oqual':{}}
                            if not samp in matrix[cycle]['rqual'] and not samp in matrix[cycle]['oqual']:
                                matrix[cycle]['rqual'][samp] = rqual
                                matrix[cycle]['oqual'][samp] = oqual                         

        printMatrix(matrix,allSamps,rqualOutFile,oqualOutFile)            

    else:
        print>>sys.stderr, "\nNo files found matching pattern entered. Exiting.\n"
        sys.exit(-1)

if __name__ == '__main__':
    makeMatrix(sys.argv[1:])
