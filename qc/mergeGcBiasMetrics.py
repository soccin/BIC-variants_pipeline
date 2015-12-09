#!/usr/bin/python

import sys
import re
import os
import fnmatch
from collections import OrderedDict


def usage():
    print "/opt/bin/python mergeGcBiasMetrics.py rootDir patternToSearch outputFileName"
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

def printMatrix(matrix,allSamps,outFile):
    """
    """

    header = "\t".join(["GC"]+allSamps)

    with open(outFile,'w') as out:
        print>>out,header
        for gc in matrix.keys():
            for samp in allSamps:
                if not samp in matrix[gc]:
                    matrix[gc][samp] = 0 
            print>>out,"\t".join([str(x) for x in [gc]+[matrix[gc][samp] for samp in allSamps]])
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
    allSamps = []

    filePattern = '*'+filePattern
    ## find all gc metrics files using pattern 
    files = findFiles(rootDir,filePattern)

    if files:
        print>>sys.stderr, "\nCombining the following files:\n"
        for file in files:
            print>>sys.stderr, file
            if os.stat(file)[6]==0: ## file is empty
                print>>sys.stderr, "WARNING: This file is empty!"
            else:
                fname = file.split("/")[-1].replace(".txt","")
                pat_idx = fname.index(filePattern.replace("*",""))
                samp = fname[pat_idx+len(filePattern.replace("*",""))+1:] ### WARNING: this is dumb as it assumes a certain pattern and should only be used with current naming convention of GCbias metrics files!
                if samp in allSamps:
                    print "ERROR: sample %s found multiple times!!" %samp
                    continue
                else:
                    allSamps.append(samp)
                with open(file,'r') as fl:
                    for line in fl:
                        if "#" in line or len(line.strip()) == 0:
                            continue
                        if line.startswith("GC"):
                            header = line.strip("\n").split("\t")
                            gc_idx = header.index("GC")
                            nc_idx = header.index("NORMALIZED_COVERAGE")
                            continue
                        cols = line.strip().split("\t")
                        gc = cols[gc_idx]
                        nc = cols[nc_idx]
                        if not gc in matrix:
                            matrix[gc] = {}
                        if not samp in matrix[gc]:
                            matrix[gc][samp] = nc                        

        printMatrix(matrix,allSamps,outFile)            

    else:
        print>>sys.stderr, "\nNo files found matching pattern entered. Exiting.\n"
        sys.exit(-1)

if __name__ == '__main__':
    makeMatrix(sys.argv[1:])
