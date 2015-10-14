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
        #else:
        #    if "Proj" in path.split("/")[-1] and "-" in path.split("/")[-1]:
        #        print>>sys.stderr, "WARNING: No files matching pattern %s found in %s" %(pattern, path)

    return filepaths


def printMatrix(matrix,outFile):
    """
    """

    header = "\t".join(["Sample"] + matrix[matrix.keys()[0]].keys())
    with open(outFile,'w') as out:
        print>>out,header
        for sample,info in matrix.items():
            print>>out,"\t".join([sample]+[str(info[x]) for x in info.keys()])    

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
                print "WARNING: This file is empty!"
            else:
                samp = file[file.find('intFiles')+9:].split("/")[0]
                fName = file.split("/")[-1]
                if not samp in matrix:
                    matrix[samp] = OrderedDict()
                    matrix[samp]['Total_Reads_Processed'] = 0	
                    matrix[samp]['Total_Bases_Processed'] = 0
                    matrix[samp]['Min_Raw_Quality'] = 1000
                    matrix[samp]['Max_Raw_Quality'] = -1
                    matrix[samp]['Quality_Scale'] = ''
                    matrix[samp]['Num_Quality_Score_Truncated'] = 0	
                    matrix[samp]['Perc_Quality_Score_Truncated'] = 0.0
                    matrix[samp]['Num_Low_Quality_Score_Truncated'] = 0	
                    matrix[samp]['Perc_Low_Quality_Score_Truncated'] = 0.0
                    matrix[samp]['Num_High_Quality_Score_Truncated'] = 0	
                    matrix[samp]['Perc_High_Quality_Score_Truncated'] = 0.0
                                    

                with open(file,'r') as fl:
                    fl.readline()
                    for line in fl:
                        f1,f2,trp,tbp,mnrq,mxrq,qs,nqst,pqst,nlqst,plqst,nhqst,phqst = line.strip().split("\t")
                        matrix[samp]['Total_Reads_Processed'] += int(trp)
                        matrix[samp]['Total_Bases_Processed'] += int(tbp)
                        matrix[samp]['Min_Raw_Quality'] = min(matrix[samp]['Min_Raw_Quality'],int(mnrq))
                        matrix[samp]['Max_Raw_Quality'] = max(matrix[samp]['Max_Raw_Quality'],int(mxrq))
                        if matrix[samp]['Quality_Scale'] and not qs == matrix[samp]['Quality_Scale']:
                            print>>sys.stderr,"ERROR: different quality scales found within sample "+samp+": "+maxtrix[samp]['Quality_Scale']+" and "+qs
                            sys.exit(-1)
                        matrix[samp]['Quality_Scale'] = qs
                        matrix[samp]['Num_Quality_Score_Truncated'] += int(nqst)
                        matrix[samp]['Perc_Quality_Score_Truncated'] = format((float(matrix[samp]['Num_Quality_Score_Truncated'])/matrix[samp]['Total_Bases_Processed'])*100,'.2f')
                        matrix[samp]['Num_Low_Quality_Score_Truncated'] += int(nlqst)
                        matrix[samp]['Perc_Low_Quality_Score_Truncated'] = format((float(matrix[samp]['Num_Low_Quality_Score_Truncated'])/matrix[samp]['Total_Bases_Processed'])*100,'.2f')
                        matrix[samp]['Num_High_Quality_Score_Truncated'] += int(nhqst)
                        matrix[samp]['Perc_High_Quality_Score_Truncated'] = format((float(matrix[samp]['Num_High_Quality_Score_Truncated'])/matrix[samp]['Total_Bases_Processed'])*100,'.2f')
 
        printMatrix(matrix,outFile)            

    else:
        print>>sys.stderr, "\nNo files found matching pattern entered. Exiting.\n"
        sys.exit(-1)

if __name__ == '__main__':
    makeMatrix(sys.argv[1:])
