#!/opt/bin/python
import os
import sys
import argparse
import fnmatch
from collections import OrderedDict

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

def getSampleName(filePath,pre):
    '''
    Extract sample name from file name
    '''
    pattern = '_indelRealigned_recal_'
    fname = filePath.rstrip('/').split('/')[-1]
    if not pattern in fname:
        print>>sys.stderr, "ERROR: search pattern '%s' is not in file name %s, so we have no way to extract sample name!" %(pattern,fname)
        return None
    start = fname.index(pattern)+len(pattern)
    end = fname.index('.bam')
    return fname[start,end]

def mergeOutput(fileList,pre):
    coverage = OrderedDict()
    samples = []
    for file in fileList:
        sample = getSampleName(file,pre)
        if not sample:
            continue
        print>>sys.stderr, "Merging coverage for sample %s" %sample
        samples.append(sample)
        try:
            with open(file,'r') as f:
                f.readline() ## skip header
                for line in f:
                    cov = line.strip().split("\t")
                    if not cov[0] in coverage:
                        coverage[cov[0]] = {}
                    if sample in coverage[cov[0]]:
                        print>>sys.stderr,"ERROR: Duplicate sample found: %s" %sample
                        continue
                    coverage[cov[0]][sample] = cov[2]
        except IOError:
            print>>sys.stderr,"Error reading file %s" %file
    return coverage,samples

def printCoverageSummary(coverage,samples,pre,outdir):
    file = os.path.join(pre+"_CanonicalExonCoverage.txt")
    samples.sort()
    try:
        with open(file,'w') as out:
            ln = 0
            print>>out,"\t".join([str(ln),"Target"]+samples)
            for i in coverage:
                ln+=1
                row = [str(ln),i]
                for s in samples:
                    row.append(coverage[i][s])
                print>>out,"\t".join(row)
    except IOError:
        print>>sys.stderr,"Error writing file %s" %file
    return

def get_args(sysargs):
    parser = argparse.ArgumentParser()

    parser.add_argument('-pattern',required=True,help='pattern to use to find DepthOfCoverage output files')
    parser.add_argument('-pre',required=True,help='project prefix')
    parser.add_argument('-dir',help='directory to search; cwd by default')
    parser.add_argument('-outdir',help='output directory; cwd by default')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args(sys.argv)
    if not args.dir:
        dir = os.getcwd()
    else:
        dir = args.dir
    if not args.outdir:
        outdir = os.getcwd()
    else:
        outdir = args.outdir
    fileList = findFiles(dir,args.pattern)
    cov,samps = mergeOutput(fileList,args.pre)
    printCoverageSummary(cov,samps,args.pre,outdir)
