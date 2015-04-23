#!/opt/common/python

import sys
import os
import fnmatch
import argparse
from collections import OrderedDict

def usage():
    print>>sys.stderr,"%prog pre patternToSearchForDOCoutputFiles fingerprintGenotypesFile"
    sys.exit(-1)

def printMajorContamination(het,allSamples,pre,outdir):
    FILE = os.path.join(outdir,pre+'_MajorContamination.txt')
    try:
        with open(FILE,'w') as out:
            print>>out,"\t".join(['Sample','PerHeterozygousPos'])
            for s in allSamples:
                print>>out,"\t".join([s,str(het[s])])
    except IOError:
        print >> sys.stderr, "Error writing file: %s" %FILE
    return

def printMinorContamination(homPos,allSamples,pre,outdir):
    FILE = os.path.join(pre+'_MinorContamination.txt')
    try:
        with open(FILE,'w') as out:
            print>>out,"\t".join(['Sample','AvgMinorHomFreq'])
            for s in allSamples:
                f = float(sum(homPos[s]))/len(homPos[s])
                print>>out,"\t".join([s,str(f)])
    except IOError:
        print >> sys.stderr, "Error writing file: %s" %FILE
    return

def printSummary(finalOut,allSamples,pre,outdir):
    FILE = os.path.join(pre+'_FingerprintSummary.txt')
    try:
        with open(FILE,'w') as out:
            header = ['Locus']
            for s in allSamples:
                header += [s+'_Counts',s+'_Genotypes',s+'_MinorAlleleFreq']
            print>>out,"\t".join(header)
            for loc in finalOut:
                row = [loc]
                for s in allSamples:
                    row += finalOut[loc][s]
                print>>out,"\t".join(row)
    except IOError:
        print >> sys.stderr, "Error writing file: %s" %FILE
    return

def indexFPgenotypes(fpFile):
    ## index fingerprint genotypes by locus
    fps = {}  ## key: locus; value: [a1,a2]
              ## e.g., fps['1:123456'] = ['C','T']
    try:
        with open(fpFile,'r') as fp:
            for line in fp:
                loc,gt = line.strip().split()
                fps[loc] = gt.split("/")
    except IOError:
        print >> sys.stderr, "Error reading file %s" %fpFile
    return fps

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
    '''
    if not pre in fname:
        print >> sys.stderr, "ERROR: Prefix %s is not in file name %s, so we have no way to extract sample name!" %(pre,fname)
        return None
    return fname[:fname.index(pre)-1] 
    '''

def printSampleErrors(contamination,pre,group,outdir):
    '''
    Print file of unexpected sample matches and a file of unexpected sample mismatches
    based on sample grouping file
    '''
    groups = {}
    # store sample grouping
    try:
        with open(group,'r') as g:
            for line in g:
                s,grp = line.strip().split("\t")
                if s in groups:
                    print>>sys.stderr,"ERROR: sample %s belongs to more than one group!" %s
                    continue
                groups[s] = grp
    except IOError:
        print>>sys.stderr, "Error reading grouping file %s" %group

    allSamps = contamination.keys()
    ueMatches = []
    ueMismatches = []
    for rSamp in allSamps:
        for qSamp in allSamps[allSamps.index(rSamp)+1:]:
            row = []
            fda = []
            #if not groups[rSamp] == groups[qSamp]:
            if not groups[rSamp[:rSamp.index('_bc')]] == groups[qSamp[:qSamp.index('_bc')]]:
                ## unexpected matches - alleleFrac < 0.05, but samples are in different groups
                if float(contamination[rSamp][qSamp]) < 0.05:
                    row += [rSamp,qSamp]
                    fda.append(str(contamination[rSamp][qSamp]))
                if float(contamination[qSamp][rSamp]) < 0.05:
                    if not row:
                        row += [qSamp,rSamp]
                    fda.append(contamination[qSamp][rSamp])
                if row and fda:
                    ueMatches.append("\t".join(row+[",".join(fda)]))
            else:
                ## unexpected mismatches - alleleFrac > 0.05, but samples are in the same group
                if float(contamination[rSamp][qSamp]) > 0.05:
                    row += [rSamp,qSamp]
                    fda.append(str(contamination[rSamp][qSamp]))
                if float(contamination[qSamp][rSamp]) > 0.05:
                    if not row:
                        row += [qSamp,rSamp]
                    fda.append(contamination[qSamp][rSamp])
                if row and fda:
                    ueMismatches.append("\t".join(row+[",".join(fda)]))
    try:
        with open(os.path.join(outdir,pre+'_UnexpectedMatches.txt'),'w') as out:
            for x in ueMatches:
                print>>out,x
    except IOError:
        print>>sys.stderr,"Error writing file %s" %s
    try:
        with open(os.path.join(pre+'_UnexpectedMismatches.txt'),'w') as out:
            for x in ueMismatches:
                print>>out,x
    except IOError:
        print>>sys.stderr,"Error writing file %s" %s

    return

def getDiscordantAlleleFractions(fpSummary,allSamples):
    '''
    Create and print a matrix of discorant allele fractions for all possible
    sample pairs
    '''
    
    #fpSummary[locus][sample] = [counts,gt,freq]
    homGTs = ['AA','GG','TT','CC']
    contamination = {}

    for rSamp in sorted(allSamples):
        contamination[rSamp] = {}
        for qSamp in sorted(allSamples):
            contamination[rSamp][qSamp] = 0           
            rCount = qCount = 0
            for loc in fpSummary:
                rGT = fpSummary[loc][rSamp][1]
                qGT = fpSummary[loc][qSamp][1]
                if rGT in homGTs:
                    rCount += 1
                    if qGT in homGTs and not qGT == rGT:
                        qCount += 1
            if rCount:
                contamination[rSamp][qSamp] = format(float(qCount)/rCount, '.2f')
    return contamination

def printDiscordantAlleleFractions(contamination,pre,outdir):
    FILE = os.path.join(pre+'_DiscordantAlleleFractions.txt')
    ## print matrix to file
    allSamples = contamination.keys()
    try:
        allSamples.sort()
        with open(FILE,'w') as out:
            header = "\t".join(['']+allSamples)
            print >> out, header
            for r in allSamples:
                row = [r]
                for q in allSamples:
                    row.append(str(contamination[r][q]))
                print >> out, "\t".join(row)
    except IOError:
        print >> sys.stderr, "Error writing file %s" %FILE
    return

def getContamination(docFiles,fpFile,pre,group):
    '''
    This script has four goals:
       1) print fingerprint summary file that contains allele counts, 
          genotype, and minor allele freq for each sample at each FP locus
       2) calculate and print avg allele freq at homozygous positions for each sample (minor contamination)
       3) calculate and print fraction of heterozygous positions for each sample
       4) calculate fraction of discordant alleles for all possible sample pairs
    '''

    fps = indexFPgenotypes(fpFile)
    homPos = {}     ## key: sample; value: list of frequencies at homozygous positions
    het = {}        ## key: sample; value: count of heterozygous positions
    finalOut = OrderedDict()   ## key: locus; value: dictionary in which key: sample; value = [counts,genotype,minAlleleFreq]
    allSamples = []    ## list of all samples
    homGTs = ['AA','CC','GG','TT']
    order = 'ACGT' ## order of bases in DepthOfCoverage output

    for doc in docFiles:
        #print >>sys.stderr,"Parsing counts file %s" %doc
        sample = getSampleName(doc,pre)
        if not sample:
            continue 
        print >>sys.stderr,"Parsing counts file for sample: %s" %sample
        if sample in allSamples:
            print>>sys.stderr,"ERROR: Duplicate sample found: %s" %sample
            continue
        allSamples.append(sample)
        homPos[sample] = []
        het[sample] = 0

        with open(doc,'r') as c:
            ## output from GATK DepthOfCoverage:
            ## 1:4849325	313	313.00	313	A:313 C:0 G:0 T:0 N:0 
            a1 = a2 = 0
            a1freq = a2freq = 0
            numHet = tot = 0.0 ## number of heterozygous positions and total num positions
            ## skip header
            c.readline()
            for line in c:
                row = [] ## store: [1:4849324	A:10 C:14	0.416666]
                locus,x,y,z,baseCounts = line.strip().split("\t")
                if not locus in fps:
                    continue
                if not locus in finalOut:
                    finalOut[locus] = {}
                tot += 1
                counts = baseCounts.split()
                ## get allele 1, allele 2, and calculate frequency of each 
                a1,a2 = fps[locus][:2]
                a1freq = float(counts[order.index(a1)].split(":")[1])
                a2freq = float(counts[order.index(a2)].split(":")[1])
                ## store base counts of the two alleles for this locus from the FP genotype file 
                row.append(counts[order.index(a1)] + " " + counts[order.index(a2)])
                gt = a1+a2 
                if a1freq == 0 and a2freq == 0:
                    frac = 0
                    gt = "--"
                elif a1freq < a2freq:
                    frac = a1freq/(a1freq+a2freq)
                    if frac < 0.1:
                        gt = a2+a2
                else:
                    frac = a2freq/(a1freq+a2freq)
                    if frac < 0.1:
                        gt = a1+a1
                ## store final genotype and fraction of minor allele
                row += [gt,str(frac)]
                ## store counts, genotype, fraction of minor allele for this sample at this locus
                finalOut[locus][sample] = row
                ## keep track of minor allele freqs at homozygous positions for minor contamination
                ## and number of heterozygous positions for major contamination
                if gt in homGTs:
                    homPos[sample].append(frac)
                elif not gt == "--":
                    numHet += 1
            ## store fraction of heterozygous positions for this sample
            het[sample] = float(numHet/tot)
    ## print output files
    printMinorContamination(homPos,allSamples,pre,outdir)
    printMajorContamination(het,allSamples,pre,outdir)
    printSummary(finalOut,allSamples,pre,outdir)

    contamination = getDiscordantAlleleFractions(finalOut,allSamples)
    printDiscordantAlleleFractions(contamination,pre,outdir)
    printSampleErrors(contamination,pre,group,outdir) 

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

def get_args(sysargs):
    parser = argparse.ArgumentParser()

    parser.add_argument('-pattern',required=True,help='pattern to use to find DepthOfCoverage output files')
    parser.add_argument('-pre',required=True,help='project prefix')
    parser.add_argument('-fp',required=True,help='fingerprint genotypes file')
    parser.add_argument('-group',required=True,help='sample grouping file')
    parser.add_argument('-dir',help='directory to search; cwd by default')
    parser.add_argument('-outdir',help='output directory')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args(sys.argv)
    if not args.dir:
        dir=os.getcwd()
    else:
        dir=args.dir
    docFiles = findFiles(dir,args.pattern)
    getContamination(docFiles,args.fp,args.pre,args.group,args.outdir)
