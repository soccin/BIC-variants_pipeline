#!/opt/common/python

from __future__ import print_function
import sys
import os
import fnmatch
import argparse
from collections import OrderedDict

def printMajorContamination(het,allSamples,pre,outdir):
    FILE = os.path.join(outdir,pre+'_MajorContamination.txt')
    try:
        with open(FILE,'w') as out:
            print("\t".join(['Sample','PerHeterozygousPos']),file=out)
            for s in allSamples:
                print("\t".join([s,str(het[s])]),file=out)
    except IOError:
        print("Error writing file: "+FILE,file=sys.stderr)
        sys.exit(1)
    return

def printMinorContamination(homPos,allSamples,pre,pairs,outdir):
    FILE = os.path.join(outdir,pre+'_MinorContamination.txt')
    try:
        with open(FILE,'w') as out:
            print("\t".join(['Sample','AvgMinorHomFreq']),file=out)
            for s in homPos:
                fracs = []
                if s in pairs and not pairs[s].lower() == 'na': ## sample is tumor and has a matched normal
                    norm = pairs[s]
                    for loc in homPos[s]:
                        if loc in homPos[norm]:# and homPos[norm][loc] < 0.10:
                            fracs.append(homPos[s][loc])
                else:
                    fracs = homPos[s].values()
                try:
                    f = float(sum(fracs))/len(fracs)
                    #f = float(sum(homPos[s]))/len(homPos[s])
                except ZeroDivisionError:
                    f = 0.00
                print("\t".join([s,str(f)]),file=out)
    except IOError:
        print("Error writing file: "+FILE,file=sys.stderr)
        sys.exit(1)
    return

def printSummary(finalOut,allSamples,pre,outdir):
    FILE = os.path.join(outdir,pre+'_FingerprintSummary.txt')
    try:
        with open(FILE,'w') as out:
            header = ['Locus']
            for s in allSamples:
                header += [s+'_Counts',s+'_Genotypes',s+'_MinorAlleleFreq']
            print("\t".join(header),file=out)
            for loc in finalOut:
                row = [loc]
                if len(finalOut[loc].keys())>0:
                    for s in allSamples:
                        row += finalOut[loc][s]
                    print("\t".join(row),file=out)
                else:
                    print("WARNING: NO SAMPLES HAD COUNTS AT LOCUS "+loc)
    except IOError:
        print("Error writing file: "+FILE,file=sys.stderr)
        sys.exit(1)
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
        print("Error reading file:",fpFile,file=sys.stderr)
        sys.exit(1)
    return fps

def getSampleName(filePath,pre):
    '''
    Extract sample name from file name

    ****WARNING: THIS FUNCTION IS SPECIFIC TO VARIANTS PIPELINE!! WILL NEED TO BE MODIFIED FOR OTHER APPLICATIONS!!***

    '''
    #pattern = '_indelRealigned_recal_'
    pattern = '__1'
    #pattern = '_bc'
    fname = filePath.rstrip('/').split('/')[-1]
    if not pattern in fname:
        print("ERROR: search pattern '"+pattern+"' is not in file name "+fname+" so we have no way to extract sample name!",file=sys.stderr)
        sys.exit(1)
    #start = fname.index(pattern)+len(pattern)
    #end = fname.index('.bam')
    start = 0
    end = fname.rindex(pattern)
    return fname[start:end]

def makeGroups(file):
    groups = {}
    try:
        with open(file,'r') as f:
            for line in f:
                s,grp = line.strip().split("\t")
                s = s.strip()
                grp = grp.strip()
                if not s in groups:
                    groups[s] = []
                groups[s].append(grp)
    except IOError:
        print("ERROR: Could not open grouping file "+file,file=sys.stderr)
        exit(1)
    return groups

def makePairs(file):
    pairs = {}
    try:
        with open(file,'r') as f:
            for line in f:
                n,t = line.strip('\n').split('\t')
                n = n.strip()
                t = t.strip()
                if t.lower() == 'na':
                    continue
                if "pool" in n.lower():
                    continue
                if not t in pairs:
                    pairs[t] = []
                pairs[t] = n
    except IOError:
        print("ERROR: Could not open pairing file "+file,file=sys.stderr)
        exit(1)
    except ValueError:
        print("ERROR: Pairing file "+file+"format invalid.",file=sys.stderr)
        exit(1)
    return pairs

def printSampleErrors(contamination,pre,groups,outdir):
    '''
    Print file of unexpected sample matches and a file of unexpected sample mismatches
    based on sample grouping file
    '''

    allSamps = contamination.keys()
    ueMatches = []
    ueMismatches = []
    for rSamp in allSamps:
        for qSamp in allSamps[allSamps.index(rSamp)+1:]:
            row = []
            fda = []
            if len(list(set(groups[rSamp]) & set(groups[qSamp]))) == 0:
                ## unexpected matches - alleleFrac < 0.05, but samples are in different groups
                if float(contamination[rSamp][qSamp]) < 0.05:
                    row += [rSamp,qSamp]
                    fda.append(str(contamination[rSamp][qSamp]))
                if float(contamination[qSamp][rSamp]) < 0.05:
                    if not row:
                        row += [qSamp,rSamp]
                    fda.append(str(contamination[qSamp][rSamp]))
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
                    fda.append(str(contamination[qSamp][rSamp]))
                if row and fda:
                    ueMismatches.append("\t".join(row+[",".join(fda)]))
    try:
        with open(os.path.join(outdir,pre+'_UnexpectedMatches.txt'),'w') as out:
            for x in ueMatches:
                print(x,file=out)
    except IOError:
        print("Error writing file "+os.path.join(outdir,pre+'_UnexpectedMatches.txt'),file=sys.stderr)
        sys.exit(1)
    try:
        with open(os.path.join(outdir,pre+'_UnexpectedMismatches.txt'),'w') as out:
            for x in ueMismatches:
                print(x,file=out)
    except IOError:
        print("Error writing file "+os.path.join(outdir,pre+'_UnexpectedMismatches.txt'),file=sys.stderr)
        sys.exit(1)
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
                if len(fpSummary[loc])>0:
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
    FILE = os.path.join(outdir,pre+'_DiscordantHomAlleleFractions.txt')
    ## print matrix to file
    allSamples = contamination.keys()
    try:
        allSamples.sort()
        with open(FILE,'w') as out:
            header = "\t".join(['']+allSamples)
            print(header,file=out)
            for r in allSamples:
                row = [r]
                for q in allSamples:
                    row.append(str(contamination[r][q]))
                print("\t".join(row),file=out)
    except IOError:
        print("Error writing file " +FILE,file=sys.stderr)
        sys.exit(1)
    return

def getContamination(docFiles,fpFile,pre,group,pairs,outdir):
    '''
    This script has four goals:
       1) print fingerprint summary file that contains allele counts, 
          genotype, and minor allele freq for each sample at each FP locus
       2) calculate and print avg allele freq at homozygous positions for each sample (minor contamination)
       3) calculate and print fraction of heterozygous positions for each sample
       4) calculate fraction of discordant alleles for all possible sample pairs
    '''

    fps = indexFPgenotypes(fpFile)
    homPos = {}     ## key: sample; value: dictionary where key: locus with hom, value:freq
    het = {}        ## key: sample; value: count of heterozygous positions
    finalOut = OrderedDict()   ## key: locus; value: dictionary in which key: sample; value = [counts,genotype,minAlleleFreq]
    allSamples = []    ## list of all samples
    homGTs = ['AA','CC','GG','TT']
    order = 'ACGT' ## order of bases in DepthOfCoverage output

    for doc in docFiles:
        sample = getSampleName(doc,pre)
        if not sample:
            continue 
        print("Parsing counts file for sample: "+sample)
        if sample in allSamples:
            print("ERROR: Duplicate sample found: "+sample,file=sys.stderr)
            continue
        allSamples.append(sample)
        homPos[sample] = {}
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
                counts = baseCounts.split()
                ## get allele 1, allele 2, and calculate frequency of each 
                a1,a2 = fps[locus][:2]
                if len(a1)>1 or len(a2)>1:
                    continue
                tot += 1
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
                    if frac < 0.10:
                        gt = a2+a2
                else:
                    frac = a2freq/(a1freq+a2freq)
                    if frac < 0.10:
                        gt = a1+a1
                ## store final genotype and fraction of minor allele
                row += [gt,str(frac)]
                ## store counts, genotype, fraction of minor allele for this sample at this locus
                finalOut[locus][sample] = row
                ## keep track of minor allele freqs at homozygous positions for minor contamination
                ## and number of heterozygous positions for major contamination
                #if gt in homGTs and frac < 0.1:
                if frac < 0.10:
                    homPos[sample][locus] = frac
                elif not gt == "--":
                    numHet += 1
            ## store fraction of heterozygous positions for this sample
            het[sample] = float(numHet/tot)
    ## print output files
    printMinorContamination(homPos,allSamples,pre,pairs,outdir)
    printMajorContamination(het,allSamples,pre,outdir)
    printSummary(finalOut,allSamples,pre,outdir)

    contamination = getDiscordantAlleleFractions(finalOut,allSamples)
    printDiscordantAlleleFractions(contamination,pre,outdir)
    printSampleErrors(contamination,pre,groups,outdir) 

def findFiles(rootDir,pattern):
    """
    return a list of full paths to files that match pattern entered
    """
    filepaths = []
    for path, dirs, files in os.walk(os.path.abspath(rootDir)):
        if fnmatch.filter(files, pattern):
            for file in fnmatch.filter(files, pattern):
                filepaths.append(os.path.join(path,file))
    if not filepaths:
        print("ERROR: No files matching pattern "+pattern+" found", file=sys.stderr)
        sys.exit(1)
    return filepaths

def get_args(sysargs):
    parser = argparse.ArgumentParser()

    parser.add_argument('-pattern',required=True,help='pattern to use to find DepthOfCoverage output files')
    parser.add_argument('-pre',required=True,help='project prefix')
    parser.add_argument('-fp',required=True,help='fingerprint genotypes file')
    parser.add_argument('-group',required=True,help='sample grouping file')
    parser.add_argument('-pair',required=True,help='sample pairing file')
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
    groups = makeGroups(args.group)
    pairs = makePairs(args.pair)
    getContamination(docFiles,args.fp,args.pre,groups,pairs,args.outdir)
