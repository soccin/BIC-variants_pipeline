#!/opt/common/CentOS_6/python/python-2.7.8/bin/python

import sys
import os
import re
import fnmatch

def indexCosmic(cosmicVcf):
    cosmicMutations = {}
    with open(cosmicVcf,'rU') as cf:
        header = cf.readline()
        for line in cf:
            if line.strip('\n') == '' or line[0] == '#':
                continue
            pts = line.strip('\n').split('\t')
            key = ':'.join([pts[0],pts[1],pts[3],pts[4]])
            cosmicMutations[key] = '|'.join([pts[2],pts[7]])
    return cosmicMutations

def indexTumors(pairingFile):
    pairs = {}
    with open(pairingFile,'r') as pairing:
        for line in pairing:
            n,t = line.strip('\n').split('\t')
            if t in pairs:
                print>>sys.stderr,"WARNING: Tumor",t,"is paired with multiple normals. Using",n+"."
            pairs[t] = n
    return pairs

def getSampleName(filePath):
    if 'intFiles' in filePath:
        samp = filePath[filePath.find('intFiles')+9:].split("/")[0] 
        return samp
    else:
        print>>sys.stderr,"ERROR: Don't know how to extract file name from path",filePath
        sys.exit(-1)

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
                print "Found file",os.path.join(path,file)
        else:
            if "Proj" in path.split("/")[-1] and "-" in path.split("/")[-1]:
                print>>sys.stderr, "WARNING: No files matching pattern %s found in %s" %(pattern, path)

    return filepaths


def writeHotspotFile(cosmicMuts,pairing,mutationFiles,outFile):
    normPos = {}
    tumorPos = {}

    for f in mutationFiles:
        sampleName = getSampleName(f)
        if sampleName in pairing:
            continue ## work with only normals for now
        if not sampleName in pairing.values():
            print>>sys.stderr, "WARNING: Sample",sampleName,"is not in pairing file. Don't know if it's a normal or tumor. Skipping."
            continue
        with open(f,'r') as fl:
            for line in fl:
                if line.startswith("Ref"):
                    continue
                pts = line.strip('\n').split('\t')
                name,x,chr,pos,ref,alt = pts[:6]
                dp,y,ad,vf = pts[6:10]
                key = ':'.join([chr,pos,ref,alt])

                if int(dp) > 0 and float(vf) > 0.01 and int(ad) >= 3 and key in cosmicMuts:
                    if not key in normPos:
                        normPos[key] = []
                    normPos[key].append('\t'.join(['NORMAL',sampleName,cosmicMuts[key],chr,pos,ref,alt,dp,ad,vf]))

    for f in mutationFiles:
        sampleName = getSampleName(f)
        if not sampleName in pairing:
            continue ## work with tumors now
        with open(f,'r') as fl:
            for line in fl:
                if line.startswith("Ref"):
                    continue
                pts = line.strip('\n').split('\t')
                name,x,chr,pos,ref,alt = pts[:6]
                dp,y,ad,vf = pts[6:10]
                key = ':'.join([chr,pos,ref,alt])

                if key in normPos:
                    for i in normPos[key]:
                        if not pairing[sampleName].lower() == 'na' and not 'pool' in pairing[sampleName].lower():
                            matchedNormal = pairing[sampleName]
                            if not key in tumorPos:
                                tumorPos[key] = []
                            tumorPos[key].append('\t'.join(['MATCHED_TM',sampleName,cosmicMuts[key],chr,pos,ref,alt,dp,ad,vf]))
                        elif int(dp) > 0 and float(vf) > 0.01 and int(ad) >= 3 :
                            if not key in tumorPos:
                                tumorPos[key] = []
                            tumorPos[key].append('\t'.join(['TM',sampleName,cosmicMuts[key],chr,pos,ref,alt,dp,ad,vf]))

    with open(outFile,'w') as out:
        if len(normPos) == 0:
            print>>out,"No hotspots detected in any normals"
        else:
            for val in normPos.values():
                for x in val:
                    print>>out,x
            for val in tumorPos.values():
                for x in val:
                    print>>out,x

    return

def usage():
    print>>sys.stderr,"%prog rootDir patternToSearch cosmicFile pairingFile outFile"
    return

if __name__ == '__main__':
    ## hotspots_in_normals.py rootDir pattern cosmicFile pairingFile outFile
    if len(sys.argv) != 6:
        usage()
        sys.exit(1)

    rootDir,pattern,cosmicFile,pairingFile,outFile = sys.argv[1:]
    cosmicMuts = indexCosmic(cosmicFile)
    pairing = indexTumors(pairingFile)
    mutationFiles = findFiles(rootDir,pattern)
    writeHotspotFile(cosmicMuts,pairing,mutationFiles,outFile)
    sys.exit(0)
