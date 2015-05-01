'''
Created on 03/10/2014
@Ronak Shah

'''

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from collections import defaultdict
import sys

def main():
   parser = argparse.ArgumentParser(prog='AnnotateSVs.py', description='Add Annotation to the Structural Variants', usage='%(prog)s [options]')
   parser.add_argument("-r", "--repeatFile", action="store", dest="rrFilename", required=True, metavar='RepeatRegionFile.tsv', help="Location of the Repeat Region Bed File") 
   parser.add_argument("-d", "--dgvFile", action="store", dest="dgvFilename", required=True, metavar='DGvFile.tsv', help="Location of the Database of Genomic Variants Bed File")
   parser.add_argument("-c", "--cosmicConsensusFile", action="store", dest="ccFilename", required=True, metavar='CosmicConsensus.tsv', help="Location of the Cosmic Consensus TSV file")
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-s", "--svFile", action="store", dest="svFilename", required=True, metavar='SVfile.txt', help="Location of the structural variant file to be annotated")
   parser.add_argument("-o", "--outputFilePrefix", action="store", dest="outFilePrefix", required=True, metavar='AnnotatedSV', help="Full path with prefix name for the output file")
   args = parser.parse_args()
   outFileTxt = args.outFilePrefix + ".txt"
   outFileExl = args.outFilePrefix + ".xlsx"
   outFileJson = args.outFilePrefix + ".json"
   if args.verbose:
        print "Reading %s..." % args.svFilename
        data = ReadSVFile(args.svFilename, args.outFilePrefix ,args.verbose)
        print "Finished Reading %s" % args.svFilename
        print "Reading %s..." % args.rrFilename
        repeatRegionDict = ReadRepeatFile(args.rrFilename, args.verbose)
        print "Finished Reading %s" % args.rrFilename
        print "Reading %s..." % args.dgvFilename
        dgvDict = ReadRepeatFile(args.dgvFilename, args.verbose)
        print "Finished Reading %s" % args.dgvFilename
        data['repName-repClass-repFamily:-site1'] = "-"
        data['repName-repClass-repFamily:-site2'] = "-"
        data['CC_Chr_Band'] = "-"
        data['CC_Tumour_Types(Somatic)'] = "-"
        data['CC_Cancer_Syndrome'] = "-"
        data['CC_Mutation_Type'] = "-"
        data['CC_Translocation_Partner'] = "-"
        data['DGv_Name-DGv_VarType-site1'] = "-"
        data['DGv_Name-DGv_VarType-site2'] = "-"
        for count, row in data.iterrows():
            sv_chr1 = row.loc['Chr1']
            sv_pos1 = row.loc['Pos1']
            sv_chr2 = row.loc['Chr2']
            sv_pos2 = row.loc['Pos2']
            sv_gene1 = row.loc['Gene1']
            sv_gene2 = row.loc['Gene2']
            sv_type = row.loc['SV_Type']
            print "Processing Record:"
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (sv_chr1, sv_pos1, sv_chr2, sv_pos2, sv_gene1, sv_gene2, sv_type)
            # Repeat Region Data
            (rr_loc1, rr_loc2) = AnnotateRepeatRegion(args.verbose, count, row, repeatRegionDict)
            data.ix[count, 'repName-repClass-repFamily:-site1'] = "<=>".join(rr_loc1) 
            data.ix[count, 'repName-repClass-repFamily:-site2'] = "<=>".join(rr_loc2) 
            # Cosmic Consensus Data
            cc_SV = ReadCosmicCensusFile(args.ccFilename, args.verbose, count, row)
            ccA, ccB, ccC, ccD , ccE= ([] for i in range(5))
            for cc in cc_SV:
                ccData = cc.split('\t')
                ccA.append(ccData[0])
                ccB.append(ccData[1])
                ccC.append(ccData[2])
                ccD.append(ccData[3])
                ccE.append(ccData[4])
            data.ix[count, 'CC_Chr_Band'] = "<=>".join(ccA) 
            data.ix[count, 'CC_Tumour_Types(Somatic)'] = "<=>".join(ccB)  
            data.ix[count, 'CC_Cancer_Syndrome'] = "<=>".join(ccC) 
            data.ix[count, 'CC_Mutation_Type'] = "<=>".join(ccD)
            data.ix[count, 'CC_Translocation_Partner'] = "<=>".join(ccE)  
            # DGvData    
            (dgv_loc1, dgv_loc2) = AnnotateDGv(args.verbose, count, row, dgvDict)
            data.ix[count, 'DGv_Name-DGv_VarType-site1'] = "<=>".join(dgv_loc1)
            data.ix[count, 'DGv_Name-DGv_VarType-site2'] = "<=>".join(dgv_loc2)
   else:
        data = ReadSVFile(args.svFilename, args.verbose)
        repeatRegionDict = ReadRepeatFile(args.rrFilename, args.verbose)
        dgvDict = ReadRepeatFile(args.dgvFilename, args.verbose)
        for count, row in data.iterrows():
            (rr_loc1, rr_loc2) = AnnotateRepeatRegion(args.verbose, count, row, repeatRegionDict)
            cc_SV = ReadCosmicCensusFile(args.ccFilename, args.verbose, count, row)
            (dgv_loc1, dgv_loc2) = AnnotateDGv(args.verbose, count, row, dgvDict)
            data.ix[count, 'repName-repClass-repFamily:-site1'] = "<=>".join(rr_loc1) 
            data.ix[count, 'repName-repClass-repFamily:-site2'] = "<=>".join(rr_loc2)
            ccA, ccB, ccC, ccD = ([] for i in range(4))
            for cc in cc_SV:
                ccData = cc.split('\t')
                ccA.append(ccData[0])
                ccB.append(ccData[1])
                ccC.append(ccData[2])
                ccD.append(ccData[3])
            data.ix[count, 'CC_Chr_Band'] = "<=>".join(ccA) 
            data.ix[count, 'CC_Tumour_Types(Somatic)'] = "<=>".join(ccB)  
            data.ix[count, 'CC_Cancer_Syndrome'] = "<=>".join(ccC) 
            data.ix[count, 'CC_Mutation_Type'] = "<=>".join(ccD)
            data.ix[count, 'CC_Translocation_Partner'] = "<=>".join(ccE)  
            data.ix[count, 'DGv_Name-DGv_VarType-site1'] = "<=>".join(dgv_loc1)
            data.ix[count, 'DGv_Name-DGv_VarType-site2'] = "<=>".join(dgv_loc2) 
    
   # Print to TSV file
   data.to_csv(outFileTxt, sep='\t', index=False)
   # Print to Json
   data.to_json(outFileJson)
   # Print to Excel
   data.to_excel(outFileExl, sheet_name='Annotated_SVs', index=False)
            
def ReadRepeatFile(file, verbose):
    if(verbose):
        print ("Reading & Storing Repeat TSV file as dictionary")
    # Initialize dictionary of lists 
    dataDict = defaultdict(list)
    with open(file, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip('\n').split('\t')
            processedData = (data[0].replace('chr', ''))
            slicedData = data[1:]
            joinedData = '\t'.join(slicedData)
            dataDict[processedData].append(joinedData)    
    return dataDict        

def AnnotateRepeatRegion (verbose, count, sv, rrDict):
    if(verbose):
        print ("Checking Entry in Repeat data")
    # Initialize List to store repeat annotation
    list_svloc1 = []
    list_svloc2 = []
    # Read SV Data
    sv_chr1 = str(sv.loc['Chr1'])
    sv_pos1 = int(sv.loc['Pos1'])
    sv_chr2 = str(sv.loc['Chr2'])
    sv_pos2 = int(sv.loc['Pos2'])
    sv_type = sv.loc['SV_Type']
    # Traverse through Repeat Data Dict
    list_loc1 = rrDict.get(sv_chr1, "None")  # Get the values for the chromosome
    if(list_loc1 != "None"):  # Check if there are no keys with a particular chromosome
        for loc in list_loc1:  # For each location in all values check the overlap
            data = loc.split('\t') 
            rr_pos1 = int(data[0])
            rr_pos2 = int(data[1])
            if (rr_pos1 <= sv_pos1 <= rr_pos2):
                slicedData = data[2:]
                joinedData = '-'.join(slicedData)
                list_svloc1.append(joinedData)
    else:
        if(verbose):
            print "Chromosome ", sv_chr1, " is not there in the repeat dictionary"
    list_loc2 = rrDict.get(sv_chr2, "None")
    if(list_loc2 != "None"):
        for loc in list_loc2:
            data = loc.split('\t') 
            rr_pos1 = int(data[0])
            rr_pos2 = int(data[1])
            if (rr_pos1 <= sv_pos2 <= rr_pos2):
                slicedData = data[2:]
                joinedData = '-'.join(slicedData)
                list_svloc2.append(joinedData)
    else:
        if(verbose):
            print "Chromosome ", sv_chr2, " is not there in the repeat dictionary" 
    return (list_svloc1, list_svloc2)  
                                  
def ReadCosmicCensusFile (file, verbose, count, sv):
    if(verbose):
        print ("Checking Entry in Cosmic data")
    # Initialize List to store repeat annotation
    list_ccData = []
    sv_chr1 = sv.loc['Chr1']
    sv_pos1 = sv.loc['Pos1']
    sv_chr2 = sv.loc['Chr2']
    sv_pos2 = sv.loc['Pos2']
    sv_gene1 = str(sv.loc['Gene1'])
    sv_gene2 = str(sv.loc['Gene2'])
    sv_type = sv.loc['SV_Type']
   
    with open(file, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip('\n').split('\t')
            if(str(data[0]) == sv_gene1): 
                slicedData = getVar(data,[4,7,9,12,13])
                slicedProcessedData = []
                for sData in slicedData:
                    if(sData):
                        sData = "site1:" + sData
                        slicedProcessedData.append(sData)
                    else:
                        slicedProcessedData.append(" ")
                joinedData = '\t'.join(slicedProcessedData)
                list_ccData.append(joinedData)  
            if(str(data[0]) == sv_gene2):
                slicedData = getVar(data,[4,7,9,12,13]) 
                slicedProcessedData = []
                for sData in slicedData:
                    if(sData):
                        sData = "site2:" + sData   
                        slicedProcessedData.append(sData)
                    else:
                        slicedProcessedData.append(" ")
                joinedData = '\t'.join(slicedProcessedData)
                list_ccData.append(joinedData)             
    return list_ccData

def ReadDGvFile(file, verbose):
    if(verbose):
        print ("Reading & Storing DGV TSV file as dictionary")
    # Initialize dictionary of lists 
    dataDict = defaultdict(list)
    with open(file, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip('\n').split('\t')
            processedData = (data[0].replace('chr', ''))
            slicedData = data[1:]
            joinedData = '\t'.join(slicedData)
            dataDict[processedData].append(joinedData)    
    return dataDict        
         
def AnnotateDGv (verbose, count, sv, dgvDict):
    if(verbose):
        print ("Checking Entry in DGv data")
    # Initialize List to store repeat annotation
    list_svloc1 = []
    list_svloc2 = []
    # Read SV Data
    sv_chr1 = str(sv.loc['Chr1'])
    sv_pos1 = int(sv.loc['Pos1'])
    sv_chr2 = str(sv.loc['Chr2'])
    sv_pos2 = int(sv.loc['Pos2'])
    sv_type = sv.loc['SV_Type']
    # Traverse through DGv Data Dict
    list_loc1 = dgvDict.get(sv_chr1, "None")  # Get the values for the chromosome
    if(list_loc1 != "None"):  # Check if there are no keys with a particular chromosome
        for loc in list_loc1:  # For each location in all values check the overlap
            data = loc.split('\t') 
            dgv_pos1 = int(data[0])
            dgv_pos2 = int(data[1])
            if (dgv_pos1 <= sv_pos1 <= dgv_pos2):
                slicedData = getVar(data, [2, 8])
                joinedData = '-'.join(slicedData)
                list_svloc1.append(joinedData)
    else:
        if(verbose):
            print "Chromosome ", sv_chr1, " is not there in the repeat dictionary"        
    list_loc2 = dgvDict.get(sv_chr2, "None")
    if(list_loc2 != "None"):
        for loc in list_loc2:
            data = loc.split('\t') 
            dgv_pos1 = int(data[0])
            dgv_pos2 = int(data[1])
            if (dgv_pos1 <= sv_pos2 <= dgv_pos2):
                slicedData = getVar(data, [2, 8])
                joinedData = '-'.join(slicedData)
                list_svloc2.append(joinedData)
    else:
        if(verbose):
            print "Chromosome ", sv_chr2, " is not there in the repeat dictionary"  
    return (list_svloc1, list_svloc2)

def ReadSVFile (file, outFilePrefix, verbose):
    if(verbose):
        print ("Reading Structural Variant File")
        count = len(open(file).readlines(  ))
        if(count > 1):
            data = pd.read_csv(file, sep='\t', header=0, keep_default_na='True')
        else:
            if(verbose):
                print "File %s doesnot have any structural variants to annotate." %(file)
            data = pd.read_csv(file, sep='\t', header=0, keep_default_na='True')
            outFileTxt = outFilePrefix + ".txt"
            outFileExl = outFilePrefix + ".xlsx"
            outFileJson = outFilePrefix + ".json"
            # Print to TSV file
            data.to_csv(outFileTxt, sep='\t', index=False)
            # Print to Excel
            data.to_excel(outFileExl, sheet_name='Annotated_SVs', index=False)
            # Print to Json
            data.to_json(outFileJson)
            sys.exit()        
    return (data)
    
# Gives elements at particular index in list
getVar = lambda searchList, ind: [searchList[i] for i in ind]
        
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
