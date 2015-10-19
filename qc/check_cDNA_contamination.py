"""
Created on 11/01/2015.
@author: Ronak H Shah
###Required Input columns in the structural variant file
TumorId    NormalId    Chr1    Pos1    Chr2    Pos2    SV_Type    Gene1    Gene2    Transcript1    Transcript2    Site1Description    Site2Description    Fusion    Confidence    Comments    Connection_Type    SV_LENGTH    MAPQ    PairEndReadSupport    SplitReadSupport    BrkptType    ConsensusSequence    TumorVariantCount    TumorSplitVariantCount    TumorReadCount    TumorGenotypeQScore    NormalVariantCount    NormalSplitVariantCount    NormalReadCount    NormalGenotypeQScore
###Output columns
"TumorId    Genes_with_cDNA_contamination"
"""

import argparse
import pandas as pd
import time
import sys

def main():
    parser = argparse.ArgumentParser(
        prog='check_cDNA_contamination.py',
        description='Calculate cDNA contamination per sample based of the Structural Variants Pipeline result',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]")
    parser.add_argument(
        "-s",
        "--svFile",
        action="store",
        dest="svFilename",
        required=True,
        metavar='SVfile.txt',
        help="Location of the structural variant file to be used")
    parser.add_argument(
        "-o",
        "--outputFileName",
        action="store",
        dest="outFileName",
        required=True,
        metavar='cDNA_contamination',
        help="Full path name for the output file")
    args = parser.parse_args()
    # Read the structural variant file
    (dataDF) = ReadSV(args.svFilename, args)
    # Process the structural variant file
    (dataDict) = ProcessData(dataDF)
    # Write results in the output file specified
    WriteResults(dataDict, args)
# Read the structural variant file
def ReadSV(file, args):
    dataDF = pd.read_csv(file, sep='\t', header=0, keep_default_na='True')
    return(dataDF)
# Process the pandas dataframe to get contaminataed genes
def ProcessData(dataDF):
    # Group the data by tumor id, gene1 name, type of sv
    gDF = dataDF.groupby(['TumorId', 'Gene1', 'SV_Type']).groups
    # initialize and empty dictionary
    resultdict = {}
    # traverse through gDF dictionary 
    for key, value in gDF.iteritems():
        # get name of the sample its meta data
        (tumorID, gene1, svtype) = key
        # check how many entries are there of sv type
        entires = len(value)
        # initialize the number of cDNA events and a list of corresponding genes
        count = 0
        geneList = []
        # run only if the event is deletion and has more then 2 entries for the same gene
        if(svtype == 'DEL' and entires >= 2):
            for idx in value:
                record = dataDF.loc[idx]
                site1 = str(record.loc['Site1Description'])
                site2 = str(record.loc['Site2Description'])
                fusion = str(record.loc['Fusion'])
                brkptType = str(record.loc['BrkptType'])
                # Skip entries that are within exon or are in-frame or out-of frame or are IMPRECISE.
                if(("Exon" in site1 and "Exon" in site2) or ("in frame" in fusion or "out of frame" in fusion) or (brkptType == "IMPRECISE")):
                    continue
                else:
                    count = count + 1
        # count the entries,genes and fill the resultdict
        if tumorID in resultdict:   
            geneList = resultdict[tumorID]
            if(count >= 2):
                geneList.append(gene1)
            resultdict[tumorID] = geneList
        else:
            if(count >= 2):
                geneList.append(gene1)
                resultdict[tumorID] = geneList
    return(resultdict)

def WriteResults(dataDict, args):
    if not dataDict:
        with open(args.outFileName,'w') as out:
            print>>out,"No contamination found"
        return()
    else:
        # Convert the dictionaries values which are as list to string
        for key, value in dataDict.iteritems():
            newvalue = ', '.join(value)
            dataDict[key] = newvalue
            
        resultDF = pd.DataFrame(dataDict.items(), columns=['TumorId', 'Genes_with_cDNA_contamination']).sort('TumorId')
        # Print to TSV file
        resultDF.to_csv(args.outFileName, sep='\t', index=False)
    return()
    
#Run the main program      
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    seconds = end_time - start_time
    print "Elapsed time was ", time.strftime("%H:%M:%S", time.gmtime(seconds))
