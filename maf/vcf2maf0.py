#!/opt/bin/python
""" convert vcf to pseudo-MAF format

"""
from __future__ import print_function
import sys
import re
import csv
##import vcf
from collections import defaultdict
import VcfParser
import argparse
from GeneralizedRead import GeneralizedRead

parser = argparse.ArgumentParser(description="VCF to MAF convertor")
parser.add_argument('-p','--pairing',action='store',
                    help="Specify a pairing file for T/N pairs")
parser.add_argument('-c','--caller',action='store',help="Specify a caller")
parser.add_argument('-i','--inputFile',action='store',help="The VCF file to be fransformed. Also used as a hint to find the .txt when using the mutect caller.")
parser.add_argument('-o','--outputFile',action='store',help='The output filename. If not specified, it is the VCF_NAME.maf')
parser.add_argument('-aF', '--additionalFile',action='store',help='If using mutect, this is the txt file. If omited, the program assumes the text file has a similar name to that maf')
parser.add_argument('-n','--normal',action='store',help='If a somatic caller produces VCFs with the column name NORMAL, renames it in the MAF with this argument')
parser.add_argument('-t','--tumor',action='store',help='If a somater callers produces VCFs with the column name TUMOR, renames it in the MAF with this argument')
parser.add_argument('-v','--verbose',action='store_true',help='If specified, create a verbose maf that includes vcf entries with no reads')
parser.add_argument('-T', '--tsvMode', action='store',help="If specified, make traditional maf0 files where indels are with vcf indexing and multiple alt alleles are returned")

callers = {"varscan": "VAR", "Varscan":"VAR", "var":"VAR", "UnifiedGenotyper":"UG", "unifiedgenotyper":"UG", "ug":"UG", "HaplotypeCaller":"HC", "haplotypecaller":"HC", "hap":"HC", "haplo":"HC", "Mutect":"MUT", "mutect":"MUT", "mut":"MUT", "SomaticSniper":"SS", "somaticsniper":"SS", "snip":"SS", "Strelka":"STR", "strelka":"STR", "str":"STR"}
somaticCallers = ["MUT",  "SS", "STR"] ##important for dealing with vcfs that just encode "NORMAL" and "TUMOR" instead of sample name
args = parser.parse_args()

if len(sys.argv) <= 1:
    print("""The program takes the following arguments:
-c, --caller: The caller used for the VCF. Currently this takes values varscan, unifiedgenotyper, haplotypecaller, mutect, somaticsniper.
-i, --inputFile: the vcf file to transform.
-o, --outputFile: the name of the output file. OPTIONAL. If not provided, it will output to VCF_NAME.maf
-p, --pairing: a tumor/normal pairing file.
-aF, --additionalFile: specify the mutect .txt file if the not default.
-n, --normal: the sample name of the NORMAL. Required for somatic callers. Only useful for certain somatic callers.
-t, --tumor: the sample name of the TUMOR. Required for somatic callers. Only useful for certain somatic callers.
-v, --verbose: If specified, create a verbose maf that includes vcf entries with no reads.
-T, --tsvMode: If specified, make traditional maf0 files where indels are with vcf indexing and multiple alt alleles are returned""", file=sys.stderr)
    sys.exit(64)
if callers.has_key(args.caller):
    1 
else:
    choices = callers.keys()
    choices.sort()
    choices = ",".join(choices) 
    print("The caller", args.caller, "is not a valid choice. Valid caller options are: ", choices, file=sys.stderr)
    sys.exit(64)

if args.inputFile:
    input=args.inputFile
else:
    print("This version of the caller requires an input file, provided with the -i flag.", file=sys.stderr)
    sys.exit(64)

if callers[args.caller] in somaticCallers and not (args.tumor and args.normal):
    print("For a somatic caller, you must give the tumor and normal with the -t and -n arguments", file=sys.stderr)
    sys.exit(64)

if args.outputFile:
    output=args.outputFile
else:
    arOutput = args.inputFile.split(".")
    arOutput[-1] = "maf"
    output = ".".join(arOutput)

pairMap=None

if args.pairing:
    pairMap=VcfParser.VcfParser.parsePairFile(args.pairing)
    if not pairMap:
      print("ERROR: Pairing file creates a null mapping")
      sys.exit(78)
tsvMode = False
if args.tsvMode:
    tsvMode = True

fo = sys.stdout

if output != "stdout":
    fo = open(output, 'w')
fi = open(input, 'r')
mafout = csv.writer(fo, delimiter="\t", lineterminator="\n")

header=["GENE","SAMPLE","CHROM","POS","REF","ALT","FILTER","QUAL","ID"]
sample_fields=["GT","GQ","ALT_FREQ","AD_REF","AD_ALT"]
portal_fields=["Mutation_Status","Validation_Status"]

parser = None


if callers[args.caller] == "MUT":
  parser = VcfParser.MutectParser(fi, pairMap, args.tumor, args.normal, tsvMode)
  mutectText = args.additionalFile
  if not mutectText:
    arInput = input.split(".")
    arInput[-1] = "txt"
    mutectText = ".".join(arInput) 
  parser.parseAdditionalFile(mutectText)
  parser.trimPairmap()
elif (callers[args.caller] == "VAR" and args.normal) or callers[args.caller] ==  "SS" or callers[args.caller] == "STR":
  pairMap = {"TUMOR": ["NORMAL"]}
  parser = VcfParser.SomaticParser(fi, pairMap, args.tumor, args.normal, args.verbose, tsvMode)
elif callers[args.caller] == "HC":
  parser =VcfParser.HapParser(fi, pairMap, args.tumor, args.normal, args.verbose, tsvMode)
else:
  parser = VcfParser.VcfParser(fi, pairMap, args.tumor, args.normal, args.verbose, tsvMode)

uniq_formats=parser.calculateCallerSpecificFormat(callers[args.caller])
uniq_infos=parser.calculateCallerSpecificInfos(callers[args.caller])
parser.setCaller(callers[args.caller])

sample_uniq_formats = []
if pairMap:
    sample_fields+=["NORM_SAMPLE","NORM_GT","NORM_GQ","NORM_ALT_FREQ","NORM_AD_REF","NORM_AD_ALT"]
    for f in uniq_formats:
        sample_uniq_formats.append("NORM_" + f)
sample_fields.append("CALLER")
mafout.writerow(header+sample_fields+uniq_infos+uniq_formats+sample_uniq_formats+portal_fields)

parser.parse(GeneralizedRead, fo)

fi.close()
if output != "stdout":
    fo.close()
