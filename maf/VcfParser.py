import sys
import re
import csv
import vcf
from GeneralizedRead import GeneralizedRead, MutectRead, GenericSomaticRead, TumorNormalPair
from collections import defaultdict

class VcfParser:
  def __init__(self, file, tumor2normals=None, tumor=None, normal=None, verbose=False):
    self.file = file
    self.pairMap = tumor2normals
    self.tumor = tumor
    self.normal = normal
    self.vcfin = vcf.Reader(file) 
    self.uis = []
    self.ufs = []
    self.extraInfo = {}
    self.caller = ""
    self.mode = None
    if(verbose):
      self.mode = "VERBOSE"
  
  @staticmethod
  def parsePairFile(file):
    pairMap=defaultdict(list)
    with open(file) as fp:
        for line in fp:
            F=line.strip().split()
            if len(F)==0:
                # Skip blank lines
                continue
            elif len(F)!=2:
                print >>sys.stderr, "Invalid Pairing File Line"
                print >>sys.stderr, "["+line[:-1]+"]"
                return None

            if re.match('NA', F[0], re.IGNORECASE) or re.match('NA', F[1], re.IGNORECASE):
                continue

            pairMap[F[1]].append(F[0])
    return pairMap

  def getDP(self, samp):
    # Make work for Haplotype caller which does not have DP in format
    # just sum AD
    if "DP" in samp.data._fields:
      return samp.data.DP
    elif "AD" in samp.data._fields:
      return sum(samp.data.AD)

  def getAD(self, sample, gtAlt):
   ##
   ## Need to deal with case of multiple variants at a locus (C G,T)
   ##
   # First just get reference allele depth, sample_AD is a temp helper variable
    sample_AD = []
    if "AD" in sample.data._fields and "RD" in sample.data._fields:
      sample_AD.append(sample.data.RD)
      sample_AD.append(sample.data.AD)
    elif "AD" in sample.data._fields:
      sample_AD=sample.data.AD
    elif  "DP4" in sample.data._fields:
      sample_AD.append(sample.data.DP4[0]+ sample.data.DP4[1])
      sample_AD.append(sample.data.DP4[2]+ sample.data.DP4[3])
    alleleDepth = [sample_AD[0]]
    if gtAlt == "0":
      variantDepth=sum(sample_AD[1:])
    else:
      variantDepth=sample_AD[int(gtAlt)]
    alleleDepth.append(variantDepth)
    nraf="%.4g" % (float(variantDepth)/float(self.getDP(sample)))
    return {"nraf":nraf, "ad":alleleDepth}


  def setCaller(self, caller):
    if caller == "MUT":
      self.caller = "mutect"
    elif caller == "VAR":
      self.caller = "varscan"
    elif caller == "UG":
      self.caller = "unifiedgenotyper"
    elif caller == "HC":
      self.caller = "haplotypecaller"
    elif caller == "SS":
      self.caller = "somaticsniper"
    elif caller == "STR":
      self.caller = "strelka"

  def calculateCallerSpecificInfos(self, caller):
    standard = []
    uniq_infos = []
    infos = self.vcfin.infos.keys()
    for info in infos:
      if info not in standard:
        uniq_infos.append(caller + "_" + info)
        self.uis.append(info)
    return uniq_infos

  def calculateCallerSpecificFormat(self, caller):
    standard = ['GT', 'AD', 'DP', 'GQ']
    uniq_formats = []
    formats = self.vcfin.formats.keys()
    for format in formats:
      if format not in standard:
        uniq_formats.append(caller + "_" + format)
        self.ufs.append(format)
    return uniq_formats

  def trimPairmap(self):
    if(self.pairMap and self.tumor and self.normal):
      if(self.pairMap[self.tumor] and self.normal in self.pairMap[self.tumor]):
        self.pairMap = {self.tumor : [self.normal]}

        

  def getFormats(self):
    return self.vcfin.formats


  def parseSample(self, rec, sample):
    gt = gq = nraf = "NA"
    if "GT" in sample.data._fields:
      gt = sample.data.GT
      if sample.data.GT == "0":    #get rid of encoding "return 00/0" as "0"
        gt = "0/0"
    if "GQ" in sample.data._fields and sample.data.GQ:
      gq = sample.data.GQ
    if "AD" in sample.data._fields and sample.data.AD and self.getDP(sample)>0:
      ##
      ## Need to deal with case of multiple variants at a locus (C G,T)
      ##
      # First just get reference allele depth, sample_AD is a temp helper variable
      resmap = self.getAD(sample, gt[-1])
      alleleDepth = resmap["ad"]
      nraf = resmap["nraf"]
    else:
      alleleDepth = ["0","0"]
    resmap = {"gene":".", "sample":sample.sample, "chrom":rec.CHROM, "pos":rec.POS, "ref":rec.REF, "alt":rec.ALT, "filter":rec.FILTER, "qual":rec.QUAL, "id":rec.ID, "gt":gt, "gq":gq, "nraf":nraf, "adRef":alleleDepth[0], "adAlt": alleleDepth[1]}
    return resmap

  def populateUniqInfos(self, rec):
    uniques = []
    for uniq_info in self.uis:
      if uniq_info in rec.INFO:
        if type(rec.INFO[uniq_info]) == list:
          uniques.append(";".join(map(str, rec.INFO[uniq_info])))
        else:
          uniques.append(rec.INFO[uniq_info])
      else:
        uniques.append(".")
    return map(str, uniques)

  def populateUniqFormats(self, sample):
    uniques = []
    for uniq_format in self.ufs:
      if uniq_format in sample.data._fields and getattr(sample.data, uniq_format):
        if type(getattr(sample.data, uniq_format)) == list:
          flatFormat = ";".join(map(str, getattr(sample.data, uniq_format)))
          uniques.append(flatFormat)
        else:
          uniques.append(getattr(sample.data, uniq_format))
      else:
        uniques.append(".")
    return map(str, uniques)


  def testParserHarness(self, readClass, fileOut):
    reads = []

    for rec in self.vcfin:
        reads = reads + (self.parseRecord(readClass, rec))
        
    return reads

  def parse(self, readClass, fileOut):
    reads = []
    for rec in self.vcfin:
        reads = self.parseRecord(readClass, rec)
        for read in reads:
           fileOut.write(str(read) + "\n")

  def parseRecord(self, readClass, rec):
    reads = []
    positionSamples = {}
    uniqInfo = self.populateUniqInfos(rec)
    for sample in rec.samples:
       ##assert isinstance(sample,vcf._Call)
      if self.mode == "VERBOSE" and "GT" in sample.data._fields and (sample.data.GT=="./." or sample.data.GT=="0/0"):
        adRef = 0
        if "AD" in sample.data._fields and isinstance(sample.data.AD, int):
          adRef = sample.data.AD
        alt = rec.ALT
        if not alt[0]:
          alt = '.'
        positionSamples[sample.sample] = readClass(".", sample.sample, rec.CHROM, rec.POS, rec.REF, alt, rec.FILTER, rec.QUAL, rec.ID, "0/0", 0, 0, adRef, 0, uniqInfo, "")
        continue
      if "GT" in sample.data._fields and (not sample.data.GT or sample.data.GT=="./."):
        positionSamples[sample.sample] = None
        continue
      parseMap = self.parseSample(rec, sample)
      uniqFormat = self.populateUniqFormats(sample)
      positionSamples[sample.sample] =  readClass(parseMap["gene"], parseMap["sample"], parseMap["chrom"], parseMap["pos"], parseMap["ref"], parseMap["alt"], parseMap["filter"], parseMap["qual"], parseMap["id"], parseMap["gt"], parseMap["gq"], parseMap["nraf"], parseMap["adRef"], parseMap["adAlt"], uniqInfo, uniqFormat)
    if self.pairMap:
      for tumor, normals in self.pairMap.iteritems():
        for normal in normals:
          if(positionSamples[tumor] and positionSamples[normal]):
            read = TumorNormalPair(positionSamples[tumor], positionSamples[normal])
            read.updateInfo(self.extraInfo)
            read.addCaller(self.caller)
            reads.append(read)
    else:
      for k, v in positionSamples.iteritems():
        if v:
          v.updateInfo(self.extraInfo)
          v.addCaller(self.caller)
          reads.append(v)
    return reads

class SomaticParser(VcfParser): #generic somatic parser that doesn't use proper sample names but instead uses thing like "NORMAL" and "TUMOR"
  def __init__(self, file, tumor2normals, tumor, normal, verbose):
    VcfParser.__init__(self, file, tumor2normals, tumor, normal, verbose)
    self.extraInfo = {'TUMOR' : tumor, "NORMAL" : normal }

  def parse(self, generalizedRead, fileOut):
      return VcfParser.parse(self, GenericSomaticRead, fileOut)


class HapParser(VcfParser):
  #right now, we are currently having Haplotype always ignore DP as we don't understand it in 3.x gVCF
  def getDP(self, samp):
    return sum(samp.data.AD)

class MutectParser(VcfParser):
  def populateUniqInfos(self, rec):
    uniques = VcfParser.populateUniqInfos(self, rec)
    uniques.append(self.extraInfo[rec.CHROM + "_" + str(rec.POS)][1])
    uniques.append(self.extraInfo[rec.CHROM + "_" + str(rec.POS)][2])
    uniques.append(self.extraInfo[rec.CHROM + "_" + str(rec.POS)][3])
    return map(str, uniques)


  def calculateCallerSpecificInfos(self, caller):
    standard = []
    uniq_infos = []
    infos = self.vcfin.infos.keys()
    for info in infos:
      if info not in standard:
        uniq_infos.append(caller + "_" + info)
        self.uis.append(info)
    uniq_infos.append("MUT_COVERED")
    uniq_infos.append("MUT_KEEP")
    uniq_infos.append("MUT_MQ0_READS")
    return uniq_infos

  def parseAdditionalFile(self, mutectText):
    mutFilter = -1
    mutKeep = -1
    mutCover = -1
    mq0 = -1
    with open(mutectText) as mt:
      mt.readline()
      header = mt.readline().strip().split("\t")
      mutFilter = header.index("failure_reasons")
      mutCover = header.index("covered")
      mutKeep = header.index("judgement")
      mq0 = header.index("map_Q0_reads")
      for line in mt:
        if line[0] == "#":
          continue
        F = line.strip().split("\t")
        key = F[0] + "_" + F[1]
        self.extraInfo[key] = (F[mutFilter], F[mutCover], F[mutKeep], str(F[mq0]))

  def parse(self, generalizedRead, fileOut):
    return VcfParser.parse(self, MutectRead, fileOut)
