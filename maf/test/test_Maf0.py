#!/opt/common/CentOS_6/python/python-2.7.7/bin/python
""" Test the VCF to maf0 converter by whole file comparisons and by unit tests
"""

import sys
import re
import csv
import vcf
import VcfParser
from GeneralizedRead import GeneralizedRead
from collections import defaultdict
import subprocess
import unittest
import shlex

class TestMaf(unittest.TestCase):  
  def comparePairToMutectVcf(self, tumorNormal):
    mRead = tumorNormal.tumor
    nRead = tumorNormal.normal
    p = subprocess.Popen(['grep', '29419659', '/ifs/work/gabow/testing/maf/test_mutect1.vcf'], stdout=subprocess.PIPE)
    o, e = p.communicate()
    lineArray = o.strip().split("\t")
    gti = 0
    adi= 1
    dpi = 3
    mutInfo = lineArray[-2].split(":")
    normInfo = lineArray[-1].split(":")
    ngt = normInfo[gti]
    if ngt == "0":
      ngt = "0/0"
    self.assertEqual(nRead.gt, ngt)
    self.assertEqual(mRead.gt, mutInfo[gti])
    self.assertEqual(str(nRead.adRef), normInfo[adi].split(",")[0])
    self.assertEqual(str(mRead.adAlt), mutInfo[adi].split(",")[1]) 
    self.assertEqual(mRead.chrom, lineArray[0])
    self.assertEqual(mRead.pos, lineArray[1])
    self.assertEqual(mRead.id, lineArray[2])
    self.assertEqual(mRead.ref, lineArray[3])
    self.assertEqual(mRead.alt, lineArray[4])
    self.assertEqual(mRead.filter, lineArray[6])

  def checkMutectText(self, mRead):
    p = subprocess.Popen(['grep', '29448410', '/ifs/work/gabow/testing/maf/test_mutect1.txt'],stdout=subprocess.PIPE)
    o, e = p.communicate()
    lineArray = o.strip().split("\t")
    reason = lineArray[-2]
    self.assertEqual(mRead.filter, 'fstar_tumor_lod,possible_contamination')
    p = subprocess.Popen(['grep', '29121285', '/ifs/work/gabow/testing/maf/test_mutect1.txt'], stdout=subprocess.PIPE)
    o, e = p.communicate()
    lineArray = o.strip().split("\t")
    mapQ0Reads = lineArray[17]
    #self.assertEqual(mRead.uniq[map_q0_reads], mapQ0Reads)

  def test_pairFileParse(self):
    pairMap = VcfParser.VcfParser.parsePairFile("/ifs/work/gabow/testing/maf/test_mutect1_sample_pairing.txt")
    self.assertEqual(pairMap["s_DS_blarad_001_P"],  ["s_DS_blarad_001_N3"])
    self.assertEqual(pairMap["s_DS_blarad_002_P"],  ["s_DS_blarad_002_N"])
    self.assertEqual(pairMap["s_DS_blarad_003_P2"],  ["s_DS_blarad_003_N"])
    self.assertEqual(pairMap["s_DS_blarad_004_M2"],  ["s_DS_blarad_004_N"])
    self.assertEqual(pairMap["s_DS_blarad_004_P"],  ["s_DS_blarad_004_N"])
    self.assertEqual(pairMap["s_DS_blarad_005_M"],  ["s_DS_blarad_005_N"])
    self.assertEqual(pairMap["s_DS_blarad_005_P"],  ["s_DS_blarad_005_N"])
    self.assertEqual(pairMap["s_DS_blarad_006_P2"],  ["s_DS_blarad_006_N"])
    self.assertEqual(pairMap["s_DS_blarad_007_P2"],  ["s_FFPE_pool_Normal"])

  def test_parseMutectPairedRead(self):
    fo = open("/ifs/work/gabow/testing/maf/test_mutect1.vcf", 'r')
    junk = open('junk.txt', 'w')
    mp = VcfParser.MutectParser(fo,  {"s_DS_blarad_005_M":["s_DS_blarad_005_N"]}) 
    mp.parseAdditionalFile("/ifs/work/gabow/testing/maf/test_mutect1.txt")
    reads = mp.testParserHarness(GeneralizedRead, junk)
    fo.close()
    junk.close()
    i = 0
    r = reads[i]
    while(r.tumor.pos != '29419659'): 
      i += 1
      r = reads[i]
    self.comparePairToMutectVcf(r)

  def test_additionalInfo(self):
    fo = open("/ifs/work/gabow/testing/maf/test_mutect1.vcf", 'r')
    mp = VcfParser.MutectParser(fo,  {"s_DS_blarad_005_M":["s_DS_blarad_005_N"]}, "s_DS_blarad_005_M", "s_DS_blarad_005_N", False)
    mp.parseAdditionalFile("/ifs/work/gabow/testing/maf/test_mutect1.txt")
    uniq_formats=mp.calculateCallerSpecificFormat("MUT")
    uniq_infos=mp.calculateCallerSpecificInfos("MUT")
    reads = []
    rec =  mp.vcfin.next()
    positionSamples = {}
    uniqInfo = mp.populateUniqInfos(rec)
    self.assertEqual(len(uniqInfo), 7)
    self.assertEqual(uniqInfo[0], "True")
    self.assertEqual(uniqInfo[1], ".")
    self.assertEqual(uniqInfo[4], "UNCOVERED")
    self.assertEqual(uniqInfo[5], "REJECT")

    for sample in rec.samples:
      if "GT" in sample.data._fields and (not sample.data.GT or sample.data.GT=="./."):
        continue
      parseMap = mp.parseSample(rec, sample)
      uniqFormat = mp.populateUniqFormats(sample)

    fo.close()
    fo = open("/ifs/work/gabow/testing/maf/test_mutect1.vcf", 'r')
    junk = open('junk.txt', 'w')
    mp = VcfParser.MutectParser(fo,  {"s_DS_blarad_005_M":["s_DS_blarad_005_N"]}, "s_DS_blarad_005_M", "s_DS_blarad_005_N")
    mp.parseAdditionalFile("/ifs/work/gabow/testing/maf/test_mutect1.txt")
    uniq_formats=mp.calculateCallerSpecificFormat("MUT")
    uniq_infos=mp.calculateCallerSpecificInfos("MUT")
    reads = mp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    i = 0
    r = reads[i]
    while(r.tumor.pos != '29419659'):
      i += 1
      r = reads[i]
    fo.close()


#  @unittest.skip("for now")
  def test_system_mutect(self):
    cmd = "/opt/common/CentOS_6/python/python-2.7.7/bin/python vcf2maf0.py -c mutect -i /ifs/work/gabow/testing/maf/test_mutect1.vcf -o test_diff.maf0 -p /ifs/work/gabow/testing/maf/test_mutect1_sample_pairing.txt -t s_DS_blarad_005_M -n s_DS_blarad_005_N"
    p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    o, e = p.communicate()
    cmd = "diff test_diff.maf0 /ifs/work/gabow/testing/maf/test_maf0_snpeff_removed_caller_added_c.txt"
    p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    o, e = p.communicate()
    self.assertEqual(o, "")

  def test_parseStrelkaPairedRead(self):
    fo =open("/ifs/work/gabow/testing/maf/test_stelka.indels.vcf", 'r')
    sp = VcfParser.SomaticParser(fo,  {"TUMOR":["NORMAL"]}, "s_DS_blarad_005_M", "s_DS_blarad_005_N", False, True)
    uniq_formats=sp.calculateCallerSpecificFormat("STR")
    uniq_infos=sp.calculateCallerSpecificInfos("STR")
    junk = open('junk.txt', 'w')
    reads = sp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    fo.close()
    i = 0
    r = reads[i]
    cp = '32007827'
    while(r.tumor.pos != cp):
      i += 1
      r = reads[i]
    nRead = r.normal
    mRead = r.tumor
    #self.assertEqual(nRead.sample,"s_DS_blarad_005_N")
    #self.assertEqual(mRead.sample,"s_DS_blarad_005_M")
    self.assertEqual(nRead.sample,"NORMAL")
    self.assertEqual(mRead.sample,"TUMOR")
    self.assertEqual(mRead.chrom, "chr6")
    self.assertEqual(mRead.pos, "32007827")
    self.assertEqual(mRead.ref, "CT")
    self.assertEqual(mRead.alt, "C")
    self.assertEqual(mRead.filter, "PASS")

  def test_parseVarscanPairedRead(self):
    fo  = open("/ifs/work/gabow/testing/maf/test_varscan.snp.vcf", 'r')
    sp = VcfParser.SomaticParser(fo,  {"TUMOR":["NORMAL"]}, "s_DS_blarad_005_M", "s_DS_blarad_005_N", False, True)
    uniq_formats=sp.calculateCallerSpecificFormat("VAR")
    uniq_infos=sp.calculateCallerSpecificInfos("VAR")
    junk = open('junk.txt', 'w')
    reads = sp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    fo.close()
    i = 0
    r = reads[i]
    cp = '16350323' 
    while(r.tumor.pos != cp):
      i += 1
      r = reads[i]
    nRead = r.normal
    mRead = r.tumor
#chr22   16350323    .   T   C   .   PASS    DP=260;SS=1;SSC=106;GPV=1E0;SPV=2.1033E-11  GT:GQ:DP:RD:AD:FREQ:DP4 1/1:.:136:28:105:78.95%:3,25,33,72  0/1:.:124:76:47:38.21%:9,67,5,42
    self.assertEqual(nRead.gt, "1/1")
    self.assertEqual(mRead.gt, "0/1")
    self.assertEqual(nRead.gq, "NA")
    self.assertEqual(mRead.gq, "NA")
    self.assertEqual(nRead.uniqF[1], "78.95%")
    self.assertEqual(mRead.uniqF[1], "38.21%")
    self.assertEqual(mRead.uniqI[0], "260")
    self.assertEqual(mRead.uniqI[1], ".")
    self.assertEqual(mRead.uniqI[2], "1")
    self.assertEqual(mRead.uniqI[3], "106")
    self.assertEqual(mRead.uniqI[4], "1.0")
    self.assertEqual(mRead.uniqI[5], "2.1033e-11")
    self.assertEqual(mRead.chrom, "chr22")
    self.assertEqual(mRead.pos, cp)
    self.assertEqual(mRead.ref, "T")
    self.assertEqual(mRead.alt, "C")
    self.assertEqual(mRead.filter, "PASS")


  def test_parseSniperPairedRead(self):
    fo  = open("/ifs/work/gabow/testing/maf/test_sniper.vcf", 'r')
    sp = VcfParser.SomaticParser(fo,  {"TUMOR":["NORMAL"]}, "s_DS_blarad_005_M", "s_DS_blarad_005_N", False, True)
    uniq_formats=sp.calculateCallerSpecificFormat("SS")
    uniq_infos=sp.calculateCallerSpecificInfos("SS")
    junk = open('junk.txt', 'w')
    reads = sp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    fo.close()
    i = 0
    r = reads[i]
    cp = '14377474'
    while(r.tumor.pos != cp):
      i += 1
      r = reads[i]
    nRead = r.normal
    mRead = r.tumor
#chr12   14377474    .   g   A   .   .   .   GT:IGT:DP:DP4:BCOUNT:GQ:JGQ:VAQ:BQ:MQ:AMQ:SS:SSC    0/0:0/0:1:0,1,0,0:0,0,1,0:30:.:0:36:60:60:0:.   1/1:1/1:3:0,1,1,1:2,0,1,0:3:.:31:32:44:36:2:24
    self.assertEqual(nRead.gt, "0/0")
    self.assertEqual(mRead.gt, "1/1")
    self.assertEqual(nRead.gq, "30")
    self.assertEqual(mRead.gq, "3")
    self.assertEqual(nRead.uniqF[-1], ".")
    self.assertEqual(mRead.uniqF[-1], "24")
    self.assertEqual(mRead.chrom, "chr12")
    self.assertEqual(mRead.pos, cp)
    self.assertEqual(mRead.ref, "g")
    self.assertEqual(mRead.alt, "A")
    self.assertEqual(mRead.filter, "PASS")


  def test_parseHaploUnpairedRead(self):
    fo  = open("/ifs/work/gabow/testing/maf/test_haplotype.vcf", 'r')
    hp = VcfParser.HapParser(fo, verbose=True)
    uniq_formats=hp.calculateCallerSpecificFormat("HAP")
    uniq_infos=hp.calculateCallerSpecificInfos("HAP")
    junk = open('junk.txt', 'w')
    reads = hp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    fo.close()
    i = 0
    r = reads[i]
    cp = '20170150'
    while(r.pos != cp):
      i += 1
      r = reads[i]
    read = r

  def test_parseHaploPairedRead(self):
    pairMap = VcfParser.VcfParser.parsePairFile("/ifs/work/gabow/testing/maf/test_mutect1_sample_pairing.txt")
    fo  = open("/ifs/work/gabow/testing/maf/test_haplotype.vcf", 'r')
    hp = VcfParser.VcfParser(fo, pairMap)
    uniq_formats=hp.calculateCallerSpecificFormat("HAP")
    uniq_infos=hp.calculateCallerSpecificInfos("HAP")
    junk = open('junk.txt', 'w')
    reads = hp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    fo.close()
    for read in reads:
      self.assertEqual(read.tumor.id, read.normal.id)
    i = 0
    r = reads[i]
    cp = '20104688'
    while(r.tumor.pos != cp):
      i += 1
      r = reads[i]
    nRead = r.normal
    mRead = r.tumor
#s_DS_blarad_001_N3  s_DS_blarad_001_P   s_DS_blarad_002_N   s_DS_blarad_002_P   s_DS_blarad_003_N   s_DS_blarad_003_P2  s_DS_blarad_004_M2  s_DS_blarad_004_N   s_DS_blarad_004_P   s_DS_blarad_005_M   s_DS_blarad_005_N   s_DS_blarad_005_P   s_DS_blarad_006_N   s_DS_blarad_006_P2  s_DS_blarad_007_P2  s_FFPE_pool_Normal  s_Frozen_pool_normal
#chr15   20104688    rs7179076   G   A   123.58  PASS    AC=6;AF=0.600;AN=10;BaseQRankSum=-1.221;ClippingRankSum=0.322;DB;DP=8;FS=0.000;MLEAC=6;MLEAF=0.600;MQ=55.06;MQ0=0;MQRankSum=0.956;QD=24.72;ReadPosRankSum=-0.322;VQSLOD=5.30;culprit=FS GT:AD:GQ:PL ./. ./. ./. ./. ./. ./. 0/0:2,0:6:0,6,75    ./. ./. ./. ./. ./. 1/1:0,3:9:85,9,0    1/1:0,1:3:36,3,0    1/1:0,1:3:37,3,0    0/0:1,0:3:0,3,40    ./.
    self.assertEqual(mRead.chrom, "chr15")
    self.assertEqual(mRead.pos, cp)
    self.assertEqual(mRead.ref, "G")
    self.assertEqual(mRead.alt, "A")
    self.assertEqual(mRead.filter, "PASS")
    self.assertEqual(nRead.sample, "s_DS_blarad_006_N")
    self.assertEqual(mRead.sample, "s_DS_blarad_006_P2")
    self.assertEqual(nRead.gt, "1/1")
    self.assertEqual(mRead.gt, "1/1")
    self.assertEqual(nRead.gq, "9")
    self.assertEqual(mRead.gq, "3")
    self.assertEqual(nRead.uniqF[-1], "85;9;0")
    self.assertEqual(mRead.uniqF[-1], "36;3;0")
    r = reads[i+1]
    nRead = r.normal
    mRead = r.tumor
    self.assertEqual(nRead.sample, "s_FFPE_pool_Normal")
    self.assertEqual(mRead.sample, "s_DS_blarad_007_P2")
    self.assertEqual(nRead.gt, "0/0")
    self.assertEqual(mRead.gt, "1/1")
    self.assertEqual(nRead.gq, "3")
    self.assertEqual(mRead.gq, "3")
    self.assertEqual(nRead.uniqF[-1], "0;3;40")
    self.assertEqual(mRead.uniqF[-1], "37;3;0")
    
  def test_VarscanUnpaired(self):
    fo  = open("/ifs/work/gabow/testing/maf/test_varscan.snp.vcf", 'r')
    sp = VcfParser.VcfParser(fo)
    uniq_formats=sp.calculateCallerSpecificFormat("VAR")
    uniq_infos=sp.calculateCallerSpecificInfos("VAR")
    junk = open('junk.txt', 'w')
    reads = sp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    fo.close()
    i = 0
    r = reads[i]
    cp = '7023565'
    while(r.pos != cp):
      i += 1
      r = reads[i]
    read = r
#OLD chr18   60795860    .   A   G   .   PASS    ADP=94144;WT=0;HET=1;HOM=0;NC=0 GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    0/1:141:94144:94144:94117:16:0.02%:6.3462E-15:27:15:41930:52187:8:8
#chr18   7023565 .   A   G   .   PASS    DP=106;SS=1;SSC=2;GPV=3.5492E-11;SPV=5.7955E-1  GT:GQ:DP:RD:AD:FREQ:DP4 0/1:.:48:34:14:29.17%:22,12,8,6 0/1:.:58:41:17:29.31%:23,18,14,3

    self.assertEqual(read.gt, "0/1")
    self.assertEqual(read.gq, ".")
    self.assertEqual(read.uniqF[2], "0.02%")
    self.assertEqual(read.uniqI[0], "94144")
    self.assertEqual(read.uniqI[1], "0")
    self.assertEqual(read.uniqI[2], "1")
    self.assertEqual(read.uniqI[3], "0")
    self.assertEqual(read.uniqI[4], "0")

 
  def test_SomaticArguments(self):
    cmd = "/opt/common/CentOS_6/python/python-2.7.7/bin/python vcf2maf0.py -c mutect -i /ifs/work/gabow/testing/maf/test_varscan.snp.vcf -o test_fail.maf0 -p /ifs/work/gabow/testing/maf/test_mutect1_sample_pairing.txt -t foo -T"
    p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = p.communicate()
    self.assertEqual(e, "For a somatic caller, you must give the tumor and normal with the -t and -n arguments\n")


  def test_VcfStyleIndel(self):
    pairMap = VcfParser.VcfParser.parsePairFile("/ifs/work/gabow/testing/maf/test_haplo_multi_sample_pairing.txt")
    fo  = open( '/ifs/work/gabow/testing/maf/test_haplo_multi.vcf', 'r')
    sp = VcfParser.VcfParser(fo, pairMap, "DS_bla_188_N", "DS_bla_188_T1")
    uniq_formats=sp.calculateCallerSpecificFormat("VAR")
    uniq_infos=sp.calculateCallerSpecificInfos("VAR")
    junk = open('junk.txt', 'w')
    reads = sp.testParserHarness(GeneralizedRead, junk)
    junk.close()
    fo.close()
    i = 0
    r = reads[i]
    cp = '12252337'
    while(r.tumor.pos != cp):
      i += 1
      r = reads[i]
    #CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT DS_bla_188_N  DS_bla_188_T1
    #chr1 12252333  . CTGTG C,CTG,CTGTGTG 3583.85 VQSRTrancheINDEL99.00to99.90  AC=12,7,17;AF=0.176,0.103,0.250;AN=68;BaseQRankSum=1.980;ClippingRankSum=0.032;DP=494;FS=14.240;MLEAC=12,6,17;MLEAF=0.176,0.088,0.250;MQ=59.90;MQRankSum=-0.218;NEGATIVE_TRAIN_SITE;QD=7.77;ReadPosRankSum=-0.576;SOR=1.632;VQSLOD=-1.612e+00;culprit=FS  GT:AD:GQ:PL 0/3:2,0,0,1:16:16,40,304,40,249,238,0,126,123,111 0/3:8,0,1,8:99:202,254,738,221,604,576,0,247,192,203
    nRead = r.normal
    tRead = r.tumor
    self.assertEqual(tRead.chrom, "chr1")
    self.assertEqual(tRead.pos, cp)
    self.assertEqual(tRead.ref, "-")
    self.assertEqual(tRead.alt, "TG")
    self.assertEqual(nRead.gt, "0/1")
    self.assertEqual(tRead.gt, "0/1")
    #self.assertEqual(nRead.uniqF[-1], "85;9;0")
    #self.assertEqual(tRead.uniqF[-1], "36;3;0")

if __name__ ==  '__main__':
    unittest.main()
