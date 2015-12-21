class GeneralizedRead:
  def __init__(self, gene,sample,chrom,pos,ref,alt,filter,qual, id, gt,gq,altFreq,adRef,adAlt, uniqInfos, uniqFormats, inTsvMode=False):
    filter = filter or []  
    self.gene = gene
    self.sample = sample
    self.chrom = chrom
    self.pos = str(pos)
    self.ref = ref
    self.alt = ";".join(map(str, alt))
    self.filter = ";".join(filter)
    self._qual = qual
    self._id = id
    self.gt = gt
    self.gtOrig = gt
    self.gq = str(gq)
    self.altFreq = altFreq
    self.adRef = adRef
    self.adAlt = adAlt
    self.uniqI = uniqInfos
    self.uniqF = uniqFormats
    self.caller = ""
    if self.filter == "":
       self.filter = "PASS"
    # if there is a multiallelic alt, use the one that corresponds to the gt. 
    if not inTsvMode and  len(alt) > 1:
       if self.gt == "0/0" or self.gt == "./.":
          self.alt = self.ref
       else:
          genoElems = self.gt.split("/")
          self.alt = str(alt[int(genoElems[1]) - 1])
          if(genoElems[0] == genoElems[1]):
            self.gt = "1/1"
          elif(genoElems[0] != "0"):
            self.gt = "1/1"
          else:
            self.gt = "0/1"
 ## map vcf indel information into maf info
    if not inTsvMode and len(ref) > len(self.alt): #deletion event, like CTG -> C, including cases that can occur with multiallelic events like CTGTGTG->CTG
       self.pos =  str(pos + len(self.alt)) #str(pos + 1)
       self.ref = ref[len(self.alt):]
       self.alt = "-"
    elif not inTsvMode and len(self.alt) > len(ref): #insertion event, like C -> CTG, including cases that can occur with multiallelic events like CTG->CTGTGTG
       self.pos = str(pos + len(ref) - 1)
       self.ref = "-"
       self.alt = self.alt[len(ref):]

  @property
  def qual(self):
    return (self._qual or ".")
 
  @qual.setter
  def qual(self, val):
    self._qual = (val or ".")

  @property
  def id(self):
    return (self._id or ".")

  @id.setter
  def id(self, val):
    self._id = (val or ".")

  def mainString(self):
    return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.gene, self.sample, self.chrom, self.pos, self.ref, self.alt, self.filter, self.qual, self.id, self.gt, self.gq, self.altFreq, self.adRef, self.adAlt)

  def __str__(self):
    s =  "%s\t%s\t%s\t%s\t.\t." % (self.mainString(), self.caller, "\t".join(self.uniqI), "\t".join(self.uniqF))   
    return s

  def updateInfo(self, infoMap):
    pass

  def addCaller(self, caller):
    self.caller = caller

class MutectRead(GeneralizedRead):
   
  def updateInfo(self, infoMap):
    key = self.chrom + "_" + self.pos
    info = infoMap[key]
    if(self.filter == "REJECT"):
      self.filter = info[0]

class GenericSomaticRead(GeneralizedRead):
  def updateInfo(self, infoMap):
    self.sample = infoMap[self.sample] #rename samples of the type "TUMOR" and "NORMAL" that Strelka et al. produces

class TumorNormalPair:
  def __init__(self, tumor, normal):
    self.tumor = tumor
    self.normal = normal
    for index in range(len(self.tumor.uniqF)):
      #assume any format that has more than four semicolons is a likelihood
      if(self.tumor.uniqF[index].count(";") > 4):
        self.tumor.uniqF[index] = self.filterLikelihoods( self.tumor.uniqF[index], self.tumor.gtOrig)
    for index in range(len(self.normal.uniqF)):
      #assume any format that has more than four semicolons is a likelihood
      if(self.normal.uniqF[index].count(";") > 4):
        self.normal.uniqF[index] = self.filterLikelihoods( self.normal.uniqF[index], self.tumor.gtOrig)

  def __str__(self):
    s = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t." % (self.tumor.mainString(), self.normal.sample, self.normal.gt, self.normal.gq, self.normal.altFreq, self.normal.adRef, self.normal.adAlt, self.normal.caller, "\t".join(self.tumor.uniqI), "\t".join(self.tumor.uniqF), "\t".join(self.normal.uniqF)) 
    return s

  def updateInfo(self, infoMap):
    self.tumor.updateInfo(infoMap)
    self.normal.updateInfo(infoMap)

  def addCaller(self, caller):
    self.normal.caller = caller
    self.tumor.caller= caller

#### FROM THE VCF SPECS
# If A is the allele in REF and B,C,... are the alleles as ordered in ALT, 
# the ordering of genotypes for the likelihoods is given by: 
#F(j/k) = (k*(k+1)/2)+j. In other words, for biallelic sites the 
#ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC
  def filterLikelihoods(self, likelihood, origGt):
    likelihoods = likelihood.split(";")
    gtElems = origGt.split("/")
    newLikelihood = likelihoods[0] + ";" + str(likelihoods[int(float(gtElems[1]) * (float(gtElems[1]) + 1.0)/2.0 + 0.0)]) + ";" + str(likelihoods[int(float(gtElems[1]) * (float(gtElems[1]) + 1.0)/2.0 + float(gtElems[1]))])
    return newLikelihood


