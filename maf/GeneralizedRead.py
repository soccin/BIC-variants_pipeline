class GeneralizedRead:
  def __init__(self, gene,sample,chrom,pos,ref,alt,filter,qual, id, gt,gq,altFreq,adRef,adAlt, uniqInfos, uniqFormats):
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
    self.gq = str(gq)
    self.altFreq = altFreq
    self.adRef = adRef
    self.adAlt = adAlt
    self.uniqI = uniqInfos
    self.uniqF = uniqFormats
    self.caller = ""
    if self.filter == "":
       self.filter = "PASS"

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

  def __str__(self):
    s = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t." % (self.tumor.mainString(), self.normal.sample, self.normal.gt, self.normal.gq, self.normal.altFreq, self.normal.adRef, self.normal.adAlt, self.normal.caller, "\t".join(self.tumor.uniqI), "\t".join(self.tumor.uniqF), "\t".join(self.normal.uniqF)) 
    return s

  def updateInfo(self, infoMap):
    self.tumor.updateInfo(infoMap)
    self.normal.updateInfo(infoMap)

  def addCaller(self, caller):
    self.normal.caller = caller
    self.tumor.caller= caller
