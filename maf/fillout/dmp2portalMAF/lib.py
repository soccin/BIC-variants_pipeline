import TCGA_MAF

class Bunch:
    def __init__(self,dictRec):
        self.__dict__.update(dictRec)
    def __str__(self):
        return str(self.__dict__)

def bunchStream(cin):
    for rec in cin:
        yield Bunch(rec)

def getVarType(s):
    ref=s.Ref.replace('-','')
    alt=s.Alt.replace('-','')
    if len(ref)==len(alt):
        if len(ref)==1:
            return "SNP"
        elif len(ref)==2:
            return "DNP"
        elif len(ref)==3:
            return "TNP"
        else:
            return "ONP"
    elif len(ref)>len(alt):
        return "DEL"
    else:
        return "INS"

def getEventSig(maf):
    return (maf.Chromosome,maf.Start_Position,maf.Reference_Allele,maf.Tumor_Seq_Allele2)

def fillSampleMAFFields(maf,si,ei,eventDb,pairs,caller,species):
    maf.Tumor_Sample_Barcode=si
    if not pairs.startswith("REF"):
        matchedNormal=pairs
        maf.n_ref_count=eventDb[ei][matchedNormal]["RD"]
        maf.n_alt_count=eventDb[ei][matchedNormal]["AD"]

    else:
        matchedNormal="REF." + species
        maf.n_ref_count="NA"
        maf.n_ref_count="NA"
    maf.Matched_Norm_Sample_Barcode= matchedNormal
    maf.t_ref_count=eventDb[ei][si]["RD"]
    maf.t_alt_count=eventDb[ei][si]["AD"]
    maf.Caller=caller
    return maf

class TCGA_MAF_Ext(TCGA_MAF.TCGA_MAF):
    pass

