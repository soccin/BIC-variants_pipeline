library(ggplot2)
library(gplots)
library(scales)
library(reshape)
library(plyr)


get.latest.qc.review <- function(proj.id){
  projdir = paste("/home/shiny/CMO/projects/Proj_",proj.id,sep="")
  file = paste(projdir,"/Proj_",proj.id,"_QC_signoff.txt",sep="")
  dir.create(projdir,showWarnings=FALSE)
  if(!file.exists(file)){ return(NULL) }
  qc.reviews = read.delim(file)
  lastDate = as.vector(qc.reviews$Date[length(qc.reviews$Date)])[1]
  lastReviewer = as.vector(qc.reviews$User[length(qc.reviews$User)])[1]
  lastReview = qc.reviews[which(qc.reviews$Date == lastDate),]
  lastFinalComment = as.vector(qc.reviews$ProjComment[length(qc.reviews$ProjComment)])[1]

  review = list(date=lastDate,user=lastReviewer,rev=lastReview,projComment=lastFinalComment)
  return(review)
}

get.title <- function(path,type){

  if(type == 'exome'){
      return(NULL)
  } else {
      path = unlist(strsplit(path,"/"))
      path = paste(path[-c(length(path),length(path)-1)],collapse="/",sep="")
      filename = dir(path)[grep("_title.txt",dir(path))]
      if(length(filename)==0){ return(NULL) }

      file = paste(path,filename,sep="/")
      if(file.exists(file)){ 
          title = read.delim(file)
          return(title)
      }
      return(NULL)
  }
}

get.num.samples <- function(path,type){
    n = 0
    if(type == 'exome'){
        path = unlist(strsplit(path,"/"))
        path = paste(c(path[-length(path)],"alignments"),collapse="/",sep="")
        if(nchar(path)>1){
            n = length(grep(".bam$",dir(path)))
        }
    } else if(type == 'targeted' & !is.null(get.title(path,type))){
        n = nrow(get.title(path,type))
    }
    return(n)
}

get.is.metrics <- function(path,type){
    if(type == 'exome'){
        filename = dir(path)[grep("_InsertSizeMetrics_Histograms.txt",dir(path))]
    } else {
        filename = dir(path)[grep("_ALL_insertsizemetrics.txt",dir(path))]
    }
    if(length(filename) == 0) { return(NULL) }

    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

    is <- read.delim(file,header=T,sep="\t")
    ## normalize header
    colnames(is)[1] = "insert_size"    
    if(type == 'targeted'){
        title = get.title(path,type)
        colnames(is)[2:length(colnames(is))] = paste(title$Barcode,title$Sample_ID,sep=": ")
    }
    ### dmp IS file has one row for every insert size 1-500, plus a line for peaks
    ### exome IS file has one row for every insert size existing in data
    if(nrow(is)>500){ is = is[1:500,] }

    return(is)
}

get.is.peaks <- function(path,type){
    is = get.is.metrics(path,type)
    if(is.null(is)){ return(NULL) }
    peaks = as.data.frame(apply(is[,2:ncol(is)],2,which.max))
    p = data.frame(Sample=rownames(peaks),Peak=peaks[,1])
    return(p)
}

get.hs.metrics <- function(path,type){
    if(type == 'exome'){
        filename = dir(path)[grep("HsMetrics.txt",dir(path))]
    } else {
        filename = dir(path)[grep("_ALL_HSmetrics.txt",dir(path))]
    }
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }
    hs = read.delim(file,header=T,sep="\t")
    return(hs)
}

get.md.metrics <- function(path,type){
    if(type == 'targeted'){ return(NULL) }
    filename = dir(path)[grep("_markDuplicatesMetrics.txt",dir(path))]
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }
    dat <- read.delim(file)
    return(dat)
}

get.fpc.sum <- function(path,type){
    if(type == 'exome'){
        path=paste(path,"fingerprint",sep="/")
        filename = dir(path)[grep("_DiscordantHomAlleleFractions.txt",dir(path))]
        if(length(filename)==0) { return(NULL) }
        file = paste(path,filename,sep="/")
        if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }
        
        dat = read.delim(file)
        rownames(dat) = dat[,1]
        dat = data.matrix(dat[,-1])
        return(dat)
    }
    title = get.title(path,type)
    if(is.null(title)) { return(NULL) }

    filename = dir(path)[grep("_ALL_FPCsummary.txt",dir(path))]
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

    xt = read.table(file, as.is=TRUE, header=TRUE, sep="\t", row.names=1)
    Period = paste(title$Barcode,title$Sample_ID,sep=": ")
    colnames(xt) = Period
    rownames(xt) = Period

    xm = as.matrix(xt)
    if(isSymmetric(xm)==TRUE){ return(NULL) }
    return(xm)
}

get.unexpected.matches <- function(path,type){
    if(type == 'exome'){ return(NULL) }
    title = get.title(path,type)
    if(is.null(title)) { return(NULL) }
    
    filename = dir(path)[grep("_ALL_FPCResultsUnMatch.txt",dir(path))]
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

    match<-read.delim(file, header=T, comment.char="#")
    if(nrow(match)>0){
        ## ignore matches with pools
        pt <- title[title$Class == "PoolTumor",]["Sample_ID"]
        pn <- title[title$Class == "PoolNormal",]["Sample_ID"]
        ntc <- title[title$Class == "NTC",]["Sample_ID"]
        mn <- title[grep("Normal_Pooled",title$Sample_ID),]["Sample_ID"]
        ignore <- rbind(pt,pn,ntc,mn)
        for (norm in ignore$Sample_ID){
            if(nrow(match)>0){
                if(length(grep(norm,match[,1]))>0){
                    match <- match[-grep(norm,match[,1]),]
                }
                if(length(grep(norm,match[,2]))>0){
                    match <- match[-grep(norm,match[,2]),]
                }
            }
        }
    }
    return(match)
}

get.unexpected.mismatches <- function(path,type){
    if(type == 'exome'){ return(NULL) }
    title = get.title(path,type)
    if(is.null(title)) { return(NULL) }
    
    filename = dir(path)[grep("_ALL_FPCResultsUnMismatch.txt",dir(path))]
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

    mismatch<-read.delim(file, header=T, comment.char="#")
    if(nrow(mismatch)>0){
        ## ignore mismatches with pools
        pt <- title[title$Class == "PoolTumor",]["Sample_ID"]
        pn <- title[title$Class == "PoolNormal",]["Sample_ID"]
        ntc <- title[title$Class == "NTC",]["Sample_ID"]
        mn <- title[grep("Normal_Pooled",title$Sample_ID),]["Sample_ID"]
        ignore <- rbind(pt,pn,ntc,mn)
        for (norm in ignore$Sample_ID){
            if(nrow(mismatch)>0){
                if(length(grep(norm,mismatch[,1]))>0){
                    mismatch <- mismatch[-grep(norm,mismatch[,1]),]
                }
                if(length(grep(norm,mismatch[,2]))>0){
                    mismatch <- mismatch[-grep(norm,mismatch[,2]),]
                }
            }
        }
    }
    return(mismatch)
}

get.major.contamination <- function(path,type){
    if(type == 'exome'){
        path = paste(path,"fingerprint",sep="/")  
        filename = dir(path)[grep("_MajorContamination.txt",dir(path))]
        if(length(filename)==0) { return(NULL) }
        file = paste(path,filename,sep="/")
        if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

        dat = read.delim(file)
        return(dat)
    }    
    title = get.title(path,type)
    if(is.null(title)) { return(NULL) }
    
    filename = dir(path)[grep("_ALL_FPhet.txt",dir(path))]
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

    dat = read.delim(file)
    dat$Sample = paste(dat$Sample, title$Sample_ID, sep=": ")

    pt <- title[title$Class == "PoolTumor",]["Sample_ID"]
    pn <- title[title$Class == "PoolNormal",]["Sample_ID"]
    ntc <- title[title$Class == "NTC",]["Sample_ID"]
    mn <- title[grep("Normal_Pooled",title$Sample_ID),]["Sample_ID"]
    ignore = rbind(pt,pn,ntc,mn)

    for(pool in ignore$Sample_ID){
        if(length(grep(pool,dat$Sample))>0){
            dat = dat[-grep(pool,dat$Sample),]
        }
    }
    return(dat)
}

get.minor.contamination <- function(path,type){
    if(type == 'exome'){ 
        path = paste(path,"fingerprint",sep="/")
        filename = dir(path)[grep("_MinorContamination.txt",dir(path))]
        if(length(filename)==0) { return(NULL) }
        file = paste(path,filename,sep="/")
        if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

        dat = read.delim(file)
        return(dat)
    }
    title = get.title(path,type)
    if(is.null(title)) { return(NULL) }
    
    filename = dir(path)[grep("_ALL_FPavgHom.txt",dir(path))]
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

    dat = read.delim(file)
    dat$Sample = paste(dat$Sample, title$Sample_ID, sep=": ")

    pt <- title[title$Class == "PoolTumor",]["Sample_ID"]
    pn <- title[title$Class == "PoolNormal",]["Sample_ID"]
    ntc <- title[title$Class == "NTC",]["Sample_ID"]
    mn <- title[grep("Normal_Pooled",title$Sample_ID),]["Sample_ID"]
    ignore = rbind(pt,pn,ntc,mn)

    for(pool in ignore$Sample_ID){
        if(length(grep(pool,dat$Sample))>0){
            dat = dat[-grep(pool,dat$Sample),]
        }
    }
    return(dat)
}


get.base.qualities <- function(path,type){
    if(type == 'exome'){
        filePreName = dir(path)[grep("_pre_recal_MeanQualityByCycle.txt",dir(path))]
        filePostName = dir(path)[grep("_post_recal_MeanQualityByCycle.txt",dir(path))]
    } else {
        filePreName = dir(path)[grep("_ALL_orgbasequalities.txt",dir(path))]
        filePostName = dir(path)[grep("_ALL_basequalities.txt",dir(path))]
    }
    if(length(filePreName)==0 || length(filePostName)==0){ return(NULL) }
    filePre = paste(path,filePreName,sep="/")
    filePost = paste(path,filePostName,sep="/")
    if(!file.exists(filePre) || file.info(filePre)$size == 0){ return(NULL) }
    if(!file.exists(filePost) || file.info(filePost)$size == 0){ return(NULL) }

    ## get pre-recalibration data
    qualPre = read.delim(filePre)[1:200,]
    samples = colnames(qualPre)[-1]
    colnames(qualPre) = c("cycle",samples)
    qualPre.m = melt(qualPre,id.vars="cycle")
    qualPre.m$type = "Pre-Recalibration"
    ## get post-recalibration data
    qualPost = read.delim(filePost)[1:200,]
    colnames(qualPost) = c("cycle",samples)
    qualPost.m = melt(qualPost,id.vars="cycle")
    qualPost.m$type = "Post-Recalibration"
    
    qual.m = rbind(qualPre.m,qualPost.m)
    qual.m$type = factor(qual.m$type, levels = c("Pre-Recalibration","Post-Recalibration"))

    return(qual.m)
}

load.metrics.db <- function(){
    paths = as.matrix(read.delim("allProjects.txt",header=F))

    projects = list()
    types = c('exome','targeted')   

    for(path in paths){
        x = unlist(strsplit(path,"/"))[-1]
        proj.idx = grep("proj",x,ignore.case=TRUE)
        proj = sub("Proj_","",x[proj.idx])
        inv = x[proj.idx-1]
        pi = x[proj.idx-2]
        if(x[proj.idx+1] == "metrics"){
            if(length(grep("rna",dir(path),ignore.case=TRUE))==0){
                type = "exome"
            } else {
                type = 'rnaseq'
            } 
        } else { ## right now there are only two 'types'
            type = "targeted"
        }
        if(type %in% types){
            indexStr = paste(proj,"(",pi,":",inv,":",type,")",sep="  ")
            projects[[indexStr]] = list('pi'=pi,'inv'=inv,'id'=proj,'type'=type,'path'=path)
        }
    }
    return(projects)
}

get.coverage <- function(path,type){
    hs <- get.hs.metrics(path,type)
    if(is.null(hs)){ return(NULL) }
    if(type == 'exome'){ 
        coverage <- data.frame(Samples=hs$SAMPLE,Cov=hs$MEAN_TARGET_COVERAGE)
        return(coverage)
    }
    if(type == 'targeted'){
        filename = dir(path)[grep("_ALL_Canonical_exoncoverage.txt",dir(path))]
    }
    if(length(filename)==0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }

    title <- get.title(path,type)
    coverage <- read.delim(file)
    var <- colnames(coverage)[2:ncol(coverage)]
    a <- colwise(median)(coverage[, c(3:ncol(coverage))])
    coverage <- melt(a)
    coverage$variable <- paste(hs$Sample, title$Sample_ID, sep = ": ")
    colnames(coverage) <- c("Samples", "Cov")
    return(coverage)
}

get.duplication <- function(path,type){
    if(type == 'targeted'){
        dat <- get.hs.metrics(path,type)
        if(is.null(dat)){ return(NULL) }
        title <- get.title(path,type)
        if(!is.null(title)){
            sample.labels = paste(dat$Sample,title$Sample_ID,sep=": ")
        } else {
            return(NULL)
        }
    } else {
        dat = get.md.metrics(path,type)
        if(is.null(dat)){ return(NULL) }
        sample.labels = dat$LIBRARY
    }
    duplication <- data.frame(Samples = sample.labels,DupRate = dat$PERCENT_DUPLICATION)
    return(duplication)
}

get.library.size <- function(path,type){
    if(type == 'targeted'){
        title = get.title(path,type)
        hs = get.hs.metrics(path,type)
        if(is.null(hs) || is.null(title)){ return(NULL) }
        librarySize = data.frame(Samples = paste(hs$Sample,title$Sample_ID,sep = ": "),Comp = hs$ESTIMATED_LIBRARY_SIZE)
    } else {
        dat = get.md.metrics(path,type)
        if(is.null(dat)){ return(NULL) }
        librarySize = data.frame(Samples = dat$LIBRARY, Comp = dat$ESTIMATED_LIBRARY_SIZE)
    }
    return(librarySize)
}

get.capture.specificity <- function(path,type){
    hs = get.hs.metrics(path,type)
    if(is.null(hs)){ return(NULL) }
    if(type == 'exome'){ 
        labels = hs$SAMPLE
    } else {
        title = get.title(path,type)
        labels = paste(hs$Sample, title$Sample_ID, sep = ": ")
    }
    dat <- data.frame(Sample = labels, OnBait = hs$ON_BAIT_BASES, NearBait = hs$NEAR_BAIT_BASES, OffBait = hs$OFF_BAIT_BASES)
    return(dat)
}

get.alignment <- function(path,type){
    dat=NULL
    if(type == 'exome'){
       md = get.md.metrics(path,type)
       if(is.null(md)){ return(NULL) }
       n = (md$UNMAPPED_READS - md$UNPAIRED_READS_EXAMINED)/2
       n[which(n<0)] = 0
       dat = data.frame(Samples=md$LIBRARY,BothAlign=md$READ_PAIRS_EXAMINED,OneAlign=md$UNPAIRED_READS_EXAMINED,NeitherAlign=n)
    } else {
       hs = get.hs.metrics(path,type)
       title = get.title(path,type)
       if(is.null(hs) || is.null(title)){ return(NULL) }
       labels = paste(title$Barcode,title$Sample_ID,sep=": ")
       dat = data.frame(Samples = labels, BothAlign=hs$both.reads.align, OneAlign=hs$one.read.aligns, NeitherAlign=hs$neither.read.aligns)
    }
    return(dat)
}

get.trimmed.reads <- function(path,type){
    reads = NULL
    if(type == 'targeted'){
        title = get.title(path,type)
        hs = get.hs.metrics(path,type)
        if(is.null(title) || is.null(hs)){ return(NULL) }
        labels = paste(title$Barcode,title$Sample_ID,sep=": ")
        reads <- data.frame(Samples = labels, Read1=hs$PerRead1Trimmed, Read2=hs$PerRead2Trimmed)
    } else {
        filename = dir(path)[grep("_CutAdaptStats.txt",dir(path))]
        if(length(filename)==0) { return(NULL) }
        file = paste(path,filename,sep="/")
        if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }
        stats <- read.delim(file)    
        reads <- data.frame(Samples = stats$Sample, Read1=stats$R1_PercentTrimmed, Read2=stats$R2_PercentTrimmed)
    }
    return(reads)    
}

get.pool.norm.genotype <- function(path,type){
    dat = NULL
    if(type == 'exome'){ return(NULL) }
    filename = dir(path)[grep("_ALL_GenotypePooledNormal.txt",dir(path))]
    if(length(filename) == 0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }
    dat <- read.delim(file,header=T,sep="\t")    
    return(dat)
}

get.gc.bias <- function(path,type){
    dat = NULL
    if(type == 'exome'){ return(NULL) }
    filename = dir(path)[grep("_ALL_gcbias.txt",dir(path))]
    if(length(filename) == 0) { return(NULL) }
    file = paste(path,filename,sep="/")
    if(!file.exists(file) || file.info(file)$size == 0){ return(NULL) }
    dat <- read.delim(file,header=T,sep="\t")
    title = get.title(path,type)
    if(is.null(title)){ return(NULL) }
    colnames(dat) = c("X",paste(title$Barcode,title$Sample_ID, sep=": "))
    return(dat)
}

get.detail.table <- function(path,type,id,user=NULL){
    mets = c("Unexpected Match(es)",
             "Unexpected Mismatch(es)",
             "Major Contamination",
             "Minor Contamination",
             "Coverage",
             "Duplication",
             "Library Size (millions)","On Bait Bases (millions)","Aligned Reads (millions)","Insert Size Peak","Percentage Trimmed Reads")

    other = c(" ","Sample")
    detail = matrix(NA,nrow=get.num.samples(path,type),ncol=length(other)+length(mets))

    colnames(detail) = c(other,mets)

    if(type == 'targeted'){
        title = get.title(path,type)
        samples = paste(title$Barcode,title$Sample_ID,sep=": ")
    } else {
        samples = as.vector(get.md.metrics(path,type)$LIBRARY)
        ## this is a hacky way of normalizing sample names, since some
        ## stats are run before library "_1" is appended to sample ID
        ## and some are run afterwards in the pipeline
        samples = sub("__1","",samples)                
    }

    unm = get.unexpected.matches(path,type)
    unmm = get.unexpected.mismatches(path,type)
    majc = get.major.contamination(path,type)
    minc = get.minor.contamination(path,type)
    cov = get.coverage(path,type)
    dup = get.duplication(path,type)
    lib = get.library.size(path,type)
    cs = get.capture.specificity(path,type)
    aln = get.alignment(path,type)
    is = get.is.peaks(path,type)
    tr = get.trimmed.reads(path,type)

    ## populate matrix
    detail[," "] = "2PASS"
    detail[,"Sample"] = samples
    detail[,"Unexpected Match(es)"] = "None"
    detail[,"Unexpected Mismatch(es)"] = "None"
    rownames(detail) = samples
    if(!is.null(unm) && nrow(unm)>0){ 
      for(x in c(1,2)){
        for(i in 1:length(unm[,x])){
            samp = unm[i,x]
            if( x == 1 ){ 
              col=2 
            } else { 
              col=1 
            }
            unms = paste(unm[which(unm[,x]==samp),col],collapse="<br>",sep="")
            if(type == 'targeted'){
                parts = unlist(strsplit(as.character(samp),"_")) 
                bc = parts[length(parts)]
                finalSamp = paste(bc,paste(parts[-length(parts)],collapse="_",sep=""),sep=": ")
            } else {
                finalSamp = sub("__1","",samp)
            }
            detail[which(rownames(detail) == finalSamp),"Unexpected Match(es)"] = unms
            detail[which(rownames(detail) == finalSamp)," "] = "0FAIL"
        }
      }
    }
    if(!is.null(unmm) && nrow(unmm)>0){
      for(x in c(1,2)){
        for(i in 1:length(unmm[,x])){
            samp = unmm[i,x]
            if( x == 1 ){
              col=2
            } else {
              col=1
            }
            unmms = paste(unmm[which(unmm[,x]==samp),col],collapse="<br>",sep="")
            if(type == 'targeted'){
                parts = unlist(strsplit(as.character(samp),"_"))
                bc = parts[length(parts)]
                finalSamp = paste(bc,paste(parts[-length(parts)],collapse="_",sep=""),sep=": ")
            } else {
                finalSamp = sub("__1","",samp)
            }
            detail[which(rownames(detail) == finalSamp),"Unexpected Mismatch(es)"] = unmms
            detail[which(rownames(detail) == finalSamp)," "] = "0FAIL"
        }
      }
    }
    if(!is.null(majc)){ detail[sub("__1","",majc[,1]),"Major Contamination"] = majc[,2] }
    if(!is.null(minc)){ detail[sub("__1","",minc[,1]),"Minor Contamination"] = minc[,2] }
    if(!is.null(cov)){ detail[sub("__1","",cov[,1]),"Coverage"] = round(cov[,2]) }
    if(!is.null(dup)){ detail[sub("__1","",dup[,1]),"Duplication"] = signif(dup[,2]*100,digits=3) }
    if(!is.null(lib)){ detail[sub("__1","",lib[,1]),"Library Size (millions)"] = round(lib[,2]/1000000) }
    if(!is.null(cs)){ detail[sub("__1","",cs[,1]),"On Bait Bases (millions)"] = round(cs[,2]/1000000) }
    if(!is.null(aln)){ detail[sub("__1","",aln[,1]),"Aligned Reads (millions)"] = round(aln[,2]/1000000) }
    if(!is.null(is)){ detail[sub("__1","",is$Sample),"Insert Size Peak"] = is$Peak }  
    if(!is.null(tr)){ detail[sub("__1","",tr[,1]),"Percentage Trimmed Reads"] = apply(tr[,c(2,3)],1,sum) }

    ## assign status to each sample based on fixed thresholds
    detail[intersect(which(as.numeric(detail[,"Duplication"])>50),which(detail[,1]=="2PASS")),1] = "1WARN" ## warn high duplication
    detail[intersect(which(as.numeric(detail[,"Coverage"])<200),which(detail[,1]=="2PASS")),1] = "1WARN" ## warn low-ish coverage
    detail[which(as.numeric(detail[,"Coverage"])<50),1] = "0FAIL"  ## fail low coverage
    detail[which(as.numeric(detail[,"Minor Contamination"])>0.02),1] = "0FAIL" ## fail minor contam
    detail[which(as.numeric(detail[,"Major Contamination"])>0.55),1] = "0FAIL" ## fail major contam

    detail = as.data.frame(detail)

    if(!is.null(user) && nchar(user)>0 && user!='Unidentified'){

        lastReview = get.latest.qc.review(id)
        failures = detail$Sample[which(detail[,1] == '0FAIL')]
        actions = c()
        comments = c()
        for(samp in detail$Sample){
            akeep = paste("<select name='action_",samp,"'>",
                            " <option value='keep' selected>Keep</option>",
                            " <option value='drop'>Drop</option>",
                            " <option value='und'>...</option>",
                            "</select>",sep="")
            adrop = paste("<select name='action_",samp,"'>",
                            " <option value='keep' >Keep</option>",
                            " <option value='drop' selected>Drop</option>",
                            " <option value='und'>...</option>",
                            "</select>",sep="")
            aund = paste("<select name='action_",samp,"'>",
                            " <option value='keep'>Keep</option>",
                            " <option value='drop'>Drop</option>",
                            " <option value='und' selected>...</option>",
                            "</select>",sep="")

            action = akeep
            comment = ""
            if(is.null(lastReview)){ ## determine actions based on qc thresholds
                if(samp %in% failures){
                    action = adrop
                }
            } else{ ## determine actions based on previous qc review
                rev = lastReview$rev
                a = as.vector(rev[which(rev$Sample == samp),"Action"])[1]
                c = as.vector(rev[which(rev$Sample == samp),"Comment"])[1]
                if(a == "drop"){
                    action = adrop
                } else if(a == "..."){
                    action = aund
                }
                if(!is.na(c) && !is.null(c)){
                    comment = c
                }
            }
            actions = c(actions,action)
            comments = c(comments,comment)
        }
        detail$Action = actions
        detail$Comment = paste("<input type='text' name='comment_",detail$Sample,"' value='",comments,"'>",sep="")

        inp = c("Action","Comment")

        colOrder = c(other,inp,mets)
        detail = detail[colOrder]
    }

    detail = detail[order(detail[," "]),]
    return(detail)
}

samples.unexpected.matches <- function(path,type){
    
    alert = list("status"="2PASS","fails"=NULL)

    unmatch <- get.unexpected.matches(path,type)
    if(is.null(unmatch) || dim(unmatch)[1] == 0){
        alert[["fails"]] = "None"
        alert[["num"]] = 0
        return(alert)
    }

    alert[["num"]] <- length(unique(c(as.vector(unmatch[,1]),as.vector(unmatch[,2]))))
    unmatch <- paste(unmatch[,1],unmatch[,2],sep="&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;")  

    if (length(unmatch) == 0){
        alert[["fails"]] = "None"
        alert[["num"]] = 0
        return(alert)
    } 
    if (length(unmatch) %in% seq(1,5)){
        alert[["status"]] = "0FAIL"
        alert[["fails"]]<- paste(unmatch,collapse="</br>")
    } else {
        alert[["status"]] <- "0FAIL"
        alert[["fails"]] <- paste("<strong>",as.character(length(unmatch))," pairs.</strong>",sep="")
    }
    return(alert)
}

samples.unexpected.mismatches <- function(path,type){

    alert = list("status"="2PASS","fails"=NULL)

    unmismatch = get.unexpected.mismatches(path,type)
    if(is.null(unmismatch) || dim(unmismatch)[1] == 0){
        alert[["fails"]] = "None"
        alert[["num"]] = 0
        return(alert)
    }

    alert[["num"]] <- length(unique(c(as.vector(unmismatch[,1]),as.vector(unmismatch[,2]))))
    unmismatch <- paste(unmismatch[,1],unmismatch[,2],sep="&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;")  

    if (length(unmismatch) == 0){
        alert[["fails"]] = "None"
        return(alert)
    } 
    if (length(unmismatch) %in% seq(1,5)){
        alert[["status"]] = "0FAIL"
        alert[["fails"]]<- paste(unmismatch,collapse="</br>")
    } else {
        alert[["status"]] <- "0FAIL"
        alert[["fails"]] <- paste("<strong>",as.character(length(unmismatch))," pairs.</strong>",sep="")
    }
    return(alert)
}

samples.small.library <- function(path,type,min.lib.size=10000000){

    hs.df = get.library.size(path,type)
    samps <- hs.df[which(hs.df$LS<min.lib.size),]["Samples"]
    if(length(samps[,1]) == 0){
        return("None")
    }
    if(length(samps[,1]) %in% seq(1,5)){
        return(paste(samps[,1],collapse="</br>",sep=""))
    } else {
        return(paste("<strong>",as.character(length(samps[,1]))," samples.</strong>",sep=""))
    }

}

samples.major.contamination <- function(path,type,max.het=0.55){

    fpHet = get.major.contamination(path,type)
    contam.samples <- fpHet[which(fpHet$PerHetrozygosPos>max.het),1]
    samps = list("fail"=c(),"warn"=c())
    samps[["fail"]] <- as.vector(contam.samples)

    return(create.alert(samps))
}

samples.minor.contamination <- function(path,type,max.hom=0.02){

    fpHom = get.minor.contamination(path,type)
    contam.samples <- fpHom[which(fpHom$AvgMinorHomFreq>max.hom),1]

    samps = list("fail"=c(),"warn"=c())
    samps[["fail"]] <- as.vector(contam.samples)

    return(create.alert(samps))
}

samples.low.specificity <- function(min.spec=0.75){
    hs = get.capture.specificity(path,type)
    totals <- hs$ON_BAIT_BASES+hs$NEAR_BAIT_BASES+hs$OFF_BAIT_BASES
    on_bait_perc <- hs$ON_BAIT_BASES/totals

    samps = list("fail"=c(),"warn"=c())
    samps[["fail"]] <- as.vector(hs[which(on_bait_perc < min.spec),1])

    return(create.alert(samps))
}

samples.low.coverage <- function(path,type,min.cov.fail=50,min.cov.warn=200){

    cov = get.coverage(path,type)
    samps = list("fail"=c(),"warn"=c())
    samps[["fail"]] <- as.vector(cov[which(cov[,2]<min.cov.fail),1])
    samps[["warn"]] <- as.vector(cov[which(cov[,2]>=min.cov.fail & cov[,2]<min.cov.warn),1])
  
    return(create.alert(samps))
}

samples.low.alignment.rate <- function(path,type,min.alignment=0.95){
    hs = get.alignment(path,type)
    totals <- hs$both.reads.align + hs$one.read.aligns + hs$neither.read.aligns
    both.perc <- hs$both.reads.align/totals

    samps = list("fail"=c(),"warn"=c())
    samps[["fail"]] <- as.vector(hs[which(both.perc<min.alignment),1])
  
    return(create.alert(samps))
}

samples.high.duplication <- function(path,type,max.dup=0.50){
    duplication = get.duplication(path,type)
    samps = list("fail"=c(),"warn"=c())
    samps[["warn"]] <- as.vector(duplication[which(duplication$DupRate>max.dup),1])
    
    return(create.alert(samps))
}

create.alert <- function(sampList){
    alerts <- list()
    status <- "2PASS"
    failedSamps <- c(sampList[["fail"]],sampList[["warn"]])

    if(length(failedSamps) == 0){
        alerts[["status"]] <- status
        alerts[["fails"]] <- "None"
        alerts[["num"]] <- 0
        return (alerts)
    }

    ## set status
    if(length(sampList[["fail"]])>0){
        status <- "0FAIL"
    } else {
        status <- "1WARN"
    }
    
    ## set failure message (either short list of samples, or num total samples)
    if(length(failedSamps) %in% seq(1,5)){
        failMsg <- paste(failedSamps,collapse="</br>",sep="")
    } else {
        failMsg <- paste("<strong>",as.character(length(failedSamps))," samples.</strong>",sep="")
    }

    alerts[["status"]] <- status
    alerts[["fails"]] <- failMsg
    alerts[["num"]] <- length(failedSamps)

    return(alerts)
}

get.summary.table <- function(path,type){

    un.match = samples.unexpected.matches(path,type)
    un.mismatch = samples.unexpected.mismatches(path,type)
    maj.contam = samples.major.contamination(path,type)
    min.contam = samples.minor.contamination(path,type)
    cov = samples.low.coverage(path,type)
    dup = samples.high.duplication(path,type)

    mets = c("Unexpected Matches", "Unexpected Mismatches", "Major Contamination", "Minor Contamination", "Low Coverage","High Duplication")
    statuses = c(un.match$status, un.mismatch$status, maj.contam$status, min.contam$status, cov$status, dup$status)
    samps = c(un.match$fails, un.mismatch$fails, maj.contam$fails, min.contam$fails, cov$fails, dup$fails)
    nums = c(un.match$num, un.mismatch$num, maj.contam$num, min.contam$num, cov$num, dup$num)

    summ <- data.frame(Status=statuses,Metric=mets,Samples=samps,x=nums)
    colnames(summ)[4] = paste("N (total=",get.num.samples(path,type),")",sep="")

    return(summ)

}
