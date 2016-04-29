## Usage: /opt/common/CentOS_6/R-3.2.0/bin/Rscript qc_images.R "pre='$PRE'" "bin='$BIN'" "path='$PATH_TO_METRICS'" "logfile='$LOGFILE'"

options(echo = FALSE)

type = "exome"

args=(commandArgs(TRUE))
for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

path = as.character(path)

source(file.path(bin,"get_metrics_from_files.R"))
source(file.path(bin,"plot_qc.R"))

print.image <- function(dat,metricType,sortOrder,plot.function,square=FALSE){
    units = "in"
    width = 11
    height = 7
    res = 600
    type = "cairo"
    if(square == TRUE){
        width = 5
        height = 5
    }

    if(!is.null(dat)){
        tryCatch({
                    #png(filename=paste(path,"/images/",pre,"_",sortOrder,"_",metricType,".png",sep=""),type=type,units=units,width=width,height=height,res=res)
                    pdf(file=paste(path,"/images/",pre,"_",sortOrder,"_",metricType,".pdf",sep=""),width=width,height=height)
                    print(plot.function(dat))
                    #plot.function(dat)
                    dev.off()
                 },
          error = function(e){
            cat(paste("ERROR: could not write ",metricType," image\n",sep=""), file=logfile, append=TRUE)
            cat(paste(e,"\n"),file=logfile,append=TRUE)
            no.fails = FALSE
          }
        )
    } else {
        cat(paste("ERROR: could not get ",metricType," metrics\n",sep=""),file=logfile, append=TRUE)
        no.fails = FALSE
    }
}

cat(paste(date(),"\n"),file=logfile, append=TRUE)
no.fails = TRUE

dir.create(paste(path,"/images",sep=""),showWarnings=FALSE)

## read in all metrics from text files in path
cat(c(path,"\n"))
is = get.is.metrics(path,type)
hs = get.hs.metrics(path,type)
bq = get.base.qualities(path,type)
dp = get.duplication(path,type)
ls = get.library.size(path,type)
al = get.alignment(path,type)
tr = get.trimmed.reads(path,type)
cs = get.capture.specificity(path,type)
cv = get.coverage(path,type)
da = get.fpc.sum(path,type)
mjc = get.major.contamination(path,type)
mnc = get.minor.contamination(path,type)
cc = get.cdna.contamination(path,type)
gc = get.gc.bias(path,type)

is.summary = NULL
dp.summary = NULL
al.summary = NULL
cs.summary = NULL
cv.summary = NULL
mjc.summary = NULL
mnc.summary = NULL

## get summary values (averages, etc)
if(!is.null(is)){ is.summary = get.mean.is.peak(type,dat=get.is.peaks(path,type)) } 
if(!is.null(dp)){ dp.summary = get.mean.duplication(type,dat=dp)*100 }
if(!is.null(ls)){ ls.summary = get.mean.library.size(type,dat=ls) }
if(!is.null(al)){ al.summary = get.alignment.totals(type,dat=al) }
if(!is.null(cs) && !is.null(hs)){ cs.summary = get.capture.specificity.summary(type,cs=cs,hs=hs) }
if(!is.null(cv)){ cv.summary = get.mean.coverage(type,cv)$All }
if(!is.null(mjc)){ mjc.summary = get.mean.frac.het.pos(type,dat=mjc) }
if(!is.null(mnc)){ mnc.summary = get.mean.minor.allele.freq(type,dat=mnc) }

## print images
print.image(al,"alignment","01",plot.alignment)
print.image(al,"alignment_percentage","02",plot.alignment.percentage)
print.image(cs,"capture_specificity","03",plot.capture.specificity)
print.image(cs,"capture_specificity_percentage","04",plot.capture.specificity.percentage)
print.image(is,"insert_size","05",plot.insert.size.distribution)
print.image(is,"insert_size_peaks","06",plot.insert.peaks)
print.image(da,"fingerprint","07",plot.fpc.sum) #,square=TRUE) }
print.image(mjc,"major_contamination","08",plot.major.contamination) 
print.image(mnc,"minor_contamination","09",plot.minor.contamination) 
print.image(cc,"cdna_contamination","10",plot.cdna.contamination) 
print.image(dp,"duplication","11",plot.duplication) 
print.image(ls,"library_size","12",plot.library.size) 
print.image(cv,"coverage","13",plot.coverage) 
print.image(tr,"trimmed_reads","14",plot.trimmed.reads) 
print.image(bq,"base_qualities","15",plot.base.qualities) 
print.image(gc,"gc_bias","16",plot.gc.bias) #,square=TRUE)

## write sample level summary table
tryCatch({
    detail = as.matrix(get.detail.table(path,type,id))
    colnames(detail)[1] = "Auto-status"
    detail[which(detail[,1]=='0FAIL'),1] = 'FAIL'
    detail[which(detail[,1]=='1WARN'),1] = 'WARN'
    detail[which(detail[,1]=='2PASS'),1] = 'PASS'

    summary.row = rep("",ncol(detail))
    summary.row[which(colnames(detail)=="Sample")] = "Project Average"
    if(!is.null(mjc.summary)){ summary.row[which(colnames(detail)=="Major Contamination")] = mjc.summary }
    if(!is.null(mnc.summary)){ summary.row[which(colnames(detail)=="Minor Contamination")] = mnc.summary }
    if(!is.null(cv.summary)){ summary.row[which(colnames(detail)=="Coverage")] = round(cv.summary) }
    if(!is.null(dp.summary)){ summary.row[which(colnames(detail)=="Duplication")] = dp.summary }
    if(!is.null(ls.summary)){ summary.row[which(colnames(detail)=="Library Size (millions)")] = ls.summary }
    if(!is.null(is.summary)){ summary.row[which(colnames(detail)=="Insert Size Peak")] = round(is.summary) }
    if(!is.null(cs.summary)){ summary.row[which(colnames(detail)=="On Bait Bases (millions)")] = round(cs.summary$meanOnBait/1000000) }
    if(!is.null(al.summary)){ summary.row[which(colnames(detail)=="Aligned Reads (millions)")] = round((al.summary$totalClusters/nrow(detail))/1000000) }
    summary.row[which(colnames(detail)=="Percentage Trimmed Reads")] = round(mean(as.numeric(detail[,which(colnames(detail)=="Percentage Trimmed Reads")])),digits=2) 

    detail = rbind(summary.row,detail)
    write.table(detail,file=file.path(path,paste(pre,"_SampleSummary.txt",sep="")),sep="\t",quote=F,row.names=F,col.names=T)
}, error = function(e){
        cat(paste("ERROR: could not write ",paste(path,"/",pre,"_SampleSummary.txt",sep=""),"\n",sep=""), file=logfile,append=TRUE)
        cat(paste(e,"\n"),file=logfile,append=TRUE)
        no.fails = FALSE
})


## write project level summary table
tryCatch({
    write.table(get.project.summary(type,path),file=file.path(path,paste(pre,"_ProjectSummary.txt",sep="")),sep="\t",quote=F,row.names=F,col.names=T)
}, error = function(e){
        cat(paste("ERROR: could not write ",paste(path,"/",pre,"_ProjectSummary.txt",sep=""),"\n",sep=""), file=logfile,append=TRUE)
        cat(paste(e,"\n"),file=logfile,append=TRUE)
        no.fails = FALSE
})

if (no.fails == FALSE){
    quit(save="no",status=15,runLast=TRUE)
}
