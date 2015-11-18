## take in path to metrics dir as argument
## take in 'pre' as argument


type = "exome"

args=(commandArgs(TRUE))
for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

source(file.path(bin,"get_metrics.R"))
source(file.path(bin,"plot_qc.R"))

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

fileName = paste(pre,"_QC.pdf",sep="")
fp = file.path(path,fileName)


pdf(file=fp, width=16,height=10)
## sample mixup heatmap
if(!is.null(da)){
  tryCatch(plot.fpc.sum(da), error = function(e){cat("ERROR: could not make fingerprint heatmap\n")})
} else {
  cat("ERROR: No fingerprint summary!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## coverage
if(!is.null(cv)){
  plot.coverage(cv)
} else {
  cat("ERROR: No coverage metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## major contamination
if(!is.null(mjc)){
  plot.major.contamination(mjc)
} else {
  cat("ERROR: No major contamination metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## minor contamination
if(!is.null(mnc)){
  plot.minor.contamination(mnc)
} else {
  cat("ERROR: No minor contamination metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## duplication
if(!is.null(dp)){
  plot.duplication(dp) 
} else {
  cat("ERROR: No duplication metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## library size
if(!is.null(ls)){
  plot.library.size(ls)
} else {
  cat("ERROR: No library size metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## alignment
if(!is.null(al)){
  plot.alignment(al)
} else {
  cat("ERROR: No alignment metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## alignment percentage
if(!is.null(al)){
  plot.alignment.percentage(al)
} else {
  cat("ERROR: No alignment metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## insert size distribution
if(!is.null(is)){
  plot.insert.size.distribution(is)
} else {
  cat("ERROR: No insert size metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## insert size peaks
if(!is.null(is)){
  plot.insert.peaks(is)
} else {
  cat("ERROR: No insert size metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## trimmed reads
if(!is.null(tr)){
  plot.trimmed.reads(tr)
} else {
  cat("ERROR: No read trimming metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## base qualitites
if(!is.null(bq)){
  plot.base.qualities(bq) 
} else {
  cat("ERROR: No base quality metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## capture specificity
if(!is.null(cs)){
  plot.capture.specificity(cs)
} else {
  cat("ERROR: No capture specificity metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
## capture specificity percentage
if(!is.null(cs)){
  plot.capture.specificity.percentage(cs)
} else {
  cat("ERROR: No capture specifity metrics!\n")
  quit(save="no",status=15,runLast=TRUE)
}
dev.off()
