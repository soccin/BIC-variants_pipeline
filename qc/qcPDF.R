## take in path to metrics dir as argument
## take in 'pre' as argument
options(echo = FALSE)

type = "exome"

args=(commandArgs(TRUE))
for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

cat(paste(date(),"\n"),file=logfile, append=TRUE)

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

no.fails = TRUE

pdf(file=fp, width=16,height=10)
## sample mixup heatmap
if(!is.null(da)){
  tryCatch(plot.fpc.sum(da), 
   error = function(e){
    cat("ERROR: could not make FINGERPRINT HEATMAP\n",file=logfile, append=TRUE)
    no.fails = FALSE
   }
  )
} else {
  cat("ERROR: No fingerprint summary!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## coverage
if(!is.null(cv)){
 tryCatch(plot.coverage(cv),
   error = function(e){
    cat("ERROR: could not make COVERAGE plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 ) 
} else {
  cat("ERROR: No coverage metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## major contamination
if(!is.null(mjc)){
  tryCatch(plot.major.contamination(mjc),
    error = function(e){
    cat("ERROR: could not make MAJOR CONTAMINATION plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No major contamination metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## minor contamination
if(!is.null(mnc)){
  tryCatch(plot.minor.contamination(mnc),
    error = function(e){
    cat("ERROR: could not make MINOR CONTAMINATION plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No minor contamination metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## duplication
if(!is.null(dp)){
  tryCatch(plot.duplication(dp),
   error = function(e){
    cat("ERROR: could not make DUPLICATION plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 ) 
} else {
  cat("ERROR: No duplication metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## library size
if(!is.null(ls)){
  tryCatch(plot.library.size(ls),
   error = function(e){
    cat("ERROR: could not make LIBRARY SIZE plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No library size metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## alignment
if(!is.null(al)){
  tryCatch(plot.alignment(al),
   error = function(e){
    cat("ERROR: could not make ALIGNMENT plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No alignment metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## alignment percentage
if(!is.null(al)){
  tryCatch(plot.alignment.percentage(al),
   error = function(e){
    cat("ERROR: could not make ALIGNMENT PERCENTAGE plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No alignment metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## insert size distribution
if(!is.null(is)){
  tryCatch(plot.insert.size.distribution(is),
   error = function(e){
    cat("ERROR: could not make INSERT SIZE DISTRIBUTION plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No insert size metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## insert size peaks
if(!is.null(is)){
  tryCatch(plot.insert.peaks(is),
   error = function(e){
    cat("ERROR: could not make INSERT SIZE plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
  )
} else {
  cat("ERROR: No insert size metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## trimmed reads
if(!is.null(tr)){
  tryCatch(plot.trimmed.reads(tr),
   error = function(e){
    cat("ERROR: could not make TRIMMED READS plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No read trimming metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## base qualitites
if(!is.null(bq)){
  tryCatch(plot.base.qualities(bq),
   error = function(e){
    cat("ERROR: could not make BASE QUALITIES plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 ) 
} else {
  cat("ERROR: No base quality metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## capture specificity
if(!is.null(cs)){
  tryCatch(plot.capture.specificity(cs),
   error = function(e){
    cat("ERROR: could not make CAPTURE SPECIFICITY plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No capture specificity metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
## capture specificity percentage
if(!is.null(cs)){
  tryCatch(plot.capture.specificity.percentage(cs),
   error = function(e){
    cat("ERROR: could not make CAPTURE SPECIFICITY PERCENTAGE plot\n",file=logfile, append=TRUE)
    no.fails = FALSE
   } 
 )
} else {
  cat("ERROR: No capture specifity metrics!\n",file=logfile, append=TRUE)
  no.fails = FALSE
}
dev.off()

if (no.fails == FALSE){
    quit(save="no",status=15,runLast=TRUE)
}
