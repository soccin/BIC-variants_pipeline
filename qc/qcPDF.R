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
isp = get.is.peaks(path,type)
hs = get.hs.metrics(path,type)
bq = get.base.qualities(path,type)
dp = get.duplication(path,type)
ls = get.library.size(path,type)
al = get.alignment(path,type)
tr = get.trimmed.reads(path,type)
cs = get.capture.specificity(path,type)

fileName = paste(pre,"_QC.pdf",sep="")
fp = file.path(path,fileName)


pdf(file=fp, width=16,height=10)
if(!is.null(dp)) plot.duplication(dp) else cat("ERROR: No duplication metrics!\n")
if(!is.null(ls)) plot.library.size(ls) else cat("ERROR: No library size metrics!\n")
if(!is.null(al)) plot.alignment(al) else cat("ERROR: No alignment metrics!\n")
if(!is.null(al)) plot.alignment.percentage(al) else cat("ERROR: No alignment metrics!\n")
if(!is.null(is)) plot.insert.size.distribution(is) else cat("ERROR: No insert size metrics!\n")
if(!is.null(is)) plot.insert.peaks(is) else cat("ERROR: No insert size metrics!\n")
if(!is.null(tr)) plot.trimmed.reads(tr) else cat("ERROR: No read trimming metrics!\n")
if(!is.null(bq)) plot.base.qualities(bq) else cat("ERROR: No base quality metrics!\n")
if(!is.null(cs)) plot.capture.specificity(cs) else cat("ERROR: No capture specificity metrics!\n")
if(!is.null(cs)) plot.capture.specificity.percentage(cs) else cat("ERROR: No capture specifity metrics!\n")
dev.off()
