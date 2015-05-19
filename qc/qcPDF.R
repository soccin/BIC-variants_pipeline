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
plot.duplication(dp)
plot.library.size(ls)
plot.alignment(al)
plot.alignment.percentage(al)
plot.insert.size.distribution(is)
plot.insert.peaks(is)
plot.trimmed.reads(tr)
plot.base.qualities(bq)
plot.capture.specificity(cs)
plot.capture.specificity.percentage(cs)
dev.off()
