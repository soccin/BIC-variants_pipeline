args<-commandArgs(TRUE)
INFILE=args[1]
print(INFILE)
dd=read.delim(INFILE,check.names=F)

jj=dd$CHROM!="chrM" & dd$CHROM!="chrY" & dd$CHROM!="chrX"
dd=dd[jj,]
qqdd=quantile(dd$INFO.DP)
dd=dd[dd$INFO.DP>qqdd[2] & dd$INFO.DP<qqdd[4],]

dx=dd[,1:8]
gt=dd[,9:ncol(dd)]
gb=gt=="0/0"

jj=apply(gb,1,mean)>.2 & apply(gb,1,mean)<.8
dx=dx[jj,]
gt=gt[jj,]
gb=gb[jj,]

ii=rev(order(dx$QUAL))[1:1000]
gbb=gb[ii,]
pdf(paste(gsub(".txt","",INFILE),"_TREE",".pdf",sep=""),height=8.5,width=9)
tree=as.dendrogram(hclust(dist(t(gbb),method="binary")))
mar.orig=par()$mar
par(mar=mar.orig+c(-1,-2,-3,6),cex=0.6)
plot(tree,horiz=T,xlim=c(1,0))
dev.off()
#dev.copy2pdf(file=cc(gsub(".txt","",INFILE),"TREE",DATE(),".pdf"),height=8.5,width=8.5)

library(gplots)
pdf(paste(gsub(".txt","",INFILE),"_HeatMap",".pdf",sep=""),height=8.5,width=9)
heatmap.2(as.matrix(dist(t(gbb),method="binary")),
    trace="none",Rowv=F,Colv=F,
    dendrogram="none",
    #colsep=1:8,rowsep=1:8,sepwidth=c(.01,.01),sepcolor=1,
    cexRow=0.4,cexCol=0.4,
    col=colorpanel(9,"#880000","#EEEEFF"),margins=c(12,12))
dev.off()
#dev.copy2pdf(file=cc(gsub(".txt","",INFILE),"HeatMap",DATE(),".pdf"),height=8.5,width=8.5)
