#projID="Proj_4340"
#manifestPairing="/home/byrne/fingerprint/tests/Proj_4340_IBS/Proj_4340_sample_pairing.txt"
#IBSFile="/home/byrne/fingerprint/tests/Proj_4340_IBS/Proj_4340_pairs.txt_IBS"

args=(commandArgs(TRUE))
for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

#cat(c(projID,"\n",manifestPairing,"\n",IBSFile,"\n"))

## pairing according to manifest
manifestPairs=as.matrix(read.delim(manifestPairing,sep="",header=F))

manifestPairs=manifestPairs[which(!is.na(manifestPairs[,1]) & !is.na(manifestPairs[,2])),]
manifestPairNames=c(paste(manifestPairs[,1],manifestPairs[,2],sep="___"),paste(manifestPairs[,2],manifestPairs[,1],sep="___"))

pairs=as.matrix(read.delim(IBSFile,sep="\t",header=T))
pairnames=paste(pairs[,1],pairs[,2],sep="___")
meanIBS=as.numeric(pairs[,4])
IBSvar=as.numeric(pairs[,5])

pdf(paste(projID,"IBS_plot.pdf",sep="_"))
par(mar=c(5,5,5,12),xpd=TRUE)
plot(meanIBS,IBSvar,type='p',xlab="IBS mean",ylab="IBS variance",cex=0.45,pch=16,col="gray",bty="n")
title(projID)
colors=as.factor(pairnames[which(pairnames %in% manifestPairNames)])
points(x=meanIBS[which(pairnames %in% manifestPairNames)],IBSvar[which(pairnames %in% manifestPairNames)],cex=0.7,col=colors,pch=16)

legendX=max(meanIBS)+0.1
legendY=max(IBSvar)
legend(legendX,legendY,pairnames[which(pairnames %in% manifestPairNames)],col=colors,cex=0.6,pch=16,bty="n")

#text(x=as.numeric(pairs[which(pairnames %in% manifestPairNames),4]),as.numeric(pairs[which(pairnames %in% manifestPairNames),5]),labels=pairnames[which(pairnames %in% manifestPairNames)],cex=0.6,pos=4,col="red")
dev.off()
