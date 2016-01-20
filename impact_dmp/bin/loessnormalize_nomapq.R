#!/ifs/opt/bin/Rscript

########################
# Loess normalize target exon coverage based on GC content.
# Modified from Rose Brannon's script: LoessNormalizeExonData_v1.1.R
# Changes:
# - Include 1000 tiling probe regions
# - Include X chromosome regions
# submit as
# R --slave --vanilla --args <prefix> <coverage-file> <gc_percent-file>
########################

# R --slave --vanilla --args <prefix> <_ALL_intervalcoverage.txt> <gc_percent-file>

rm(list=ls(all=T));

#source("textplot.R")
#setwd('/Volumes/phoenix-home/chengd1/LabPipeline/CopyNumberWork/CMO');
#prefix <- 'MPNST';
#outdir <- '/Volumes/phoenix-home/chengd1/LabPipeline/CopyNumberWork/CMO';
#gcfile <- '/Volumes/phoenix-home/chengd1/LabPipeline/ImpactPanelDesigns/cv3/cv3_hg19_all_GC150bp.txt';

#source("/dmp/resources/prod/tools/bio/misc/textplot.R")
source('textplot.R');
args <- commandArgs();
prefix <- args[5];
outdir <- args[6];   # Example: "/ifs/e63data/bergerm1/Analysis/Projects/Production/VEM/Vem3_TestCopyNumber";
gcfile <- args[7];   # Example: "v4_hg19_all_GC200bp.txt";

exon_file <- file.path(outdir,paste(prefix,"_ALL_exonnomapqcoverage.txt",sep=""));
tiling_file <- file.path(outdir,paste(prefix,"_ALL_tilingnomapqcoverage.txt",sep=""));
interval_file <- file.path(outdir,paste(prefix,"_ALL_intervalnomapqcoverage.txt",sep=""));
d_exon <- read.delim(exon_file,sep="\t",header=T,as.is=T,fill=T);
d_tiling <- read.delim(tiling_file,sep="\t",header=T,as.is=T,fill=T);
d <- rbind(d_exon, d_tiling);
write.table(d,interval_file,sep="\t",row.names=F,col.names=T,quote=F);

tiling.probes <- d_tiling[,'Target'];

outtablename <- file.path(outdir,paste(prefix,"_ALL_intervalnomapqcoverage_loess.txt",sep=""));
tar <- read.delim(gcfile,sep="\t",header=T,as.is=T,fill=T);

cat(dim(tar),"\t",dim(d),"\n");

#Order by target
chrom <- c(seq(1,22,1),'X','Y');
order.genomic <- function(dat){
	chr.pos <- do.call('rbind',lapply(dat[,'Target'], function(x){
		f <- unlist(strsplit(x,'\\:'));
		g <- unlist(strsplit(f[2],'\\-'));
		return(c('chr0'=f[1],'pos0'=g[1]));
	}));
	dat <- cbind(dat,chr.pos);
	ret <- do.call('rbind',lapply(chrom,function(chr.i){
		dat. <- dat[which(dat[,'chr0'] == chr.i),];
		dat. <- dat.[order(as.numeric(as.character(dat.[,'pos0'])),decreasing=F),];
		return(dat.);
	}));
	return(ret[,-c(ncol(ret)-1,ncol(ret))]);
}	
tar.sort <- order.genomic(tar);
d.sort <- order.genomic(d);

cat(dim(tar.sort),"\n");
cat(dim(d),"\n");
cat(dim(d.sort),"\n");

a <- tar.sort[,'Target'] != d.sort[,'Target'];
cat("Match test1 ",sum(a)," ",dim(tar.sort)," ",dim(d.sort),"\n");


targets <- tar.sort[which(!(tar.sort[,'Chr'] %in% c('chrY'))),];
gc2 <- d.sort[match(targets[,'Target'],d.sort[,'Target']),];


cat(gc2[1:10,'Target'],"\n");
cat(targets[1:10,'Target'],"\n");
cat(dim(gc2),"\n");
cat(dim(targets),"\n");
a <- gc2[,'Target'] != targets[,'Target'];
cat("Match test2  ",sum(a)," ",dim(gc2)," ",dim(targets),"\n");

index <- which(is.na(gc2[,'Target']));
cat(length(index),"\n");
index2 <- which(is.na(targets[,'Target']));
cat(length(index2),"\n");


if(sum(gc2[,'Target'] != targets[,'Target']) != 0){
	stop("Order of targets do not match\n");
}
gc2 <- gc2[,-c(1,2)];

GC<-as.numeric(targets[,grep("GC",colnames(targets))]);  
           

# New implementation - use optimize to search for function minimum - DC June 11th 2013
index <- which(targets[,'Target'] %in% tiling.probes);
span.fits <- do.call('rbind',list(apply(gc2,2,function(column){	
	column_sqrt<-sqrt(column[index]);
	gc.bias <- GC[index];

	testspan <- function(spanvalue){
		temp<-loess(column_sqrt~gc.bias,span=spanvalue);
		temp2<-predict(temp) #Calculation of the loess fit for each spanvalue
		normalized<-column_sqrt-temp2+median(column_sqrt) #Data normalized for each spanvalue
		
		fit<-loess(normalized~gc.bias);
		fit2<-predict(fit) #The "fit" of each normalized data point - the result gives the flat-ish line
		spanvar=var(fit2,na.rm=TRUE) #Calculate the variance to find the flattest line after fitting
		return(round(spanvar,5));
	}
	optimize.obj <-	optimize(testspan,interval=c(0.3,0.75),maximum=F);
	return(c('min'=optimize.obj$minimum,'obj'=optimize.obj$objective));
})));
span.fits <- t(span.fits);

out_file= paste(prefix,"_loessnorm.pdf", sep="")
pdf(out_file, height=8.5, width=11)
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(11, 10, "Memorial Sloan Kettering",adj=1)
text(11, 9.5, "Berger Lab",adj=1)
text(11, 9, date(),adj=1)
text(1, 8, "Loess Normalization results",adj=0,cex=1.5,font=2)

targets.inclY <- tar.sort;
gc2.inclY <- d.sort[match(targets.inclY[,'Target'],d.sort[,'Target']),];
gc2.inclY <- gc2.inclY[,-c(1,2)];

index <- which(targets[,'Target'] %in% tiling.probes);
index.inclY <- which(targets.inclY[,'Target'] %in% tiling.probes);

GC.inclY <-as.numeric(targets.inclY[,grep("GC",colnames(targets.inclY))]); 

norm_rt <- do.call('cbind',lapply(seq(1,ncol(gc2),1),function(i){
	column_sqrt<-sqrt(gc2.inclY[-index.inclY,i]);
	gc.bias <- GC.inclY[-index.inclY];
	loess.obj <-loess(column_sqrt~gc.bias,span=span.fits[i,'min']);
	temp2<-predict(loess.obj);
	normalized.filt<-(column_sqrt-temp2+median(column_sqrt))/(median(column_sqrt[which(column_sqrt != 0)]));
		
	column_sqrt.tiling <- sqrt(gc2[index,i]); # exclude Y
	gc.tiling <- GC[index]; 				  # exclude Y
	
	column_sqrt.tiling.inclY <- sqrt(gc2.inclY[index.inclY,i]); # exclude Y
	gc.tiling.inclY <- GC.inclY[index.inclY];
	
	loess.obj.tiling <-loess(column_sqrt.tiling~gc.tiling,span=span.fits[i,'min']); # Learn model excluding Y, but apply on newdata which has Y.
	temp2.tiling<-predict(loess.obj.tiling,newdata=gc.tiling.inclY);
	
	normalized.tiling<-(column_sqrt.tiling.inclY-temp2.tiling+median(column_sqrt))/(median(column_sqrt[which(column_sqrt != 0)]));

	temp2.all <- GC.inclY;
	temp2.all[which(targets.inclY[,'Target'] %in% tiling.probes)] <- temp2.tiling;
	temp2.all[which(!(targets.inclY[,'Target'] %in% tiling.probes))] <- temp2;

	normalized <-GC.inclY;
	normalized[which(targets.inclY[,'Target'] %in% tiling.probes)] <- normalized.tiling;
	normalized[which(!(targets.inclY[,'Target'] %in% tiling.probes))] <- normalized.filt;
	
	## Debug: Plot
	column_sqrt.all <- sqrt(gc2.inclY[,i]);
	par(mfrow=c(2,2))
	plot(GC.inclY[-index.inclY],column_sqrt.all[-index.inclY],ylim=c(0,50),main=paste("SqRt_",colnames(gc2)[i],sep=""),col='black',xlim=c(0.2,0.9),xlab='pGC',ylab='sqrt_cov');
	par(new=T);
	plot(GC.inclY[index.inclY],column_sqrt.all[index.inclY],ylim=c(0,50),col='red',xlim=c(0.2,0.9),xlab='',ylab='');	
	legend(x='topright',col=c('black','red'),legend=c('Target','Tiling'),pch=1);
	par(new=F);
	
	plot(GC.inclY,temp2.all,ylim=c(0,50),main=paste("Loess fit. Span:",span.fits[i,'min'],sep=""));
	##Normalize the data for the purpose of the graph.
	
        #normalized<-column_sqrt.all-temp2.all+median(column_sqrt);

	plot(GC.inclY,normalized,ylim=c(0,2),main="Normalized")
	fit<-loess(normalized~GC.inclY);
	fit2<-predict(fit,newdata=GC.inclY);
	plot(GC.inclY,fit2,ylim=c(0,2),main="Normalized fit")	
	return(normalized);
}));
dev.off();

norm=(norm_rt^2);
colnames(norm)<-colnames(gc2.inclY);
normtable=data.frame('Order'=targets.inclY[,"Order"],'Target'=targets.inclY[,"Target"],'genes'=as.character(targets.inclY[,"Gene"]),norm);

write.table(normtable,outtablename,sep="\t",row.names=FALSE,col.names=TRUE)
