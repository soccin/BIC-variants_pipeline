# Donavan Cheng
# 2013 July 24th
# -- 
# 1. Use tiling probes
# 2. Compare vs. best normal method vs. all other normals within pool method
# 3. Compare vs. list of other unmatched normals

# 4. Normalization across other samples within pool - current: normalized wrt sample median, consider: quantile normalization

rm(list=ls(all=T));
library('MASS');
library('gplots');
library('RColorBrewer');
library('DNAcopy');
library('Ckmeans.1d.dp');
library('rjson');

#source('/dmp/resources/prod/tools/bio/misc/textplot.R');
source('textplot.R');
args <- commandArgs();
prefix <- args[5];
stdnormal_loess <- args[6];
gene_intervals.ann <- args[7];
tiling_intervals.ann <- args[8];
tested.class.str <- args[9];
doFull <- args[10];
if(doFull != "FULL" & doFull != "MIN"){
    stop("doFull should be either FULL or MIN\n");
}     

#setwd('/Volumes/phoenix-home/chengd1/LabPipeline/CopyNumberWork');
#source('textplot.R');
#prefix <- 'IMPACTv3-CLIN-20140010-Rerun';
#stdnormal_loess <- 'CombinedStdNormals';
#gene_intervals.ann <- 'gene_intervals.list.annotated';
#tiling_intervals.ann <- 'tiling_intervals.list.annotated';
#tested.class.str <- 'NormalGL,Normal';
#doFull <- 'MIN';

tested.class <- unlist(strsplit(tested.class.str,'\\,'));
	
#collecting the barcodes from the title slide
titlefile <- paste(prefix,"_title.txt",sep="")
title <- read.delim(titlefile,sep="\t",header=T,as.is=T,fill=T);

barcode <- title$Sample_ID;
covfile_tiling <-"_ALL_tilingcoverage.txt";
covfile_loess <-"_ALL_intervalcoverage_loess.txt";

gc_tiling <- read.delim(paste(prefix,covfile_tiling,sep=""),sep="\t",header=T,as.is=T,fill=T);
tiling_probe_ids <- gc_tiling[,'Target'];

gc<-read.delim(paste(prefix,covfile_loess,sep=""),sep="\t",header=T,as.is=T,fill=T);	
stdn_cov <- read.delim(paste(stdnormal_loess,"_ALL_intervalcoverage_loess.txt",sep=""),sep="\t",header=T,as.is=T,fill=T);
stdn_cov.title <- read.delim(paste(stdnormal_loess,"_title.txt",sep=""),sep="\t",header=T,as.is=T,fill=T);

gene.ann <- read.delim(gene_intervals.ann,sep='\t',header=F,as.is=T,fill=T);
colnames(gene.ann) <- c('Target','Cyt','GeneExon');
tiling.ann <- read.delim(tiling_intervals.ann,sep='\t',header=F,as.is=T,fill=T); 
colnames(tiling.ann) <- c('Target','Cyt');

stdn_cov <- stdn_cov[,c(1:3,which(colnames(stdn_cov) %in% stdn_cov.title[which(stdn_cov.title[,'Class'] == 'Normal'),'Barcode']))];

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
gc.sort <- order.genomic(gc);
gc <- gc.sort;
stdn_cov.sort <- order.genomic(stdn_cov);
stdn_cov <- stdn_cov.sort;

chr <- do.call('c',lapply(gc[,'Target'],function(x){
	f <- unlist(strsplit(x,'\\:'));
	return(f[1]);
}));
chr <- factor(chr,levels=c(seq(1,22,1),'X'));
chr.counts <- table(chr);
chr.midpt <- sapply(chr.counts,function(x){return(round(x/2,0));});
numberpos <- chr.midpt + c(0,cumsum(chr.counts)[-length(chr.counts)]);
linepos <- cumsum(chr.counts)[-length(chr.counts)];

d <-gc[,-c(1:3)];
colnames(d) <- title[match(colnames(d),title[,'Barcode']),'Sample_ID'];

panel.genes <- unique(gc[which(!grepl('Tiling',gc[,3]) & !grepl('FP_',gc[,3])),3]);
panel.genes <- panel.genes[which(!(panel.genes %in% c("","CDKN2A")))];
panel.genes <- c(panel.genes,'CDKN2Ap16INK4A','CDKN2Ap14ARF');
panel.genes <- unique(panel.genes);
gene.ann <- cbind(gene.ann,'gene'=do.call('c',lapply(gene.ann[,'GeneExon'],function(x){
			     	                      if(x != ""){
			     	            		  		f <- unlist(strsplit(x,'\\,'));
			     	            		  		genes <- do.call('c',lapply(f,function(f.i){
			     	            		  		                            g <- unlist(strsplit(f.i,'\\:'));			     	            		  			
			     	            		  		                            return(g[1]);
			     	            		  		                 }));
												idx <- which(genes %in% panel.genes);
												if(length(idx) >= 1){
													return(paste(unique(genes[idx]),collapse=','));
												}else{
													return("");
												}			     	            		  		                 
			     	            		  }else{
			     	            		  		return("");
			     	            		  }
			     	            		  })));
all.ann <- rbind(gene.ann,cbind(tiling.ann,'GeneExon'=rep('-',nrow(tiling.ann)),'gene'=rep('',nrow(tiling.ann))));

gc.target <- as.character(all.ann[match(gc[,'Target'],all.ann[,'Target']),'gene']);
genelevel<-cbind(gc[,1:2],'genes'=gc.target,d) #cbinds annotation information (chr, gene, etc.) to the barcode informatio

common.targets <- genelevel[which(genelevel[,'Target'] %in% stdn_cov[,'Target']),'Target'];
genelevel <- genelevel[match(common.targets,genelevel[,'Target']),];

d.stdn <- stdn_cov[match(common.targets,stdn_cov[,'Target']),-c(1:3)];		
colnames(d.stdn) <- paste('STDN_',colnames(d.stdn),sep="");

#TestedClass
index <- which(title$Class %in% tested.class); #selects the tumor barcodes 
if(length(index)==0){
	res0 <- paste('No samples matching tested classes in pool:',tested.class.str,collapse='');
	write.table(res0, file=paste(prefix,"_nvn_results.txt", sep=""),sep="\t",row.names=F,col.names=T,quote=F); 
	stop('No tested classes in pool\n');
}
test <- do.call('cbind',list(apply(genelevel[,title[index,'Sample_ID'],drop=F],2,function(x){
      ret <- as.numeric(x);
      ret[which(is.infinite(ret) | is.nan(ret))] <- NA;
      return(ret);
})));
gender <- rep('-',ncol(test));
if(any(grepl('Gender',colnames(title)))){
	gender <- as.character(title[index,'Gender',drop=F][,1]);
	idx <- which(gender %in% c('Male','M','male'));
	gender[idx] <- 'M';
	gender[-idx] <- 'F';
}
names(gender) <- colnames(test);

index <- which(title$Class=="PoolNormal"); #selects the normal barcodes 
PoolNormalPresent <- 0;
if(length(index)>=1){
	PoolNormalPresent <- 1;	
	if(length(index)>1){
		index <- index[1]; 
	}
	d.pooln <- do.call('cbind',list(apply(genelevel[,title[index,'Sample_ID'],drop=F],2,as.numeric)));	
	d.stdn <- cbind(d.stdn, d.pooln);
}

#Make a matrix of normal samples, then replace values less than .1*median value with NA
#This is to filter out positions with low reads that could skew the resulting ratios after division.
doMedianOnCols <- function(dat){
	normal.median <- apply(dat,2,median,na.rm=T);
	ret <- do.call('cbind',lapply(seq(1,ncol(dat),1),function(i){
		col.i <- dat[,i];
		col.i[which(col.i < 0.1*normal.median[i])] <- NA;  # Difference from previous
		col.i[which(is.infinite(col.i) | is.nan(col.i))] <- NA;  
		return(col.i);
	}));
	colnames(ret) <- names(normal.median);
	return(ret);
}
d.stdn <- doMedianOnCols(d.stdn);

low.cov.probes <- do.call('c',list(apply(d.stdn[,grepl('STDN',colnames(d.stdn))],1,function(cov){
	return(sum(is.na(cov))/length(cov));
})));
low.cov.probes <- cbind('Target'=genelevel[,'Target'],'pct_low'=low.cov.probes,all.ann[match(genelevel[,'Target'],all.ann[,'Target']),]);
low.cov.probes <- low.cov.probes[which(low.cov.probes[,'pct_low'] > 0.2 & !grepl('Y\\:',low.cov.probes[,'Target'])),];
idx <- which(genelevel[,'Target'] %in% low.cov.probes[,'Target']);
d.stdn[idx,] <- matrix(rep(NA,length(idx)*ncol(d.stdn)),nrow=length(idx),ncol=ncol(d.stdn));

y.probes <- genelevel[grep('^Y\\:',genelevel[,'Target']),'Target'];
y.probes <- y.probes[which(!(y.probes %in% c("Y:6111337-6111456","Y:6736141-6736260","Y:25310891-25311010","Y:9238599-9238718")))];
y.cov <- do.call('c',list(apply(d.stdn[which(genelevel[,'Target'] %in% y.probes),colnames(d.stdn)],2,function(col.i){
	col.i[which(is.na(col.i))] <- 0;
	return(median(as.numeric(col.i)));
})));
y.cov_grp <- abs(tapply(y.cov,stdn_cov.title[match(gsub('STDN_','',names(y.cov)),stdn_cov.title[,'Barcode']),'Sample_type'],mean,na.rm=T));	

calc.fc <- function(tm.vec, nm.vec){
	ratio <- tm.vec/nm.vec;
	lr <- log2(tm.vec) - log2(nm.vec);
	fc <- sapply(lr,function(x){
		if(is.na(x)){return(NA);}
		if(x > 0){return(2^x);}
		else{return(-2^(-x));}
	});
	return(cbind('ratio'=ratio,
				 'lr'=lr,
				 'fc'=fc));
}

levelplot <- function(yvalues,grouping,title){
	idx <- order(yvalues,decreasing=F);
	yvalues <- yvalues[idx];
	grouping <- grouping[idx];
	
	grouping[which(is.na(grouping))] <- 0;
	num.clus <- as.numeric(unique(grouping));
	num.clus <- num.clus[order(num.clus,decreasing=F)];
	my.col <- rainbow(length(num.clus)-1);
	my.col <- c('black',my.col);
	legend.txt <- do.call('c',lapply(num.clus,function(k){
		return(paste('Cluster ',k,sep=''));
	}));
	for(k in 1:length(num.clus)){
		a1 <- yvalues;
		a1[which(grouping!=num.clus[k])] <- NA;	
		plot(x=1:length(yvalues),a1,col=my.col[k],
		     xlim=c(1,length(yvalues)),xlab="",
		     ylim=c(min(yvalues),max(yvalues)),ylab="");
		par(new=T);
	}
	legend(x='bottomright',col=my.col,legend=legend.txt,lty=1,lwd=1.5,cex=0.8);
	par(new=F);
}

rough.clus <- function(yvalues,sep=0.1,num.steps=3,make.plot=F){
	seg.levels <- table(yvalues);
	seg.levels <- seg.levels[order(as.numeric(names(seg.levels)),decreasing=F)];
	dc.clus <- do.call('c',lapply(2:length(seg.levels),function(i){
		i.prev1 <- i-1;
		if(i.prev1 <= 0){
			i.prev1 <- NA;
		} 
		i.prev1.flag <- F;
		if(!is.na(i.prev1)){
			if(abs(as.numeric(names(seg.levels)[i])-as.numeric(names(seg.levels[i.prev1]))) > sep){
				i.prev1.flag <- T;
			}
		}				
		if(i.prev1.flag){				
			return(1);
		}else{
			return(0);
		}
	}));
	dc.clus <- c(0,dc.clus);
	
	clus <- 1;
	dc.clus.names <- c();
	for(i in 1:length(seg.levels)){		
		if(dc.clus[i] == 1){
			flag <- 1;
			if(flag){
				clus <- clus+1;
			}
		}
		dc.clus.names <- c(dc.clus.names,clus);		
	}
	small.clusters <- names(table(dc.clus.names)[which(table(dc.clus.names) < num.steps)]);
	dc.clus.names[which(dc.clus.names %in% small.clusters)] <- 0;
	names(dc.clus.names) <- names(seg.levels);
		
	if(make.plot){
		par(mfrow=c(3,1))
		plot(as.numeric(names(seg.levels)),main='Unclustered seg values',
 	    	 xlim=c(1,length(seg.levels)),xlab="",
 		 	 ylim=c(min(yvalues),max(yvalues)),ylab="");
		plot(dc.clus,xlim=c(1,length(seg.levels)),ylim=c(0,2))
		plot(dc.clus.names,xlim=c(1,length(seg.levels)),ylim=c(0,max(dc.clus.names)));
	}
		
	#browser();
	ret <- do.call('c',lapply(yvalues,function(y.i){
		return(dc.clus.names[as.character(y.i)]);
	}));	
	names(ret) <- names(yvalues);
	return(ret);
}		
				
seg.clus.p <- function(logratio, genelabel, make.plots=F,sample.id,bc){
	
	# CBS
	cna.obj.mat <- do.call('rbind',lapply(seq(1,length(logratio),1),function(j){
		f <- unlist(strsplit(names(logratio)[j],'\\:'));
		chr <- f[1];
		g <- unlist(strsplit(f[2],'\\-'));
		start <- as.numeric(g[1]); stop <- as.numeric(g[2]);
		if((stop-start)%%2==0){
			midpt <- start + (stop-start)/2;
		}else{
			midpt <- start + (stop-start+1)/2;
		}
		return(c('lr'=as.numeric(logratio[j]),'chr'=chr,'pos'=midpt));
	}));
	cna.obj <- CNA(as.numeric(cna.obj.mat[,'lr']),cna.obj.mat[,'chr'],as.numeric(cna.obj.mat[,'pos']),data.type='logratio',sampleid=sample.id);
	names(cna.obj)[3] <- sample.id;
	smoothed.cna.obj <- smooth.CNA(cna.obj);
	segment.smoothed.cna.obj  <- segment(smoothed.cna.obj, undo.splits="sdundo", undo.SD=2, verbose=1);
	segment.smoothed.cna.obj.srt <- do.call('rbind',lapply(c(1:22,"X"),function(chr){
		ret <- segment.smoothed.cna.obj$output[which(segment.smoothed.cna.obj$output[,'chrom'] == chr),,drop=F];
		return(ret[order(as.numeric(ret[,'loc.start']),decreasing=F),]);
	}));
	write.table(segment.smoothed.cna.obj.srt, 
			    paste(sample.id,"_",bc,"_",prefix,"_copynumber.nvn.seg",sep=""),
			    col.names=T,row.names=F,quote=F);
	data.srt <- do.call('rbind',lapply(c(1:22,"X"),function(chr){
		idx <- which(segment.smoothed.cna.obj$data$chrom == chr);		
		return(cbind('chr'=segment.smoothed.cna.obj$data$chrom[idx],
					 'maploc'=segment.smoothed.cna.obj$data$maploc[idx],
					 'values'=segment.smoothed.cna.obj$data[[length(segment.smoothed.cna.obj$data)]][idx]));
	}));
	segment.smoothed.cna.obj$data$chrom <- data.srt[,'chr'];
	segment.smoothed.cna.obj$data$maploc <- data.srt[,'maploc'];
	segment.smoothed.cna.obj$data[[3]] <- as.numeric(data.srt[,'values']);
	
	segment.smoothed.cna.obj$output <- segment.smoothed.cna.obj.srt;
	seg.out <- segment.smoothed.cna.obj$output;

	seg.probes <- do.call('c',lapply(names(logratio),function(targ.i){
		f <- unlist(strsplit(targ.i,'\\:'));
		chr <- f[1];
		g <- unlist(strsplit(f[2],'\\-'));
		start <- as.numeric(g[1]);
		end <- as.numeric(g[2]);
		
		X <- seg.out[which(seg.out[,'chrom'] == chr & !(end < as.numeric(seg.out[,'loc.start']) | start > as.numeric(seg.out[,'loc.end']))),,drop=F];
		if(nrow(X) == 0){
			return(NA);
		}else if(nrow(X) == 1){
			return(X[1,'seg.mean']);
		}else{
			cat("Multiple seg\t",sample.i,"\t",targ.i,"\n");
			return(mean(as.numeric(X[,'seg.mean'])));
		}
	}));
	names(seg.probes) <- names(logratio);
	not_na.idx <- which(!is.na(seg.probes));
	seg.probes <- seg.probes[not_na.idx];

	all.genes <- panel.genes;
	idx <- which(all.genes == 'PDPK1');
	if(length(idx)==1){
		all.genes <- all.genes[-idx];
	}

	seg.genes <- do.call('rbind',lapply(all.genes,function(gene.i){
		flag <- genelabel == gene.i | grepl(paste('^',gene.i,',',sep=''),genelabel) | grepl(paste(',',gene.i,'$',sep=''),genelabel) | grepl(paste(',',gene.i,',',sep=''),genelabel); 
		seg.gene <- seg.probes[names(genelabel)[which(flag)]];
		seg.gene <- seg.gene[which(!is.na(seg.gene))];
		f <- unlist(strsplit(names(seg.gene)[1],'\\:'));
		return(c('chr'=f[1],'num.exons'=length(seg.gene),'seg'=mean(seg.gene)));
	}));
	rownames(seg.genes) <- all.genes;			

	if(make.plots){
		a<-max(abs(logratio),na.rm=TRUE)
		if (a > 4){ a=signif(a, digits = 1)}else {a= 4}
		
		layout(matrix(c(1,1,2,2),nrow=2,ncol=2,byrow=T),heights=c(1,1),widths=c(1,1));	
		plot(segment.smoothed.cna.obj, plot.type="w",xlab="");
		abline(v=linepos,lty=3);
		abline(h=-1,lty=2);
		abline(h=1,lty=2);
		mtext(c(1:22,'X'),side=1,line=0,at=numberpos,cex=0.7);
		
		plot(seg.probes,xlab="Targeted exons per chromosome", ylab="LOG2 Tumor/Normal Ratio",main=paste("Segmented IMPACT LogRatio",sample.id,sep="\n"),xlim=c(0,length(seg.probes)),ylim=c(-a,a),xaxt="n",pch=20,col="black",mgp=c(2,1,0),mar=c(3,5,2,5));
		abline(v=linepos,lty=3);
		abline(h=-1,lty=2);
		abline(h=1,lty=2);
		mtext(c(1:22,'X'),side=1,line=0,at=numberpos,cex=0.7);
	}
		    
	all.clus <- rough.clus(seg.probes,sep=0.08,num.steps=3);
	all.clus.means <- tapply(seg.probes,all.clus,mean,na.rm=T);
	all.clus.p <- table(all.clus)/sum(table(all.clus));
	
	lr.clus.means <- tapply(logratio[not_na.idx],all.clus,mean,na.rm=T);
	lr.clus.sds <- tapply(logratio[not_na.idx],all.clus,sd,na.rm=T);
	idx <- which(names(lr.clus.means) != "0");
	lr.clus.means <- lr.clus.means[idx];
	lr.clus.sds <- lr.clus.sds[idx];
	cl0 <- names(lr.clus.means)[which(abs(lr.clus.means) == min(abs(lr.clus.means)))];
	cl0.mean <- lr.clus.means[cl0];
	cl0.sd <- lr.clus.sds[cl0];

	if(make.plots){
		par(mfrow=c(2,2));
		plot(seg.probes[order(seg.probes,decreasing=F)],main='Unclustered seg values',xlim=c(1,length(seg.probes)),xlab="",ylim=c(min(seg.probes),max(seg.probes)),ylab="");
		levelplot(seg.probes, all.clus, "Clustered seg values");
		for(clus.mean.i in all.clus.means[which(names(all.clus.means) != "0")]){
		    abline(h=clus.mean.i,lty=2);
		}	
		plot(density(logratio,na.rm=T),xlim=c(-2,2),ylim=c(0,1),main='Empirical LogRatio - Probes',xlab="Log2 T/N Ratio");
		com.col <- rainbow(length(lr.clus.means));		
		plot(density(logratio,na.rm=T),xlim=c(-2,2),ylim=c(0,1),main='Mixture Model - Probes',xlab="Log2 T/N Ratio");
		all.clus.p <- all.clus.p[names(lr.clus.means)];
		legend.txt <- do.call('c',lapply(1:length(lr.clus.means),function(k){
			lines(x=seq(-2,2,0.02),y=all.clus.p[k]*dnorm(seq(-2,2,0.02),mean=lr.clus.means[k],sd=lr.clus.sds[k]),xlim=c(-2,2),ylim=c(0,1),col=com.col[k],xlab="",ylab="",type='l');
			if(names(lr.clus.means)[k] == cl0){
				return(paste('Fit ',names(lr.clus.means)[k],' - NULL',sep=''));
			}else{
				return(paste('Fit ',names(lr.clus.means)[k],sep=''));
			}
		}));
		legend(x='bottomright',col=com.col,legend=legend.txt,lty=1,lwd=1.5,cex=0.8);
	}

	fc <- sapply(logratio,function(x){
		if(is.na(x)){return(NA);}
		if(x > 0){return(2^x);}
		else{return(-2^(-x));}
	});
	z <- (logratio-cl0.mean)/cl0.sd;
	p <- sapply(z,function(z.i){
		if(is.na(z.i)){return(NA);}
		if(z.i <= 0){ return(2*pnorm(z.i,mean=0,sd=1));}
	     else{ return(2*(1-pnorm(z.i,mean=0,sd=1)));}		
	});	
	p.adj <- p.adjust(p,method='BH');
	sig <- (fc >= 1.3 | fc <= -1.7) & p.adj < 0.05;
	sig[which(is.na(sig))] <- 0;
	
	fc.genes <- sapply(as.numeric(seg.genes[,'seg']),function(x){
		if(is.na(x)){return(NA);}
		if(x > 0){return(2^x);}
		else{return(-2^(-x));}
	});
	zg <- (as.numeric(seg.genes[,'seg'])-cl0.mean)/cl0.sd;
	pg <- sapply(zg,function(zg.i){
		if(is.na(zg.i)){return(NA);}
		if(zg.i <= 0){ return(2*pnorm(zg.i,mean=0,sd=1));}
	     else{ return(2*(1-pnorm(zg.i,mean=0,sd=1)));}		
	});	
	pg.adj <- p.adjust(pg,method='BH');
	sig.g <- (fc.genes >= 1.3 | fc.genes <= -1.7) & pg.adj < 0.05;
	sig.g[which(is.na(sig.g))] <- 0;

	ret.intragenic <- do.call('rbind',lapply(all.genes,function(gene.i){
		if(grepl('Tiling',gene.i)){return(NULL);}
		flag <- genelabel == gene.i | grepl(paste('^',gene.i,',',sep=''),genelabel) | grepl(paste(',',gene.i,'$',sep=''),genelabel) | grepl(paste(',',gene.i,',',sep=''),genelabel);  
		lr.x <- logratio[names(genelabel)[which(flag)]];
		if(grepl("^X",names(lr.x)[1])){return(NULL);}

		fc.x <- fc[names(genelabel)[which(flag)]];
		p.adj.x <- p.adj[names(genelabel)[which(flag)]];
		idx <- which(is.na(lr.x));
		if(length(idx) > 0){
			lr.x <- lr.x[-idx];
			fc.x <- fc.x[-idx];
			p.adj.x <- p.adj.x[-idx];
		}
				
		ann <- do.call('rbind',lapply(names(lr.x),function(targ.i){
			if(targ.i %in% gene.ann[,'Target']){
				return(gene.ann[which(gene.ann[,'Target']==targ.i),c('Cyt','GeneExon')]);
			}else{
				return(c('Cyt'=NA,'GeneExon'=NA));				
			}
		}));
		
		# Test for intragenic loss
		sig.x <- fc.x <= -1.7 & p.adj.x < 0.05;
		ret <- c();		
		if(length(lr.x) >= 3 & length(unique(lr.x)) > 2 & length(table(sig.x)) > 1){ 
			k2 <- Ckmeans.1d.dp(lr.x,2);
			k1 <- Ckmeans.1d.dp(lr.x,1);
			ratio <- k1$withinss/sum(k2$withinss);
									
			if(ratio >= 3.5){
				clus.del <- which(k2$centers == min(k2$centers));
				clus.not_del <- which(k2$centers != min(k2$centers));
				clus.diff <- k2$centers[clus.not_del] - k2$centers[clus.del];
				
				clus.del.members <- as.numeric(k2$cluster == clus.del);
				runsum <- 0;
				runsum_arr <- c();
				for(i in seq(1,length(clus.del.members),1)){
					if(clus.del.members[i]==0){
						runsum <- 0;
					}else{
						runsum <- runsum + clus.del.members[i];
					}
					runsum_arr <- c(runsum_arr,runsum);
				}
				max.runsum <- max(runsum_arr);
				not.in.stretch <- k2$size[clus.del] - max.runsum;
				pct.stretch <- max.runsum/k2$size[clus.del];
		
				if((k2$size[clus.del] == 1 & clus.diff >= 0.8) |
				   (k2$size[clus.del] > 1 & k2$size[clus.del] < 6 & clus.diff >= 0.5 & pct.stretch > 0.5) | 
				   (k2$size[clus.del] >= 6 & clus.diff >= 0.5 & pct.stretch >= 0.5)){								
				   	
				   	cluster <- k2$cluster;
				   	names(cluster) <- names(lr.x);
				   	probes <- names(genelabel)[which(flag)];
				   	lr.x.i <- logratio[probes];
				   	fc.x.i <- fc[names(genelabel)[which(flag)]];
				   	p.adj.x.i <- p.adj[names(genelabel)[which(flag)]];
				   	p.cluster <- rep(NA, length(lr.x.i));				   	
				   	names(p.cluster) <- probes;
				   	p.cluster[names(cluster)] <- cluster;
				   	ann <- do.call('rbind',lapply(names(lr.x.i),function(targ.i){
						if(targ.i %in% gene.ann[,'Target']){
							return(gene.ann[which(gene.ann[,'Target']==targ.i),c('Cyt','GeneExon')]);
						}else{
							return(c('Cyt'=NA,'GeneExon'=NA));				
						}
					}));
					ret <- rbind(ret,cbind('table'=rep('Intragenic_loss',length(lr.x.i)),'lr'=lr.x.i,'NExons'=rep(1,length(lr.x.i)),
					       	                        'chr'=do.call('c',lapply(names(lr.x.i),function(x){
									f <-  unlist(strsplit(x,'\\:'));
									return(f[1]);
			    				 })),'fc'=fc.x.i,'p.adj'=p.adj.x.i,'sig'=as.character(p.cluster),'sig.intragenic'=rep('-',length(lr.x.i)),ann));
					return(ret);
			    }				
			}
		}
				
		# Test for intragenic gain
		sig.x <- fc.x >= 1.3 & p.adj.x < 0.05;
		if(length(lr.x) >= 3 & length(unique(lr.x)) > 2 & length(table(sig.x)) > 1){ 
			k2 <- Ckmeans.1d.dp(lr.x,2);
			k1 <- Ckmeans.1d.dp(lr.x,1);
			ratio <- k1$withinss/sum(k2$withinss);
									
			if(ratio >= 3.5){
				clus.gain <- which(k2$centers == max(k2$centers));
				clus.not_gain <- which(k2$centers != max(k2$centers));
				clus.diff <- k2$centers[clus.gain] - k2$centers[clus.not_gain];
				
				clus.gain.members <- as.numeric(k2$cluster == clus.gain);
				runsum <- 0;
				runsum_arr <- c();
				for(i in seq(1,length(clus.gain.members),1)){
					if(clus.gain.members[i]==0){
						runsum <- 0;
					}else{
						runsum <- runsum + clus.gain.members[i];
					}
					runsum_arr <- c(runsum_arr,runsum);
				}
				max.runsum <- max(runsum_arr);
				not.in.stretch <- k2$size[clus.gain] - max.runsum;
				pct.stretch <- max.runsum/k2$size[clus.gain];
		
				if((k2$size[clus.gain] == 1 & clus.diff >= 0.5) |
				   (k2$size[clus.gain] > 1 & k2$size[clus.gain] < 6  & clus.diff >= 0.35 & pct.stretch > 0.5) |
				   (k2$size[clus.gain] >=6 & clus.diff >= 0.35 & pct.stretch >= 0.5)){								
				   	
				   	cluster <- k2$cluster;
				   	names(cluster) <- names(lr.x);
				   	probes <- names(genelabel)[which(flag)];
				   	lr.x.i <- logratio[probes];
				   	fc.x.i <- fc[names(genelabel)[which(flag)]];
				   	p.adj.x.i <- p.adj[names(genelabel)[which(flag)]];
				   	p.cluster <- rep(NA, length(lr.x.i));				   	
				   	names(p.cluster) <- probes;
				   	p.cluster[names(cluster)] <- cluster;
				   	ann <- do.call('rbind',lapply(names(lr.x.i),function(targ.i){
						if(targ.i %in% gene.ann[,'Target']){
							return(gene.ann[which(gene.ann[,'Target']==targ.i),c('Cyt','GeneExon')]);
						}else{
							return(c('Cyt'=NA,'GeneExon'=NA));				
						}
					}));
					ret <- rbind(ret,cbind('table'=rep('Intragenic_gain',length(lr.x.i)),'lr'=lr.x.i,'NExons'=rep(1,length(lr.x.i)),
					    	                        'chr'=do.call('c',lapply(names(lr.x.i),function(x){
									f <-  unlist(strsplit(x,'\\:'));
									return(f[1]);
			    				 })),'fc'=fc.x.i,'p.adj'=p.adj.x.i,'sig'=as.character(p.cluster),'sig.intragenic'=rep('-',length(lr.x.i)),ann));
			    }				
			}
		}
		return(ret);
	}));

	ann <- do.call('rbind',lapply(names(logratio),function(targ.i){
		if(targ.i %in% gene.ann[,'Target']){
			return(gene.ann[which(gene.ann[,'Target']==targ.i),c('Cyt','GeneExon')]);
		}else if(targ.i %in% tiling.ann[,'Target']){
			return(c('Cyt'=tiling.ann[which(gene.ann[,'Target']==targ.i),'Cyt'],'GeneExon'='Tiling'));
		}else{
			stop(paste('Error: Annotation not found ',targ.i,"\n")); 
		}
	}));				
	ret.probes <- cbind('table'=rep('Probe',length(logratio)),'lr'=logratio,'NExons'=rep(1,length(logratio)),
                            'chr'=do.call('c',lapply(names(logratio),function(x){
				f <- unlist(strsplit(x,'\\:'));
				return(f[1]);
			    })),'fc'=fc,'p.adj'=p.adj,'sig'=sig);

	sig.intragenic <- rep(0,nrow(ret.probes));		
 	names(sig.intragenic) <- rownames(ret.probes);
	if(!is.null(ret.intragenic)){
		sig.intragenic[rownames(ret.intragenic)[which(is.na(ret.intragenic[,'sig']))]] <- NA;
 		sig.intragenic[rownames(ret.intragenic)[which(ret.intragenic[,'sig']==1)]] <- 1;
 		sig.intragenic[rownames(ret.intragenic)[which(ret.intragenic[,'sig']==2)]] <- 2;
 	}
 	ret.probes <- cbind(ret.probes,'sig.intragenic'=sig.intragenic,ann);
 	
	gene.cyt <- do.call('c',lapply(rownames(seg.genes),function(gene){
		cyt <- gene.ann[which(gene.ann[,'gene'] %in% gene),'Cyt'];
		cyt.table <- table(cyt);
		cyt <- names(cyt.table)[which(cyt.table == max(cyt.table))];
		f <- unlist(strsplit(cyt,'\\.'));		
		return(f[1]);
	}));
 	ann <- cbind('Cyt'=gene.cyt, 'GeneExon'=rep(NA,nrow(seg.genes))); 				        
   	ret.genes <- cbind('table'=rep('Gene',nrow(seg.genes)),'lr'=seg.genes[,'seg'],'NExons'=seg.genes[,'num.exons'],
                       'chr'=seg.genes[,'chr'],'fc'=fc.genes,
			           'p.adj'=pg.adj,'sig'=sig.g);
	sig.intragenic <- rep(0,nrow(ret.genes));		
 	names(sig.intragenic) <- rownames(ret.genes);
	if(!is.null(ret.intragenic)){
		intragenic.genes <- do.call('c',lapply(as.character(genelabel[rownames(ret.intragenic)]),function(x){
			f <- unlist(strsplit(x,'\\,'));
			return(f);
		}));			   
		sig.intragenic[unique(intragenic.genes)] <- 1;
	}
	ret.genes <- cbind(ret.genes,'sig.intragenic'=sig.intragenic,ann);	

	if(is.null(ret.intragenic)){			   
		ret <- rbind(ret.probes,ret.genes);
	}else{
		ret <- rbind(ret.probes,ret.genes,ret.intragenic);
	}
	return(ret);
}

# Compute log ratio given tumor, selected_normal, matched_normal	
cn.analysis <- function(tm.vec, auto_nm.vec, x_nm.vec, tm.label, bc.label, genelabel,make.plot=F){
	# exclude Y
	tm.vec <- tm.vec[which(!grepl('^Y\\:',names(tm.vec)))];
	auto_nm.vec <- auto_nm.vec[which(!grepl('^Y\\:',names(auto_nm.vec)))];
	x_nm.vec <- x_nm.vec[which(!grepl('^Y\\:',names(x_nm.vec)))];

	x.idx <- grep('^X\\:',names(tm.vec));
	lr <- log2(tm.vec) - log2(auto_nm.vec);
	lr[x.idx] <- log2(tm.vec)[x.idx] - log2(x_nm.vec)[x.idx];
	return(seg.clus.p(lr,genelabel,make.plots=make.plot,sample.id=tm.label,bc.label));
}
	
if(doFull == "FULL"){	
    out_file= paste(prefix,"_copynumber_segclusp.nvn.full.pdf", sep="");
}else{
    out_file= paste(prefix,"_copynumber_segclusp.nvn.pdf", sep="");
}
pdf(out_file, height=8.5, width=11)
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(11, 10, "Memorial Sloan Kettering",adj=1)
text(11, 9.5, "Molecular Diagnostics Service",adj=1)
text(11, 9, date(),adj=1)
text(1, 8, "Copy number segmentation analysis - NVN",adj=0,cex=1.5,font=2)
text(1,7.5, paste("Project pool: ", title$Pool[1],sep=""),adj=0,cex=1.5,font=2) #Add project name from title file data
text(1,7, paste("Baits version: ", title$Bait_version[1],sep=""),adj=0,cex=1.5,font=2) #Add baits version from title file data
	
test.samples <- colnames(test);

analysis.out <- do.call('rbind',lapply(test.samples,function(tst){
	tst.probes <- test[,tst];
	names(tst.probes) <- genelevel[,'Target'];
	genelabel <- genelevel[,'genes'];
	names(genelabel) <- genelevel[,'Target'];
	test.bc <- title[which(title[,'Sample_ID'] == tst & title[,'Class'] %in% tested.class),'Barcode'];
	
	cat(tst,"\n");
	
	cn.probes <- do.call('cbind',lapply(colnames(d.stdn),function(snm){		
		nm.probes.i <- d.stdn[,snm];
		ret.probes <- calc.fc(tst.probes,nm.probes.i);
		colnames(ret.probes) <- paste(snm,"|",colnames(ret.probes),sep="");
		return(ret.probes);
	}));
	cn.probes <- cbind(genelevel[,1:3],cn.probes);

	cn.probes.autosomes <- cn.probes[which(!grepl('^X\\:',cn.probes[,'Target']) & !grepl('^Y\\:',cn.probes[,'Target'])),];
	auto.noise <- do.call('c',list(apply(cn.probes.autosomes[,grep('lr',colnames(cn.probes.autosomes)),drop=F],2,function(lr){
		if(sum(is.na(lr)) > 100){return(NA);} 
		noise <- as.numeric(lr[-1] - lr[-length(lr)]);	
		#return(sum(noise^2,na.rm=T));
		return(sum(lr^2,na.rm=T));
	})));
	names(auto.noise) <- gsub('\\|lr','',names(auto.noise));
	auto.noise <- auto.noise[which(!is.na(auto.noise))];
	best.norm.auto <- names(auto.noise)[which(auto.noise == min(auto.noise))][1];		
	best.nm.probes.auto <- d.stdn[,best.norm.auto];
	names(best.nm.probes.auto) <- genelevel[,'Target'];
	
	y.cov.pt <- median(as.numeric(tst.probes[grep('^Y\\:',genelevel[,'Target'])]),na.rm=T);
	if(gender[tst] != '-'){
		pt.gender <- gender[tst];
	}else{
		if(y.cov.pt < mean(y.cov_grp)){
			pt.gender <- "F";
		}else{
			pt.gender <- "M";
		}
	}
	cn.probes.X <- cn.probes[which(grepl('^X\\:',cn.probes[,'Target'])),grep('lr',colnames(cn.probes)),drop=F];
	colnames(cn.probes.X) <- gsub('\\|lr','',colnames(cn.probes.X));
	cn.probes.X <- cn.probes.X[,which(colnames(cn.probes.X) %in% c(paste('STDN_',stdn_cov.title[which(stdn_cov.title[,'Sample_type'] == pt.gender),'Barcode'],sep='')))];				
	X.noise <- do.call('c',list(apply(cn.probes.X,2,function(lr){
		if(sum(is.na(lr)) > 100){return(NA);} 
		noise <- as.numeric(lr[-1] - lr[-length(lr)]);	
		#return(sum(noise^2,na.rm=T));
		return(sum(lr^2,na.rm=T));
	})));
	X.noise <- X.noise[which(!is.na(X.noise))];
	best.norm.X <- names(X.noise)[which(X.noise == min(X.noise))][1];		
	best.nm.probes.X <- d.stdn[,best.norm.X];
	names(best.nm.probes.X) <- genelevel[,'Target'];
		
	# Best-normal check			
	ratiolabel <- paste(tst," vs. ",best.norm.auto,"\nSex=",pt.gender,", Norm for X=",best.norm.X,sep="");

    sig.list <- cn.analysis(tst.probes, best.nm.probes.auto, best.nm.probes.X, tst, test.bc, genelabel,make.plot=F);
	sig.list[,'table'] <- as.character(sig.list[,'table']);
	sig.list[,'lr'] <- as.numeric(as.character(sig.list[,'lr']));
	sig.list[,'NExons'] <- as.numeric(as.character(sig.list[,'NExons']));
	sig.list[,'fc'] <- as.numeric(as.character(sig.list[,'fc']));
	sig.list[,'p.adj'] <- as.numeric(as.character(sig.list[,'p.adj']));
	sig.list[,'sig'] <- as.character(sig.list[,'sig']);
	sig.list[,'sig.intragenic'] <- as.character(sig.list[,'sig.intragenic']);
	sig.list[,'Cyt'] <- as.character(sig.list[,'Cyt']);
	sig.list[,'GeneExon'] <- as.character(sig.list[,'GeneExon']);
	
	sig.probes <- sig.list[which(grepl('\\:',rownames(sig.list)) & sig.list[,'table'] == 'Probe'),];
	sig.genes <- sig.list[which(!grepl('\\:',rownames(sig.list)) & sig.list[,'table'] == 'Gene'),];
		
	sig.probes <- cbind('genes'=genelevel[match(rownames(sig.probes),genelevel[,'Target']),'genes'],
	                     'loc'=rownames(sig.probes),
	                     sig.probes);	
	sig.genes <- cbind('genes'=rownames(sig.genes),sig.genes);
	sig.genes[,'chr'] <- sig.genes[,'Cyt'];	


	# Plot 1
	layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T),heights=c(1,1),widths=c(1,1));
	lr <- as.numeric(as.character(sig.probes[,'lr']));	
	a<-max(abs(lr),na.rm=TRUE)
	if (a > 4){ a=signif(a, digits = 1)}else {a= 4}
	lr.plot_tiling <- lr;
	lr.plot_exon <- lr;
	lr.plot_tiling[-which(cn.probes[,'Target'] %in% tiling_probe_ids)] <- NA;
	lr.plot_exon[which(cn.probes[,'Target'] %in% tiling_probe_ids)] <- NA;
	plot(lr.plot_tiling,xlab="", ylab="",main="",xlim=c(0,nrow(d.stdn)),ylim=c(-a,a),xaxt="n",pch=20,col="brown",mgp=c(2,1,0),mar=c(3,5,2,5));
	par(new=T);
	plot(lr.plot_exon,xlab="Targeted exons per chromosome", ylab="LOG2 Tumor/Normal Ratio",main=ratiolabel,xlim=c(0,nrow(d.stdn)),ylim=c(-a,a),xaxt="n",pch=20,col="navy",mgp=c(2,1,0),mar=c(3,5,2,5))
	par(new=F);
	abline(v=linepos,lty=3);
	abline(h=-1,lty=2);
	abline(h=1,lty=2);
	mtext(c(1:22,'X'),side=1,line=0,at=numberpos,cex=0.7);
	legend(x='topright',pch=20,col=c('navy','brown'),legend=c('Exonic probe','Tiling probe'));

	sig.genes <- sig.genes[order(as.numeric(sig.genes[,'fc']),decreasing=T),,drop=F];
	sig.gene.gains <- sig.genes;
	sig.gene.gains[,'fc'] <- formatC(as.numeric(sig.gene.gains[,'fc',drop=F][[1]]),digits=2,format='f');
	idx <- which(!grepl('Tiling',sig.gene.gains[,'genes']) & as.numeric(sig.gene.gains[,'fc']) >= 1.3 & as.numeric(sig.gene.gains[,'p.adj']) < 0.05);	   
	sig.gene.gains[idx,'fc'] <- paste(sig.gene.gains[idx,'fc'],"*",sep="");
	sig.gene.gains <- sig.gene.gains[which(!grepl('Tiling',sig.genes[,'genes']) & !grepl('FP_',sig.genes[,'genes'])),];
				
	sig.gene.gains[,'p.adj'] <- formatC(as.numeric(sig.gene.gains[,'p.adj',drop=F][[1]]),digits=1,format='e');
	sig.gene.gains <- sig.gene.gains[,c('chr','genes','NExons','fc','p.adj'),drop=F];
		
	if(nrow(sig.gene.gains)>20){
		textplot(sig.gene.gains[1:20,], show.rownames = FALSE, hadj=0,valign="top",cex=1,mar=c(2,2,2,2),main='Top 20 Gains (gene-level) ');
	}else{
		textplot(sig.gene.gains, show.rownames = FALSE, hadj=0,valign="top",cex=1,mar=c(2,2,2,2),main='Top 20 Gains (gene-level) ');
	}
	
	sig.genes <- sig.genes[order(as.numeric(sig.genes[,'fc']),decreasing=F),,drop=F];	
	sig.gene.losses <- sig.genes;
	sig.gene.losses[,'fc'] <- formatC(as.numeric(sig.gene.losses[,'fc',drop=F][[1]]),digits=2,format='f');
	idx <- which(!grepl('Tiling',sig.gene.losses[,'genes']) & as.numeric(sig.gene.losses[,'fc']) <= -1.7 & as.numeric(sig.gene.losses[,'p.adj']) < 0.05);	   
	sig.gene.losses[idx,'fc'] <- paste(sig.gene.losses[idx,'fc'],"*",sep="");
	sig.gene.losses <- sig.gene.losses[which(!grepl('Tiling',sig.genes[,'genes']) & !grepl('FP_',sig.genes[,'genes'])),];

	sig.gene.losses[,'p.adj'] <- formatC(as.numeric(sig.gene.losses[,'p.adj',drop=F][[1]]),digits=1,format='e');
	sig.gene.losses <- sig.gene.losses[,c('chr','genes','NExons','fc','p.adj'),drop=F];
	
	if(nrow(sig.gene.losses)>20){
		textplot(sig.gene.losses[1:20,], show.rownames = FALSE, hadj=0,valign="top",cex=1,mar=c(2,2,2,2),main='Top 20 Losses (gene-level) ');
	}else{
		textplot(sig.gene.losses, show.rownames = FALSE, hadj=0,valign="top",cex=1,mar=c(2,2,2,2),main='Top 20 Losses (gene-level) ');
	}	
	
	tiling <- rep('Gene',nrow(sig.probes));
	tiling[grep('Tiling',sig.probes[,'GeneExon'])] <- 'Tiling';
	gene.refseq.exon.boundaries <- do.call('rbind',lapply(as.character(sig.probes[,'GeneExon']),function(gene_exon){
	  	f <- unlist(strsplit(gene_exon,'\\:'));
	  	if(grepl('Tiling',gene_exon) | gene_exon == ""){
	  		return(c('Gene'='','Refseq'='','Exon'='','Boundaries'=''));
	  	}else{	  			
	  		return(c('Gene'=f[1],'Refseq'=f[2],'Exon'=f[3],'Boundaries'=paste(f[4:length(f)],collapse=':')));
	  	}
	}));
	
	ret <- cbind('sample'=rep(tst,nrow(sig.list)),
                 'norm_used_auto'=rep(best.norm.auto,nrow(sig.list)),
                 'norm_used_X'=rep(best.norm.X,nrow(sig.list)),
		     'region'=rownames(sig.list),
		     sig.list);
	return(ret);
}));
dev.off();

probe.table <- analysis.out[which(analysis.out[,'table'] == 'Probe'),];
write.table(probe.table[,-which(colnames(probe.table) == 'table')],paste(prefix,"_copynumber_segclusp.nvn.probes.txt",sep=''),sep='\t',row.names=F,col.names=T,quote=F);
	  
gene.table <- analysis.out[which(analysis.out[,'table'] == 'Gene'),];
write.table(gene.table[,-which(colnames(gene.table) == 'table')],paste(prefix,"_copynumber_segclusp.nvn.genes.txt",sep=''),sep='\t',row.names=F,col.names=T,quote=F);	  
	  
intragenic.table <- c();
idx <- which(analysis.out[,'table'] == 'Intragenic_loss');
if(length(idx) > 0){
	intragenic.table <- analysis.out[idx,];
}
idx <- which(analysis.out[,'table'] == 'Intragenic_gain');
if(length(idx) > 0){
	intragenic.table <- rbind(intragenic.table,analysis.out[idx,]);
}
if(!is.null(intragenic.table)){
	write.table(intragenic.table[,-which(colnames(intragenic.table) == 'table')],paste(prefix,"_copynumber_segclusp.nvn.intragenic.txt",sep=''),sep='\t',row.names=F,col.names=T,quote=F);
}	




