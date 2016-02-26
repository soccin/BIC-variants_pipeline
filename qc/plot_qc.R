library(corrplot)
library(ggplot2)
library(gplots)
library(scales)
library(reshape)
library(plyr)

library(RColorBrewer)

c3 <- c("#5395b4", "#f1592a", "#85c440")
c48 <- c("#1d915c","#5395b4",
                "#964a48",
                "#2e3b42",
                "#b14e72",
                "#402630","#f1592a",
                "#81aa90","#f79a70", # lt pink
                "#b5ddc2",
                "#8fcc8b", # lt purple
                "#9f1f63", # lt orange
                "#865444", "#a7a9ac",
                "#d0e088","#7c885c","#d22628","#343822","#231f20",
                "#f5ee31","#a99fce","#54525e","#b0accc",
                "#5e5b73","#efcd9f", "#68705d", "#f8f391", "#faf7b6", "#c4be5d", "#764c29", "#c7ac74", "#8fa7aa", "#c8e7dd", "#766a4d", "#e3a291", "#5d777a", "#299c39", "#4055a5", "#b96bac", "#d97646", "#cebb2d", "#bf1e2e", "#d89028", "#85c440", "#36c1ce", "#574a9e")

plot.fpc.sum <- function(fpc.sum){
    if(is.null(fpc.sum)){ return(NULL) }
    heatmap.2(fpc.sum,
              Rowv="NA",
              Colv="NA",
              dendrogram="none",
              distfun="dist",
              hclustfun="hclust",
              labRow=rownames(fpc.sum),
              labCol=rownames(fpc.sum),
              #key=TRUE,
              keysize=0.75,
              trace="none",
              density.info=c("none"),
              margins=c(12,12),
              col=brewer.pal(5,"Spectral"),
              main="Sample Mix-Ups",
              sepcolor="black",
              colsep=c(1:nrow(fpc.sum)),
              rowsep=c(1:nrow(fpc.sum)),
              #cexRow = 0.7,
              #cexCol = 0.7,
              cexRow = 1.0 - (0.01*(nrow(fpc.sum)-25)),
              cexCol = 1.0 - (0.01*(nrow(fpc.sum)-25)),
              sepwidth=c(0.005, 0.005)
              )
}

plot.cdna.contamination <- function(cdna.contamination){
    if(is.null(cdna.contamination)){ return(NULL) }
    corrplot(cdna.contamination,
             method="circle",
             col=c("white","red"),
             cl.pos="n",
             title="cDNA Contamination",
             pch.col="black", 
             tl.col="black",
             tl.srt=90,
             tl.cex=0.9,
             mar=c(0,5,5,5),
             is.corr=FALSE)
}

plot.major.contamination <- function(dat,sort.by='name'){
    if(is.null(dat)){ return(NULL) }
    colnames(dat) = c("Sample","PerHeterozygousPos") ## colnames are different between exome and dmp pipelines
    if(sort.by == 'badToGood'){
       dat$Sample <- factor(dat$Sample, levels=dat$Sample[order(dat$PerHeterozygousPos)])
    } else if(sort.by == 'name'){
        dat$Sample <- factor(dat$Sample, levels=dat$Sample[order(dat$Sample)])
    }

    ggplot(dat, aes(x = Sample, y = PerHeterozygousPos)) +
        geom_bar(stat="identity") +
        ylim(0, 1) +
        theme(legend.title=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1,size=12,color="black"),
            axis.text.y=element_text(size=12,color="black")
        ) +
        xlab("") +
        ylab("Fraction of position that are hetrozygous") +
        labs(title="Major Contamination Check") +
        geom_hline(aes(yintercept=0.55), color = "red", size=0.75) #+
        #coord_flip()
}

plot.minor.contamination <- function(dat,sort.by='name'){
    if(is.null(dat)){ return(NULL) }
    if(sort.by == 'badToGood'){
        dat$Sample <- factor(dat$Sample, levels=dat$Sample[order(dat$AvgMinorHomFreq)])
    } else if(sort.by == 'name'){
        dat$Sample <- factor(dat$Sample, levels=dat$Sample[order(dat$Sample)])
    }

    maximum<-1.5*(as.numeric(max(dat$AvgMinorHomFreq, na.rm=TRUE)))

    ggplot(dat, aes(x = Sample, y = AvgMinorHomFreq)) +
      geom_bar(stat="identity") +
      ylim(0, maximum) +
      theme(legend.title=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1,size=9,color="black"),
            axis.text.y=element_text(size=12,color="black")
        ) +
      xlab("") +
      ylab("Avg. Minor Allele Frequency at Homozygous Position") +
      labs(title="Minor Contamination Check") +
      geom_hline(aes(yintercept=0.02), color = "red", size=0.5) +
      geom_hline(aes(yintercept=0.01), color = "yellow", size=0.5) +
      scale_color_hue(name="Samples") #+
      #coord_flip()
}

plot.coverage <- function(dat,sort.by='name'){
    if(is.null(dat)){ return(NULL) }
    if(sort.by == "badToGood"){
        dat$Samples = factor(dat$Samples, levels=dat$Samples[order(-dat$Cov)])
    } else if(sort.by == "name"){
        dat$Samples = factor(dat$Samples, levels=dat$Samples[order(dat$Samples)])
    }
    ggplot(dat, aes(x = Samples,y = Cov)) +
       geom_bar(stat = "identity") +
       theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          legend.position="bottom",
          legend.title = element_blank()
        ) +
       labs(title="Mean Target Coverage") +
       xlab("") +
       ylab("Mean Coverage") +
       geom_hline(aes(yintercept=200),color="yellow", size=0.75) +
       geom_hline(aes(yintercept=50), color="red", size=0.75) #+
       #coord_flip()
}

plot.duplication <- function(duplication,sort.by='name'){
    if(is.null(duplication)){ return(NULL) }
    if(sort.by == 'badToGood'){
        duplication$Samples <- factor(duplication$Samples, levels=duplication$Samples[order(duplication$DupRate)])
    } else if(sort.by == 'name'){
        duplication$Samples <- factor(duplication$Samples, levels=duplication$Samples[order(duplication$Samples)])
    }
    ggplot(duplication, aes(x = Samples, y = DupRate)) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          legend.position="bottom",
          legend.title = element_blank()
        ) +
        labs(title="Estimated Duplication Rate") +
        xlab("") +
        ylab("Duplication Rate") +
        scale_y_continuous(labels=percent) +
        geom_hline(aes(yintercept=0.5),color="yellow", size=1.0) #+
        #coord_flip()
}

plot.library.size <- function(librarySize,sort.by='name'){
    if(is.null(librarySize)) { return(NULL) }
    if(sort.by == 'badToGood'){
        librarySize$Samples <- factor(librarySize$Samples, levels=librarySize$Samples[order(-librarySize$Comp)])
    } else if(sort.by == 'name'){
        librarySize$Samples <- factor(librarySize$Samples, levels=librarySize$Samples[order(librarySize$Samples)])
    }
    ggplot(librarySize, aes(x = Samples, y = Comp/1000000)) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
          axis.text.y=element_text(size=12,color="black")#,
          #legend.position="bottom",
          #legend.title = element_blank()
        ) +
        labs(title="Estimated Library Size") +
        xlab("") +
        ylab("xMillion") #+
        #coord_flip()
}

plot.capture.specificity <- function(captureSpecificity,sort.by='name'){
    if(is.null(captureSpecificity)){ return(NULL) }
    cs.m = melt(captureSpecificity, id.vars="Sample")

    if(sort.by == 'badToGood'){
        cs.m$Sample <- factor(cs.m$Sample, levels=cs.m$Sample[order(-cs.m$value)])
    } else if(sort.by == 'name'){
        cs.m$Sample <- factor(cs.m$Sample, levels=cs.m$Sample[order(cs.m$Sample)])
    }

    ggplot(cs.m, aes(x = Sample, y = value/1000000, fill = variable)) +
       geom_bar(stat="identity", width=0.7, color="black") +
       theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          legend.position="bottom",
          legend.title = element_blank()
       ) +
       scale_fill_brewer(palette="Set2",labels=c("On Bait Bases","Near Bait Bases","Off Bait Bases")) +
       labs(title="Capture Specificity") +
       xlab("") +
       ylab("Total Bases (millions)") #+
       #coord_flip()
}

plot.capture.specificity.percentage <- function(captureSpecificity,sort.by='name'){
    if(is.null(captureSpecificity)){ return(NULL) }
    cs = captureSpecificity
    x=cs[,c("OnBait","NearBait","OffBait")]
    on.bait.per <- cs$OnBait/apply(x,1,sum)
    on.bait.order <- order(-on.bait.per)
    cs.m <- melt(cs, id.vars="Sample")

    if(sort.by == 'badToGood'){
        cs.m$Sample <- factor(cs.m$Sample, levels=cs.m$Sample[on.bait.order])
    } else if(sort.by == 'name'){
        cs.m$Sample <- factor(cs.m$Sample, levels=cs.m$Sample[order(cs.m$Sample)])
    }

    ggplot(cs.m, aes(x = Sample, y = value/1000000, fill = variable)) +
        geom_bar(stat="identity", position="fill",width=0.7, color="black") +
        theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          legend.position="bottom",
          legend.title = element_blank()
        ) +
        scale_fill_brewer(palette="Set2", labels = c("On Bait Bases", "Near Bait Bases", "Off Bait Bases")) +
        labs(title="Capture Specificity") +
        xlab("") +
        ylab("Percent bases") +
        scale_y_continuous(labels=percent) #+
        #coord_flip()
}

plot.alignment <- function(alignment,sort.by='name'){
    if(is.null(alignment)){ return(NULL) }
    density.m = melt(alignment)
    density.m$value <- density.m$value/1000000

    if(sort.by == 'badToGood'){
        density.m$Samples <- factor(density.m$Samples, levels=density.m$Samples[order(-density.m$value)])
    } else if(sort.by == 'name'){
        density.m$Samples <- factor(density.m$Samples, levels=density.m$Samples[order(density.m$Samples)])
    }

    ggplot(density.m, aes(x = Samples, y = value, fill = variable)) +
        geom_bar(stat="identity", width=0.7, color="black")+
        theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          legend.position="bottom",
          legend.title = element_blank()
        ) +
        scale_fill_manual(name="Type", values = c3, labels = c("Both Reads Aligned", "One Read Aligned", "Neither Read Aligned")) +
        labs(title="Cluster Density & Alignment Rate") +
        xlab("") +
        ylab("xMillion") #+
        #coord_flip()
}

plot.alignment.percentage <- function(alignment,sort.by='name'){
    density = alignment
    if(is.null(density)){ return(NULL) }
    x=density[,c("BothAlign","OneAlign","NeitherAlign")]
    both.per <- density$BothAlign/apply(x,1,sum)
    both.order <- order(-both.per)
    density.m <- melt(density)
    density.m$value <- density.m$value/1000000
    density.m2 <- ddply(density.m, .(Samples), mutate, perc = value/sum(value))

    y <- min(density.m2[density.m2$variable == "BothAlign", ]["perc"])
    y <- y*0.90

    if(sort.by == 'badToGood'){
        density.m2$Samples <- factor(density.m2$Samples, levels=density.m$Samples[both.order])
    } else if(sort.by == 'name'){
        density.m2$Samples <- factor(density.m2$Samples, levels=density.m2$Samples[order(density.m2$Samples)])
    }

    ggplot(density.m2, aes(x = Samples, y = perc, fill = variable)) +
        geom_bar(stat="identity", position="fill", width=0.7, color="black")+
        theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          legend.position="bottom",
          legend.title = element_blank()
        ) +
        scale_fill_manual(name="Type",values = c3, labels = c("Both Reads Aligned", "One Read Aligned", "Neither Read Aligned")) +
        labs(title="Cluster Density & Alignment Rate") +
        xlab("") +
        ylab("") +
        scale_y_continuous(label=percent, limits = c(y,1), oob =rescale_none) #+
        #coord_flip()
}

plot.insert.size.distribution <- function(is.metrics){
    if(is.null(is.metrics)){ return(NULL) }
    insert_label <- colnames(is.metrics)[2:ncol(is.metrics)]
    insert.m <- melt(is.metrics, id.vars="insert_size")
    insert.m$insert_size <- as.integer(as.character(insert.m$insert_size))
    levels(insert.m$variable) <- insert_label

    legend.position = "right"
    if(length(insert_label)>=20){ legend.position = "none" }
 
    ggplot(insert.m,
      aes(x = insert_size, y = value, color = variable)) +
      geom_line(size=0.75) +
      theme(legend.position=legend.position) + #"right",
            #legend.text = element_text(size=9),
            #legend.title = element_blank()
            #) +
      labs(title="Insert Size Distribution") +
      xlab("Insert size") +
      ylab("") +
      scale_color_manual(name="Samples", values = rep(c48, 34)) #+
      #guides(colour = guide_legend(override.aes = list(size=5), ncol=ceiling(ncol(is.metrics)/20)))
}

plot.insert.peaks <- function(is.metrics){
    if(is.null(is.metrics)){ return(NULL) }
    peaks <- as.data.frame(t(apply(is.metrics[,2:ncol(is.metrics)],2,which.max)))
    peaks_label <- colnames(peaks)
    peaks.m <- melt(peaks)
    levels(peaks.m$variable) <- peaks_label

    ggplot(peaks.m, aes(x = variable, y = value)) +
        geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle=45,size=12,hjust=1,color="black"),
              axis.text.y = element_text(size=12,color="black")) +
        xlab("") +
        ylab("Insert Size") +
        labs(title="Peak Insert Size Values") #+
        #coord_flip()
}

plot.trimmed.reads <- function(reads,sort.by='name'){
    reads.m <- melt(reads)
    reads.m$value <- reads.m$value/100

    if(sort.by == 'badToGood'){
        reads.m$Samples <- factor(reads.m$Samples, levels=reads.m$Samples[order(reads.m$value)])
    } else if(sort.by == 'name'){
        reads.m$Samples <- factor(reads.m$Samples, levels=reads.m$Samples[order(reads.m$Samples)])
    }


    ggplot(reads.m, aes(x = Samples, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = position_dodge(width=0.7), width = 0.7, color = "black") +
      theme(legend.title=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1,size=12,color="black"),
            axis.text.y=element_text(size=12,color="black")
           ) +
      scale_fill_brewer(palette="Set2") +
      labs(title="Percentage of Reads Trimmed") +
      xlab("") +
      ylab("") +
      scale_y_continuous(labels=percent) #+
      #coord_flip()
}

plot.base.qualities <- function(base.qualities){
    ggplot(base.qualities, aes(x = cycle, y = value, color = variable)) +
        geom_line(size=0.5) +
        theme(legend.position="none") + #"right",
              #legend.title=element_blank(),
              #axis.text.x = element_text(angle=0)
             #) +
        facet_wrap(~type, ncol = 1) +
        scale_color_manual(name="Samples", values = rep(c48, 34)) +
        labs(title="Pre- & Post-Recalibration Quality Scores") +
        xlab("Quality Score") +
        ylab("") #+
        #guides(colour = guide_legend(override.aes = list(size=2),ncol=ceiling(ncol(base.qualities)/20)))
}

plot.pool.norm.genotype <- function(pool.norm.genotype){
    if(is.null(pool.norm.genotype)){ return(NULL) }
    poolNormal = pool.norm.genotype
    corTest <- with(poolNormal, cor.test(ExpectedVF, ObservedVF))
    ggplot(poolNormal, aes(x = ExpectedVF, y = ObservedVF)) +
        geom_point(size=2, alpha=0.7) +
        theme(legend.position="right",
              legend.title=element_blank(),
              axis.text.x = element_text(angle=0)
             ) +
        labs(title=paste("Expected and Observed variant frequencies\n for SNPs present in pool normal. Cor:", round(corTest$estimate, digits=3), sep = " "))

}


plot.gc.bias <- function(gc.bias){
    if(is.null(gc.bias)){ return(NULL) }
    xt.m <- melt(gc.bias, id.vars="X")

    legend.position = "right"
    if(ncol(gc.bias)>=20){ legend.position = "none" }

    ggplot(xt.m, aes(x = X, y = value, color = variable)) +
       geom_line(size=0.5) +
       theme(legend.position=legend.position) + #"right",
              #legend.title=element_blank(),
              #axis.text.x = element_text(angle=0)
             #) +
       labs(title="Normalized Coverage vs GC-Content") +
       xlab("GC-Content") +
       ylab("Normalized Coverage") +
       scale_color_manual(name="Samples", values = rep(c48, 34)) +
       guides(colour = guide_legend(override.aes = list(size=5), ncol=ceiling(ncol(gc.bias)/20)))

}

