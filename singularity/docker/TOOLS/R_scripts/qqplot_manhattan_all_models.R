##
##  Copyright 2002-2014 Barcelona Supercomputing Center (www.bsc.es)
##  Life Science Department, 
##  Computational Genomics Group (http://www.bsc.es/life-sciences/computational-genomics)
##
##  Licensed under the Apache License, Version 2.0 (the "License");
##  you may not use this file except in compliance with the License.
##  You may obtain a copy of the License at
##
##      http://www.apache.org/licenses/LICENSE-2.0
##
##  Unless required by applicable law or agreed to in writing, software
##  distributed under the License is distributed on an "AS IS" BASIS,
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##  See the License for the specific language governing permissions and
##  limitations under the License.
##
## Last update: $LastChangedDate: 2015-01-08 11:04:27 +0100 (Thu, 08 Jan 2015) $
## Revision Number: $Revision: 14 $
## Last revision  : $LastChangedRevision: 14 $
## Written by     : XXXXX YYYYYY.
##                  xxxxx.yyyyy@gmail.com
## Modified by    : Friman Sanchez C.
##               friman.sanchez@gmail.com
## GWImp-COMPSs web page: http://cg.bsc.es/gwimp-compss/
##

qq.plot <- function(tab,lambda=T,stat,BA=F,plot=T,mod,
                    pch.col="red",max.yaxis=10,max.xaxis=10,scale.cex=T,scale.fact=0.25,dens=T,...) {

    #Calculating LAMBDA
    chi_list <- qchisq(na.omit(tab[,col_number]),df=1,lower.tail=F);
    l.fac <- median(chi_list)/0.456;
    print(paste(c("LAMBDA = ",l.fac),collapse=""));

    #Setting format for plot
    pvals <- -log10((na.omit(tab[,col_number])));
    tab$col <- "ivory3"
    tab$col[tab[,col_number] <= pval_threshold] <- "deeppink"
    col <- tab$col[!is.na(tab[,col_number])]
    col <- col[order(pvals)]
    pvals <- pvals[order(pvals)]
    tab$log10 <- pvals
    n <- length(pvals)

    #Parameters for density plot
    frac <- 0.01;
    alpha <- 0.05
    index <- seq(1,n)
    nsmall <- round(n*frac)
    temp <- log(seq(2,n-nsmall,100))
    indlow <- round(temp/max(temp)*(n-nsmall))
    index <- c(indlow,seq(n-nsmall,n))
    A <- qbeta(alpha/2,shape1=index,shape2=n-index+1,lower.tail=T)
    B <- qbeta(alpha/2,shape1=index,shape2=n-index+1,lower.tail=F)
    lower <- qexp(A,rate=log(10))
    upper <- qexp(B,rate=log(10))
    p.exp <- sort(-log10( c(1:length(pvals))/(length(pvals)+1) ))
    
    #Scaling Size dots
    scale.fact <- 0.25
    my.cex <- scale.fact * pvals * 0.9
    ylab <- expression(-log[10]~italic(p)[Obs])
    xlab <- expression(-log[10]~italic(p)[Exp])
    max.xaxis <- max(p.exp) + 1
    max.yaxis <- max(pvals) + 6
    
    #Creating Plot
    plot(p.exp,pvals, xlab=xlab, ylab=ylab,
           xlim=c(0,max.xaxis), ylim=c(0,max.yaxis), type="n", xaxs="i", yaxs="i",bty="l")
    #plot the confidence band for the data
    smallexp <- p.exp[index]
    y <- c(rev(lower),upper)
    polygon(c(rev(smallexp),smallexp),y,col="grey80",border=NA)

    ## rescaled density plots in bg
    d.vals <- density(pvals)
    lines(d.vals$x, (max.yaxis/2)+((max.yaxis/2)*(d.vals$y/max(d.vals$y))),
 	   col="cadetblue4",lty="dotted",lwd=2)
    # add density axis
    axis(side=4, at=c(max.yaxis/2,max.yaxis), labels=c(0,round(max(d.vals$y),1)))
    mtext("Data density",side=4,at=max.yaxis*0.75,line=1,adj=0.5,cex=0.9)
    lines(c(0,-log10(1/length(pvals))),c(0,-log10(1/length(pvals))), col="lightslategrey", lwd=2,lty=2)
    points(p.exp[col == "deeppink"], pvals[col == "deeppink"], pch=18, cex=my.cex[col == "deeppink"], col="deeppink")
    points(p.exp[col == "ivory3"], pvals[col=="ivory3"], pch=18, cex=my.cex[col=="ivory3"], col="ivory3")

    #legend files

    mtext(substitute(paste(lambda~" = "~lfac),list(lfac=round(l.fac,digits=3))),line=1,adj=0.5,cex=0.9)
}


manhattan <- function(dataframe, colors=c("lightsteelblue4","lightyellow2"), 
                      ymin=0, ymax="max", limitchromosomes=1:25, 
                      suggestiveline=-log10(pval_threshold/0.005), 
                      genomewideline=-log10(pval_threshold), 
                      annotate=NULL, ...) {
   d=dataframe
   #Validate input files
   if (!("chr" %in% names(d) & "position" %in% names(d) & "pvalue" %in% names(d))) stop("Make sure your data frame contains columns chr, position, and pvalue")

   if (any(limitchromosomes)) d=d[d$chr %in% limitchromosomes, ]
   d=subset(na.omit(d[order(d$chr, d$position), ]), (pvalue>0 & pvalue<=1)) # remove na's, sort, and keep only 0<pvalue<=1
   d$logp = -log10(d$pvalue)
   d$pos=NA
   ticks=NULL
   lastbase=0
   colors <- c(rep(colors,25)[1:25],"lightsteelblue4")

   if (ymax=="max") {ymax<-ceiling(max(d$logp))}
   	ymax<-as.numeric(ymax)
        if (ymax < 8) {ymax<-8}

        #Print Max and Min P-values
        print(paste("YMAX:",ymax))
        print(paste("YMIN:",ymin))

        numchroms = length(unique(d$chr))

        if (numchroms==1) {
                d$pos=d$position
                ticks=floor(length(d$pos))/2+1
        } else {
                for (i in unique(d$chr)) {
                        if (i==as.vector(unique(d$chr))[1]) {
                                d[d$chr==i, ]$pos=d[d$chr==i, ]$position
                        } else {
                                lastbase=lastbase+tail(subset(d,chr==as.vector(unique(d$chr))[which(as.vector(unique(d$chr))==i)-1])$position, 1)
                                d[d$chr==i, ]$pos=d[d$chr==i, ]$position+lastbase
                        }
                        ticks=c(ticks, d[d$chr==i, ]$pos[floor(length(d[d$chr==i, ]$pos)/2)+1])
                }
        }

		chr <- unique(d$chr)
		if (23 %in% chr){ 	
			ticks[(length(ticks))]
			ticks[(length(ticks)-2)]
			ticks2 <- c(ticks[(length(ticks))],ticks[(length(ticks)-2)])

		    d$chr[d$chr=="23"] <- as.character("F")		
		    d$chr[d$chr=="24"] <- as.character("M")
            d$chr[d$chr=="25"] <- as.character("A")
		}

        if (numchroms==1) {
                plot(d$pos, d$logp, ylim=c(ymin,ymax), ylab=expression(-log[10](italic(p))),
                        xlab=paste("Chromosome",unique(d$chr),"position"))
        }else {
                plot(d$pos,d$logp, ylim=c(ymin,ymax), ylab=expression(-log[10](italic(p))),
                        xlab="", xaxt="n", type="n")
				mtext(text = "Chromosome",
      				side = 1,
      				line = 4)
			axis(1, at=ticks, lab=unique(d$chr))
				if (23 %in% chr){ 
                	axis(1, at=ticks2, line=2.5,tick=T,lab=rep("",2),lwd.ticks=0)
					axis(1, at=(as.numeric(ticks2[1])+as.numeric(ticks2[2]))/2, lab=c("X"),line=2,tick=F)

				}

                icol=1

                for (i in unique(d$chr)) {
                #with(d[d$chr==i, ],points(pos[d$pvalue > 5e-8], logp[d$pvalue > 5e-8], col=colors[icol], ...))
                        points(d$pos[d$chr == i], d$logp[d$chr == i], col=colors[icol],pch=".")
                        icol=icol+1
                }
                for (i in unique(d$chr)) {
                        #with(d[d$chr==i, ],points(pos[d$pvalue > 5e-8], logp[d$pvalue > 5e-8], col=colors[icol], ...))
                        points(d$pos[d$chr == i & d$logp >= -log10(pval_threshold)], 
                        d$logp[d$chr == i & d$logp >= -log10(pval_threshold)], 
                        col="lightcoral",cex=0.6,pch=16)
                }
        }


        genomewideline <- -log10(pval_threshold)
        abline(h=genomewideline, col="indianred2", lty=2)


}

################################################################################################################
library(gap)
library(sfsmisc)

	args <- commandArgs(TRUE)

	tab_file <- args[1] #tab file with chr, position , pvalue
	out_qqplot <- args[2]
	out_manhattan <- args[3]
	out_qqplot_tiff <- args[4]
	out_manhattan_tiff <- args[5]
	inh_model <- args[6]
	pval_threshold <- as.numeric(args[7])

    tab_file_data<-read.delim(tab_file,na.strings="NA")
	
	col_number <- which(colnames(tab_file_data)==inh_model)

    tab_file_data <- tab_file_data[tab_file_data[,col_number]!=0,]
    tab_file_data <- tab_file_data[!is.na(tab_file_data[,col_number]),]
	tab_file_data_qq <- tab_file_data[tab_file_data$chr!="23_males" & tab_file_data$chr!="23_females",] 

    p <- tab_file_data_qq[,col_number]
    p <- p[!is.na(p)]
    n <- length(p)
    x2obs <- qchisq(p, 1, lower.tail = FALSE)
    x2exp <- qchisq(1:n/n, 1, lower.tail = FALSE)

#   makeqqplot
    pdf(out_qqplot)
        qq.plot(tab_file_data, lambda=F, stat="pvalue" ,scale.cex=T);

    dev.off()
    tiff(out_qqplot_tiff,width=5600,height=5600,units = "px", res = 800,compression="zip")
            qq.plot(tab_file_data, lambda=F, stat="pvalue" ,scale.cex=T);
    dev.off()
    cat("Q-Q Plot Assoc successfully completed!\n")
    
    #make Manhattan on pdf

    tab_file_data$chr <- as.character(tab_file_data$chr)
    tab_file_data$chr[tab_file_data$chr=="23"] <- as.character("25")
    tab_file_data$chr[tab_file_data$chr=="23_females"] <- as.character("23")
    tab_file_data$chr[tab_file_data$chr=="23_males"] <- as.character("24")

    tab_file_data$chr <- as.numeric(tab_file_data$chr)

    tab_file_data$SNP <- paste(paste("chr",tab_file_data$chr,sep=""),tab_file_data$position,sep=":")
    tab_man <- data.frame(SNP=tab_file_data$SNP,
                          chr=tab_file_data$chr,
                          position=tab_file_data$position,
                          pvalue=tab_file_data[,col_number])

    tab_man <- tab_man[tab_man$pvalue <= 0.05,]
	dim(tab_man)

    title <- "Manhattan-plot"
    YMIN <- 0
    YMAX <- "max"

    pdf(out_manhattan,width=24,height=14.5)
        manhattan(tab_man, pch=16,cex=0.7,main=title,colors=c("lightsteelblue4","ivory3"), suggestiveline=-log10(pval_threshold), ymax=YMAX, ymin=YMIN)
    dev.off()

    #make Manhattan on tiff
    tiff(out_manhattan_tiff,width=27,height=17.5,units = "in", res = 800,compression="zip")
        manhattan(tab_man, pch=16,cex=0.70,main=title,colors=c("lightsteelblue4","ivory3"), suggestiveline=-log10(pval_threshold), ymax=YMAX, ymin=YMIN)
    dev.off()

#    write.table(tab_man,out_corrected_pvals,row.names=F,sep="\t",quote=F)

    cat("MANHATTAN PLOT successfully completed!\n")
