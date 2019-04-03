.libPaths("/gpfs/projects/bsc05/ramon/R_libs")

library(data.table)
library(plyr)
library(dplyr)
library(reshape)
library(IRanges)

args <- commandArgs(TRUE)

tophits <- args[1]
output_crosspheno <- args[2]
output_ranges <- args[3]
output_summary <- args[4]
pval <- as.numeric(args[5])

tophits <- unlist(strsplit(args[1],","))

import.list <- llply(tophits, read.delim, sep="", colClasses=c("character","numeric",
                                                             "character",
                                                             "character","character",
                                                             "numeric","numeric","numeric",
                                                             "numeric","numeric"))

cross_pheno_all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, 
													 by = c("chr","position",
															"rs_id_all","alleleA","alleleB",
															"info_all","all_maf"), 
													 all = TRUE),import.list)

write.table(cross_pheno_all,output_crosspheno,col.names=T,row.names=F,quote=F,sep="\t")

cross_pheno_all$start <- cross_pheno_all$position - 250000
chr <- as.vector(unique(cross_pheno_all$chr))

ranges <- NULL
    for(i in chr){
        dat <- cross_pheno_all[cross_pheno_all$chr == i,]
        dat <- IRanges(start=dat$start,width=500000)
        dat <- reduce(dat)
        dat <- as.data.frame(dat)
        dat$chr <- i
        ranges <- rbind(ranges,dat)
    }

if (nrow(cross_pheno_all)!=0){
    ranges$Range <- paste("chr",ranges$chr,":",ranges$start,"-",ranges$end,sep="")
} else {
    ranges <- data.frame(matrix(ncol = 4, nrow = 0))
    x <- c("chr","start","end","width")
    colnames(ranges) <- x
}

ranges

write.table(ranges,output_ranges,col.names=T,row.names=F,quote=F,sep="\t")

cat("#############################################################################\n")
cat("#####                     CROSS-PHENOTYPE SUMMARY                       #####\n")
cat("#############################################################################\n")

diseases <- names(cross_pheno_all)[which(grepl("frequentist_add_pvalue_",names(cross_pheno_all)))]

if (nrow(cross_pheno_all)!=0){
    for (i in 1:nrow(cross_pheno_all)){
            r <- which(ranges$start <= cross_pheno_all$position[i] & cross_pheno_all$position[i] <= ranges$end & ranges$chr==cross_pheno_all$chr[i])
            cross_pheno_all$RANGE[i] <- paste("chr",ranges$chr[r],":",ranges$start[r],"-",ranges$end[r],sep="")
    }

    pval_threshold <- 0.05/length(unique(cross_pheno_all$RANGE))/length(diseases)

    crosspheno_variants <- NULL
    for (i in 1:nrow(cross_pheno_all)){
        for (n in diseases){
            if (cross_pheno_all[i,n]<=pval & !is.na(cross_pheno_all[i,n])){
                for (e in diseases){
                    if (e != n){
                    if (cross_pheno_all[i,e]<=pval_threshold & !is.na(cross_pheno_all[i,e])){
                        crosspheno_variants$range[i] <- cross_pheno_all$RANGE[i]
                        crosspheno_variants$rsid[i] <- as.character(cross_pheno_all$rs_id_all[i])
                        crosspheno_variants$disease_A[i] <- as.character(lapply(strsplit(n,"frequentist_add_pvalue_"), "[", 2))
                        crosspheno_variants$frequentist_add_pvalue_A[i] <- cross_pheno_all[i,n]
                        crosspheno_variants$disease_B[i] <- as.character(lapply(strsplit(e,"frequentist_add_pvalue_"), "[", 2))
                        crosspheno_variants$frequentist_add_pvalue_B[i] <- cross_pheno_all[i,e]}}}
}}}

    crosspheno <- as.data.frame(crosspheno_variants)
    crosspheno_2 <- crosspheno[!is.na(crosspheno$disease_B),]
    crosspheno <- crosspheno_2[order(crosspheno_2$range),]

    crosspheno$disease_A_vs_disease_B <- NA
    crosspheno$disease_A_vs_disease_B <- paste(crosspheno$disease_A,"-",crosspheno$disease_B,sep="")

    associations <- as.character(unique(crosspheno$disease_A_vs_disease_B))
    range <- as.character(unique(crosspheno$range))

    crosspheno_min_pval <- NULL
    for (i in 1:length(range)){
        for (n in 1:length(associations)){
             all <- crosspheno[crosspheno$range==range[i] & crosspheno$disease_A_vs_disease_B==associations[n],]
             crosspheno_min_pval <- rbind(crosspheno_min_pval,all[which.min(all$frequentist_add_pvalue_A),])}
    }

} else {
    crosspheno_min_pval <- data.frame(matrix(ncol = 6, nrow = 0))
    x <- c("range","rsid","disease_A","frequentist_add_pvalue_A","disease_B","frequentist_add_pvalue_B")
    colnames(crosspheno_min_pval) <- x
}

write.table(crosspheno_min_pval,output_summary,col.names=T,row.names=F,quote=F,sep="\t")

cat("#############################################################################\n")
cat("#####                  CROSS-PHENOTYPE SUMMARY DONE                     #####\n")
cat("#############################################################################\n")




