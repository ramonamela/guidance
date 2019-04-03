.libPaths("/gpfs/projects/bsc05/ramon/R_libs")

library(data.table)
library(plyr)
library(dplyr)
library(reshape)


args <- commandArgs(TRUE)

cat("#############################################################################\n")
cat("#####                           TOP HITS LIST                           #####\n")
cat("#############################################################################\n")


tophits <- unlist(strsplit(args[1],","))
tophits_final <- args[2]

import.list <- llply(tophits, read.delim, sep="", colClasses=c("character","numeric",
                                                             "character","numeric",
                                                             "character","character",
                                                             "numeric","numeric",
															 "numeric","numeric"))
tophits_all <- rbindlist(import.list)

tophits_all <- as.data.frame(tophits_all)
head(tophits_all)

cat("#####                 Checking for duplicated variants				 #####\n")
table(duplicated(tophits_all$rs_id_all))

cat("#####                	Removing Duplicated Variant         	     #####\n")

tophits_all <- distinct(tophits_all, rs_id_all, .keep_all = TRUE)
table(duplicated(tophits_all$rs_id_all))

colnames(tophits_all)
head(tophits_all)
tophits_all_final <- tophits_all[,c("chr","position","rs_id_all","alleleA","alleleB")]

write.table(tophits_all_final,tophits_final,col.names=T,row.names=F,quote=F,sep="\t")

cat("#####            		  			DONE   				             #####\n")
