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
#models <- unlist(strsplit(args[3],","))

#classes <- c("character","numeric","character","numeric","character","character","numeric","charater")
#for (m in models){
#    if (m !="gen"){
#        columns <- c("numeric","numeric","numeric")
#    } else {
#        columns <- c("numeric","numeric","numeric","numeric")
#        #columns <- c("numeric","numeric","numeric")
#    }
#    classes <- c(classes,columns)
#}
#
#classes

#if (grepl("23",tophits[1])){
#	columns <- c("numeric","numeric","numeric","numeric")
#	classes <- c(classes,columns)
#}

#classes

#import.list <- llply(tophits, read.delim, sep="", colClasses=classes)

import.list <- llply(tophits, read.delim, sep="")
tophits_all <- rbindlist(import.list)

tophits_all <- as.data.frame(tophits_all)
head(tophits_all)

cat("#####                 Checking for duplicated variants				 #####\n")
table(duplicated(tophits_all$rs_id_all))

cat("#####                	Removing Duplicated Variant         	     #####\n")

tophits_all <- distinct(tophits_all,rs_id_all, .keep_all = TRUE)
table(duplicated(tophits_all$rs_id_all))

tophits_all_final <- tophits_all[,c("chr","position","rs_id_all","alleleA","alleleB","all_maf","refpanel","best_model")]

write.table(tophits_all_final,tophits_final,col.names=T,row.names=F,quote=F,sep="\t")

cat("#####            		  			DONE   				             #####\n")
