library(data.table)

args <- commandArgs(TRUE)

tophits <- args[1]
all_pheno <- args[2]
output <- args[3]
pheno <- args[4]


tophits <- read.delim(tophits,header=T,sep="") 
data <- read.delim(all_pheno,header=T,sep="")

data_merge <- merge(tophits,data,by.x=c("chr","position","rs_id_all","alleleA","alleleB"),
                                 by.y=c("chr","position","rs_id_all","alleleA","alleleB"))


colnames(data_merge)[8:10] <- paste(colnames(data_merge)[8:10],pheno,sep="_")
write.table(data_merge,output,col.names=T,row.names=F,quote=F,sep="\t")
