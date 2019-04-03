args <- commandArgs(TRUE)

filtered <- args[1]
filtered_males <- args[2]
filtered_females <- args[3]
output_condensed <- args[4]
output_tophits <- args[5]
pval_threshold <- args[6]

pval <- as.numeric(pval_threshold)

if (filtered==filtered_males){
    
    filtered <- read.table(filtered,header=T,stringsAsFactors=F,na.strings="-")
    chr <- 23
    
    if (nrow(filtered)!=0){
        chr <- filtered$chr[1]
    }
    
    if (chr!=23){
        condensedA <- filtered[,c("chr",
                                  "position",
                                  "rs_id_all",
                                  "info_all",
                                  "alleleA",
                                  "alleleB",
                                  "all_maf",
                                  "frequentist_add_se_1",
                                  "frequentist_add_beta_1",
                                  "frequentist_add_pvalue")]
        
        names(condensedA)[c(8,9)] <- c("frequentist_add_se","frequentist_add_beta")
        condensed <- condensedA
        
    } else {
        filtered_males <- read.table(filtered_males,header=T,stringsAsFactors=F,na.strings="-")
        filtered_females <- read.table(filtered_females,header=T,stringsAsFactors=F,na.strings="-")
        
        condensed_males <- filtered_males[,c("chr",
                                             "position",
                                             "rs_id_all",
                                             "info_all",
                                             "alleleA",
                                             "alleleB",
                                             "all_maf",
                                             "frequentist_add_se_1.genotype.sex.1",
                                             "frequentist_add_beta_1.genotype.sex.1",
                                             "frequentist_add_pvalue")]
        if (nrow(filtered_males)!=0){
            condensed_males$chr <- as.character("23_males")
        }
        
        names(condensed_males)[c(8,9)] <- c("frequentist_add_se","frequentist_add_beta")
        condensed_females <- filtered_females[,c("chr",
                                                 "position",
                                                 "rs_id_all",
                                                 "info_all",
                                                 "alleleA",
                                                 "alleleB",
                                                 "all_maf",
                                                 "frequentist_add_se_1.genotype.sex.2",
                                                 "frequentist_add_beta_1.genotype.sex.2",
                                                 "frequentist_add_pvalue")]
        if (nrow(filtered_females)!=0){
            condensed_females$chr <- as.character("23_females")
        }
        
        names(condensed_females)[c(8,9)] <- c("frequentist_add_se","frequentist_add_beta")
        condensed <- rbind(condensed_males,condensed_females)
    }
} else {
    
    filtered <- read.table(filtered,header=T,stringsAsFactors=F,na.strings="-")
    filtered_males <- read.table(filtered_males,header=T,stringsAsFactors=F,na.strings="-")
    filtered_females <- read.table(filtered_females,header=T,stringsAsFactors=F,na.strings="-")
    
    condensed <- filtered[,c("chr",
                             "position",
                             "rs_id_all",
                             "info_all",
                             "alleleA",
                             "alleleB",
                             "all_maf",
                             "frequentist_add_se_1",
                             "frequentist_add_beta_1",
                             "frequentist_add_pvalue")]
    
    condensed$info_all <- as.numeric(condensed$info_all)
    names(condensed)[c(8,9)] <- c("frequentist_add_se","frequentist_add_beta")
    
    condensed_males <- filtered_males[,c("chr",
                                         "position",
                                         "rs_id_all",
                                         "info_all",
                                         "alleleA",
                                         "alleleB",
                                         "all_maf",
                                         "frequentist_add_se_1.genotype.sex.1",
                                         "frequentist_add_beta_1.genotype.sex.1",
                                         "frequentist_add_pvalue")]
    
    if (nrow(filtered_males)!=0){
        condensed_males$chr <- as.character("23_males")
    }
    names(condensed_males)[c(8,9)] <- c("frequentist_add_se","frequentist_add_beta")
    
    condensed_females <- filtered_females[,c("chr",
                                             "position",
                                             "rs_id_all",
                                             "info_all",
                                             "alleleA",
                                             "alleleB",
                                             "all_maf",
                                             "frequentist_add_se_1.genotype.sex.2",
                                             "frequentist_add_beta_1.genotype.sex.2",
                                             "frequentist_add_pvalue")]
    if (nrow(filtered_females)!=0){
        condensed_females$chr <- as.character("23_females")
    }
    names(condensed_females)[c(8,9)] <- c("frequentist_add_se","frequentist_add_beta")
    
    condensed <- rbind(condensed,condensed_males,condensed_females)
    condensed <- condensed[order(condensed$chr),]
}

tophits <- condensed[condensed$frequentist_add_pvalue<=pval,]

write.table(condensed,output_condensed,col.names=T,row.names=F,quote=F,sep="\t")
write.table(tophits,output_tophits,col.names=T,row.names=F,quote=F,sep="\t")
