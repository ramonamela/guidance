library(plyr)
library(IRanges)

args <- commandArgs(TRUE)

filtered <- args[1]
filtered_males <- args[2]
filtered_females <- args[3]
output_condensed <- args[4]
output_tophits <- args[5]
output_ranges <- args[6]
pval_threshold <- args[7]

inheritance_models <- c("add","dom","rec","gen","het")

pval <- as.numeric(pval_threshold)

	if (filtered==filtered_males){

		filtered <- read.table(filtered,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		chr <- 23

        if (nrow(filtered)!=0){
			chr <- filtered$chr[1]
		}

		if (chr!=23){
			condensed <- filtered[,c("chr",
						"position",
						"rs_id_all",
						"info_all",
						"alleleA",
						"alleleB",
						"all_maf",
						"frequentist_add_se_1",
						"frequentist_add_beta_1",
						"frequentist_add_pvalue",
                        "frequentist_dom_se_1",
						"frequentist_dom_beta_1",
						"frequentist_dom_pvalue",
						"frequentist_rec_se_1",
						"frequentist_rec_beta_1",
						"frequentist_rec_pvalue",
						"frequentist_gen_se_1",
						"frequentist_gen_beta_1",
						"frequentist_gen_pvalue",
						"frequentist_het_se_1",
						"frequentist_het_beta_1",
						"frequentist_het_pvalue"),]

			names(condensed)[c(8,9,11,12,14,15,17,18,20,21)] <- c("frequentist_add_se","frequentist_add_beta",
																	"frequentist_dom_se","frequentist_dom_beta",
																	"frequentist_rec_se","frequentist_rec_beta",
																	"frequentist_gen_se","frequentist_gen_beta",
																	"frequentist_het_se","frequentist_het_beta")

		} else {
	        filtered_males <- read.table(filtered_males,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
			filtered_females <- read.table(filtered_females,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		
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

			condensed <- rbind.fill(condensed_males,condensed_females)

	}

	
	} else {

		filtered <- read.table(filtered,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		filtered_males <- read.table(filtered_males,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		filtered_females <- read.table(filtered_females,header=T,stringsAsFactors=F,na.strings=c("-","NA"))

	    condensed <- filtered[,c("chr",
    	 				           "position",
                       				"rs_id_all",
                        			"info_all",
                        			"alleleA",
                        			"alleleB",
                        			"all_maf",
                        			"frequentist_add_se_1",
                        			"frequentist_add_beta_1",
                        			"frequentist_add_pvalue",
                        			"frequentist_dom_se_1",
                        			"frequentist_dom_beta_1",
                        			"frequentist_dom_pvalue",
                        			"frequentist_rec_se_1",
                        			"frequentist_rec_beta_1",
                        			"frequentist_rec_pvalue",
                        			"frequentist_gen_se_1",
                        			"frequentist_gen_beta_1",
                        			"frequentist_gen_pvalue",
                        			"frequentist_het_se_1",
                        			"frequentist_het_beta_1",
                        			"frequentist_het_pvalue"),]

        names(condensed)[c(8,9,11,12,14,15,17,18,20,21)] <- c("frequentist_add_se","frequentist_add_beta",
                                                               "frequentist_dom_se","frequentist_dom_beta",
                                                               "frequentist_rec_se","frequentist_rec_beta",
                                                               "frequentist_gen_se","frequentist_gen_beta",
                                                               "frequentist_het_se","frequentist_het_beta")

		condensed$info_all <- as.numeric(condensed$info_all)

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

		condensed <- rbind.fill(condensed,condensed_males,condensed_females)
		condensed <- condensed[order(condensed$chr),]
	}

tophits <- condensed[condensed$frequentist_add_pvalue<=pval & !is.na(condensed$frequentist_add_pvalue) |
					 condensed$frequentist_dom_pvalue<=pval & !is.na(condensed$frequentist_dom_pvalue) | 
					 condensed$frequentist_rec_pvalue<=pval & !is.na(condensed$frequentist_rec_pvalue) |
					 condensed$frequentist_gen_pvalue<=pval & !is.na(condensed$frequentist_gen_pvalue) |
					 condensed$frequentist_het_pvalue<=pval & !is.na(condensed$frequentist_het_pvalue),]

write.table(condensed,output_condensed,col.names=T,row.names=F,quote=F,sep="\t")
write.table(tophits,output_tophits,col.names=T,row.names=F,quote=F,sep="\t")

tophits$start <- tophits$position - 250000
chr <- as.vector(unique(tophits$chr))

ranges <- NULL
    for(i in chr){
        dat <- tophits[tophits$chr == i,]
        dat <- IRanges(start=dat$start,width=500000)
        dat <- reduce(dat)
        dat <- as.data.frame(dat)
        dat$chr <- i
        dat$inheritance_models <- NA
        for (n in 1:nrow(dat)){
            dat$num_variants[n] <- nrow(tophits[(tophits$position>=dat$start[n] & tophits$position<=dat$end[n] & tophits$chr<=dat$chr),])
            for (i in inheritance_models){
                if (any(tophits$position>=dat$start[n] & tophits$position<=dat$end[n] &
                        tophits$chr<=dat$chr & tophits[,c(paste("frequentist_",i,"_pvalue",sep=""))]<=pval &
                        !(is.na(tophits[,c(paste("frequentist_",i,"_pvalue",sep=""))])))){
                    if (is.na(dat$inheritance_models[n])){
                        dat$inheritance_models[n] <- paste(i)
                    } else {
                        dat$inheritance_models[n] <- paste(dat$inheritance_models[n],i,sep=",")
                    }
                }
            }
        }
        ranges <- rbind(ranges,dat)
    }

if (nrow(tophits)!=0){
    ranges$Range <- paste("chr",ranges$chr,":",ranges$start,"-",ranges$end,sep="")
} else {
    ranges <- data.frame(matrix(ncol = 4, nrow = 0))
    x <- c("chr","start","end","width")
    colnames(ranges) <- x
}

write.table(ranges,output_ranges,col.names=T,row.names=F,quote=F,sep="\t")
