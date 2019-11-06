library(plyr)
library(IRanges)
library(data.table)

args <- commandArgs(TRUE)

filtered <- args[1]
filtered_males <- args[2]
filtered_females <- args[3]
filtered_x <- args[4]
output_condensed <- args[5]
output_tophits <- args[6]
output_ranges <- args[7]
pval_threshold <- args[8]
inheritance_models <- unlist(strsplit(args[9],","))

condensed_models_all <- NULL
pval <- as.numeric(pval_threshold)

	if(filtered==filtered_males & filtered_males==filtered_females & filtered_females==filtered_x){
        filtered_auto <- read.table(filtered,header=T,stringsAsFactors=F,na.strings=c("-","NA"))

        condensed <- filtered_auto[,c("chr",
			                    "position",
            		            "rs_id_all",
                    		    "info_all",
                       			 "alleleA",
                        		"alleleB",
                        		"all_maf",
                        		"refpanel")]
        for (m in inheritance_models){
            if (m !="gen"){
                condensed_models <- filtered_auto[,c(
                                    paste("frequentist_",m,"_se_1",sep=""),
                                    paste("frequentist_",m,"_beta_1",sep=""),
                                    paste("frequentist_",m,"_pvalue",sep=""))]
                condensed <- cbind(condensed,condensed_models)
            } else {
                condensed_models_gen <- filtered_auto[,c("frequentist_gen_se_1",
                                                "frequentist_gen_beta_1",
                                                "frequentist_gen_beta_2",
                                                "frequentist_gen_pvalue")]
                condensed <- cbind(condensed,condensed_models_gen)
            }
        }
		condensed$info_all <- as.numeric(condensed$info_all)
	}


#	if (filtered==filtered_males & filtered_males!=filtered_females & filtered_females!=filtered_x){
#
#		filtered <- read.table(filtered,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
#		chr <- 23
#
#        if (nrow(filtered)!=0){
#			chr <- filtered$chr[1]
#		}
#
#		if (chr!=23){
#			condensed <- filtered[,c("chr",
#						"position",
#						"rs_id_all",
#						"info_all",
#						"alleleA",
#						"alleleB",
#						"all_maf",
# 						"refpanel")]
#			for (m in inheritance_models){
#				if (m !="gen"){
#					condensed_models <- filtered[,c(
#										paste("frequentist_",m,"_se_1",sep=""),
#										paste("frequentist_",m,"_beta_1",sep=""),
#										paste("frequentist_",m,"_pvalue",sep=""))]
#					condensed <- cbind(condensed,condensed_models)
#				} else {
#					condensed_models_gen <- filtered[,c("frequentist_gen_se_1",
#													"frequentist_gen_beta_1",
#                        							"frequentist_gen_beta_2",
#													"frequentist_gen_pvalue")]
#					condensed <- cbind(condensed,condensed_models_gen)
#				}
#			}
#		} 

		if (filtered==filtered_males & filtered_males!=filtered_females & filtered_females!=filtered_x) {
	        filtered_x_males <- read.table(filtered_males,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
			head(filtered_x_males)
			filtered_x_females <- read.table(filtered_females,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
			head(filtered_x_females)
			filtered_x_all <- read.table(filtered_x,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		
	        condensed_males <- filtered_x_males[,c("chr",
        	                                "position",
                	                        "rs_id_all",
                        	                "info_all",
                                	        "alleleA",
											"alleleB",
        	                                "all_maf",
                	                        "frequentist_add_se_1.genotype.sex.1",
                        	                "frequentist_add_beta_1.genotype.sex.1",
                                	        "frequentist_add_pvalue",
											"refpanel")]

			if (nrow(filtered_x_males)!=0){
		        	condensed_males$chr <- as.character("23_males")
			}

	        names(condensed_males)[c(8,9)] <- c("frequentist_add_se_sex.1","frequentist_add_beta_sex.1")
			head(filtered_x_females)
	        condensed_females <- filtered_x_females[,c("chr",
	                                                "position",
        	                                        "rs_id_all",
                	                                "info_all",
                        	                        "alleleA",
                                	                "alleleB",
                                        	        "all_maf",
                                 	                "frequentist_add_se_1.genotype.sex.2",
                                        	        "frequentist_add_beta_1.genotype.sex.2",
                                                	"frequentist_add_pvalue",
													"refpanel")]

			if (nrow(filtered_x_females)!=0){
		        condensed_females$chr <- as.character("23_females")
			}

            names(condensed_females)[c(8,9)] <- c("frequentist_add_se_sex.2","frequentist_add_beta_sex.2")

            condensed_x <- filtered_x_all[,c("chr",
                                             "position",
                                             "rs_id_all",
                                             "info_all",
                                             "alleleA",
                                             "alleleB",
                                             "all_maf",
                                             "frequentist_add_se_1.genotype.sex.1",
                                             "frequentist_add_se_2.genotype.sex.2",
                                             "frequentist_add_beta_1.genotype.sex.1",
                                             "frequentist_add_beta_2.genotype.sex.2",
                                             "frequentist_add_pvalue",
                                             "refpanel")]

            if (nrow(filtered_x_all)!=0){
                    condensed_x$chr <- as.character("23")
            }

            names(condensed_x)[c(8,9,10,11)] <- c("frequentist_add_se_sex.1","frequentist_add_se_sex.2","frequentist_add_beta_sex.1","frequentist_add_beta_sex.2")

			condensed <- rbind.fill(condensed_males,condensed_females,condensed_x)

	}

	
	if (filtered!=filtered_males & filtered_males!=filtered_females & filtered_females!=filtered_x){

		filtered_auto <- read.table(filtered,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		filtered_x_males <- read.table(filtered_males,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		head(filtered_x_males)
		filtered_x_females <- read.table(filtered_females,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		filtered_x_all <- read.table(filtered_x,header=T,stringsAsFactors=F,na.strings=c("-","NA"))
		head(filtered_x_females)

	    condensed <- filtered_auto[,c("chr",
    	 				           "position",
                       				"rs_id_all",
                        			"info_all",
                        			"alleleA",
                        			"alleleB",
                        			"all_maf",
									"refpanel"),]
        	for (m in inheritance_models){
                if (m !="gen"){
                    condensed_models <- filtered_auto[,c(
                                        paste("frequentist_",m,"_se_1",sep=""),
                                        paste("frequentist_",m,"_beta_1",sep=""),
                                        paste("frequentist_",m,"_pvalue",sep=""))]
                    condensed <- cbind(condensed,condensed_models)
                } else {
                    condensed_models_gen <- filtered_auto[,c("frequentist_gen_se_1",
                                                    "frequentist_gen_beta_1",
                                                    "frequentist_gen_beta_2",
                                                    "frequentist_gen_pvalue")]
                    condensed <- cbind(condensed,condensed_models_gen)
                }
            }
		condensed$info_all <- as.numeric(condensed$info_all)

		condensed_males <- filtered_x_males[,c("chr",
											"position",
											"rs_id_all",
			                         		"info_all",
           				        	   		"alleleA",
                      		 	    		"alleleB",
                      						"all_maf",
                        			   		"frequentist_add_se_1.genotype.sex.1",
                			          		"frequentist_add_beta_1.genotype.sex.1",
           					   				"frequentist_add_pvalue",
											"refpanel")]

		if (nrow(filtered_x_males)!=0){
			condensed_males$chr <- as.character("23_males")
		}
		names(condensed_males)[c(8,9)] <- c("frequentist_add_se_sex.1","frequentist_add_beta_sex.1")

		condensed_females <- filtered_x_females[,c("chr",
			        	                		"position",
				   				    	        "rs_id_all",
           			    				        "info_all",
												"alleleA",
                      			      			"alleleB",
               				           			"all_maf",
                        				 		"frequentist_add_se_1.genotype.sex.2",
                      			     			"frequentist_add_beta_1.genotype.sex.2",
               									"frequentist_add_pvalue",
												"refpanel")]

		if (nrow(filtered_x_females)!=0){
			condensed_females$chr <- as.character("23_females")
		}
		names(condensed_females)[c(8,9)] <- c("frequentist_add_se_sex.2","frequentist_add_beta_sex.2")


        condensed_x <- filtered_x_all[,c("chr",
                                     "position",
                                     "rs_id_all",
                                     "info_all",
                                     "alleleA",
                                     "alleleB",
                                     "all_maf",
                                     "frequentist_add_se_1.genotype.sex.1",
                                     "frequentist_add_se_2.genotype.sex.2",
                                     "frequentist_add_beta_1.genotype.sex.1",
                                     "frequentist_add_beta_2.genotype.sex.2",
                                     "frequentist_add_pvalue",
                                     "refpanel")]

        if (nrow(filtered_x_all)!=0){
                condensed_x$chr <- as.character("23")
        }

        names(condensed_x)[c(8,9,10,11)] <- c("frequentist_add_se_sex.1","frequentist_add_se_sex.2","frequentist_add_beta_sex.1","frequentist_add_beta_sex.2")

		condensed <- rbind.fill(condensed,condensed_males,condensed_females,condensed_x)
		condensed <- condensed[order(condensed$chr),]
	}

#filtering info because the old GERA analysis was made for a info_score=0.5
condensed <- condensed[condensed$info_all>=0.7,]

#for (n in 1:nrow(condensed)){
#    condensed$best_model <- lapply(strsplit(strsplit(names(condensed)[grep("_pvalue",names(condensed))[which.min(condensed[,grep("_pvalue",names(condensed))])]],"_pvalue")[[1]][1],"_")[[1]][2])
#}

names(condensed)

write.table(condensed,output_condensed,col.names=T,row.names=F,quote=F,sep="\t")

condensed$best_model <- NA
tophits <- NULL
tophits_models <- NULL
tophits_models_gen <- NULL
for (m in inheritance_models){
			tophits_models <- condensed[which(as.numeric(condensed[,c(paste("frequentist_",m,"_pvalue",sep=""))])<=pval & 
							!is.na(condensed[,c(paste("frequentist_",m,"_pvalue",sep=""))])),]
			tophits <- rbind.fill(tophits,tophits_models)
}

head(tophits)
if (nrow(tophits)!=0){
	for (n in 1:nrow(tophits)){
		tophits$best_model[n] <- strsplit(strsplit(names(tophits)[grep("_pvalue",names(tophits))[which.min(tophits[n,grep("_pvalue",names(tophits))])]],"_pvalue")[[1]][1],"_")[[1]][2]
	}
}


tophits <- tophits[!duplicated(tophits),]
table(tophits$chr)
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
            dat$num_variants[n] <- dim(tophits[(tophits$position>=dat$start[n] & tophits$position<=dat$end[n] & tophits$chr==i),])[1]
			if (i !="23" & i!="23_males" & i!="23_females"){
            	for (inh in inheritance_models){
                	if (any(tophits$position>=dat$start[n] & tophits$position<=dat$end[n] &
                        	tophits$chr==i & tophits[,c(paste("frequentist_",inh,"_pvalue",sep=""))]<=pval &
                        	!(is.na(tophits[,c(paste("frequentist_",inh,"_pvalue",sep=""))])))){
                    	if (is.na(dat$inheritance_models[n])){
                        	dat$inheritance_models[n] <- paste(inh)
                    	} else {
                        	dat$inheritance_models[n] <- paste(dat$inheritance_models[n],inh,sep=",")
                    	}
                	}
            	}
			} else {
              	dat$inheritance_models[n] <- paste("add")
			}
		}
        ranges <- rbind(ranges,dat)
    }

if (nrow(tophits)!=0){
    ranges$range <- paste("chr",ranges$chr,":",ranges$start,"-",ranges$end,sep="")
} else {
    ranges <- data.frame(matrix(ncol = 4, nrow = 0))
    x <- c("chr","start","end","width")
    colnames(ranges) <- x
}

head(ranges)
dim(ranges)

tophits_all <- tophits
tophits_ranges <- ranges
topvariant <- NULL
topvariants <- NULL
tophits_topvariants_auto <- NULL
topvariant_final <- NULL

if (nrow(tophits_ranges)!=0){
	for (i in 1:length(tophits_ranges$range)){
		if (tophits_ranges$chr[i]!="23" & tophits_ranges$chr[i]!="23_males" & tophits_ranges$chr[i]!="23_females"){
    		tophits_set <- tophits_all[(tophits_all$chr==tophits_ranges$chr[i] & 
						   tophits_all$position>=tophits_ranges$start[i] & 
						   tophits_all$position<=tophits_ranges$end[i]),]
			print(head(tophits_set))
			print(table(duplicated(tophits_set$rs_id_all)))
    		for (m in unlist(strsplit(as.character(tophits_ranges$inheritance_models[i]),","))){
				if (m!="gen"){
        			if (any(!is.na((tophits_set[,c(paste("frequentist_",m,"_pvalue",sep=""))])))){
            			topvariant <- as.data.frame(tophits_set[which.min(tophits_set[,c(paste("frequentist_",m,"_pvalue",sep=""))]),
															c("position","rs_id_all","alleleA","alleleB","all_maf",
															paste("frequentist_",m,"_pvalue",sep=""),
                            	                            paste("frequentist_",m,"_se_1",sep=""),
															paste("frequentist_",m,"_beta_1",sep=""),
															"info_all","refpanel")])
            			names(topvariant)[6] <- "pvalue"
            			names(topvariant)[7] <- "se"
						names(topvariant)[8] <- "beta_1"
            			topvariant$model <- m
						topvariant$range <- i
						if ("add" %in% inheritance_models){
	            			topvariant$add_pvalue <- tophits_set[topvariant$rs_id_all==tophits_set$rs_id_all,c("frequentist_add_pvalue")]
						}
					}
        		} else {
						if (any(!is.na((tophits_set[,c(paste("frequentist_",m,"_pvalue",sep=""))])))){
                        	topvariant <- as.data.frame(tophits_set[which.min(tophits_set[,c(paste("frequentist_",m,"_pvalue",sep=""))]),
                                                            c("position","rs_id_all","alleleA","alleleB","all_maf",
                                                            paste("frequentist_",m,"_pvalue",sep=""),
                                                            paste("frequentist_",m,"_se_1",sep=""),
                                                            paste("frequentist_",m,"_beta_1",sep=""),
                                                            paste("frequentist_",m,"_beta_2",sep=""),
                                                            "info_all","refpanel")])
                    		names(topvariant)[6] <- "pvalue"
                    		names(topvariant)[7] <- "se"
                    		names(topvariant)[8] <- "beta_1"
                   			names(topvariant)[9] <- "beta_2"
                    		topvariant$model <- m
							topvariant$range <- i
							if ("add" %in% inheritance_models){
	                    		topvariant$add_pvalue <- tophits_set[topvariant$rs_id_all==tophits_set$rs_id_all,c("frequentist_add_pvalue")]
							}
						}			
				}	
				topvariants <- rbind.fill(topvariants,topvariant)
			}
			topvariants$TOPHIT <- "NO"
 			for (r in unique(topvariants$range)){
					topvariants$TOPHIT[which(topvariants$range==r)[which.min(topvariants[topvariants$range==r,c("pvalue")])]] <- "YES"
				
			}
			tophits_topvariants_auto <- topvariants[topvariants$TOPHIT=="YES",]
		}
	}		
}

topvariant <- NULL
topvariants <- NULL
tophits_topvariants_x <- NULL
topvariant_final <- NULL

if (nrow(tophits_ranges)!=0){
    for (i in 1:length(tophits_ranges$range)){
        if (tophits_ranges$chr[i]=="23" | tophits_ranges$chr[i]=="23_males" | tophits_ranges$chr[i]=="23_females"){
            tophits_set <- tophits_all[(as.character(tophits_all$chr)==tophits_ranges$chr[i] &
                           tophits_all$position>=tophits_ranges$start[i] &
                           tophits_all$position<=tophits_ranges$end[i]),]
            for (m in unlist(strsplit(as.character(tophits_ranges$inheritance_models[i]),","))){
                if (any(!is.na((tophits_set[,c(paste("frequentist_",m,"_pvalue",sep=""))])))){
                    topvariant <- as.data.frame(tophits_set[which.min(tophits_set[,c(paste("frequentist_",m,"_pvalue",sep=""))]),
                                                            c("position","rs_id_all","alleleA","alleleB","all_maf",
                                                            paste("frequentist_",m,"_pvalue",sep=""),
                                                            paste("frequentist_",m,"_se_sex.1",sep=""),
                                                            paste("frequentist_",m,"_beta_sex.1",sep=""),
                                                            paste("frequentist_",m,"_se_sex.2",sep=""),
                                                            paste("frequentist_",m,"_beta_sex.2",sep=""),
                                                            "info_all","refpanel")])
                names(topvariant)[6] <- "pvalue"
                names(topvariant)[7] <- "se_males"
		names(topvariant)[8] <- "beta_males"
                names(topvariant)[9] <- "se_females"
                names(topvariant)[10] <- "beta_females"
		topvariant$range <- i
                topvariant$model <- m
               	topvariant$add_pvalue <- tophits_set[topvariant$rs_id_all==tophits_set$rs_id_all,c("frequentist_add_pvalue")]
                }
           		topvariants <- rbind.fill(topvariants,topvariant)
			}
			topvariants$TOPHIT <- "NO"
            for (r in unique(topvariants$range)){
                    topvariants$TOPHIT[which(topvariants$range==r)[which.min(topvariants[topvariants$range==r,c("pvalue")])]] <- "YES"
                
            }
            tophits_topvariants_x <- topvariants[topvariants$TOPHIT=="YES",]
		}	   
	}
}

if (is.null(tophits_topvariants_auto)){
	tophits_topvariants_auto  <- data.frame(matrix(ncol = 14, nrow = 0))
	names(tophits_topvariants_auto) <- c("position","rs_id_all","alleleA","alleleB","all_maf",
                                    "pvalue","se_1","se_2","beta_1","beta_2","model","add_pvalue","info_all","refpanel")
}

if (is.null(tophits_topvariants_x)){
	tophits_topvariants_x  <- data.frame(matrix(ncol = 14, nrow = 0))
	names(tophits_topvariants_x) <- c("position","rs_id_all","alleleA","alleleB","all_maf",
					"pvalue","se_males","beta_males","se_females","beta_females","model","add_pvalue","info_all","refpanel")
}

if (nrow(tophits_ranges)==0){
	tophits_ranges <- NULL
	tophits_ranges  <- data.frame(matrix(ncol = 7, nrow = 0))
	names(tophits_ranges) <- c("start","end","width","chr","inheritance_models","num_variants","range")
}
	
tophits_topvariants <- rbind.fill(tophits_topvariants_auto,tophits_topvariants_x)
tophits_final <- cbind(tophits_ranges,tophits_topvariants)

if (nrow(tophits_ranges)!=0){
	tophits_final <- tophits_final[,-c(which(names(tophits_final)=="TOPHIT"),which(names(tophits_final)=="range"))]
}

#if (nrow(tophits_ranges)==0){
#	tophits_final <- data.frame(matrix(ncol = 24, nrow = 0))
#			tophits_final_names <- c("start","end","width","chr","inheritance_models",
#            	        	           "num_variants","Range","position","rs_id_all","alleleA","alleleB",
#                		               "all_maf","pvalue","se","beta_1","beta_2","se_males","beta_males","se_females","beta_females","model","add_pvalue","info_all","refpanel")
#			colnames(tophits_final) <- tophits_final_names
#		} else {
#    		tophits_final <- tophits_final[,-c(which(names(tophits_final)=="TOPHIT"),which(names(tophits_final)=="range"))]
#		}
#
#} else {
#	tophits_topvariants_auto <- data.frame(matrix(ncol = 10, nrow = 0))
#	tophits_topvariants_auto <- c("chr","position","rs_id_all","alleleA","alleleB","all_maf","pvalue","se","beta_1","beta_2")
#	tophits_topvariants <- rbind.fill(tophits_topvariants_auto,tophits_topvariants_x)
#	tophits_final <- cbind(tophits_ranges,tophits_topvariants)
#        if (nrow(tophits_ranges)==0){
#            tophits_final <- data.frame(matrix(ncol = 24, nrow = 0))
#            tophits_final_names <- c("start","end","width","chr","inheritance_models",
#                                       "num_variants","Range","position","rs_id_all","alleleA","alleleB",
#                                       "all_maf","pvalue","se","beta_1","beta_2","se_males","beta_males","se_females","beta_females","model","add_pvalue","info_all","refpanel")
#            colnames(tophits_final) <- tophits_final_names
#        } else {
#            tophits_final <- tophits_final[,-c(which(names(tophits_final)=="TOPHIT"),which(names(tophits_final)=="range"))]
#        }
#}

write.table(tophits_final,output_ranges,col.names=T,row.names=F,quote=F,sep="\t")
