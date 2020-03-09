init_time <- Sys.time()

library(data.table)
library(plyr)
library(dplyr)
library(reshape)
library(IRanges)
library(parallel)

args <- commandArgs(TRUE)

tophits <- unlist(strsplit(args[1], ","))
output_summary <- args[2]
pval <- as.numeric(args[3])
models <- unlist(strsplit(args[4], ","))

init_t0 <- Sys.time()

import.list <- llply(tophits, read.delim, sep = "")

end_t0 <- Sys.time()
print(paste("Read input files:" , end_t0 - init_t0))

#for (i in 1:length(import.list)){
#	disease <- strsplit(names(import.list[[i]])[8],"frequentist_add_se_")[[1]][2]
#	names(import.list[[i]])[11:22] <- paste(names(import.list[[i]])[11:22],"_",disease,sep="")
#}

init_t1 <- Sys.time()
cross_pheno_all <- Reduce(function(dtf1, dtf2)
  merge(
    dtf1,
    dtf2,
    by = c(
      "chr",
      "position",
      "rs_id_all",
      "alleleA",
      "alleleB",
      "info_all",
      "all_maf",
      "refpanel",
      "best_model"
    ),
    all = TRUE
  ), import.list)

#cross_pheno_all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2,
#                                                     by = c("chr","position",
#                                                            "rs_id_all","alleleA","alleleB",
#                                                            "info_all","all_maf","refpanel"),
#                                                     all = TRUE),import.list)

cross_pheno_all$start <- cross_pheno_all$position - 250000
chr <- as.vector(unique(cross_pheno_all$chr))

ranges <- NULL
for (i in chr) {
  dat <- cross_pheno_all[cross_pheno_all$chr == i, ]
  dat <- IRanges(start = dat$start, width = 500000)
  dat <- reduce(dat)
  dat <- as.data.frame(dat)
  dat$chr <- i
  ranges <- rbind(ranges, dat)
}

if (nrow(cross_pheno_all) != 0) {
  ranges$Range <-
    paste("chr", ranges$chr, ":", ranges$start, "-", ranges$end, sep = "")
} else {
  ranges <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("chr", "start", "end", "width")
  colnames(ranges) <- x
}

# new comment #ranges

# new comment #cat("#############################################################################\n")
# new comment #cat("#####                     CROSS-PHENOTYPE SUMMARY                       #####\n")
# new comment #cat("#############################################################################\n")

cross_all <- list()

i = 1
for (m in models) {
  cross_all[[i]] <-
    cross_pheno_all[, c(1:8, grep(m, names(cross_pheno_all)))]
  i = i + 1
}

#head(cross_all,n=10)

end_t1 <- Sys.time()
print(paste("First zone of code:" , end_t1 - init_t1))

add_cross_pheno_column <- function(position, chr) {
  r <-
    which(ranges$start <= position &
            position <= ranges$end &
            ranges$chr == chr)
  return(paste("chr",
               ranges$chr[r],
               ":",
               ranges$start[r],
               "-",
               ranges$end[r],
               sep = ""))
}

crosspheno_all <- NULL
for (l in 1:length(cross_all)) {
  cross_pheno_all <- cross_all[[l]]
  if (nrow(cross_pheno_all) != 0) {
    init_t2 <- Sys.time()
    amountOfAvailableCores = detectCores()
    cross_pheno_all$RANGE <-
      mcmapply(
        add_cross_pheno_column,
        cross_pheno_all$position,
        cross_pheno_all$chr,
        mc.cores = amountOfAvailableCores
      )
    end_t2 <- Sys.time()
    print(paste("Add cross pheno column:" , end_t2 - init_t2))
    #init_t2 <- Sys.time()
    #for (i in 1:nrow(cross_pheno_all)) {
    #  r <-
    #    which(
    #      ranges$start <= cross_pheno_all$position[i] &
    #        cross_pheno_all$position[i] <= ranges$end &
    #        ranges$chr == cross_pheno_all$chr[i]
    #    )
    #  cross_pheno_all$RANGE2[i] <-
    #    paste("chr",
    #          ranges$chr[r],
    #          ":",
    #          ranges$start[r],
    #          "-",
    #          ranges$end[r],
    #          sep = "")
    #}
    #end_t2 <- Sys.time()
    #print(paste("Add cross pheno column:" , end_t2 - init_t2))
    #all.equal(cross_pheno_all$RANGE,cross_pheno_all$RANGE2)
    init_t3 <- Sys.time()
    diseases <-
      names(cross_pheno_all)[which(grepl("pvalue", names(cross_pheno_all)))]
    #    pval_threshold <- 0.05/length(unique(cross_pheno_all$RANGE))/length(diseases)
    pval_threshold <- 0.05
    end_t3 <- Sys.time()
    print(paste("Grep pvalue:", end_t3 - init_t3))
    
    init_t4 <- Sys.time()
    crosspheno_variants <- NULL
    #line <- 1
    #for (i in 1:nrow(cross_pheno_all)) {
    #  for (n in diseases) {
    #    if (cross_pheno_all[i, n] <= pval & !is.na(cross_pheno_all[i, n])) {
    #      for (e in diseases) {
    #        if (e != n) {
    #          if (cross_pheno_all[i, e] <= pval_threshold &
    #              !is.na(cross_pheno_all[i, e])) {
    #            crosspheno_variants$range[line] <- cross_pheno_all$RANGE[i]
    #            crosspheno_variants$rsid[line] <-
    #              as.character(cross_pheno_all$rs_id_all[i])
    #            crosspheno_variants$disease_A[line] <-
    #              as.character(lapply(strsplit(n, "pvalue_"), "[", 2))
    #            crosspheno_variants$pvalue_A[line] <-
    #              cross_pheno_all[i, n]
    #            crosspheno_variants$disease_B[line] <-
    #              as.character(lapply(strsplit(e, "pvalue_"), "[", 2))
    #            crosspheno_variants$pvalue_B[line] <-
    #              cross_pheno_all[i, e]
    #            line <- line + 1
    #          }
    #        }
    #      }
    #    }
    #  }
    #}
    for (n in diseases) {
      for (e in diseases[diseases != n]) {
        #n <- "frequentist_rec_pvalue_ALLERGIC_RHINITIS"
        #e <- "frequentist_rec_pvalue_ASTHMA"
        first_mask <-
          unlist(lapply(cross_pheno_all[, n], function(x)
            x <= pval & !is.na(x)))
        first_true_indices <- which(first_mask == TRUE)
        first_true_elements <- cross_pheno_all[first_true_indices, ]
        second_mask <-
          unlist(lapply(first_true_elements[, e], function(x)
            x <= pval_threshold & !is.na(x)))
        second_true_indices <- which(second_mask == TRUE)
        second_true_elements <-
          first_true_elements[second_true_indices, ]
        
        #additional_crosspheno_variants_3 <- data.frame(matrix(ncol=5,nrow=nrow(second_true_elements)))
        #colnames(additional_crosspheno_variants_3) <- c("rsid", "disease_A", "pvalue_A", "disease_B", "pvalue_B")
        range_col <- as.character(second_true_elements$RANGE)
        rsid_col <- as.character(second_true_elements$rs_id_all)
        disease_A_col <-
          rep(as.character(lapply(strsplit(
            n, "pvalue_"
          ), "[", 2)),
          length(rsid_col))
        pvalue_A_col <- second_true_elements[, n]
        disease_B_col <-
          rep(as.character(lapply(strsplit(
            e, "pvalue_"
          ), "[", 2)),
          length(rsid_col))
        pvalue_B_col <- second_true_elements[, e]
        additional_crosspheno_variants_3 <-
          data.frame(
            range = range_col,
            rsid = rsid_col,
            disease_A = disease_A_col,
            pvalue_A = pvalue_A_col,
            disease_B = disease_B_col,
            pvalue_B = pvalue_B_col,
            stringsAsFactors = FALSE
          )
        
        crosspheno_variants <-
          rbind(crosspheno_variants,
                additional_crosspheno_variants_3)
      }
    }
    crosspheno_variants <-
      crosspheno_variants[order(crosspheno_variants$rsid),]
    row.names(crosspheno_variants) <- 1:nrow(crosspheno_variants)
    crosspheno_variants <- as.data.frame(crosspheno_variants)
    end_t4 <- Sys.time()
    print(paste("Crosspheno variants:", end_t4 - init_t4))
    
    if (length(crosspheno_variants) != 0) {
      crosspheno <- as.data.frame(crosspheno_variants)
      crosspheno_2 <- crosspheno[!is.na(crosspheno$disease_B), ]
      crosspheno <- crosspheno_2[order(crosspheno_2$range), ]
      
      crosspheno$disease_A_vs_disease_B <- NA
      crosspheno$disease_A_vs_disease_B <-
        paste(crosspheno$disease_A, "-", crosspheno$disease_B, sep = "")
      
      associations <-
        as.character(unique(crosspheno$disease_A_vs_disease_B))
      range <- as.character(unique(crosspheno$range))
      
      
      crosspheno_min_pval <- NULL
      for (i in 1:length(range)) {
        for (n in 1:length(associations)) {
          all <-
            crosspheno[crosspheno$range == range[i] &
                         crosspheno$disease_A_vs_disease_B == associations[n], ]
          crosspheno_min_pval <-
            rbind(crosspheno_min_pval, all[which.min(all$pvalue_A), ])
        }
      }
      model <- strsplit(diseases[1], "_")[[1]][2]
      crosspheno_min_pval$model <- model
      crosspheno_all <-
        rbind.fill(crosspheno_all, crosspheno_min_pval)
      
    } else {
      crosspheno_all <- data.frame(matrix(ncol = 2, nrow = 0))
      x <-
        c(
          "range",
          "rsid",
          "disease_A",
          "pvalue_A",
          "disease_B",
          "pvalue_B",
          "disease_A_vs_disease_B",
          "model"
        )
      colnames(crosspheno_all) <- x
    }
    
  } else {
    crosspheno_all <- data.frame(matrix(ncol = 2, nrow = 0))
    x <-
      c(
        "range",
        "rsid",
        "disease_A",
        "pvalue_A",
        "disease_B",
        "pvalue_B",
        "disease_A_vs_disease_B",
        "model"
      )
    colnames(crosspheno_min_pval) <- x
  }

  write.table(
    crosspheno_all,
    output_summary,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
  )
}

# new comment #cat("#############################################################################\n")
# new comment #cat("#####                  CROSS-PHENOTYPE SUMMARY DONE                     #####\n")
# new comment #cat("#############################################################################\n")

print(paste("Total amount of time:", Sys.time() - init_time))
