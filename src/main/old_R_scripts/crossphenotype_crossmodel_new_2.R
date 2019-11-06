library(data.table)
library(plyr)
library(dplyr)
library(reshape)
library(IRanges)

args <- commandArgs(TRUE)

tophits <- unlist(strsplit(args[1], ","))
output_summary <- args[2]
pval <- as.numeric(args[3])
models <- unlist(strsplit(args[4], ","))

import.list <- llply(tophits, read.delim, sep = "")
#for (i in 1:length(import.list)){
#	disease <- strsplit(names(import.list[[i]])[8],"frequentist_add_se_")[[1]][2]
#	names(import.list[[i]])[11:22] <- paste(names(import.list[[i]])[11:22],"_",disease,sep="")
#}

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

#ranges

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

crosspheno_all <- NULL
for (l in 1:length(cross_all)) {
  cross_pheno_all <- cross_all[[l]]
  if (nrow(cross_pheno_all) != 0) {
    for (i in 1:nrow(cross_pheno_all)) {
      r <-
        which(
          ranges$start <= cross_pheno_all$position[i] &
            cross_pheno_all$position[i] <= ranges$end &
            ranges$chr == cross_pheno_all$chr[i]
        )
      cross_pheno_all$RANGE[i] <-
        paste("chr",
              ranges$chr[r],
              ":",
              ranges$start[r],
              "-",
              ranges$end[r],
              sep = "")
    }
    diseases <-
      names(cross_pheno_all)[which(grepl("pvalue", names(cross_pheno_all)))]
    #    pval_threshold <- 0.05/length(unique(cross_pheno_all$RANGE))/length(diseases)
    pval_threshold <- 0.05
    
    crosspheno_variants <- NULL
    line <- 1
    for (i in 1:nrow(cross_pheno_all)) {
      for (n in diseases) {
        if (cross_pheno_all[i, n] <= pval & !is.na(cross_pheno_all[i, n])) {
          for (e in diseases) {
            if (e != n) {
              if (cross_pheno_all[i, e] <= pval_threshold &
                  !is.na(cross_pheno_all[i, e])) {
                crosspheno_variants$range[line] <- cross_pheno_all$RANGE[i]
                crosspheno_variants$rsid[line] <-
                  as.character(cross_pheno_all$rs_id_all[i])
                crosspheno_variants$disease_A[line] <-
                  as.character(lapply(strsplit(n, "pvalue_"), "[", 2))
                crosspheno_variants$pvalue_A[line] <-
                  cross_pheno_all[i, n]
                crosspheno_variants$beta_A[line] <-
                  cross_pheno_all[i, paste(as.character(lapply(
                    strsplit(n, "pvalue_"), "[", 1
                  )),
                  "beta_",
                  as.character(lapply(
                    strsplit(n, "pvalue_"), "[", 2
                  )),
                  sep = "")]
                crosspheno_variants$se_A[line] <-
                  cross_pheno_all[i, paste(as.character(lapply(
                    strsplit(n, "pvalue_"), "[", 1
                  )),
                  "se_",
                  as.character(lapply(
                    strsplit(n, "pvalue_"), "[", 2
                  )),
                  sep = "")]
                crosspheno_variants$disease_B[line] <-
                  as.character(lapply(strsplit(e, "pvalue_"), "[", 2))
                crosspheno_variants$pvalue_B[line] <-
                  cross_pheno_all[i, e]
                crosspheno_variants$beta_B[line] <-
                  cross_pheno_all[i, paste(as.character(lapply(
                    strsplit(e, "pvalue_"), "[", 1
                  )),
                  "beta_",
                  as.character(lapply(
                    strsplit(e, "pvalue_"), "[", 2
                  )),
                  sep = "")]
                crosspheno_variants$se_B[line] <-
                  cross_pheno_all[i, paste(as.character(lapply(
                    strsplit(e, "pvalue_"), "[", 1
                  )),
                  "se_",
                  as.character(lapply(
                    strsplit(e, "pvalue_"), "[", 2
                  )),
                  sep = "")]
                line <- line + 1
              }
            }
          }
        }
      }
    }
    
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
      crosspheno_all <- rbind.fill(crosspheno_all, crosspheno_min_pval)
      
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
