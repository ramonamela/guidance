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
inheritance_models <- unlist(strsplit(args[9], ","))

#filtered <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_21_to_22.txt.gz"
#filtered_males <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_23_males.txt.gz"
#filtered_females <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_23_females.txt.gz"
#filtered_x <- "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_23.txt.gz"
#output_condensed <-
#  "ALLERGIC_RHINITIS_uk10k_condensed_chr_21_to_23_ref.txt.gz"
#output_tophits <-
#  "tophits_ALLERGIC_RHINITIS_GERA_300_uk10k_ref.txt.gz"
#output_ranges <-
#  "tophits_ALLERGIC_RHINITIS_GERA_300_uk10k_crossmodel_ranges_ref.txt.gz"
#pval_threshold <- "0.05"
#inheritance_models <- unlist(strsplit("add,rec,dom", ","))


#filtered <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_21_to_22.txt.gz"
#filtered_males <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_21_to_22.txt.gz"
#filtered_females <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_21_to_22.txt.gz"
#filtered_x <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_21_to_22.txt.gz"
#output_condensed <-
#  "ALLERGIC_RHINITIS_uk10k_condensed_chr_21_to_23_ref.txt.gz"
#output_tophits <-
#  "tophits_ALLERGIC_RHINITIS_GERA_300_uk10k_ref.txt.gz"
#output_ranges <-
#  "tophits_ALLERGIC_RHINITIS_GERA_300_uk10k_crossmodel_ranges_ref.txt.gz"
#pval_threshold <- "0.05"
#inheritance_models <- unlist(strsplit("add,rec,dom", ","))

#filtered <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_23_males.txt.gz"
#filtered_males <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_23_males.txt.gz"
#filtered_females <-
#  "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_23_females.txt.gz"
#filtered_x <- "ALLERGIC_RHINITIS_uk10k_filteredByAll_chr_23.txt.gz"
#output_condensed <-
#  "ALLERGIC_RHINITIS_uk10k_condensed_chr_21_to_23_ref.txt.gz"
#output_tophits <-
#  "tophits_ALLERGIC_RHINITIS_GERA_300_uk10k_ref.txt.gz"
#output_ranges <-
#  "tophits_ALLERGIC_RHINITIS_GERA_300_uk10k_crossmodel_ranges_ref.txt.gz"
#pval_threshold <- "0.05"
#inheritance_models <- unlist(strsplit("add,rec,dom", ","))

condensed_models_all <- NULL
pval <- as.numeric(pval_threshold)

init_con_time <- Sys.time()

if (filtered == filtered_males &
    filtered_males == filtered_females &
    filtered_females == filtered_x) {
  filtered_auto <-
    read.table(
      filtered,
      header = T,
      stringsAsFactors = F,
      na.strings = c("-", "NA")
    )
  selected_columns <- list(
    c(
      "chr",
      "position",
      "rs_id_all",
      "info_all",
      "alleleA",
      "alleleB",
      "all_maf",
      "refpanel"
    )
  )
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_auto)[which(apply(data.frame(colnames(filtered_auto)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
    
  }
  condensed <- filtered_auto[, selected_columns]
  condensed$info_all <- as.numeric(condensed$info_all)
}

if (filtered == filtered_males &
    filtered_males != filtered_females &
    filtered_females != filtered_x) {
  files_to_load <-
    c(filtered, filtered_males, filtered_females, filtered_x)
  
  read_file_function <- function(x) {
    read.table(
      x,
      header = T,
      stringsAsFactors = F,
      na.strings = c("-", "NA")
    )
  }
  
  amountOfAvailableCores = detectCores()
  files_read <- mcmapply(read_file_function,
                         files_to_load,
                         mc.cores = amountOfAvailableCores)
  
  filtered_x_males <- files_read[[filtered_males]]
  filtered_x_females <- files_read[[filtered_females]]
  filtered_x_all <- files_read[[filtered_x]]
  
  changeName <- function(condensed) {
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_se_", x) & grepl("sex.1", x)) == TRUE)] <- paste("frequentist_",m,"_se_sex.1", sep="")
    }
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_se_", x) & grepl("sex.2", x)) == TRUE)] <- paste("frequentist_", m, "_se_sex.2", sep="")
    }
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_beta_", x) & grepl("sex.1", x)) == TRUE)] <- paste("frequentist_", m, "_beta_sex.1", sep="")
    }
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_beta_", x) & grepl("sex.2", x)) == TRUE)] <- paste("frequentist_", m, "_beta_sex.2", sep="")
    }
    return(colnames(condensed))
  }
  
  colnames(filtered_x_males) <- changeName(filtered_x_males)
  colnames(filtered_x_females) <- changeName(filtered_x_females)
  colnames(filtered_x_all) <- changeName(filtered_x_all)
  
  selected_columns <- list(
    c(
      "chr",
      "position",
      "rs_id_all",
      "info_all",
      "alleleA",
      "alleleB",
      "all_maf",
      "refpanel"
    )
  )
  
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_x_males)[which(apply(data.frame(colnames(filtered_x_males)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
  }
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_x_females)[which(apply(data.frame(colnames(filtered_x_females)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
  }
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_x_all)[which(apply(data.frame(colnames(filtered_x_all)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
  }
  selected_columns <- unique(selected_columns)
  if (nrow(filtered_x_males) != 0) {
    filtered_x_males$chr <- as.character("23_males")
  }
  if (nrow(filtered_x_females) != 0) {
    filtered_x_females$chr <- as.character("23_females")
  }
  if (nrow(filtered_x_all) != 0) {
    filtered_x_all$chr <- as.character("23")
  }
  condensed <- rbind.fill(filtered_x_males[, intersect(selected_columns, colnames(filtered_x_males))],
                          filtered_x_females[, intersect(selected_columns, colnames(filtered_x_females))],
                          filtered_x_all[, intersect(selected_columns, colnames(filtered_x_all))])
  
  condensed$info_all <- as.numeric(condensed$info_all)
  condensed <- condensed[order(condensed$chr), ]
  
}


if (filtered != filtered_males &
    filtered_males != filtered_females &
    filtered_females != filtered_x) {
  files_to_load <-
    c(filtered, filtered_males, filtered_females, filtered_x)
  
  read_file_function <- function(x) {
    read.table(
      x,
      header = T,
      stringsAsFactors = F,
      na.strings = c("-", "NA")
    )
  }
  
  amountOfAvailableCores = detectCores()
  files_read <- mcmapply(read_file_function,
                         files_to_load,
                         mc.cores = amountOfAvailableCores)
  
  filtered_auto <- files_read[[filtered]]
  filtered_x_males <- files_read[[filtered_males]]
  filtered_x_females <- files_read[[filtered_females]]
  filtered_x_all <- files_read[[filtered_x]]
  
  changeName <- function(condensed) {
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_se_", x) & grepl("sex.1", x)) == TRUE)] <- paste("frequentist_",m,"_se_sex.1", sep="")
    }
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_se_", x) & grepl("sex.2", x)) == TRUE)] <- paste("frequentist_", m, "_se_sex.2", sep="")
    }
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_beta_", x) & grepl("sex.1", x)) == TRUE)] <- paste("frequentist_", m, "_beta_sex.1", sep="")
    }
    for (m in inheritance_models) {
      colnames(condensed)[which(apply(data.frame(colnames(condensed)), 1, function(x)
        grepl(m, x) & grepl("_beta_", x) & grepl("sex.2", x)) == TRUE)] <- paste("frequentist_", m, "_beta_sex.2", sep="")
    }
    return(colnames(condensed))
  }
  
  colnames(filtered_auto) <- changeName(filtered_auto)
  colnames(filtered_x_males) <- changeName(filtered_x_males)
  colnames(filtered_x_females) <- changeName(filtered_x_females)
  colnames(filtered_x_all) <- changeName(filtered_x_all)
  
  selected_columns <- list(
    c(
      "chr",
      "position",
      "rs_id_all",
      "info_all",
      "alleleA",
      "alleleB",
      "all_maf",
      "refpanel"
    )
  )
  
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_auto)[which(apply(data.frame(colnames(filtered_auto)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
  }
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_x_males)[which(apply(data.frame(colnames(filtered_x_males)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
  }
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_x_females)[which(apply(data.frame(colnames(filtered_x_females)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
  }
  for (m in inheritance_models) {
    selected_columns_new_cols <-
      colnames(filtered_x_all)[which(apply(data.frame(colnames(filtered_x_all)), 1, function(x)
        (grepl("_se_", x) |
           grepl("pvalue", x) |
           grepl("_beta_", x)) & grepl(m, x)) == TRUE)]
    selected_columns <-
      unlist(c(selected_columns, selected_columns_new_cols))
  }
  selected_columns <- unique(selected_columns)
  if (nrow(filtered_x_males) != 0) {
    filtered_x_males$chr <- as.character("23_males")
  }
  if (nrow(filtered_x_females) != 0) {
    filtered_x_females$chr <- as.character("23_females")
  }
  if (nrow(filtered_x_all) != 0) {
    filtered_x_all$chr <- as.character("23")
  }
  condensed <- rbind.fill(filtered_auto[, intersect(selected_columns, colnames(filtered_auto))],
                     filtered_x_males[, intersect(selected_columns, colnames(filtered_x_males))],
                     filtered_x_females[, intersect(selected_columns, colnames(filtered_x_females))],
                     filtered_x_all[, intersect(selected_columns, colnames(filtered_x_all))])

  condensed$info_all <- as.numeric(condensed$info_all)
  condensed <- condensed[order(condensed$chr), ]
}

#filtering info because the old GERA analysis was made for a info_score=0.5
#condensed <- condensed[condensed$info_all >= 0.7, ]

#for (n in 1:nrow(condensed)){
#    condensed$best_model <- lapply(strsplit(strsplit(names(condensed)[grep("_pvalue",names(condensed))[which.min(condensed[,grep("_pvalue",names(condensed))])]],"_pvalue")[[1]][1],"_")[[1]][2])
#}

# names(condensed)

write.table(
  condensed,
  output_condensed,
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

end_con_time <- Sys.time()
print(paste("Condensed file generation:" , end_con_time - init_con_time))

init_top_time <- Sys.time()

condensed$best_model <- NA
tophits <- NULL
tophits_models <- NULL
tophits_models_gen <- NULL
for (m in inheritance_models) {
  tophits_models <-
    condensed[which(as.numeric(condensed[, c(paste("frequentist_", m, "_pvalue", sep =
                                                     ""))]) <= pval &
                      !is.na(condensed[, c(paste("frequentist_", m, "_pvalue", sep = ""))])), ]
  tophits <- rbind.fill(tophits, tophits_models)
}

#head(tophits)
pvalues_index <- grep("pvalue", names(tophits))
pvalues_names <- names(tophits)[pvalues_index]
pvalues_as_columns <- tophits[,pvalues_index]
index_min_value <- unlist(lapply(transpose(pvalues_as_columns), function(x) which.min(x)))
amountOfAvailableCores = detectCores()
best_model_column <- unlist(mclapply(index_min_value, function(x) strsplit(pvalues_names[x], "_")[[1]][2], mc.cores=amountOfAvailableCores))
tophits$best_model <- best_model_column

#if (nrow(tophits) != 0) {
#  for (n in 1:nrow(tophits)) {
#    tophits$best_model[n] <-
#      strsplit(strsplit(names(tophits)[grep("_pvalue", names(tophits))[which.min(tophits[n, grep("_pvalue", names(tophits))])]], "_pvalue")[[1]][1], "_")[[1]][2]
#  }
#}


tophits <- tophits[!duplicated(tophits), ]
#table(tophits$chr)

write.table(
  tophits,
  output_tophits,
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

tophits$start <- tophits$position - 250000

end_top_time <- Sys.time()
print(paste("Top hits file generation:" , end_top_time - init_top_time))

init_cross_time <- Sys.time()

chr <- as.vector(unique(tophits$chr))

ranges <- NULL
for (i in chr) {
  dat <- tophits[tophits$chr == i, ]
  dat <- IRanges(start = dat$start, width = 500000)
  dat <- reduce(dat)
  dat <- as.data.frame(dat)
  dat$chr <- i
  dat$inheritance_models <- NA
  for (n in 1:nrow(dat)) {
    dat$num_variants[n] <-
      dim(tophits[(tophits$position >= dat$start[n] &
                     tophits$position <= dat$end[n] &
                     tophits$chr == i), ])[1]
    if (i != "23" & i != "23_males" & i != "23_females") {
      for (inh in inheritance_models) {
        if (any(
          tophits$position >= dat$start[n] & tophits$position <= dat$end[n] &
          tophits$chr == i &
          tophits[, c(paste("frequentist_", inh, "_pvalue", sep = ""))] <= pval &
          !(is.na(tophits[, c(paste("frequentist_", inh, "_pvalue", sep =
                                    ""))]))
        )) {
          if (is.na(dat$inheritance_models[n])) {
            dat$inheritance_models[n] <- paste(inh)
          } else {
            dat$inheritance_models[n] <-
              paste(dat$inheritance_models[n], inh, sep = ",")
          }
        }
      }
    } else {
      dat$inheritance_models[n] <- paste("add")
    }
  }
  ranges <- rbind(ranges, dat)
}

if (nrow(tophits) != 0) {
  ranges$range <-
    paste("chr", ranges$chr, ":", ranges$start, "-", ranges$end, sep = "")
} else {
  ranges <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("chr", "start", "end", "width")
  colnames(ranges) <- x
}

#head(ranges)
#dim(ranges)

tophits_all <- tophits
tophits_ranges <- ranges
topvariant <- NULL
topvariants <- NULL
tophits_topvariants_auto <- NULL
topvariant_final <- NULL

if (nrow(tophits_ranges) != 0) {
  for (i in 1:length(tophits_ranges$range)) {
    if (tophits_ranges$chr[i] != "23" &
        tophits_ranges$chr[i] != "23_males" &
        tophits_ranges$chr[i] != "23_females") {
      tophits_set <-
        tophits_all[(
          tophits_all$chr == tophits_ranges$chr[i] &
            tophits_all$position >= tophits_ranges$start[i] &
            tophits_all$position <= tophits_ranges$end[i]
        ), ]
      #print(head(tophits_set))
      #print(table(duplicated(tophits_set$rs_id_all)))
      for (m in unlist(strsplit(as.character(
        tophits_ranges$inheritance_models[i]
      ), ","))) {
        if (m != "gen") {
          if (any(!is.na((tophits_set[, c(paste("frequentist_", m, "_pvalue", sep =
                                                ""))])))) {
            topvariant <-
              as.data.frame(tophits_set[which.min(tophits_set[, c(paste("frequentist_", m, "_pvalue", sep =
                                                                          ""))]),
                                        c(
                                          "position",
                                          "rs_id_all",
                                          "alleleA",
                                          "alleleB",
                                          "all_maf",
                                          paste("frequentist_", m, "_pvalue", sep = ""),
                                          paste("frequentist_", m, "_se_1", sep =
                                                  ""),
                                          paste("frequentist_", m, "_beta_1", sep = ""),
                                          "info_all",
                                          "refpanel"
                                        )])
            names(topvariant)[6] <- "pvalue"
            names(topvariant)[7] <- "se"
            names(topvariant)[8] <- "beta_1"
            topvariant$model <- m
            topvariant$range <- i
            if ("add" %in% inheritance_models) {
              topvariant$add_pvalue <-
                tophits_set[topvariant$rs_id_all == tophits_set$rs_id_all, c("frequentist_add_pvalue")]
            }
          }
        } else {
          if (any(!is.na((tophits_set[, c(paste("frequentist_", m, "_pvalue", sep =
                                                ""))])))) {
            topvariant <-
              as.data.frame(tophits_set[which.min(tophits_set[, c(paste("frequentist_", m, "_pvalue", sep =
                                                                          ""))]),
                                        c(
                                          "position",
                                          "rs_id_all",
                                          "alleleA",
                                          "alleleB",
                                          "all_maf",
                                          paste("frequentist_", m, "_pvalue", sep =
                                                  ""),
                                          paste("frequentist_", m, "_se_1", sep =
                                                  ""),
                                          paste("frequentist_", m, "_beta_1", sep =
                                                  ""),
                                          paste("frequentist_", m, "_beta_2", sep =
                                                  ""),
                                          "info_all",
                                          "refpanel"
                                        )])
            names(topvariant)[6] <- "pvalue"
            names(topvariant)[7] <- "se"
            names(topvariant)[8] <- "beta_1"
            names(topvariant)[9] <- "beta_2"
            topvariant$model <- m
            topvariant$range <- i
            if ("add" %in% inheritance_models) {
              topvariant$add_pvalue <-
                tophits_set[topvariant$rs_id_all == tophits_set$rs_id_all, c("frequentist_add_pvalue")]
            }
          }
        }
        topvariants <- rbind.fill(topvariants, topvariant)
      }
      topvariants$TOPHIT <- "NO"
      for (r in unique(topvariants$range)) {
        topvariants$TOPHIT[which(topvariants$range == r)[which.min(topvariants[topvariants$range ==
                                                                                 r, c("pvalue")])]] <-
          "YES"
        
      }
      tophits_topvariants_auto <-
        topvariants[topvariants$TOPHIT == "YES", ]
    }
  }
}

topvariant <- NULL
topvariants <- NULL
tophits_topvariants_x <- NULL
topvariant_final <- NULL

if (nrow(tophits_ranges) != 0) {
  for (i in 1:length(tophits_ranges$range)) {
    if (tophits_ranges$chr[i] == "23" |
        tophits_ranges$chr[i] == "23_males" |
        tophits_ranges$chr[i] == "23_females") {
      tophits_set <-
        tophits_all[(
          as.character(tophits_all$chr) == tophits_ranges$chr[i] &
            tophits_all$position >= tophits_ranges$start[i] &
            tophits_all$position <= tophits_ranges$end[i]
        ), ]
      for (m in unlist(strsplit(as.character(
        tophits_ranges$inheritance_models[i]
      ), ","))) {
        if (any(!is.na((tophits_set[, c(paste("frequentist_", m, "_pvalue", sep =
                                              ""))])))) {
          topvariant <-
            as.data.frame(tophits_set[which.min(tophits_set[, c(paste("frequentist_", m, "_pvalue", sep =
                                                                        ""))]),
                                      c(
                                        "position",
                                        "rs_id_all",
                                        "alleleA",
                                        "alleleB",
                                        "all_maf",
                                        paste("frequentist_", m, "_pvalue", sep =
                                                ""),
                                        paste("frequentist_", m, "_se_sex.1", sep =
                                                ""),
                                        paste("frequentist_", m, "_beta_sex.1", sep =
                                                ""),
                                        paste("frequentist_", m, "_se_sex.2", sep =
                                                ""),
                                        paste("frequentist_", m, "_beta_sex.2", sep =
                                                ""),
                                        "info_all",
                                        "refpanel"
                                      )])
          names(topvariant)[6] <- "pvalue"
          names(topvariant)[7] <- "se_males"
          names(topvariant)[8] <- "beta_males"
          names(topvariant)[9] <- "se_females"
          names(topvariant)[10] <- "beta_females"
          topvariant$range <- i
          topvariant$model <- m
          topvariant$add_pvalue <-
            tophits_set[topvariant$rs_id_all == tophits_set$rs_id_all, c("frequentist_add_pvalue")]
        }
        topvariants <- rbind.fill(topvariants, topvariant)
      }
      topvariants$TOPHIT <- "NO"
      for (r in unique(topvariants$range)) {
        topvariants$TOPHIT[which(topvariants$range == r)[which.min(topvariants[topvariants$range ==
                                                                                 r, c("pvalue")])]] <-
          "YES"
        
      }
      tophits_topvariants_x <-
        topvariants[topvariants$TOPHIT == "YES", ]
    }
  }
}

if (is.null(tophits_topvariants_auto)) {
  tophits_topvariants_auto  <- data.frame(matrix(ncol = 14, nrow = 0))
  names(tophits_topvariants_auto) <-
    c(
      "position",
      "rs_id_all",
      "alleleA",
      "alleleB",
      "all_maf",
      "pvalue",
      "se_1",
      "se_2",
      "beta_1",
      "beta_2",
      "model",
      "add_pvalue",
      "info_all",
      "refpanel"
    )
}

if (is.null(tophits_topvariants_x)) {
  tophits_topvariants_x  <- data.frame(matrix(ncol = 14, nrow = 0))
  names(tophits_topvariants_x) <-
    c(
      "position",
      "rs_id_all",
      "alleleA",
      "alleleB",
      "all_maf",
      "pvalue",
      "se_males",
      "beta_males",
      "se_females",
      "beta_females",
      "model",
      "add_pvalue",
      "info_all",
      "refpanel"
    )
}

if (nrow(tophits_ranges) == 0) {
  tophits_ranges <- NULL
  tophits_ranges  <- data.frame(matrix(ncol = 7, nrow = 0))
  names(tophits_ranges) <-
    c("start",
      "end",
      "width",
      "chr",
      "inheritance_models",
      "num_variants",
      "range")
}

tophits_topvariants <-
  rbind.fill(tophits_topvariants_auto, tophits_topvariants_x)
tophits_final <- cbind(tophits_ranges, tophits_topvariants)

if (nrow(tophits_ranges) != 0) {
  tophits_final <-
    tophits_final[, -c(which(names(tophits_final) == "TOPHIT"), which(names(tophits_final) ==
                                                                        "range"))]
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

end_cross_time <- Sys.time()
print(paste("Cross model file generation:" , end_cross_time - init_cross_time))

start_write_time <- Sys.time()

write.table(
  tophits_final,
  output_ranges,
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

end_write_time <- Sys.time()
print(paste("Write time:" , end_write_time - start_write_time))