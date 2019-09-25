manhattan <-
  function(dataframe,
           colors = c("lightsteelblue4", "lightyellow2"),
           ymin = 0,
           ymax = "max",
           limitchromosomes = 1:25,
           suggestiveline = -log10(1e-5),
           genomewideline = -log10(5e-8),
           annotate = NULL,
           ...) {
    d = dataframe
    #Validate input files
    if (!("chr" %in% names(d) &
          "position" %in% names(d) &
          "pvalue" %in% names(d)))
      stop("Make sure your data frame contains columns chr, position, and pvalue")
    
    if (any(limitchromosomes))
      d = d[d$chr %in% limitchromosomes,]
    d = subset(na.omit(d[order(d$chr, d$position),]), (pvalue > 0 &
                                                         pvalue <= 1)) # remove na's, sort, and keep only 0<pvalue<=1
    d$logp = -log10(d$pvalue)
    d$pos = NA
    ticks = NULL
    lastbase = 0
    colors <- c(rep(colors, 24)[1:24], "lightsteelblue4")
    
    if (ymax == "max") {
      ymax <- ceiling(max(d$logp))
    }
    ymax <- as.numeric(ymax)
    if (ymax < 8) {
      ymax <- 8
    }
    
    #Print Max and Min P-values
    print(paste("YMAX:", ymax))
    print(paste("YMIN:", ymin))
    
    numchroms = length(unique(d$chr))
    
    if (numchroms == 1) {
      d$pos = d$position
      ticks = floor(length(d$pos)) / 2 + 1
    } else {
      for (i in unique(d$chr)) {
        if (i == as.vector(unique(d$chr))[1]) {
          d[d$chr == i,]$pos = d[d$chr == i,]$position
        } else {
          lastbase = lastbase + tail(subset(d, chr == i - 1)$position, 1)
          d[d$chr == i,]$pos = d[d$chr == i,]$position +
            lastbase
        }
        ticks = c(ticks, d[d$chr == i,]$pos[floor(length(d[d$chr ==
                                                             i,]$pos) / 2) + 1])
      }
    }
    
    chr <- unique(d$chr)
    if (23 %in% chr) {
      ticks[(length(ticks))]
      ticks[(length(ticks) - 1)]
      ticks2 <- c(ticks[(length(ticks))], ticks[(length(ticks) - 1)])
      
      d$chr[d$chr == "23"] <- as.character("F")
      d$chr[d$chr == "24"] <- as.character("M")
    }
    
    if (numchroms == 1) {
      plot(
        d$pos,
        d$logp,
        ylim = c(ymin, ymax),
        ylab = expression(-log[10](italic(p))),
        xlab = paste("Chromosome", unique(d$chr), "position")
      )
    } else {
      plot(
        d$pos,
        d$logp,
        ylim = c(ymin, ymax),
        ylab = expression(-log[10](italic(p))),
        xlab = "Chromosome",
        xaxt = "n",
        type = "n"
      )
      axis(1, at = ticks, lab = unique(d$chr))
      
      if (23 %in% chr) {
        axis(
          1,
          at = ticks2,
          line = 2.5,
          tick = T,
          lab = rep("", 2),
          lwd.ticks = 0
        )
        axis(
          1,
          at = (as.numeric(ticks2[1]) + as.numeric(ticks2[2])) / 2,
          lab = c("X"),
          line = 2,
          tick = F
        )
      }
      
      icol = 1
      
      for (i in unique(d$chr)) {
        #with(d[d$chr==i, ],points(pos[d$pvalue > 5e-8], logp[d$pvalue > 5e-8], col=colors[icol], ...))
        points(d$pos[d$chr == i], d$logp[d$chr == i], col =
                 colors[icol], pch = ".")
        icol = icol + 1
      }
      for (i in unique(d$chr)) {
        #with(d[d$chr==i, ],points(pos[d$pvalue > 5e-8], logp[d$pvalue > 5e-8], col=colors[icol], ...))
        points(d$pos[d$chr == i &
                       d$logp >= -log10(5e-8)],
               d$logp[d$chr == i &
                        d$logp >= -log10(5e-8)],
               col = "lightcoral",
               cex = 0.6,
               pch = 16)
      }
    }
    
    
    genomewideline <- -log10(5e-8)
    abline(h = genomewideline, col = "indianred2", lty = 2)
    
    
  }

################################################################################################################

.libPaths("/gpfs/projects/bsc05/silvia/R_INSTALL/R_libs")

library(gap)
library(sfsmisc)


args <- commandArgs(TRUE)

#tab_file <- args[1] #tab file with chr, position , pvalue
tab_file <-
  "/home/ramela/git/guidance/CANCER_uk10k_condensed_chr_1_to_23.txt.gz"

#out_qqplot <- args[2] #name out qqplot qqplot_namestudy.pdf
#out_manhattan <- args[3] #name out manhattan manhattan_namestudy.pdf

out_qqplot <- "/home/ramela/git/guidance/qqplot_namestudy.pdf"
out_manhattan <-
  "/home/ramela/git/guidance/manhattan_namestudy.pdf"

#out_qqplot_tiff <- args[4] #Name output file of qqplot on tiff.
#out_manhattan_tiff <- args[5] #name out manhattan manhattan_namestudy.pdf

out_qqplot_tiff <-
  "/home/ramela/git/guidance/qqplot_namestudy.tiff"
out_manhattan_tiff <-
  "/home/ramela/git/guidance/manhattan_namestudy.tiff"

tab_file_data <- read.delim(tab_file)
tab_file_data <-
  tab_file_data[tab_file_data$frequentist_add_pvalue != 0, ]

p <- tab_file_data$frequentist_add_pvalue
p <- p[!is.na(p)]
n <- length(p)
x2obs <- qchisq(p, 1, lower.tail = FALSE)
x2exp <- qchisq(1:n / n, 1, lower.tail = FALSE)

#make Manhattan on pdf

tab_file_data$chr <- as.character(tab_file_data$chr)
tab_file_data$chr[tab_file_data$chr == "23_females"] <-
  as.character("23")
tab_file_data$chr[tab_file_data$chr == "23_males"] <-
  as.character("24")
tab_file_data$chr <- as.numeric(tab_file_data$chr)

tab_file_data$SNP <-
  paste(paste("chr", tab_file_data$chr, sep = ""),
        tab_file_data$position,
        sep = ":")
tab_man <- data.frame(
  SNP = tab_file_data$SNP,
  chr = tab_file_data$chr,
  position = tab_file_data$position,
  pvalue = tab_file_data$frequentist_add_pvalue
)

tab_man <- tab_man[tab_man$pvalue <= 0.05, ]
dim(tab_man)

title <- "Manhattan-plot"
YMIN <- 0
YMAX <- "max"





############################# INIT FUNCTION

#dataframe, colors=c("lightsteelblue4","lightyellow2"),
#ymin=0, ymax="max", limitchromosomes=1:25,
#suggestiveline=-log10(1e-5),
#genomewideline=-log10(5e-8),
#annotate=NULL, ...) {

#manhattan(tab_man, pch=16,cex=0.70,main=title,colors=c("lightsteelblue4","ivory3"), suggestiveline=-log10(5e-8), ymax=YMAX, ymin=YMIN)

d = tab_man
pch = 16
cex = 0.7
suggestiveline = -log10(1e-5)
genomewideline = -log10(5e-8)
limitchromosomes = 1:25
ymax = YMAX
ymin = YMIN
colors = c("lightsteelblue4", "ivory3")

#####################333


#Validate input files
if (!("chr" %in% names(d) &
      "position" %in% names(d) &
      "pvalue" %in% names(d)))
  stop("Make sure your data frame contains columns chr, position, and pvalue")

if (any(limitchromosomes))
  d = d[d$chr %in% limitchromosomes,]
d = subset(na.omit(d[order(d$chr, d$position),]), (pvalue > 0 &
                                                     pvalue <= 1)) # remove na's, sort, and keep only 0<pvalue<=1
d$logp = -log10(d$pvalue)
d$pos = NA
ticks = NULL
lastbase = 0
colors <- c(rep(colors, 24)[1:24], "lightsteelblue4")

if (ymax == "max") {
  ymax <- ceiling(max(d$logp))
}
ymax <- as.numeric(ymax)
if (ymax < 8) {
  ymax <- 8
}

#Print Max and Min P-values
print(paste("YMAX:", ymax))
print(paste("YMIN:", ymin))

numchroms = length(unique(d$chr))

if (numchroms == 1) {
  d$pos = d$position
  ticks = floor(length(d$pos)) / 2 + 1
} else {
  unique_chr = as.vector(unique(d$chr))
  for (i in 1:numchroms) {
    if (i == 1) {
      d[d$chr == unique_chr[i], ]$pos = d[d$chr == unique_chr[i], ]$position
    } else {
      lastbase = lastbase + tail(subset(d, chr == unique_chr[i - 1])$position, 1)
      d[d$chr == unique_chr[i], ]$pos = d[d$chr == unique_chr[i], ]$position +
        lastbase
    }
    ticks = c(ticks, d[d$chr == unique_chr[i], ]$pos[floor(length(d[d$chr ==
                                                                      unique_chr[i], ]$pos) / 2) + 1])
  }
}

chr <- unique(d$chr)
if (23 %in% chr) {
  #ticks[(length(ticks))]
  #ticks[(length(ticks)-1)]
  ticks2 <- c(ticks[(length(ticks))], ticks[(length(ticks) - 1)])
  
  d$chr[d$chr == "23"] <- as.character("F")
  d$chr[d$chr == "24"] <- as.character("M")
}

if (numchroms == 1) {
  plot(
    d$pos,
    d$logp,
    ylim = c(ymin, ymax),
    ylab = expression(-log[10](italic(p))),
    xlab = paste("Chromosome", unique(d$chr), "position")
  )
} else {
  plot(
    d$pos,
    d$logp,
    ylim = c(ymin, ymax),
    ylab = expression(-log[10](italic(p))),
    xlab = "Chromosome",
    xaxt = "n",
    type = "n"
  )
  axis(1, at = ticks, lab = unique(d$chr))
  
  if (23 %in% chr) {
    axis(
      1,
      at = ticks2,
      line = 2.5,
      tick = T,
      lab = rep("", 2),
      lwd.ticks = 0
    )
    axis(
      1,
      at = (as.numeric(ticks2[1]) + as.numeric(ticks2[2])) / 2,
      lab = c("X"),
      line = 2,
      tick = F
    )
  }
  
  icol = 1
  
  for (i in unique(d$chr)) {
    #with(d[d$chr==i, ],points(pos[d$pvalue > 5e-8], logp[d$pvalue > 5e-8], col=colors[icol], ...))
    points(d$pos[d$chr == i], d$logp[d$chr == i], col = colors[icol], pch =
             ".")
    icol = icol + 1
  }
  for (i in unique(d$chr)) {
    #with(d[d$chr==i, ],points(pos[d$pvalue > 5e-8], logp[d$pvalue > 5e-8], col=colors[icol], ...))
    points(d$pos[d$chr == i & d$logp >= -log10(5e-8)],
           d$logp[d$chr == i & d$logp >= -log10(5e-8)],
           col = "lightcoral",
           cex = 0.6,
           pch = 16)
  }
}


genomewideline <- -log10(5e-8)
abline(h = genomewideline, col = "indianred2", lty = 2)

#################### END FUNCTION



























pdf(out_manhattan, width = 14, height = 7.5)
manhattan(
  tab_man,
  pch = 16,
  cex = 0.70,
  main = title,
  colors = c("lightsteelblue4", "ivory3"),
  suggestiveline = -log10(5e-8),
  ymax = YMAX,
  ymin = YMIN
)
dev.off()

#make Manhattan on tiff
tiff(
  out_manhattan_tiff,
  width = 9600,
  height = 5600,
  units = "px",
  res = 800,
  compression = "lzw"
)
manhattan(
  tab_man,
  pch = 16,
  cex = 0.70,
  main = title,
  colors = c("lightsteelblue4", "ivory3"),
  suggestiveline = -log10(5e-8),
  ymax = YMAX,
  ymin = YMIN
)
dev.off()

#    write.table(tab_man,out_corrected_pvals,row.names=F,sep="\t",quote=F)

cat("MANHATTAN PLOT successfully completed!\n")
