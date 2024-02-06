library(tximport)
library(GenomicFeatures)

BiocManager::install("tximport")


  setwd("~/PreProcSEQ-main/5-expressionMatrix/tximport")
dirquant <- "~/PreProcSEQ-main/4-quantification/salmon/quant_salmon/"
files <- list.files(dirquant)
files_import <- paste0(dirquant, files[-1], "/quant.sf")
all(file.exists(files_import))