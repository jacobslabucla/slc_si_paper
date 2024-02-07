## The following code is a combination of work from https://bookdown.org/jean_souza/PreProcSEQ/annotation-of-transcripts.html
## and work from https://github.com/hbctraining/Intro-to-rnaseq-hpc-gt/blob/master/lessons/DE_analysis.md
## with the exception of additional code that enables us to use the latest release of Ensembl Mus Musculus C57Bl6 Transcriptome 

BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")

#need development version of biocfilecache due to issues with dbplyr
remove.packages("BiocFileCache")
devtools::install_github("Bioconductor/BiocFileCache")


library(tximport)
library(GenomicFeatures)
library(biomaRt)
library(BiocFileCache)
library(dplyr)
library(tidyr)

setwd("/home/julianne/Documents/slc_si_paper/RNAseq/quant_salmon/")
## List all directories containing data  
samples <- list.files(path = ".", full.names = F, pattern="_quant$")
files <- file.path(samples, "quant.sf")
names(files) <-  samples
# Create a character vector of Ensembl IDs		

## List all directories containing data  
## Since all quant files have the same name it is useful to have names for each element
ids <- read.delim("trim_JJ1715_101_S25_R1_001.fastq_paired.fq_quant/quant.sf", sep="\t",header=T)
ids <- as.character(ids[,1])
head(ids)
require(stringr)
ids.strip <- str_replace(ids, "([.][0-9])", "")
head(ids.strip)
# Create a mart object

mart <- useDataset("mmc57bl6nj_gene_ensembl", useMart("ENSEMBL_MART_MOUSE", host="www.ensembl.org"))

# Get official gene symbol and Ensembl gene IDs
tx2gene <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id", "external_gene_name"),
  values= ids.strip,
  mart= mart)

?tximport()
df <- tx2gene%>% mutate_all(~na_if(., ""))
tx2gene_noNA <- df %>% drop_na(external_gene_name)
txi <- tximport(files, type="salmon", txIn = TRUE, txOut = FALSE, tx2gene=tx2gene_noNA, ignoreTxVersion=TRUE)

# Save files 
save(txi, file = "../matrix_tximport_salmon.Rdata") # objeto R
write.csv(txi$abundance, file = "../matrix_salmon_tximport_abundance.csv") # TPM
write.csv(txi$counts, file = "../matrix_salmon_tximport_counts.csv") # counts

df <- txi$abundance
summary(colSums(df))
