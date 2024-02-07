BiocManager::install("tximeta")

library(tximeta)
library(SummarizedExperiment)
library(readxl)
library(GenomicFeatures)
library(BiocFileCache)

setwd("/home/julianne/Documents/slc_si_paper/")

## add annotation for the geneIDs
dirquant <- "RNAseq/quant_salmon/"

coldata <- read.csv("RNAseq/starting_files/Metadata.csv")
coldata$names <- coldata$SampleID
coldata$files <- paste0(dirquant,coldata$Filepath,"/","quant.sf")
all(file.exists(coldata$files))
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$names

# point to a Salmon quantification file with an additional artificial transcript
dir <- "RNAseq/quant_salmon/trim_JJ1715_101_S25_R1_001.fastq_paired.fq_quant/"
file <- file.path(dir, "quant.sf")
coldata <- data.frame(files=file, names="JJ1715_101_S25", sample="1",
                      stringsAsFactors=FALSE)

# now point to the Salmon index itself to create a linkedTxome
# as the index will not match a known txome
indexDir <- file.path("/home/julianne/Documents/transcriptomicsonhoffman/Mus_musculus_c57bl6nj_index/")

# point to the source FASTA and GTF:
fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-111/fasta/mus_musculus_c57bl6nj/cdna/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.cdna.all.fa.gz")
gtfPath <- file.path("/home/julianne/Documents/transcriptomicsonhoffman/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.111.gtf.gz")

# now create a linkedTxome, linking the Salmon index to its FASTA and GTF sources
makeLinkedTxome(indexDir=indexDir, source="Ensembl", organism="Mus musculus",
                release="111", genome="C57BL_6NJ_v1.111", fasta=fastaFTP, gtf=gtfPath, write=FALSE)

# add annotations
se <- tximeta(coldata, useHub = F)
gse <- summarizeToGene(se)
save(gse, file="RNAseq/matrix_gse_salmon_tximeta.RData")

se(tximeta())
