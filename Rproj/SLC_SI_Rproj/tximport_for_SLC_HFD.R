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

setwd("/home/julianne/Documents/slc_si_paper/")

## List all directories containing data  
samples <- read.delim("RNAseq/quant_salmon/SLC_HFD_Files.tsv",header=FALSE)
samples <- samples$V1
files <- file.path("RNAseq/quant_salmon",samples, "quant.sf")
names <-  gsub("_R1.*$","",samples)
names <- gsub("^.*JJ1715_", "", names)
names <-  gsub("_.*$","",names)
names(files) <- names

## Since all quant files have the same name it is useful to have names for each element
ids <- read.delim("RNAseq/quant_salmon/trim_JJ1715_153_S36_R1_001.fastq_paired.fq_quant/quant.sf", sep="\t",header=T)
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
txiSLCHFD <- tximport(files, type="salmon", txIn = TRUE, txOut = FALSE, tx2gene=tx2gene_noNA, ignoreTxVersion=TRUE)

# Save files 
save(txiSLCHFD, file = "RNAseq/starting_files/SLC_HFD_matrix_tximport_salmon.Rdata") # objeto R
write.csv(txiSLCHFD$abundance, file = "RNAseq/starting_files/SLC_HFD_matrix_salmon_tximport_abundance.csv") # TPM
write.csv(txiSLCHFD$counts, file = "RNAseq/starting_files/SLC_HFD_matrix_salmon_tximport_counts.csv") # counts

## Run DESeq2
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(DESeq2)
metadata <- read.delim("RNAseq/starting_files/All_Metadata.tsv", sep="\t")
metadata <- metadata %>% filter(Model=="HFD")
row.names(metadata) <- metadata$SampleID
colnames(txiSLCHFD$counts)
row.names(metadata) == colnames(txiSLCHFD$counts)


dds = DESeqDataSetFromTximport(txiSLCHFD, colData = metadata, design = ~Sex + Genotype)
dds <- DESeq(dds)

CDvsNorm=results(dds, contrast=c("Genotype", "MUT", "WT"))
# for continuous variable: CDvsNorm = results(diagadds, name="RQ_T2R138")
head(CDvsNorm)
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix <- cbind(as(CDvsNorm, "data.frame"))
head(CDvsNormMatrix)

write.csv(CDvsNormMatrix,"RNAseq/DESEQ2_HFD_MUT_vs_WT_results.csv")
plot <- CDvsNormMatrix# %>% filter(padj>1e-10)
hfd_plot <- EnhancedVolcano(plot,
                            lab = rownames(plot),
                            x = 'log2FoldChange',
                            y = 'padj',
                            ylab = bquote(~Log[10]~ '(p-adjusted)'),
                            pCutoff = 0.05,
                            title = "SLC HFD MUT vs WT",
                            subtitle = "Gene~ Sex + Genotype")
save(hfd_plot, file="RNAseq/hfd_plot.png")


CDvsNormMatrix <- CDvsNormMatrix %>% filter(padj<0.05)
write.csv(CDvsNormMatrix,"RNAseq/significant_DESEQ2_HFD_MUT_vs_WT_results.csv")


