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
samples <- read.delim("RNAseq/quant_salmon/SLC_HFD_Females.tsv",header=FALSE)
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
mart <- useDataset("mmc57bl6nj_gene_ensembl", useMart("ENSEMBL_MART_MOUSE", host="https://www.ensembl.org"))

# Get official gene symbol and Ensembl gene IDs
tx2gene <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id", "external_gene_name"),
                 values= ids.strip,
                 mart= mart)


df <- tx2gene%>% mutate_all(~na_if(., ""))
tx2gene_noNA <- df %>% drop_na(external_gene_name)
txiSLCHFD <- tximport(files, type="salmon", txIn = TRUE, txOut = FALSE, tx2gene=tx2gene_noNA, ignoreTxVersion=TRUE)

# Save files 
save(txiSLCHFD, file = "RNAseq/starting_files/females_SLC_HFD_matrix_tximport_salmon.Rdata") # objeto R
write.csv(txiSLCHFD$abundance, file = "RNAseq/starting_files/females_SLC_HFD_matrix_salmon_tximport_abundance.csv") # TPM
write.csv(txiSLCHFD$counts, file = "RNAseq/starting_files/females_SLC_HFD_matrix_salmon_tximport_counts.csv") # counts

## Run DESeq2
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(DESeq2)
library(plotly)
library(cowplot)

metadata <- read.delim("RNAseq/starting_files/All_Metadata.tsv", sep="\t")
metadata <- metadata %>% filter(Model=="HFD") %>% filter(Sex=="Female")
row.names(metadata) <- metadata$SampleID
colnames(txiSLCHFD$counts)
row.names(metadata) == colnames(txiSLCHFD$counts)


dds = DESeqDataSetFromTximport(txiSLCHFD, colData = metadata, design = ~Genotype)
dds <- DESeq(dds)

CDvsNorm=results(dds, contrast=c("Genotype", "MUT", "WT"))

# for continuous variable: CDvsNorm = results(diagadds, name="RQ_T2R138")
head(CDvsNorm)
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix <- cbind(as(CDvsNorm, "data.frame"))
head(CDvsNormMatrix)

write.csv(CDvsNormMatrix,"RNAseq/females_DESEQ2_HFD_MUT_vs_WT_results.csv")

plot <- read.csv("RNAseq/females_DESEQ2_HFD_MUT_vs_WT_results.csv",row.names = 1)
summary(plot$padj)
plot <- plot #remove unannotated outlier
female_hfd_plot <- EnhancedVolcano(plot,
                                 lab = rownames(plot),
                                 labSize = 3,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 ylab = bquote(~-Log[10]~ '(p-adjusted)'),
                                 pCutoff = 0.05,
                                 title = "Females: SLC HFD MUT vs WT",
                                 subtitle = "Gene~ Genotype",
                                 FCcutoff = 1)+
  ylim(0,2.5)+
  theme_cowplot(12)+
  theme(legend.position = "top")
female_hfd_plot
ggsave(female_hfd_plot, file="RNAseq/female_hfd_plot.png",width = 7,height=7)

CDvsNormMatrix <- CDvsNormMatrix %>% filter(padj<0.05)
write.csv(CDvsNormMatrix,"RNAseq/females_significant_DESEQ2_HFD_MUT_vs_WT_results.csv")

### HETs
CDvsNorm=results(dds, contrast=c("Genotype", "HET", "WT"))

# for continuous variable: CDvsNorm = results(diagadds, name="RQ_T2R138")
head(CDvsNorm)
CDvsNorm = CDvsNorm[order(CDvsNorm$padj, na.last = NA), ]
CDvsNormMatrix <- cbind(as(CDvsNorm, "data.frame"))
head(CDvsNormMatrix)

write.csv(CDvsNormMatrix,"RNAseq/females_DESEQ2_HFD_HET_vs_WT_results.csv")

plot <- read.csv("RNAseq/females_DESEQ2_HFD_HET_vs_WT_results.csv",row.names = 1)
summary(plot$padj)
plot <- plot #remove unannotated outlier
het_female_hfd_plot <- EnhancedVolcano(plot,
                                   lab = rownames(plot),
                                   labSize = 3,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   ylab = bquote(~-Log[10]~ '(p-adjusted)'),
                                   pCutoff = 0.05,
                                   title = "Females: SLC HFD HET vs WT",
                                   subtitle = "Gene~ Genotype",
                                   FCcutoff = 1)+
  ylim(0,2.5)+
  theme_cowplot(12)+
  theme(legend.position = "top")
het_female_hfd_plot
ggsave(het_female_hfd_plot, file="RNAseq/het_female_hfd_plot.png",width = 7,height=7)

CDvsNormMatrix <- CDvsNormMatrix %>% filter(padj<0.05)
write.csv(CDvsNormMatrix,"RNAseq/females_significant_DESEQ2_HFD_HET_vs_WT_results.csv")



