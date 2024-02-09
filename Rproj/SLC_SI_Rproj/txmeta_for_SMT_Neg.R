BiocManager::install("tximeta")

library(tximeta)
library(SummarizedExperiment)
library(readxl)
library(GenomicFeatures)
library(BiocFileCache)
library(dplyr)

setwd("/home/julianne/Documents/slc_si_paper/")

## add annotation for the geneIDs
dirquant <- "RNAseq/quant_salmon/"

coldata <- read.delim("RNAseq/starting_files/All_Metadata.tsv",sep="\t")
coldata <- coldata %>% filter(Model=="SMT_Neg")
coldata$names <- coldata$SampleID
coldata$files <- paste0(dirquant,coldata$Filepath,"/","quant.sf")
all(file.exists(coldata$files))
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$names


# now point to the Salmon index itself to create a linkedTxome
# as the index will not match a known txome
indexDir <- file.path("/home/julianne/Documents/transcriptomicsonhoffman/Mus_musculus_c57bl6nj_index/")

# point to the source FASTA and GTF:
fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-111/fasta/mus_musculus_c57bl6nj/cdna/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.cdna.all.fa.gz")
gtfPath <- file.path("/home/julianne/Documents/transcriptomicsonhoffman/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.111.gtf.gz")

# now create a linkedTxome, linking the Salmon index to its FASTA and GTF sources
makeLinkedTxome(indexDir=indexDir, source="Ensembl", organism="Mus musculus",
                release="111", genome="C57BL_6NJ_v1.111", fasta=fastaFTP, gtf=gtfPath, write=FALSE)

# summarize transcripts to gene level 
se <- tximeta(coldata, useHub = F)
gse <- summarizeToGene(se)
save(gse, file="RNAseq/SMT_Neg_matrix_gse_salmon_tximeta.RData")
save(se, file = "RNAseq/SMT_Neg_matrix_se_salmon_tximeta.RData")

## Annotating the transcripts further(?)
library(AnnotationHub)
library(GenomicFeatures)
library(ensembldb)
library(SummarizedExperiment)
ah <- AnnotationHub()
ah <- subset(ah, species == "Mus musculus")

edb <- query(ah, "C57BL_6NJ")

gns <- ensembldb::genes(edb)

EnsDbAnnotation <- as.data.frame(gns)
EnsDbAnnotation <- EnsDbAnnotation[,c("gene_id","symbol","gene_biotype","entrezid")]
dim(EnsDbAnnotation)
colnames(EnsDbAnnotation) <- c("ensemblid","symbol","gene_biotype","entrezid")

?query()
  load("../5-expressionMatrix/tximeta/matrix_gse_salmon_tximeta.RData")

nrow(gse)
gseAnnotation <- rowData(gse)
## Run DESeq2 --
BiocManager::install("DESeq2")
library(DESeq2)
load(file="RNAseq/SMT_Neg_matrix_gse_salmon_tximeta.RData")
?DESeqDataSet()
colData(gse)
dds <- DESeqDataSet(gse, ~ Sex + Genotype)
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table
res <- res[order(res$padj),]
head(res)
