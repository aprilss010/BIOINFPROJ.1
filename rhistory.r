#first download the bioconductor package , this allows us to run DESEQ2 for RNA-SEQ data downloaded
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("rnaseqGene")
#these are the packages I downloaded 
library(tidyverse)
library(DESeq2)
library(htmltools)
#Extraction of Downloaded RNA-SEQ data
zipF<- "C:\\path\\to\\my\\zipfile\\https://GSE52778_All_Sample_FPKM_Matrix.txt.gz"
outDir<-"C:\\Users\\Name\\Documents\\unzipfolder"
unzip(zipF,exdir=outDir)
countData <- read.csv(airway_scaledcounts.csv,header=TRUE,sep=",")
countData <- read.csv('airway_scaledcounts.csv',header=TRUE,sep=",")
head(countData)
`metadata<- read.csv('airway matrix metadata.csv', header = TRUE,sep = ",")
metadata<- read.csv('airway_matrix_metadata.csv', header = TRUE,sep = ",")
metaData<- read.csv('airway_matrix_metadata.csv', header = TRUE,sep = ",")
metaData
#creation of DESEQ2 environment using dds command
dds <- DESeqDataSetFromMatrix(countData=countData,
colData=metaData,
design=~dex, tidy = TRUE)
dds
dds <- DESeq(dds)
# commands used to reveal results from DESEQ2 
res <- results(dds)
res
head(results(dds, tidy=TRUE))
res <- res[order(res$padj),]
head(res)
#Command used to create plot Map for the 6 genes
plotCounts
par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000179094", intgroup="dex")
plotCounts(dds, gene="ENSG00000148175", intgroup="dex")
plotCounts(dds, gene="ENSG00000152583", intgroup="dex")
plotCounts(dds, gene="ENSG00000120129", intgroup="dex")
plotCounts(dds, gene="ENSG00000189221", intgroup="dex")
plotCounts(dds, gene="ENSG00000116584", intgroup="dex")
#Command used to create the Log-fold change Map 
plotMA(res)
#final product links
link to gene dot plot :https://ibb.co/h2pVHSy
link to log fold count map:https://ibb.co/jDnC2kR
link to p<0.05 results:https://ibb.co/bswZK8f
link to p-value and log fold count assembly:https://ibb.co/MCkNx27
