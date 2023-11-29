if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("rnaseqGene")
library(tidyverse)
library(DESeq2)
library(htmltools)
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
dds <- DESeqDataSetFromMatrix(countData=countData,
colData=metaData,
design=~dex, tidy = TRUE)
dds
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds)
res
head(results(dds, tidy=TRUE))
res <- res[order(res$padj),]
head(res)
plotCounts
par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000179094", intgroup="dex")
plotCounts(dds, gene="ENSG00000148175", intgroup="dex")
plotCounts(dds, gene="ENSG00000152583", intgroup="dex")
plotCounts(dds, gene="ENSG00000120129", intgroup="dex")
plotCounts(dds, gene="ENSG00000189221", intgroup="dex")
plotCounts(dds, gene="ENSG00000116584", intgroup="dex")
savehistory("~/gene expression assembly/rhistory.r")
