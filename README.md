# Transcriptomics
Scripts to analyze RNA-Seq data. During my work at FPGMX, I collaborate in different projects related to the analysis of RNA-Seq data.

## Libraries
We have used DESeq2 package (https://doi.org/doi:10.18129/B9.bioc.DESeq2), the most popular and widespread R package among the scientific community for the analysis of RNA-Seq data. Furthermore, for the normalization of RNA-Seq read counts between samples we have used RUVSeq (https://doi.org/doi:10.18129/B9.bioc.RUVSeq), a R package that removes the unwanted variation.

```
library(DESeq2)
library(RUVSeq)
```

## Data
One of the key points of writing a script is the possibility to reutilice it for different analysis without changing many things. Because of that, here we designate different variables that can allow us to work with our data:

```
pheno <- read.table() of the pheno data.
counts <- read.table() of the count data.
```
One of the requirements of DESeq2 is to order the sample data in the same way in both files (pheno and counts). Sometimes, both files have different order so we have to re-order it.

```
samples <- colnames(counts)
pheno <- subset(pheno, rownames(pheno) %in% samples)
counts <- counts[, rownames(pheno)]
all(rownames(pheno) == colnames(counts))
```

## Normalization of data
It is recommendable to remove unwanted variation (batch, library preparation, and other nuisance) to improve our analysis.
```
RUVgApply <- function(DESEQcounts, controls, cases)
batchRemoval <- RUVgApply(counts, CONTROLS, CASES)  
set2 <- batchRemoval[[3]]
```

## DESeq object
1. **countData**. We specified the normalized counts after applying RUV with the command `counts()`.
2. **colData**. We specified the dataframe with the phenotype and the covariables of our study. It is recommendable to factorize those that can be factorize before create the DESeq object, although they will be factorize if we did not do it.
3. **design**. We specify the columns to run the differential expression analysis, putting the one that we want to analyze as the last one. In the RUV normalization we set the *empirical (k)* as 3, so we have to add to the design formula the 3 factors that accounts for the effects.
```
dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = covs,
                              design = ~ W_1 + W_2 + W_3 + covs + Pheno)
dds$Pheno <- relevel(dds$Pheno, ref = "Control")
```

## DE analysis
Once the DESeq object is created, we can filtered by removing all the genes/transcripts that have no counts or that not passed a minimun number of counts. These filters should be adapted to our data. These filter will allow us to increase the speed of our analysis.
```
ds <- estimateSizeFactors(dds)
keep <- rowSums(counts(ds, normalized=T) >= nlect) >= npat
ds <- ds[keep,]
```
Once we have filtered the low counts, we can run the analysis.
```
cds <- DESeq(ds, parallel = T)
```
We can extract the results and store it as dataframe. 
```
res <- na.omit(data.frame(results(cds, contrast = c("Pheno", "Case", "Control"), alpha = 0.05)))
res <- res[order(res$padj),]
```

## Plot the results

### Counts plot
DESeq has inserted a function 



padj<-read.delim("sinCOV_caseVScontrol.tsv",sep = "\t",header = TRUE, na.strings = "NA", dec = ".")
rownames(padj)<-padj$X
EnhancedVolcano(padj,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano plot - CaCo',
                ylab = bquote(~-Log[10] ''~italic(P)~'-value'),
                xlab = bquote(~Log[2] italic('fold change')),
                pointSize = 1.0,
                gridlines.minor = FALSE,
                col = c('gray', 'skyblue', 'pink1', 'purple'),
                pCutoff = 5.5e-2)
rm(list=ls())
