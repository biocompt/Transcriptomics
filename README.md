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

```
dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=T)), decreasing=TRUE)[1:500]
vsd <- vst(dds, blind=FALSE)
groups <- as.data.frame(colData(dds)[,c(cov, Pheno)])
```

### Counts plot
DESeq has a function that allowed us to plot counts in each of the groups that we want to compare, to see if there are a particular distribution in each one.
```
plotCounts(cds, gene = "Gene_name", intgroup = "Pheno")
```

### PCA
PCA can show differences between the groups that we want to compare.
```
plotPCA(vsd, intgroup="Pheno")
```
![Counts' PCA](https://i.stack.imgur.com/txw4G.png)
### Volcano plot


### Heatmap
```
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F,
         cluster_cols=T, annotation_col=groups)
```
![Heatmap_DESEq2](https://bioinformatics.stackexchange.com/questions/13976/rnaseq-biological-replicates-not-clustering-in-pca-plots)
