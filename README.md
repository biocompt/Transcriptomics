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
To run the differential expression analysis, we need to create the DESeq object where we specify what our counts are, the phenotype data, and the design of our analysis. Due to normalization, we indicate in countData the normalized counts of the set2 object. In colData, we have to specify the data frame with our phenotype and covariates. It is advisable to factor those that can be factored. Finally, in the design we specify the effect factors (W_1, W_2 and W_3 in our case because our k = 3).

1. Prueba 1
2. Prueba 2
3. Prueba 3

```
dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = covs,
                              design = ~ W1 + W2 + W3 + covs + Pheno)
dds$Pheno <- relevel(dds$Pheno, ref = "Control")
```
