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
