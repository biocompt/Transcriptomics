# Transcriptomics
Scripts to analyze RNA-Seq data. During my work at FPGMX, I collaborate in different projects related to the analysis of RNA-Seq data.

## Libraries
We have used DESeq2 package (https://doi.org/doi:10.18129/B9.bioc.DESeq2), the most popular and widespread R package among the scientific community for the analysis of RNA-Seq data. Furthermore, for the normalization of RNA-Seq read counts between samples we have used RUVSeq (https://doi.org/doi:10.18129/B9.bioc.RUVSeq), a R package that removes the unwanted variation.

```
library(DESeq2)
library(RUVSeq)
```
