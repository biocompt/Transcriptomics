# Transcriptomics
Scripts to analyze RNA-Seq data. During my work at FPGMX, I collaborate in different projects related to the analysis of RNA-Seq data.

## Libraries
We have used DESeq2 package (https://doi.org/doi:10.18129/B9.bioc.DESeq2), the most popular and widespread R package among the scientific community for the analysis of RNA-Seq data. Furthermore, for the normalization of RNA-Seq read counts between samples we have used RUVSeq (https://doi.org/doi:10.18129/B9.bioc.RUVSeq), a R package that removes the unwanted variation.

```
library(DESeq2)
library(RUVSeq)
```

## The variables
One of the key points of writing a script is the possibility to reutilice it for different analysis without changing many things. Because of that, here we designate different variables that can allow us to work with our data:

```
pheno <- Here we indicate the name of the phenotype that we are analyzing. It is recommendable to set the same name that the file that we want to upload.
counts <- HTSEQ/RSEM. Here we specify the tool that we have used to extract the counts. HTSEQ is to map counts to Genes, while RSEM is to extract transcripts.
```
