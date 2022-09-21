library(DESeq2)
colnames(cds) <- rep(c(paste("Control", 1:CONTROLS), paste("Caso", 1:CASES)))

## PCA of object ##
vsd <- vst(ds, blind = F)
jpeg(paste("plots/PCA_", tipo, "_", feno, "_", tecnica, ".jpeg", sep = ""),
     quality = 100)
plotPCA(vsd, intgroup = 'Feno')
dev.off()

## Heatmap ##
library(pheatmap)
res <- res[res$pvalue < 0.05,]
de_genes <- rownames(res)
select <- subset(counts(ds), rownames(ds) %in% de_genes)
select <- order(rowMeans(select), decreasing=T)
df <- as.data.frame(colData(ds)[,c("cov1","Feno")])
jpeg(paste("plots/Heatmap_", tipo, "_", feno, "_", tecnica, ".jpeg", sep = ""),
     quality = 100)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T, annotation_col=df)
dev.off()

## Sample dists ##
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Feno
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
jpeg(paste("plots/Sample_Dist_", tipo, "_", feno, "_", tecnica, ".jpeg", sep = ""),
     quality = 100)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, col=colors)
dev.off()

## Dispersion ##
jpeg(paste("plots/Dispersion_", tipo, "_", feno, "_", tecnica, ".jpeg", sep = ""),
     quality = 100)
plotDispEsts(cds)
dev.off()

## Box-plot ##
jpeg(paste("plots/BoxPlot_", tipo, "_", feno, "_", tecnica, ".jpeg", sep = ""),
     quality = 100)
boxplot(log10(assays(cds)[["cooks"]]), range=0, las=2)
dev.off()
