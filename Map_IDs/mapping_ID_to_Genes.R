library(tidyr)
library(BUSpaRse)

genes <- data.frame(ID = ids)
genes <- genes %>% separate(ID, c("id", "identifier", "extra"))

mapping <- tr2g_ensembl("Homo sapiens", ensembl_version = NULL, 
                        other_attrs = NULL)

if (id == "gene"){
  mapping <- mapping %>% separate(gene, c("id", "identifier"))
  mapping <- mapping[,c(2,4,5)]
  mapping <- unique(mapping)
  
  gene_mapped <- merge(genes, mapping, by = "id")
  gene_mapped$ID <- paste(gene_mapped$id, gene_mapped$identifier, sep = ".")
  gene_mapped$ID <- paste(gene_mapped$ID, gene_mapped$extra, sep = "_")

  gene_mapped <- gene_mapped[,c("ID", "gene_name", "chromosome_name")]
  res$ID <- rownames(res)
  res <- merge(res, gene_mapped)
}
if (id == "transcript") {
  mapping <- mapping %>% separate(transcript, c("id", "identifier"))
  mapping <- mapping[,c(1,4,5)]
  mapping <- unique(mapping)
  
  gene_mapped <- merge(genes, mapping, by = "id")
  gene_mapped <- na.omit(gene_mapped)
  gene_mapped$ID <- paste(gene_mapped$id, gene_mapped$identifier, sep = ".")
  gene_mapped$ID <- paste(gene_mapped$ID, gene_mapped$extra, sep = "_")

  gene_mapped <- gene_mapped[,c("ID", "gene_name", "chromosome_name")]
  res$ID <- rownames(res)
  res <- merge(res, gene_mapped)
}

