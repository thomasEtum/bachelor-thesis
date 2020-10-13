# create a gene space from genes shared in TFs and drugs
# measure their distance

# needs all_knock_fc from knockTF_process.R and long_cmap from combat_foldchanges.R

library(data.table)
library(dplyr)

cm_genes <- long_cmap$genes %>% unique() # all genes always appear
small_knock_fc <- all_knock_fc[,mean(Log2FC), by=list(TF, Gene)]
colnames(small_knock_fc) <- c("TF", "Gene", "FC")
TF_most_wanted <- small_knock_fc[, .N, by=Gene]
tf_genes <- TF_most_wanted[N==308]$Gene # genes that appear in all 308 TFs in KnockTF

gene_space <- intersect(cm_genes, tf_genes) # 1485 genes
TF_space <- small_knock_fc[Gene %in% gene_space]
CM_space <- long_cmap[genes %in% gene_space]

# distance functions
euclidean_dist <- function(vec1, vec2){
  d <- vec1 - vec2
  di <- d^2
  dis <- sum(di)
  dist <- dis^(1/2)
  return(dist)
}
distant_space <- function(tra_space, drug_space){
  dr <- drug_space$Drug %>% unique()
  fac <- tra_space$TF %>% unique()
  results <- data.table()
  for (drug in dr) {
    for (tf in fac){
      s <- tra_space[TF==tf]
      s <- s[order(Gene)]
      t <- drug_space[Drug==drug]
      t <- t[order(Gene)]
      distance <- euclidean_dist(s$FC, t$FC)
      results <- rbind(results, data.table(tf, drug, distance))
    }
    print(paste("drug", drug, "done"))
  }
  return(results)
}

# prepare for standardization in python using sklearn StandardScaler
pyCM_space <- dcast(CM_space, Drug~Gene, value.var = "FC")
pyTF_space <- dcast(TF_space, TF~Gene, value.var = "FC")
write.table(pyCM_space, file="path/to/cmap/space/values.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")
write.table(pyTF_space, file="path/to/knocktf/space/values.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")
# normalized dataset
pyCM_space <- fread("path/to/cmap/space/normalized/values.tsv")
pyTF_space <- fread("path/to/knocktf/space/normalized/values.tsv")

pyCM_space <- melt(pyCM_space, id.vars = "Drug", variable.name = "Gene", value.name = "FC")
pyTF_space <- melt(pyTF_space, id.vars = "TF", variable.name = "Gene", value.name = "FC")

distances <- distant_space(pyTF_space, pyCM_space)

# distances can be compared to network proximity values
