# comparison to msigdb and omnipathdb

#1 omnipathdb
# TF interactions from omnipath compared to grn edges
library(data.table)
library(dplyr)
library(rlist)
library(OmnipathR)

t <- import_TFregulons_Interactions(confidence_level = c("A", "B", "C", "D"))
dt <- as.data.table(dt)
omniedges <- dt[, c(3, 4)]
omniedges[, forward := paste0(source_genesymbol, "_", target_genesymbol)]
omniedges[, backward := paste0(target_genesymbol, "_", source_genesymbol)]

compare_edges <- function(grntbl, omniedgetbl=omniedges){
  if(length(grntbl)>1){
    grnedges <- grntbl[, c(1, 2)]
    grnedges[, combined:= paste0(regulatoryGene, "_", targetGene)]
    commonedges <- as.data.table(intersect(omniedgetbl$forward, grnedges$combined))
    commonedges <- rbind(commonedges, as.data.table(intersect(omniedgetbl$backward, grnedges$combined)))
    commonedges <- unique(commonedges)
    return(commonedges)
  } else return(data.table())
}
shared_omniedges <- lapply(comp, function(ls){
  retls <- lapply(ls, compare_edges)
})
shared_omni_table <- lapply(shared_omniedges, function(x){return(rbindlist(x, idcol="diseases"))}) %>% rbindlist()

#2 MSigDB
## hubgenes from grns
## input: alldt from grn_pairwise_comparison.R

 path_all_edges_file <- "path/to/fused/edges/file.tsv"
 alldt <- fread(path_all_edges_file)
 srcs <- unique(alldt$src)

## sort genes by number of connections
hub_list <- lapply(srcs, function(s){
   rel <- alldt[src==s]
   rel$weight <- NULL
   rel <- melt(rel, id.vars="src", variable.name="type", value.name="gene")
   rN <- rel[, .N, by=gene]
   rN <- rN[order(-N)]
   return(rN)
 })
names(hub_list) <- srcs

hub_tbl <- lapply(hub_list, function(mat){
   return(mat[1:30])
})
hub_tbl <- rbindlist(hub_tbl, idcol = "src") %>% as.data.table()

## MSIGDB comparison
path_geneset <- "path/to/msigdb/geneset"
gset <- fread(path_geneset)
gset <- gset[-1,]

dt <- hub_tbl[src=="source_icd10"] 

shared_genes <- intersect(gset, dt$gene)






