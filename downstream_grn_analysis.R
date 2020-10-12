# comparison to msigdb and omnipathdb

#1 omnipathdb
# TF interactions from omnipath compared to grn edges
library(data.table)
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