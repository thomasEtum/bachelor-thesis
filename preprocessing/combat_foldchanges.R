# 1. correct the data for batch effects 
# 2. change affy probes to gene symbols
# 3. calculate fold canges

library(data.table)
library(dplyr)
library(sva)
library(rlist)
library(hthgu133a.db)

# input
path_to_norm_expr <- "path/to/normalized/expression/profiles.tsv"
path_to_subset_instances <- "path/to/subset/instances.tsv"
# output
path_to_combat_expr <- "path/to/combat/corrected/expression/values.tsv"
path_to_cmap_fc <- "path/to/cmap/foldchanges.tsv"

raw_expr <- fread(path_to_norm_expr)
colnames(raw_expr)[colnames(raw_expr)=="V1"] <- "probes"
colnames(raw_expr) <- sub(".CEL", "", colnames(raw_expr))
#write.table(raw_expr, path_to_norm_expr, quote = F, col.names = T, row.names = F, sep="\t")
instances <- fread(path_to_subset_instances)

# create lookup tables
batch_lookup <- instances[, .SD, .SDcols=c(2, 15, 14)]
batch_lookup <- unique(batch_lookup)
drug_lookup <- instances[, .SD, .SDcols=c(3, 5, 14, 15)] %>% unique()
drug_lookup$drug <- paste(drug_lookup$cmap_name, drug_lookup$`concentration (M)`, sep = "_")
drug_lookup$cmap_name <- NULL
drug_lookup$`concentration (M)` <-NULL
setcolorder(drug_lookup, neworder = c(3, 2, 1))

# batch correction
# prepare for combat
batch_vec <- sapply(colnames(raw_expr), function(x, lookup){
  return(lookup[scan_id==x,]$batch_id)
}, lookup = batch_lookup) %>% unlist()
batch_lookup[type=="perturbation_scan_id", type:= "1"]
batch_lookup[type=="vehicle_scan_id4", type:= "0"]
batch_lookup$type <- as.numeric(batch_lookup$type)
batch_lookup <- batch_lookup[match(colnames(raw_expr[, 2:1871]), scan_id),]
# combat
mod <- model.matrix(~as.factor(type), data = batch_lookup)
cor_expr <- ComBat(as.matrix(raw_expr[, 2:1871]), batch_vec, mod=mod)
# corrected expression table
cor_expr <- as.data.table(cor_expr)
cor_expr$probes <- raw_expr$probes
setcolorder(cor_expr, c("probes", setdiff(names(cor_expr), "probes")))
# write.table(cor_expr, path_to_combat_expr, quote = F, col.names = T, row.names = F, sep="\t")

# fold changes
druglist <- drug_lookup$drug %>% unique()
cmap_fc <- lapply(druglist, function(drugname, lookup, expr){
  perturbcel <- lookup[drug==drugname & type=="perturbation_scan_id"]$scan_id
  vehiclecel <- lookup[drug==drugname & type=="vehicle_scan_id4"]$scan_id
  means <- data.table()
  means$perturb <- expr[, ..perturbcel,] %>% apply(MARGIN = 1, mean)
  means$vehicle <- expr[, ..vehiclecel,] %>% apply(MARGIN = 1, mean)
  means$fc <- means[, log2(perturb/vehicle)]
  return(means$fc)
}, lookup=drug_lookup, expr=cor_expr)
cmap_fc <- list.cbind(cmap_fc)
cmap_fc <- as.data.table(cmap_fc)
colnames(cmap_fc) <- druglist
cmap_fc$probes <- cor_expr$probes
setcolorder(cmap_fc, c("probes", setdiff(names(cmap_fc), "probes")))

# change affy ids to gene symbols (hgnc)
hgnc_lookup <- AnnotationDbi::select(hthgu133a.db, keys = cmap_fc$probes, columns = c("SYMBOL"))
hgnc_lookup <- as.data.table(hgnc_lookup)

cmap_fc$genes <- sapply(cmap_fc$probes, function(probe, lookup){
  return(lookup[PROBEID==probe,]$SYMBOL[1])
}, lookup=hgnc_lookup)
setcolorder(cmap_fc, c("genes", setdiff(names(cmap_fc), "genes")))
cmap_fc$probes <- NULL

# drop rows with NA, average duplicate genes
cmap_fc <- cmap_fc[!is.na(genes)]
cmap_fc <- cmap_fc[, lapply(.SD, mean, na.rm=F), by=genes]


# write.table(cmap_fc, file=path_to_cmap_fc, quote = F, col.names = T, row.names = F, sep = "\t")
# long_cmap <- melt(cmap_fc, id.vars = "genes", variable.name = "drug", value.name = "fc")
