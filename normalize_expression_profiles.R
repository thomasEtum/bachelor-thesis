# read and normalize (gcrma) data from Affymetrix GeneChips (CEL files)
# careful: generates large files; might want to break it up

library(data.table)
library(dplyr)
library(affy)
library(gcrma)

# input
path_to_subinstances <- "path/to/subset/instances/file.tsv"
path_to_CEL_dir <- "path/to/CELfile/directory/"
# output
path_to_raw_expr <- "path/to/normalized/expression/output/file.tsv"

instances <- fread(path_to_subinstances)
CELnames <- instances$scan_id %>% unique()
CELnames <- paste0(path_to_CEL_dir, CELnames, ".CEL")

rawexpr <- ReadAffy(filenames=CELnames)
normed <- gcrma(rawexpr)
exprval <- exprs(normed)
write.table(exprval, file = path_to_raw_expr,
            quote = F, row.names = T, col.names = NA, sep="\t")
