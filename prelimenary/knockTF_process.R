# KnockTF 

# takes expression profiles from knockTF as input and 
#   - calculates fold changes
#   - prepares gene sets for gsea in gmt format
#   - prepares ranked gene lists for gsea in rnk format
# User input:
#   - change input and output paths to correct file locations
#   - change cutoff for gene sets if necessary


# libraries
library(data.table)
library(dplyr)

# input files (download from knockTF)
path_to_TFinfo <- "path/to/TFinfo.csv"
path_to_TFexpr_dir <- "path/to/knockTF/expression/directory"

path_to_all_diff <- "path/to/KnockTF_alldatasets_diffexpr.txt"

# output files
path_to_TFfc_matrix <- "path/to/calculated/fold/changes/file.tsv"
path_to_geneSets <- "path/to/created/genesets.gmt"
path_to_rankedGene_dir <- "path/to/ranking/file/directory"


# all datasets fc (downloaded directly)
all_knock_data <- fread(path_to_all_diff)
all_knock_fc <- all_knock_data[, .SD, .SDcols=c("TF", "Gene", "Log2FC")]
# or calculated

# read expression data
TFinfo <- fread(path_to_TFinfo, sep = ",")
profilist <- list.files(path_to_TFexpr_dir, pattern = ".txt", full.names = F)
profilist <- lapply(profilist, function(x){
  path <- paste0(path_to_TFexpr_dir, x)
  id <- sub(".txt", "", x)
  tbl <- fread(path, sep="\t", header = T)
  tbl$src <- id
  return(tbl)
})

# calculate fold changes from expression values
fc <- function(dt){
  newdt <- data.table()
  newdt$Gene <- dt$X
  newdt$src <- dt$src
  newdt$meank <- dt[,grep("k[[:digit:]]", names(dt)), with=F ] %>% apply( MARGIN = 1, mean)
  newdt$meanc <- dt[,grep("c[[:digit:]]", names(dt)), with=F ] %>% apply( MARGIN = 1, mean)
  newdt$fc <- newdt[, log2(meank/meanc)]
  return(newdt)
}

profilist <- lapply(profilist, as.data.table)

# removed DataSet_01_363 because of formatting issues
profilist <- profilist[-1]

TFfolds <- lapply(profilist, fc) %>% rbindlist()

# get TF name instead of dataset ID
getTFName <- function(dtID){
  return(TFinfo[TFinfo$DatasetID == dtID,]$TF)
}
TFfolds$TF <- lapply(TFfolds$src, getTFName)

# fold change matrix: gene vs TF
TFtable <- dcast(TFfolds, Gene ~ TF, value.var = "fc", fun.aggregate = mean)
write.table(TFtable, file=path_to_TFfc_matrix, sep = "\t", quote = F, row.names = F, col.names = T)


### prepare gmt for gsea

#alltable <- dcast(all_knock_fc, Gene ~ TF, value.var = "Log2FC", fun.aggregate = mean)

# if it needs to be reloaded
#TFtable <- fread(path_to_TFfc_matrix, sep="\t")
# long data format
#TFfolds <- melt(TFtable, id.vars = "Gene", variable.name = "TF", value.name = "FC")

# change cutoff as necessary
cutoff <- 1.5

subTF <- TFfolds[FC < -cutoff | FC > cutoff,]
TFnames <- unique(subTF$TF)
TFlist <- lapply(TFnames, function(x){
  place <- data.frame(subTF[TF==x]$Gene) %>% t() %>% as.data.frame()
  return(cbind(x, "describe_me", place))
})
names(TFlist) <- TFnames
dt <- rbindlist(TFlist, use.names = T, fill = T)
write.table(dt, file = path_to_geneSets,
            quote=F, col.names = F, row.names = F, sep="\t")
# replace all \tNA with "" in notepad
###

### prepare ranked gene lists for gsea

lapply(seq_along(TFtable), function(i, dt, n){
  if(i!=1){
    towrite <- dt[,c(1, ..i)]
    name <- paste0(path_to_rankedGene_dir, n[[i]], ".rnk")
    write.table(towrite, file = name, quote = F, sep = "\t", row.names = F, col.names = F)
  }
}, dt=TFtable, n=colnames(TFtable))


### all TF gene lists

allTFs <- all_knock_fc$TF %>% unique()
lapply(allTFs, function(x, dt){
  if (x!= 1){
    towrite <- dt[TF==x]
    towrite <- towrite[, c(2, 3)]
    name <- paste0(path_to_rankedGene_dir, x, ".rnk")
    write.table(towrite, file = name, quote = F, row.names = F, col.names = F, sep = "\t")
  }
}, dt = all_knock_fc)



