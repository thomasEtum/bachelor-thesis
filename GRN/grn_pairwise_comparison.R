#1
## get edges above a certain weight from all grns in a directory
## expects tab separated grn (format: gene1 | gene2 | weight)
## output: gene1 | gene2 | weight | src


library(data.table)

path_grn_dir <- "path/to/grn/directory/"
weight_threshold <- 0.1
path_all_edges_file <- "path/to/fused/edges/file.tsv"

lsf <- list.files(path_grn_dir, pattern = ".tsv")
alldt <- data.table()
for(i in lsf){
  dt <- fread(paste0(path_grn_dir, i))
  dt <- dt[weight>0.1]
  name <- sub("_file_ending.tsv", "", i)
  dt$src <- rep(name, nrow(dt))
  alldt <- rbind(alldt, dt)
}
write.table(alldt, file=path_all_edges_file, sep="\t",
            quote=FALSE, row.names = FALSE, col.names = TRUE)

#2
## get edges common between different grns
## input == alldt from #1

library(data.table)
library(dplyr)
library(ggplot2)

# alldt <- fread(path_all_edges_file)
allt$weight <- NULL
srcs <- unique(alldt$src)

# common edges list
comp <- lapply(srcs, function(s){
  dt1 <- dt[src==s]
  dt1$src <- NULL
  setkey(dt1, regulatoryGene, targetGene)
  x <- lapply(srcs, function(s2){
    if (s == s2)return(NULL)
    dt2 <- dt[src==s2]
    dt2$src <- NULL
    setkey(dt2, regulatoryGene, targetGene)
    newdt <- merge(dt1, dt2, all=FALSE)
    return(newdt)
  }dt=alldt)
  return(x)
},dt=alldt)

# identifiers icd10_icd10
y <- list()
for (i in 1:33){
  v <- c()
  for (j in 1:33){
    x <- paste(srcs[i], srcs[j], sep="_")
    comp[[i]][[j]]$src <- x
    v <- append(v, x)
  }
  y[[i]] <- v
}
for (i in 1:33){
  names(comp[[i]]) <- y[[i]]
}

# save it
path_common_edges <- "path/to/common/edge/files/directory/"
for (i in 1:33){
  n <- names(comp[[i]])
  for (j in 1:33){
    fn <- paste0(path_common_edges, "common_edges_", n[[j]], ".tsv")
    write.table(comp[[i]][[j]], fn, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
  }
}

# number of shared edges per disease pair
path_shared_num <- "path/to/number/of/shared/edges/file.tsv"
shared_num <- data.table()
for (i in 1:33){
  for (j in 1:33){
    if(i != j){
      x <- comp[[i]][[j]]
      shared_num <- rbind(shared_num, data.table(srcs[i], srcs[j], nrow(x)))
    }
  }
}
colnames(shared_num) <- c("disease_1", "disease_2", "num_shared_edges")
write.table(shared_num, path_shared_num,
            quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

## plot it
ggplot(shared_num, aes(disease_1, disease_2, fill=num_shared_edges))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_continuous("shared TFs",type = "viridis")

