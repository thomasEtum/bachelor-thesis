# infer GRNs from gene-expression profiles using either aracne or mrnet
# input: tab separated disease-specific gene-expression files (format: genes | profile1 | profile2 | ...)
# all input files (and nothing else) must be in the specified folder

# uses slurm!

library(data.table)
library(rslurm)
library(minet)

# input
path_expr_files <- "/path/to/expression/directory/"
# output
path_grn_files <- "/path/to/GRN/directory/"

expr_files <- list.files(path_expr_files, pattern=".tsv", full.names=F)
expr_files <- as.data.frame(expr_files)
colnames(expr_files) <- c("this_file")

calculate_grn <- function(this_file){
  this_path <- paste0(path_expr_files, this_file)
  tst_mat <- fread(this_path)
  tst_mat <- as.data.frame(tst_mat)
  rownames(tst_mat) <- tst_mat$genes
  tst_mat$genes <- NULL
  tst_mat <- as.data.frame(t(tst_mat))
  mim <- build.mim(tst_mat)
  mim <- mim[rowSums(is.na(mim)) < ncol(mim)-1,]
  mim <- mim[,colSums(is.na(mim)) < nrow(mim)]
  ## change this to switch algorithms
  tst_grn <- aracne(mim)
  #  tst_grn <- mrnet(mim)
  ##
  this_name <- sub("_file_ending.tsv", "", this_file)
  path_this_grn <- paste0(path_grn_files, this_name, "_new_file_ending.tsv")
  write.table(tst_grn, file=path_this_grn, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
}

sjob <- slurm_apply(calculate_grn, expr_files, jobname="grn_inference",
                    add_objects=c("path_expr_files", "path_grn_files"),
                    slurm_options=list(ntasks=1, "cpus-per-task"=4, share=TRUE), submit = TRUE)
