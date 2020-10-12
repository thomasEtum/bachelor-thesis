# compare grns from ARACNE and GENIE3 graphically
# expects grns in long format (wide format can be converted using the format_grn() function in cellline_comparison.R)
# grns should be generated using the same celllines (mixed or separate)

# uses slurm!

library(ggplot2)
library(data.table)
library(rslurm)

icd <- read.table("path/to/relevant_diseases_icd10.txt")
colnames(icd) <- c("icdCode")

compare_genie3 <- function(icdCode){
  path_aracne <- "path/to/aracne/grn/directory/"
  path_genie <- "path/to/genie3/grn/directory/"
  path_plot <- "path/to/output/plot/directory"
  
	aracne <- fread(paste0(path_aracne, icdCode, "_grn_aracne.tsv"))
	genie <- fread(paste0(path_genie, icdCode, "drugs_linkList.txt"))
	genie$V1 <- NULL
	aracne$src <- NULL
	dt <- merge(aracne, genie, by=c("regulatoryGene", "targetGene"), suffixes=c(".aracne", ".genie3"), all=TRUE)
	rm(aracne)
	rm(genie)
	ggplot(dt[!(weight.aracne==0&weight.genie3==0)&weight.aracne<=1], aes(weight.aracne, weight.genie3))+
	       geom_point(shape=1)+labs(title=paste0("Disease ", icdCode, " ARACNE vs GENIE3"))
	ggsave(paste0(path_plot, icdCode, "_aracneVSgenie3.png"))
}

sjob <- slurm_apply(compare_genie3, icd, slurm_options = list(ntasks=1, "cpus-per-task"=10, share=TRUE),
		    jobname="aracneVSg3", submit=TRUE)

