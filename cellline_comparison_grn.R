# compare GRNs from ARACNE and MRNET in celllines (graphically)
# if necessary use format_grn() to convert a grn from wide to long format
# merges grns to one file for comparison
# for first usage uncomment lines 34-68 (replace paths accordingly)

# uses slurm!

library(ggplot2)
library(ggpubr)
library(data.table)
library(rslurm)


pref <- read.table("path/to/relevant_diseases_icd10.txt")
colnames(pref) <- c("icdCode")

format_grn <- function(grnfile){
  grnmat <- fread(grnfile)
if(ncol(grnmat)>5){
  rown <- grnmat[[1]]
  grnmat <- grnmat[,-1]
  grnmat[lower.tri(grnmat)] <- NA
  grnmat$regulatoryGene <- rown
  grnmat <- melt(grnmat, id.vars = "regulatoryGene", variable.name = "targetGene", value.name = "weight")
  grnmat <- grnmat[!is.na(weight)]
  write.table(grnmat, grnfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")}
  return(grnmat)
}

analize_grn <- function(icdCode){
  path_comp <- "/path/to/merged/grns/directory/"
  path_plot <- "/path/to/plot/directory/" 
  
#  path_grn_p1 <- "/path/to/aracne/pc3/grn/directory/"
#  path_grn_p2 <- "/path/to/mrnet/pc3/grn/directory/"
#  path_grn_m1 <- "/path/to/aracne/mcf7/grn/directory/"
#  path_grn_m2 <- "/path/to/mrnet/mcf7/grn/directory/"
#  path_grn_h1 <- "/path/to/aracne/hl60/grn/directory/"
#  path_grn_h2 <- "/path/to/mrnet/hl60/grn/directory/"
#  suf_p1 <- "_grn_pc3_aracne.tsv"
#  suf_p2 <- "_grn_pc3_mrnet.tsv"
#  suf_m1 <- "_grn_mcf7_aracne.tsv"
#  suf_m2 <- "_grn_mcf7_mrnet.tsv"
#  suf_h1 <- "_grn_hl60_aracne.tsv"
#  suf_h2 <- "_grn_hl60_mrnet.tsv"

#  p1 <- format_grn(paste0(path_grn_p1, icdCode, suf_p1))
#  p2 <- format_grn(paste0(path_grn_p2, icdCode, suf_p2))
#  p <- merge(p1, p2, by=c("regulatoryGene", "targetGene"), suffixes = c(".aracne", ".mrnet"))
#  p1 <- NULL
#  p2 <- NULL
#  m1 <- format_grn(paste0(path_grn_m1, icdCode, suf_m1))
#  m2 <- format_grn(paste0(path_grn_m2, icdCode, suf_m2))
#  m <- merge(m1, m2, by=c("regulatoryGene", "targetGene"), suffixes = c(".aracne", ".mrnet"))
#  m1 <- NULL
#  m2 <- NULL
#  h1 <- format_grn(paste0(path_grn_h1, icdCode, suf_h1))
#  h2 <- format_grn(paste0(path_grn_h2, icdCode, suf_h2))
#  h <- merge(h1, h2, by=c("regulatoryGene", "targetGene"), suffixes = c(".aracne", ".mrnet"))
#  h1 <- NULL
#  h2 <- NULL    
#  mph = Reduce(function(...) merge(..., all = TRUE, by=c("regulatoryGene", "targetGene"), suffixes=c(".mcf7", ".pc3", ".hl60")), list(m, p, h))
#  colnames(mph) <- c("regulatoryGene", "targetGene", "weight.aracne.mcf7", "weight.mrnet.mcf7", "weight.aracne.pc3", "weight.mrnet.pc3", "weight.aracne.hl60", "weight.mrnet.hl60")
#  m <- NULL
#  p <- NULL
#  h <- NULL
#  write.table(mph, paste0(path_comp, icdCode, "_all_grns.tsv"),
#              quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t"

  mph <- fread(paste0(path_comp, icdCode, "_all_grns.tsv"))

 # aracne
  p1 <- ggplot(mph[!(weight.aracne.mcf7==0&weight.aracne.pc3==0)&weight.aracne.mcf7<=1&weight.aracne.pc3<=1],
               aes(weight.aracne.mcf7, weight.aracne.pc3))+geom_point(shape=1)+
               labs(title = paste0(icdCode, "   MCF7 vs PC3"))
  
  p2 <- ggplot(mph[!(weight.aracne.mcf7==0&weight.aracne.hl60==0)&weight.aracne.mcf7<=1&weight.aracne.hl60<=1],
         aes(weight.aracne.mcf7, weight.aracne.hl60))+geom_point(shape=1)+
         labs(title = paste0("MCF7 vs HL60"))
  
  p3 <- ggplot(mph[!(weight.aracne.pc3==0&weight.aracne.hl60==0)&weight.aracne.pc3<=1&weight.aracne.hl60<=1],
         aes(weight.aracne.pc3, weight.aracne.hl60))+geom_point(shape=1)+
         labs(title = paste0("PC3 vs HL60"))
  
  ls1 <- list(p1, p2, p3)
  ggexport(plotlist=ls1, filename=paste0(path_plot, icdCode, "_aracne.png"), ncol=3, nrow=1, width=900, height=300)
  
  # mrnet
  p4 <- ggplot(mph[!(weight.mrnet.mcf7==0&weight.mrnet.pc3==0)&weight.mrnet.mcf7<2&weight.mrnet.pc3<2],
         aes(weight.mrnet.mcf7, weight.mrnet.pc3))+geom_point(shape=1)+
         labs(title = paste0(icdCode, "   MCF7 vs PC3"))
  
  p5 <- ggplot(mph[!(weight.mrnet.mcf7==0&weight.mrnet.hl60==0)&weight.mrnet.mcf7<2&weight.mrnet.hl60<2],
         aes(weight.mrnet.mcf7, weight.mrnet.hl60))+geom_point(shape=1)+
         labs(title = paste0("MCF7 vs HL60"))
  
  p6 <- ggplot(mph[!(weight.mrnet.pc3==0&weight.mrnet.hl60==0)&weight.mrnet.pc3<2&weight.mrnet.hl60<2],
         aes(weight.mrnet.pc3, weight.mrnet.hl60))+geom_point(shape=1)+
         labs(title = paste0("PC3 vs HL60"))

  ls2 <- list(p4, p5, p6)
  ggexport(plotlist=ls2, filename=paste0(path_plot, icdCode, "_mrnet.png"), ncol=3, nrow=1, width=900, height=300)
  
  # compare methods
  p7 <- ggplot(mph[!(weight.aracne.mcf7==0&weight.mrnet.mcf7==0)&weight.aracne.mcf7<=1&weight.mrnet.mcf7<2],
         aes(weight.aracne.mcf7, weight.mrnet.mcf7))+geom_point(shape=1)+
    labs(title = paste0(icdCode, "   MCF7"))
  
  p8 <- ggplot(mph[!(weight.aracne.pc3==0&weight.mrnet.pc3==0)&weight.aracne.pc3<=1&weight.mrnet.pc3<2],
         aes(weight.aracne.pc3, weight.mrnet.pc3))+geom_point(shape=1)+
    labs(title = paste0("PC3"))
 
  p9 <- ggplot(mph[!(weight.aracne.hl60==0&weight.mrnet.hl60==0)&weight.aracne.hl60<=1&weight.mrnet.hl60<2],
         aes(weight.aracne.hl60, weight.mrnet.hl60))+geom_point(shape=1)+
    labs(title = paste0("HL60"))

  ls3 <- list(p7, p8, p9)
  ggexport(plotlist=ls3, filename=paste0(path_plot, icdCode, "_aracnevsmrnet.png"), ncol=3, nrow=1, width =900, height=300)
	   

## all plots:
# pltls <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9)
#  ggexport(plotlist = pltls, filename = paste0(path_plot, icdCode, "_all_plots.png"), ncol=3, nrow=3) 
}

sjob <- slurm_apply(analize_grn, pref, jobname = "cmap_grn_plot",
                    add_objects = c("format_grn"), slurm_options = list(ntasks=1, "cpus-per-task"=10, share=TRUE),
                    submit = TRUE)
