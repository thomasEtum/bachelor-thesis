import gseapy
import os
import pandas as pd

gene_sets = "path/to/genesets.gmt"
rnk_dir = "path/to/ranked/foldfiles/directory"
pval_out = "/path/to/pvalue/output/file.tsv"
gsea_outdir = "path/to/full/gsea/output/directory"

rnks = os.listdir(rnk_dir)
df = pd.DataFrame()

for file in rnks:
    full_path = rnk_dir + "/" + file
    out = gsea_outdir + "/" + os.path.splitext(file)[0]
    gs = gseapy.prerank(rnk=full_path, gene_sets=gene_sets, outdir=out, min_size=0, max_size=15000, no_plot=True)
    pval = gs.res2d.pval
    df[file] = pval
    df.to_csv(pval_out, sep="\t")