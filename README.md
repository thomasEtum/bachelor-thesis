# Scripts, supplementary results and data for my bachelor's thesis.

! All scripts contain pathnames that have to be changed before execution !

Files:

 BA_thesis_Thomas_Eska.pdf  
  
  the thesis  
 README.md  
 
  you are here  

  preprocessing:  
    scripts to preprocess the CMAP dataset  
    - correct_cmap_instances_file.py  
      corrects mistakes and shorthand in the instance file downloaded from CMAP  
    - subset_cmap_instances.R  
      select the right expression profiles from the data  
    - normalize_expression_profiles.R  
      read the CEL files and apply gcrma  
    - combat_foldchanges.R  
      apply combat batch correction and calculate foldchanges  

  prelimenary:  
    scripts for GSEA, genespace  
    - knockTF_process.R  
      some preprocessing to prepare knockTF downloads for GSEA and genespace  
    - gsea.py  
      applies GSEA  
    - genespace.R  
      calculates distances between drugs and TFs in genespace  
    - network_proximity.py  
      calculates network proximity between drugs and TFs  
    - dict_tf_targets.pickle  
      contains TF target data as a python dictionary  
    - drugname_uniprot_dict.pickle  
      contains associations of drugs and genes as a python dictionary  

GRN:  
    scripts used in GRN analysis  
    - grn_inference  
      use ARACNE or MRNET to infer GRNs from gene expression data  
    - cellline_comparison_grn.R  
      compare ARACNE and MRNET results between celllines  
    - genie3_comparison.R  
      compare results to GENIE3  
    - grn_pairwise_comparison.R  
      compare GRN features between diseases  
    - downstream_grn_analysis.R  
      compare results to MSigDB and OmnipathDB  
    - dict_diseaseicd10_drug_ctd_repodb2020_complete.pickle  
      drug-disease associations as a python dictionary  
    - cellline_comparison_plots  
      supplementary plots from celline_comparison_grn.R  
    - genie3_comparison_plots  
      supplementary plots from genie3_comparison.R  

