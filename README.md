# README for code samples

This repository contains examples of code I generated to analyze data from different projects.

Descriptions of the code snippets are below:

Visualization of genomic interaction datasets:
- TMEM161B_interactions.Rmd (code with the actual visualization exploration)
(https://oldmantae.github.io/Code_Sample/TMEM161B_Genome_Interaction_Visualizations/TMEM161B_interactions.html)
- Genome_Viewer_Functions_hg38.R (functions used for creating visualizations)
- Interaction_Analysis_Functions.R (functions for formatting different datatypes and databases to play nicely)
- File_Conversions.Rmd (file conversions and liftover for different reference genomes)

A simple pipeline for pre-processing raw data from a molecular tool we developed (capture-based Massively parallel reporter assay):
- caMPRA_tag_seq.umi_tools.10nt.Snakefile.yaml
- EXP_00007_Process_Tag_Seq_umi_tools_10nt_snakemake.sh

Script for pre-processing targeted sequencing data based on custom designed molecular inversion probes:
- MIPs_non_MIPgen_v8.sh

More lengthy analysis of random mutagenesis + capture-based Massively parallel reporter assay (.Rmd files as well as html output for graphics):
- EXP_00007_HAR_caMPRA_Tag_Seq_Analysis_Mutagenesis_Analysis_1_QC_and_filtering.Rmd
(https://oldmantae.github.io/Code_Sample/random_mut_caMPRA_analysis/EXP_00007_HAR_caMPRA_Tag_Seq_Analysis_Mutagenesis_Analysis_1_QC_and_filtering.html)
- EXP_00007_HAR_caMPRA_Tag_Seq_Analysis_Mutagenesis_Analysis_2_Prelim_variant_effects.Rmd
(https://oldmantae.github.io/Code_Sample/random_mut_caMPRA_analysis/EXP_00007_HAR_caMPRA_Tag_Seq_Analysis_Mutagenesis_Analysis_2_Prelim_variant_effects.html)
- EXP_00007_HAR_caMPRA_Tag_Seq_Analysis_Mutagenesis_Analysis_3_Masked_variant_effects.Rmd
(https://oldmantae.github.io/Code_Sample/random_mut_caMPRA_analysis/EXP_00007_HAR_caMPRA_Tag_Seq_Analysis_Mutagenesis_Analysis_3_Masked_variant_effects.html)

On-going analysis of variant effects of patient variants (small snippet):
- EXP_00010_Variant_TF_Analysis_Sample.Rmd

Design of MPRA constructs for patient variant validation (designing controls as well as python script/jupyter notebook to design the actual probes):
- EXP_00027_SCZ_Variant_Validation.Rmd
- FASTA_extraction_functions.R
- EXP_00027_SCZ_pwmscan_check.sh (database wide TFBS motif scanning - for speed - necessary for checking changing RE sites without impacting potential TF binding sites)
- EXP_00027_SCZ_Variant_Validation.ipynb (for designing the full oligo to be printed (including check and removal of sequences that contain restriction enzymes needed for experiment))
