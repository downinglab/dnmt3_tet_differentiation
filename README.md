# General

This repository contains Jupyter Notebooks and scripts for processing data and generating figures for "DNMT3s and TETs adjust CpG methylation canyon width to regulate gene expression and cell fate."

Some processed data is provided, such as BASiCS output and UMRs generated with PopCanyon.

# Dependencies

1. scRNA-seq and methylation data must be downloaded from GEO using get_data.sh
2. Python dependencies must be installed with requirements.txt. Development was done in Python 3.12.
3. Jupyter must be installed
4. For gene pair correlations, a GPU with CUDA is required.
5. For BASiCS processing, R 4.4.1 was used. Requirements are not provided, but processed data is stored in the basics directory

# Workflow

## scRNA-seq Preprocessing

(Figure 1)

Download data, filter/downsample counts, generates BED files for plotting

1. get_data.sh (Download data)
2. preprocessing/downsample_and_normalize.py (Generates multiple different normalizatin h5ads and parquet files)
3. preprocessing/bootstrap_stats.py (Bootstraps multiple statistics)
4. basics/basics_pipeline_wt_tko_prefilt.R (Optional) (BASiCS processing)
5. basics/BASiCS.ipynb (Plotting BASiCS output and creating BED file with statistics)
6. sc_rna_sq/Statistics.ipynb (Combines ad hoc statistics and BASiCS statistics and outputs BED file)

## Gene Pair Correlations

(Figure 1)

1. gene_pair_correlations/find_interactions (executable, requires CUDA-enabled GPU, outputs CSV of Pearson or Spearman correlation of all gene pairs in a parquet.gz file)
2. gene_pair_correlations/Gene Pair Correlation.ipynb

## Methylation Processing

(Figures 2 & 3)

UMRs generated via PopCanyon are already created

1. methylation/create_cpg_differences.py (Creates differential methylation BED files)
2. methylation/create_dmrs.py (Creates DMRs)
3. methylation/Log Odds Heatmaps.ipynb (Creates log odds intersection heatmaps of UMRs & DMRs)
4. methylation/Methylation Canyon Analysis.ipynb (Creates UMR plots)

## Lowess Plotting

(Figures 2 & 3)

1. plotting/lowess_bedgraph.py
2. plotting/Lowess Plotting.ipynb

## UMAP Plotting

(Figures 4 & 5)

1. sc_rna_seq/UMAP.ipynb (Integrates datasets and creates combined UMAP plots)
2. sc_rna_seq/scVelo.ipynb (Requires running velocyto on TKO data to generate loom files, separate virtual environment required)

## smFISH Plotting

(Figure 4)

1. fish/Merge FISH.ipynb (Convenience GUI for merging smFISH TIFF files)