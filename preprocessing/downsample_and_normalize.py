import sys, os, argparse

import numpy as np
from tqdm import tqdm
import seaborn as sns
import pandas as pd

import random

import scipy.sparse as sp_sparse
from scipy.stats import median_abs_deviation

import anndata as ad
import scanpy as sc

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial','Roboto']
rcParams['text.usetex'] = False
rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt

############
# This script performs different normalizations and downsamplings of 
# transcript counts for testing and downstream processing

# Input: design csv file, includes multiple h5 files from cellranger
# Output: normalized transcript counts (csv, h5)

# Pseudocode:
# 1. Downsample cells in datasets to ensure same number of cells
# 2. Downsample transcripts to ensure same distribution of transcript counts per cell
# 3. Normalize per library size
# 4. Plot and save output
############

# Usage:

#python preprocessing/downsample_and_normalize.py -path_input_csv preprocessing/input_tko.csv -path_out_dir transcript_counts/wt_tko
#python preprocessing/downsample_and_normalize.py -path_input_csv preprocessing/input_dko.csv -path_out_dir transcript_counts/wt_dko

def sc_load(path: str, sample: str):

	if os.path.isdir(path):

		if os.path.exists(f'{path}/matrix.mtx.gz'):

			 # this is basically only for the weird data... it is actually a large dataset, kinda nice to keep it gzipped
		
			adata = sc.read_mtx(f'{path}/matrix.mtx.gz')
			adata= adata.T
			
			adata_bc=pd.read_csv(f'{path}/barcodes.tsv.gz',header=None,sep='\t')
			adata_features=pd.read_csv(f'{path}/features.tsv.gz',header=None,sep='\t')
			
			adata.obs.index = adata_bc[0].tolist()
			adata.var.index = adata_features[0].tolist()
	
			adata.var['gene_ids'] = ''
			adata.var['feature_types'] = 'Gene Expression'
			adata.var['genome'] = 'hg38'
			
			path_met = f'{path}/cell_metadata.txt.gz'
			df_metadata = pd.read_csv(path_met, sep='\t')
			adata.obs['Reference_Annotation'] = df_metadata['Reference_Annotation']
			adata.obs['Seurat_Clusters_res.0.1'] = df_metadata['Seurat_Clusters_res.0.1']

		else:

			adata = sc.read_10x_mtx(path)
			adata.obs['Reference_Annotation'] = sample
			adata.obs['Seurat_Clusters_res.0.1'] = sample
			
	elif path.endswith('.h5ad'):
		adata = sc.read_h5ad(path)
		adata.obs['Reference_Annotation'] = sample
		adata.obs['Seurat_Clusters_res.0.1'] = sample
	else:
		adata = sc.read_10x_h5(path)
		adata.obs['Reference_Annotation'] = sample
		adata.obs['Seurat_Clusters_res.0.1'] = sample

	print(adata.obs)
	print(adata.var)
	
	adata.var_names_make_unique(join='_dup_')
	
	return adata
	
def sc_filter(adata, sample, path_out_dir):

	# Generate qc metrics
	adata.var["mt"] = adata.var_names.str.startswith("MT-")

	sc.pp.calculate_qc_metrics(
		adata, 
		qc_vars=["mt"], 
		inplace=True, percent_top=[20], 
		log1p=True
		)
		
	sc.pl.violin(
		adata,
		["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes", "pct_counts_mt"],
		jitter=0.4,
		multi_panel=True,
	)
	plt.savefig(os.path.join(path_out_dir, f'{sample}.qc.png'))
	plt.close()
	
	# Identify outliers per qc metric
	adata.obs["counts_outlier"] = is_outlier(adata, "log1p_total_counts", 5)
	adata.obs["genes_outlier"] = is_outlier(adata, "log1p_n_genes_by_counts", 5)
	adata.obs["topgenes_outlier"] = is_outlier(adata, "pct_counts_in_top_20_genes", 5)
	adata.obs["mito_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] <= 1)

	# Identify outliers based on ALL qc metrics 
	adata.obs["outlier"] = (
		adata.obs["counts_outlier"]  
		| adata.obs["genes_outlier"] 
		| adata.obs["topgenes_outlier"] 
		| adata.obs["mito_outlier"]
	)
	# Remove outliers
	adata = adata[~adata.obs["outlier"]].copy()
	
	# Remove non-developmental genes
	mt_genes = adata.var_names[adata.var_names.str.startswith('MT-')]
	rpl_genes = adata.var_names[adata.var_names.str.startswith('RPL')]
	rps_genes = adata.var_names[adata.var_names.str.startswith('RPS')]
	hist_genes = adata.var_names[adata.var_names.str.startswith('HIST')]
	rna_genes = adata.var_names[adata.var_names.str.startswith('RNA')]
	
	remove_genes = list(mt_genes) + list(rpl_genes) + list(rps_genes) + list(hist_genes) + list(rna_genes)
	adata = adata[:, ~adata.var_names.isin(remove_genes)]

	# Perform gene filtering
	sc.pp.filter_genes(adata, min_cells=3)

	# doublets
	sc.pp.scrublet(adata)
	
	adata.obs.drop(columns=['predicted_doublet'], inplace=True)

	return adata

# our outlier detection function from before
def is_outlier(adata, metric: str, nmads: int):
	M = adata.obs[metric]
	outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
			np.median(M) + nmads * median_abs_deviation(M) < M
	)
	
	print(f'{metric}: {np.median(M) - nmads * median_abs_deviation(M)} - {np.median(M) + nmads * median_abs_deviation(M)}')
	
	return outlier

def plot_sequencing_depth(df, filter_type, name, path_out_dir):
	########
	# PLOT DIFFERENT GENES vs. SUM TO CONFIRM NORMALIZATION SUCCESS
	########
	path_png = os.path.join(path_out_dir, f'sequencing_depth.{filter_type}.{name}.png')
	sums = df.sum(axis=1)
	
	genes_to_plot = ['GAPDH', 'TPM3', 'SPCS1', 'SESN3', 'RAD9A']
	
	for gene in genes_to_plot:
		plt.scatter(np.log(sums), np.log(df[gene] + 1), label=gene, s=2)
		b, a = np.polyfit(np.log(sums), np.log(df[gene] + 1), deg=1)
		xseq = np.linspace(min(np.log(sums)), max(np.log(sums)), num=100)
		plt.plot(xseq, a + b * xseq, color="k", lw=2.5)
		
	plt.xlabel('log(Sequencing Depth) (Sum counts in cell)')
	plt.ylabel('log(Gene Count in Cell + 1)')
	plt.legend()
	plt.title(f'{name} Sequencing Depth {filter_type}')
	plt.savefig(path_png, dpi=200)
	plt.close()

def qc_plots(dfs, names, filter_type, path_out_dir):

	path_plot = os.path.join(path_out_dir, filter_type)
	os.makedirs(path_plot, exist_ok=True)
	
	for df, name in zip(dfs, names):
		plot_sequencing_depth(df, filter_type, name, path_plot)
	
	sum_data = []
	for df, name in zip(dfs, names):
		cell_sums = df.sum(axis=1)
		temp_df = pd.DataFrame({
			'sum': cell_sums,
			'sample': name
		})
		sum_data.append(temp_df)
	
	dfs_concat = pd.concat(sum_data, ignore_index=True)
	
	sns.kdeplot(data=dfs_concat, x='sum', hue='sample', common_norm=False)
	
	plt.xlabel('Sum Transcripts in Cell')
	plt.ylabel('Density')
	plt.title(f'Sum of Transcripts in Cell {filter_type}')
	
	path_png = os.path.join(path_plot, f'transcripts_per_cell.{filter_type}.png')
	plt.savefig(path_png, dpi=200)
	plt.close()

def normalize(df, grand_median=None):
	# normalize library size as per the first part of https://kb.10xgenomics.com/hc/en-us/articles/115004583806-How-are-the-UMI-counts-normalized-before-PCA-and-differential-expression-
	
	sums = df.sum(axis=1)
	if grand_median is None:
		grand_median = np.median(sums)
		print(f'grand median: {grand_median}')
	
	scaling_factors = grand_median / sums
	
	df_normalized = df.mul(scaling_factors, axis=0)
	
	return df_normalized, grand_median

def normalize_and_save(df, filter_type, name, path_out_dir):
	
	df_normalized, grand_median = normalize(df)
	df_normalized_10k, grand_median_10k = normalize(df, grand_median=10000)

	os.makedirs(os.path.join(path_out_dir, f'{filter_type}.normalized.gmauto'), exist_ok=True)
	os.makedirs(os.path.join(path_out_dir, f'{filter_type}.normalized.gm10000'), exist_ok=True)
	os.makedirs(os.path.join(path_out_dir, f'{filter_type}'), exist_ok=True)
	
	path_csv = os.path.join(path_out_dir, f'{filter_type}.normalized.gmauto', f'{name}.{filter_type}.normalized.gmauto.parquet.gz')
	df_normalized.to_parquet(path_csv, compression='gzip')
	
	path_csv = os.path.join(path_out_dir, f'{filter_type}.normalized.gm10000', f'{name}.{filter_type}.normalized.gm10000.parquet.gz')
	df_normalized_10k.to_parquet(path_csv, compression='gzip')
	
	path_csv = os.path.join(path_out_dir, f'{filter_type}', f'{name}.{filter_type}.parquet.gz')
	df.to_parquet(path_csv, compression='gzip')
	
	adata = ad.AnnData(df)
	
	path_h5ad = os.path.join(path_out_dir, f'{filter_type}', f'{name}.{filter_type}.h5ad')
	adata.write_h5ad(path_h5ad, compression='gzip')	
	
def downsample_transcripts(downsampled_dfs, names, min_cell_sums):

	dfs_return = []
	
	for i_df, df in enumerate(downsampled_dfs):
		outfile = names[i_df]
		sums = df._sums.values
		
		data = df.values.astype(np.int32)
		n_cells, n_genes = data.shape
		
		# Find cells that need downsampling
		needs_downsampling = sums > min_cell_sums
		cells_to_process = np.where(needs_downsampling)[0]
		
		print(f'Processing {len(cells_to_process)} cells that need downsampling...')
		
		for cell_idx in tqdm(cells_to_process):
			counts = data[cell_idx]
			target_sum = min_cell_sums[cell_idx]
			current_sum = sums[cell_idx]
			
			# Skip if somehow already correct
			if current_sum <= target_sum:
				continue
				
			# Create gene index array more efficiently
			# Use cumsum to create boundaries for each gene
			cumsum = np.cumsum(counts)
			n_to_remove = int(current_sum - target_sum)
			
			# Random selection of transcript positions
			remove_positions = np.sort(np.random.choice(int(current_sum), n_to_remove, replace=False))
			
			# Use searchsorted to find which gene each position belongs to
			gene_assignments = np.searchsorted(cumsum, remove_positions, side='right')
			
			# Count removals per gene using bincount
			remove_counts = np.bincount(gene_assignments, minlength=n_genes)
			
			# Apply removals
			data[cell_idx] -= remove_counts.astype(np.int32)
		
		# Create result dataframe
		result_df = pd.DataFrame(data=data, index=df.index, columns=df.columns)
		dfs_return.append(result_df)

	return dfs_return

parser = argparse.ArgumentParser()

parser.add_argument("-path_input_csv", required = True)
parser.add_argument("-path_out_dir", required = True)

# get arguments
args = parser.parse_args()
path_input_csv = args.path_input_csv
path_out_dir = args.path_out_dir

df_input = pd.read_csv(path_input_csv)

paths = list(df_input['path'])
names = list(df_input['name'])

#####
# RAW
#####
adatas = [sc_load(path, name) for path, name in zip(paths, names)]

for adata, name in zip(adatas, names):
	normalize_and_save(adata.to_df(), 'raw', name, path_out_dir)
    
qc_plots([adata.to_df() for adata in adatas], names, 'raw', path_out_dir)
	
###########
# FILTERED
###########
adatas = [sc_filter(adata, name, path_out_dir) for adata, name in zip(adatas, names)]

# do separate filtered plotting
for adata, name in zip(adatas, names):
	normalize_and_save(adata.to_df(), 'filtered', name, path_out_dir)

qc_plots([adata.to_df() for adata in adatas], names, 'filtered', path_out_dir)
	
################
# SAME CELL/GENE - OPTIMIZED
################
print('downsampling cells...')
count_dfs = [adata.to_df() for adata in adatas]
min_cell_count = min(len(df) for df in count_dfs)

# Sample cells more efficiently
downsampled_dfs = [df.sample(n=min_cell_count, random_state=42) for df in count_dfs]

print('getting common genes...')
# More efficient intersection using sets
common_genes = set(downsampled_dfs[0].columns)
for df in downsampled_dfs[1:]:  # Skip first df since we already have its columns
	common_genes.intersection_update(df.columns)
common_genes = list(common_genes)

# Filter to common genes in-place
for i, df in enumerate(downsampled_dfs):
	downsampled_dfs[i] = df[common_genes]

# Filtered plotting
for df, name in zip(downsampled_dfs, names):
	normalize_and_save(df, 'filtered.samecg', name, path_out_dir)

qc_plots(downsampled_dfs, names, 'filtered.samecg', path_out_dir)

##############
# DOWNSAMPLING TRANSCRIPTS
##############
print('calculating sums and sorting...')
for i, df in enumerate(downsampled_dfs):
	df_sums = df.sum(axis=1)
	# Sort both dataframe and sums together
	sort_idx = df_sums.argsort()
	downsampled_dfs[i] = df.iloc[sort_idx].copy()
	df_sums = df_sums.iloc[sort_idx]
	# Store sums as attribute to avoid recalculation
	downsampled_dfs[i]._sums = df_sums

# Calculate minimum sums for each cell position more efficiently
min_cell_sums = np.minimum.reduce([df._sums.values for df in downsampled_dfs])

print('downsampling transcripts...')
downsampled_dfs = downsample_transcripts(downsampled_dfs, names, min_cell_sums)

for df, name in zip(downsampled_dfs, names):
	normalize_and_save(df, 'filtered.samecg.downsampled', name, path_out_dir)

qc_plots(downsampled_dfs, names, 'filtered.samecg.downsampled', path_out_dir)
