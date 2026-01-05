import sys, os, re, argparse, pickle

import numpy as np
import scipy as sp
import math

import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr

import statsmodels.api as sm

import random
import itertools
import time

import matplotlib.pyplot as plt

def nb_zero_prob(mu, r):
	return (r / (r + mu)) ** r

def compute_diff_stats(mean_1, mean_2, var_1, var_2, size_1, size_2, r_1, r_2, countnonzero_1, countnonzero_2):
	
	mean_diff = mean_1 - mean_2
	lfc_mean = np.log2(mean_1 / mean_2)
	fano_diff = (var_1 / mean_1) - (var_2 / mean_2)
	lfc_fano = np.log2((var_1 / mean_1) / (var_2 / mean_2))
	
	zero_ratio_1 = (size_1 - countnonzero_1) / size_1
	zero_ratio_2 = (size_2 - countnonzero_2) / size_2
	
	zero_ratio_diff = zero_ratio_1 - zero_ratio_2
	lfc_zero_ratio = np.log2(zero_ratio_1 / zero_ratio_2)
	
	epsilon_1 = zero_ratio_1 - nb_zero_prob(mean_1, r_1)
	epsilon_2 = zero_ratio_2 - nb_zero_prob(mean_2, r_2)
	
	epsilon_diff = epsilon_1 - epsilon_2
	
	return np.array([mean_diff, lfc_mean, fano_diff, lfc_fano, zero_ratio_diff, lfc_zero_ratio, epsilon_diff])

def compute_all_stats(group1, group2, axis=-1):
	
	# calculates sum of 
	sum_1 = np.sum(group1, axis=0, keepdims=True)
	sum_2 = np.sum(group2, axis=0, keepdims=True)
	mean_1 = np.mean(group1, axis)
	mean_2 = np.mean(group2, axis)
	var_1 = np.var(group1, axis)
	var_2 = np.var(group2, axis)
	
	group1_norm = group1 * (np.median(sum_1) / sum_1)
	group2_norm = group2 * (np.median(sum_2) / sum_2)
	
	group1_norm_10k = group1 * (10000 / sum_1)
	group2_norm_10k = group2 * (10000 / sum_2)
	
	mean_1_norm = np.mean(group1_norm, axis)
	mean_2_norm = np.mean(group2_norm, axis)
	var_1_norm = np.var(group1_norm, axis)
	var_2_norm = np.var(group2_norm, axis)
	
	mean_1_norm_10k = np.mean(group1_norm_10k, axis)
	mean_2_norm_10k = np.mean(group2_norm_10k, axis)
	var_1_norm_10k = np.var(group1_norm_10k, axis)
	var_2_norm_10k = np.var(group2_norm_10k, axis)
	
	n_cells_1 = group1.shape[axis]
	n_cells_2 = group2.shape[axis]
	
	n_transcripts_1 = group1.sum().sum()
	n_transcripts_2 = group2.sum().sum()
	
	r_1 = np.sqrt(2) * np.sqrt(n_transcripts_1 / n_cells_1)
	r_2 = np.sqrt(2) * np.sqrt(n_transcripts_2 / n_cells_2)
	
	stats = compute_diff_stats(mean_1, mean_2, var_1, var_2, n_cells_1, n_cells_2, r_1, r_2, np.count_nonzero(group1, axis=axis), np.count_nonzero(group2, axis=axis))
	stats_norm = compute_diff_stats(mean_1_norm, mean_2_norm, var_1_norm, var_2_norm, n_cells_1, n_cells_2, r_1, r_2, np.count_nonzero(group1_norm, axis=axis), np.count_nonzero(group2_norm, axis=axis))
	stats_norm_10k = compute_diff_stats(mean_1_norm_10k, mean_2_norm_10k, var_1_norm_10k, var_2_norm_10k, n_cells_1, n_cells_2, r_1, r_2, np.count_nonzero(group1_norm_10k, axis=axis), np.count_nonzero(group2_norm_10k, axis=axis))
	
	return np.concatenate([stats, stats_norm, stats_norm_10k])

def get_stats(df, suffix=''):
	
	df_stats = pd.DataFrame()
	df_stats.index = df.index
	df_stats[f'mean{suffix}'] = df.mean(axis=1)
	df_stats[f'var{suffix}'] = df.var(axis=1)
	df_stats[f'std{suffix}'] = df.std(axis=1)
	df_stats[f'dispersion{suffix}'] = [df_stats[f'var{suffix}'][i]/df_stats[f'mean{suffix}'][i] if df_stats[f'mean{suffix}'][i] > 0 else np.nan for i in range(len(df_stats[f'var{suffix}']))]
	df_stats[f'coeff_var{suffix}'] = [df_stats[f'std{suffix}'][i]/df_stats[f'mean{suffix}'][i] if df_stats[f'mean{suffix}'][i] > 0 else np.nan for i in range(len(df_stats[f'std{suffix}']))]
	
	for column in df_stats.columns:
		df_stats[f'log_{column}'] = np.log(df_stats[column])
	
	return df_stats

def compare_count_changes(df_counts_1, df_counts_2, counts_name_1, counts_name_2, full_comparison_name, path_out_dir):
		
	df_stats_1 = get_stats(df_counts_1, suffix='')
	df_stats_2 = get_stats(df_counts_2, suffix='')
	
	df_stats_1.sort_index(inplace=True)
	df_stats_2.sort_index(inplace=True)
	
	df_diff = df_stats_2 - df_stats_1
	
	# rename columns with diff
	df_diff.columns = [f'{column}_diff' for column in df_diff.columns]
	
	##########
	# EXPORT DIFF CSV
	##########
	path_out = os.path.join(path_out_dir, 'diff_expanded_data.csv')
	df_diff.to_csv(path_out)
	
	return df_diff
	
def plot_results(df_return, statistic_functions, counts_name_1, normalization, comparison, path_out_dir):
	
	for i_stat, statistic_function in enumerate(statistic_functions):
		
		if 'zero_ratio' in statistic_function:
			df_return[f'valid_{statistic_function}'] = np.where(
				df_return[f'abs_value_{statistic_function}'] > df_return[f'ci_range_{statistic_function}'],
				'red',
				'blue'
			)
	
		path_bootstrap_dir_plot = os.path.join(path_out_dir, normalization)
			
		if 'norm10k' in statistic_function:
			path_bootstrap_dir_plot += '.normalized.gm10000'
		elif 'norm' in statistic_function:
			path_bootstrap_dir_plot += '.normalized.gmauto'
			
		path_bootstrap_dir_plot = os.path.join(path_bootstrap_dir_plot, comparison)
		
		# some values might be inf, so drop these for plotting
		df_clean = df_return[[f'ci_range_{statistic_function}', f'value_{statistic_function}', 'xvalue', f'valid_{statistic_function}']].copy()
		df_clean.replace([np.inf, -np.inf], np.nan, inplace=True)
		df_clean.dropna(inplace=True)
		
		df_clean = df_clean.loc[df_clean['xvalue'] > 0]
		
		ci_ranges_clean = df_clean[f'ci_range_{statistic_function}']
		xvalues_clean = df_clean['xvalue']
		real_values_clean = df_clean[f'value_{statistic_function}']
		colors_clean = df_clean[f'valid_{statistic_function}']
		
		path_out_png = os.path.join(path_bootstrap_dir_plot, f'{statistic_function}.png')
		plt.scatter(np.log(xvalues_clean), real_values_clean, s=1, c=colors_clean)
		plt.plot([min(np.log(xvalues_clean)), max(np.log(xvalues_clean))], [0, 0], color='red')
		plt.title(f'{counts_name_1} Expression Mean vs.\n{comparison} {statistic_function}')
		plt.xlabel(f'ln({counts_name_1} Expression Mean)')
		plt.ylabel(f'Non-bootstrapped {comparison}\n{statistic_function}')
		plt.tight_layout()
		plt.savefig(path_out_png)
		plt.close()
		
		path_out_png = os.path.join(path_bootstrap_dir_plot, f'{statistic_function}_95ci.png')
		plt.scatter(np.log(xvalues_clean), ci_ranges_clean, s=1)
		plt.title(f'{counts_name_1} Expression Mean vs.\n{statistic_function} 95% CI')
		plt.xlabel(f'ln({counts_name_1} Expression Mean)')
		plt.ylabel(f'{comparison}\n{statistic_function} 95% CI')
		plt.tight_layout()
		plt.savefig(path_out_png)
		plt.close()
	
		
def perform_bootstrap(np_counts_1, np_counts_2, statistic_functions, counts_name_1, normalization, comparison, path_out_dir):
	
	xvalues = np.mean(np_counts_1, axis=1)
	all_stats = compute_all_stats(np_counts_2, np_counts_1, axis=1)
		
	path_pkl = os.path.join(path_out_dir, normalization, comparison, f'bootstrap.pkl')
	if not os.path.exists(path_pkl) or True:
		t1 = time.time()
		bootstrap_stats = sp.stats.bootstrap((np_counts_2, np_counts_1), statistic=compute_all_stats, method='basic', n_resamples=10, vectorized=True, axis=1)
		t2 = time.time()
		print(f'bootstrapping took {t2-t1} seconds')
		with open(path_pkl, 'wb') as handle:
			pickle.dump(bootstrap_stats, handle, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		with open(path_pkl, 'rb') as handle:
			bootstrap_stats = pickle.load(handle)
		
	df_return = pd.DataFrame()
	df_return['xvalue'] = xvalues
		
	for i_stat, statistic_function in enumerate(statistic_functions):
			
		real_values = all_stats[i_stat]
		
		ci_lows = bootstrap_stats.confidence_interval[0][i_stat]
		ci_highs = bootstrap_stats.confidence_interval[1][i_stat]
		ci_ranges = ci_highs - ci_lows
		real_values = all_stats[i_stat]
		
		df_return[f'ci_range_{statistic_function}'] = ci_ranges
		df_return[f'value_{statistic_function}'] = real_values
		df_return[f'ci_high_{statistic_function}'] = ci_highs
		df_return[f'ci_low_{statistic_function}'] = ci_lows
		df_return[f'abs_value_{statistic_function}'] = np.abs(df_return[f'value_{statistic_function}'])
		
		df_return[f'valid_{statistic_function}'] = np.where(
			df_return[f'abs_value_{statistic_function}'] > df_return[f'ci_range_{statistic_function}'],
			'red',
			'blue'
		)
	
	return df_return

# Example
#python ./bootstrap_stats.py -path_counts_dir ./transcript_counts -path_out_dir ./bootstrap

parser = argparse.ArgumentParser()

parser.add_argument("-path_counts_dir", required = True)
parser.add_argument("-path_out_dir", required = True)

# get arguments
args = parser.parse_args()

path_counts_dir = args.path_counts_dir
path_out_dir = args.path_out_dir

os.makedirs(path_out_dir, exist_ok=True)

# get files that have been filtered & normalized in the same way
regex_method = '(?<=\.).*(?=\.parquet\.gz)'
input_files = os.listdir(path_counts_dir)
input_files = [file for file in input_files if '.parquet.gz' in file]

normalization_methods = list(set([re.search(regex_method, file).group(0) for file in input_files]))

stat_functions = [	
	'linear_diff_mean', 'lfc_mean', 'linear_diff_fano', 'lfc_fano', 'linear_diff_zero_ratio', 'lfc_zero_ratio', 'epsilon_diff',
	'linear_diff_mean_norm', 'lfc_mean_norm', 'linear_diff_fano_norm', 'lfc_fano_norm', 'linear_diff_zero_ratio_norm', 'lfc_zero_ratio_norm', 'epsilon_diff_norm',
	'linear_diff_mean_norm10k', 'lfc_mean_norm10k', 'linear_diff_fano_norm10k', 'lfc_fano_norm10k', 'linear_diff_zero_ratio_norm10k', 'lfc_zero_ratio_norm10k', 'epsilon_diff_norm10k',
]
	
for normalization_method in normalization_methods:
		
	# this script is a bit wonky... we actually bootstrap counts first and
	# only then normalize counts... so we skip the actually normalized ones
	# we instead look at downsampled, non-downsampled things, etc.
	if 'raw' in normalization_method or 'normalized' in normalization_method:
		continue
	
	path_normalization_dir = os.path.join(path_out_dir, normalization_method)
	
	files = [input_file for input_file in input_files if f'{normalization_method}.parquet.gz' in input_file]
	
	filepaths = [os.path.join(path_counts_dir, file) for file in files]
	filenames = [file.replace(f'.{normalization_method}.parquet.gz', '') for file in files]
	
	file_names_paths = zip(filenames, filepaths)
	
	file_pair_combinations = itertools.combinations(file_names_paths, 2)
	
	for file_pair_combination in file_pair_combinations:
		
		counts_name_1 = file_pair_combination[0][0]
		path_counts_1 = file_pair_combination[0][1]
		
		counts_name_2 = file_pair_combination[1][0]
		path_counts_2 = file_pair_combination[1][1]
		
		names_list = sorted([counts_name_1, counts_name_2], reverse=True)
		paths_list = [x for _, x in sorted(zip([counts_name_1, counts_name_2], [path_counts_1, path_counts_2]), reverse=True)]
		
		counts_name_1 = names_list[0]
		counts_name_2 = names_list[1]
		
		path_counts_1 = paths_list[0]
		path_counts_2 = paths_list[1]
		
		comparison = f'{counts_name_2}-{counts_name_1}'
		path_comparison_dir = os.path.join(path_normalization_dir, comparison)
		os.makedirs(path_comparison_dir, exist_ok=True)
		
		print(f'processing for {comparison} {normalization_method}...')
		
		full_comparison_name = f'{normalization_method} {comparison}'
		
		path_bootstrap_csv = os.path.join(path_comparison_dir, f'{normalization_method}.{comparison}.bootstrapstats.csv')
		
		if not os.path.exists(path_bootstrap_csv) or True:
		
			#######
			# GET STATISTICS
			#######
			
			print('loading counts...')
			df_counts_1 = pd.read_parquet(path_counts_1)
			df_counts_2 = pd.read_parquet(path_counts_2)
			
			df_counts_1 = df_counts_1.T
			df_counts_2 = df_counts_2.T
			
			common_genes = list(set(df_counts_1.index).intersection(df_counts_2.index))
			num_genes = len(common_genes)
			
			df_counts_1.drop([gene for gene in df_counts_1.index if gene not in common_genes], inplace=True)
			df_counts_2.drop([gene for gene in df_counts_2.index if gene not in common_genes], inplace=True)
			
			df_diff = compare_count_changes(df_counts_1, df_counts_2, counts_name_1, counts_name_2, full_comparison_name, path_comparison_dir)
			
			print('bootstrapping')
			
			df_counts_1.sort_index(inplace=True)
			df_counts_2.sort_index(inplace=True)
			
			# DO NOT REMOVE THIS, NEED FOR PROPER INDEXING!!!!
			common_genes = df_counts_1.index
			
			# explicitly round floats resulting from normalization into ints 
			np_counts_1 = df_counts_1.to_numpy()
			np_counts_2 = df_counts_2.to_numpy()
			
			print(np.shape(np_counts_1))
			print(np.shape(np_counts_2))
			
			df_bootstrap_stats = perform_bootstrap(np_counts_1, np_counts_2, stat_functions, counts_name_1, normalization_method, comparison, path_out_dir)
			
			df_bootstrap_stats.index = common_genes
			path_bootstrap_csv = os.path.join(path_comparison_dir, f'{normalization_method}.{comparison}.bootstrapstats.csv')
			df_bootstrap_stats.to_csv(path_bootstrap_csv)
			
		else:
		
			df_bootstrap_stats = pd.read_csv(path_bootstrap_csv)
		
		plot_results(df_bootstrap_stats, stat_functions, counts_name_1, normalization_method, comparison, path_out_dir)