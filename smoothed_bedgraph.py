import os
import sys
import argparse
import random
import re

import pybedtools

import tempfile
import uuid

import time

from tqdm import tqdm
from joblib import Parallel, delayed, dump, load, parallel_backend

import statsmodels.api as sm
from fractions import Fraction

import numpy as np
import pandas as pd

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	'''
	alist.sort(key=natural_keys) sorts in human order
	http://nedbatchelder.com/blog/200712/human_sorting.html
	(See Toothy's implementation in the comments)
	'''
	return [ atoi(c) for c in re.split(r'(\d+)', text) ]
	
def bed_file(filename):
	ext = [".bed"]
	return filename.endswith(tuple(ext))

def get_df(path, category_field_name):

	if os.path.isdir(path):
	
		filenames = os.listdir(path)
		filenames = list(filter(bed_file, filenames))
		filenames.sort(key=natural_keys)
		
		filepaths = [os.path.join(path, filename) for filename in filenames]
		
		num_categories = len(filepaths)
		
		dfs = []
		
		strand_column = 4
		group_column = 5
		
		for i, filepath in enumerate(filepaths):
			
			cols_to_use = [0, 1, 2, int(strand_column), int(group_column)]
			df = pd.read_csv(filepath, sep='\t', header=None, usecols=cols_to_use, comment='#')[cols_to_use]
			
			df.rename(columns={int(strand_column):3, int(group_column):5}, inplace=True)
			
			df[4] = filenames[i].replace('.bed', '')
			
			df = df[[0,1,2,3,4,5]]
			
			# ensure consistent naming convention
			df[0] = [f'chr{i}' if 'chr' not in i else i for i in list(df[0])]
			
			dfs.append(df)
			
		df_merged = pd.concat(dfs)
		
		return df_merged, num_categories
		
	else:
		
		# we could evaluate which field num the column name is in
		df = pd.read_csv(path, sep='\t')
		columns = list(df.columns)
		category_field_num = columns.index(category_field_name)
		group_column = columns.index('gene')
		strand_column = columns.index('strand')
		
		if category_field_num is not None:
			
			cols_to_use = [0, 1, 2, strand_column, category_field_num, group_column]
			df = pd.read_csv(path, sep='\t', header=None, usecols=cols_to_use, comment='#')[cols_to_use]
			
			df.rename(columns={int(strand_column):3, int(category_field_num):4, int(group_column):5}, inplace=True)
			
			df = df[[0,1,2,3,4,5]]
			
			num_categories = len(list(set(list(df[4]))))
			
			# ensure consistent naming convention
			df[0] = [f'chr{i}' if 'chr' not in i else i for i in list(df[0])]
			
			return df, num_categories
			
		else:
			
			print(f'Please provide a field number for identifying categories in {path}')
			sys.exit()

def contiguous_bedgraph_file(filename):
	ext = [".cbg"]
	return filename.endswith(tuple(ext))

def bedgraph_file(filename):
	ext = [".bedGraph", ".bedgraph", ".bg", ".bedgraph.gz", ".bedGraph.gz", ".bg.gz"]
	return filename.endswith(tuple(ext))
	
def make_contiguous(path):
	
	last_end = 0
	last_chrom = ''
	newpath = path + '.cbg'
	
	with open(path, 'r') as infile, open(newpath, 'w') as outfile:
		for i_line, lineread in enumerate(infile):
		
			line = lineread.split('\t')
			
			chrom = line[0]
			bedgraph_start = int(line[1])
			bedgraph_end = int(line[2])
			
			if 'chr' not in line[0]:
				line[0] = f'chr{line[0]}'
				chrom =  f'chr{chrom}'
				
			lineread = '\t'.join(line)
			
			if chrom == last_chrom and bedgraph_start != last_end:				
				outfile.write('\t'.join([chrom, str(last_end), str(bedgraph_start), '0.0']) + '\n')
				
			outfile.write(lineread)
				
			last_end = bedgraph_end
			last_chrom = chrom
			
	return newpath
	
def get_smoothed_values(np_positions, np_xvals, np_values, num_lowess_bins, lowess_float, bedgraph, region):
	
	smoothed_cpg_values = sm.nonparametric.lowess(np_values, np_positions, xvals=np_xvals, frac=lowess_float, it=0)
	smoothed_cpg_values_list = [[bedgraph, region, xval, smoothed_cpg_values[i]] for i, xval in enumerate(np_xvals)]
		
	return smoothed_cpg_values_list

def compute_lowess(i_slice, slice, num_lowess_bins, lowess_float, metadata, value_counters, refpoint):
	
	np_positions = slice[0:value_counters[i_slice], 0]
	np_values = slice[0:value_counters[i_slice], 1]
	
	region = metadata[i_slice]['region']
	bedgraph = metadata[i_slice]['bedgraph']
	
	# get bounds for lowess
	min_bound = min(np_positions)
	max_bound = max(np_positions)
	np_xvals = np.linspace(min_bound, max_bound, num_lowess_bins)
	
	lowess_data_list = get_smoothed_values(np_positions, np_xvals, np_values, num_lowess_bins, lowess_float, bedgraph, region)
	
	df = pd.DataFrame.from_records(lowess_data_list, columns=['bedgraph', 'region', 'position', 'value'])
	
	return df

def intersect_regions(file_data, path_temp_pybedtools):
	
	# have to set this again for some reason? weird
	pybedtools.helpers.set_tempdir(path_temp_pybedtools)
	
	file = file_data['file']
	region = file_data['region']
	
	path_feature = file_data['full_path']
	feature_name = file.replace('.bed', '')
	
	filename_bedgraph = file_data['filename_bedgraph'].replace('.cbg', '')
	path_bedgraph = file_data['path_bedgraph']
	
	bedtool_feature = pybedtools.BedTool(path_feature)
	num_fields = bedtool_feature.field_count()
	
	bedtool_bedgraph = pybedtools.BedTool(path_bedgraph)
	
	# have to move since pybedtools tries to delete
	moveto_path = os.path.join(path_temp_pybedtools, f'{region}.{filename_bedgraph}.bed')
	
	intersection = bedtool_feature.intersect(bedtool_bedgraph, wa=True, wb=True).moveto(moveto_path)
	
	num_total_columns = intersection.field_count()
	
	return [intersection, region, filename_bedgraph, num_fields, num_total_columns]
	
def get_relative_position(df, pos_col):
	# Limit pos to region limits
	pos = df[pos_col].clip(lower=df[1], upper=df[2])
	
	# Calculate relative position based on strand
	relpos = np.where(
		df[3] == '+',
		pos - df['refpos'],      # positive strand
		df['refpos'] - pos        # negative strand
	)
	
	return relpos

def downsample_intersection(path_intersection, i_bedgraph, bedgraph_file, regions, num_downsampled_values, refpoint, reflimit, group_dict, output):
	
	# pseudocode:
	# 1. each intersection is a bedgraph intersected with all regions
	# 2. downsample values
	#		using all possible values for lowess is slow and unnecessary so downsample
	# 		downsampling after loading all values is memory intensive, so values are populated by probability
	#		downsampling occurs at bp-resolution value so does not bias for long or short regions
	# 3. fill memmap
	
	# fields:
	# 0:chrom, 1:region_start, 2:region_end, 3:region_strand, 4:category, 5:gene, 6:chrom, 7:bedgraph_start, 8:bedgraph_end, 9:bedgraph_value
	
	df = pd.read_csv(path_intersection, header=None, sep='\t')
	
	print(df)
	
	count_dict = {bedgraph_file: {region:0 for region in regions}}
	
	# for now we just use midpoint
	if refpoint == 'midpoint':
		df['refpos'] = (df[2] + df[1])/2
	elif refpoint == 'tss':
		df['refpos'] = np.where(df[3] == '+', df[1], df[2])
	else:
		print('please select a valid refpoint, ("midpoint" or "tss")')
		sys.exit()
		
	df['relative_start'] = get_relative_position(df, 7)
	df['relative_end'] = get_relative_position(df, 8)
	
	start = df['relative_start']
	end = df['relative_end']
	
	# for negative strands, relative_end ends up being smaller so we must take the minimum
	df['relative_start'] = np.minimum(start, end)
	df['relative_end'] = np.maximum(start, end)
	
	df['length'] = df['relative_end'] - df['relative_start']
	
	# remove rows where the start of the intersection is greater than our limit
	df = df[df['relative_start'] >= -reflimit]
	df = df[df['relative_end'] <= reflimit]
	
	print(df)
	
	# now we have relative positions of bedgraph starts and ends
	# need to get lengths of values contained
	for i_region, region in enumerate(regions):
		
		df_region = df.loc[df[4] == region]
		
		count = 0
		
		for record in df_region.to_dict('records'):
			
			bedgraph_value = record[9]
			gene_length = np.minimum(record[2] - record[1], reflimit)
			
			rel_start = int(record['relative_start'])
			rel_end = int(record['relative_end'])
			
			for position in range(rel_start, rel_end):
				
				if count < num_downsampled_values:
					
					count_dict[bedgraph_file][region] += 1
					
					if refpoint == 'tss':
						output[i_region][i_bedgraph][count, 0] = position / gene_length
					else:
						output[i_region][i_bedgraph][count, 0] = position
						
					output[i_region][i_bedgraph][count, 1] = bedgraph_value
					
				else:
					
					# replace with decreasing probability
						
					j = random.randint(0, count - 1)
					
					if j < num_downsampled_values:
						
						if refpoint == 'tss':
							output[i_region][i_bedgraph][j, 0] = position / gene_length
						else:
							output[i_region][i_bedgraph][j, 0] = position
							
						output[i_region][i_bedgraph][j, 1] = bedgraph_value
						
				count += 1
	
	count_dict[bedgraph_file][region] = count
	
	# we can delete intersection file here
	os.remove(path_intersection)
	del df
	
	return count_dict


# Example
#python ./smoothed_bedgraph.py -path_bedgraph ./GSM3618718 -path_regions ./regions.bed -field_name "Combined TKO ZR Change Class" -path_out_csv ./smooth.csv -lowess_frac 1/8 -num_lowess_bins 250 -num_downsampled_values 1000000 -refpoint midpoint -reflimit 50000

# construct arguments
parser = argparse.ArgumentParser()

parser.add_argument("-path_bedgraph", required = True)
parser.add_argument("-path_regions", required = True)
parser.add_argument("-lowess_frac", required = True)
parser.add_argument("-num_lowess_bins", required = True)
parser.add_argument("-num_downsampled_values", required = True)
parser.add_argument("-refpoint", required = True)
parser.add_argument("-reflimit", required = True)
parser.add_argument("-path_out_csv", required = True)
parser.add_argument("-field_name", required = False)
parser.add_argument("--na_as_zero", action='store_true')

# get arguments
args = parser.parse_args()

path_bedgraph = args.path_bedgraph
path_regions = args.path_regions
lowess_frac = args.lowess_frac
num_lowess_bins = int(args.num_lowess_bins)
num_downsampled_values = int(args.num_downsampled_values)
refpoint = args.refpoint
reflimit = int(args.reflimit)
path_out_csv = args.path_out_csv
field_name = args.field_name
na_as_zero = args.na_as_zero

lowess_float = float(0) + Fraction(lowess_frac)

##################
# DIRECTORIES
##################
temp_dir = tempfile.gettempdir()

##################
# LOAD DATA
##################
path_features = os.path.join(temp_dir, f'{str(uuid.uuid4())}.bed')
df_regions, num_categories = get_df(path_regions, field_name)
df_regions = df_regions[df_regions[1] < df_regions[2]]
df_regions.to_csv(path_features, index=False, header=False, sep='\t')

group_dict = {group:i for i, group in enumerate(list(df_regions[5]))}

print(df_regions)

print('Counts of bed regions:')
regions = list(set(list(df_regions[4])))
counts = df_regions[4].value_counts()
print(counts)

feature_bedtool = pybedtools.BedTool(path_features)

all_files = os.listdir(path_bedgraph)
bedgraph_files = list(filter(bedgraph_file, all_files))
bedgraph_files.sort(key=natural_keys)
bedgraph_paths = [os.path.join(path_bedgraph, file) for file in bedgraph_files]
	
##################
# MAKE CONTIGUOUS IF NA_AS_ZERO
##################
if na_as_zero:
	print('Making bedgraphs contiguous...')
	with parallel_backend("loky", inner_max_num_threads=2):
		bedgraph_paths = Parallel(n_jobs=len(bedgraph_paths), pre_dispatch="all", verbose=100)(delayed(make_contiguous)(path) for path in bedgraph_paths)
	
##################
# INTERSECT IN SERIES
##################
intersection_paths = []
for path_bedgraph in bedgraph_paths:
	
	filename = os.path.basename(path_bedgraph)
	bedgraph_bedtool = pybedtools.BedTool(path_bedgraph)
	
	print(f'Intersecting {filename}...')
	t1 = time.time()
	path_bed = os.path.join(temp_dir, f'{str(uuid.uuid4())}.bed')
	intersection = feature_bedtool.intersect(bedgraph_bedtool, wa=True, wb=True)
	intersection.saveas(path_bed)
	intersection_paths.append(path_bed)
	t2 = time.time()
	print(f'Intersection took {t2-t1} seconds')
	
##################
# PROCESS AND DOWNSAMPLE
##################
print(f'Downsampling and processing...')
# output holds all data
# shape = num_regions x num_bedgraphs x num_downsampled_values x (position & value)
# this is gated by how quickly we can write to memory... need to find max processes that make sense. more = inefficient disk space usage
# also gated by how much disk space we have
path_output_memmap = os.path.join(temp_dir, f'{str(uuid.uuid4())}.memmap')
output = np.memmap(path_output_memmap, dtype=float, shape=(num_categories, len(bedgraph_paths), num_downsampled_values, 2), mode='w+')
with parallel_backend("loky", inner_max_num_threads=2):
	count_dicts = Parallel(n_jobs=min(24, len(bedgraph_files)), pre_dispatch="all", verbose=100)(delayed(downsample_intersection)(path, i_bedgraph, bedgraph_files[i_bedgraph], regions, num_downsampled_values, refpoint, reflimit, group_dict, output) for i_bedgraph, path in enumerate(intersection_paths))
	
count_dicts = {k: v for d in count_dicts for k, v in d.items()}
df_counts = pd.DataFrame(count_dicts)
print('Number of basepairs in each category:')
print(df_counts)

##################
# COMPUTE LOWESS
##################
print(f'Computing lowess...')
t1 = time.time()
# shape = region x bedgraph x parameter x level x num_downsampled_cpgs x 3
slices = []
metadata = []
counters = []
for i_region, region in enumerate(regions):
	for i_bedgraph, bedgraph in enumerate(bedgraph_files):
		slices.append(output[i_region, i_bedgraph])
		metadata.append({'region': region, 'bedgraph': bedgraph})
		counters.append(count_dicts[bedgraph][region])

# can parallelize infinitely
num_processes = 90
with parallel_backend("loky", inner_max_num_threads=2):
	frames = Parallel(n_jobs=num_processes, pre_dispatch="all")(delayed(compute_lowess)(i_sl, sl, num_lowess_bins, lowess_float, metadata, counters, refpoint) for i_sl, sl in enumerate(slices))
t2 = time.time()
print(f'Lowess took {t2-t1} seconds...')

df_main = pd.concat(frames)
df_main = df_main.reset_index()
df_main.to_csv(path_out_csv)

print(df_main)

##################
# CLEANUP
##################
if os.path.exists(path_output_memmap):
	os.remove(path_output_memmap)
if os.path.exists(path_features):
	os.remove(path_features)
for intersection_path in intersection_paths:
	if os.path.exists(intersection_path):
		os.remove(intersection_path)

	