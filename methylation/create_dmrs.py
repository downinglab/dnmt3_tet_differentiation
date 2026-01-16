import os
import sys
import pandas as pd
from tqdm import tqdm
import argparse

# HUES64 WT and HUES64 DKO
# python /home/data/Shared/shared_datasets/wgbs/nlaszik/code/create_dmrs.py -path_bed /home/data/Shared/shared_datasets/wgbs/nlaszik/hues64wt_dko_wgbs_intersect/GSM1112840.HUES64.WGBS_GSM3662265.HUES64_DKO.WGBS_difference.bed -value_thresholds=-1,-0.8,0.8,1 -distance_threshold 1000 -path_out_dir /home/data/Shared/shared_datasets/wgbs/nlaszik/hues64wt_dko_wgbs_intersect/dmr_0.8 -suffix _dko -min_num_values 3

# HUES8 WT and HUES8 TKO
# python /home/data/Shared/shared_datasets/wgbs/nlaszik/code/create_dmrs.py -path_bed /home/data/Shared/shared_datasets/wgbs/nlaszik/hues8wt_tko_wgbs_intersect/GSM3618720_HUES8_TKO.GSM3618718_HUES8_WT.difference.bedGraph -value_thresholds="-1,-0.3,0.3,1" -distance_threshold 1000 -path_out_dir /home/data/Shared/shared_datasets/wgbs/nlaszik/hues8wt_tko_wgbs_intersect/dmr_0.3_len1000 -suffix _tko -min_num_values 3

# Neural DMRs
# python /home/data/Shared/shared_datasets/wgbs/nlaszik/code/create_dmrs.py -path_bed /home/data/Shared/shared_datasets/wgbs/data/GSE90553/difference/D14_R1_D0_R2_difference.bed -value_thresholds="-1,-0.3,0.3,1" -distance_threshold 1000 -path_out_dir /home/data/Shared/shared_datasets/wgbs/data/GSE90553/popcanyon/dmr_SRR8611939_D14_0.3_len1000 -suffix _nscd14 -min_num_values 5

def get_clusters(df, value_threshold_pair, distance_threshold, min_num_values, suffix=None):
	
	##########
	# 1. Filter differentially methylated CpGs
	##########
	df_filt = df.loc[(df['value'] >= value_threshold_pair[0]) & (df['value'] <= value_threshold_pair[1])]
	
	##########
	# 2. Find DMRs by clustering adjacent CpGs
	##########
	clustered_regions = cluster_regions(df_filt, distance_threshold)
	
	# Convert DMRs to DataFrame
	df_clustered_regions = pd.DataFrame(clustered_regions, columns=['chromosome', 'start', 'end', 'total_value', 'num_values'])
	df_clustered_regions['category'] = f'{value_threshold_pair[0]}_{value_threshold_pair[1]}{suffix}'
	
	##########
	# 3. Filter by number of CpGs in DMR
	##########
	df_clustered_regions = df_clustered_regions.loc[df_clustered_regions['num_values'] > min_num_values]
	return df_clustered_regions

def cluster_regions(df, distance_threshold):
	
	regions = []
	current_region = None
	previous_end = None
	previous_chrom = None
	value_counter = 0
	
	df_records = df.to_records('dict')
	
	previous_chrom = df_records[0]['chromosome']
	
	for row in tqdm(df_records):
		# Check if the current CpG is within the distance cutoff from the previous one
		
		this_chrom = row['chromosome']
		
		if current_region is None or (previous_end is not None and row['start'] - previous_end > distance_threshold) or this_chrom != previous_chrom:
			if current_region is not None:
				current_region.append(value_counter)
				regions.append(current_region)
			current_region = [row['chromosome'], row['start'], row['end'], row['value']]
			value_counter = 1
		else:
			current_region[2] = row['end']  # Extend the current DMR
			current_region[3] += row['value']  # Accumulate methylation diff
			value_counter += 1
			
		previous_end = row['end']
		previous_chrom = row['chromosome']
	
	# Add the last DMR if it exists
	if current_region is not None:
		current_region.append(value_counter)
		regions.append(current_region)
		
	return regions
	
	
def parse_arguments():
	"""
	Parses command-line arguments.
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-path_bed", required = True)
	parser.add_argument("-value_thresholds", required = True)
	parser.add_argument("-distance_threshold", required = True)
	parser.add_argument("-path_out_dir", required = True)
	parser.add_argument("-suffix", required = False)
	parser.add_argument("-min_num_values", required = False, default=1)
	
	# get arguments
	return parser.parse_args()

def main():
	
	"""
	Main function for command-line usage.
	"""
	args = parse_arguments()
	
	path_bed = args.path_bed
	value_thresholds = args.value_thresholds
	distance_threshold = int(args.distance_threshold)
	path_out_dir = args.path_out_dir
	suffix = args.suffix
	min_num_values = int(args.min_num_values)
	
	if suffix is None:
		suffix=''
	
	os.makedirs(path_out_dir, exist_ok=True)
	
	value_thresholds = [float(val) for val in value_thresholds.split(',')]
	value_threshold_pairs = [[value_thresholds[i], value_thresholds[i+1]] for i in range(len(value_thresholds) - 1)]
	
	print(value_threshold_pairs)
	
	df = pd.read_csv(path_bed, sep='\t', header=None, names=['chromosome', 'start', 'end', 'value'], usecols=[0,1,2,3])
	
	cluster_dfs = [get_clusters(df, value_threshold_pair, distance_threshold, min_num_values, suffix) for value_threshold_pair in value_threshold_pairs]
	
	for i, dmr_df in enumerate(cluster_dfs):
		
		path_out = os.path.join(path_out_dir, f'{value_threshold_pairs[i][0]}_{value_threshold_pairs[i][1]}{suffix}.bed')
		dmr_df.to_csv(path_out, index=False, header=False, sep='\t')

if __name__ == "__main__":
	
	main()






