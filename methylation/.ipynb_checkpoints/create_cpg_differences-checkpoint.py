import pandas as pd
import pybedtools

#python ./create_cpg_differences.py

def load_meth(path, precision, skiprows):
    
    df = pd.read_csv(path, usecols=[0,1,2,3,4], names=['chr', 'start', 'stop', 'methylation', 'num_cpgs'], sep='\t', skiprows=1)
    df.drop(columns=['num_cpgs'], inplace=True)
    df['methylation'] = df['methylation'] * precision
    df[['methylation']] = df[['methylation']].astype(int)
    df['start'] = df['start'] + 1

    return df

def subtract_meth(df_1, df_2, precision, path_out):
	
	df_difference = df_1 - df_2
	df_difference.reset_index(inplace=True)
	df_difference.dropna(inplace=True)
	df_difference['methylation'] = df_difference['methylation'] / precision
	
	print(df_difference)
	
	df_difference.to_csv(path_out, index=False, header=False, sep='\t')

# converting to int for subtraction prevents floating point errors
precision = 1e6

df_cpg_wt_hues64 = load_meth('/home/data/Shared/shared_datasets/wgbs/data/GSE126958/GSM4458668/GSM4458668_WGBS_HUES64_WT_ESCs.bed', precision, 1)
print(df_cpg_wt_hues64)

df_cpg_wt_hues8 = load_meth('/home/data/Shared/shared_datasets/wgbs/data/GSE126958/GSM3618718/GSM3618718_HUES8_WT_WGBS.bed', precision, 0)
print(df_cpg_wt_hues8)

df_cpg_dko = load_meth('/home/data/Shared/shared_datasets/wgbs/data/GSE126958/GSM3662265/GSM3662265_HUES64_DKO_P28_WGBS.bed', precision, 0)
print(df_cpg_dko)
	
df_cpg_tko = load_meth('/home/data/Shared/shared_datasets/wgbs/data/GSE126958/GSM3618720/GSM3618720_HUES8_TKO_WGBS.bed', precision, 0)
print(df_cpg_tko)

# set common index
df_cpg_wt_hues64.set_index(['chr', 'start', 'stop'], inplace=True)
df_cpg_wt_hues8.set_index(['chr', 'start', 'stop'], inplace=True)
df_cpg_dko.set_index(['chr', 'start', 'stop'], inplace=True)
df_cpg_tko.set_index(['chr', 'start', 'stop'], inplace=True)

########
# subtract dko
########
print('sub DKO')
path_dko_out = '/home/data/Shared/shared_datasets/wgbs/nlaszik/wt_dko_wgbs_intersect/GSM4458668.HUES64_GSM3662265.HUES64_DKO_WGBS_difference.bed'
subtract_meth(df_cpg_dko, df_cpg_wt_hues64, precision, path_dko_out)

########
# subtract tko
########
print('sub TKO')
path_tko_out = '/home/data/Shared/shared_datasets/wgbs/nlaszik/wt_tko_intersect/GSM3618718.HUES8_GSM3618720.HUES8_TKO_WGBS_difference.bed'
subtract_meth(df_cpg_tko, df_cpg_wt_hues8, precision, path_tko_out)
	

# Neural Dataset

# both these files are stranded, with strand information in subsequent CpGs
df_cpg_nsc_d0 = pd.read_csv('/home/data/Shared/shared_datasets/wgbs/data/GSE90553/data/GSM3039356_BiSeq_WT_D0_R2.bedGraph', names=['chr', 'start', 'stop', 'methylation'], sep='\t')
df_cpg_nsc_d0['methylation'] = df_cpg_nsc_d0['methylation'] * precision
df_cpg_nsc_d0[['methylation']] = df_cpg_nsc_d0[['methylation']].astype(int)

df_cpg_nsc_d14 = pd.read_csv('/home/data/Shared/shared_datasets/wgbs/data/GSE90553/data/GSM3039348_BiSeq_WT_D14_R1.sorted.bedGraph', names=['chr', 'start', 'stop', 'methylation'], sep='\t')
df_cpg_nsc_d14['methylation'] = df_cpg_nsc_d14['methylation'] * precision
df_cpg_nsc_d14[['methylation']] = df_cpg_nsc_d14[['methylation']].astype(int)

df_cpg_nsc_d0.set_index(['chr', 'start', 'stop'], inplace=True)
df_cpg_nsc_d14.set_index(['chr', 'start', 'stop'], inplace=True)

path_out = '/home/data/Shared/shared_datasets/wgbs/data/GSE90553/difference/D14_R1_D0_R2_difference.bed'
subtract_meth(df_cpg_nsc_d14, df_cpg_nsc_d0, precision, path_out)




