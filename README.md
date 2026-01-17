# General

This repository contains Jupyter Notebooks and Python scripts for processing data and generating figures for "DNMT3s and TETs adjust CpG methylation canyon width to regulate gene expression and cell fate."

Some processed data is provided, such as BASiCS output and UMRs generated with PopCanyon.

# Dependencies

1. scRNA-seq and methylation data must be downloaded from GEO using get_data.sh
2. Python dependencies must be installed with requirements.txt. Development was done in Python 3.12.
3. Jupyter must be installed
4. For gene pair correlations, a GPU with CUDA is required.
5. For BASiCS processing, R 4.4.1 was used. Requirements are not provided, but processed data is stored in the baiscs directory

