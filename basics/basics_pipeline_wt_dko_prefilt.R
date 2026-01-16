.libPaths(c("/home/nlaszik/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

library(Seurat)
library("anndata")
library(BASiCS)

# Define file paths
base_path <- "../preprocessing/filtered_transcript_counts"
samples <- c("wt_rep1", "wt_rep2", "dko_rep1", "dko_rep2")

# Convert h5ad files to h5seurat
h5ad_paths <- file.path(base_path, paste0(samples, ".h5ad"))
h5seurat_paths <- file.path(base_path, paste0(samples, ".h5seurat"))
names(h5ad_paths) <- samples
names(h5seurat_paths) <- samples

# Read anndata
data_list <- lapply(samples, function(sample) {
	read_h5ad(h5ad_paths[sample])
})

# Create Seurat objects
seurat_list <- lapply(data_list, function(data) {
	CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
})
names(seurat_list) <- samples

# Extract count matrices
counts_list <- lapply(seurat_list, function(seurat_obj) {
  counts_sparse <- GetAssayData(seurat_obj, layer = "counts")
  as.matrix(counts_sparse)
})

# Apply BASiCS_Filter (tech parameter not needed for filtering)
# MinAvCountsPerCellsWithExpression = 1 filters out genes/cells with very low counts
filtered_list <- lapply(counts_list, function(counts) {
  BASiCS_Filter(counts, Tech = rep(FALSE, nrow(counts)), 
				MinAvCountsPerCellsWithExpression = 1)$Counts
})

# Find common genes across all samples
expressed_genes <- Reduce(union, lapply(filtered_list, rownames))
cat("Number of expressed genes:", length(expressed_genes), "\n")

# Get genes present in all samples
common_genes <- Reduce(intersect, lapply(filtered_list, rownames))
cat("Number of common genes:", length(common_genes), "\n")

# Subset all matrices to common genes
filtered_list <- lapply(filtered_list, function(mat) {
  mat[common_genes, , drop = FALSE]
})

# Merge replicates by condition
filtered_wt <- cbind(filtered_list$wt_rep1, filtered_list$wt_rep2)
filtered_dko <- cbind(filtered_list$dko_rep1, filtered_list$dko_rep2)

# Create BASiCS data objects with batch information
basics_data_wt <- newBASiCS_Data(
  filtered_wt, 
  BatchInfo = rep(c(1, 2), c(ncol(filtered_list$wt_rep1), ncol(filtered_list$wt_rep2)))
)

basics_data_dko <- newBASiCS_Data(
  filtered_dko, 
  BatchInfo = rep(c(1, 2), c(ncol(filtered_list$dko_rep1), ncol(filtered_list$dko_rep2)))
)

# Run BASiCS MCMC
# Note: Regression = FALSE means residual over-dispersion may be confounded by mean expression
# Consider setting Regression = TRUE for better estimates (recommended)
basics_results_wt <- BASiCS_MCMC(
  basics_data_wt, 
  Threads = 96, 
  N = 4000, 
  Thin = 10, 
  Burn = 2000, 
  WithSpikes = FALSE, 
  Regression = FALSE, 
  PrintProgress = TRUE, 
  PriorParam = BASiCS_PriorParam(basics_data_wt, PriorDelta = 'log-normal')
)

basics_results_dko <- BASiCS_MCMC(
  basics_data_dko, 
  Threads = 96, 
  N = 4000, 
  Thin = 10, 
  Burn = 2000, 
  WithSpikes = FALSE, 
  Regression = FALSE, 
  PrintProgress = TRUE, 
  PriorParam = BASiCS_PriorParam(basics_data_dko, PriorDelta = 'log-normal')
)

# Differential expression testing
de_test <- BASiCS_TestDE(
  Chain1 = basics_results_dko, 
  Chain2 = basics_results_wt, 
  GroupLabel1 = "DKO", 
  GroupLabel2 = "WT", 
  EpsilonM = log2(1.5), 
  EpsilonD = log2(1.5), 
  ProbThresholdM = 0.85, 
  ProbThresholdD = 0.85, 
  Plot = TRUE
)

# Merge results
mean_df <- as.data.frame(de_test@Results$Mean@Table)
disp_df <- as.data.frame(de_test@Results$Disp@Table)
merged_df <- merge(mean_df, disp_df, by = "GeneName")

# Save results
path_out <- "dko_wt_basics_results.csv"
write.csv(merged_df, file = path_out, row.names = FALSE)

cat("Analysis complete. Results saved to:", path_out, "\n")