#===============================
# Setup and Configuration
#===============================
# Install required packages if needed
required_packages <- c("BiocManager", "data.table", "ggplot2", "viridis", "reshape2")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

if (!require("qvalue")) {
  BiocManager::install("qvalue")
}

if (!require("HiClimR")) {
  BiocManager::install("HiClimR")
}

#===============================
# Source Function Scripts
#===============================
# Source the main analysis functions
source("eqtl_gcn_analysis_optim.R")  # Contains the core analysis functions
source("eqtl_gcn_plots.R")  # Contains the plotting functions

#===============================
# Set Analysis Parameters
#===============================
params <- list(
  data_file = "data/MergedTagseqQTL_GENOCORRECTED.csv",
  gene_map_file = "data/MapENSGACT_ENSGACG.csv",
  output_dir = "results",
  batch_size = 50,
  min_reads = 30,
  min_presence = 0.8,
  eqtl_significance = 0.05,
  correlation_threshold = 0.7
)

#===============================
# Create Output Directories
#===============================
dir.create(params$output_dir, showWarnings = FALSE)
dir.create(file.path(params$output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(params$output_dir, "output"), showWarnings = FALSE)

#===============================
# Run Analysis Pipeline
#===============================
# Run main analysis
cat("Starting eQTL analysis pipeline...\n")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

# Run the analysis
results <- run_analysis(
  data_file = params$data_file,
  gene_map_file = params$gene_map_file,
  batch_size = params$batch_size
)

# Save results
cat("Saving analysis results...\n")
saveRDS(results, file.path(params$output_dir, "output", 
                           paste0("eqtl_results_", timestamp, ".rds")))

#===============================
# Generate Visualizations
#===============================
cat("Generating plots...\n")

# 1. eQTL Manhattan and Cis/Trans Plots
plot_eqtl_results(
  eqtl_results = results$eqtl,
  output_dir = file.path(params$output_dir, "plots"),
  significance_threshold = params$eqtl_significance
)

# 2. Co-expression Heatmap and Network
plot_coexpression(
  cor_matrix = results$coexpression$correlations,
  gene_map = gene_map,
  output_dir = file.path(params$output_dir, "plots"),
  correlation_threshold = params$correlation_threshold
)

# 3. Chromosome Analysis Plots
plot_chromosome_analysis(
  chr_correlations = results$coexpression,
  output_dir = file.path(params$output_dir, "plots")
)

#===============================
# Generate Summary Report
#===============================
cat("Generating summary report...\n")

# Create summary statistics
summary_stats <- list(
  total_genes = ncol(results$eqtl),
  significant_eqtls = sum(results$eqtl$p_value < params$eqtl_significance),
  mean_correlation = mean(abs(results$coexpression$correlations), na.rm = TRUE),
  chromosome_test_pvalue = results$coexpression$chromosome_test$p.value
)

# Write summary report
report <- c(
  "eQTL Analysis Summary Report",
  "========================",
  paste("Analysis Date:", Sys.Date()),
  paste("Total Genes Analyzed:", summary_stats$total_genes),
  paste("Significant eQTLs found:", summary_stats$significant_eqtls),
  paste("Mean absolute correlation:", round(summary_stats$mean_correlation, 3)),
  paste("Chromosome test p-value:", format(summary_stats$chromosome_test_pvalue, scientific = TRUE)),
  "\nAnalysis Parameters:",
  paste("Minimum reads:", params$min_reads),
  paste("Minimum presence:", params$min_presence),
  paste("eQTL significance threshold:", params$eqtl_significance),
  paste("Correlation threshold:", params$correlation_threshold)
)

writeLines(report, file.path(params$output_dir, 
                             paste0("analysis_summary_", timestamp, ".txt")))

cat("Analysis pipeline completed!\n")
cat("Results saved in:", params$output_dir, "\n")