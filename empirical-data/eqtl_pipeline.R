#===============================
# Setup and Configuration
#===============================
rm(list = ls())

# Remove old checkpoint files (if they exist)
unlink("checkpoints/*")
unlink("results/*")

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
source("eqtl_analysis_checkpoints.R")  # Contains the core analysis functions
source("eqtl_gcn_plots.R")  # Contains the plotting functions

#===============================
# Set Analysis Parameters
#===============================
params <- list(
  data_file = "data/MergedTagseqQTL_GENOCORRECTED.csv",
  gene_map_file = "data/MapENSGACT_ENSGACG.csv",
  output_dir = "results",
  checkpoint_dir = "checkpoints",
  batch_size = 50,
  min_reads = 30,
  min_presence = 0.8,
  eqtl_significance = 0.05,
  correlation_threshold = 0.7,
  start_from = "beginning"
)

#===============================
# Create Output Directories
#===============================
dir.create(params$output_dir, showWarnings = FALSE)
dir.create(file.path(params$output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(params$output_dir, "output"), showWarnings = FALSE)
dir.create(params$checkpoint_dir, showWarnings = FALSE)

#===============================
# Run Analysis Pipeline
#===============================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

tryCatch({
  cat("Starting eQTL analysis pipeline...\n")
  
  # Run the analysis with checkpoints
  results <- run_analysis(
    data_file = params$data_file,
    gene_map_file = params$gene_map_file,
    batch_size = params$batch_size,
    checkpoint_dir = params$checkpoint_dir,
    start_from = params$start_from
  )
  
  # Save final results
  cat("Saving analysis results...\n")
  saveRDS(results, file.path(params$output_dir, "output", 
                             paste0("eqtl_results_", timestamp, ".rds")))
  
  #===============================
  # Generate Visualizations
  #===============================
  cat("Generating plots...\n")
  
  tryCatch({
    # 1. eQTL Manhattan Plot
    cat("Creating Manhattan plot...\n")
    plot_eqtl_results(
      eqtl_results = results$eqtl,
      output_dir = file.path(params$output_dir, "plots"),
      significance_threshold = params$eqtl_significance
    )
    
    # 2. Co-expression and Chromosome Analysis
    if (!is.null(results$coexpression)) {
      cat("Creating co-expression plots...\n")
      plot_coexpression(
        cor_matrix = results$coexpression$correlations,
        gene_map = results$gene_map,
        output_dir = file.path(params$output_dir, "plots"),
        correlation_threshold = params$correlation_threshold
      )
      
      cat("Creating chromosome analysis plots...\n")
      plot_chromosome_analysis(
        chr_correlations = results$coexpression,
        output_dir = file.path(params$output_dir, "plots")
      )
      
      # Update summary statistics to include chromosome analysis
      summary_stats$chr_stats <- if (!is.null(results$coexpression$chr_stats)) {
        do.call(rbind, lapply(results$coexpression$chr_stats, 
                              function(x) data.frame(chromosome = x$chromosome,
                                                     mean_correlation = x$mean_correlation,
                                                     n_genes = x$n_genes)))
      } else NULL
      
      # Add chromosome statistics to report
      if (!is.null(summary_stats$chr_stats)) {
        report <- c(report,
                    "\nChromosome Analysis Results:",
                    paste("Number of chromosomes analyzed:", nrow(summary_stats$chr_stats)),
                    paste("Overall chromosome effect p-value:", 
                          format(results$coexpression$chromosome_test$p_value, 
                                 scientific = TRUE)))
      }
    } else {
      cat("Skipping co-expression and chromosome analysis: no results available\n")
    }
    
  }, error = function(e) {
    cat("Error in plotting:", e$message, "\n")
    print(traceback())
  })
  
  #===============================
  # Generate Summary Report
  #===============================
  cat("Generating summary report...\n")
  
  # Create summary statistics
  summary_stats <- list(
    total_genes = nrow(results$eqtl),
    significant_eqtls = sum(results$eqtl$p_value < params$eqtl_significance),
    analysis_start = params$start_from
  )
  
  # Write summary report
  report <- c(
    "eQTL Analysis Summary Report",
    "========================",
    paste("Analysis Date:", Sys.Date()),
    paste("Analysis Start Point:", summary_stats$analysis_start),
    paste("Total Genes Analyzed:", summary_stats$total_genes),
    paste("Significant eQTLs found:", summary_stats$significant_eqtls),
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
  
}, error = function(e) {
  cat("Error in analysis pipeline:", e$message, "\n")
  cat("Last checkpoint can be found in:", params$checkpoint_dir, "\n")
})