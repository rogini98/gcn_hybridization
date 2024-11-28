#########################################
# Master script to run both hybridization and introgression simulations
# Author: Rogini Runghen
# Last updated: November 28, 2024
#########################################

# Record start time
start_time <- Sys.time()
cat("Starting simulations at:", format(start_time), "\n\n")

#########################################
# Install and load required packages
#########################################
required_packages <- c(
  "igraph", "ggplot2", "dplyr", "tidyr", "patchwork", 
  "parallel", "doParallel", "foreach", "pheatmap", 
  "RColorBrewer", "pbapply", "Matrix", "tidyverse"
)

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#########################################
# Create master output directory
#########################################
main_output_dir <- "simulation_results"
dir.create(main_output_dir, showWarnings = FALSE)

# Create log file
log_file <- file(file.path(main_output_dir, "complete_simulation_log.txt"), open = "w")
sink(log_file, type = "output")
sink(log_file, type = "message")

#########################################
# Run Hybridization Simulation
#########################################
cat("\nStarting Hybridization Simulation...\n")
cat("=====================================\n\n")

tryCatch({
  source("01-run_hybridization_simulation.R")
  cat("Hybridization simulation completed successfully!\n")
}, error = function(e) {
  cat("Error in hybridization simulation:", conditionMessage(e), "\n")
}, warning = function(w) {
  cat("Warning in hybridization simulation:", conditionMessage(w), "\n")
})

#########################################
# Run Introgression Simulation
#########################################
cat("\nStarting Introgression Simulation...\n")
cat("====================================\n\n")

tryCatch({
  source("02-run_introgression_simulation.R")
  cat("Introgression simulation completed successfully!\n")
}, error = function(e) {
  cat("Error in introgression simulation:", conditionMessage(e), "\n")
}, warning = function(w) {
  cat("Warning in introgression simulation:", conditionMessage(w), "\n")
})

#########################################
# Create Combined Summary
#########################################
cat("\nGenerating Combined Summary...\n")
cat("=============================\n\n")

# Calculate total runtime
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

# Create combined summary
sink(file.path(main_output_dir, "combined_summary.txt"))

cat("Combined Simulation Summary\n")
cat("=========================\n\n")
cat("Run Date:", format(start_time, "%Y-%m-%d"), "\n")
cat("Start Time:", format(start_time, "%H:%M:%S"), "\n")
cat("End Time:", format(end_time, "%H:%M:%S"), "\n")
cat("Total Runtime:", round(total_time, 2), "minutes\n\n")

cat("Output Directories:\n")
cat("------------------\n")
cat("1. Hybridization Results: output-hybridization/\n")
cat("2. Introgression Results: output-introgression/\n\n")

cat("Simulation Parameters:\n")
cat("--------------------\n")
cat("Chromosome lengths:", paste(Chr.lengths, collapse = ", "), "\n")
cat("Population size:", PopSize, "\n")
cat("Number of replicates:", n_replicates, "\n")
cat("Correlation threshold:", threshold, "\n\n")

cat("Files Generated:\n")
cat("---------------\n")
cat("1. Hybridization:\n")
cat("   - Network plots and heatmaps\n")
cat("   - Evolution and relationship plots\n")
cat("   - Complete analysis data\n\n")
cat("2. Introgression:\n")
cat("   - Network analysis results\n")
cat("   - Migration rate effects\n")
cat("   - Tracking simulation results\n\n")

sink()

# Close log file
sink(type = "output")
sink(type = "message")
close(log_file)

cat("\nAll simulations completed successfully!\n")
cat("Total runtime:", round(total_time, 2), "minutes\n")
cat("Results saved in:", main_output_dir, "\n")
cat("See 'combined_summary.txt' for detailed information.\n")