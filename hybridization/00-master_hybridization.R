#########################################
# Master script for network evolution analysis (Scenario 1)
# This script runs the entire workflow:
# 1. Loads functions
# 2. Runs simulations
# 3. Creates visualizations
#########################################

# Clear workspace and set random seed for reproducibility
rm(list = ls())
set.seed(123)

# Check and install required packages if needed
required_packages <- c("igraph", "ggplot2", "dplyr", "tidyr", 
                       "patchwork", "parallel", "doParallel", "foreach")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Create output directory if it doesn't exist
if(!dir.exists("output")) {
  dir.create("output")
}

# Create plots directory if it doesn't exist
if(!dir.exists("output/plots")) {
  dir.create("output/plots")
}

# Start logging
log_file <- file("output/analysis_log.txt", open = "w")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("Starting analysis:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

#########################################
# 1. Source functions
#########################################
cat("Loading functions...\n")
tryCatch({
  source("01-functions_hybridization.R")
  cat("Successfully loaded network functions\n")
}, error = function(e) {
  cat("Error loading network functions:", conditionMessage(e), "\n")
  quit(status = 1)
})

#########################################
# 2. Run simulation
#########################################
cat("\nRunning simulation...\n")
tryCatch({
  source("02-simulation_hybridization.R")
  cat("Successfully completed simulation\n")
}, error = function(e) {
  cat("Error in simulation:", conditionMessage(e), "\n")
  quit(status = 1)
})

#########################################
# 3. Create plots
#########################################
cat("\nGenerating plots...\n")
tryCatch({
  source("03-visualization_hybridization.R")
  cat("Successfully created plots\n")
}, error = function(e) {
  cat("Error creating plots:", conditionMessage(e), "\n")
  quit(status = 1)
})

#########################################
# 4. Save workspace and summary
#########################################

# Save complete workspace
save.image("output/complete_analysis.RData")

# Create summary statistics
summary_stats_final <- data.frame(
  Metric = c("Mean Degree", "Modularity", "Transitivity"),
  Parent_AA = c(
    parent_metrics$AA$mean_degree,
    parent_metrics$AA$modularity,
    parent_metrics$AA$transitivity
  ),
  Parent_BB = c(
    parent_metrics$BB$mean_degree,
    parent_metrics$BB$modularity,
    parent_metrics$BB$transitivity
  ),
  Final_Mean = c(
    tail(summary_stats$mean_degree_mean, 1),
    tail(summary_stats$modularity_mean, 1),
    tail(summary_stats$transitivity_mean, 1)
  ),
  Final_SE = c(
    tail(summary_stats$mean_degree_se, 1),
    tail(summary_stats$modularity_se, 1),
    tail(summary_stats$transitivity_se, 1)
  )
)

# Save summary statistics
write.csv(summary_stats_final, 
          "output/final_summary_statistics.csv", 
          row.names = FALSE)

# Print final summary
cat("\nFinal Summary Statistics:\n")
print(summary_stats_final)

# Print correlation matrix
cat("\nCorrelation Matrix:\n")
print(round(cor_matrix, 3))

# Print analysis parameters
cat("\nAnalysis Parameters:\n")
cat("Number of replicates:", n_replicates, "\n")
cat("Population size:", PopSize, "\n")
cat("Number of generations:", max(summary_stats$generation), "\n")
cat("Chromosome lengths:", paste(Chr.lengths, collapse = ", "), "\n")
cat("Network threshold:", threshold, "\n")

# End logging
cat("\nAnalysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink(type = "output")
sink(type = "message")
close(log_file)

# Print success message to console
cat("\nAnalysis complete! Results saved in output directory.\n")
cat("See 'output/analysis_log.txt' for detailed log.\n")
