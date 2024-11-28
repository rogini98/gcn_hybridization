# run_simulation.R
# Enhanced script to run introgression simulations with network analysis
# Author: Original by Dan Bolnick, modified by Rogini Runghen
# Last updated: November 25, 2024

{
# Clear workspace and set random seed for reproducibility
rm(list = ls())
set.seed(123)
  
#############################################
### Setup
#############################################
# Load required libraries
library(ggplot2)
library(patchwork)
library(parallel)
library(pbapply)
library(igraph)
library(Matrix)

# Source required functions
source("01-network_functions_introgressionSim.R")
source("02-visualization_functions_introgression.R")
source("03-simulation_runners_introgression.R")
source("04-network_analysis_introgression.R")

# Create output directories
output_dir <- "output-introgression/"
for(dir in c("figures", "data", "network_analysis",
             "network_analysis/averaged_heatmaps", 
             "network_analysis/averaged_networks",
             "network_analysis/metrics")) {
  dir.create(file.path(output_dir, dir), recursive = TRUE, showWarnings = FALSE)
}

#############################################
### Run Simulations
#############################################
{
# Initialize parameters
params <- initialize_parameters()
saveRDS(params, file.path(output_dir, "data", "simulation_parameters.rds"))

# Print parameter summary
cat("\nSimulation Parameters:\n")
cat("Number of genes:", sum(params$Chr.lengths), "\n")
cat("Population size:", params$PopSize, "\n")
cat("Generations:", params$Generations, "\n")
cat("Network density (00):", params$network_params$density00, "\n")
cat("Network density (11):", params$network_params$density11, "\n")
cat("NetQTL locus:", params$neteQTLlocus, "\n")


#############################################
### Run Network Simulations
#############################################
cat("\nRunning network simulations...\n")

# Initialize network metrics storage
network_metrics <- data.frame(
  migration_rate = double(),
  replicate = integer(),
  mean_degree = double(),
  transitivity = double(),
  modularity = double()
)

# Run simulations with network analysis
migration_rates_network <- exp(seq(-10, 0, by = 0.01))
#migration_rates_network <- exp(seq(-10, 0, by = 1))
threshold <- 0.2

network_results <- run_network_simulations(
  migration_rates = migration_rates_network,
  params = params,
  n_replicates = 1,
  output_dir = file.path(output_dir, "network_analysis")
)

#############################################
### Process Network Analysis Results
#############################################
cat("\nProcessing network analysis results...\n")

# Process each migration rate
for(i in seq_along(network_results)) {
  rate <- migration_rates_network[i]
  result <- network_results[[i]]
  
  # Calculate average correlation matrix across replicates
  avg_cor <- calculate_average_correlation(result)
  
  if(!is.null(avg_cor)) {
    # Generate and save averaged heatmap
    png(file.path(output_dir, "network_analysis/averaged_heatmaps", 
                  sprintf("avg_heatmap_mig%.6f.png", rate)),
        width = 800, height = 800, res = 100)
    plot_averaged_heatmap(avg_cor, params$Chr.lengths, rate)
    dev.off()
    
    # Generate and save averaged network plot
    png(file.path(output_dir, "network_analysis/averaged_networks", 
                  sprintf("avg_network_mig%.6f.png", rate)),
        width = 800, height = 800, res = 100)
    par(mar = c(1,1,3,1))
    plot_averaged_network(avg_cor, params$Chr.lengths, threshold)
    dev.off()
  }
  
  # Calculate and store network metrics for each replicate
  for(rep in 1:length(result$replicates)) {
    if(!is.null(result$replicates[[rep]]$transcriptomes)) {
      cor_matrix <- cor(result$replicates[[rep]]$transcriptomes$r)
      net_stats <- calculate_network_stats(cor_matrix, threshold)
      
      network_metrics <- rbind(network_metrics, data.frame(
        migration_rate = rate,
        replicate = rep,
        mean_degree = net_stats$mean_degree,
        transitivity = net_stats$transitivity,
        modularity = net_stats$modularity
      ))
    }
  }
}

# Save network metrics
write.csv(network_metrics,
          file.path(output_dir, "network_analysis/metrics", "network_metrics.csv"),
          row.names = FALSE)

# Create network metrics plots
cat("\nGenerating network metrics visualizations...\n")
network_plots <- plot_network_metrics_with_interpretation(
  metrics_df = network_metrics,
  output_dir = file.path(output_dir, "network_analysis")
)

# Create interpretation summary
trends <- create_interpretation_summary(
  metrics_df = network_metrics,
  output_dir = file.path(output_dir, "network_analysis")
)


#################################################################
### Run Introgression with Allele frequency Tracking Simulations
#################################################################
cat("\nRunning tracking simulations...\n")
tracking_results <- run_tracked_simulations(
  migration_rates = c(0, 0.05, 0.5, 1),
  params = params,
  n_replicates = 10,
  output_dir = file.path(output_dir, "tracking")
)

#############################################
### Generate SI Figures
#############################################
# Figure S4: Grid plot of all migration rates
cat("\nGenerating grid plot...\n")
grid_plot <- plot_grid_introgression(
  tracking_results$results,
  params$Chr.lengths,
  n_loci_ci = 5
)

ggsave(
  filename = file.path(output_dir, "figures", "migration_rates_grid.pdf"),
  plot = grid_plot,
  width = 15,
  height = 12
)

#############################################
### Generate Summary Report
#############################################

cat("\nGenerating summary report...\n")
sink(file.path(output_dir, "simulation_summary.txt"))

cat("Introgression Simulation Summary\n")
cat("==============================\n\n")

cat("Parameters:\n")
cat("-----------\n")
cat("Number of genes:", sum(params$Chr.lengths), "\n")
cat("Population size:", params$PopSize, "\n")
cat("Generations:", params$Generations, "\n")
cat("Chromosome lengths:", paste(params$Chr.lengths, collapse = ", "), "\n")
cat("Network parameters:\n")
cat("  Density (00):", params$network_params$density00, "\n")
cat("  Effect (00):", params$network_params$effect00, "\n")
cat("  Density (11):", params$network_params$density11, "\n")
cat("  Effect (11):", params$network_params$effect11, "\n")
cat("NetQTL locus:", params$neteQTLlocus, "\n\n")

# Modified network analysis section
cat("\nNetwork Analysis Results:\n")
cat("----------------------\n")
cat("Number of migration rates tested:", length(migration_rates_network), "\n")
cat("Migration rate range:", sprintf("%.6f to %.6f\n", 
                                     min(migration_rates_network), 
                                     max(migration_rates_network)))
cat("Network metrics tracked: mean degree, transitivity, modularity\n")
cat("Correlation threshold used:", threshold, "\n\n")

cat("Network Analysis Files:\n")
cat("   - Average correlation heatmaps (avg_heatmap_mig*.png)\n")
cat("   - Average network plots (avg_network_mig*.png)\n")
cat("   - Network metrics data (network_metrics.csv)\n")
cat("   - Network metrics summary plot (network_metrics_summary.pdf)\n")
cat("   - Network trends analysis (network_trends.csv)\n")
cat("   - Network interpretation report (network_analysis_report.txt)\n")

sink()


cat("Simulation Results:\n")
cat("-----------------\n")
cat("Network simulations:\n")
cat("  Migration rates tested:", length(network_results), "\n")
cat("  Replicates per rate:", network_results[[1]]$parameters$n_replicates, "\n\n")

cat("Tracking simulations:\n")
cat("  Migration rates tested:", length(tracking_results$results), "\n")
cat("  Replicates per rate:", tracking_results$results[[1]]$parameters$n_replicates, "\n\n")

cat("Output Files:\n")
cat("-------------\n")
cat("1. Figures:\n")
cat("   - Combined visualization (combined_migration_rates.pdf)\n")
cat("   - Heatmap panel (heatmap_panel.pdf)\n")
cat("   - Trajectories panel (trajectories_panel.pdf)\n")
cat("   - Migration rates grid (migration_rates_grid.pdf)\n")
cat("   - Individual replicate comparisons (replicates_rate_*.pdf)\n\n")

cat("2. Data:\n")
cat("   - Simulation parameters (simulation_parameters.rds)\n")
cat("   - Network simulation results (network_simulation_results.rds)\n")
cat("   - Tracking simulation results (all_results.rds)\n")

sink()

cat("\nSimulation completed successfully!\n")
cat("Results saved in:", output_dir, "\n")
cat("See simulation_summary.txt for detailed information.\n")
}
}
