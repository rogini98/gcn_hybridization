# master_script.R
# Master script to run the entire simulation pipeline
# Authors: Original by Dan Bolnick, reorganized by Claude
# Last updated: November 26, 2024

library(tictoc)
library(ggplot2)
library(patchwork)
library(parallel)
library(pbapply)
library(igraph)
library(Matrix)

# Source required scripts
source("1-base_functions_introgressionSim.R")
source("2-run_introgressionSim.R")
source("3-generate_introgression_plots.R")


tic()

# Simulation Parameters
config <- list(
  # Output directory
  output_dir = "results/introgression_simulation",
  
  # Basic simulation parameters
  chromosome_lengths = c(10, 10, 5, 5, 10),  # Length of each chromosome
  population_size = 500,                     # Size of population
  generations = 100,                         # Number of generations to simulate
  n_replicates = 1,                        # Number of replicates per migration rate
  
  # Network parameters
  network = list(
    density00 = 0.02,        # Network density for parent population A
    effect00 = 2,            # Network effect size for parent population A
    density11 = 0.03,        # Network density for parent population B
    effect11 = 1.5           # Network effect size for parent population B
  ),
  
  # Expression parameters
  expression = list(
    mean = 5,                # Mean expression level
    sd = 2,                  # Standard deviation of expression
    beta_sd = 3              # Standard deviation of genetic effects
  ),
  
  # Migration rates
  migration_rates = list(
    # Rates for network analysis (log scale)
    network = exp(seq(-10, 0, by = 1)),
    
    # Rates for tracking analysis (specific values)
    tracking = c(0, 0.05, 0.5, 1)
  ),
  
  # Visualization parameters
  visualization = list(
    correlation_threshold = 0.2,  # Threshold for network visualization
    n_loci_ci = 5                # Number of loci for confidence intervals
  )
)

# Initialize parameters with configuration
params <- initialize_parameters(
  Chr.lengths = config$chromosome_lengths,
  PopSize = config$population_size,
  Generations = config$generations,
  network_params = config$network,
  expression_params = config$expression
)

# Run full simulation
cat("\nStarting simulation...\n")
results <- run_full_simulation(
  output_dir = config$output_dir,
  params = params,
  n_replicates = config$n_replicates,
  network_migration_rates = config$migration_rates$network,
  tracking_migration_rates = config$migration_rates$tracking
)


# Create visualizations
viz_results <- create_all_visualizations(results, output_dir)

# Generate summary report
cat("\nGenerating summary report...\n")
generate_summary_report(
  results = results,
  config = config,
  output_dir = config$output_dir
)
toc()

cat("\nSimulation completed successfully!\n")
cat("Results saved in:", output_dir, "\n")
cat("See simulation_summary.txt for detailed information.\n")