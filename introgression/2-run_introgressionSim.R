# simulation_runner.R
# Script to run the introgression simulations
# Authors: Original by Dan Bolnick, reorganized by Claude
# Last updated: November 26, 2024

run_full_simulation <- function(output_dir, params, n_replicates = 10,
                                network_migration_rates, tracking_migration_rates) {
  # Create output directories
  for(dir in c("figures", "data", "network_analysis", 
               "network_analysis/averaged_heatmaps",
               "network_analysis/averaged_networks",
               "network_analysis/metrics")) {
    dir.create(file.path(output_dir, dir), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save parameters
  saveRDS(params, file.path(output_dir, "data", "simulation_parameters.rds"))
  
  # Print simulation parameters
  cat("\nStarting simulation with parameters:\n")
  cat("Number of genes:", sum(params$Chr.lengths), "\n")
  cat("Population size:", params$PopSize, "\n")
  cat("Generations:", params$Generations, "\n")
  cat("Network analysis migration rates:", length(network_migration_rates), "\n")
  cat("Tracking analysis migration rates:", length(tracking_migration_rates), "\n")
  
  # Run network simulations
  cat("\nRunning network simulations...\n")
  network_results <- run_network_simulations(
    migration_rates = network_migration_rates,
    params = params,
    n_replicates = n_replicates,
    output_dir = file.path(output_dir, "network_analysis")
  )
  
  # Run tracking simulations
  cat("\nRunning tracking simulations...\n")
  tracking_results <- run_tracked_simulations(
    migration_rates = tracking_migration_rates,
    params = params,
    n_replicates = n_replicates,
    output_dir = file.path(output_dir, "tracking")
  )
  
  # Extract and compile network metrics
  network_metrics <- data.frame(
    migration_rate = double(),
    replicate = integer(),
    mean_degree = double(),
    transitivity = double(),
    modularity = double()
  )
  
  # Process network results and collect metrics
  for(i in seq_along(network_results)) {
    rate <- network_migration_rates[i]
    result <- network_results[[as.character(rate)]]
    
    for(rep in seq_along(result$replicates)) {
      if(!is.null(result$replicates[[rep]]$transcriptomes)) {
        cor_matrix <- cor(result$replicates[[rep]]$transcriptomes$r)
        net_stats <- calculate_network_stats(cor_matrix, threshold = 0.2)
        
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
  write.csv(
    network_metrics,
    file.path(output_dir, "network_analysis", "metrics", "network_metrics.csv"),
    row.names = FALSE
  )
  
  return(list(
    params = params,
    network_results = network_results,
    tracking_results = tracking_results,
    network_metrics = network_metrics
  ))
}