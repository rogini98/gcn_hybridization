# 04-network_analysis_introgression.R
# Enhanced network analysis tracking and visualization
# Author: Modified based on original code
# Last updated: November 25, 2024

library(igraph)
library(Matrix)
library(ggplot2)
library(patchwork)

#############################################
### Network Analysis Functions
#############################################

#' Calculate network metrics for simulation output
#' @param co_expr Correlation matrix
#' @param threshold Correlation threshold
#' @return List of network metrics
calculate_network_stats <- function(co_expr, threshold = 0.2) {
  # Calculate adjacency matrix
  adjacency <- (cor(co_expr) > threshold) * 1
  diag(adjacency) <- 0
  
  # Create igraph object
  g <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
  
  # Calculate metrics
  degree_vector <- degree(g)
  mean_degree <- mean(degree_vector)
  transitivity <- transitivity(g, type = "global")
  
  # Calculate modularity
  comm <- cluster_louvain(g)
  modularity <- modularity(comm)
  
  # Handle NaN values
  if(is.nan(transitivity)) transitivity <- 0
  if(is.nan(modularity)) modularity <- 0
  
  return(list(
    mean_degree = mean_degree,
    transitivity = transitivity,
    modularity = modularity,
    graph = g
  ))
}

#' Plot gene network
#' @param co_expr Correlation matrix
#' @param Chr.lengths Vector of chromosome lengths
#' @param threshold Correlation threshold
#' @return Plot of gene network
plot_gene_network <- function(co_expr, Chr.lengths, threshold = 0.2) {
  # Create adjacency matrix
  adjacency <- (cor(co_expr) > threshold) * 1
  diag(adjacency) <- 0
  
  # Create igraph object
  g <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
  
  # Calculate node colors based on chromosomes
  chr_ends <- cumsum(Chr.lengths)
  chr_starts <- c(1, chr_ends[-length(chr_ends)] + 1)
  node_colors <- rep(NA, vcount(g))
  chr_colors <- rainbow(length(Chr.lengths))
  
  for(i in 1:length(Chr.lengths)) {
    node_colors[chr_starts[i]:chr_ends[i]] <- chr_colors[i]
  }
  
  # Create plot
  plot(g,
       layout = layout_in_circle(g),
       vertex.size = 5,
       vertex.label = NA,
       vertex.color = node_colors,
       edge.width = 1,
       edge.color = "gray70")
  
  return(g)
}

#############################################
### Main Simulation Analysis
#############################################

run_network_analysis <- function(migration_rates = exp(seq(-10, 0, by = 0.01)),
                                 Nreps = 1,
                                 Fgen = 100,
                                 threshold = 0.2,
                                 output_dir = "results/network_analysis") {
  
  # Create output directories
  for(dir in c("heatmaps", "networks", "data")) {
    dir.create(file.path(output_dir, dir), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Initialize metrics storage
  network_metrics <- data.frame(
    migration_rate = double(),
    replicate = integer(),
    mean_degree = double(),
    transitivity = double(),
    modularity = double()
  )
  
  # Initialize parameters 
  params <- initialize_parameters()
  
  # Run simulations for each migration rate
  for(mig_rate in migration_rates) {
    cat("\nProcessing migration rate:", mig_rate, "\n")
    
    for(rep in 1:Nreps) {
      cat("Replicate:", rep, "\n")
      
      # Run simulation
      output <- IntrogressionSim(
        Generations = Fgen,
        Chr.lengths = params$Chr.lengths,
        PopSize = 500,
        intercepts = params$intercepts,
        betas = params$betas,
        effect_matrices = params$effect_matrices,
        neteQTLlocus = params$neteQTLlocus,
        migration_rate = mig_rate
      )
      
      # Calculate network metrics
      net_stats <- calculate_network_stats(output$co_expr, threshold)
      
      # Store metrics
      network_metrics <- rbind(network_metrics, data.frame(
        migration_rate = mig_rate,
        replicate = rep,
        mean_degree = net_stats$mean_degree,
        transitivity = net_stats$transitivity,
        modularity = net_stats$modularity
      ))
      
      # Generate and save heatmap
      png(file.path(output_dir, "heatmaps", sprintf("heatmap_mig%.6f_rep%d.png", mig_rate, rep)),
          width = 500, height = 500)
      heatmap(output$co_expr,
              Rowv = NA,
              Colv = NA,
              main = sprintf("Migration Rate: %.6f", mig_rate))
      dev.off()
      
      # Generate and save network plot
      png(file.path(output_dir, "networks", sprintf("network_mig%.6f_rep%d.png", mig_rate, rep)),
          width = 500, height = 500)
      par(mar = c(0,0,0,0))
      plot_gene_network(output$co_expr, params$Chr.lengths, threshold)
      dev.off()
      
      # Save simulation output
      saveRDS(output, 
              file.path(output_dir, "data", sprintf("simulation_mig%.6f_rep%d.rds", mig_rate, rep)))
    }
  }
  
  # Save network metrics
  write.csv(network_metrics,
            file.path(output_dir, "data", "network_metrics.csv"),
            row.names = FALSE)
  
  # Create summary plots
  create_summary_plots(network_metrics, output_dir)
  
  return(network_metrics)
}

#' Create summary plots of network metrics
#' @param network_metrics Data frame of network metrics
#' @param output_dir Output directory
create_summary_plots <- function(network_metrics, output_dir) {
  # Create long format for plotting
  metrics_long <- network_metrics %>%
    pivot_longer(
      cols = c(mean_degree, transitivity, modularity),
      names_to = "metric",
      values_to = "value"
    )
  
  # Create plot
  p <- ggplot(metrics_long, aes(x = migration_rate, y = value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    scale_x_log10() +
    labs(x = "Migration Rate (log scale)",
         y = "Metric Value",
         title = "Network Metrics vs Migration Rate") +
    theme_minimal()
  
  # Save plot
  ggsave(
    file.path(output_dir, "network_metrics_summary.pdf"),
    p,
    width = 10,
    height = 12
  )
}
