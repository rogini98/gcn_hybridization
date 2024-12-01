# visualization_functions.R
# Functions for visualizing introgression simulation results
# Authors: Original by Dan Bolnick, reorganized by Rogini Runghen
# Last updated: November 25, 2024

library(ggplot2)
library(patchwork)
library(gridExtra)

#' Plot replicated heatmap of allele frequencies
#' @param sim_results Results from run_replicated_introgression
#' @param Chr.lengths Vector of chromosome lengths
#' @return ggplot object
plot_replicated_heatmap <- function(sim_results, Chr.lengths = NULL) {
  freq_matrix <- sim_results$summary$mean_frequencies
  
  df <- data.frame(
    Generation = rep(1:nrow(freq_matrix), ncol(freq_matrix)),
    Locus = rep(1:ncol(freq_matrix), each = nrow(freq_matrix)),
    Frequency = as.vector(freq_matrix)
  )
  
  mig_rate_formatted <- sprintf("%.3f", sim_results$parameters$migration_rate)
  
  p <- ggplot(df, aes(x = Locus, y = Generation, fill = Frequency)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "white",
      mid = "skyblue",
      high = "darkblue",
      midpoint = 0.5,
      limits = c(0, 1)
    ) +
    labs(
      title = paste0("Mean Allele Frequency Evolution\n",
                     sim_results$parameters$n_replicates,
                     " Replicates, Migration Rate = ",
                     mig_rate_formatted),
      x = "Locus",
      y = "Generation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  if (!is.null(Chr.lengths)) {
    chr_boundaries <- cumsum(Chr.lengths)
    chr_centers <- c(0, chr_boundaries[-length(chr_boundaries)]) + Chr.lengths/2
    
    p <- p +
      geom_vline(xintercept = chr_boundaries,
                 linetype = "dashed",
                 alpha = 0.5) +
      scale_x_continuous(
        breaks = chr_centers,
        labels = paste("Chr", 1:length(Chr.lengths))
      )
  }
  
  return(p)
}

#' Plot confidence intervals for selected loci
#' @param sim_results Results from run_replicated_introgression
#' @param n_loci Number of loci to plot
#' @return ggplot object
plot_confidence_intervals <- function(sim_results, n_loci = 5) {
  total_loci <- ncol(sim_results$summary$mean_frequencies)
  if (n_loci > total_loci) n_loci <- total_loci
  selected_loci <- round(seq(1, total_loci, length.out = n_loci))
  
  plot_data <- data.frame()
  generations <- 1:nrow(sim_results$summary$mean_frequencies)
  
  for(locus in selected_loci) {
    plot_data <- rbind(plot_data, data.frame(
      Generation = generations,
      Locus = paste("Locus", locus),
      Mean = sim_results$summary$mean_frequencies[, locus],
      Lower = sim_results$summary$ci_lower[, locus],
      Upper = sim_results$summary$ci_upper[, locus]
    ))
  }
  
  mig_rate_formatted <- sprintf("%.3f", sim_results$parameters$migration_rate)
  blue_palette <- colorRampPalette(c("lightblue", "darkblue"))(n_loci)
  
  p <- ggplot(plot_data, aes(x = Generation, y = Mean,
                             color = Locus, fill = Locus)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper),
                alpha = 0.2, color = NA) +
    scale_color_manual(values = blue_palette) +
    scale_fill_manual(values = blue_palette) +
    labs(
      title = paste0("Allele Frequency Trajectories with 95% CI\n",
                     sim_results$parameters$n_replicates,
                     " Replicates, Migration Rate = ",
                     mig_rate_formatted),
      x = "Generation",
      y = "Allele Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, 1))
  
  return(p)
}

#' Create combined visualization of introgression results
#' @param sim_results Results from run_replicated_introgression
#' @param Chr.lengths Vector of chromosome lengths
#' @param n_loci_ci Number of loci to show in confidence interval plot
#' @return Combined ggplot object
plot_introgression_summary <- function(sim_results, Chr.lengths = NULL, n_loci_ci = 5) {
  p1 <- plot_replicated_heatmap(sim_results, Chr.lengths)
  p2 <- plot_confidence_intervals(sim_results, n_loci_ci)
  
  mig_rate_formatted <- sprintf("%.3f", sim_results$parameters$migration_rate)
  
  combined_plot <- p1 / p2 +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = paste0("Introgression Simulation Results (Migration Rate = ",
                     mig_rate_formatted, ")"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    )
  
  return(combined_plot)
}


#' Create grid plot for all migration rates
#' @param results_list List of results for different migration rates
#' @param Chr.lengths Vector of chromosome lengths
#' @param n_loci_ci Number of loci for confidence interval plots
#' @return Combined ggplot object
plot_grid_introgression <- function(results_list, Chr.lengths, n_loci_ci = 5) {
  plots <- list()
  
  for(rate in names(results_list)) {
    p <- plot_introgression_summary(
      results_list[[rate]], 
      Chr.lengths, 
      n_loci_ci = n_loci_ci
    ) +
      ggtitle(paste("Migration rate =", rate)) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plots[[rate]] <- p
  }
  
  combined_plot <- wrap_plots(plots, ncol = 2, nrow = 2) +
    plot_annotation(
      title = "Introgression Results Across Migration Rates",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  return(combined_plot)
}

#' Create combined visualization for multiple migration rates
#' @param results_list List of simulation results
#' @param output_dir Output directory
#' @return List of plots
create_combined_visualizations <- function(results_list, output_dir) {
  plot_data <- data.frame()
  ci_data <- data.frame()
  
  for(rate in names(results_list)) {
    result <- results_list[[rate]]
    
    # Add mean frequencies to plot data
    freq_matrix <- result$summary$mean_frequencies
    temp_data <- data.frame(
      Generation = rep(1:nrow(freq_matrix), ncol(freq_matrix)),
      Locus = rep(1:ncol(freq_matrix), each = nrow(freq_matrix)),
      Frequency = as.vector(freq_matrix),
      Migration_Rate = rate
    )
    plot_data <- rbind(plot_data, temp_data)
    
    # Add confidence intervals for selected loci
    selected_loci <- round(seq(1, ncol(freq_matrix), length.out = 5))
    for(locus in selected_loci) {
      ci_data <- rbind(ci_data, data.frame(
        Generation = 1:nrow(freq_matrix),
        Locus = paste("Locus", locus),
        Mean = result$summary$mean_frequencies[, locus],
        Lower = result$summary$ci_lower[, locus],
        Upper = result$summary$ci_upper[, locus],
        Migration_Rate = rate
      ))
    }
  }
  
  # Create heatmap plot
  p1 <- ggplot(plot_data, aes(x = Locus, y = Generation, fill = Frequency)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "white",
      mid = "skyblue",
      high = "darkblue",
      midpoint = 0.5,
      limits = c(0, 1)
    ) +
    facet_wrap(~Migration_Rate, ncol = 2) +
    labs(
      title = "Allele Frequency Evolution Across Migration Rates",
      x = "Locus",
      y = "Generation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  # Create trajectory plot
  p2 <- ggplot(ci_data, aes(x = Generation, y = Mean, color = Locus, group = interaction(Locus, Migration_Rate))) +
    geom_line() +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Locus), alpha = 0.2, color = NA) +
    facet_wrap(~Migration_Rate, ncol = 2) +
    labs(
      title = "Allele Frequency Trajectories with 95% CI",
      x = "Generation",
      y = "Allele Frequency"
    ) +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  # Combine plots
  combined_plot <- p1 / p2 +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = "Introgression Patterns Across Different Migration Rates",
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )
  
  # Save plots
  ggsave(
    filename = file.path(output_dir, "plots", "combined_migration_rates.pdf"),
    plot = combined_plot,
    width = 15,
    height = 20
  )
  
  return(list(heatmap = p1, trajectories = p2, combined = combined_plot))
}

#' Create an enhanced visualization of network metrics with interpretation zones
#' @param metrics_df Data frame containing network metrics
#' @param add_interpretation Logical; whether to add interpretation zones
plot_network_metrics_with_interpretation <- function(metrics_df, output_dir) {
  # Create output directory if it doesn't exist
  dir.create(file.path(output_dir, "figures", "network_metrics", "interpretation"), 
             recursive = TRUE, showWarnings = FALSE)
  
  # Calculate means and confidence intervals
  summary_stats <- metrics_df %>%
    pivot_longer(
      cols = c(mean_degree, transitivity, modularity),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(migration_rate, metric) %>%
    summarize(
      mean_val = mean(value),
      ci_lower = mean_val - qt(0.975, n() - 1) * sd(value) / sqrt(n()),
      ci_upper = mean_val + qt(0.975, n() - 1) * sd(value) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Create individual plots with enhanced styling
  plots <- list()
  
  # Mean Degree Plot
  plots$mean_degree <- ggplot(subset(summary_stats, metric == "mean_degree"), 
                              aes(x = migration_rate, y = mean_val)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "blue") +
    geom_line(color = "blue") +
    geom_point() +
    scale_x_log10() +
    labs(x = "Migration Rate (log scale)", 
         y = "Mean Degree",
         title = "Network Connectivity") +
    theme_minimal() +
    annotate("text", x = 1e-4, y = max(subset(summary_stats, metric == "mean_degree")$ci_upper),
             hjust = 0, vjust = 1,
             label = "Higher values indicate\nmore connected networks", 
             size = 3)
  
  # Transitivity Plot
  plots$transitivity <- ggplot(subset(summary_stats, metric == "transitivity"), 
                               aes(x = migration_rate, y = mean_val)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "red") +
    geom_line(color = "red") +
    geom_point() +
    scale_x_log10() +
    labs(x = "Migration Rate (log scale)", 
         y = "Transitivity",
         title = "Local Clustering") +
    theme_minimal() +
    annotate("text", x = 1e-4, y = max(subset(summary_stats, metric == "transitivity")$ci_upper),
             hjust = 0, vjust = 1,
             label = "Higher values indicate\nmore local clustering", 
             size = 3)
  
  # Modularity Plot
  plots$modularity <- ggplot(subset(summary_stats, metric == "modularity"), 
                             aes(x = migration_rate, y = mean_val)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "green4") +
    geom_line(color = "green4") +
    geom_point() +
    scale_x_log10() +
    labs(x = "Migration Rate (log scale)", 
         y = "Modularity",
         title = "Network Modularity") +
    theme_minimal() +
    annotate("text", x = 1e-4, y = max(subset(summary_stats, metric == "modularity")$ci_upper),
             hjust = 0, vjust = 1,
             label = "Higher values indicate\nmore distinct communities", 
             size = 3)
  
  # Combine plots
  combined_plot <- (plots$mean_degree / plots$transitivity / plots$modularity) +
    plot_layout(guides = 'collect') +
    plot_annotation(
      title = "Gene Regulatory Network Structure Changes with Migration Rate",
      subtitle = "Impact of introgression on network topology",
      caption = "Shaded areas represent 95% confidence intervals",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 10)
      )
    )
  
  # Save individual plots
  for (metric in names(plots)) {
    ggsave(
      filename = file.path(output_dir, "figures", "network_metrics", "interpretation", 
                           paste0(metric, "_plot.png")),
      plot = plots[[metric]],
      width = 8,
      height = 6
    )
  }
  
  # Save combined plot
  ggsave(
    filename = file.path(output_dir, "figures", "network_metrics", "interpretation", 
                         "combined_metrics.png"),
    plot = combined_plot,
    width = 10,
    height = 12
  )
  
  return(combined_plot)
}

create_interpretation_summary <- function(metrics_df, output_dir) {
  # Create output directory if it doesn't exist
  dir.create(file.path(output_dir, "results"), recursive = TRUE, showWarnings = FALSE)
  
  # Calculate overall trends
  trends <- metrics_df %>%
    pivot_longer(
      cols = c(mean_degree, transitivity, modularity),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(metric) %>%
    summarize(
      correlation = cor(log10(migration_rate), value),
      trend = case_when(
        correlation > 0.5 ~ "Strong positive",
        correlation > 0.3 ~ "Moderate positive",
        correlation > 0.1 ~ "Weak positive",
        correlation < -0.5 ~ "Strong negative",
        correlation < -0.3 ~ "Moderate negative",
        correlation < -0.1 ~ "Weak negative",
        TRUE ~ "No clear trend"
      )
    )
  
  # Save trends to CSV
  write.csv(
    trends,
    file = file.path(output_dir, "results", "network_metrics_trends.csv"),
    row.names = FALSE
  )
  
  return(trends)
}
#############################################
### Network Averaging Functions
#############################################

#' Create averaged correlation matrix from multiple replicates
#' @param results List of simulation results for a given migration rate
#' @return Average correlation matrix
calculate_average_correlation <- function(results) {
  # Extract correlation matrices from all replicates
  cor_matrices <- lapply(results$replicates, function(rep) {
    if(!is.null(rep$transcriptomes)) {
      return(cor(rep$transcriptomes$r))
    }
    return(NULL)
  })
  
  # Remove NULL entries
  cor_matrices <- cor_matrices[!sapply(cor_matrices, is.null)]
  
  # Calculate average correlation matrix
  if(length(cor_matrices) > 0) {
    avg_cor <- Reduce('+', cor_matrices) / length(cor_matrices)
    return(avg_cor)
  }
  return(NULL)
}

# Update these functions in 2-visualization_functions.R

#' #' Plot averaged correlation heatmap with hybridization style
#' #' @param cor_matrix Average correlation matrix
#' #' @param Chr.lengths Vector of chromosome lengths
#' #' @param migration_rate Migration rate
#' plot_averaged_heatmap <- function(cor_matrix, Chr.lengths, migration_rate) {
#'   # Set margins
#'   par(mar = c(5, 5, 4, 2))
#'   
#'   # Create heatmap using base R to match hybridization style
#'   heatmap(cor_matrix,
#'           Rowv = NA,  # No row dendrogram
#'           Colv = NA,  # No column dendrogram
#'           col = colorRampPalette(c("#FFFFD4", "#FED98E", "#FE9929", "#CC4C02"))(100),  # Yellow to brown palette
#'           main = sprintf("Generation %d Average Correlation", migration_rate),
#'           labRow = paste0("X", 1:nrow(cor_matrix)),
#'           labCol = paste0("X", 1:ncol(cor_matrix)),
#'           scale = "none",  # No scaling
#'           margins = c(5, 5))  # Adjust margins as needed
#' }
#' 
#' #' Plot averaged gene network with hybridization style
#' #' @param cor_matrix Average correlation matrix
#' #' @param Chr.lengths Vector of chromosome lengths
#' #' @param threshold Correlation threshold
#' #' @return Plot of gene network
#' plot_averaged_network <- function(cor_matrix, Chr.lengths, threshold = 0.2) {
#'   # Create adjacency matrix
#'   adjacency <- (abs(cor_matrix) > threshold) * 1
#'   diag(adjacency) <- 0
#'   
#'   # Create igraph object
#'   g <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
#'   
#'   # Calculate node colors based on chromosomes
#'   node_colors <- rep(NA, vcount(g))
#'   chr_start <- 1
#'   for(i in 1:length(Chr.lengths)) {
#'     node_colors[chr_start:(chr_start + Chr.lengths[i] - 1)] <- i
#'     chr_start <- chr_start + Chr.lengths[i]
#'   }
#'   
#'   # Set margins for plot
#'   par(mar = c(1, 1, 1, 1))
#'   
#'   # Plot network with hybridization style
#'   plot(g,
#'        layout = layout_in_circle(g),
#'        vertex.size = 6,  # Match hybridization size
#'        vertex.label = NA,  # No labels
#'        vertex.color = node_colors,  # Color by chromosome
#'        edge.arrow.size = 0.2,
#'        edge.width = 2,
#'        edge.color = "gray70")  # Gray edges to match hybridization style
#' }
#' 


#' Plot averaged gene network
#' @param cor_matrix Average correlation matrix
#' @param Chr.lengths Vector of chromosome lengths
#' @param threshold Correlation threshold
#' @return Plot of gene network
plot_averaged_network <- function(cor_matrix, Chr.lengths, threshold = 0.2) {
  # Create adjacency matrix
  adjacency <- (abs(cor_matrix) > threshold) * 1
  diag(adjacency) <- 0

  # Create igraph object
  g <- graph_from_adjacency_matrix(adjacency, mode = "undirected")

  # Calculate node colors based on chromosomes
  ## Calculate node colors using the simpler method from plot.gene.net
  Chr.col <- c()
  for(i in 1:length(Chr.lengths)) {
    Chr.col <- c(Chr.col, rep(i, Chr.lengths[i]))
  }
  V(g)$color <- Chr.col
  
  # Set other properties
  V(g)$size <- 6
  V(g)$label <- NA
  E(g)$arrow.size <- .2
  E(g)$width <- 2
  E(g)$color <- "gray70"
  
  # Create plot
  plot(g,
       layout = layout_in_circle(g),
       vertex.size = 6,
       vertex.label = NA,
       vertex.color = Chr.col,
       edge.width = 2,
       edge.color = "gray70")

  # # Add legend for chromosomes
  # legend("topright",
  #        legend = paste("Chr", 1:length(Chr.lengths)),
  #        fill = Chr.col,
  #        title = "Chromosome",
  #        cex = 0.8)

  return(g)
}

#' Plot averaged correlation heatmap
#' @param cor_matrix Average correlation matrix
#' @param Chr.lengths Vector of chromosome lengths
#' @param migration_rate Migration rate
plot_averaged_heatmap <- function(cor_matrix, Chr.lengths, migration_rate) {
  # Add chromosome boundaries
  chr_ends <- cumsum(Chr.lengths)

  # Create heatmap
  heatmap(cor_matrix,
          Rowv = NA,
          Colv = NA,
          col = colorRampPalette(c("#FFFFD4", "#FED98E", "#FE9929", "#CC4C02"))(100)#,
          #main = sprintf("Average Gene Expression Correlation\nMigration Rate: %.6f", migration_rate)
          )

  # Add chromosome boundary lines
  for(pos in chr_ends) {
    abline(h = pos + 0.5, col = "black", lwd = 2)
    abline(v = pos + 0.5, col = "black", lwd = 2)
  }
}

#' Generate network visualizations for all migration rates
#' @param results_list List of simulation results
#' @param Chr.lengths Vector of chromosome lengths
#' @param threshold Correlation threshold
#' @param output_dir Output directory
generate_network_visualizations <- function(results_list, Chr.lengths, threshold = 0.2,
                                            output_dir = "results/network_analysis") {
  # Create directories if they don't exist
  dir.create(file.path(output_dir, "averaged_heatmaps"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "averaged_networks"), recursive = TRUE, showWarnings = FALSE)
  
  # Process each migration rate
  for(rate in names(results_list)) {
    cat("\nProcessing migration rate:", rate, "\n")
    
    # Calculate average correlation matrix
    avg_cor <- calculate_average_correlation(results_list[[rate]])
    
    if(!is.null(avg_cor)) {
      # Generate and save heatmap
      png(file.path(output_dir, "averaged_heatmaps", 
                    sprintf("avg_heatmap_mig%.6f.png", as.numeric(rate))),
          width = 800, height = 800, res = 100)
      plot_averaged_heatmap(avg_cor, Chr.lengths, as.numeric(rate))
      dev.off()
      
      # Generate and save network plot
      png(file.path(output_dir, "averaged_networks", 
                    sprintf("avg_network_mig%.6f.png", as.numeric(rate))),
          width = 800, height = 800, res = 100)
      par(mar = c(1,1,3,1))
      plot_averaged_network(avg_cor, Chr.lengths, threshold)
      dev.off()
    }
  }
}