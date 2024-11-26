# Script to generate plots for network metrics evolution and relationships

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(igraph)

# Load results if not in environment
if(!exists("results")) {
  load("simulation_results.RData")
}

# Calculate parent metrics if not already done
if(!exists("parent_metrics")) {
  parent_metrics <- list(
    AA = {
      g_AA <- graph_from_adjacency_matrix(effect.matrix00 > threshold, mode="undirected", diag=FALSE)
      comm_AA <- cluster_fast_greedy(g_AA)
      list(
        mean_degree = mean(degree(g_AA)),
        modularity = modularity(g_AA, membership(comm_AA)),
        transitivity = transitivity(g_AA)
      )
    },
    BB = {
      g_BB <- graph_from_adjacency_matrix(effect.matrix11 > threshold, mode="undirected", diag=FALSE)
      comm_BB <- cluster_fast_greedy(g_BB)
      list(
        mean_degree = mean(degree(g_BB)),
        modularity = modularity(g_BB, membership(comm_BB)),
        transitivity = transitivity(g_BB)
      )
    }
  )
}

###########################################
# 1. Create Evolution Plots
###########################################

# Function to create evolution plots with parent references
create_evolution_plots <- function(summary_stats, parent_metrics) {
  # Mean Degree Plot
  p1 <- ggplot(summary_stats, aes(x = generation)) +
    geom_ribbon(aes(ymin = mean_degree_mean - 1.96*mean_degree_se, 
                    ymax = mean_degree_mean + 1.96*mean_degree_se),
                alpha = 0.2, fill = "#2C3E50") +
    geom_line(aes(y = mean_degree_mean), linewidth = 0.8, color = "#2C3E50") +
    geom_point(aes(y = mean_degree_mean), size = 2, color = "#2C3E50") +
    geom_hline(yintercept = parent_metrics$AA$mean_degree, 
               linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = parent_metrics$BB$mean_degree, 
               linetype = "dashed", color = "blue", linewidth = 0.8) +
    annotate("text", x = -1, y = parent_metrics$AA$mean_degree, 
             label = "AA", color = "red", hjust = 1, vjust = -0.5) +
    annotate("text", x = -1, y = parent_metrics$BB$mean_degree, 
             label = "BB", color = "blue", hjust = 1, vjust = -0.5) +
    labs(title = "Mean Degree",
         x = "Generation",
         y = "Mean Degree") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  # Modularity Plot
  p2 <- ggplot(summary_stats, aes(x = generation)) +
    geom_ribbon(aes(ymin = modularity_mean - 1.96*modularity_se, 
                    ymax = modularity_mean + 1.96*modularity_se),
                alpha = 0.2, fill = "#2C3E50") +
    geom_line(aes(y = modularity_mean), linewidth = 0.8, color = "#2C3E50") +
    geom_point(aes(y = modularity_mean), size = 2, color = "#2C3E50") +
    geom_hline(yintercept = parent_metrics$AA$modularity, 
               linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = parent_metrics$BB$modularity, 
               linetype = "dashed", color = "blue", linewidth = 0.8) +
    annotate("text", x = -1, y = parent_metrics$AA$modularity, 
             label = "AA", color = "red", hjust = 1, vjust = -0.5) +
    annotate("text", x = -1, y = parent_metrics$BB$modularity, 
             label = "BB", color = "blue", hjust = 1, vjust = -0.5) +
    labs(title = "Modularity",
         x = "Generation",
         y = "Modularity") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  # Transitivity Plot
  p3 <- ggplot(summary_stats, aes(x = generation)) +
    geom_ribbon(aes(ymin = transitivity_mean - 1.96*transitivity_se, 
                    ymax = transitivity_mean + 1.96*transitivity_se),
                alpha = 0.2, fill = "#2C3E50") +
    geom_line(aes(y = transitivity_mean), linewidth = 0.8, color = "#2C3E50") +
    geom_point(aes(y = transitivity_mean), size = 2, color = "#2C3E50") +
    geom_hline(yintercept = parent_metrics$AA$transitivity, 
               linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = parent_metrics$BB$transitivity, 
               linetype = "dashed", color = "blue", linewidth = 0.8) +
    annotate("text", x = -1, y = parent_metrics$AA$transitivity, 
             label = "AA", color = "red", hjust = 1, vjust = -0.5) +
    annotate("text", x = -1, y = parent_metrics$BB$transitivity, 
             label = "BB", color = "blue", hjust = 1, vjust = -0.5) +
    labs(title = "Transitivity",
         x = "Generation",
         y = "Transitivity") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  # Comparative plot
  metrics_long <- summary_stats %>%
    select(generation, 
           mean_degree = mean_degree_mean, 
           modularity = modularity_mean, 
           transitivity = transitivity_mean) %>%
    pivot_longer(
      -generation,
      names_to = "metric",
      values_to = "value"
    )
  
  p4 <- ggplot(metrics_long, aes(x = generation, y = value, color = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#E74C3C", "#3498DB", "#2ECC71")) +
    labs(title = "Comparative Network Metrics",
         x = "Generation",
         y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80"),
      legend.position = "bottom"
    )
  
  # Combine all plots
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Network Metrics Evolution",
      subtitle = "Dashed lines indicate parental genotype values (AA: red, BB: blue)",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40")
      )
    )
  
  return(combined_plot)
}

###########################################
# 2. Create Relationship Plots
###########################################

# Function to create relationship plots
create_relationship_plots <- function(results) {
  # Calculate correlations
  correlations <- list(
    deg_mod = cor(results$mean_degree, results$modularity),
    deg_trans = cor(results$mean_degree, results$transitivity),
    mod_trans = cor(results$modularity, results$transitivity)
  )
  
  # Mean Degree vs Modularity
  p1 <- ggplot(results, aes(x = mean_degree, y = modularity)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "red") +
    annotate("text", x = min(results$mean_degree), y = max(results$modularity),
             label = sprintf("r = %.3f", correlations$deg_mod),
             hjust = 0, vjust = 1) +
    labs(title = "Mean Degree vs Modularity",
         x = "Mean_Degree",
         y = "Modularity") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  # Mean Degree vs Transitivity
  p2 <- ggplot(results, aes(x = mean_degree, y = transitivity)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "red") +
    annotate("text", x = min(results$mean_degree), y = max(results$transitivity),
             label = sprintf("r = %.3f", correlations$deg_trans),
             hjust = 0, vjust = 1) +
    labs(title = "Mean Degree vs Transitivity",
         x = "Mean_Degree",
         y = "Transitivity") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  # Modularity vs Transitivity
  p3 <- ggplot(results, aes(x = modularity, y = transitivity)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "red") +
    annotate("text", x = min(results$modularity), y = max(results$transitivity),
             label = sprintf("r = %.3f", correlations$mod_trans),
             hjust = 0, vjust = 1) +
    labs(title = "Modularity vs Transitivity",
         x = "Modularity",
         y = "Transitivity") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  # Combined plot
  relationship_plot <- (p1 + p2 + p3) +
    plot_annotation(
      title = "Network Metrics Relationships",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold")
      )
    )
  
  return(relationship_plot)
}

# Generate and save plots
evolution_plot <- create_evolution_plots(summary_stats, parent_metrics)
relationship_plot <- create_relationship_plots(results)

# Save plots
ggsave("network_evolution.pdf", evolution_plot, width = 12, height = 10)
ggsave("network_relationships.pdf", relationship_plot, width = 15, height = 5)

# Print correlation matrix
cor_matrix <- cor(results[c("mean_degree", "modularity", "transitivity")])
print("Correlation Matrix:")
print(round(cor_matrix, 3))