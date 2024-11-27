# 1.Basic correlation calculation used by both simulations #####
coexpr <- function(genoreads, threshold) {
  cormat <- cor(genoreads$r)
  if(!is.na(threshold)) {
    cormat[abs(cormat) < threshold] <- 0
    cormat[abs(cormat) > threshold] <- 1
  }
  return(cormat)
}


#' Calculate summary statistics for network metrics
#' @param results Data frame of network metrics
#' @param group_var Variable to group by (generation or migration_rate)
calculate_summary_stats <- function(results, group_var) {
  results %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      mean_degree_mean = mean(mean_degree),
      mean_degree_se = sd(mean_degree)/sqrt(n()),
      modularity_mean = mean(modularity),
      modularity_se = sd(modularity)/sqrt(n()),
      transitivity_mean = mean(transitivity),
      transitivity_se = sd(transitivity)/sqrt(n()),
      .groups = 'drop'
    )
}

#' Save simulation results and parameters
#' @param results List of results
#' @param params List of parameters
#' @param output_file Output file path
save_simulation_results <- function(results, params, output_file) {
  saveRDS(
    list(
      results = results,
      parameters = params
    ),
    output_file
  )
}

#' Create metric visualization with confidence intervals
#' @param plot_data Data frame containing metric data
#' @param metric_name Name of metric to plot
#' @param color Color for plot elements
#' @param x_var Name of x-axis variable
#' @param title Plot title
#' @return ggplot object
create_metric_plot <- function(plot_data, metric_name, color,
                               x_var = "generation", title = NULL) {
  # Ensure required columns exist
  required_cols <- c("mean_val", "ci_lower", "ci_upper")
  if(!all(required_cols %in% names(plot_data))) {
    stop("plot_data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = mean_val)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.2, fill = color) +
    geom_line(color = color, size = 1) +
    geom_point(size = 2) +
    labs(
      title = title %||% stringr::str_to_title(metric_name),
      x = stringr::str_to_title(x_var),
      y = stringr::str_to_title(metric_name)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  return(p)
}

#' Plot averaged network heatmap
#' @param cor_matrix Correlation matrix
#' @param Chr.lengths Vector of chromosome lengths
#' @param title Optional plot title
#' @return NULL (creates plot)
plot_averaged_heatmap <- function(cor_matrix, Chr.lengths, title = NULL) {
  # Add chromosome boundaries
  chr_ends <- cumsum(Chr.lengths)
  
  # Create custom color palette
  heatmap_colors <- colorRampPalette(c(
    "#FFFFD4", # light yellow
    "#FED98E", # medium yellow
    "#FE9929", # orange
    "#CC4C02"  # dark orange/red
  ))(100)
  
  # Create heatmap
  heatmap(cor_matrix,
          Rowv = NA,
          Colv = NA,
          col = heatmap_colors,
          main = title)
  
  # Add chromosome boundary lines
  for(pos in chr_ends) {
    abline(h = pos + 0.5, col = "blue", lwd = 1)
    abline(v = pos + 0.5, col = "blue", lwd = 1)
  }
}
