# 1.Basic correlation calculation used by both simulations #####
coexpr <- function(genoreads, threshold) {
  cormat <- cor(genoreads$r)
  if(!is.na(threshold)) {
    cormat[abs(cormat) < threshold] <- 0
    cormat[abs(cormat) > threshold] <- 1
  }
  return(cormat)
}

# 2. Create averaged correlation matrix from multiple replicates #######
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