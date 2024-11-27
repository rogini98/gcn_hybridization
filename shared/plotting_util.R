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
