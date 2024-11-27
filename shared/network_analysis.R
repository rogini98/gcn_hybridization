# shared network functions
# script contains all functions handling network 
# creation, analysis and metrics

# 1. Network architecture #######################
#' Create network effect matrix with specified architecture
#' @param N_genes Number of genes in network
#' @param network_density Probability of edge presence
#' @param network_effect Standard deviation of edge effects
#' @param type Network type: "random" (default) or "BarabasiAlbert"
#' @param power_law_exp Power law exponent for BA networks (default: 1)
#' @return Matrix of network effects
get_effect_matrix <- function(N_genes, 
                              network_density, 
                              network_effect,
                              type = c("random", "BarabasiAlbert"),
                              power_law_exp = 1) {
  
  # Input validation
  type <- match.arg(type)
  if(N_genes < 1) stop("N_genes must be positive")
  if(network_density < 0 || network_density > 1) stop("network_density must be between 0 and 1")
  if(network_effect < 0) stop("network_effect must be non-negative")
  
  # Generate base incidence matrix
  incidences <- rbinom(n = N_genes^2, size = 1, prob = network_density)
  
  # Generate effect sizes
  effects <- rnorm(n = N_genes^2, mean = 0, sd = network_effect)
  
  if(type == "random") {
    # Simple random network
    effect_matrix <- matrix(incidences * effects, 
                            nrow = N_genes, 
                            ncol = N_genes)
  } else {
    # Barabási-Albert preferential attachment network
    ba_network <- igraph::sample_pa(n = N_genes, 
                                    power = power_law_exp, 
                                    m = 1, 
                                    directed = FALSE)
    ba_matrix <- as.matrix(igraph::as_adj(ba_network))
    
    # Combine BA structure with random effects
    effect_matrix <- matrix(incidences * effects, nrow = N_genes) * ba_matrix
  }
  
  # Ensure matrix is symmetric (undirected network)
  effect_matrix[lower.tri(effect_matrix)] <- t(effect_matrix)[lower.tri(effect_matrix)]
  
  return(effect_matrix)
}


# 2. Calculate network metrics for simulation output ###########
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


# 3. Population transcriptomes ####################
#' Get transcriptomes for a population
#' @param Population List of individual genotypes
#' @param N_genes Number of genes
#' @param intercepts Expression intercepts
#' @param betas Genetic effect sizes
#' @param effect_matrices List of network effect matrices
#' @param neteQTLlocus Network QTL locus
#' @return List containing genotypes and readcounts
getTranscriptomes <- function(Population, N_genes, intercepts, betas, effect_matrices, neteQTLlocus) {
  PopSize <- length(Population)
  genotypes <- data.frame(matrix(NA, nrow = PopSize, ncol = N_genes))
  readcounts <- data.frame(matrix(NA, nrow = PopSize, ncol = N_genes))
  
  for(i in 1:PopSize) {
    ind1 <- Population[[i]]
    NChr <- length(ind1)
    genome_val <- c()
    for(j in seq(1, NChr-1, by = 2)) {
      genome_val <- c(genome_val, ind1[[j]] + ind1[[j+1]])
    }
    genotypes[i,] <- genome_val
    readcounts[i,] <- genotype_to_netexpr(ind1, intercepts, betas, effect_matrices, neteQTLlocus)
  }
  
  return(list(g = genotypes, r = readcounts))
}

# 4. Network expression ########################
#' Generate transcriptome with network effects
#' @param genotype Individual genotype
#' @param intercepts Expression intercepts
#' @param betas Genetic effect sizes
#' @param effect_matrices List of network effect matrices
#' @param neteQTLlocus Network QTL locus
#' @return Expression values
genotype_to_netexpr <- function(genotype, intercepts, betas, effect_matrices, neteQTLlocus) {
  NChr <- length(genotype)
  genome_val <- c()
  for(i in seq(1, NChr-1, by = 2)) {
    genome_val <- c(genome_val, genotype[[i]] + genotype[[i+1]])
  }
  
  lambda <- exp(intercepts + genome_val*betas)
  null_expression <- rpois(n = length(lambda), lambda)
  N_genes <- length(null_expression)
  effect.matrix <- effect_matrices[[genome_val[neteQTLlocus] + 1]]
  
  for(i in 1:(N_genes-1)) {
    for(j in i:N_genes) {
      null_expression[j] <- null_expression[j] + 
        effect.matrix[i,j] * (null_expression[i] - lambda[i])
    }
  }
  
  return(null_expression)
}