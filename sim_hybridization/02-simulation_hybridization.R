# Script to run multiple replicates of hybridization simulation
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)
source("./sim_hybridization/01-functions_hybridization.R")

# Set up initial parameters to match documentation
Chr.lengths <- c(10, 10, 5, 5, 10)
N_genes <- sum(Chr.lengths)
meanexpr <- 5          # μ_α = 5 from text
sd.expr <- 2           # σ_α = 2 from text
sd.beta <- 2          # Changed from 3 to 2 to match σ_β = 2 from text
threshold <- 0.2

# Generate transcriptome parameters
set.seed(123)
intercepts <- rnorm(N_genes, meanexpr, sd.expr)
betas <- rnorm(N_genes, mean = 0, sd.beta)

# Network parameters updated to match text
network_density00 <- 0.2    # Changed from 0.03 to 0.2 to match text
network_effect00 <- 2       # σ_γ = 2 from text
network_density11 <- 0.1    # Changed from 0.01 to 0.1 to match text
network_effect11 <- 2.5     # From text


# Create network architecture
neteQTLlocus <- sample(1:N_genes, 1)
effect.matrix00 <- get_effect.matrix(N_genes, network_density00, network_effect00)
effect.matrix11 <- get_effect.matrix(N_genes, network_density11, network_effect11)
effect.matrix01.add <- (effect.matrix00 + effect.matrix11) / 2
effect.matrices <- list(AA = effect.matrix00, AB = effect.matrix01.add, BB = effect.matrix11)

# Function to run a single replicate
run_single_replicate <- function(rep_num) {
  cat(sprintf("Running replicate %d...\n", rep_num))
  
  metrics <- data.frame(
    replicate = rep_num,
    generation = 0:20,
    mean_degree = NA,
    modularity = NA,
    transitivity = NA
  )
  
  for(gen in 0:20) {
    sim_results <- netQTLsim(
      Fgen = gen,
      Chr.lengths = Chr.lengths,
      PopSize = 500,
      intercepts = intercepts,
      betas = betas,
      effect.matrices = effect.matrices,
      neteQTLlocus = neteQTLlocus
    )
    
    # Calculate network metrics
    adj_matrix <- sim_results$co_expr > threshold
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    comm <- cluster_fast_greedy(g)
    
    metrics$mean_degree[gen + 1] <- mean(degree(g))
    metrics$modularity[gen + 1] <- modularity(g, membership(comm))
    metrics$transitivity[gen + 1] <- transitivity(g, type = "global")
  }
  
  return(metrics)
}

# Set up parallel processing
n_replicates <- 10
n_cores <- min(detectCores() - 1, n_replicates)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary objects to the cluster
clusterExport(cl, c("Chr.lengths", "intercepts", "betas", 
                    "effect.matrices", "neteQTLlocus", "threshold",
                    "netQTLsim"))

# Run parallel simulations
results <- foreach(rep = 1:n_replicates, 
                   .combine = rbind,
                   .packages = c("igraph")) %dopar% {
                     run_single_replicate(rep)
                   }

# Stop cluster
stopCluster(cl)

# Calculate summary statistics
summary_stats <- results %>%
  group_by(generation) %>%
  summarise(
    mean_degree_mean = mean(mean_degree),
    mean_degree_se = sd(mean_degree)/sqrt(n()),
    modularity_mean = mean(modularity),
    modularity_se = sd(modularity)/sqrt(n()),
    transitivity_mean = mean(transitivity),
    transitivity_se = sd(transitivity)/sqrt(n())
  )

# Save results
saveRDS(list(
  results = results,
  summary_stats = summary_stats,
  parameters = list(
    Chr.lengths = Chr.lengths,
    network_params = list(
      density00 = network_density00,
      effect00 = network_effect00,
      density11 = network_density11,
      effect11 = network_effect11
    ),
    n_replicates = n_replicates,
    threshold = threshold
  )
), "simulation_results.rds")
