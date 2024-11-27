# Script to run multiple replicates of hybridization simulation
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)
source("1-network_functions_hybridizationSim.R")

# Set up initial parameters
Chr.lengths <- c(10, 10, 5, 5, 10)
N_genes <- sum(Chr.lengths)
meanexpr <- 5
sd.expr <- 2
sd.beta <- 3
threshold <- 0.2

# Generate transcriptome parameters
set.seed(123)
intercepts <- rnorm(N_genes, meanexpr, sd.expr)
betas <- rnorm(N_genes, mean = 0, sd.beta)

# Network parameters
network_density00 <- 0.03
network_effect00 <- 2
network_density11 <- 0.01
network_effect11 <- 2.5

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
