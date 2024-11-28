#########################################
# Master script for network evolution analysis (Scenario 1)
# This script runs the entire workflow:
# 1. Loads functions
# 2. Runs multiple replicates of simulations in parallel
# 3. Creates visualizations including averaged heatmaps and networks
#########################################

# Clear workspace and set random seed for reproducibility
rm(list = ls())
set.seed(123)

# Check and install required packages if needed
required_packages <- c("igraph", "ggplot2", "dplyr", "tidyr", 
                       "patchwork", "parallel", "doParallel", "foreach",
                       "pheatmap", "RColorBrewer")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Create output directory structure
dirs <- c("output-hybridization",
          "output-hybridization/plots",
          "output-hybridization/plots/heatmaps",
          "output-hybridization/plots/networks")
for(dir in dirs) {
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Start logging
log_file <- file("output-hybridization/analysis_log.txt", open = "w")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("Starting analysis:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

#########################################
# 1. Source functions and set parameters
#########################################
cat("Loading functions...\n")
source("./hybridization/01-network_functions_hybridizationSim.R")

# Define plot.gene.net function
plot.gene.net <- function(cormat, Chr.lengths, layout = "circle"){
  # Ensure matrix symmetry for undirected graph
  adj_matrix <- as.matrix(as.dist(cormat))
  # Make symmetric by taking maximum of corresponding elements
  adj_matrix[upper.tri(adj_matrix)] <- t(adj_matrix)[upper.tri(adj_matrix)]
  
  gr <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  V(gr)$size <- 6
  Chr.col <- c()
  for(i in 1:length(Chr.lengths)){
    Chr.col <- c(Chr.col, rep(i, Chr.lengths[i]))
  }
  V(gr)$color <- Chr.col
  V(gr)$label <- NA
  E(gr)$arrow.size <- .2
  E(gr)$width <- 2
  if(layout == "net"){
    l = layout_with_fr(gr)
  }else(if(layout == "circle"){l = layout_in_circle(gr)}else(l = layout_randomly(gr)))
  plot(gr, layout=l)
}

# Set up initial parameters
Chr.lengths <- c(10, 10, 5, 5, 10)
N_genes <- sum(Chr.lengths)
PopSize <- 500
n_replicates <- 10
threshold <- 0.2
generations <- 0:20

# Generate transcriptome parameters
set.seed(123)
meanexpr <- 5
sd.expr <- 2
sd.beta <- 3
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

#########################################
# 2. Run single replicate function
#########################################

run_single_replicate <- function(rep_num) {
  cat(sprintf("Running replicate %d...\n", rep_num))
  
  # Store results for each generation
  gen_results <- list()
  
  for(gen in generations) {
    # Run simulation
    sim_results <- netQTLsim(
      Fgen = gen,
      Chr.lengths = Chr.lengths,
      PopSize = PopSize,
      intercepts = intercepts,
      betas = betas,
      effect.matrices = effect.matrices,
      neteQTLlocus = neteQTLlocus
    )
    
    # Calculate correlation matrix
    cor_matrix <- cor(sim_results$genoreads$r)
    
    # Calculate network metrics
    adj_matrix <- cor_matrix > threshold
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    comm <- cluster_fast_greedy(g)
    
    # Store results
    gen_results[[gen + 1]] <- list(
      generation = gen,
      replicate = rep_num,
      cor_matrix = cor_matrix,
      mean_degree = mean(degree(g)),
      modularity = modularity(g, membership(comm)),
      transitivity = transitivity(g, type = "global")
    )
  }
  
  return(gen_results)
}

#########################################
# 3. Run parallel simulations
#########################################
cat("\nStarting parallel simulations...\n")

# Set up parallel processing
n_cores <- min(detectCores() - 1, n_replicates)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary objects to the cluster
clusterExport(cl, c("Chr.lengths", "PopSize", "intercepts", "betas", 
                    "effect.matrices", "neteQTLlocus", "threshold",
                    "generations", "netQTLsim", "N_genes"))

# Run parallel simulations
all_results <- foreach(rep = 1:n_replicates, 
                       .packages = c("igraph"),
                       .errorhandling = "pass") %dopar% {
                         run_single_replicate(rep)
                       }

stopCluster(cl)

#########################################
# 4. Process results and create visualizations
#########################################
cat("\nProcessing results and creating visualizations...\n")

# Initialize storage for averaged results
avg_cor_matrices <- vector("list", length(generations))
metrics_df <- data.frame(
  generation = numeric(),
  replicate = numeric(),
  mean_degree = numeric(),
  modularity = numeric(),
  transitivity = numeric()
)

# Process results
for(gen in generations) {
  gen_idx <- gen + 1
  
  # Initialize average correlation matrix
  avg_cor_matrices[[gen_idx]] <- matrix(0, N_genes, N_genes)
  
  # Sum up results across replicates
  for(rep in 1:n_replicates) {
    rep_results <- all_results[[rep]][[gen_idx]]
    
    # Add to average correlation matrix
    avg_cor_matrices[[gen_idx]] <- avg_cor_matrices[[gen_idx]] + rep_results$cor_matrix
    
    # Store metrics in data frame
    metrics_df <- rbind(metrics_df, data.frame(
      generation = gen,
      replicate = rep,
      mean_degree = rep_results$mean_degree,
      modularity = rep_results$modularity,
      transitivity = rep_results$transitivity
    ))
  }
  
  # Calculate final average correlation matrix
  avg_cor_matrices[[gen_idx]] <- avg_cor_matrices[[gen_idx]] / n_replicates

  # Create and save averaged heatmap
  png(sprintf("output-hybridization/plots/heatmaps/heatmap_gen_%d.png", gen),
      width = 800, height = 800, res = 100)
  
  # Set margins
  par(mar = c(5, 5, 4, 2))
  
  # Create heatmap using base R heatmap function
  heatmap(avg_cor_matrices[[gen_idx]], 
          Rowv = NA,  # No row dendrogram
          Colv = NA,  # No column dendrogram
          scale = "none",  # No scaling
          margins = c(5,5),  # Adjust margins
          #col = colorRampPalette(c("#FFFFD4", "#FED98E", "#FE9929", "#CC4C02"))(100),
          # main = sprintf("Generation %d Average Correlation", gen),
          # labRow = paste0("X", 1:N_genes),
          # labCol = paste0("X", 1:N_genes)
          )
  dev.off()
  
  # Create and save averaged network plot
  png(sprintf("output-hybridization/plots/networks/network_gen_%d.png", gen),
      width = 800, height = 800, res = 100)
  par(mar = c(1, 1, 1, 1))
  # Use threshold to create binary adjacency matrix
  adj_matrix <- avg_cor_matrices[[gen_idx]] > threshold
  plot.gene.net(adj_matrix, Chr.lengths, layout = "circle")
  dev.off()
}

# Calculate summary statistics
summary_stats <- metrics_df %>%
  group_by(generation) %>%
  summarise(
    mean_degree_mean = mean(mean_degree),
    mean_degree_se = sd(mean_degree)/sqrt(n()),
    modularity_mean = mean(modularity),
    modularity_se = sd(modularity)/sqrt(n()),
    transitivity_mean = mean(transitivity),
    transitivity_se = sd(transitivity)/sqrt(n())
  )

# Calculate correlations between metrics
cor_matrix <- cor(metrics_df[c("mean_degree", "modularity", "transitivity")])

# Print correlation matrix
cat("\nCorrelation Matrix:\n")
print(round(cor_matrix, 3))

#########################################
# 5. Save results
#########################################
results <- metrics_df  # Convert the metrics dataframe to match expected format

# Save all results including the reformatted 'results' variable
save(results, metrics_df, summary_stats, avg_cor_matrices, cor_matrix,
     file = "output-hybridization/simulation_results.RData")

# Save all results
# saveRDS(list(metrics = metrics_df, 
#              summary_stats = summary_stats,
#              avg_cor_matrices = avg_cor_matrices,
#              cor_matrix = cor_matrix),
#         "output-hybridization/simulation_results.rds")

#########################################
# 6. create final plots
#########################################

# Create final summary plots
source("./hybridization/03-generate_hybridization_plots.R")

# End logging
cat("\nAnalysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink(type = "output")
sink(type = "message")
close(log_file)

# Print success message to console
cat("\nAnalysis complete! Results saved in output directory.\n")
cat("See 'output/analysis_log.txt' for detailed log.\n")