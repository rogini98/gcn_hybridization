# Core introgression simulation functions
# Authors: Original by Dan Bolnick, enhanced & reorganized by Rogini Runghen
# Last updated: November 25, 2024

#############################################
### Core Functions for introgression simulation
#############################################

# 1. Calculate allele frequencies #################
#' @param population List of individual genotypes
#' @return Vector of allele frequencies
get_allele_freq <- function(population) {
  first_ind <- population[[1]]
  chr_lengths <- sapply(seq(1, length(first_ind), by=2), function(i) {
    length(first_ind[[i]])
  })
  
  allelecounts <- matrix(0, nrow = 1, ncol = sum(chr_lengths))
  
  for(i in 1:length(population)) {
    genome <- c()
    for(chr in seq(1, 2*length(chr_lengths), by = 2)) {
      Chromatid1 <- population[[i]][[chr]]
      Chromatid2 <- population[[i]][[chr + 1]]
      Genotypes <- Chromatid1 + Chromatid2
      genome <- c(genome, Genotypes)
    }
    allelecounts <- allelecounts + genome
  }
  
  return(allelecounts / (2*length(population)))
}

#############################################
### Expression and Network Functions
#############################################

#' Generate transcriptome genetic architecture
#' @param N_genes Number of genes
#' @param meanexpr Mean expression level
#' @param sd.expr Standard deviation of expression
#' @param sd.beta Standard deviation of genetic effects
#' @return List of transcriptome parameters
generate_transcriptome_parameters <- function(N_genes, meanexpr = 5, sd.expr = 2, sd.beta = 3) {
  intercepts <- rnorm(N_genes, meanexpr, sd.expr)
  betas <- rnorm(N_genes, mean = 0, sd.beta)
  neteQTLlocus <- sample(1:N_genes, 1)
  
  return(list(
    intercepts = intercepts,
    betas = betas,
    neteQTLlocus = neteQTLlocus
  ))
}






#############################################
### Network Analysis Functions
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


#############################################
### Parameter Initialization
#############################################

#' Initialize all simulation parameters
#' @return List of all parameters needed for simulation
initialize_parameters <- function(Chr.lengths, PopSize, Generations,
                                  network_params, expression_params) {
  # Basic parameters validation
  stopifnot(
    is.numeric(Chr.lengths),
    is.numeric(PopSize),
    is.numeric(Generations),
    PopSize > 0,
    Generations > 0,
    all(Chr.lengths > 0)
  )
  
  N_genes <- sum(Chr.lengths)
  
  # Generate transcriptome parameters
  trans_params <- generate_transcriptome_parameters(
    N_genes = N_genes,
    meanexpr = expression_params$mean,
    sd.expr = expression_params$sd,
    sd.beta = expression_params$beta_sd
  )
  
  # Generate network architectures
  effect.matrix00 <- get_effect_matrix(
    N_genes,
    network_params$density00,
    network_params$effect00
  )
  
  effect.matrix11 <- get_effect_matrix(
    N_genes,
    network_params$density11,
    network_params$effect11
  )
  
  effect.matrix01.add <- (effect.matrix00 + effect.matrix11) / 2
  
  return(list(
    Generations = Generations,
    PopSize = PopSize,
    Chr.lengths = Chr.lengths,
    intercepts = trans_params$intercepts,
    betas = trans_params$betas,
    neteQTLlocus = trans_params$neteQTLlocus,
    effect_matrices = list(effect.matrix00, effect.matrix01.add, effect.matrix11),
    network_params = network_params
  ))
}

########################################
# Simulation runner functions
########################################
#' Run introgression simulation with tracking
#' @param Generations Number of generations
#' @param PopSize Population size
#' @param Chr.lengths Vector of chromosome lengths
#' @param migration_rate Rate of migration
#' @param params List of simulation parameters
#' @return List containing simulation results
introgressing_with_tracking <- function(Generations, PopSize, Chr.lengths, migration_rate, params) {
  # Input validation
  if (!is.numeric(Generations) || Generations < 1) stop("Generations must be a positive integer")
  if (!is.numeric(PopSize) || PopSize < 1) stop("PopSize must be a positive integer")
  if (!is.numeric(migration_rate) || migration_rate < 0) stop("migration_rate must be non-negative")
  if (!is.numeric(Chr.lengths) || any(Chr.lengths < 1)) stop("Chr.lengths must be positive integers")
  
  # Create initial parents and population
  parents <- makeparents(Chr.lengths)
  PopulationA <- vector("list", PopSize)
  names(PopulationA) <- paste0("Ind", 1:PopSize)
  PopulationA[] <- list(parents$pA)
  
  # Initialize frequency tracking
  frequency_history <- matrix(0, nrow = Generations, ncol = sum(Chr.lengths))
  frequency_history[1,] <- get_allele_freq(PopulationA)
  
  # Run generations
  for(gen in 2:Generations) {
    # Add migrants
    Nmigrants <- rpois(1, migration_rate)
    total_pop <- PopSize + Nmigrants
    
    PopulationA_after_migration <- vector("list", total_pop)
    names(PopulationA_after_migration) <- paste0("Ind", 1:total_pop)
    PopulationA_after_migration[1:PopSize] <- PopulationA
    if(Nmigrants > 0) {
      PopulationA_after_migration[(PopSize+1):total_pop] <- list(parents$pB)
    }
    
    # Create next generation
    PopulationA_nextgen <- vector("list", PopSize)
    names(PopulationA_nextgen) <- paste0("Ind", 1:PopSize)
    
    for(i in 1:PopSize) {
      dam_sire <- sample(total_pop, 2)
      PopulationA_nextgen[[i]] <- mate(
        PopulationA_after_migration[[dam_sire[1]]],
        PopulationA_after_migration[[dam_sire[2]]]
      )
    }
    
    PopulationA <- PopulationA_nextgen
    frequency_history[gen,] <- get_allele_freq(PopulationA)
  }
  
  # Get final transcriptomes if parameters are provided
  transcriptomes <- NULL
  if (!is.null(params)) {
    transcriptomes <- getTranscriptomes(
      Population = PopulationA,
      N_genes = sum(Chr.lengths),
      intercepts = params$intercepts,
      betas = params$betas,
      effect_matrices = params$effect_matrices,
      neteQTLlocus = params$neteQTLlocus
    )
  }
  
  return(list(
    final_population = PopulationA,
    frequency_history = frequency_history,
    transcriptomes = transcriptomes
  ))
}

#' Run replicated introgression simulations
#' @param Generations Number of generations
#' @param PopSize Population size
#' @param Chr.lengths Vector of chromosome lengths
#' @param migration_rate Rate of migration
#' @param n_replicates Number of replicates
#' @param parallel Whether to use parallel processing
#' @param params List of simulation parameters
#' @param n_cores Number of cores for parallel processing
#' @return List containing simulation results and summary statistics
run_replicated_introgression <- function(Generations, PopSize, Chr.lengths, migration_rate,
                                         n_replicates = 10, parallel = TRUE, params = NULL,
                                         n_cores = parallel::detectCores() - 1) {
  
  if (parallel) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    # Export all required core functions to the cluster
    core_functions <- c(
      # Core simulation functions
      "makeparents", "get.gamete", "mate", "get_allele_freq",
      # Network and expression functions
      "generate_transcriptome_parameters", "get_effect_matrix",
      "genotype_to_netexpr", "getTranscriptomes",
      # Main simulation function
      "introgressing_with_tracking"
    )
    
    # Export functions and required objects
    parallel::clusterExport(cl, core_functions, envir = environment())
    if (!is.null(params)) {
      parallel::clusterExport(cl, "params", envir = environment())
    }
    
    # Load required packages on each cluster node
    parallel::clusterEvalQ(cl, {
      library(stats)  # For rnorm, rpois, etc.
    })
    
    results <- pbapply::pblapply(1:n_replicates, function(rep) {
      introgressing_with_tracking(
        Generations = Generations,
        PopSize = PopSize,
        Chr.lengths = Chr.lengths,
        migration_rate = migration_rate,
        params = params
      )
    }, cl = cl)
    
  } else {
    results <- pbapply::pblapply(1:n_replicates, function(rep) {
      introgressing_with_tracking(
        Generations = Generations,
        PopSize = PopSize,
        Chr.lengths = Chr.lengths,
        migration_rate = migration_rate,
        params = params
      )
    })
  }
  
  # Calculate summary statistics
  n_loci <- ncol(results[[1]]$frequency_history)
  n_generations <- nrow(results[[1]]$frequency_history)
  
  mean_frequencies <- sd_frequencies <- ci_lower <- ci_upper <- 
    array(0, dim = c(n_generations, n_loci))
  
  for(gen in 1:n_generations) {
    for(locus in 1:n_loci) {
      freq_values <- sapply(results, function(x) x$frequency_history[gen, locus])
      mean_frequencies[gen, locus] <- mean(freq_values)
      sd_frequencies[gen, locus] <- sd(freq_values)
      ci <- quantile(freq_values, c(0.025, 0.975))
      ci_lower[gen, locus] <- ci[1]
      ci_upper[gen, locus] <- ci[2]
    }
  }
  
  return(list(
    replicates = results,
    summary = list(
      mean_frequencies = mean_frequencies,
      sd_frequencies = sd_frequencies,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    ),
    parameters = list(
      Generations = Generations,
      PopSize = PopSize,
      Chr.lengths = Chr.lengths,
      migration_rate = migration_rate,
      n_replicates = n_replicates
    )
  ))
}

#' Run network simulations
#' @param migration_rates Vector of migration rates to test
#' @param params List of simulation parameters
#' @param n_replicates Number of replicates per migration rate
#' @param output_dir Output directory
#' @return List of results for each migration rate
run_network_simulations <- function(migration_rates, params, n_replicates = 10,
                                    output_dir = "results/network_analysis") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_list <- list()
  
  for(rate in migration_rates) {
    cat("\nRunning simulations for migration rate:", rate, "\n")
    
    replicate_results <- run_replicated_introgression(
      Generations = params$Generations,
      PopSize = params$PopSize,
      Chr.lengths = params$Chr.lengths,
      migration_rate = rate,
      n_replicates = n_replicates,
      parallel = TRUE,
      params = params
    )
    
    results_list[[as.character(rate)]] <- replicate_results
  }
  
  saveRDS(results_list, file.path(output_dir, "network_simulation_results.rds"))
  
  return(results_list)
}

#' Run tracked simulations
#' @param migration_rates Vector of migration rates to test
#' @param params List of simulation parameters
#' @param n_replicates Number of replicates per migration rate
#' @param output_dir Output directory
#' @return List of results for each migration rate
run_tracked_simulations <- function(migration_rates, params, n_replicates = 10,
                                    output_dir = "results/introgression_tracking") {
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
  
  results_list <- list()
  
  for(rate in migration_rates) {
    cat("\nRunning simulations for migration rate:", rate, "\n")
    
    replicate_results <- run_replicated_introgression(
      Generations = params$Generations,
      PopSize = params$PopSize,
      Chr.lengths = params$Chr.lengths,
      migration_rate = rate,
      n_replicates = n_replicates,
      parallel = TRUE,
      params = params
    )
    
    results_list[[as.character(rate)]] <- replicate_results
    
    # Generate and save individual plots
    p_individual <- plot_introgression_summary(replicate_results, params$Chr.lengths)
    ggsave(
      filename = file.path(output_dir, "plots", sprintf("migration_rate_%s_summary.pdf", rate)),
      plot = p_individual,
      width = 12,
      height = 10
    )
  }
  
  # Create combined visualizations
  plots <- create_combined_visualizations(results_list, output_dir)
  
  # Save complete results
  saveRDS(results_list, file.path(output_dir, "data", "all_results.rds"))
  
  return(list(
    results = results_list,
    plots = plots
  ))
}

# Generate a summary report
generate_summary_report <- function(results, config, output_dir) {
  # Create report file
  report_file <- file.path(output_dir, "simulation_summary.txt")
  
  sink(report_file)
  
  cat("Introgression Simulation Summary\n")
  cat("==============================\n\n")
  
  # Parameters section
  cat("Simulation Parameters\n")
  cat("--------------------\n")
  cat("Basic Parameters:\n")
  cat(sprintf("- Number of genes: %d\n", sum(results$params$Chr.lengths)))
  cat(sprintf("- Population size: %d\n", results$params$PopSize))
  cat(sprintf("- Generations: %d\n", results$params$Generations))
  cat(sprintf("- Number of replicates: %d\n", config$n_replicates))
  cat(sprintf("- Chromosome lengths: %s\n", paste(results$params$Chr.lengths, collapse = ", ")))
  
  cat("\nNetwork Parameters:\n")
  cat(sprintf("- Density (00): %.3f\n", config$network$density00))
  cat(sprintf("- Effect (00): %.3f\n", config$network$effect00))
  cat(sprintf("- Density (11): %.3f\n", config$network$density11))
  cat(sprintf("- Effect (11): %.3f\n", config$network$effect11))
  
  cat("\nExpression Parameters:\n")
  cat(sprintf("- Mean expression: %.2f\n", config$expression$mean))
  cat(sprintf("- Expression SD: %.2f\n", config$expression$sd))
  cat(sprintf("- Beta SD: %.2f\n", config$expression$beta_sd))
  
  # Migration rates
  cat("\nMigration Rates\n")
  cat("--------------\n")
  cat("Network Analysis Rates:\n")
  cat(sprintf("- Range: %.6f to %.6f\n", 
              min(config$migration_rates$network),
              max(config$migration_rates$network)))
  cat("Tracking Analysis Rates:\n")
  cat(sprintf("- Values: %s\n", 
              paste(config$migration_rates$tracking, collapse = ", ")))
  
  # Results summary
  cat("\nSimulation Results\n")
  cat("-----------------\n")
  
  # Network results
  if(!is.null(results$network_metrics)) {
    cat("\nNetwork Analysis Summary:\n")
    metrics_summary <- aggregate(
      . ~ migration_rate,
      data = results$network_metrics[c("migration_rate", "mean_degree", "transitivity", "modularity")],
      FUN = mean
    )
    cat(sprintf("- Number of migration rates analyzed: %d\n", nrow(metrics_summary)))
    cat(sprintf("- Average mean degree range: %.2f to %.2f\n",
                min(metrics_summary$mean_degree),
                max(metrics_summary$mean_degree)))
    cat(sprintf("- Average transitivity range: %.2f to %.2f\n",
                min(metrics_summary$transitivity),
                max(metrics_summary$transitivity)))
    cat(sprintf("- Average modularity range: %.2f to %.2f\n",
                min(metrics_summary$modularity),
                max(metrics_summary$modularity)))
  }
  
  # Tracking results
  if(!is.null(results$tracking_results)) {
    cat("\nTracking Analysis Summary:\n")
    cat(sprintf("- Number of migration rates analyzed: %d\n",
                length(results$tracking_results$results)))
    cat(sprintf("- Number of generations tracked: %d\n",
                nrow(results$tracking_results$results[[1]]$summary$mean_frequencies)))
  }
  
  # Output files
  cat("\nOutput Files\n")
  cat("------------\n")
  cat("1. Data Files:\n")
  cat("   - simulation_parameters.rds\n")
  cat("   - network_metrics.csv\n")
  cat("   - simulation_config.rds\n")
  
  cat("\n2. Visualization Files:\n")
  cat("   - migration_rates_grid.pdf\n")
  cat("   - network_metrics_summary.pdf\n")
  cat("   - averaged_heatmaps/*.png\n")
  cat("   - averaged_networks/*.png\n")
  
  # Close the report
  sink()
  
  # Also print a brief summary to console
  cat("\nSummary report generated at:", report_file, "\n")
  cat("Number of migration rates analyzed:\n")
  cat("- Network analysis:", length(config$migration_rates$network), "\n")
  cat("- Tracking analysis:", length(config$migration_rates$tracking), "\n")
}

