# simulation_runners.R
# Functions for running introgression simulations with network analysis and tracking
# Authors: Original by Dan Bolnick, reorganized by Rogini Runghen
# Last updated: November 25, 2024
library(parallel)
library(pbapply)
library(ggplot2)
library(patchwork)

# Source core functions
source("./introgression/01-network_functions_introgressionSim.R")
source("./introgression/02-visualization_functions_introgression.R")

#############################################
### Parameter Initialization
#############################################

#' Initialize all simulation parameters
#' @return List of all parameters needed for simulation
initialize_parameters <- function() {
  # Basic parameters
  Chr.lengths <- c(10, 10, 5, 5, 10)
  N_genes <- sum(Chr.lengths)
  
  # Define transcriptome parameters
  meanexpr <- 5
  sd.expr <- 2
  sd.beta <- 3
  
  # Generate transcriptome parameters
  trans_params <- generate_transcriptome_parameters(
    N_genes = N_genes,
    meanexpr = meanexpr,
    sd.expr = sd.expr,
    sd.beta = sd.beta
  )
  
  # Define network parameters
  network_density00 <- 0.02
  network_effect00 <- 2
  network_density11 <- 0.03
  network_effect11 <- 1.5
  
  # Generate network architectures
  effect.matrix00 <- get_effect_matrix(N_genes, network_density00, network_effect00)
  effect.matrix11 <- get_effect_matrix(N_genes, network_density11, network_effect11)
  effect.matrix01.add <- (effect.matrix00 + effect.matrix11) / 2
  
  return(list(
    Generations = 100,
    PopSize = 500,
    Chr.lengths = Chr.lengths,
    intercepts = trans_params$intercepts,
    betas = trans_params$betas,
    neteQTLlocus = trans_params$neteQTLlocus,
    effect_matrices = list(effect.matrix00, effect.matrix01.add, effect.matrix11),
    network_params = list(
      density00 = network_density00,
      effect00 = network_effect00,
      density11 = network_density11,
      effect11 = network_effect11
    )
  ))
}

#############################################
### Core Simulation Functions
#############################################

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

#' #############################################
#' ### Replicated Simulation Functions
#' #############################################
#' 
#' #' Run replicated introgression simulations
#' #' @param Generations Number of generations
#' #' @param PopSize Population size
#' #' @param Chr.lengths Vector of chromosome lengths
#' #' @param migration_rate Rate of migration
#' #' @param n_replicates Number of replicates
#' #' @param parallel Whether to use parallel processing
#' #' @param params List of simulation parameters
#' #' @param n_cores Number of cores for parallel processing
#' #' @return List containing simulation results and summary statistics
#' run_replicated_introgression <- function(Generations, PopSize, Chr.lengths, migration_rate,
#'                                          n_replicates = 10, parallel = TRUE, params = NULL,
#'                                          n_cores = parallel::detectCores() - 1) {
#'   
#'   if (parallel) {
#'     cl <- parallel::makeCluster(n_cores)
#'     on.exit(parallel::stopCluster(cl))
#'     
#'     parallel::clusterExport(cl, 
#'                             c("introgressing_with_tracking", "makeparents", "mate", 
#'                               "get.gamete", "get_allele_freq", "getTranscriptomes"),
#'                             envir = environment())
#'     
#'     if (!is.null(params)) {
#'       parallel::clusterExport(cl, "params", envir = environment())
#'     }
#'     
#'     results <- pbapply::pblapply(1:n_replicates, function(rep) {
#'       introgressing_with_tracking(
#'         Generations = Generations,
#'         PopSize = PopSize,
#'         Chr.lengths = Chr.lengths,
#'         migration_rate = migration_rate,
#'         params = params
#'       )
#'     }, cl = cl)
#'     
#'   } else {
#'     results <- pbapply::pblapply(1:n_replicates, function(rep) {
#'       introgressing_with_tracking(
#'         Generations = Generations,
#'         PopSize = PopSize,
#'         Chr.lengths = Chr.lengths,
#'         migration_rate = migration_rate,
#'         params = params
#'       )
#'     })
#'   }
#'   
#'   # Calculate summary statistics
#'   n_loci <- ncol(results[[1]]$frequency_history)
#'   n_generations <- nrow(results[[1]]$frequency_history)
#'   
#'   mean_frequencies <- sd_frequencies <- ci_lower <- ci_upper <- 
#'     array(0, dim = c(n_generations, n_loci))
#'   
#'   for(gen in 1:n_generations) {
#'     for(locus in 1:n_loci) {
#'       freq_values <- sapply(results, function(x) x$frequency_history[gen, locus])
#'       mean_frequencies[gen, locus] <- mean(freq_values)
#'       sd_frequencies[gen, locus] <- sd(freq_values)
#'       ci <- quantile(freq_values, c(0.025, 0.975))
#'       ci_lower[gen, locus] <- ci[1]
#'       ci_upper[gen, locus] <- ci[2]
#'     }
#'   }
#'   
#'   return(list(
#'     replicates = results,
#'     summary = list(
#'       mean_frequencies = mean_frequencies,
#'       sd_frequencies = sd_frequencies,
#'       ci_lower = ci_lower,
#'       ci_upper = ci_upper
#'     ),
#'     parameters = list(
#'       Generations = Generations,
#'       PopSize = PopSize,
#'       Chr.lengths = Chr.lengths,
#'       migration_rate = migration_rate,
#'       n_replicates = n_replicates
#'     )
#'   ))
#' }
#' 
#' #' Run network simulations
#' #' @param migration_rates Vector of migration rates to test
#' #' @param params List of simulation parameters
#' #' @param n_replicates Number of replicates per migration rate
#' #' @param output_dir Output directory
#' #' @return List of results for each migration rate
#' run_network_simulations <- function(migration_rates, params, n_replicates = 10,
#'                                     output_dir = "results/network_analysis") {
#'   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
#'   
#'   results_list <- list()
#'   
#'   for(rate in migration_rates) {
#'     cat("\nRunning simulations for migration rate:", rate, "\n")
#'     
#'     replicate_results <- run_replicated_introgression(
#'       Generations = params$Generations,
#'       PopSize = params$PopSize,
#'       Chr.lengths = params$Chr.lengths,
#'       migration_rate = rate,
#'       n_replicates = n_replicates,
#'       parallel = TRUE,
#'       params = params
#'     )
#'     
#'     results_list[[as.character(rate)]] <- replicate_results
#'   }
#'   
#'   saveRDS(results_list, file.path(output_dir, "network_simulation_results.rds"))
#'   
#'   return(results_list)
#' }
#' 
#' #' Run tracked simulations
#' #' @param migration_rates Vector of migration rates to test
#' #' @param params List of simulation parameters
#' #' @param n_replicates Number of replicates per migration rate
#' #' @param output_dir Output directory
#' #' @return List of results for each migration rate
#' run_tracked_simulations <- function(migration_rates, params, n_replicates = 10,
#'                                     output_dir = "results/introgression_tracking") {
#'   dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
#'   dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
#'   
#'   results_list <- list()
#'   
#'   for(rate in migration_rates) {
#'     cat("\nRunning simulations for migration rate:", rate, "\n")
#'     
#'     replicate_results <- run_replicated_introgression(
#'       Generations = params$Generations,
#'       PopSize = params$PopSize,
#'       Chr.lengths = params$Chr.lengths,
#'       migration_rate = rate,
#'       n_replicates = n_replicates,
#'       parallel = TRUE,
#'       params = params
#'     )
#'     
#'     results_list[[as.character(rate)]] <- replicate_results
#'     
#'     # Generate and save individual plots
#'     p_individual <- plot_introgression_summary(replicate_results, params$Chr.lengths)
#'     ggsave(
#'       filename = file.path(output_dir, "plots", sprintf("migration_rate_%s_summary.pdf", rate)),
#'       plot = p_individual,
#'       width = 12,
#'       height = 10
#'     )
#'   }
#'   
#'   # Create combined visualizations
#'   plots <- create_combined_visualizations(results_list, output_dir)
#'   
#'   # Save complete results
#'   saveRDS(results_list, file.path(output_dir, "data", "all_results.rds"))
#'   
#'   return(list(
#'     results = results_list,
#'     plots = plots
#'   ))
#' }
#' 
#' # # Example usage
#' # if (FALSE) {
#' #   # Initialize parameters
#' #   params <- initialize_parameters()
#' #   
#' #   # Run network simulations
#' #   network_results <- run_network_simulations(
#' #     migration_rates = exp(seq(-10, 0, by = 2)),
#' #     params = params,
#' #     n_replicates = 10
#' #   )
#' #   
#' #   # Run tracking simulations
#' #   tracking_results <- run_tracked_simulations(
#' #     migration_rates = c(0, 0.05, 0.5, 1),
#' #     params = params,
#' #     n_replicates = 10
#' #   )
#' # }