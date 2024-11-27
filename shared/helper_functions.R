# Core genetics and population operations shared between
# both the hybridization and introgression simulations

# 1. Function to create parents ####################
#' Create parent populations with specified chromosome lengths
#' @param Chr.lengths Vector specifying length of each chromosome
#' @return List containing parent populations A and B
makeparents <- function(Chr.lengths) {
  Chr.names <- paste(rep(paste("Chr", c(1:length(Chr.lengths)), sep = ""), each = 2), c(".1", ".2"), sep = "")
  
  ParentA <- sapply(Chr.names, function(x) NULL)
  ParentB <- sapply(Chr.names, function(x) NULL)
  
  for(i in seq(2, length(Chr.names), by = 2)) {
    ParentA[[i-1]] <- ParentA[[i]] <- rep(0, Chr.lengths[i/2])
    ParentB[[i-1]] <- ParentB[[i]] <- rep(1, Chr.lengths[i/2])
  }
  
  return(list(pA = ParentA, pB = ParentB))
}

# 2. Function for recombination and mating ##########
#' Generate a gamete through recombination
#' @param parent Parent genotype
#' @return Recombined gamete
get.gamete <- function(parent) {
  Chr.names.p <- names(parent)
  Chr.names <- Chr.names.p[seq(1, length(Chr.names.p)-1, by = 2)]
  gamete <- sapply(Chr.names, function(x) NULL)
  
  for(i in seq(2, length(Chr.names.p), by = 2)) {
    chromatid1 <- parent[[i-1]]
    chromatid2 <- parent[[i]]
    
    if(runif(1) < 0.5) {
      chromatid1.1 <- chromatid2
      chromatid2.1 <- chromatid1
      chromatid1 <- chromatid1.1
      chromatid2 <- chromatid2.1
    }
    
    breakpoint <- sample(c(2:(length(chromatid1)-1)), 1)
    after_crossingover <- c(chromatid1[1:breakpoint], chromatid2[(breakpoint+1):length(chromatid2)])
    gamete[[i/2]] <- after_crossingover
  }
  
  return(gamete)
}

# 3. Mating function #####################
#' Mate two parents to produce offspring
#' @param parent1,parent2 Parent genotypes
#' @return Offspring genotype
mate <- function(parent1, parent2) {
  gamete1 <- get.gamete(parent1)
  gamete2 <- get.gamete(parent2)
  Chr.names.p <- names(parent1)
  offspring <- sapply(Chr.names.p, function(x) NULL)
  
  for(i in seq(2, length(Chr.names.p), by = 2)) {
    offspring[[i-1]] <- gamete1[[i/2]]
    offspring[[i]] <- gamete2[[i/2]]
  }
  
  return(offspring)
}

# 4. Expression functions ###########################
#' Generates base expression without network effects
#' Calculates baseline gene expression levels for an individual based on their genotype,
#' before any network effects are applied. Expression follows a Poisson distribution
#' with mean determined by genotype effects.
#' @param genotype List of chromosomes, where each chromosome contains allele values (0/1)
#' @param intercepts Numeric vector of baseline expression levels for each gene
#' @param betas Numeric vector of genetic effect sizes
#' 
#' @return A list containing: 
#' expression: baseline expression values,  
#' genome_val: Numeric vector of diploid genotype values (0,1,2) for each gene
#' 

genotype_to_nullexpr <- function(genotype, intercepts, betas) {
  NChr <- length(genotype) 
  genome_val <- c()
  for(i in seq(1, NChr-1, by = 2)) {
    genome_val <- c(genome_val, genotype[[i]] + genotype[[i+1]])
  }
  N_genes <- length(genome_val)
  lambda <- exp(intercepts + genome_val*betas)
  expression <- rpois(n = N_genes, lambda)
  return(list(expression = expression, genome_val = genome_val))
}

#' Set up parallel processing cluster
#' @param n_replicates Number of replicates
#' @param required_objects Vector of object names to export
#' @param required_packages Vector of package names to load
setup_parallel_cluster <- function(n_replicates, required_objects, required_packages) {
  n_cores <- min(detectCores() - 1, n_replicates)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to cluster
  clusterExport(cl, required_objects)
  
  # Load required packages on each node
  if(!is.null(required_packages)) {
    clusterEvalQ(cl, {
      for(pkg in required_packages) {
        library(pkg, character.only = TRUE)
      }
    })
  }
  
  return(cl)
}

#' Create output directories
#' @param base_dir Base directory path
#' @param subdirs Vector of subdirectory names
create_output_dirs <- function(base_dir, subdirs) {
  for(dir in subdirs) {
    dir.create(
      file.path(base_dir, dir), 
      recursive = TRUE, 
      showWarnings = FALSE
    )
  }
}