# Core introgression simulation functions
# Authors: Original by Dan Bolnick, reorganized by Rogini Runghen
# Last updated: November 25, 2024

#############################################
### Core Functions for introgression simulation
#############################################

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

#' Calculate allele frequencies
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

#' Create gene coexpression network architecture
#' @param N_genes Number of genes
#' @param network_density Network density
#' @param network_effect Network effect size
#' @return Adjacency matrix
get_effect_matrix <- function(N_genes, network_density, network_effect) {
  incidences <- rbinom(n = N_genes^2, size = 1, prob = network_density)
  effects <- rnorm(n = N_genes^2, mean = 0, sd = network_effect)
  effect.matrix <- matrix(incidences * effects, nrow = N_genes, ncol = N_genes)
  return(effect.matrix)
}

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