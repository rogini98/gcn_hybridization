# Core functions for network simulation and analysis
# scenario 1: simulate hybridization
library(igraph)

#############################################
### 1. Function to create parents
makeparents <- function(Chr.lengths) {
  Chr.names <- paste(rep(paste("Chr", c(1:length(Chr.lengths)), sep = ""), each = 2), c(".1", ".2"), sep = "")
  
  ParentA <- sapply(Chr.names, function(x) NULL)
  ParentB <- sapply(Chr.names, function(x) NULL)
  for(i in seq(2, length(Chr.names), by = 2)) {
    ParentA[[i-1]] <- rep(0, Chr.lengths[i/2])
    ParentA[[i]] <- rep(0, Chr.lengths[i/2])
    ParentB[[i-1]] <- rep(1, Chr.lengths[i/2])
    ParentB[[i]] <- rep(1, Chr.lengths[i/2])
  }
  return(list(pA = ParentA, pB = ParentB))
}

#############################################
### 2. Functions for recombination and mating
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

#############################################
### 3. Mating function
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

#############################################
### 4. Population breeding
breeding <- function(Fgen, PopSize, Parent1, Parent2) {
  Population_parents <- sapply(paste("Ind", c(1:PopSize), sep = ""), function(x) NULL)
  
  if(Fgen == 0) {
    for(i in 1:PopSize) {
      Population_parents[[i]] <- Parent1
    }
  } else {
    Population_parents <- list(Parent1, Parent2)
    Population_nextgen <- sapply(paste("Ind", c(1:PopSize), sep = ""), function(x) NULL)
    
    for(gen in 1:Fgen) {
      for(i in 1:PopSize) {
        Nparents <- length(Population_parents)
        dam_sire <- sample(c(1:Nparents), 2)
        par1 <- Population_parents[[dam_sire[1]]]
        par2 <- Population_parents[[dam_sire[2]]]
        Population_nextgen[[i]] <- mate(par1, par2)
      }
      Population_parents <- Population_nextgen
    }
  }
  return(Population_parents)
}

#############################################
### 5. Expression functions
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

#############################################
### 6. Network architecture
get_effect.matrix <- function(N_genes, network_density, network_effect, type = "BarabasiAlbert", powerval = 1) {
  if(type == "binom") {
    incidences <- rbinom(n = N_genes^2, size = 1, prob = network_density)
    effects <- rnorm(n = N_genes^2, mean = 0, sd = network_effect)
    effect.matrix <- matrix(incidences*effects, nrow = N_genes)
  }
  if(type == "BarabasiAlbert") {
    netw <- sample_pa(n=N_genes, powerval, m=1, directed=F)
    netw.matrix <- as.matrix(as_adj(netw))
    incidences <- rbinom(n = N_genes^2, size = 1, prob = network_density)
    effects <- rnorm(n = N_genes^2, mean = 0, sd = network_effect)
    subsample.matrix <- matrix(incidences*effects, nrow = N_genes)
    effect.matrix <- incidences*subsample.matrix
  }
  return(effect.matrix)
}

#############################################
### 7. Network expression
genotype_to_netexpr <- function(genotype, intercepts, betas, effect.matrices, neteQTLlocus) {
  placeholder <- genotype_to_nullexpr(genotype, intercepts, betas)
  genome_val <- placeholder$genome_val
  null_expression <- placeholder$expression
  lambda <- exp(intercepts + betas)
  N_genes <- length(null_expression)
  effect.matrix <- effect.matrices[[genome_val[neteQTLlocus] + 1]]
  
  for(i in 1:(N_genes-1)) {
    for(j in i:N_genes) {
      null_expression[j] <- null_expression[j] + effect.matrix[i,j]*(null_expression[i] - lambda[i])
    }
  }
  return(null_expression)
}

#############################################
### 8. Population transcriptomes
getTranscriptomes <- function(Population, N_genes, intercepts, betas, effect.matrices, neteQTLlocus) {
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
    readcounts[i,] <- genotype_to_netexpr(ind1, intercepts, betas, effect.matrices, neteQTLlocus)
  }
  return(list(g = genotypes, r = readcounts))
}

#############################################
### 9. Network analysis functions
coexpr <- function(genoreads, threshold) {
  cormat <- cor(genoreads$r)
  if(!is.na(threshold)) {
    cormat[abs(cormat) < threshold] <- 0
    cormat[abs(cormat) > threshold] <- 1
  }
  return(cormat)
}

#############################################
### 10. Main simulation function
netQTLsim <- function(Fgen, Chr.lengths, PopSize, intercepts, betas, effect.matrices, neteQTLlocus) {
  N_genes <- sum(Chr.lengths)
  parents <- makeparents(Chr.lengths)
  F_i_population <- breeding(Fgen, PopSize, parents$pA, parents$pB)
  genoreads <- getTranscriptomes(F_i_population, N_genes, intercepts, betas, effect.matrices, neteQTLlocus)
  co_expr <- coexpr(genoreads, threshold = 0.2)
  
  return(list(
    genoreads = genoreads,
    F_i_population = F_i_population,
    co_expr = co_expr
  ))
}