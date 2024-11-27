# Core functions for network simulation and analysis
# scenario 1: simulate hybridization
library(igraph)


# 1.Population breeding #####################
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

########################################
# Simulation runner functions
########################################
# 2. Main simulation function ###################
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