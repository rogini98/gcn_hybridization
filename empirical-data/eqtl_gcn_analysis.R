# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("qvalue")

# Required packages
library(qvalue) # pkg to estimate q-vlaues and false discovery rates
library(HiClimR)

#' Filter and normalize gene expression data
#' @param tagseq Expression data matrix
#' @param metadata Sample metadata
#' @param min_reads Minimum mean read count (default 30)
#' @param min_presence Minimum proportion of samples with expression (default 0.8)
#' @return List containing filtered and normalized expression data
process_expression_data <- function(tagseq, metadata, min_reads = 30, min_presence = 0.8) {
  # Calculate mean reads per gene
  meanreads <- colSums(tagseq, na.rm = TRUE)/nrow(tagseq)
  
  # Calculate proportion of samples with expression
  presence <- colMeans(tagseq > 0, na.rm = TRUE)
  
  # Filter genes based on criteria
  keep_genes <- meanreads > min_reads & presence > min_presence
  tagseq_filtered <- tagseq[, keep_genes]
  
  # Calculate total reads per sample for normalization
  total_reads <- rowSums(tagseq_filtered, na.rm = TRUE)
  
  # Normalize to proportions
  tagseq_norm <- sweep(tagseq_filtered, 1, total_reads, '/')
  
  return(list(
    filtered = tagseq_filtered,
    normalized = tagseq_norm,
    total_reads = total_reads
  ))
}

#' Perform eQTL analysis for all marker-gene pairs
#' @param expr_data Processed expression data
#' @param marker_data Genetic marker data
#' @param metadata Sample metadata
#' @return Data frame of eQTL results
perform_eqtl_analysis <- function(expr_data, marker_data, metadata) {
  results <- data.frame()
  
  for(snp_i in 1:ncol(marker_data)) {
    for(gene_i in 1:ncol(expr_data$filtered)) {
      # Prepare response variable (read counts)
      Y <- cbind(
        expr_data$filtered[, gene_i],
        expr_data$total_reads - expr_data$filtered[, gene_i]
      )
      
      # Fit quasibinomial GLM
      tryCatch({
        model <- glm(
          Y ~ as.factor(marker_data[, snp_i]),
          family = "quasibinomial"
        )
        
        summ <- summary(model)
        
        # Store results
        results <- rbind(results, data.frame(
          marker = colnames(marker_data)[snp_i],
          gene = colnames(expr_data$filtered)[gene_i],
          coef = coef(model)[2],
          se = summ$coefficients[2,2],
          p_value = summ$coefficients[2,4]
        ))
      }, error = function(e) {
        message(sprintf("Error in marker %d, gene %d: %s", snp_i, gene_i, e$message))
      })
    }
  }
  
  return(results)
}

#' Analyze gene co-expression patterns
#' @param expr_data Normalized expression data
#' @param metadata Sample metadata
#' @param gene_map Gene position information
#' @return List containing correlation results and statistics
analyze_coexpression <- function(expr_data, metadata, gene_map) {
  # Get residuals controlling for covariates
  expr_resid <- matrix(NA, nrow = nrow(expr_data$normalized), 
                       ncol = ncol(expr_data$normalized))
  
  for(i in 1:ncol(expr_data$normalized)) {
    model <- lm(expr_data$normalized[,i] ~ 
                  metadata$sex + 
                  metadata$mass + 
                  metadata$parasite_infection)
    expr_resid[,i] <- residid(model)
  }
  
  # Calculate pairwise correlations
  cor_matrix <- fastCor(t(expr_resid), upperTri = TRUE)
  
  # Create chromosome comparison matrix
  same_chr <- matrix(FALSE, nrow = ncol(expr_resid), ncol = ncol(expr_resid))
  for(i in 1:nrow(same_chr)) {
    for(j in 1:ncol(same_chr)) {
      same_chr[i,j] <- gene_map$chromosome[i] == gene_map$chromosome[j]
    }
  }
  
  # Test for correlation differences between same/different chromosomes
  same_chr_cor <- abs(cor_matrix[same_chr])
  diff_chr_cor <- abs(cor_matrix[!same_chr])
  chr_test <- t.test(same_chr_cor, diff_chr_cor)
  
  return(list(
    correlations = cor_matrix,
    chromosome_test = chr_test,
    same_chr_mean = mean(same_chr_cor, na.rm = TRUE),
    diff_chr_mean = mean(diff_chr_cor, na.rm = TRUE)
  ))
}

# Main analysis pipeline
run_analysis <- function(data_file, gene_map_file) {
  # Load data
  data <- read.csv(data_file)
  gene_map <- read.csv(gene_map_file)
  
  # Split data
  metadata <- data[,1:112]
  markers <- data[,113:346]
  expression <- data[,441:ncol(data)]
  
  # Process expression data
  expr_data <- process_expression_data(expression, metadata)
  
  # Perform eQTL analysis
  eqtl_results <- perform_eqtl_analysis(expr_data, markers, metadata)
  
  # Analyze co-expression
  coexp_results <- analyze_coexpression(expr_data, metadata, gene_map)
  
  return(list(
    eqtl = eqtl_results,
    coexpression = coexp_results
  ))
}

# Run analysis on data
results <- run_analysis("./data/MergedTagseqQTL_GENOCORRECTED.csv", "./data/MapENSGACT_ENSGACG.csv")
print(paste("Number of genes analyzed:", ncol(results$eqtl)))
print(results$coexpression$chromosome_test)