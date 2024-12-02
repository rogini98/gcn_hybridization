# Required packages
library(qvalue)
library(HiClimR)
library(data.table) # For more efficient data handling

#' Process expression data with memory-efficient operations
#' @param tagseq Expression data matrix
#' @param metadata Sample metadata
#' @param min_reads Minimum mean read count (default 30)
#' @param min_presence Minimum proportion of samples with expression (default 0.8)
process_expression_data <- function(tagseq, metadata, min_reads = 30, min_presence = 0.8) {
  # Convert to data.table for more efficient operations
  tagseq_dt <- as.data.table(tagseq)
  
  # Calculate metrics using data.table operations
  meanreads <- colSums(tagseq_dt, na.rm = TRUE)/nrow(tagseq_dt)
  presence <- colMeans(tagseq_dt > 0, na.rm = TRUE)
  
  # Filter genes
  keep_genes <- which(meanreads > min_reads & presence > min_presence)
  tagseq_filtered <- tagseq_dt[, ..keep_genes]
  
  # Normalize (using matrix operations for efficiency)
  total_reads <- rowSums(tagseq_filtered, na.rm = TRUE)
  tagseq_norm <- sweep(as.matrix(tagseq_filtered), 1, total_reads, '/')
  
  return(list(
    filtered = tagseq_filtered,
    normalized = tagseq_norm,
    total_reads = total_reads
  ))
}

#' Perform eQTL analysis in batches
#' @param expr_data Processed expression data
#' @param marker_data Genetic marker data
#' @param batch_size Number of genes to process at once
perform_eqtl_analysis <- function(expr_data, marker_data, metadata, batch_size = 100) {
  # Pre-allocate results list
  total_genes <- ncol(expr_data$filtered)
  total_markers <- ncol(marker_data)
  n_batches <- ceiling(total_genes / batch_size)
  
  all_results <- vector("list", n_batches)
  
  # Process in batches
  for(batch in 1:n_batches) {
    start_idx <- (batch-1) * batch_size + 1
    end_idx <- min(batch * batch_size, total_genes)
    
    # Create batch results
    batch_results <- data.frame(
      marker = character(),
      gene = character(),
      coef = numeric(),
      se = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Process each gene in batch
    for(gene_i in start_idx:end_idx) {
      Y <- cbind(
        expr_data$filtered[[gene_i]],
        expr_data$total_reads - expr_data$filtered[[gene_i]]
      )
      
      # Process each marker for this gene
      for(snp_i in 1:total_markers) {
        tryCatch({
          model <- glm(
            Y ~ as.factor(marker_data[[snp_i]]),
            family = "quasibinomial"
          )
          
          summ <- summary(model)
          
          batch_results <- rbind(batch_results, data.frame(
            marker = names(marker_data)[snp_i],
            gene = names(expr_data$filtered)[gene_i],
            coef = coef(model)[2],
            se = summ$coefficients[2,2],
            p_value = summ$coefficients[2,4]
          ))
        }, error = function(e) {
          message(sprintf("Error in batch %d, gene %d, marker %d: %s", 
                          batch, gene_i, snp_i, e$message))
        })
      }
    }
    
    all_results[[batch]] <- batch_results
    
    # Clean up memory
    gc()
    
    # Progress update
    cat(sprintf("Completed batch %d of %d\n", batch, n_batches))
  }
  
  # Combine results
  final_results <- do.call(rbind, all_results)
  return(final_results)
}

#' Memory-efficient co-expression analysis
#' @param expr_data Normalized expression data
#' @param metadata Sample metadata
#' @param gene_map Gene position information
#' @param batch_size Size of gene batches for correlation
analyze_coexpression <- function(expr_data, metadata, gene_map, batch_size = 500) {
  # Calculate residuals in batches
  n_genes <- ncol(expr_data$normalized)
  n_batches <- ceiling(n_genes / batch_size)
  expr_resid <- matrix(NA, nrow = nrow(expr_data$normalized), 
                       ncol = ncol(expr_data$normalized))
  
  for(batch in 1:n_batches) {
    start_idx <- (batch-1) * batch_size + 1
    end_idx <- min(batch * batch_size, n_genes)
    
    for(i in start_idx:end_idx) {
      model <- lm(expr_data$normalized[,i] ~ 
                    metadata$sex + 
                    metadata$mass + 
                    metadata$parasite_infection)
      expr_resid[,i] <- residuals(model)
    }
    
    cat(sprintf("Processed residuals batch %d of %d\n", batch, n_batches))
  }
  
  # Calculate correlations using fastCor (already optimized)
  cor_matrix <- fastCor(t(expr_resid), upperTri = TRUE)
  
  # Efficient chromosome comparison
  same_chr <- outer(gene_map$chromosome, gene_map$chromosome, "==")
  
  # Calculate statistics
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

# Main analysis function with progress reporting
run_analysis <- function(data_file, gene_map_file, batch_size = 100) {
  cat("Loading data...\n")
  data <- fread(data_file)  # Using data.table's faster reader
  gene_map <- fread(gene_map_file)
  
  cat("Splitting data...\n")
  metadata <- data[, 1:112]
  markers <- data[, 113:346]
  expression <- data[, 441:ncol(data)]
  
  cat("Processing expression data...\n")
  expr_data <- process_expression_data(expression, metadata)
  
  cat("Performing eQTL analysis...\n")
  eqtl_results <- perform_eqtl_analysis(expr_data, markers, metadata, batch_size)
  
  cat("Analyzing co-expression...\n")
  coexp_results <- analyze_coexpression(expr_data, metadata, gene_map, batch_size)
  
  return(list(
    eqtl = eqtl_results,
    coexpression = coexp_results
  ))
}

# Example usage with smaller batch size for laptop
results <- run_analysis(
  "./data/MergedTagseqQTL_GENOCORRECTED.csv",
  "./data/MapENSGACT_ENSGACG.csv",
  batch_size = 50  # Adjust based on available memory
)