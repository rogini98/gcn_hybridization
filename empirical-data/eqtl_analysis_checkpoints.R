# Script containing all of the core functions to run the eQTL analysis 
# and correlation test

# required libraries
library(qvalue)
library(HiClimR)
library(data.table)

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
  total_genes <- ncol(expr_data$filtered)
  total_markers <- ncol(marker_data)
  n_batches <- ceiling(total_genes / batch_size)
  
  all_results <- vector("list", n_batches)
  
  for(batch in 1:n_batches) {
    start_idx <- (batch-1) * batch_size + 1
    end_idx <- min(batch * batch_size, total_genes)
    
    batch_results <- data.frame(
      marker = character(),
      gene = character(),
      coef = numeric(),
      se = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for(gene_i in start_idx:end_idx) {
      Y <- cbind(
        expr_data$filtered[[gene_i]],
        expr_data$total_reads - expr_data$filtered[[gene_i]]
      )
      
      for(snp_i in 1:total_markers) {
        tryCatch({
          model <- glm(Y ~ as.factor(marker_data[[snp_i]]), family = "quasibinomial")
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
    gc()
    cat(sprintf("Completed batch %d of %d\n", batch, n_batches))
  }
  
  final_results <- do.call(rbind, all_results)
  return(final_results)
}

#' Analyze co-expression patterns
#' @param expr_data Normalized expression data
#' @param gene_map Gene mapping information
analyze_coexpression <- function(expr_data, gene_map) {
  # Clean gene names
  expr_genes <- colnames(expr_data$normalized)
  if(any(grepl("\\.\\d+$", expr_genes))) {
    cat("Removing version numbers from gene names...\n")
    expr_genes_clean <- sub("\\.\\d+$", "", expr_genes)
    colnames(expr_data$normalized) <- expr_genes_clean
  }
  
  # Find matching genes
  common_genes <- intersect(colnames(expr_data$normalized), gene_map$ENSGACT)
  cat("Number of matching genes:", length(common_genes), "\n")
  
  if(length(common_genes) == 0) {
    stop("No matching genes found between expression data and gene map")
  }
  
  # Subset data to matching genes
  expr_matrix <- expr_data$normalized[, common_genes, drop = FALSE]
  matched_gene_map <- gene_map[match(common_genes, gene_map$ENSGACT), ]
  
  # Calculate correlations
  cat("Calculating correlations...\n")
  cor_matrix <- fastCor(t(expr_matrix), upperTri = TRUE)
  
  # Process chromosome data
  cat("Processing chromosome data...\n")
  chr_info <- matched_gene_map$LG
  
  # Print chromosome information
  cat("Unique chromosomes found:", paste(sort(unique(chr_info)), collapse=", "), "\n")
  cat("Number of genes per chromosome:\n")
  print(table(chr_info))
  
  # Remove NA chromosomes if any
  valid_idx <- !is.na(chr_info)
  if(!all(valid_idx)) {
    cat("Removing", sum(!valid_idx), "genes with missing chromosome information\n")
    cor_matrix <- cor_matrix[valid_idx, valid_idx]
    chr_info <- chr_info[valid_idx]
    common_genes <- common_genes[valid_idx]
  }
  
  # Convert correlation matrix to regular matrix if it isn't already
  cor_matrix <- as.matrix(cor_matrix)
  
  cat("Dimensions of correlation matrix after filtering:", dim(cor_matrix), "\n")
  
  # Process correlations
  unique_chrs <- sort(unique(chr_info))
  chr_stats <- list()
  same_chr_correlations <- numeric()
  diff_chr_correlations <- numeric()
  
  cat("Processing correlations by chromosome...\n")
  
  # Create matrix to store chromosome relationships
  n_genes <- length(chr_info)
  same_chr_matrix <- matrix(FALSE, nrow = n_genes, ncol = n_genes)
  
  # Fill the chromosome relationship matrix
  for(i in 1:n_genes) {
    for(j in i:n_genes) {
      same_chr_matrix[i,j] <- same_chr_matrix[j,i] <- (chr_info[i] == chr_info[j])
    }
  }
  
  # Get upper triangle indices
  upper_tri <- upper.tri(cor_matrix)
  
  # Extract correlations
  same_chr_cor <- cor_matrix[upper_tri & same_chr_matrix]
  diff_chr_cor <- cor_matrix[upper_tri & !same_chr_matrix]
  
  # Calculate chromosome-specific statistics
  for(chr in unique_chrs) {
    chr_idx <- chr_info == chr
    if(sum(chr_idx) > 1) {
      sub_matrix <- cor_matrix[chr_idx, chr_idx]
      chr_stats[[as.character(chr)]] <- list(
        chromosome = chr,
        mean_correlation = mean(sub_matrix[upper.tri(sub_matrix)], na.rm = TRUE),
        n_genes = sum(chr_idx)
      )
    }
  }
  
  # Perform t-test
  chr_test <- tryCatch({
    if(length(same_chr_cor) > 0 && length(diff_chr_cor) > 0) {
      test_result <- t.test(same_chr_cor, diff_chr_cor)
      list(
        statistic = test_result$statistic,
        p_value = test_result$p.value,
        mean_same_chr = mean(same_chr_cor, na.rm = TRUE),
        mean_diff_chr = mean(diff_chr_cor, na.rm = TRUE)
      )
    } else {
      NULL
    }
  }, error = function(e) {
    cat("Warning: Could not perform chromosome correlation t-test:", e$message, "\n")
    NULL
  })
  
  cat("Analysis complete.\n")
  
  return(list(
    correlations = cor_matrix,
    same_chr_cor = same_chr_cor,
    diff_chr_cor = diff_chr_cor,
    chromosome_test = chr_test,
    chr_stats = chr_stats,
    chr_info = chr_info,
    same_chr_mean = mean(abs(same_chr_cor), na.rm = TRUE),
    diff_chr_mean = mean(abs(diff_chr_cor), na.rm = TRUE)
  ))
}

#' Run main analysis with checkpoints
#' @param data_file Path to data file
#' @param gene_map_file Path to gene map file
#' @param batch_size Batch size for processing
#' @param checkpoint_dir Directory for checkpoints
#' @param start_from Starting point for analysis
run_analysis <- function(data_file, gene_map_file, batch_size = 100,
                         checkpoint_dir = "checkpoints",
                         start_from = "beginning") {
  
  # Step 1: Load data
  if(start_from == "beginning") {
    cat("Loading data...\n")
    data <- fread(data_file)
    gene_map <- fread(gene_map_file)
    
    metadata <- data[, 1:112]
    markers <- data[, 113:346]
    expression <- data[, 441:ncol(data)]
    
    split_data <- list(
      metadata = metadata,
      markers = markers,
      expression = expression,
      gene_map = gene_map
    )
    saveRDS(split_data, file.path(checkpoint_dir, "split_data.rds"))
  } else {
    split_data <- readRDS(file.path(checkpoint_dir, "split_data.rds"))
    metadata <- split_data$metadata
    markers <- split_data$markers
    expression <- split_data$expression
    gene_map <- split_data$gene_map
  }
  
  # Step 2: Process expression data
  if(start_from %in% c("beginning", "expression")) {
    cat("Processing expression data...\n")
    expr_data <- process_expression_data(expression, metadata)
    saveRDS(expr_data, file.path(checkpoint_dir, "processed_expression.rds"))
  } else {
    expr_data <- readRDS(file.path(checkpoint_dir, "processed_expression.rds"))
  }
  
  # Step 3: Perform eQTL analysis
  if(start_from %in% c("beginning", "expression", "eqtl")) {
    cat("Performing eQTL analysis...\n")
    eqtl_results <- perform_eqtl_analysis(expr_data, markers, metadata, batch_size)
    saveRDS(eqtl_results, file.path(checkpoint_dir, "eqtl_results.rds"))
  } else {
    eqtl_results <- readRDS(file.path(checkpoint_dir, "eqtl_results.rds"))
  }
  
  # Step 4: Analyze co-expression
  if(start_from %in% c("beginning", "expression", "eqtl", "coexpression")) {
    cat("Performing co-expression analysis...\n")
    coexp_results <- analyze_coexpression(expr_data, gene_map)
    saveRDS(coexp_results, file.path(checkpoint_dir, "coexp_results.rds"))
  } else {
    coexp_results <- readRDS(file.path(checkpoint_dir, "coexp_results.rds"))
  }
  
  return(list(
    eqtl = eqtl_results,
    coexpression = coexp_results,
    metadata = metadata,
    gene_map = gene_map,
    expr_data = expr_data
  ))
}