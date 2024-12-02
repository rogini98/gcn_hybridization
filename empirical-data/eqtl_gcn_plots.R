# Required packages
library(ggplot2)
library(reshape2)
library(viridis)

#' 1. eQTL Manhattan Plot and Heatmap
#' @param eqtl_results Data frame with columns: marker_pos, gene_pos, p_value
plot_eqtl_results <- function(eqtl_results) {
  # Manhattan plot: Shows significance of all eQTL associations
  png("eqtl_manhattan.png", width = 1000, height = 600)
  ggplot(eqtl_results, aes(x = marker_pos, y = -log10(p_value))) +
    geom_point(alpha = 0.6, size = 0.8) +
    facet_grid(~chromosome, scales = "free_x", space = "free_x") +
    theme_minimal() +
    labs(x = "Marker Position", y = "-log10(p-value)", 
         title = "Manhattan Plot of eQTL Associations")
  dev.off()
  
  # Cis/Trans Plot: Shows relationship between gene and marker positions
  png("eqtl_cistrans.png", width = 800, height = 800)
  ggplot(eqtl_results[eqtl_results$p_value < 0.05,], 
         aes(x = marker_pos, y = gene_pos)) +
    geom_point(aes(color = -log10(p_value)), alpha = 0.7) +
    scale_color_viridis() +
    theme_minimal() +
    labs(x = "Marker Position", y = "Gene Position",
         title = "Cis/Trans eQTL Plot",
         color = "-log10(p-value)")
  dev.off()
}

#' 2. Gene Co-expression Network and Heatmap
#' @param cor_matrix Correlation matrix of gene expression
#' @param gene_info Data frame with gene information
plot_coexpression <- function(cor_matrix, gene_info) {
  # Correlation Heatmap
  png("coexpression_heatmap.png", width = 800, height = 800)
  heatmap(cor_matrix,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          main = "Gene Co-expression Heatmap")
  dev.off()
  
  # Network visualization (for highly correlated genes only)
  high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
  network_data <- data.frame(
    gene1 = rownames(cor_matrix)[high_cor[,1]],
    gene2 = colnames(cor_matrix)[high_cor[,2]],
    correlation = cor_matrix[high_cor]
  )
  
  png("coexpression_network.png", width = 1000, height = 1000)
  ggplot(network_data, aes(x = gene1, y = gene2)) +
    geom_tile(aes(fill = correlation)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Gene Co-expression Network")
  dev.off()
}

#' 3. Chromosome Analysis Visualization
#' @param chr_correlations List with same/different chromosome correlation results
plot_chromosome_analysis <- function(chr_correlations) {
  # Boxplot comparing within/between chromosome correlations
  png("chromosome_correlations.png", width = 600, height = 800)
  boxplot(list(
    "Same Chromosome" = chr_correlations$same_chr_cor,
    "Different Chromosomes" = chr_correlations$diff_chr_cor
  ),
  main = "Distribution of Gene Correlations by Chromosomal Location",
  ylab = "Absolute Correlation",
  col = c("lightblue", "lightgreen"))
  dev.off()
  
  # Chromosome-wise correlation heatmap
  png("chromosome_heatmap.png", width = 800, height = 800)
  heatmap(chr_correlations$chr_cor_matrix,
          col = viridis(100),
          main = "Chromosome-wise Correlation Heatmap")
  dev.off()
}