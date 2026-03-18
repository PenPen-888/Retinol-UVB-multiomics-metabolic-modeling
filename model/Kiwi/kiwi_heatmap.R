library(ComplexHeatmap)
library(circlize)
library(grid)
library(stringr)

kiwi_heatmap <- function(df,
                         title = '',
                         width = 12,
                         height = 8,
                         save_prefix = 'kiwi_plot',
                         output_dir = 'kiwi_output',
                         filter_gene = FALSE,
                         filter_metabo = FALSE,            
                         filter_method = "combined",
                         top_n = 50,
                         top_n_metab = 50,                
                         save_score_table = FALSE
) {
  

  rownames(df) <- df[, 1]
  df_clean <- df[, -1]
  df_clean <- as.data.frame(lapply(df_clean, as.numeric), row.names = rownames(df_clean))
  df_clean[is.na(df_clean)] <- 0

  if (filter_metabo) {
    metab_nonzero_freq <- apply(df_clean, 1, function(x) sum(x != 0) / length(x))
    metab_max_score <- apply(df_clean, 1, function(x) max(abs(x)))
    metab_mean_score <- apply(df_clean, 1, function(x) mean(abs(x)))
    
    freq_scaled <- (metab_nonzero_freq - min(metab_nonzero_freq)) / (max(metab_nonzero_freq) - min(metab_nonzero_freq))
    max_scaled <- (metab_max_score - min(metab_max_score)) / (max(metab_max_score) - min(metab_max_score))
    mean_scaled <- (metab_mean_score - min(metab_mean_score)) / (max(metab_mean_score) - min(metab_mean_score))
    
    if (filter_method == "combined") {
      final_score <- 0.5 * freq_scaled + 0.5 * max_scaled
    } else if (filter_method == "max") {
      final_score <- max_scaled
    } else if (filter_method == "mean") {
      final_score <- mean_scaled
    } else {
      stop("filter_method only 'max'、'mean' 或 'combined'")
    }
    
    top_metabs <- names(sort(final_score, decreasing = TRUE))[1:min(top_n_metab, length(final_score))]
    df_clean <- df_clean[top_metabs, , drop = FALSE]
  
  }

  if (filter_gene) {
    gene_nonzero_freq <- apply(df_clean, 2, function(x) sum(x != 0) / length(x))
    gene_max_score <- apply(df_clean, 2, function(x) max(abs(x)))
    gene_mean_score <- apply(df_clean, 2, function(x) mean(abs(x)))
    
    freq_scaled <- (gene_nonzero_freq - min(gene_nonzero_freq)) / (max(gene_nonzero_freq) - min(gene_nonzero_freq))
    max_scaled <- (gene_max_score - min(gene_max_score)) / (max(gene_max_score) - min(gene_max_score))
    mean_scaled <- (gene_mean_score - min(gene_mean_score)) / (max(gene_mean_score) - min(gene_mean_score))
    
    if (filter_method == "combined") {
      final_score <- 0.5 * freq_scaled + 0.5 * max_scaled
    } else if (filter_method == "max") {
      final_score <- max_scaled
    } else if (filter_method == "mean") {
      final_score <- mean_scaled
    } else {
      stop("filter_method only 'max'、'mean' 或 'combined'")
    }
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    if (save_score_table) {
      score_table <- data.frame(
        Gene = names(final_score),
        NonzeroFreq = gene_nonzero_freq[names(final_score)],
        MaxAbsScore = gene_max_score[names(final_score)],
        MeanAbsScore = gene_mean_score[names(final_score)],
        ScaledFreq = freq_scaled[names(final_score)],
        ScaledMax = max_scaled[names(final_score)],
        ScaledMean = mean_scaled[names(final_score)],
        FinalScore = final_score[names(final_score)]
      )
      write.csv(score_table, file = file.path(output_dir, paste0(save_prefix, "_score_table.csv")), row.names = FALSE)
    }
    
    top_genes <- names(sort(final_score, decreasing = TRUE))[1:min(top_n, length(final_score))]
    df_clean <- df_clean[, top_genes, drop = FALSE]
    title <- paste(title, "(Top", top_n, "Genes)", sep = " ")
  }
  
 
  

  rownames(df_clean) <- stringr::str_wrap(rownames(df_clean), width = 30)
  

  max_abs <- max(abs(df_clean), na.rm = TRUE)
  
  col_fun <- circlize::colorRamp2(
    c(-max_abs, 0, max_abs),
    c("blue", "white", "red")
  )
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  ht <- ComplexHeatmap::Heatmap(
    df_clean,
    name = "Score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    show_row_dend = TRUE,
    heatmap_legend_param = list(
      title = "Score",
      title_position = "topcenter",
      legend_direction = "vertical"
    ),
    column_names_side = "bottom",
    row_names_side = "left",
    column_names_gp = grid::gpar(fontsize = 14),
    row_names_gp = grid::gpar(fontsize = 14),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 14)
  )
  
  png(file.path(output_dir, paste0(save_prefix, ".png")), width = width, height = height, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  pdf(file.path(output_dir, paste0(save_prefix, ".pdf")), width = width, height = height)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  return(df_clean)
}
