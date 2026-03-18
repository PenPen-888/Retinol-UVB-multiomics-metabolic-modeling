library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

fc_function <- function(df, fc, p = NULL, padj = NULL, kegg = TRUE, GO = TRUE) {
  
  df_name <- deparse(substitute(df))
  folder_name <- paste0(df_name, "_fc", fc)
  dir.create(folder_name, showWarnings = FALSE)
  

  df_filtered <- df %>%
    filter(2^logFC >= fc | 2^logFC <= 1/fc)
  
  if (!is.null(p)) {
    df_filtered <- df_filtered %>% filter(pvalue <= p)
  }
  if (!is.null(padj)) {
    df_filtered <- df_filtered %>% filter(padj <= padj)
  }
  
  write.csv(df_filtered, file = file.path(folder_name, "filtered_diff_genes.csv"), row.names = FALSE)
  

  all_genes <- na.omit(df_filtered$geneID)
  up_genes <- df_filtered %>% filter(logFC > 0) %>% pull(geneID) %>% na.omit()
  down_genes <- df_filtered %>% filter(logFC < 0) %>% pull(geneID) %>% na.omit()
  

  convert_ids_to_symbols <- function(enrich_df) {
    enrich_df <- enrich_df %>%
      mutate(geneSymbol = sapply(strsplit(geneID, "/"), function(ids) {
        symbols <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
        matched_symbols <- symbols$SYMBOL[match(ids, symbols$ENTREZID)]
        matched_symbols[is.na(matched_symbols)] <- ids[is.na(matched_symbols)]  # fallback
        paste(matched_symbols, collapse = "/")
      }))
    return(enrich_df)
  }
  

  enrich_and_save <- function(gene_set, type = c("ALL", "UP", "DOWN")) {
    if (length(gene_set) < 10) return(NULL)
    

    id_symbol_map <- df %>% dplyr::select(geneID, symbol) %>% distinct() %>% na.omit()
    

    process_result <- function(result_df) {
      result_df <- result_df %>%
        mutate(GeneRatio = as.character(GeneRatio)) %>%
        mutate(symbol_list = sapply(strsplit(geneID, "/"), function(ids) {
          mapped <- id_symbol_map$symbol[match(ids, id_symbol_map$geneID)]
          paste(na.omit(mapped), collapse = "/")
        }))
      return(result_df)
    }
    
    # KEGG
    if (kegg) {
      kegg_result <- enrichKEGG(
        gene = gene_set,
        organism = "hsa",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
      if (!is.null(kegg_result)) {
        kegg_df <- as.data.frame(kegg_result) |> process_result()
        write.csv(kegg_df, file = file.path(folder_name, paste0("KEGG_", type, "_all.csv")), row.names = FALSE)
        sig_kegg <- kegg_df %>% filter(pvalue < 0.05)
        write.csv(sig_kegg, file = file.path(folder_name, paste0("KEGG_", type, "_p0.05.csv")), row.names = FALSE)
      }
    }
    
    # GO
    if (GO) {
      go_result <- enrichGO(
        gene = gene_set,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "ALL",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
      if (!is.null(go_result)) {
        go_df <- as.data.frame(go_result) |> process_result()
        write.csv(go_df, file = file.path(folder_name, paste0("GO_", type, "_all.csv")), row.names = FALSE)
        sig_go <- go_df %>% filter(pvalue < 0.05)
        write.csv(sig_go, file = file.path(folder_name, paste0("GO_", type, "_p0.05.csv")), row.names = FALSE)
      }
    }
  }
  
  

  enrich_and_save(all_genes, "ALL")
  enrich_and_save(up_genes, "UP")
  enrich_and_save(down_genes, "DOWN")
  
  message("✅ Saved in", folder_name)
}