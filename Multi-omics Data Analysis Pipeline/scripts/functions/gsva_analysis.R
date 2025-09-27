prepare_module_expression_data <- function(module_genes_data, filtered_datasets, target_modules) {
  module_expression <- list()
  
  for (module in target_modules) {
    module_data <- module_genes_data[[module]]$gene_data
    
    module_expr <- list()
    
    for (dataset_name in names(filtered_datasets)) {
      dataset_samples <- rownames(filtered_datasets[[dataset_name]])
      expr_data <- as.matrix(module_data[, colnames(module_data) %in% dataset_samples])
      rownames(expr_data) <- module_data$gene_name
      
      module_expr[[dataset_name]] <- expr_data
    }
    
    module_expression[[module]] <- module_expr
  }
  
  return(module_expression)
}

combine_module_expression <- function(module_expression, target_modules) {
  combined_expression <- list()
  
  for (dataset_name in names(module_expression[[1]])) {
    dataset_expr_list <- list()
    
    for (module in target_modules) {
      dataset_expr_list[[module]] <- module_expression[[module]][[dataset_name]]
    }
    
    combined_expr <- do.call(rbind, dataset_expr_list)
    combined_expression[[dataset_name]] <- combined_expr
  }
  
  return(combined_expression)
}

filter_significant_genesets <- function(enrichment_results, geneset_databases, pvalue_threshold = 0.05) {
  significant_genesets <- unique(enrichment_results$ID[enrichment_results$p.adjust < pvalue_threshold])
  
  filtered_genesets <- list()
  
  for (db_name in names(geneset_databases)) {
    db_genesets <- geneset_databases[[db_name]]
    significant_in_db <- significant_genesets[significant_genesets %in% db_genesets$gs_name]
    
    if (length(significant_in_db) > 0) {
      filtered_db <- db_genesets[db_genesets$gs_name %in% significant_in_db, ]
      geneset_list <- split(filtered_db$gene_symbol, filtered_db$gs_name)
      filtered_genesets <- c(filtered_genesets, geneset_list)
    }
  }
  
  return(filtered_genesets)
}

perform_gsva_analysis <- function(expression_data, geneset_list, dataset_name) {
  library(GSVA)
  
  message(paste("Performing GSVA analysis for", dataset_name))
  
  gsva_result <- gsva(
    expr = expression_data,
    gset.idx.list = geneset_list,
    kcdf = GSVA_PARAMS$kcdf,
    min.sz = GSVA_PARAMS$min_sz,
    max.sz = GSVA_PARAMS$max_sz,
    parallel.sz = GSVA_PARAMS$parallel_cores,
    verbose = TRUE
  )
  
  return(gsva_result)
}

associate_gsva_with_features <- function(gsva_scores, clinical_features, dataset_name) {
  association_results <- list()
  
  for (feature_name in colnames(clinical_features)) {
    feature_vector <- clinical_features[[feature_name]]
    
    if (is.numeric(feature_vector)) {
      correlations <- apply(gsva_scores, 1, function(gsva_scores_row) {
        cor_result <- cor.test(gsva_scores_row, feature_vector, method = "spearman")
        return(c(correlation = cor_result$estimate, pvalue = cor_result$p.value))
      })
      
      correlations_df <- as.data.frame(t(correlations))
      correlations_df$pathway <- rownames(correlations_df)
      correlations_df$feature <- feature_name
      correlations_df$dataset <- dataset_name
      
      association_results[[feature_name]] <- correlations_df
    }
  }
  
  all_associations <- do.call(rbind, association_results)
  rownames(all_associations) <- NULL
  
  return(all_associations)
}

identify_significant_associations <- function(association_results, pvalue_threshold = 0.05) {
  significant_associations <- association_results[association_results$pvalue < pvalue_threshold, ]
  significant_associations <- significant_associations[order(significant_associations$pvalue), ]
  
  return(significant_associations)
}

plot_gsva_heatmap <- function(gsva_scores, clinical_features, dataset_name, output_dir) {
  library(pheatmap)
  
  significant_associations <- identify_significant_associations(
    associate_gsva_with_features(gsva_scores, clinical_features, dataset_name)
  )
  
  if (nrow(significant_associations) == 0) {
    message(paste("No significant associations found for", dataset_name))
    return(NULL)
  }
  
  top_pathways <- head(unique(significant_associations$pathway), 20)
  pathway_scores <- gsva_scores[top_pathways, ]
  
  annotation_df <- clinical_features[, sapply(clinical_features, is.numeric)]
  annotation_df <- annotation_df[, colSums(!is.na(annotation_df)) > 0]
  
  output_file <- file.path(output_dir, paste0("gsva_heatmap_", tolower(dataset_name), ".pdf"))
  
  pdf(output_file, width = 12, height = 8)
  
  pheatmap_result <- pheatmap(
    pathway_scores,
    annotation_col = annotation_df,
    scale = "row",
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = paste("GSVA Pathway Activity -", dataset_name)
  )
  
  dev.off()
  
  message(paste("GSVA heatmap saved:", output_file))
  
  return(pheatmap_result)
}

save_gsva_results <- function(gsva_results, association_results, output_dir) {
  for (dataset_name in names(gsva_results)) {
    gsva_scores <- gsva_results[[dataset_name]]
    associations <- association_results[[dataset_name]]
    
    gsva_file <- file.path(output_dir, paste0("gsva_scores_", tolower(dataset_name), ".csv"))
    assoc_file <- file.path(output_dir, paste0("gsva_associations_", tolower(dataset_name), ".csv"))
    
    write.csv(gsva_scores, gsva_file)
    write.csv(associations, assoc_file, row.names = FALSE)
  }
}

generate_gsva_report <- function(gsva_results, association_results) {
  report_file <- file.path(output_dirs$rdata, "gsva_analysis_report.txt")
  
  sink(report_file)
  
  cat("GSVA Analysis Report\n")
  cat("===================\n\n")
  cat("Generated on:", date(), "\n\n")
  
  cat("GSVA Scores Summary:\n")
  cat("-------------------\n")
  
  for (dataset_name in names(gsva_results)) {
    gsva_data <- gsva_results[[dataset_name]]
    
    cat(paste("\nDataset:", dataset_name, "\n"))
    cat(paste("Number of pathways:", nrow(gsva_data), "\n"))
    cat(paste("Number of samples:", ncol(gsva_data), "\n"))
    cat(paste("Score range:", round(min(gsva_data), 3), "to", round(max(gsva_data), 3), "\n"))
  }
  
  cat("\nSignificant Associations:\n")
  cat("------------------------\n")
  
  for (dataset_name in names(association_results)) {
    associations <- association_results[[dataset_name]]
    significant <- associations[associations$pvalue < GSVA_PARAMS$significance_threshold, ]
    
    cat(paste("\nDataset:", dataset_name, "\n"))
    cat(paste("Total significant associations:", nrow(significant), "\n"))
    
    if (nrow(significant) > 0) {
      top_associations <- head(significant[order(significant$pvalue), ], 5)
      cat("Top associations:\n")
      print(top_associations[, c("pathway", "feature", "correlation", "pvalue")])
    }
  }
  
  sink()
  
  message(paste("GSVA report saved:", report_file))
}