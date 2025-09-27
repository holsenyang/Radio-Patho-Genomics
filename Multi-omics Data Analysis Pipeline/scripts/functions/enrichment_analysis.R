extract_module_genes <- function(module_colors_data, gene_mapping, target_modules) {
  module_genes <- list()
  
  for (module in target_modules) {
    module_data <- module_colors_data[module_colors_data$mergeModule == module, ]
    module_genes[[module]] <- list(
      all_genes = module_data$gene_name,
      gene_data = module_data
    )
  }
  
  return(module_genes)
}

prepare_geneset_databases <- function(species = "Homo sapiens") {
  library(msigdbr)
  
  genesets <- list()
  
  if ("GO" %in% ENRICHMENT_PARAMS$databases) {
    go_df <- msigdbr(species = species, category = "C5")
    go_df <- go_df[go_df$gs_subcat != "HPO", ]
    genesets$GO <- dplyr::select(go_df, gs_name, gene_symbol)
  }
  
  if ("KEGG" %in% ENRICHMENT_PARAMS$databases) {
    kegg_df <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG")
    genesets$KEGG <- dplyr::select(kegg_df, gs_name, gene_symbol)
  }
  
  if ("PID" %in% ENRICHMENT_PARAMS$databases) {
    pid_df <- msigdbr(species = species, category = "C2", subcategory = "CP:PID")
    genesets$PID <- dplyr::select(pid_df, gs_name, gene_symbol)
  }
  
  if ("REACTOME" %in% ENRICHMENT_PARAMS$databases) {
    reactome_df <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
    genesets$REACTOME <- dplyr::select(reactome_df, gs_name, gene_symbol)
  }
  
  if ("WP" %in% ENRICHMENT_PARAMS$databases) {
    wp_df <- msigdbr(species = species, category = "C2", subcategory = "CP:WIKIPATHWAYS")
    genesets$WP <- dplyr::select(wp_df, gs_name, gene_symbol)
  }
  
  if ("BIOCARTA" %in% ENRICHMENT_PARAMS$databases) {
    biocarta_df <- msigdbr(species = species, category = "C2", subcategory = "CP:BIOCARTA")
    genesets$BIOCARTA <- dplyr::select(biocarta_df, gs_name, gene_symbol)
  }
  
  if ("HALLMARK" %in% ENRICHMENT_PARAMS$databases) {
    hallmark_df <- msigdbr(species = species, category = "H")
    genesets$HALLMARK <- dplyr::select(hallmark_df, gs_name, gene_symbol)
  }
  
  return(genesets)
}

perform_enrichment_analysis <- function(gene_list, geneset_database) {
  library(clusterProfiler)
  
  enrichment_results <- list()
  
  for (db_name in names(geneset_database)) {
    term2gene <- geneset_database[[db_name]]
    
    enrich_result <- enricher(
      gene_list,
      pAdjustMethod = ENRICHMENT_PARAMS$p_adjust_method,
      pvalueCutoff = ENRICHMENT_PARAMS$pvalue_cutoff,
      qvalueCutoff = ENRICHMENT_PARAMS$qvalue_cutoff,
      TERM2GENE = term2gene
    )
    
    enrichment_results[[db_name]] <- enrich_result
  }
  
  return(enrichment_results)
}

extract_correlated_genes <- function(expression_data, module_eigengenes, module_colors, 
                                     target_module, correlation_threshold = 0) {
  
  module_genes <- names(module_colors)[module_colors == target_module]
  module_eigengene_name <- paste0("ME", target_module)
  
  if (!module_eigengene_name %in% colnames(module_eigengenes)) {
    return(list(positive = character(), negative = character()))
  }
  
  correlations <- cor(expression_data[, module_genes], 
                      module_eigengenes[, module_eigengene_name])
  
  positive_genes <- rownames(correlations)[correlations > correlation_threshold]
  negative_genes <- rownames(correlations)[correlations < -correlation_threshold]
  
  return(list(positive = positive_genes, negative = negative_genes))
}

combine_enrichment_results <- function(enrichment_list, dataset_name) {
  combined_results <- data.frame()
  
  for (i in seq_along(enrichment_list)) {
    enrich_name <- names(enrichment_list)[i]
    enrich_result <- enrichment_list[[i]]
    
    if (!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
      result_df <- enrich_result@result
      result_df$dataset <- dataset_name
      result_df$database <- gsub("_.*", "", enrich_name)
      result_df$module <- gsub(".*(M[0-9]).*", "\\1", enrich_name)
      result_df$correlation <- gsub(".*M[0-9]([A-Z]+)_.*", "\\1", enrich_name)
      
      combined_results <- rbind(combined_results, result_df)
    }
  }
  
  return(combined_results)
}

assess_reproducibility <- function(results_list, datasets) {
  reproducibility_summary <- list()
  
  for (module in ENRICHMENT_PARAMS$target_modules) {
    module_results <- list()
    
    for (dataset in datasets) {
      dataset_results <- results_list[[dataset]]
      module_data <- dataset_results[dataset_results$module == module, ]
      module_results[[dataset]] <- unique(module_data$ID)
    }
    
    common_pathways <- Reduce(intersect, module_results)
    reproducibility_summary[[module]] <- list(
      total_common = length(common_pathways),
      common_pathways = common_pathways,
      dataset_specific = lapply(module_results, function(x) setdiff(x, common_pathways))
    )
  }
  
  return(reproducibility_summary)
}

save_enrichment_results <- function(enrichment_results, output_dir) {
  for (dataset_name in names(enrichment_results)) {
    dataset_results <- enrichment_results[[dataset_name]]
    
    output_file <- file.path(output_dir, 
                             paste0("enrichment_results_", tolower(dataset_name), ".csv"))
    write.csv(dataset_results, output_file, row.names = FALSE)
  }
}

generate_enrichment_report <- function(enrichment_results, reproducibility_summary) {
  report_file <- file.path(output_dirs$rdata, "enrichment_analysis_report.txt")
  
  sink(report_file)
  
  cat("Enrichment Analysis Report\n")
  cat("==========================\n\n")
  cat("Generated on:", date(), "\n\n")
  
  cat("Summary Statistics:\n")
  cat("-------------------\n")
  
  for (dataset_name in names(enrichment_results)) {
    dataset_data <- enrichment_results[[dataset_name]]
    
    cat(paste("\nDataset:", dataset_name, "\n"))
    cat(paste("Total enriched pathways:", nrow(dataset_data), "\n"))
    
    for (module in ENRICHMENT_PARAMS$target_modules) {
      module_data <- dataset_data[dataset_data$module == module, ]
      cat(paste("  Module", module, ":", nrow(module_data), "pathways\n"))
    }
  }
  
  cat("\nReproducibility Analysis:\n")
  cat("-------------------------\n")
  
  for (module in names(reproducibility_summary)) {
    module_summary <- reproducibility_summary[[module]]
    
    cat(paste("\nModule:", module, "\n"))
    cat(paste("Common pathways across datasets:", module_summary$total_common, "\n"))
    
    for (dataset_name in names(module_summary$dataset_specific)) {
      specific_count <- length(module_summary$dataset_specific[[dataset_name]])
      cat(paste("  ", dataset_name, "specific pathways:", specific_count, "\n"))
    }
  }
  
  sink()
}