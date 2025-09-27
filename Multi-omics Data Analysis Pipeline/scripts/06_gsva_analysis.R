library(GSVA)
library(pheatmap)
library(tidyverse)

source("config/paths.R")
source("config/parameters.R")
source("scripts/functions/gsva_analysis.R")

perform_gsva_analysis <- function() {
  message("Starting GSVA analysis...")
  
  load(file.path(output_dirs$rdata, "enrichment_results.RData"))
  load(file.path(output_dirs$rdata, "processed_data.RData"))
  
  significant_genesets <- extract_significant_genesets(combined_results)
  
  module_expression <- prepare_module_expression_data(module_genes, filtered_datasets, 
                                                      ENRICHMENT_PARAMS$target_modules)
  
  combined_expression <- combine_module_expression(module_expression, 
                                                   ENRICHMENT_PARAMS$target_modules)
  
  gsva_results <- perform_gsva_all_datasets(combined_expression, significant_genesets)
  
  clinical_features <- prepare_clinical_features()
  
  association_results <- analyze_gsva_associations(gsva_results, clinical_features)
  
  plot_gsva_results(gsva_results, clinical_features)
  
  save_gsva_results(gsva_results, association_results, output_dirs$tables)
  
  generate_gsva_report(gsva_results, association_results)
  
  save(gsva_results, association_results, significant_genesets,
       file = file.path(output_dirs$rdata, "gsva_results.RData"))
  
  message("GSVA analysis completed successfully!")
  
  return(list(
    gsva_results = gsva_results,
    association_results = association_results,
    significant_genesets = significant_genesets
  ))
}

extract_significant_genesets <- function(enrichment_results) {
  significant_ids <- unique(enrichment_results$ID[enrichment_results$p.adjust < 0.05])
  
  load(file.path(output_dirs$rdata, "enrichment_results.RData"))
  
  geneset_databases <- list()
  for (db_name in ENRICHMENT_PARAMS$databases) {
    db_file <- file.path(output_dirs$rdata, paste0(tolower(db_name), "_geneset.RData"))
    if (file.exists(db_file)) {
      load(db_file)
      geneset_databases[[db_name]] <- get(paste0(tolower(db_name), "_dt"))
    }
  }
  
  significant_genesets <- list()
  for (db_name in names(geneset_databases)) {
    db_data <- geneset_databases[[db_name]]
    significant_in_db <- significant_ids[significant_ids %in% db_data$gs_name]
    
    if (length(significant_in_db) > 0) {
      geneset_list <- split(db_data$gene_symbol, db_data$gs_name)
      significant_genesets <- c(significant_genesets, geneset_list[significant_in_db])
    }
  }
  
  return(significant_genesets)
}

prepare_clinical_features <- function() {
  clinical_features <- list()
  
  for (dataset_name in DATASET_CONFIG) {
    features <- data.frame(
      sample_id = rownames(filtered_datasets[[dataset_name]]),
      simulated_feature1 = rnorm(nrow(filtered_datasets[[dataset_name]])),
      simulated_feature2 = sample(c("A", "B", "C"), nrow(filtered_datasets[[dataset_name]]), replace = TRUE),
      simulated_feature3 = runif(nrow(filtered_datasets[[dataset_name]]))
    )
    rownames(features) <- features$sample_id
    clinical_features[[dataset_name]] <- features
  }
  
  return(clinical_features)
}

perform_gsva_all_datasets <- function(combined_expression, significant_genesets) {
  gsva_results <- list()
  
  for (dataset_name in names(combined_expression)) {
    expr_data <- combined_expression[[dataset_name]]
    
    if (nrow(expr_data) > 0 && length(significant_genesets) > 0) {
      gsva_result <- perform_gsva_analysis(expr_data, significant_genesets, dataset_name)
      gsva_results[[dataset_name]] <- gsva_result
    } else {
      warning(paste("Skipping GSVA for", dataset_name, "- no expression data or genesets"))
    }
  }
  
  return(gsva_results)
}

analyze_gsva_associations <- function(gsva_results, clinical_features) {
  association_results <- list()
  
  for (dataset_name in names(gsva_results)) {
    gsva_scores <- gsva_results[[dataset_name]]
    clinical_data <- clinical_features[[dataset_name]]
    
    numeric_features <- clinical_data[, sapply(clinical_data, is.numeric)]
    
    if (ncol(numeric_features) > 0) {
      associations <- associate_gsva_with_features(gsva_scores, numeric_features, dataset_name)
      association_results[[dataset_name]] <- associations
    }
  }
  
  return(association_results)
}

plot_gsva_results <- function(gsva_results, clinical_features) {
  for (dataset_name in names(gsva_results)) {
    gsva_scores <- gsva_results[[dataset_name]]
    clinical_data <- clinical_features[[dataset_name]]
    
    numeric_features <- clinical_data[, sapply(clinical_data, is.numeric)]
    
    if (ncol(numeric_features) > 0) {
      plot_gsva_heatmap(gsva_scores, numeric_features, dataset_name, output_dirs$figures)
    }
  }
}

if (!interactive()) {
  gsva_results <- perform_gsva_analysis()
}