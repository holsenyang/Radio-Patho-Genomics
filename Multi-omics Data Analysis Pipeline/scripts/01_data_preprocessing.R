library(WGCNA)
library(tidyverse)
library(impute)

source("config/paths.R")
source("config/parameters.R")
source("scripts/functions/data_processing.R")

set.seed(12345)
enableWGCNAThreads()

perform_data_preprocessing <- function() {
  load("data/raw/expression_data.RData")
  
  processed_data <- preprocess_expression_data(expression_data)
  
  datasets <- prepare_datasets(processed_data, DATASET_CONFIG)
  
  filtered_datasets <- perform_quality_control(datasets)
  
  common_genes <- find_common_genes(filtered_datasets)
  
  for (name in names(filtered_datasets)) {
    filtered_datasets[[name]] <- filtered_datasets[[name]][, common_genes]
  }
  
  highly_variable_genes <- select_variable_genes(filtered_datasets, ANALYSIS_PARAMS$top_genes)
  
  gene_id_name <- create_gene_mapping(common_genes)
  
  save(filtered_datasets, highly_variable_genes, gene_id_name,
       file = file.path(output_dirs$rdata, "processed_data.RData"))
  
  return(list(filtered_datasets = filtered_datasets, 
              highly_variable_genes = highly_variable_genes,
              gene_id_name = gene_id_name))
}

create_gene_mapping <- function(gene_ids) {
  gene_mapping <- data.frame(
    gene_id = gene_ids,
    gene_name = paste0("GENE_", seq_along(gene_ids)),
    stringsAsFactors = FALSE
  )
  
  return(gene_mapping)
}

if (!interactive()) {
  results <- perform_data_preprocessing()
}