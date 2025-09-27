library(WGCNA)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)

source("config/paths.R")
source("config/parameters.R")
source("scripts/functions/enrichment_analysis.R")

perform_enrichment_analysis <- function() {
  message("Starting enrichment analysis...")
  
  load(file.path(output_dirs$rdata, "network_results.RData"))
  load(file.path(output_dirs$rdata, "processed_data.RData"))
  
  geneset_databases <- prepare_geneset_databases()
  
  module_genes <- extract_module_genes_results()
  
  correlated_genes <- extract_correlated_genes_all_datasets(network_results, filtered_datasets)
  
  enrichment_results <- perform_all_enrichment_analyses(correlated_genes, geneset_databases)
  
  combined_results <- combine_all_enrichment_results(enrichment_results)
  
  reproducibility_summary <- assess_reproducibility(combined_results, 
                                                    names(filtered_datasets))
  
  save_enrichment_results(combined_results, output_dirs$tables)
  save_geneset_databases(geneset_databases)
  generate_enrichment_report(combined_results, reproducibility_summary)
  
  save(module_genes, correlated_genes, enrichment_results, combined_results,
       reproducibility_summary, file = file.path(output_dirs$rdata, "enrichment_results.RData"))
  
  message("Enrichment analysis completed successfully!")
  
  return(list(
    module_genes = module_genes,
    correlated_genes = correlated_genes,
    enrichment_results = enrichment_results,
    combined_results = combined_results,
    reproducibility_summary = reproducibility_summary
  ))
}
save_geneset_databases <- function(geneset_databases) {
  for (db_name in names(geneset_databases)) {
    assign(paste0(tolower(db_name), "_dt"), geneset_databases[[db_name]])
    save_file <- file.path(output_dirs$rdata, paste0(tolower(db_name), "_geneset.RData"))
    save(list = paste0(tolower(db_name), "_dt"), file = save_file)
  }
}
extract_module_genes_results <- function() {
  load(file.path(output_dirs$rdata, "module_preservation_results.RData"))
  
  module_genes <- list()
  
  for (module in ENRICHMENT_PARAMS$target_modules) {
    module_data <- preservation_results$batch$preservation_result$multiColor$Reference
    module_genes[[module]] <- names(module_data)[module_data == module]
  }
  
  return(module_genes)
}

extract_correlated_genes_all_datasets <- function(network_results, filtered_datasets) {
  correlated_genes <- list()
  
  for (dataset_name in names(filtered_datasets)) {
    dataset_correlated <- list()
    
    for (module in ENRICHMENT_PARAMS$target_modules) {
      module_eigengenes <- network_results[[dataset_name]]$module_trait$module_eigengenes
      module_colors <- network_results[[dataset_name]]$network$colors
      
      correlated <- extract_correlated_genes(
        filtered_datasets[[dataset_name]],
        module_eigengenes,
        module_colors,
        module
      )
      
      dataset_correlated[[module]] <- correlated
    }
    
    correlated_genes[[dataset_name]] <- dataset_correlated
  }
  
  return(correlated_genes)
}

perform_all_enrichment_analyses <- function(correlated_genes, geneset_databases) {
  enrichment_results <- list()
  
  for (dataset_name in names(correlated_genes)) {
    dataset_enrichment <- list()
    
    for (module in ENRICHMENT_PARAMS$target_modules) {
      module_correlated <- correlated_genes[[dataset_name]][[module]]
      
      for (cor_type in c("positive", "negative")) {
        genes <- module_correlated[[cor_type]]
        
        if (length(genes) > 0) {
          enrichment <- perform_enrichment_analysis(genes, geneset_databases)
          
          for (db_name in names(enrichment)) {
            result_name <- paste(db_name, "M", 
                                 gsub("green", "1", gsub("yellow", "2", module)),
                                 toupper(substr(cor_type, 1, 1)), "cor",
                                 dataset_name, sep = "_")
            
            dataset_enrichment[[result_name]] <- enrichment[[db_name]]
          }
        }
      }
    }
    
    enrichment_results[[dataset_name]] <- dataset_enrichment
  }
  
  return(enrichment_results)
}

combine_all_enrichment_results <- function(enrichment_results) {
  combined_results <- list()
  
  for (dataset_name in names(enrichment_results)) {
    dataset_combined <- combine_enrichment_results(enrichment_results[[dataset_name]], 
                                                   dataset_name)
    combined_results[[dataset_name]] <- dataset_combined
  }
  
  all_combined <- do.call(rbind, combined_results)
  rownames(all_combined) <- NULL
  
  return(all_combined)
}

if (!interactive()) {
  enrichment_results <- perform_enrichment_analysis()
}