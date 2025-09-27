library(WGCNA)
library(tidyverse)

source("config/paths.R")
source("config/parameters.R")
source("scripts/functions/network_analysis.R")
source("scripts/functions/visualization.R")

perform_network_analysis <- function() {
  load(file.path(output_dirs$rdata, "processed_data.RData"))
  
  trait_data <- simulate_trait_data(filtered_datasets$batch)
  
  results <- list()
  
  for (dataset_name in names(filtered_datasets)) {
    expression_data <- filtered_datasets[[dataset_name]]
    traits <- trait_data[[dataset_name]]
    
    plot_sample_dendrogram(expression_data, traits, dataset_name, 
                           output_dirs$figures, ANALYSIS_PARAMS$cluster_cut_height)
    
    threshold_result <- select_soft_threshold(
      expression_data, 
      ANALYSIS_PARAMS$network_type,
      ANALYSIS_PARAMS$r_sq_cutoff
    )
    
    plot_soft_threshold(threshold_result$sft, output_dirs$figures, 
                        ANALYSIS_PARAMS$r_sq_cutoff)
    
    network <- build_coexpression_network(
      expression_data,
      threshold_result$power,
      ANALYSIS_PARAMS$min_module_size,
      ANALYSIS_PARAMS$merge_cut_height,
      ANALYSIS_PARAMS$network_type
    )
    
    plot_module_dendrogram(network, output_dirs$figures)
    
    module_trait_results <- plot_module_trait_heatmap(
      expression_data, traits, network$colors,
      dataset_name, output_dirs$figures
    )
    
    results[[dataset_name]] <- list(
      soft_threshold = threshold_result,
      network = network,
      module_trait = module_trait_results
    )
  }
  
  save(results, file = file.path(output_dirs$rdata, "network_results.RData"))
  
  return(results)
}

if (!interactive()) {
  network_results <- perform_network_analysis()
}