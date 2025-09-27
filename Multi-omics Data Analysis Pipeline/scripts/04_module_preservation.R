library(WGCNA)
library(tidyverse)

source("config/paths.R")
source("config/parameters.R")
source("scripts/functions/module_preservation.R")

perform_preservation_analysis <- function() {
  load(file.path(output_dirs$rdata, "processed_data.RData"))
  load(file.path(output_dirs$rdata, "network_results.RData"))
  
  module_colors <- network_results$batch$network$colors
  
  preservation_results <- list()
  
  test_datasets <- list(
    tcga = list(name = "TCGA", expression = "exp_dat_TCGA18"),
    cptac = list(name = "CPTAC", expression = "exp_dat_CPTAC")
  )
  
  for (dataset_key in names(test_datasets)) {
    dataset_config <- test_datasets[[dataset_key]]
    
    dataset_result <- analyze_dataset_preservation(
      dataset_config, module_colors
    )
    
    preservation_results[[dataset_key]] <- dataset_result
  }
  
  save(preservation_results, 
       file = file.path(output_dirs$rdata, "preservation_results.RData"))
  
  generate_preservation_report(preservation_results)
  
  return(preservation_results)
}

if (!interactive()) {
  preservation_results <- perform_preservation_analysis()
}