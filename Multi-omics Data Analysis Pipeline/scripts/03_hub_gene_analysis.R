library(WGCNA)
library(tidyverse)
library(ggrepel)

source("config/paths.R")
source("config/parameters.R")
source("scripts/functions/hub_gene_analysis.R")

perform_hub_gene_analysis <- function() {
  load(file.path(output_dirs$rdata, "network_results.RData"))
  load(file.path(output_dirs$rdata, "processed_data.RData"))
  
  target_modules <- c("green", "yellow")
  all_results <- list()
  hub_gene_summary <- data.frame()
  
  for (module in target_modules) {
    module_results <- analyze_module_across_datasets(
      module, network_results, filtered_datasets
    )
    
    all_results[[module]] <- module_results
    
    if (!is.null(module_results$hub_gene_summary)) {
      module_results$hub_gene_summary$module <- module
      hub_gene_summary <- rbind(hub_gene_summary, module_results$hub_gene_summary)
    }
  }
  
  save(all_results, hub_gene_summary,
       file = file.path(output_dirs$rdata, "hub_gene_results.RData"))
  
  generate_hub_gene_report(all_results, hub_gene_summary)
  
  return(list(all_results = all_results, hub_gene_summary = hub_gene_summary))
}

if (!interactive()) {
  hub_gene_results <- perform_hub_gene_analysis()
}