select_soft_threshold <- function(expression_data, 
                                  network_type = "unsigned",
                                  r_sq_cutoff = 0.85,
                                  power_range = c(seq(1, 10, by = 1), 
                                                  seq(12, 21, by = 2))) {
  
  sft <- pickSoftThreshold(
    expression_data,
    networkType = network_type,
    powerVector = power_range,
    RsquaredCut = r_sq_cutoff,
    verbose = 5
  )
  
  power_estimate <- sft$powerEstimate
  
  if (is.na(power_estimate)) {
    n_samples <- nrow(expression_data)
    
    power_estimate <- ifelse(n_samples < 20, 9,
                             ifelse(n_samples < 30, 8,
                                    ifelse(n_samples < 40, 7, 6)))
    
    if (network_type != "unsigned") {
      power_estimate <- power_estimate * 2
    }
  }
  
  return(list(
    sft = sft,
    power = power_estimate
  ))
}

build_coexpression_network <- function(expression_data, 
                                       power,
                                       min_module_size = 30,
                                       merge_cut_height = 0.35,
                                       network_type = "unsigned") {
  
  net <- blockwiseModules(
    expression_data,
    power = power,
    maxBlockSize = ncol(expression_data),
    corType = "pearson",
    networkType = network_type,
    TOMType = network_type,
    minModuleSize = min_module_size,
    mergeCutHeight = merge_cut_height,
    numericLabels = FALSE,
    saveTOMs = FALSE,
    verbose = 3
  )
  
  module_table <- table(net$colors)
  
  return(net)
}