preprocess_expression_data <- function(raw_expression) {
  gene_names <- raw_expression$gene_id
  expression_values <- raw_expression[, -1]
  
  imputed_data <- impute.knn(
    as.matrix(expression_values),
    k = ANALYSIS_PARAMS$impute_k,
    rowmax = 0.5,
    colmax = 0.8,
    rng.seed = 362436069
  )
  
  processed_data <- as.data.frame(imputed_data$data)
  processed_data <- cbind(geneNames = gene_names, processed_data)
  
  return(processed_data)
}

prepare_datasets <- function(processed_data, dataset_config) {
  datasets <- list()
  
  for (config_name in names(dataset_config)) {
    config <- dataset_config[[config_name]]
    sample_indices <- config$samples
    
    if (is.null(sample_indices[2])) {
      sample_indices <- sample_indices[1]:ncol(processed_data)
    }
    
    dataset <- t(processed_data[, sample_indices])
    colnames(dataset) <- processed_data$geneNames
    datasets[[config_name]] <- dataset
  }
  
  return(datasets)
}

perform_quality_control <- function(datasets) {
  filtered_datasets <- list()
  
  for (dataset_name in names(datasets)) {
    expression_matrix <- datasets[[dataset_name]]
    
    gsg <- goodSamplesGenes(expression_matrix, verbose = 3)
    
    if (!gsg$allOK) {
      if (sum(!gsg$goodGenes) > 0) {
        message(paste("Removed genes:", sum(!gsg$goodGenes)))
      }
      if (sum(!gsg$goodSamples) > 0) {
        message(paste("Removed samples:", sum(!gsg$goodSamples)))
      }
      
      expression_matrix <- expression_matrix[gsg$goodSamples, gsg$goodGenes]
    }
    
    sample_tree <- hclust(dist(expression_matrix), method = "average")
    clusters <- cutreeStatic(sample_tree, 
                             cutHeight = ANALYSIS_PARAMS$cluster_cut_height,
                             minSize = ANALYSIS_PARAMS$min_cluster_size)
    
    keep_samples <- (clusters == 1)
    filtered_matrix <- expression_matrix[keep_samples, ]
    
    message(paste("Final dimensions for", dataset_name, ":", 
                  nrow(filtered_matrix), "samples ¡Á", ncol(filtered_matrix), "genes"))
    
    filtered_datasets[[dataset_name]] <- filtered_matrix
  }
  
  return(filtered_datasets)
}

find_common_genes <- function(datasets) {
  common_genes <- Reduce(intersect, lapply(datasets, colnames))
  return(common_genes)
}

select_variable_genes <- function(datasets, n_genes) {
  combined_expression <- do.call(rbind, datasets)
  mad_values <- apply(combined_expression, 2, mad)
  highly_variable <- combined_expression[, order(mad_values, decreasing = TRUE)[1:n_genes]]
  
  return(highly_variable)
}

simulate_trait_data <- function(expression_data) {
  trait_data <- list()
  
  for (dataset_name in names(expression_data)) {
    n_samples <- nrow(expression_data[[dataset_name]])
    
    traits <- data.frame(
      label = sample(c("Group_A", "Group_B", "Group_C"), n_samples, replace = TRUE),
      age = rnorm(n_samples, 60, 10),
      gender = sample(c("Male", "Female"), n_samples, replace = TRUE),
      stage = sample(paste0("Stage_", 1:4), n_samples, replace = TRUE)
    )
    
    trait_data[[dataset_name]] <- traits
  }
  
  return(trait_data)
}