plot_sample_dendrogram <- function(expression_data, trait_data, dataset_name, 
                                   output_dir, cut_height = 140) {
  
  sample_tree <- hclust(dist(expression_data), method = "average")
  
  if (!is.null(trait_data$label)) {
    trait_factor <- as.numeric(factor(trait_data$label))
    sample_colors <- numbers2colors(
      trait_factor,
      colors = rainbow(length(unique(trait_data$label))),
      signed = FALSE
    )
    group_label <- "Phenotype"
  } else {
    sample_colors <- NULL
    group_label <- NULL
  }
  
  output_file <- file.path(output_dir, 
                           paste0("sample_dendrogram_", tolower(dataset_name), ".pdf"))
  
  pdf(output_file, width = 8, height = 6)
  par(mar = c(1, 4, 3, 1), cex = 0.8)
  
  plotDendroAndColors(
    sample_tree, 
    sample_colors,
    groupLabels = group_label,
    cex.dendroLabels = 0.8,
    marAll = c(1, 4, 3, 1),
    cex.rowText = 0.01,
    main = paste("Sample dendrogram -", dataset_name)
  )
  
  if (!is.null(cut_height)) {
    abline(h = cut_height, col = "red")
  }
  
  dev.off()
}

plot_soft_threshold <- function(sft, output_dir, r_sq_cutoff = 0.85) {
  output_file <- file.path(output_dir, "soft_threshold_analysis.pdf")
  
  pdf(output_file, width = 10, height = 5)
  par(mfrow = c(1, 2))
  
  plot(sft$fitIndices[, 1], 
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n",
       main = "Scale Free Topology Model Fit")
  
  text(sft$fitIndices[, 1], 
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = sft$fitIndices[, 1], 
       cex = 0.7, 
       col = "red")
  
  abline(h = r_sq_cutoff, col = "red", lty = 2)
  
  plot(sft$fitIndices[, 1], 
       sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", 
       ylab = "Mean Connectivity", 
       type = "n",
       main = "Mean Connectivity")
  
  text(sft$fitIndices[, 1], 
       sft$fitIndices[, 5], 
       labels = sft$fitIndices[, 1], 
       cex = 0.7, 
       col = "red")
  
  abline(h = 100, col = "blue", lty = 2)
  
  dev.off()
}

plot_module_dendrogram <- function(net, output_dir) {
  output_file <- file.path(output_dir, "module_cluster_dendrogram.pdf")
  
  module_colors <- net$colors
  
  pdf(output_file, width = 10, height = 6)
  plotDendroAndColors(
    net$dendrograms[[1]], 
    module_colors[net$blockGenes[[1]]],
    "Module Colors",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05,
    main = "Gene Clustering and Module Colors"
  )
  dev.off()
}

plot_module_trait_heatmap <- function(expression_data, trait_data, module_colors, 
                                      dataset_name, output_dir, 
                                      excluded_modules = c("gold", "grey", "magenta", "pink")) {
  
  me_list <- moduleEigengenes(expression_data, module_colors)
  me_data <- me_list$eigengenes
  
  excluded_pattern <- paste0("ME", excluded_modules)
  me_data_filtered <- me_data[, !colnames(me_data) %in% excluded_pattern]
  
  trait_numeric <- trait_data
  for (col_name in colnames(trait_numeric)) {
    if (!is.numeric(trait_numeric[[col_name]])) {
      trait_numeric[[col_name]] <- as.numeric(factor(trait_numeric[[col_name]]))
    }
  }
  
  module_trait_cor <- cor(me_data_filtered, trait_numeric, 
                          use = "p", method = "spearman")
  module_trait_pvalue <- corPvalueStudent(module_trait_cor, 
                                          nrow(me_data_filtered))
  
  significance_matrix <- matrix("", nrow = nrow(module_trait_pvalue), 
                                ncol = ncol(module_trait_pvalue))
  
  significance_matrix[module_trait_pvalue < 0.001] <- "***"
  significance_matrix[module_trait_pvalue < 0.01 & module_trait_pvalue >= 0.001] <- "**"
  significance_matrix[module_trait_pvalue < 0.05 & module_trait_pvalue >= 0.01] <- "*"
  
  output_file <- file.path(output_dir, 
                           paste0("module_trait_heatmap_", tolower(dataset_name), ".pdf"))
  
  pdf(output_file, 
      width = max(6, 1 * ncol(trait_data)), 
      height = max(5, 0.6 * nrow(module_trait_cor)))
  
  par(mar = c(5, 9, 3, 3))
  
  labeledHeatmap(
    Matrix = module_trait_cor,
    xLabels = colnames(trait_data),
    yLabels = colnames(me_data_filtered),
    ySymbols = colnames(me_data_filtered),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = significance_matrix,
    setStdMargins = FALSE,
    cex.text = 1.0,
    zlim = c(-0.7, 0.7),
    main = paste("Module-Trait Relationships -", dataset_name)
  )
  
  dev.off()
  
  return(list(
    correlation = module_trait_cor,
    pvalue = module_trait_pvalue,
    module_eigengenes = me_data_filtered
  ))
}