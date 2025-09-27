perform_module_preservation <- function(multi_expr, multi_color, 
                                        n_permutations = 200, 
                                        random_seed = 1,
                                        reference_network = 1) {
  
  system.time({
    preservation_result <- modulePreservation(
      multiExpr = multi_expr,
      multiColor = multi_color,
      referenceNetworks = reference_network,
      nPermutations = n_permutations,
      randomSeed = random_seed,
      quickCor = 0,
      verbose = 3
    )
  })
  
  return(preservation_result)
}

extract_preservation_stats <- function(preservation_result, ref_index = 1, test_index = 2) {
  stats_observed <- cbind(
    preservation_result$quality$observed[[ref_index]][[test_index]][, -1],
    preservation_result$preservation$observed[[ref_index]][[test_index]][, -1]
  )
  
  stats_z <- cbind(
    preservation_result$quality$Z[[ref_index]][[test_index]][, -1],
    preservation_result$preservation$Z[[ref_index]][[test_index]][, -1]
  )
  
  result_df <- data.frame(
    mergeModule = rownames(stats_observed),
    medianRank.pres = stats_observed[, "medianRank.pres"],
    medianRank.qual = stats_observed[, "medianRank.qual"],
    Zsummary.pres = signif(stats_z[, "Zsummary.pres"], 3),
    Zsummary.qual = signif(stats_z[, "Zsummary.qual"], 3),
    stringsAsFactors = FALSE
  )
  
  return(list(
    observed = stats_observed,
    z = stats_z,
    summary = result_df
  ))
}

filter_preserved_modules <- function(preservation_stats, z_threshold = 10, 
                                     excluded_modules = c("grey", "gold")) {
  
  preserved_modules <- preservation_stats$summary[
    preservation_stats$summary$Zsummary.pres > z_threshold & 
      !preservation_stats$summary$mergeModule %in% excluded_modules,
  ]
  
  return(preserved_modules)
}

plot_preservation_summary <- function(preservation_result, ref_index = 1, test_index = 2,
                                      dataset_name = "Test", output_dir) {
  
  mod_colors <- rownames(preservation_result$preservation$observed[[ref_index]][[test_index]])
  module_sizes <- preservation_result$preservation$Z[[ref_index]][[test_index]][, 1]
  
  plot_mods <- !(mod_colors %in% c("grey", "gold"))
  text_labels <- mod_colors[plot_mods]
  
  plot_data <- cbind(
    preservation_result$preservation$observed[[ref_index]][[test_index]][, 2],
    preservation_result$preservation$Z[[ref_index]][[test_index]][, 2]
  )
  
  output_file <- file.path(output_dir, 
                           paste0("module_preservation_summary_", tolower(dataset_name), ".pdf"))
  
  pdf(output_file, width = 10, height = 5)
  par(mfrow = c(1, 2))
  par(mar = c(4.5, 4.5, 2.5, 1))
  
  mains <- c("Preservation Median rank", 
             paste("Preservation Zsummary (", dataset_name, ")", sep = ""))
  
  for (p in 1:2) {
    y_data <- plot_data[plot_mods, p]
    min_val <- min(y_data, na.rm = TRUE)
    max_val <- max(y_data, na.rm = TRUE)
    
    if (p %in% c(1, 2)) {
      if (min_val > -max_val/10) min_val <- -max_val/10
      y_lim <- c(min_val - 0.1 * (max_val - min_val), 
                 max_val + 0.1 * (max_val - min_val))
    } else {
      y_lim <- c(max_val + 0.1 * (max_val - min_val), 
                 min_val - 0.1 * (max_val - min_val))
    }
    
    plot(module_sizes[plot_mods], y_data,
         col = 1, bg = mod_colors[plot_mods], pch = 21,
         main = mains[p],
         cex = 2.4,
         ylab = mains[p], xlab = "Module size", log = "x",
         ylim = y_lim,
         xlim = c(10, 2000), 
         cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2)
    
    labelPoints(module_sizes[plot_mods], y_data, text_labels, 
                cex = 0.8, offs = 0.2)
    
    if (p == 2) {
      abline(h = 0)
      abline(h = 2, col = "blue", lty = 2)
      abline(h = 10, col = "darkgreen", lty = 2)
    }
  }
  
  dev.off()
}

plot_detailed_zstats <- function(preservation_stats, preservation_result, 
                                 ref_index = 1, test_index = 2,
                                 dataset_name = "Test", output_dir) {
  
  stats_z <- preservation_stats$z
  mod_colors <- rownames(stats_z)
  module_sizes <- preservation_result$quality$Z[[ref_index]][[test_index]][, 1]
  
  plot_mods <- !(mod_colors %in% c("grey", "gold"))
  
  labs <- match(mod_colors[plot_mods], standardColors(50))
  
  output_file <- file.path(output_dir, 
                           paste0("module_preservation_zstats_", tolower(dataset_name), ".pdf"))
  
  pdf(output_file, width = 12, height = 9)
  par(mfrow = c(4, 5))
  par(mar = c(3, 3, 2, 1))
  par(mgp = c(1.6, 0.4, 0))
  
  for (s in 1:ncol(stats_z)) {
    y_data <- stats_z[plot_mods, s]
    min_val <- min(y_data, na.rm = TRUE)
    max_val <- max(y_data, na.rm = TRUE)
    
    if (min_val > -max_val/5) min_val <- -max_val/5
    
    plot(module_sizes[plot_mods], y_data,
         col = 1, bg = mod_colors[plot_mods], pch = 21,
         main = colnames(stats_z)[s],
         cex = 1.7,
         ylab = colnames(stats_z)[s], xlab = "Module size", log = "x",
         ylim = c(min_val - 0.1 * (max_val - min_val), 
                  max_val + 0.1 * (max_val - min_val)),
         xlim = c(20, 1000))
    
    labelPoints(module_sizes[plot_mods], y_data, labs, 
                cex = 0.7, offs = 0.04)
    
    abline(h = 0)
    abline(h = 2, col = "blue", lty = 2)
    abline(h = 10, col = "darkgreen", lty = 2)
  }
  
  dev.off()
}

save_preservation_results <- function(preservation_stats, module_counts, 
                                      dataset_name = "Test", output_dir) {
  
  result_df <- merge(preservation_stats$summary, module_counts, 
                     by = "mergeModule", all.x = TRUE)
  
  result_df <- result_df[!result_df$mergeModule %in% c('grey', 'gold'), ]
  
  output_file <- file.path(output_dir, 
                           paste0("module_preservation_results_", tolower(dataset_name), ".csv"))
  
  write.csv(result_df, file = output_file, row.names = FALSE)
  
  return(result_df)
}

analyze_dataset_preservation <- function(dataset_config, module_colors) {
  reference_expr <- get("exp_dat", envir = .GlobalEnv)
  test_expr <- get(dataset_config$expression, envir = .GlobalEnv)
  
  multi_expr <- list(
    Reference = list(data = reference_expr),
    Test = list(data = test_expr)
  )
  
  multi_color <- list(Reference = module_colors)
  
  preservation_result <- perform_module_preservation(
    multi_expr = multi_expr,
    multi_color = multi_color,
    n_permutations = ANALYSIS_PARAMS$n_permutations,
    random_seed = ANALYSIS_PARAMS$random_seed,
    reference_network = 1
  )
  
  preservation_stats <- extract_preservation_stats(
    preservation_result = preservation_result,
    ref_index = 1,
    test_index = 2
  )
  
  preserved_modules <- filter_preserved_modules(
    preservation_stats = preservation_stats,
    z_threshold = ANALYSIS_PARAMS$z_threshold,
    excluded_modules = c("grey", "gold")
  )
  
  module_counts <- as.data.frame(table(module_colors))
  colnames(module_counts) <- c("mergeModule", "Number_of_genes")
  
  result_table <- save_preservation_results(
    preservation_stats = preservation_stats,
    module_counts = module_counts,
    dataset_name = dataset_config$name,
    output_dir = output_dirs$tables
  )
  
  plot_preservation_summary(
    preservation_result = preservation_result,
    ref_index = 1,
    test_index = 2,
    dataset_name = dataset_config$name,
    output_dir = output_dirs$figures
  )
  
  plot_detailed_zstats(
    preservation_stats = preservation_stats,
    preservation_result = preservation_result,
    ref_index = 1,
    test_index = 2,
    dataset_name = dataset_config$name,
    output_dir = output_dirs$figures
  )
  
  return(list(
    preservation_result = preservation_result,
    preservation_stats = preservation_stats,
    preserved_modules = preserved_modules,
    result_table = result_table
  ))
}

generate_preservation_report <- function(preservation_results) {
  report_file <- file.path(output_dirs$rdata, "module_preservation_report.txt")
  
  sink(report_file)
  
  cat("Module Preservation Analysis Report\n")
  cat("===================================\n\n")
  cat("Generated on:", date(), "\n\n")
  
  cat("Preservation Summary:\n")
  cat("---------------------\n\n")
  
  for (dataset_key in names(preservation_results)) {
    result <- preservation_results[[dataset_key]]
    
    cat(paste("Dataset:", toupper(dataset_key), "\n"))
    cat(paste("Number of preserved modules (Z > 10):", 
              nrow(result$preserved_modules), "\n"))
    
    if (nrow(result$preserved_modules) > 0) {
      cat("Preserved modules:\n")
      print(result$preserved_modules$mergeModule)
    }
    
    cat("\n")
  }
  
  cat("Zsummary Statistics:\n")
  cat("--------------------\n\n")
  
  for (dataset_key in names(preservation_results)) {
    result <- preservation_results[[dataset_key]]
    
    cat(paste("Dataset:", toupper(dataset_key), "\n"))
    
    z_summary <- result$preservation_stats$summary
    z_summary <- z_summary[!z_summary$mergeModule %in% c('grey', 'gold'), ]
    
    z_summary <- z_summary[order(-z_summary$Zsummary.pres), ]
    
    print(z_summary[, c("mergeModule", "Zsummary.pres", "Number_of_genes")])
    cat("\n")
  }
  
  sink()
}