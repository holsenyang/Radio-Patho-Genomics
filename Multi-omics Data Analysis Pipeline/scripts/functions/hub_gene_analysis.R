calculate_module_membership <- function(expression_data, module_eigengenes, 
                                        trait_data, module_colors, 
                                        target_module, dataset_name) {
  
  trait_numeric <- as.numeric(trait_data$label)
  trait_df <- data.frame(label = trait_numeric)
  colnames(trait_df) <- "label"
  
  mod_names <- substring(colnames(module_eigengenes), 3)
  
  gene_module_membership <- as.data.frame(
    cor(expression_data, module_eigengenes, use = "p")
  )
  
  gene_trait_significance <- as.data.frame(
    cor(expression_data, trait_df, use = "p", method = "spearman")
  )
  
  module_genes <- module_colors == target_module
  module_gene_names <- colnames(expression_data)[module_genes]
  
  module_column <- match(target_module, mod_names)
  mm_values <- abs(gene_module_membership[module_genes, module_column])
  gs_values <- abs(gene_trait_significance[module_genes, 1])
  
  result_df <- data.frame(
    gene_id = module_gene_names,
    MM = mm_values,
    GS = gs_values,
    stringsAsFactors = FALSE
  )
  
  if (exists("gene_id_name")) {
    result_df <- merge(result_df, gene_id_name, by = "gene_id", all.x = TRUE)
  }
  
  rownames(result_df) <- result_df$gene_id
  
  return(list(
    mm_gs_data = result_df,
    module_genes = module_gene_names,
    dataset = dataset_name
  ))
}

identify_hub_genes <- function(mm_gs_data, mm_threshold = 0.8, gs_threshold = 0.2) {
  mm_gs_data$is_hub <- (abs(mm_gs_data$MM) > mm_threshold) & 
    (abs(mm_gs_data$GS) > gs_threshold)
  
  hub_genes <- mm_gs_data[mm_gs_data$is_hub, ]
  
  return(hub_genes)
}

plot_mm_gs_scatter <- function(mm_gs_data, hub_genes, target_module, 
                               dataset_name, output_dir, label_genes = NULL) {
  
  plot_data <- mm_gs_data
  plot_data$group <- ifelse(plot_data$is_hub, "Hub", "Non-hub")
  
  if (!is.null(label_genes)) {
    label_data <- plot_data[plot_data$gene_id %in% label_genes, ]
  } else {
    label_data <- hub_genes
  }
  
  plot <- ggplot(data = plot_data, aes(x = MM, y = GS, color = group)) + 
    geom_point(size = 2) +
    scale_colour_manual(values = c("Hub" = "#DE6757", "Non-hub" = "grey60")) +
    theme_bw() +  
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = 'none'
    ) +
    labs(
      x = paste("Module Membership in", target_module, "module"),
      y = paste("Gene significance for", dataset_name),
      title = "Module membership vs. gene significance"
    ) +
    geom_hline(yintercept = 0.2, colour = "#5B9BD5", linewidth = 0.8, linetype = 3) + 
    geom_vline(xintercept = 0.8, colour = "#5B9BD5", linewidth = 0.8, linetype = 3)
  
  if (nrow(label_data) > 0) {
    plot <- plot +
      geom_point(data = label_data, aes(x = MM, y = GS), 
                 color = target_module, size = 3, shape = 16) +
      geom_text_repel(
        data = label_data,
        aes(label = ifelse(exists("gene_name"), gene_name, gene_id)),
        color = "black", 
        size = 4, 
        fontface = "italic",
        box.padding = 0.2,
        point.padding = 0.3,
        segment.color = 'black', 
        segment.size = 0.3,
        max.overlaps = 20
      )
  }
  
  output_file <- file.path(
    output_dir, 
    paste0("mm_gs_", target_module, "_", tolower(dataset_name), ".pdf")
  )
  
  ggsave(
    filename = output_file,
    plot = plot,
    width = 8,
    height = 6
  )
  
  return(plot)
}

find_common_hub_genes <- function(hub_gene_list, dataset_names) {
  hub_gene_ids <- lapply(hub_gene_list, function(x) x$gene_id)
  names(hub_gene_ids) <- dataset_names
  
  common_genes <- Reduce(intersect, hub_gene_ids)
  
  if (length(common_genes) > 0) {
    common_df <- data.frame(
      gene_id = common_genes,
      stringsAsFactors = FALSE
    )
    
    if (exists("gene_id_name")) {
      common_df <- merge(common_df, gene_id_name, by = "gene_id", all.x = TRUE)
    }
  } else {
    common_df <- data.frame(
      gene_id = character(),
      gene_name = character(),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    common_genes = common_df,
    all_hub_genes = hub_gene_ids
  ))
}

create_hub_gene_summary <- function(mm_gs_results, common_hub_genes, target_module) {
  common_genes <- common_hub_genes$common_genes$gene_id
  
  if (length(common_genes) == 0) {
    return(data.frame())
  }
  
  summary_data <- list()
  
  for (i in seq_along(mm_gs_results)) {
    dataset_name <- names(mm_gs_results)[i]
    dataset_data <- mm_gs_results[[i]]$mm_gs_data
    
    common_data <- dataset_data[dataset_data$gene_id %in% common_genes, 
                                c("MM", "GS")]
    colnames(common_data) <- paste0(dataset_name, "_", c("MM", "GS"))
    
    summary_data[[dataset_name]] <- common_data
  }
  
  if (length(summary_data) > 0) {
    summary_df <- do.call(cbind, summary_data)
    summary_df$gene_id <- common_genes
    
    if (exists("gene_id_name")) {
      summary_df <- merge(summary_df, gene_id_name, by = "gene_id", all.x = TRUE)
    }
    
    gene_info_cols <- c("gene_id", "gene_name")
    data_cols <- setdiff(colnames(summary_df), gene_info_cols)
    summary_df <- summary_df[, c(gene_info_cols, data_cols)]
    
    rownames(summary_df) <- summary_df$gene_id
    
    return(summary_df)
  } else {
    return(data.frame())
  }
}

analyze_module_across_datasets <- function(target_module, network_results, 
                                           filtered_datasets) {
  
  mm_gs_results <- list()
  hub_genes_list <- list()
  
  for (dataset_key in names(filtered_datasets)) {
    expression_data <- filtered_datasets[[dataset_key]]
    trait_data <- simulate_trait_data(filtered_datasets)[[dataset_key]]
    module_colors <- network_results[[dataset_key]]$network$colors
    
    mm_gs_result <- calculate_module_membership(
      expression_data = expression_data,
      module_eigengenes = network_results[[dataset_key]]$module_trait$module_eigengenes,
      trait_data = trait_data,
      module_colors = module_colors,
      target_module = target_module,
      dataset_name = dataset_key
    )
    
    mm_gs_results[[dataset_key]] <- mm_gs_result
    
    hub_genes <- identify_hub_genes(
      mm_gs_result$mm_gs_data,
      mm_threshold = ANALYSIS_PARAMS$mm_threshold,
      gs_threshold = ANALYSIS_PARAMS$gs_threshold
    )
    
    hub_genes_list[[dataset_key]] <- hub_genes
  }
  
  common_hub_genes <- find_common_hub_genes(
    hub_genes_list,
    dataset_names = names(hub_genes_list)
  )
  
  hub_gene_summary <- create_hub_gene_summary(
    mm_gs_results,
    common_hub_genes,
    target_module
  )
  
  plot_results <- list()
  for (dataset_name in names(mm_gs_results)) {
    plot_result <- plot_mm_gs_scatter(
      mm_gs_data = mm_gs_results[[dataset_name]]$mm_gs_data,
      hub_genes = hub_genes_list[[dataset_name]],
      target_module = target_module,
      dataset_name = dataset_name,
      output_dir = output_dirs$figures,
      label_genes = common_hub_genes$common_genes$gene_id
    )
    
    plot_results[[dataset_name]] <- plot_result
  }
  
  return(list(
    mm_gs_results = mm_gs_results,
    hub_genes = hub_genes_list,
    common_hub_genes = common_hub_genes,
    hub_gene_summary = hub_gene_summary,
    plots = plot_results
  ))
}

generate_hub_gene_report <- function(all_results, hub_gene_summary) {
  report_file <- file.path(output_dirs$rdata, "hub_gene_analysis_report.txt")
  
  sink(report_file)
  
  cat("Hub Gene Analysis Report\n")
  cat("========================\n\n")
  cat("Generated on:", date(), "\n\n")
  
  cat("Summary Statistics:\n")
  cat("-------------------\n")
  
  for (module_name in names(all_results)) {
    module_result <- all_results[[module_name]]
    
    cat(paste("\nModule:", module_name, "\n"))
    
    for (dataset_name in names(module_result$hub_genes)) {
      n_hub_genes <- nrow(module_result$hub_genes[[dataset_name]])
      n_total_genes <- nrow(module_result$mm_gs_results[[dataset_name]]$mm_gs_data)
      
      cat(paste("  ", dataset_name, ":", n_hub_genes, "/", n_total_genes, 
                "genes (", round(n_hub_genes/n_total_genes*100, 1), "%)\n"))
    }
    
    n_common <- nrow(module_result$common_hub_genes$common_genes)
    cat(paste("  Common hub genes:", n_common, "\n"))
  }
  
  if (nrow(hub_gene_summary) > 0) {
    cat("\nCommon Hub Genes Across All Modules:\n")
    cat("------------------------------------\n")
    print(hub_gene_summary)
  }
  
  sink()
}