ANALYSIS_PARAMS <- list(
  top_genes = 5000,
  cluster_cut_height = 150,
  min_cluster_size = 10,
  impute_k = 10,
  r_sq_cutoff = 0.85,
  min_module_size = 30,
  merge_cut_height = 0.35,
  network_type = "unsigned",
  mm_threshold = 0.8,
  gs_threshold = 0.2,
  n_permutations = 200,
  z_threshold = 10
)

DATASET_CONFIG <- list(
  batch = list(name = "Batch", samples = 2:36),
  tcga = list(name = "TCGA", samples = 37:443),
  cptac = list(name = "CPTAC", samples = 444:ncol(expr_matrix))
)

ENRICHMENT_PARAMS <- list(
  target_modules = c("green", "yellow"),
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  p_adjust_method = "BH",
  top_pathways = 150,
  databases = c("GO", "KEGG", "PID", "REACTOME", "WP", "BIOCARTA", "HALLMARK")
)

GSVA_PARAMS <- list(
  kcdf = "Gaussian",
  min_sz = 1,
  max_sz = Inf,
  parallel_cores = parallel::detectCores(),
  significance_threshold = 0.05
)