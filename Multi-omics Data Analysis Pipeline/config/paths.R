# rm(list = ls())
# rstudioapi::getActiveDocumentContext()$path
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# project_root=getwd()

project_root <- here::here()
data_dirs <- list(
  raw = file.path(project_root, "data/raw"),
  processed = file.path(project_root, "data/processed")
)

output_dirs <- list(
  figures = file.path(project_root, "output/figures"),
  tables = file.path(project_root, "output/tables"),
  rdata = file.path(project_root, "output/rdata")
)

lapply(c(data_dirs, output_dirs), function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
})
