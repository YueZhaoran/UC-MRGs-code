# ============================ R environment inventory ============================
# Purpose: record the R version and package versions used in this project
# Usage: Rscript 01_env_versions.R
# Output: R_environment.txt (R version and package versions)

rm(list = ls())
gc()

# Output directory (script directory; works with both Rscript and RStudio)
out_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
)
if (!dir.exists(out_dir)) out_dir <- getwd()
setwd(out_dir)
cat("Working directory: ", getwd(), "\n")

# R packages used across the analysis modules
pkg_list <- c(
  "AnnotationDbi", "AUCell", "Augur", "BayesSpace", "BiocParallel", "biomaRt",
  "caret", "CellChat", "celldex", "clusterProfiler", "clustree", "ComplexHeatmap",
  "data.table", "dplyr", "e1071", "enrichplot", "forcats", "GEOquery",
  "ggh4x", "ggplot2", "ggpubr", "ggrepel", "ggsci", "ggvenn",
  "Giotto", "glmnet", "GSEABase", "GseaVis", "GSVA", "gtable",
  "gwasvcf", "harmony", "httr", "igraph", "jsonlite", "KEGG.db",
  "limma", "Matrix", "mclust", "miloR", "monocle3", "msigdbr",
  "multiMiR", "neuralnet", "org.Hs.eg.db", "paletteer", "patchwork", "pdftools",
  "pheatmap", "png", "pROC", "psych", "randomForest", "RColorBrewer",
  "ReactomeGSA", "readxl", "reshape2", "reticulate", "rstatix", "rstudioapi",
  "rtracklayer", "scales", "scater", "scCustomize", "Scissor", "scPagwas",
  "scran", "scTenifoldKnk", "scTenifoldNet", "Seurat", "SingleCellExperiment", "SingleR",
  "singscore", "SPOTlight", "stringr", "tibble", "tidyr", "tidyverse",
  "UCell", "VariantAnnotation", "VennDiagram", "viridis", "xCell"
)

# Record R version and platform
r_version <- R.version.string
platform  <- R.version$platform
cat("R version: ", r_version, "\n")
cat("Platform: ", platform, "\n")

# Check and record each package version
ver_df <- data.frame(package = character(), version = character(),
                     installed = logical(), stringsAsFactors = FALSE)
for (pkg in pkg_list) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (ok) as.character(packageVersion(pkg)) else NA
  ver_df <- rbind(ver_df, data.frame(package = pkg, version = ver,
                                     installed = ok, stringsAsFactors = FALSE))
  cat(sprintf("%-20s %s\n", pkg, ifelse(ok, ver, "not installed")))
}

# Write the combined environment inventory
sink("R_environment.txt")
cat("===== R environment inventory =====\n")
cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version: ", r_version, "\n")
cat("Platform: ", platform, "\n\n")
cat("===== Package versions used in this project =====\n")
print(ver_df[, c("package", "version")], row.names = FALSE)
sink()

cat("\nDone. R_environment.txt has been generated.\n")
