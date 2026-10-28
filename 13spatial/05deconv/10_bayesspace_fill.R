library(Seurat)
library(BayesSpace)
library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
key_genes <- read.csv(file.path(project_dir, "..", "05DEGs", "hub_gene.csv"))$symbol

done <- sub("^09_bayesspace_(.*)\\.rds$", "\\1",
            list.files(".", pattern = "^09_bayesspace_.*\\.rds$"))
todo <- setdiff(sheet$section, done)

for (sec in todo) {
  i <- match(sec, names(seu_list))
  seu <- seu_list[[i]]
  counts <- GetAssayData(seu, assay = "Spatial", layer = "counts")

  pos <- read.csv(file.path(project_dir, "00raw", paste0(sheet$gsm[i], "_", sec),
                            "spatial", "tissue_positions_list.csv"), header = FALSE,
                  col.names = c("barcode", "in_tissue", "array_row", "array_col",
                                "pxl_row", "pxl_col"))
  pos <- pos[match(colnames(counts), pos$barcode), ]

  sce <- SingleCellExperiment(assays = list(counts = counts))
  colData(sce)$array_row <- pos$array_row
  colData(sce)$array_col <- pos$array_col
  colData(sce)$pxl_row_in_fullres <- pos$pxl_row
  colData(sce)$pxl_col_in_fullres <- pos$pxl_col

  sce <- logNormCounts(sce)
  sce <- scater::runPCA(sce, subset_row = getTopHVGs(modelGeneVar(sce), n = 2000), ncomponents = 15)

  sce <- qTune(sce, qs = seq(4, 12, by = 2), burn.in = 100, nrep = 1000)
  sce <- spatialCluster(sce, q = 8, use.dimred = "PCA", d = 15, nrep = 10000, burn.in = 500)

  sce_enh <- spatialEnhance(sce, q = 8, use.dimred = "PCA", d = 15, nrep = 10000, burn.in = 500)
  enh_mat <- enhanceFeatures(sce_enh, sce, feature_names = key_genes, model = "lm",
                             feature.matrix = as.matrix(assay(sce, "logcounts")))

  enh <- data.frame(
    section = sec, group = seu$group[1],
    x = colData(sce_enh)$pxl_col_in_fullres, y = colData(sce_enh)$pxl_row_in_fullres,
    cluster = sce_enh$spatial.cluster,
    t(enh_mat[key_genes, , drop = FALSE]))
  write.csv(enh, paste0("09_enhanced_subspot_", sec, ".csv"), row.names = FALSE)

  saveRDS(list(sce = sce, sce_enhanced = sce_enh), paste0("09_bayesspace_", sec, ".rds"))

  p_q <- qPlot(sce) + labs(title = sec)
  ggsave(paste0("09_bayesspace_", sec, "_qtune.pdf"), p_q, width = 8, height = 6)
  ggsave(paste0("09_bayesspace_", sec, "_qtune.png"), p_q, width = 8, height = 6, dpi = 300)

}
