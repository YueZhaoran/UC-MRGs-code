rm(list = ls()); gc()

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
project_dir <- dirname(getwd())
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(jsonlite)
})

section_keep <- NULL
out_mtx     <- "01_palantir_counts.mtx"
out_genes   <- "01_palantir_genes.txt"
out_barcodes<- "01_palantir_barcodes.txt"
out_meta    <- "01_palantir_metadata.csv"

seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))
if (!is.null(section_keep)) {
  seu_list <- seu_list[intersect(section_keep, names(seu_list))]
}

counts_all <- list()
meta_all <- data.frame()

for (i in 1:length(seu_list)) {
  sec <- names(seu_list)[i]
  seu <- seu_list[[i]]
  cnt <- GetAssayData(seu, assay = "Spatial", layer = "counts")
  colnames(cnt) <- paste0(sec, "_", colnames(cnt))

  coords <- GetTissueCoordinates(seu)
  counts_all[[sec]] <- cnt
  meta_all <- rbind(meta_all, data.frame(
    cell = colnames(cnt), section = sec, group = seu$group[1],
    x = coords$x, y = coords$y, cluster = as.character(seu$seurat_clusters),
    nCount = seu$nCount_Spatial, nFeature = seu$nFeature_Spatial))
}

genes <- Reduce(intersect, lapply(counts_all, rownames))
counts_all <- lapply(counts_all, function(x) x[genes, , drop = FALSE])
mat <- do.call(cbind, counts_all)

prop <- readRDS(file.path(project_dir, "05deconv", "06_spotlight_prop.rds"))
cell_cols <- c("B/Plasma", "T cell", "Epithelial", "Myeloid", "Stromal", "NK/ILC", "Unassigned")
prop$cell <- paste0(prop$sample, "_", prop$barcode)
prop <- prop[match(meta_all$cell, prop$cell), ]
meta_all$celltype <- cell_cols[max.col(as.matrix(prop[, cell_cols]), ties.method = "first")]

writeMM(mat, out_mtx)
writeLines(rownames(mat), out_genes)
writeLines(colnames(mat), out_barcodes)
write.csv(meta_all, out_meta, row.names = FALSE)

kg <- intersect(c("ADAMTS1", "FBN1", "SPARC", "VCAM1"), rownames(mat))
