library(Seurat)
library(SPOTlight)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

sc <- readRDS(file.path(project_dir, "..", "01singleCell", "05_scData_annotated.rds"))
sc <- subset(sc, cells = colnames(sc)[!is.na(sc$external_annotation)])
DefaultAssay(sc) <- "RNA"
sc <- NormalizeData(sc, verbose = FALSE)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sc <- ScaleData(sc, verbose = FALSE)

keep_cells <- character()
for (ct in sort(unique(sc$external_annotation))) {
  cells <- colnames(sc)[sc$external_annotation == ct]
  keep_cells <- c(keep_cells, sample(cells, min(length(cells), 200)))
}
sc <- subset(sc, cells = keep_cells)

mgs <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, verbose = FALSE)
mgs$gene <- rownames(mgs)

seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))
celltypes <- sort(unique(sc$external_annotation))
prop <- data.frame()

for (i in 1:length(seu_list)) {
  seu <- seu_list[[i]]
  DefaultAssay(seu) <- "Spatial"

  res <- SPOTlight(
    x = as.matrix(GetAssayData(sc, assay = "RNA", layer = "counts")),
    y = as.matrix(GetAssayData(seu, assay = "Spatial", layer = "counts")),
    groups = as.character(sc$external_annotation),
    mgs = mgs,
    n_top = 1000,
    weight_id = "avg_log2FC",
    group_id = "cluster",
    gene_id = "gene",
    verbose = FALSE)

  mat <- as.data.frame(res$mat)
  mat$barcode <- rownames(mat)
  mat$sample <- seu$section[1]
  mat$group <- seu$group[1]
  prop <- rbind(prop, mat)
}
prop$dominant <- celltypes[max.col(as.matrix(prop[, celltypes]), ties.method = "first")]
write.csv(prop, "06_spotlight_prop.csv", row.names = FALSE)
saveRDS(prop, "06_spotlight_prop.rds")

prop_sum <- data.frame()
for (s in unique(prop$sample)) {
  sub <- prop[prop$sample == s, ]
  for (ct in celltypes) {
    prop_sum <- rbind(prop_sum, data.frame(
      sample = s, group = sub$group[1], celltype = ct,
      mean_prop = round(mean(sub[[ct]]), 4), median_prop = round(median(sub[[ct]]), 4)))
  }
}
write.csv(prop_sum, "06_prop_summary.csv", row.names = FALSE)
