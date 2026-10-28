rm(list = ls()); gc()

suppressMessages({
  library(Seurat)
  library(Matrix)
})

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

target_celltype   <- "Stromal"
compare_celltypes <- c("B/Plasma", "Epithelial", "Myeloid", "NK/ILC", "T cell")
group_keep        <- "UC"
n_per_type        <- 300
min_cells_expr    <- 10

sc <- readRDS(file.path(project_dir, "..", "01singleCell", "05_scData_annotated.rds"))
sc <- subset(sc, subset = group == group_keep &
               external_annotation %in% c(target_celltype, compare_celltypes))
sc$celltype_group <- as.character(sc$external_annotation)

keep_cells <- unlist(lapply(split(colnames(sc), sc$celltype_group), function(x) {
  sample(x, min(length(x), n_per_type))
}))
sc <- subset(sc, cells = keep_cells)

DefaultAssay(sc) <- "RNA"
counts <- GetAssayData(sc, assay = "RNA", layer = "counts")
counts <- counts[Matrix::rowSums(counts > 0) >= min_cells_expr, ]

write.csv(as.matrix(counts), "09_scenic_counts.csv")
meta_out <- data.frame(cell = colnames(sc), celltype_group = sc$celltype_group,
                       sample = sc$sample_name, group = sc$group)
write.csv(meta_out, "09_scenic_cells.csv", row.names = FALSE)
