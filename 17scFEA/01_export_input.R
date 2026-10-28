rm(list = ls()); gc()

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
project_dir <- file.path(getwd(), "..")

suppressMessages({
  library(Seurat)
  library(Matrix)
})

target_celltype <- "Stromal"
out_csv         <- "01_input_for_scfea.csv"
n_per_type      <- 400

scData <- readRDS(file.path(project_dir, "01singleCell", "05_scData_annotated.rds"))

DefaultAssay(scData) <- "RNA"
if ("JoinLayers" %in% ls("package:Seurat")) {
  scData <- JoinLayers(scData, assay = "RNA")
}

set.seed(123)
cell_keep <- unlist(lapply(split(colnames(scData), scData$external_annotation), function(x) {
  sample(x, min(length(x), n_per_type))
}))
scData <- subset(scData, cells = cell_keep)

scData <- NormalizeData(scData, assay = "RNA", normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)

expr <- as.matrix(GetAssayData(scData, assay = "RNA", layer = "data"))
expr <- round(expr, 4)
expr_out <- data.frame(gene = rownames(expr), expr, check.names = FALSE)
write.csv(expr_out, out_csv, row.names = FALSE)

meta <- data.frame(cell = colnames(scData),
                   celltype = scData$external_annotation,
                   group = scData$group,
                   keycell = ifelse(scData$external_annotation == target_celltype, "Yes", "No"))
write.csv(meta, "01_cell_metadata.csv", row.names = FALSE)
