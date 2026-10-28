work_dir <- getwd()
setwd(work_dir)

suppressPackageStartupMessages({
  library(reticulate)
  library(Seurat)
  library(Matrix)
})

key_genes     <- c("ADAMTS1", "FBN1", "SPARC", "VCAM1")
cell_type_col <- "external_annotation"
cell_keep     <- c("Stromal", "Epithelial", "Myeloid")
drop_subtypes <- c("Endothelial_cells", "Neurons")
group_keep    <- "UC"
seurat_rds    <- "05_scData_annotated.rds"
n_hvg         <- 3000
output_h5ad   <- "03_celloracle_input.h5ad"
tf_file       <- "01_tf_candidates_all.csv"

conda_env   <- Sys.getenv("CELLORACLE_ENV", "celloracle")
python_path <- Sys.getenv("CELLORACLE_PYTHON", "")
if (nzchar(python_path)) {
  use_python(python_path, required = TRUE)
} else {
  use_condaenv(conda_env, required = TRUE)
}

seurat_obj <- readRDS(seurat_rds)
DefaultAssay(seurat_obj) <- "RNA"
if ("JoinLayers" %in% ls("package:Seurat")) seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

Idents(seurat_obj) <- cell_type_col
keep_cell <- seurat_obj[[cell_type_col]][, 1] != "Stromal" |
  !(as.character(seurat_obj$labels) %in% drop_subtypes)
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_cell])
seurat_sub <- subset(seurat_obj, idents = cell_keep, subset = group == group_keep)

seurat_sub <- NormalizeData(seurat_sub, verbose = FALSE)
seurat_sub <- FindVariableFeatures(seurat_sub, nfeatures = n_hvg, verbose = FALSE)
hvg_genes <- VariableFeatures(seurat_sub)

seurat_sub <- ScaleData(seurat_sub, verbose = FALSE)
seurat_sub <- RunPCA(seurat_sub, npcs = 30, verbose = FALSE)
seurat_sub <- RunUMAP(seurat_sub, dims = 1:30, seed.use = 123, verbose = FALSE)

tf_data <- read.csv(tf_file, stringsAsFactors = FALSE)
key_tf <- unique(toupper(trimws(tf_data$TF)))

extra_tf_file <- "scenic_candidate_TFs.txt"
extra_tf <- if (file.exists(extra_tf_file)) {
  unique(toupper(trimws(readLines(extra_tf_file))))
} else {
  character(0)
}

all_genes_of_interest <- unique(c(key_genes, key_tf, extra_tf))
gene_use <- unique(c(head(hvg_genes, n_hvg), all_genes_of_interest))
gene_use <- intersect(gene_use, rownames(seurat_sub))
seurat_sub <- seurat_sub[gene_use, ]

missing_genes <- setdiff(all_genes_of_interest, rownames(seurat_obj))
if (length(missing_genes) > 0) cat("表达矩阵中缺失：", paste(missing_genes, collapse = ", "), "\n")
if (!"umap" %in% names(seurat_sub@reductions)) stop("单细胞对象没有 umap 降维")

ad <- import("anndata", convert = FALSE)
np <- import("numpy", convert = FALSE)

counts_mat <- GetAssayData(seurat_sub, layer = "counts", assay = "RNA")
data_mat   <- GetAssayData(seurat_sub, layer = "data", assay = "RNA")

meta_data <- seurat_sub@meta.data
for (col in colnames(meta_data)) {
  if (is.factor(meta_data[[col]])) meta_data[[col]] <- as.character(meta_data[[col]])
}
meta_data$celltype <- as.character(meta_data[[cell_type_col]])
meta_data$celltype_group <- meta_data$celltype

to_np <- function(m) np$array(t(as.matrix(m)))
umap_coords <- as.matrix(Embeddings(seurat_sub, reduction = "umap"))

adata <- ad$AnnData(
  X = to_np(counts_mat),
  obs = meta_data,
  var = data.frame(gene = rownames(seurat_sub), row.names = rownames(seurat_sub)),
  layers = list(data = to_np(data_mat)),
  obsm = list(X_umap = np$array(umap_coords))
)

output_path <- file.path(work_dir, output_h5ad)
adata$write_h5ad(output_path)
