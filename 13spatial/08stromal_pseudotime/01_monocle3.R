rm(list = ls()); gc()

library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(igraph)
library(patchwork)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

target_celltype <- "Stromal"
rest_markers <- c("RGS5", "NOTCH3", "PDGFRB")
n_stages <- 3

sc <- readRDS(file.path(project_dir, "..", "01singleCell", "05_scData_annotated.rds"))
key_genes <- read.csv(file.path(project_dir, "..", "05DEGs", "hub_gene.csv"))$symbol

pbmc <- subset(sc, subset = external_annotation == target_celltype)
drop_subtypes <- c("Endothelial_cells", "Neurons")
pbmc <- subset(pbmc, subset = !labels %in% drop_subtypes)
n_sample <- length(unique(pbmc$sample_name))
rm(sc); gc()

DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, verbose = FALSE)

data <- GetAssayData(pbmc, assay = "RNA", layer = "counts")
cell_metadata <- pbmc@meta.data[, setdiff(colnames(pbmc@meta.data), "sample_name"), drop = FALSE]
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP), ]
colnames(int.embed) <- c("umap_1", "umap_2")
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds, resolution = 1e-3)

count_branch_points <- function(cds) {
  g <- principal_graph(cds)[["UMAP"]]
  c(branch_points = sum(igraph::degree(g) >= 3), leaves = sum(igraph::degree(g) == 1))
}

bp_ok <- FALSE
for (nc in c(200, 160, 150, 120, 75, 55, 20, 15)) {
  cds_try <- learn_graph(cds, use_partition = FALSE,
                         learn_graph_control = list(euclidean_distance_ratio = 0.5, ncenter = nc))
  st <- count_branch_points(cds_try)
  if (st["branch_points"] <= 1) {
    cds <- cds_try
    bp_ok <- TRUE
    break
  }
}
if (!bp_ok) {
  cds <- cds_try
}
st <- count_branch_points(cds)

act_markers <- intersect(c("COL1A1", "COL1A2", "MMP1", "MMP3", "IL11", "POSTN", "INHBA"),
                         rownames(data))
cnt_all <- counts(cds)
cpm <- log1p(sweep(as.matrix(cnt_all), 2, Matrix::colSums(cnt_all), "/") * 1e4)
rest_score_all <- colMeans(cpm[rest_markers, , drop = FALSE])
act_score_all <- colMeans(cpm[act_markers, , drop = FALSE])
axis_score <- scale(rest_score_all) - scale(act_score_all)
axis_score <- as.numeric(axis_score)

g <- principal_graph(cds)[["UMAP"]]
closest_vertex <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds), ])
node_axis <- tapply(axis_score, closest_vertex[, 1], mean)
node_n <- tapply(axis_score, closest_vertex[, 1], length)
root_idx <- as.numeric(names(node_axis)[which.max(node_axis)])
root_node <- igraph::V(g)$name[root_idx]
cds <- order_cells(cds, root_pr_nodes = root_node)

pt <- pseudotime(cds)
rho_rest <- cor(rest_score_all, pt, method = "spearman")
rho_act <- cor(act_score_all, pt, method = "spearman")
if (rho_rest < 0) {
  pt <- max(pt, na.rm = TRUE) - pt
  cds@principal_graph_aux[["UMAP"]]$pseudotime <- pt
  rho_rest <- cor(rest_score_all, pt, method = "spearman")
  rho_act <- cor(act_score_all, pt, method = "spearman")
}

saveRDS(cds, "01_cds_pseudotime.rds")
pt_tbl <- data.frame(cell = colnames(cds), pseudotime = pt,
                     subtype = colData(cds)$labels, group = colData(cds)$group)
write.csv(pt_tbl, "01_pseudotime.csv", row.names = FALSE)

umap_coords <- reducedDims(cds)[["UMAP"]]
brks <- unique(quantile(pt, probs = seq(0, 1, length.out = n_stages + 1), na.rm = TRUE))
stage_bins <- cut(pt, breaks = brks, include.lowest = TRUE, labels = seq_len(length(brks) - 1))
stage_pos <- data.frame(stage = seq_len(length(brks) - 1), umap_1 = NA, umap_2 = NA,
                        pt_mid = NA, n_cells = NA)
for (s in seq_len(nrow(stage_pos))) {
  cells_s <- names(pt)[which(stage_bins == as.character(s))]
  stage_pos$umap_1[s] <- median(umap_coords[cells_s, 1])
  stage_pos$umap_2[s] <- median(umap_coords[cells_s, 2])
  stage_pos$pt_mid[s] <- median(pt[cells_s], na.rm = TRUE)
  stage_pos$n_cells[s] <- length(cells_s)
}
write.csv(stage_pos, "01_stage_position.csv", row.names = FALSE)

p_pt <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) +
  annotate("label", x = stage_pos$umap_1, y = stage_pos$umap_2, label = stage_pos$stage,
           size = 3, fontface = "bold", fill = "white", colour = "black",
           linewidth = 0.3, alpha = 0.9) +
  labs(title = paste0(target_celltype, " pseudotime")) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
ggsave("02_trajectory_pseudotime.pdf", p_pt, width = 8, height = 7, units = "cm", dpi = 600)
ggsave("02_trajectory_pseudotime.png", p_pt, width = 8, height = 7, units = "cm", dpi = 600)

subtype_levels <- setdiff(sort(unique(as.character(colData(cds)$labels))),
                          c("Neurons", "Endothelial_cells"))
subtype_cols <- c("Fibroblasts" = "#F4A582", "Tissue_stem_cells" = "#C0392B",
                  "Smooth_muscle_cells" = "#7B3294", "Chondrocytes" = "#4D9221",
                  "MSC" = "#00A087", "Neurons" = "#BDBDBD", "Endothelial_cells" = "#2166AC")
p_ct <- plot_cells(cds, color_cells_by = "labels", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) +
  scale_color_manual(values = subtype_cols, name = "Stromal subtype") +
  labs(title = "Stromal subtypes along pseudotime") +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 6), legend.title = element_text(size = 7))
ggsave("03_trajectory_subtype.pdf", p_ct, width = 9, height = 7, units = "cm", dpi = 600)
ggsave("03_trajectory_subtype.png", p_ct, width = 9, height = 7, units = "cm", dpi = 600)

avail_genes <- intersect(key_genes, rownames(cds))
p_gene <- plot_genes_in_pseudotime(cds[avail_genes, ], color_cells_by = "pseudotime",
                                   min_expr = 0.1, ncol = 2)
ggsave("04_keygene_pseudotime.pdf", p_gene, width = 12, height = 9, units = "cm", dpi = 600)
ggsave("04_keygene_pseudotime.png", p_gene, width = 12, height = 9, units = "cm", dpi = 600)

plate <- (p_pt | p_ct) / p_gene +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("05_plate_7jk.pdf", plate, width = 15.5, height = 16, units = "cm", dpi = 600)
ggsave("05_plate_7jk.png", plate, width = 15.5, height = 16, units = "cm", dpi = 600)


