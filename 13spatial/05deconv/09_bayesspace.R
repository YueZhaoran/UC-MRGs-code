library(Seurat)
library(BayesSpace)
library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)
library(patchwork)
library(mclust)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
key_genes <- read.csv(file.path(project_dir, "..", "05DEGs", "hub_gene.csv"))$symbol
clu_tab <- data.frame()
spot_tab <- data.frame()

for (i in 1:length(seu_list)) {
  sec <- names(seu_list)[i]
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

  clu_tab <- rbind(clu_tab, data.frame(
    section = sec, group = seu$group[1],
    cluster = names(table(sce$spatial.cluster)),
    n_spot = as.integer(table(sce$spatial.cluster))))
  spot_tab <- rbind(spot_tab, data.frame(
    barcode = colnames(sce), section = sec, group = seu$group[1],
    x = colData(sce)$pxl_col_in_fullres, y = colData(sce)$pxl_row_in_fullres,
    cluster = sce$spatial.cluster))

  enh <- data.frame(
    section = sec, group = seu$group[1],
    x = colData(sce_enh)$pxl_col_in_fullres, y = colData(sce_enh)$pxl_row_in_fullres,
    cluster = sce_enh$spatial.cluster,
    t(enh_mat[key_genes, , drop = FALSE]))
  write.csv(enh, paste0("09_enhanced_subspot_", sec, ".csv"), row.names = FALSE)

  saveRDS(list(sce = sce, sce_enhanced = sce_enh), paste0("09_bayesspace_", sec, ".rds"))
}

write.csv(clu_tab, "09_cluster_count.csv", row.names = FALSE)
write.csv(spot_tab, "09_spot_cluster.csv", row.names = FALSE)

clu_cols <- c("#2166AC",
              "#C0392B",
              "#4D9221",
              "#7B3294",
              "#F4A582",
              "#00A087",
              "#E0C097",
              "#3C5488")
spot_tab$section <- factor(spot_tab$section, levels = sheet$section)
spot_tab$cluster <- factor(spot_tab$cluster)
spot_tab$facet <- factor(paste0(spot_tab$section, "\n", spot_tab$group),
                         levels = paste0(sheet$section, "\n", sheet$group))

p_clu <- ggplot(spot_tab, aes(x, y, colour = cluster)) +
  geom_point(size = 0.4) +
  scale_colour_manual(values = clu_cols, name = "BayesSpace\ncluster") +
  facet_wrap(~facet, scales = "free", nrow = 2) +
  scale_y_reverse() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(size = 7, face = "bold", lineheight = 0.9, margin = margin(b = 2)),

        legend.key.size = unit(0.28, "cm"), panel.spacing = unit(1.5, "pt"))
ggsave("09_cluster_spatial.pdf", p_clu, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("09_cluster_spatial.png", p_clu, width = 15.5, height = 12, units = "cm", dpi = 600)

seurat_clu <- read.csv(file.path(project_dir, "02cluster", "02_spot_cluster.csv"))
merged <- merge(spot_tab[, c("barcode", "section", "cluster")],
                seurat_clu[, c("barcode", "section", "cluster")],
                by = c("barcode", "section"), suffixes = c("_bs", "_seurat"))
ari <- mclust::adjustedRandIndex(merged$cluster_bs, merged$cluster_seurat)

