library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

seu_list <- readRDS(file.path(project_dir, "01load_qc", "01_seu_list_qc.rds"))

for (i in 1:length(seu_list)) {
  seu <- seu_list[[i]]
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  seu_list[[i]] <- seu
}
saveRDS(seu_list, "02_seu_list_cluster.rds")

sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
clu_count <- data.frame()
spat <- data.frame()
umap <- data.frame()

for (i in 1:length(seu_list)) {
  seu <- seu_list[[i]]
  coords <- GetTissueCoordinates(seu)
  emb <- Embeddings(seu, "umap")
  cluster <- as.character(seu$seurat_clusters)

  clu_count <- rbind(clu_count, data.frame(
    sample = seu$section[1], group = seu$group[1],
    cluster = names(table(cluster)), n_spot = as.integer(table(cluster))))

  spat <- rbind(spat, data.frame(barcode = colnames(seu), section = seu$section[1],
                                 group = seu$group[1], x = coords$x, y = coords$y,
                                 cluster = cluster))
  umap <- rbind(umap, data.frame(barcode = colnames(seu), section = seu$section[1],
                                 group = seu$group[1], umap1 = emb[, 1], umap2 = emb[, 2],
                                 cluster = cluster))
}
write.csv(clu_count, "02_cluster_count.csv", row.names = FALSE)
write.csv(spat, "02_spot_cluster.csv", row.names = FALSE)
write.csv(umap, "02_umap_coords.csv", row.names = FALSE)

clu_lab <- sort(unique(as.integer(spat$cluster)))
spat$section <- factor(spat$section, levels = sheet$section)
spat$cluster <- factor(spat$cluster, levels = clu_lab)
umap$section <- factor(umap$section, levels = sheet$section)
umap$cluster <- factor(umap$cluster, levels = clu_lab)

p_spat <- ggplot(spat, aes(x, y, colour = cluster)) +
  geom_point(size = 0.4) +
  scale_colour_manual(name = "Cluster",
                      values = c("#2166AC", "#C0392B", "#4D9221", "#7B3294", "#F4A582",
                                 "#00A087", "#E0C097", "#3C5488", "#B2182B", "#8AB17D",
                                 "#E76F51", "#BDBDBD")) +
  facet_wrap(~section, scales = "free", nrow = 2) +
  scale_y_reverse() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(face = "bold", margin = margin(b = 1)),
        legend.key.size = unit(0.28, "cm"),
        panel.spacing = unit(1.5, "pt"))
ggsave("02_cluster_spatial.pdf", p_spat, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("02_cluster_spatial.png", p_spat, width = 15.5, height = 12, units = "cm", dpi = 600)
saveRDS(p_spat, "02_p_spat.rds")

p_umap <- ggplot(umap, aes(umap1, umap2, colour = cluster)) +
  geom_point(size = 0.35) +
  scale_colour_manual(name = "Cluster",
                      values = c("#2166AC", "#C0392B", "#4D9221", "#7B3294", "#F4A582",
                                 "#00A087", "#E0C097", "#3C5488", "#B2182B", "#8AB17D",
                                 "#E76F51", "#BDBDBD")) +
  facet_wrap(~section, scales = "free", nrow = 2) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"))
ggsave("02_cluster_umap.pdf", p_umap, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("02_cluster_umap.png", p_umap, width = 15.5, height = 12, units = "cm", dpi = 600)

plate <- p_spat / p_umap +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("02_cluster_plate.pdf", plate, width = 15.5, height = 24, units = "cm", dpi = 600)
ggsave("02_cluster_plate.png", plate, width = 15.5, height = 24, units = "cm", dpi = 600)

