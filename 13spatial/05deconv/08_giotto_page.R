library(Seurat)
library(Giotto)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

options(giotto.no_python_warn = TRUE)
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

sc <- readRDS(file.path(project_dir, "..", "01singleCell", "05_scData_annotated.rds"))
sc <- subset(sc, cells = colnames(sc)[!is.na(sc$external_annotation)])
DefaultAssay(sc) <- "RNA"
sc <- NormalizeData(sc, verbose = FALSE)
Idents(sc) <- "external_annotation"

markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.5, verbose = FALSE)
marker_list <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  summarise(genes = list(gene)) %>%
  tibble::deframe()

sign_mat <- makeSignMatrixPAGE(sign_names = names(marker_list), sign_list = marker_list)

seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))
score_tab <- data.frame()
spat <- data.frame()

for (i in 1:length(seu_list)) {
  sec <- names(seu_list)[i]
  seu <- seu_list[[i]]
  counts <- GetAssayData(seu, assay = "Spatial", layer = "counts")

  coords <- GetTissueCoordinates(seu)
  coords <- as.data.frame(coords[colnames(counts), 1:2])
  colnames(coords) <- c("sdimx", "sdimy")
  coords$cell_ID <- colnames(counts)

  gobj <- createGiottoObject(expression = list("raw" = counts), spatial_locs = coords)
  gobj <- normalizeGiotto(gobj, norm_methods = "standard", verbose = FALSE)
  gobj <- createSpatialNetwork(gobject = gobj, method = "kNN", k = 6,
                               name = "knn_network", verbose = FALSE)
  gobj <- runPAGEEnrich(gobject = gobj, sign_matrix = sign_mat,
                        expression_values = "normalized", min_overlap_genes = 5)

  enr <- as.data.frame(getSpatialEnrichment(gobj, name = "PAGE", output = "data.table"))
  rownames(enr) <- enr$cell_ID
  enr$cell_ID <- NULL

  for (ct in colnames(enr)) {
    score_tab <- rbind(score_tab, data.frame(
      section = sec, group = seu$group[1], celltype = ct,
      mean_score = round(mean(enr[[ct]]), 4), median_score = round(median(enr[[ct]]), 4)))
  }
  spat <- rbind(spat, data.frame(barcode = rownames(enr), section = sec,
                                 group = seu$group[1],
                                 x = coords$sdimx, y = coords$sdimy, enr))
}
write.csv(score_tab, "08_page_score_summary.csv", row.names = FALSE)
write.csv(spat, "08_page_score_spot.csv", row.names = FALSE)

sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
celltypes <- sort(unique(score_tab$celltype))
score_tab$section <- factor(score_tab$section, levels = sheet$section)
score_tab$celltype <- factor(score_tab$celltype, levels = celltypes)
max_score <- max(score_tab$mean_score)
score_tab$lab_col <- ifelse(score_tab$mean_score > max_score * 0.55, "white", "grey20")

p_heat <- ggplot(score_tab, aes(celltype, section, fill = mean_score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", mean_score), colour = lab_col), size = 7 / .pt) +
  scale_colour_identity() +
  scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                       limits = c(min(score_tab$mean_score), max_score),
                       oob = squish, name = "PAGE score") +
  scale_y_discrete(limits = rev(levels(score_tab$section))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.6, "cm"))
ggsave("08_page_heatmap.pdf", p_heat, width = 11, height = 7, units = "cm", dpi = 600)
ggsave("08_page_heatmap.png", p_heat, width = 11, height = 7, units = "cm", dpi = 600)

best_ct <- as.character(score_tab$celltype[which.max(score_tab$mean_score)])
spat$section <- factor(spat$section, levels = sheet$section)
spat$facet <- factor(paste0(spat$section, "\n", spat$group),
                     levels = paste0(sheet$section, "\n", sheet$group))

p_spat <- ggplot(spat, aes(x, y, colour = .data[[best_ct]])) +
  geom_point(size = 0.4) +
  scale_colour_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                         limits = c(min(spat[[best_ct]]), max(spat[[best_ct]])),
                         oob = squish, name = "PAGE score") +
  facet_wrap(~facet, scales = "free", nrow = 2) +
  scale_y_reverse() +
  labs(x = NULL, y = NULL, title = paste("Giotto PAGE enrichment:", best_ct)) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(size = 7, face = "bold", lineheight = 0.9, margin = margin(b = 2)),

        plot.title = element_text(size = 9, face = "bold"),
        legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.6, "cm"),
        panel.spacing = unit(1.5, "pt"))
ggsave("08_page_spatial.pdf", p_spat, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("08_page_spatial.png", p_spat, width = 15.5, height = 12, units = "cm", dpi = 600)
saveRDS(p_heat, "08_p_heat.rds")
saveRDS(p_spat, "08_p_spat.rds")

plate <- p_spat / p_heat +
  plot_layout(heights = c(1.6, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"),
        plot.tag.position = c(0, 1))
ggsave("08_page_plate.pdf", plate, width = 15.5, height = 19, units = "cm", dpi = 600)
ggsave("08_page_plate.png", plate, width = 15.5, height = 19, units = "cm", dpi = 600)

