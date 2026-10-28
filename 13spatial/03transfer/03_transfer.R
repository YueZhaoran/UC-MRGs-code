library(Seurat)
library(ggplot2)
library(patchwork)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

sc <- readRDS(file.path(project_dir, "..", "01singleCell", "05_scData_annotated.rds"))

sc <- subset(sc, cells = colnames(sc)[!is.na(sc$external_annotation)])

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose = FALSE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

sc <- ScaleData(sc, verbose = FALSE)

seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))

for (i in 1:length(seu_list)) {
  seu <- seu_list[[i]]
  DefaultAssay(seu) <- "Spatial"

  anchors <- FindTransferAnchors(reference = sc, query = seu, normalization.method = "LogNormalize",
                                 reference.assay = "RNA", query.assay = "Spatial", verbose = FALSE)
  pred <- TransferData(anchorset = anchors, refdata = sc$external_annotation,
                       weight.reduction = seu[["pca"]], dims = 1:30, verbose = FALSE)
  seu <- AddMetaData(seu, metadata = pred[colnames(seu), , drop = FALSE])
  seu_list[[i]] <- seu
}
saveRDS(seu_list, "03_seu_list_transfer.rds")

sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
score_cols <- grep("^prediction.score\\.", colnames(seu_list[[1]]@meta.data), value = TRUE)

lab_count <- data.frame()
spat <- data.frame()

for (i in 1:length(seu_list)) {
  seu <- seu_list[[i]]
  coords <- GetTissueCoordinates(seu)
  score <- seu@meta.data[, score_cols]
  best <- max.col(score, ties.method = "first")

  lab_count <- rbind(lab_count, data.frame(
    sample = seu$section[1], group = seu$group[1],
    celltype = names(table(seu$predicted.id)), n_spot = as.integer(table(seu$predicted.id))))

  spat <- rbind(spat, data.frame(
    barcode = colnames(seu), section = seu$section[1], group = seu$group[1],
    x = coords$x, y = coords$y, celltype = seu$predicted.id,
    score = score[cbind(1:nrow(score), best)]))
}
write.csv(lab_count, "03_label_count.csv", row.names = FALSE)
write.csv(spat, "03_spot_celltype.csv", row.names = FALSE)

cell_cols <- c("B/Plasma" = "#2166AC", "T cell" = "#C0392B", "Epithelial" = "#4D9221",
               "Myeloid" = "#7B3294", "Stromal" = "#F4A582", "NK/ILC" = "#00A087",
               "Unassigned" = "#BDBDBD")
spat$section <- factor(spat$section, levels = sheet$section)
spat$celltype <- factor(spat$celltype, levels = names(cell_cols))
lab_count$section <- factor(lab_count$sample, levels = sheet$section)
lab_count$celltype <- factor(lab_count$celltype, levels = names(cell_cols))

p_label <- ggplot(spat, aes(x, y, colour = celltype)) +
  geom_point(size = 0.4) +
  scale_colour_manual(values = cell_cols, name = "Cell type", drop = FALSE) +
  facet_wrap(~section, scales = "free", nrow = 2) +
  scale_y_reverse() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(face = "bold", margin = margin(b = 1)),
        legend.key.size = unit(0.28, "cm"),
        panel.spacing = unit(1.5, "pt"))
ggsave("03_label_spatial.pdf", p_label, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("03_label_spatial.png", p_label, width = 15.5, height = 12, units = "cm", dpi = 600)

p_bar <- ggplot(lab_count, aes(x = section, y = n_spot, fill = celltype)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cell_cols, name = "Cell type", drop = FALSE) +
  labs(x = NULL, y = "Number of spots") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("03_label_bar.pdf", p_bar, width = 8, height = 7, units = "cm", dpi = 600)
ggsave("03_label_bar.png", p_bar, width = 8, height = 7, units = "cm", dpi = 600)

p_score <- ggplot(spat, aes(x = section, y = score, fill = group)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  scale_fill_manual(values = c(UC = "#C0392B", Normal = "#2166AC"), name = "Group") +
  labs(x = NULL, y = "Max prediction score") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("03_label_score.pdf", p_score, width = 8, height = 7, units = "cm", dpi = 600)
ggsave("03_label_score.png", p_score, width = 8, height = 7, units = "cm", dpi = 600)

plate <- p_label / (p_bar | p_score) +
  plot_layout(heights = c(1.5, 1), widths = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("03_label_plate.pdf", plate, width = 15.5, height = 18, units = "cm", dpi = 600)
ggsave("03_label_plate.png", plate, width = 15.5, height = 18, units = "cm", dpi = 600)

