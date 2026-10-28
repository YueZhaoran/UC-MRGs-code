library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
seu_list <- list()
qc_tab <- data.frame()

for (i in 1:nrow(sheet)) {
  sid <- paste0(sheet$gsm[i], "_", sheet$section[i])
  seu <- Load10X_Spatial(file.path(project_dir, "00raw", sid),
                         filename = "filtered_feature_bc_matrix.h5")
  seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$section <- sheet$section[i]
  seu$group <- sheet$group[i]

  keep <- seu$nFeature_Spatial >= 200 & seu$nFeature_Spatial <= 6000 &
    seu$nCount_Spatial >= 500 & seu$nCount_Spatial <= 30000 &
    seu$percent.mt <= 20
  n_raw <- ncol(seu)
  seu <- subset(seu, cells = colnames(seu)[keep])

  qc_tab <- rbind(qc_tab, data.frame(
    sample = sid, section = sheet$section[i], group = sheet$group[i],
    n_spot_raw = n_raw, n_spot_keep = ncol(seu),
    keep_rate = round(ncol(seu) / n_raw * 100, 1),
    n_feature_median = median(seu$nFeature_Spatial),
    n_count_median = median(seu$nCount_Spatial),
    pct_mt_median = round(median(seu$percent.mt), 2)
  ))

  seu_list[[sheet$section[i]]] <- seu
}
write.csv(qc_tab, "01_qc_summary.csv", row.names = FALSE)

plot_tab <- data.frame()
for (i in 1:nrow(sheet)) {
  raw <- Load10X_Spatial(file.path(project_dir, "00raw",
                                   paste0(sheet$gsm[i], "_", sheet$section[i])),
                         filename = "filtered_feature_bc_matrix.h5")
  raw$percent.mt <- PercentageFeatureSet(raw, pattern = "^MT-")

  plot_tab <- rbind(plot_tab,
    data.frame(section = sheet$section[i], group = sheet$group[i], qc = "Before QC",
               nFeature = raw$nFeature_Spatial, nCount = raw$nCount_Spatial,
               percent.mt = raw$percent.mt),
    data.frame(section = sheet$section[i], group = sheet$group[i], qc = "After QC",
               nFeature = seu_list[[sheet$section[i]]]$nFeature_Spatial,
               nCount = seu_list[[sheet$section[i]]]$nCount_Spatial,
               percent.mt = seu_list[[sheet$section[i]]]$percent.mt)
  )
}
plot_tab$section <- factor(plot_tab$section, levels = sheet$section)
plot_tab$qc <- factor(plot_tab$qc, levels = c("Before QC", "After QC"))

plot_long <- pivot_longer(plot_tab, c(nFeature, nCount, percent.mt),
                          names_to = "metric", values_to = "value")
plot_long$metric <- factor(plot_long$metric,
                           levels = c("nFeature", "nCount", "percent.mt"),
                           labels = c("nFeature per spot", "nCount per spot", "Mitochondrial %"))

p_dist <- ggplot(plot_long, aes(section, value, fill = qc)) +
  geom_violin(scale = "width", position = position_dodge(0.8), linewidth = 0.2) +
  scale_fill_manual(values = c("Before QC" = "#999999", "After QC" = "#7B3294"), name = NULL) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "top", legend.margin = margin(0, 0, -4, 0),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"))
ggsave("01_qc_metric_violin.pdf", p_dist, width = 15.5, height = 7, units = "cm", dpi = 600)
ggsave("01_qc_metric_violin.png", p_dist, width = 15.5, height = 7, units = "cm", dpi = 600)
saveRDS(p_dist, "01_p_dist.rds")

spat_tab <- data.frame()
for (i in 1:nrow(sheet)) {
  coords <- GetTissueCoordinates(seu_list[[sheet$section[i]]])
  spat_tab <- rbind(spat_tab, data.frame(
    section = sheet$section[i], group = sheet$group[i],
    x = coords$x, y = coords$y,
    nFeature = seu_list[[sheet$section[i]]]$nFeature_Spatial))
}
spat_tab$section <- factor(spat_tab$section, levels = sheet$section)
spat_tab$facet <- factor(paste0(spat_tab$section, "\n", spat_tab$group),
                         levels = paste0(sheet$section, "\n", sheet$group))
write.csv(spat_tab, "01_spot_metadata.csv", row.names = FALSE)

p_spat <- ggplot(spat_tab, aes(x, y, colour = nFeature)) +
  geom_point(size = 0.35) +
  scale_colour_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                         limits = c(0, 3500), oob = squish, name = "nFeature") +
  facet_wrap(~facet, scales = "free", nrow = 2) +
  scale_y_reverse() +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(face = "bold", lineheight = 0.9),
        legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.6, "cm"),
        panel.spacing = unit(1.5, "pt"))
ggsave("01_qc_spatial_scatter.pdf", p_spat, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("01_qc_spatial_scatter.png", p_spat, width = 15.5, height = 12, units = "cm", dpi = 600)

bar_tab <- pivot_longer(qc_tab, c(n_spot_raw, n_spot_keep),
                        names_to = "qc", values_to = "n")
bar_tab$qc <- factor(bar_tab$qc, levels = c("n_spot_raw", "n_spot_keep"),
                     labels = c("Before QC", "After QC"))
bar_tab$section <- factor(bar_tab$section, levels = sheet$section)

p_bar <- ggplot(bar_tab, aes(section, n, fill = qc)) +
  geom_col(position = position_dodge(0.75), width = 0.68) +
  geom_text(aes(label = n, y = n + 100), position = position_dodge(0.75),
            angle = 90, hjust = 0, size = 7 / .pt, show.legend = FALSE) +
  scale_fill_manual(values = c("Before QC" = "#999999", "After QC" = "#7B3294"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of spots") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.margin = margin(0, 0, -4, 0))
ggsave("01_qc_spot_count.pdf", p_bar, width = 7.5, height = 7, units = "cm", dpi = 600)
ggsave("01_qc_spot_count.png", p_bar, width = 7.5, height = 7, units = "cm", dpi = 600)
saveRDS(p_bar, "01_p_bar.rds")

plate <- p_dist / (p_spat | p_bar) +
  plot_layout(heights = c(1, 1.4), widths = c(1.7, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"),
        plot.tag.position = c(0, 1))
ggsave("01_qc_plate.pdf", plate, width = 15.5, height = 16, units = "cm", dpi = 600)
ggsave("01_qc_plate.png", plate, width = 15.5, height = 16, units = "cm", dpi = 600)

saveRDS(seu_list, "01_seu_list_qc.rds")

