library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

prop <- readRDS(file.path(project_dir, "05deconv", "06_spotlight_prop.rds"))
prop_sum <- read.csv(file.path(project_dir, "05deconv", "06_prop_summary.csv"))
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))

cell_cols <- c("B/Plasma" = "#2166AC", "T cell" = "#C0392B", "Epithelial" = "#4D9221",
               "Myeloid" = "#7B3294", "Stromal" = "#F4A582", "NK/ILC" = "#00A087",
               "Unassigned" = "#BDBDBD")
celltypes <- names(cell_cols)

coords <- data.frame()
for (i in 1:length(seu_list)) {
  co <- GetTissueCoordinates(seu_list[[i]])
  coords <- rbind(coords, data.frame(sample = seu_list[[i]]$section[1],
                                     barcode = colnames(seu_list[[i]]), x = co$x, y = co$y))
}
spat <- merge(prop, coords, by = c("sample", "barcode"))
spat$section <- factor(spat$sample, levels = sheet$section)
spat$facet <- factor(paste0(spat$sample, "\n", spat$group),
                     levels = paste0(sheet$section, "\n", sheet$group))
spat$dominant <- factor(spat$dominant, levels = celltypes)
prop_sum$section <- factor(prop_sum$sample, levels = sheet$section)
prop_sum$celltype <- factor(prop_sum$celltype, levels = celltypes)

p_dom <- ggplot(spat, aes(x, y, colour = dominant)) +
  geom_point(size = 0.4) +
  scale_colour_manual(values = cell_cols, name = "Dominant\ncell type", drop = FALSE) +
  facet_wrap(~facet, scales = "free", nrow = 2) +
  scale_y_reverse() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(size = 7, face = "bold", lineheight = 0.9, margin = margin(b = 2)),

        legend.key.size = unit(0.28, "cm"), panel.spacing = unit(1.5, "pt"))
ggsave("06_deconv_dominant.pdf", p_dom, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("06_deconv_dominant.png", p_dom, width = 15.5, height = 12, units = "cm", dpi = 600)

p_stack <- ggplot(prop_sum, aes(section, mean_prop, fill = celltype)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = cell_cols, name = "Cell type", drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean proportion per spot") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.28, "cm"))
ggsave("06_deconv_stack.pdf", p_stack, width = 8.5, height = 7, units = "cm", dpi = 600)
ggsave("06_deconv_stack.png", p_stack, width = 8.5, height = 7, units = "cm", dpi = 600)

prop_sum$lab_col <- ifelse(prop_sum$mean_prop > max(prop_sum$mean_prop) * 0.62, "white", "grey20")

p_heat <- ggplot(prop_sum, aes(celltype, section, fill = mean_prop)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", mean_prop), colour = lab_col), size = 7 / .pt) +
  scale_colour_identity() +
  scale_fill_gradientn(colours = c("#F7F7F7", "#F4A582", "#B2182B"),
                       limits = c(0, max(prop_sum$mean_prop)), oob = squish, name = "Proportion") +
  scale_y_discrete(limits = rev(levels(prop_sum$section))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.6, "cm"))
ggsave("06_deconv_heatmap.pdf", p_heat, width = 9.5, height = 7, units = "cm", dpi = 600)
ggsave("06_deconv_heatmap.png", p_heat, width = 9.5, height = 7, units = "cm", dpi = 600)

saveRDS(p_dom, "06_p_dom.rds")
saveRDS(p_stack, "06_p_stack.rds")
saveRDS(p_heat, "06_p_heat.rds")

plate <- p_dom / (p_stack | p_heat) +
  plot_layout(heights = c(1.5, 1), widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"),
        plot.tag.position = c(0, 1))
ggsave("06_deconv_plate.pdf", plate, width = 15.5, height = 19, units = "cm", dpi = 600)
ggsave("06_deconv_plate.png", plate, width = 15.5, height = 19, units = "cm", dpi = 600)

top_ct <- prop_sum[prop_sum$mean_prop == max(prop_sum$mean_prop), ]

