library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

key_genes <- read.csv(file.path(project_dir, "..", "05DEGs", "hub_gene.csv"))$symbol
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))

stat <- data.frame()
expr <- data.frame()

for (i in 1:length(seu_list)) {
  seu <- seu_list[[i]]
  DefaultAssay(seu) <- "Spatial"
  seu <- NormalizeData(seu, verbose = FALSE)
  coords <- GetTissueCoordinates(seu)
  dat <- GetAssayData(seu, assay = "Spatial", layer = "data")

  for (g in key_genes) {
    v <- as.numeric(dat[g, ])
    stat <- rbind(stat, data.frame(
      sample = seu$section[1], group = seu$group[1], gene = g,
      n_spot = ncol(seu), n_spot_detected = sum(v > 0),
      pct_detected = round(sum(v > 0) / ncol(seu) * 100, 1),
      mean_expr = round(mean(v), 3), max_expr = round(max(v), 3)))
    expr <- rbind(expr, data.frame(
      sample = seu$section[1], group = seu$group[1], gene = g,
      x = coords$x, y = coords$y, expr = v))
  }
}
write.csv(stat, "04_keygene_expr_stat.csv", row.names = FALSE)
write.csv(expr, "04_keygene_expr_spot.csv", row.names = FALSE)

prop <- read.csv(file.path(project_dir, "05deconv", "06_prop_summary.csv"))
prop_wide <- reshape(prop[, c("sample", "celltype", "mean_prop")],
                     idvar = "sample", timevar = "celltype", direction = "wide")
names(prop_wide) <- sub("^mean_prop\\.", "", names(prop_wide))
celltypes <- setdiff(names(prop_wide), "sample")

cor_tab <- data.frame()
for (g in key_genes) {
  d <- merge(stat[stat$gene == g, ], prop_wide, by = "sample")
  for (ct in celltypes) {
    test <- cor.test(d$pct_detected, d[[ct]], method = "spearman")
    cor_tab <- rbind(cor_tab, data.frame(
      gene = g, celltype = ct, n_section = nrow(d),
      spearman_rho = round(unname(test$estimate), 3),
      p_value = signif(test$p.value, 3)))
  }
}
write.csv(cor_tab, "04_keygene_celltype_cor.csv", row.names = FALSE)

plot_genes <- c("SPARC", "FBN1")
plot_ct <- c("Epithelial", "Stromal")
cor_sel <- cor_tab[cor_tab$gene %in% plot_genes & cor_tab$celltype %in% plot_ct, ]

cor_long <- merge(stat[stat$gene %in% plot_genes, c("sample", "group", "gene", "pct_detected")],
                  prop[prop$celltype %in% plot_ct, c("sample", "celltype", "mean_prop")],
                  by = "sample")
cor_long <- merge(cor_long, cor_sel, by = c("gene", "celltype"))
cor_long$gene <- factor(cor_long$gene, levels = plot_genes)
cor_long$group <- factor(cor_long$group, levels = c("UC", "Normal"))
cor_long$panel <- factor(sprintf("%s\nrho = %.2f, P = %.3f", cor_long$celltype,
                                 cor_long$spearman_rho, cor_long$p_value),
                         levels = unique(sprintf("%s\nrho = %.2f, P = %.3f",
                                                 cor_sel$celltype, cor_sel$spearman_rho, cor_sel$p_value)))

expr$section <- factor(expr$sample, levels = sheet$section)
expr$gene <- factor(expr$gene, levels = key_genes)
expr$facet <- factor(paste0(expr$section, "\n", expr$group),
                     levels = paste0(sheet$section, "\n", sheet$group))

p_expr <- ggplot(expr, aes(x, y, colour = expr)) +
  geom_point(size = 0.3) +
  scale_colour_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                         limits = c(0, 3), oob = squish, breaks = c(0, 1, 2, 3),
                         name = "Expression") +
  facet_grid(gene ~ facet, scales = "free", switch = "y") +
  scale_y_reverse() +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text.x = element_text(size = 7, face = "bold", lineheight = 0.9,
                                    margin = margin(b = 2)),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.placement = "outside",
        legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.55, "cm"),
        panel.spacing = unit(1.2, "pt"), plot.margin = margin(1, 1, 1, 1))
ggsave("04_keygene_spatial.pdf", p_expr, width = 15.5, height = 14, units = "cm", dpi = 600)
ggsave("04_keygene_spatial.png", p_expr, width = 15.5, height = 14, units = "cm", dpi = 600)

stat$gene <- factor(stat$gene, levels = key_genes)
stat$group <- factor(stat$group, levels = c("UC", "Normal"))
stat_mean <- aggregate(pct_detected ~ gene, stat, mean)
stat_med <- aggregate(pct_detected ~ gene + group, stat, median)

p_det <- ggplot(stat, aes(gene, pct_detected)) +
  geom_col(data = stat_mean, aes(gene, pct_detected), inherit.aes = FALSE,
           fill = "grey93", width = 0.6) +
  geom_point(aes(fill = group), shape = 21, size = 2.1, colour = "white", stroke = 0.4,
             position = position_jitter(width = 0.11, height = 0, seed = 123)) +
  geom_segment(data = stat_med,
               aes(x = as.numeric(gene) - 0.22, xend = as.numeric(gene) + 0.22,
                   y = pct_detected, yend = pct_detected, colour = group),
               linewidth = 0.55, inherit.aes = FALSE) +
  scale_fill_manual(values = c(UC = "#C0392B", Normal = "#2166AC"), name = "Group") +
  scale_colour_manual(values = c(UC = "#C0392B", Normal = "#2166AC"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Spots with detectable expression (%)") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold"),
        legend.position = "top", legend.margin = margin(0, 0, -4, 0))
ggsave("04_keygene_detection.pdf", p_det, width = 8, height = 7, units = "cm", dpi = 600)
ggsave("04_keygene_detection.png", p_det, width = 8, height = 7, units = "cm", dpi = 600)

p_cor <- ggplot(cor_long, aes(mean_prop, pct_detected)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_point(aes(fill = group), shape = 21, size = 2, colour = "white", stroke = 0.35) +
  facet_grid(gene ~ panel, scales = "free") +
  scale_fill_manual(values = c(UC = "#C0392B", Normal = "#2166AC"), name = "Group") +
  scale_x_continuous(labels = percent_format(accuracy = 1), n.breaks = 3) +
  labs(x = "Mean cell-type proportion per spot", y = "Spots with detectable expression (%)") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7), axis.title = element_text(size = 9),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 7, face = "bold", lineheight = 1.05,
                                  margin = margin(2, 1, 2, 1)),
        strip.text.y = element_text(angle = 0, face = "bold"),
        panel.spacing = unit(2.5, "pt"),
        legend.position = "top", legend.margin = margin(0, 0, -4, 0))
ggsave("04_keygene_celltype_cor.pdf", p_cor, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("04_keygene_celltype_cor.png", p_cor, width = 15.5, height = 12, units = "cm", dpi = 600)

saveRDS(p_expr, "04_p_expr.rds")
saveRDS(p_det, "04_p_det.rds")
saveRDS(p_cor, "04_p_cor.rds")

plate <- p_expr / (p_det | p_cor) +
  plot_layout(heights = c(2.4, 1.5), widths = c(1, 1.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"),
        plot.tag.position = c(0, 1))
ggsave("04_keygene_plate.pdf", plate, width = 15.5, height = 21, units = "cm", dpi = 600)
ggsave("04_keygene_plate.png", plate, width = 15.5, height = 21, units = "cm", dpi = 600)

top_cor <- cor_sel[order(-abs(cor_sel$spearman_rho)), ]
cor_txt <- sprintf("%s vs %s rho = %.2f, P = %.3f",
                   top_cor$gene, top_cor$celltype, top_cor$spearman_rho, top_cor$p_value)

