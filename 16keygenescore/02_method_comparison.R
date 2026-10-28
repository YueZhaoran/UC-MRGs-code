library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

long <- read.csv("01_CSV/02_Violin_Combined_Methods_Data.csv")
summary_ct <- read.csv("01_CSV/Summary_Mean_Scores_by_Celltype.csv")

methods <- c("AUCell", "UCell", "singscore", "ssGSEA", "AddModuleScore")
mean_tab <- long %>%
  group_by(external_annotation, Method) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Method, values_from = mean_score) %>%
  as.data.frame()
rownames(mean_tab) <- mean_tab$external_annotation

z_tab <- as.data.frame(t(scale(t(mean_tab[, methods]))))
z_tab$celltype <- rownames(mean_tab)
norm_long <- pivot_longer(z_tab, all_of(methods), names_to = "Method", values_to = "score_z")
norm_long$Method <- factor(norm_long$Method, levels = methods)
write.csv(norm_long, "01_CSV/05_Method_Comparison_Zscore.csv", row.names = FALSE)
order_ct <- summary_ct$external_annotation[order(-summary_ct$KeyGenes_Combined_Mean)]

rank_tab <- apply(mean_tab[, methods], 2, rank)
rho_tab <- cor(rank_tab, method = "spearman")
write.csv(round(rho_tab, 3), "01_CSV/05_Method_Rank_Correlation.csv")

rank_tab <- apply(mean_tab[, methods], 2, rank)
rank_long <- as.data.frame(as.table(rank_tab))
names(rank_long) <- c("celltype", "Method", "rank")
rank_long$celltype <- factor(rank_long$celltype, levels = order_ct)
rank_long$Method <- factor(rank_long$Method, levels = methods)

p_line <- ggplot(rank_long, aes(celltype, Method, fill = rank)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = as.integer(rank)), size = 6.5 / .pt) +
  scale_fill_gradientn(colours = c("#B2182B", "#F7F7F7", "#2166AC"), breaks = c(1, 4, 7),
                       name = "Rank within\nmethod\n(1 = highest)") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.5, "cm"))

p_bar <- ggplot(summary_ct, aes(reorder(external_annotation, KeyGenes_Combined_Mean),
                                KeyGenes_Combined_Mean, fill = KeyGenes_Combined_Mean)) +
  geom_col(width = 0.68) +
  scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), guide = "none") +
  labs(x = NULL, y = "Combined score (mean)") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

rho_long <- as.data.frame(as.table(rho_tab))
p_rho <- ggplot(rho_long, aes(Var2, Var1, fill = Freq)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 6 / .pt) +
  scale_fill_gradientn(colours = c("#F7F7F7", "#F4A582", "#B2182B"), limits = c(0.5, 1),
                       oob = squish, name = "Spearman\nrho") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.5, "cm"))

plate <- p_line / (p_bar | p_rho) +
  plot_layout(heights = c(1.25, 1), widths = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10, family = "sans"),
        plot.tag.position = c(0, 1))
ggsave("02_Plots/05_Method_Comparison_Plate.pdf",
       plate,
       width = 17,
       height = 13,
       units = "cm",
       dpi = 600)
ggsave("02_Plots/05_Method_Comparison_Plate.png",
       plate,
       width = 17,
       height = 13,
       units = "cm",
       dpi = 600)
saveRDS(plate, "02_Plots/05_p_method_plate.rds")

best_tab <- data.frame(Method = methods,
                       best_celltype = apply(mean_tab[, methods],
                                        2,
                                        function(v) rownames(mean_tab)[which.max(v)]),
                       best_mean = round(apply(mean_tab[, methods], 2, max), 4))
write.csv(best_tab, "01_CSV/05_Method_Best_Celltype.csv", row.names = FALSE)


