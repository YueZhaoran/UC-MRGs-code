library(ggplot2)
library(patchwork)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

p_expr <- readRDS(file.path(project_dir, "04keygene", "04_p_expr.rds"))
p_det <- readRDS(file.path(project_dir, "04keygene", "04_p_det.rds"))
p_cor <- readRDS(file.path(project_dir, "04keygene", "04_p_cor.rds"))
p_dom <- readRDS(file.path(project_dir, "05deconv", "06_p_dom.rds"))
p_qc <- readRDS(file.path(project_dir, "01load_qc", "01_p_bar.rds"))
p_clu <- readRDS(file.path(project_dir, "02cluster", "02_p_spat.rds"))

qc <- read.csv(file.path(project_dir, "01load_qc", "01_qc_summary.csv"))
stat <- read.csv(file.path(project_dir, "04keygene", "04_keygene_expr_stat.csv"))
cor_tab <- read.csv(file.path(project_dir, "04keygene", "04_keygene_celltype_cor.csv"))
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
key_genes <- unique(stat$gene)

fig10b <- p_expr / (p_det | p_cor) +
  plot_layout(heights = c(2.4, 1.5), widths = c(1, 1.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"),
        plot.tag.position = c(0, 1))
ggsave("05_keygene_plate.pdf", fig10b, width = 15.5, height = 21, units = "cm", dpi = 600)
ggsave("05_keygene_plate.png", fig10b, width = 15.5, height = 21, units = "cm", dpi = 600)

cor_sel <- cor_tab[cor_tab$gene %in% c("SPARC", "FBN1") &
                     cor_tab$celltype %in% c("Epithelial", "Stromal"), ]
cor_txt <- sprintf("%s vs %s rho = %.2f, P = %.3f",
                   cor_sel$gene, cor_sel$celltype, cor_sel$spearman_rho, cor_sel$p_value)

supp <- p_clu / p_dom / p_qc +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("06_supp_plate.pdf", supp, width = 15.5, height = 26, units = "cm", dpi = 600)
ggsave("06_supp_plate.png", supp, width = 15.5, height = 26, units = "cm", dpi = 600)

