rm(list = ls()); gc()

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

key_genes <- c("ADAMTS1", "FBN1", "SPARC", "VCAM1")
top_lr    <- 6
top_path  <- 8

meta <- read.table("03_run_meta.txt", sep = "=", stringsAsFactors = FALSE, strip.white = TRUE)
mv   <- setNames(meta[[2]], meta[[1]])

gene_cols <- c(ADAMTS1 = "#2166AC", FBN1 = "#4D9221", SPARC = "#C0392B", VCAM1 = "#7B3294")

plate_theme <- theme_classic(base_size = 8, base_family = "sans") +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.key.size = unit(0.28, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

sig_tab <- data.frame()
plots <- list()
z_rng <- NULL

for (g in key_genes) {
  d <- read.csv(paste0(g, "_OE_diffRegulation.csv"))
  d <- d[!is.na(d$logFC) & is.finite(d$logFC) & !is.na(d$neg_log10_pval), ]
  z_rng <- rbind(z_rng,
                 data.frame(gene = g, z_min = min(d$Z, na.rm = TRUE), z_max = max(d$Z, na.rm = TRUE)))

  d$pass_z    <- abs(d$Z) > 2
  d$pass_padj <- d$p.adj < 0.01
  d$pass_both <- d$pass_z & d$pass_padj
  d$pass_p    <- d$p.value < 0.05

  sig_tab <- rbind(sig_tab, data.frame(
    gene = g, n_gene = nrow(d),
    n_pass_Z2 = sum(d$pass_z, na.rm = TRUE),
    n_pass_padj001 = sum(d$pass_padj, na.rm = TRUE),
    n_pass_both = sum(d$pass_both, na.rm = TRUE),
    n_pass_p005 = sum(d$pass_p, na.rm = TRUE)))

  top_lab <- d %>% dplyr::arrange(p.value) %>% head(top_lr)

  p <- ggplot(d, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(colour = pass_p), size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = c(`TRUE` = gene_cols[[g]], `FALSE` = "#BDBDBD"),
                        labels = c(`TRUE` = "P < 0.05", `FALSE` = "NS"),
                        name = NULL) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = top_lab, aes(label = gene), size = 2.2,
                             segment.size = 0.2, max.overlaps = 20, min.segment.length = 0) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                       breaks = c(-80, -40, -10, -1, 0, 1, 5),
                       labels = c("-80", "-40", "-10", "-1", "0", "1", "5")) +
    labs(x = expression(log[2]~"(likelihood ratio)"), y = expression(-log[10]~"(P)"),
         title = paste0(g, " overexpression in Stromal (UC)")) +
    plate_theme
  plots[[g]] <- p
}

write.csv(sig_tab, "01_drg_threshold_summary.csv", row.names = FALSE)

plate_8d <- wrap_plots(plots, ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("01_plate_8d_volcano.pdf", plate_8d, width = 15.5, height = 13, units = "cm", dpi = 600)
ggsave("01_plate_8d_volcano.png", plate_8d, width = 15.5, height = 13, units = "cm", dpi = 600)

kegg_all <- data.frame()
for (g in key_genes) {
  k <- read.csv(paste0(g, "_OE_KEGG_results.csv"))
  k <- k %>% dplyr::arrange(p.adjust) %>% head(top_path) %>%
    mutate(gene = g,
           GeneRatio_Num = sapply(strsplit(GeneRatio, "/"),
                                  function(x) as.numeric(x[1]) / as.numeric(x[2])),
           label = paste0(g, ": ", Description))
  kegg_all <- rbind(kegg_all, k[, c("gene", "Description", "label", "GeneRatio_Num",
                                    "p.adjust", "Count", "pvalue")])
}
write.csv(kegg_all, "02_kegg_top_pathways.csv", row.names = FALSE)

kegg_all$label <- factor(kegg_all$label,
                         levels = rev(unique(kegg_all$label[order(kegg_all$gene, kegg_all$p.adjust)])))

p_kegg <- ggplot(kegg_all, aes(x = GeneRatio_Num, y = label)) +
  geom_point(aes(size = Count, colour = -log10(p.adjust))) +
  scale_colour_gradient(low = "#FDDBC7", high = "#B2182B", name = expression(-log[10]~"(p.adj)")) +
  scale_size_continuous(range = c(1, 4), name = "Count") +
  labs(x = "Gene ratio", y = NULL) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6))
ggsave("02_plate_8e_kegg.pdf", p_kegg, width = 15.5, height = 14, units = "cm", dpi = 600)
ggsave("02_plate_8e_kegg.png", p_kegg, width = 15.5, height = 14, units = "cm", dpi = 600)

z_sparc <- z_rng[z_rng$gene == "SPARC", ]

