rm(list = ls()); gc()

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  .args <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .args[grep("^--file=", .args)])
  setwd(if (length(.f) > 0) dirname(normalizePath(.f)) else getwd())
}

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

key_genes <- c("ADAMTS1", "FBN1", "SPARC", "VCAM1")
top_lr    <- 6
top_path  <- 8

meta_lines <- readLines("03_run_meta.txt", encoding = "UTF-8")
mv <- setNames(sub("^[^=]+=", "", meta_lines), sub("=.*$", "", meta_lines))

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
p_floor_used <- NA_real_

for (g in key_genes) {
  d <- read.csv(paste0(g, "_KO_all_genes_sorted_by_pvalue.csv"))
  d$logFC <- ifelse(d$FC > 0, log2(d$FC), NA)
  p_pos <- d$p.value[d$p.value > 0]
  p_floor <- min(p_pos, na.rm = TRUE) / 10
  p_floor_used <- if (is.na(p_floor_used)) p_floor else min(p_floor_used, p_floor)
  d$neg_log10_pval <- -log10(pmax(d$p.value, p_floor))
  d$clipped <- d$p.value <= p_floor

  d$is_drg <- !is.na(d$logFC) & d$p.value < 0.05 & abs(d$logFC) > 0.5

  sig_tab <- rbind(sig_tab, data.frame(
    gene = g, n_gene = nrow(d),
    n_pass_p005 = sum(d$p.value < 0.05, na.rm = TRUE),
    n_pass_drg = sum(d$is_drg, na.rm = TRUE),
    n_pass_padj005 = sum(d$p.adj < 0.05, na.rm = TRUE),
    n_pass_padj005_logfc05 = sum(d$p.adj < 0.05 & abs(d$logFC) > 0.5, na.rm = TRUE),
    n_pass_Z2 = sum(abs(d$Z) > 2, na.rm = TRUE),
    n_pass_padj001 = sum(d$p.adj < 0.01, na.rm = TRUE),
    n_pass_both = sum(abs(d$Z) > 2 & d$p.adj < 0.01, na.rm = TRUE),
    n_clipped = sum(d$clipped, na.rm = TRUE)))

  z_rng <- rbind(z_rng,
                 data.frame(gene = g, z_min = min(d$Z, na.rm = TRUE), z_max = max(d$Z, na.rm = TRUE)))

  d_plot <- d[!is.na(d$logFC) & is.finite(d$logFC), ]
  top_lab <- d_plot %>% dplyr::arrange(p.value) %>% head(top_lr)

  p <- ggplot(d_plot, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(colour = is_drg), size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = c(`TRUE` = gene_cols[[g]], `FALSE` = "#BDBDBD"),
                        labels = c(`TRUE` = "P < 0.05", `FALSE` = "NS"),
                        name = NULL) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = top_lab, aes(label = gene), size = 2.2,
                             segment.size = 0.2, max.overlaps = 20, min.segment.length = 0) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                       breaks = c(-20, -5, -1, 0, 1, 5, 20),
                       labels = c("-20", "-5", "-1", "0", "1", "5", "20")) +
    labs(x = expression(log[2]~"(likelihood ratio)"), y = expression(-log[10]~"(P)"),
         title = paste0(g, " knockout in Stromal (UC)")) +
    plate_theme
  plots[[g]] <- p
}

write.csv(sig_tab, "04_ko_threshold_summary.csv", row.names = FALSE)

plate_8a <- wrap_plots(plots, ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("01_plate_8a_volcano.pdf", plate_8a, width = 15.5, height = 13, units = "cm", dpi = 600)
ggsave("01_plate_8a_volcano.png", plate_8a, width = 15.5, height = 13, units = "cm", dpi = 600)

enrich_plate <- function(file_fun, label_fun) {
  all <- data.frame()
  for (g in key_genes) {
    f <- file_fun(g)
    if (!file.exists(f)) next
    e <- read.csv(f)
    if (nrow(e) == 0) next
    e <- e %>% dplyr::arrange(p.adjust) %>% head(top_path) %>%
      mutate(gene = g,
             GeneRatio_Num = sapply(strsplit(GeneRatio, "/"),
                                    function(x) as.numeric(x[1]) / as.numeric(x[2])),
             label = label_fun(g, .))
    all <- rbind(all,
                 e[, c("gene", "Description", "label", "GeneRatio_Num", "p.adjust", "Count", "pvalue")])
  }
  if (nrow(all) == 0) return(NULL)
  all$label <- factor(all$label, levels = rev(unique(all$label[order(all$gene, all$p.adjust)])))
  ggplot(all, aes(x = GeneRatio_Num, y = label)) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    scale_colour_gradient(low = "#FDDBC7",
                          high = "#B2182B",
                          name = expression(-log[10]~"(p.adj)")) +
    scale_size_continuous(range = c(1, 4), name = "Count") +
    labs(x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6))
}

p_go <- enrich_plate(function(g) paste0(g, "_GO_ALL_results.csv"),
                     function(g, e) paste0(g, ": ", e$ONTOLOGY, " - ", e$Description))
if (!is.null(p_go)) {
  ggsave("02_plate_8b_go.pdf", p_go, width = 15.5, height = 14, units = "cm", dpi = 600)
  ggsave("02_plate_8b_go.png", p_go, width = 15.5, height = 14, units = "cm", dpi = 600)
}

p_kegg <- enrich_plate(function(g) paste0(g, "_KEGG_results.csv"),
                       function(g, e) paste0(g, ": ", e$Description))
if (!is.null(p_kegg)) {
  ggsave("03_plate_8c_kegg.pdf", p_kegg, width = 15.5, height = 14, units = "cm", dpi = 600)
  ggsave("03_plate_8c_kegg.png", p_kegg, width = 15.5, height = 14, units = "cm", dpi = 600)
}


