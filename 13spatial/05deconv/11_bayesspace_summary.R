library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)
library(mclust)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
files <- paste0("09_bayesspace_", sheet$section, ".rds")
missing <- sheet$section[!file.exists(files)]
if (length(missing) > 0) stop("缺切片结果: ", paste(missing, collapse = ", "))

clu_tab <- data.frame()
spot_tab <- data.frame()
for (i in 1:length(files)) {
  o <- readRDS(files[i])
  sec <- sheet$section[i]
  clu_tab <- rbind(clu_tab, data.frame(
    section = sec, group = sheet$group[i],
    cluster = names(table(o$sce$spatial.cluster)),
    n_spot = as.integer(table(o$sce$spatial.cluster))))
  spot_tab <- rbind(spot_tab, data.frame(
    barcode = colnames(o$sce), section = sec, group = sheet$group[i],
    x = colData(o$sce)$pxl_col_in_fullres, y = colData(o$sce)$pxl_row_in_fullres,
    cluster = o$sce$spatial.cluster))
}

write.csv(clu_tab, "09_cluster_count.csv", row.names = FALSE)
write.csv(spot_tab, "09_spot_cluster.csv", row.names = FALSE)

clu_cols <- c("#2166AC",
              "#C0392B",
              "#4D9221",
              "#7B3294",
              "#F4A582",
              "#00A087",
              "#E0C097",
              "#3C5488")
spot_tab$section <- factor(spot_tab$section, levels = sheet$section)
spot_tab$cluster <- factor(spot_tab$cluster)
spot_tab$facet <- factor(paste0(spot_tab$section, "\n", spot_tab$group),
                         levels = paste0(sheet$section, "\n", sheet$group))

p_clu <- ggplot(spot_tab, aes(x, y, colour = cluster)) +
  geom_point(size = 0.4) +
  scale_colour_manual(values = clu_cols, name = "BayesSpace\ncluster") +
  facet_wrap(~facet, scales = "free", nrow = 2) +
  scale_y_reverse() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(size = 7, face = "bold", lineheight = 0.9, margin = margin(b = 2)),

        legend.key.size = unit(0.28, "cm"), panel.spacing = unit(1.5, "pt"))
ggsave("09_cluster_spatial.pdf", p_clu, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("09_cluster_spatial.png", p_clu, width = 15.5, height = 12, units = "cm", dpi = 600)

seurat_clu <- read.csv(file.path(project_dir, "02cluster", "02_spot_cluster.csv"))
merged <- merge(spot_tab[, c("barcode", "section", "cluster")],
                seurat_clu[, c("barcode", "section", "cluster")],
                by = c("barcode", "section"), suffixes = c("_bs", "_seurat"))
ari_by_section <- sapply(split(merged, merged$section), function(sub) {
  mclust::adjustedRandIndex(sub$cluster_bs, sub$cluster_seurat)
})
ari_median <- median(ari_by_section)
ari_min <- min(ari_by_section)
ari_max <- max(ari_by_section)

ari_tab <- data.frame(
  section = names(ari_by_section),
  group = sheet$group[match(names(ari_by_section), sheet$section)],
  n_spot = as.integer(table(merged$section)),
  n_cluster_bayesspace = as.integer(sapply(split(merged$cluster_bs, merged$section), function(x) length(unique(x)))),
  n_cluster_seurat = as.integer(sapply(split(merged$cluster_seurat, merged$section), function(x) length(unique(x)))),
  ARI = round(as.numeric(ari_by_section), 3))
write.csv(ari_tab, "09_cluster_ari.csv", row.names = FALSE)

