library(Seurat)
library(ggplot2)
library(patchwork)
library(grid)
library(scales)

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else if (file.exists(file.path("13spatial", "06plate", "slice_plot.R"))) {
  setwd(file.path("13spatial", "06plate"))
}
project_dir <- dirname(getwd())
source(file.path(dirname(getwd()), "06plate", "slice_plot.R"))

sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
sec_lab <- sprintf("%s\n%s", sheet$section, sheet$group)
names(sec_lab) <- sheet$section
panel_cm <- 4.3
fig_w <- 3 * panel_cm + 3.2
fig_h <- 3 * panel_cm + 0.5

draw_slices <- function(dat, fig_name, col_name, pal, legend_name, col_lim = NULL,
                        point_size = 0.4, legend_ncol = 1) {
  for (s in sheet$section) {
    d <- dat[dat$section == s, ]
    d$section <- factor(d$section)
    save_slice(d, s, fig_name, col_name = col_name, pal = pal, legend_name = legend_name,
               col_lim = col_lim, point_size = point_size, legend_ncol = legend_ncol, label = "2 mm")
  }
  vals <- if (is.null(col_lim)) levels(dat[[col_name]]) else col_lim
  col_cm <- legend_col_cm(if (is.null(col_lim)) vals else sprintf("%g", col_lim),
                          title = legend_name, bar = !is.null(col_lim))
  g <- combine_slices(dat, sheet$section, sec_lab, "2 mm", col_name = col_name, pal = pal,
                      legend_name = legend_name, col_lim = col_lim, point_size = point_size,
                      legend_ncol = legend_ncol, panel_cm = panel_cm, col_cm = col_cm)
  save_gtable(g, fig_name, width = 3 * panel_cm + col_cm + 1.0, height = fig_h)
  saveRDS(g, sprintf("%s_obj.rds", fig_name))
  g
}

seu_list <- readRDS(file.path(project_dir, "01load_qc", "01_seu_list_qc.rds"))
qc <- data.frame()
for (seu in seu_list) {
  coords <- GetTissueCoordinates(seu)
  qc <- rbind(qc, data.frame(section = seu$section[1], group = seu$group[1],
                             x = coords$x, y = coords$y, nFeature = seu$nFeature_Spatial))
}
qc$section <- factor(qc$section, levels = sheet$section)
qc <- slice_um(qc)

size_tab <- do.call(rbind, lapply(sheet$section, function(s) {
  d <- qc[qc$section == s, ]
  data.frame(section = s, group = sheet$group[sheet$section == s], n_spot = nrow(d),
             width_um = round(diff(range(d$x))), height_um = round(diff(range(d$y))),
             median_nFeature = median(d$nFeature))
}))
write.csv(size_tab, "07_slice_qc_size.csv", row.names = FALSE)

draw_slices(qc, "07_slice_qc", "nFeature", grad_cols, "nFeature", col_lim = c(0, 3500))

clu <- read.csv(file.path(project_dir, "02cluster", "02_spot_cluster.csv"))
clu$section <- factor(clu$section, levels = sheet$section)
clu$cluster <- factor(clu$cluster)
clu <- slice_um(clu)
n_clu <- sapply(sheet$section, function(s) length(unique(clu$cluster[clu$section == s])))

draw_slices(clu, "08_slice_cluster", "cluster", cluster_cols, "Cluster")

lab <- read.csv(file.path(project_dir, "03transfer", "03_spot_celltype.csv"))
lab$section <- factor(lab$section, levels = sheet$section)
lab$celltype <- factor(lab$celltype, levels = names(cell_cols))
lab <- slice_um(lab)

draw_slices(lab, "09_slice_label", "celltype", cell_cols, "Cell type")

key <- read.csv(file.path(project_dir, "04keygene", "04_keygene_expr_spot.csv"))
names(key)[names(key) == "sample"] <- "section"
key$section <- factor(key$section, levels = sheet$section)
key$gene <- factor(key$gene, levels = unique(key$gene))
key <- slice_um(key)

for (g in levels(key$gene)) {
  for (s in sheet$section) {
    d <- key[key$gene == g & key$section == s, ]
    d$section <- factor(d$section)
    save_slice(d, s, sprintf("10_slice_keygene_%s", g), col_name = "expr",
               pal = grad_cols, legend_name = "Expression", col_lim = c(0, 3),
               point_size = 0.3, label = "2 mm")
  }
}
for (g in levels(key$gene)) {
  col_cm <- legend_col_cm(sprintf("%g", c(0, 3)), title = "Expression", bar = TRUE)
  p <- combine_slices(key[key$gene == g, ], sheet$section, sec_lab, "2 mm", col_name = "expr",
                      pal = grad_cols, legend_name = "Expression", col_lim = c(0, 3),
                      point_size = 0.3, panel_cm = panel_cm, col_cm = col_cm)
  save_gtable(p, sprintf("10_slice_keygene_%s", g),
              width = 3 * panel_cm + col_cm + 1.0, height = fig_h)
  saveRDS(p, sprintf("10_p_keygene_%s.rds", g))
}

prop <- read.csv(file.path(project_dir, "05deconv", "06_spotlight_prop.csv"))
spat <- read.csv(file.path(project_dir, "02cluster", "02_spot_cluster.csv"))[,
                 c("barcode", "section", "x", "y")]
dom <- merge(prop[, c("sample", "barcode", "dominant")], spat,
             by.x = c("sample", "barcode"), by.y = c("section", "barcode"))
names(dom)[1] <- "section"
dom$section <- factor(dom$section, levels = sheet$section)
dom$dominant <- factor(dom$dominant, levels = names(cell_cols))
dom <- slice_um(dom)

draw_slices(dom, "11_slice_deconv", "dominant", cell_cols, "Dominant\ncell type")

page <- read.csv(file.path(project_dir, "05deconv", "08_page_score_spot.csv"))
names(page) <- gsub("\\.", "/", names(page))
page$section <- factor(page$section, levels = sheet$section)
page_ct <- intersect(names(cell_cols), names(page))
best_ct <- page_ct[which.max(colMeans(page[, page_ct]))]
page <- slice_um(page)

draw_slices(page, "12_slice_page", best_ct, grad_cols, "PAGE score",
            col_lim = range(page[[best_ct]]))

bs <- read.csv(file.path(project_dir, "05deconv", "09_spot_cluster.csv"))
bs$section <- factor(bs$section, levels = sheet$section)
bs$cluster <- factor(bs$cluster)
bs <- slice_um(bs)

draw_slices(bs, "13_slice_bayesspace", "cluster", bayesspace_cols, "BayesSpace\ncluster")
