library(Seurat)
library(CellChat)
library(ggplot2)
library(ggrepel)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

cell_cols <- c("B/Plasma" = "#2166AC", "T cell" = "#C0392B", "Epithelial" = "#4D9221",
               "Myeloid" = "#7B3294", "Stromal" = "#F4A582", "NK/ILC" = "#00A087",
               "Unassigned" = "#BDBDBD")

st_list <- readRDS("00_st_list.rds")
cellchat_spatial <- readRDS("01_cellchat_spatial_list.rds")
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))

snap_node_coords <- function(coords, labels, celltypes) {
  node_coords <- matrix(NA, nrow = length(celltypes), ncol = 2,
                        dimnames = list(celltypes, c("x_raw", "y_raw")))
  for (ct in celltypes) {
    idx <- which(labels == ct)
    if (length(idx) == 0) next
    med_x <- median(coords[idx, 1])
    med_y <- median(coords[idx, 2])
    dist2 <- (coords[idx, 1] - med_x)^2 + (coords[idx, 2] - med_y)^2
    best <- idx[which.min(dist2)]
    node_coords[ct, ] <- as.numeric(coords[best, 1:2])
  }
  node_coords
}

plot_snap <- function(net, coords_raw, labels, title) {
  celltypes <- rownames(net)
  node_coords_raw <- snap_node_coords(coords_raw, labels, celltypes)

  spot_plot <- data.frame(x_cent = coords_raw[, 2], y_cent = coords_raw[, 1],
                          labels = labels)
  node_plot <- data.frame(X1 = node_coords_raw[, 2], X2 = node_coords_raw[, 1],
                          node_idents = factor(celltypes, levels = celltypes))
  node_plot <- node_plot[complete.cases(node_plot), , drop = FALSE]

  diag(net) <- 0
  nonzero <- which(net > 0, arr.ind = TRUE)
  edges <- NULL
  if (nrow(nonzero) > 0) {
    edges <- data.frame(
      X1 = node_coords_raw[rownames(net)[nonzero[, 1]], 2],
      Y1 = node_coords_raw[rownames(net)[nonzero[, 1]], 1],
      X2 = node_coords_raw[colnames(net)[nonzero[, 2]], 2],
      Y2 = node_coords_raw[colnames(net)[nonzero[, 2]], 1],
      source = rownames(net)[nonzero[, 1]],
      width = as.numeric(net[nonzero]))
    edges <- edges[complete.cases(edges), , drop = FALSE]
  }

  p <- ggplot(spot_plot, aes(x_cent, y_cent)) +
    geom_point(aes(colour = labels), alpha = 0.15, size = 1.5, show.legend = FALSE)

  if (!is.null(edges) && nrow(edges) > 0) {
    p <- p + geom_curve(data = edges,
                        aes(x = X1, y = Y1, xend = X2, yend = Y2,
                            colour = source, size = width),
                        curvature = 0.2, alpha = 0.7,
                        arrow = arrow(angle = 5, type = "closed", length = unit(0.1, "inches")),
                        show.legend = FALSE) +
      scale_size_continuous(range = c(0.3, 3), guide = "none")
  }
  p +
    geom_point(data = node_plot, aes(X1, X2, colour = node_idents), size = 2.2) +
    ggrepel::geom_text_repel(data = node_plot, aes(X1, X2, label = node_idents),
                             colour = "black", size = 2.2, segment.size = 0.2,
                             max.overlaps = Inf) +
    scale_colour_manual(values = cell_cols) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title) +
    theme_void(base_size = 8, base_family = "sans") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
          legend.position = "none")
}

for (sec in names(cellchat_spatial)) {
  st <- st_list[[sec]]
  coords_raw <- as.matrix(GetTissueCoordinates(st)[colnames(st), 1:2])
  rownames(coords_raw) <- colnames(st)
  colnames(coords_raw) <- c("x_raw", "y_raw")
  labels <- as.character(st$labels)
  names(labels) <- colnames(st)

  net_count <- as.matrix(cellchat_spatial[[sec]]@net$count)
  net_weight <- as.matrix(cellchat_spatial[[sec]]@net$weight)

  p <- plot_snap(net_count, coords_raw, labels,
                 paste0("Spatial cell-cell communication (count) - ", sec))
  ggsave(paste0("02_snap_count_", sec, ".pdf"), p, width = 15.5, height = 12,
         units = "cm", dpi = 600)
  ggsave(paste0("02_snap_count_", sec, ".png"), p, width = 15.5, height = 12,
         units = "cm", dpi = 600)

  p <- plot_snap(net_weight, coords_raw, labels,
                 paste0("Spatial cell-cell communication (weight) - ", sec))
  ggsave(paste0("02_snap_weight_", sec, ".pdf"), p, width = 15.5, height = 12,
         units = "cm", dpi = 600)
  ggsave(paste0("02_snap_weight_", sec, ".png"), p, width = 15.5, height = 12,
         units = "cm", dpi = 600)

}
