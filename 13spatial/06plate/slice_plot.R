library(grid)
library(gtable)

px_per_um <- c(B5 = 2.25, B4 = 2.25, B8 = 1.106, B9 = 1.106, B12 = 2.26,
               B13 = 2.26, C2 = 1.92, B10 = 2.26, C5 = 1.92)
span_um <- 7600
bar_um <- 2000

cell_cols <- c("B/Plasma" = "#2166AC", "T cell" = "#C0392B", "Epithelial" = "#4D9221",
               "Myeloid" = "#7B3294", "Stromal" = "#F4A582", "NK/ILC" = "#00A087",
               "Unassigned" = "#BDBDBD")
cluster_cols <- c("#2166AC", "#C0392B", "#4D9221", "#7B3294", "#F4A582",
                  "#00A087", "#E0C097", "#3C5488", "#B2182B", "#8AB17D",
                  "#E76F51", "#BDBDBD")
bayesspace_cols <- c("#2166AC", "#C0392B", "#4D9221", "#7B3294", "#F4A582",
                     "#00A087", "#E0C097", "#3C5488")
grad_cols <- c("#2166AC", "#F7F7F7", "#B2182B")

slice_um <- function(dat) {
  dat$x <- dat$x / px_per_um[as.character(dat$section)]
  dat$y <- dat$y / px_per_um[as.character(dat$section)]
  dat
}

legend_size <- function(base_size) {
  theme(legend.title = element_text(size = base_size * 1.05),
        legend.text = element_text(size = base_size * 0.95),
        legend.key.size = unit(base_size * 0.05, "cm"),
        legend.spacing.y = unit(0.03, "cm"),
        legend.margin = margin(0, 0, 0, 0))
}

slice_plot <- function(dat, col_name, pal, legend_name, title = NULL, point_size = 0.4,
                       col_lim = NULL, legend_ncol = 1, label = NULL, base_size = 8) {
  cx <- mean(range(dat$x, na.rm = TRUE))
  cy <- mean(range(dat$y, na.rm = TRUE))
  half <- max(span_um, diff(range(dat$x)), diff(range(dat$y)), na.rm = TRUE) / 2 * 1.04

  p <- ggplot(dat, aes(x, y, colour = .data[[col_name]])) +
    geom_point(size = point_size) +
    scale_x_continuous(limits = cx + c(-half, half)) +
    scale_y_continuous(limits = cy + c(-half, half)) +
    scale_y_reverse() +
    coord_fixed(1) +
    labs(x = NULL, y = NULL, title = title) +
    theme_void(base_size = base_size, base_family = "sans") +
    legend_size(base_size) +
    theme(plot.title = element_text(size = base_size, face = "bold", hjust = 0.5, margin = margin(b = 2)),

          plot.margin = margin(1, 1, 1, 1), panel.spacing = unit(1.5, "pt"))

  if (is.null(col_lim)) {
    p <- p + scale_colour_manual(values = pal, name = legend_name, drop = FALSE,
                                 guide = guide_legend(override.aes = list(size = 2),
                                        ncol = legend_ncol))
  } else {
    p <- p + scale_colour_gradientn(colours = pal, limits = col_lim, oob = scales::squish,
                                    name = legend_name,
                                    guide = guide_colourbar(barwidth = unit(base_size * 0.035, "cm"),

                                                            barheight = unit(base_size * 0.22,
                                        "cm")))
  }
  bar_dat <- data.frame(x = cx - half * 0.92, xend = cx - half * 0.92 + bar_um,
                        y = cy + half * 0.88, lab = label)
  p + geom_segment(data = bar_dat, aes(x = x, xend = xend, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.45, colour = "black") +
    geom_text(data = bar_dat, aes(x = (x + xend) / 2, y = y, label = lab),
              inherit.aes = FALSE, vjust = -0.35, size = base_size * 0.75 / .pt, colour = "black")
}

legend_col_cm <- function(labels, title = NULL, base_size = 8, bar = FALSE) {
  txt <- max(nchar(c(labels, title))) * base_size * 0.48 / 72 * 2.54
  min(3.2, (if (bar) 0.30 else 0.75) + txt)
}

wrap_label <- function(x, width) {
  if (nchar(x) <= width) return(x)
  paste(strwrap(x, width = width), collapse = "\n")
}

legend_grob <- function(vals, cols, legend_name, col_lim = NULL, base_size = 8,
                        col_cm = 2.6, height_cm = 12.4) {
  nx <- function(v) unit(v, "native")
  vp <- viewport(width = unit(col_cm, "cm"), height = unit(height_cm, "cm"),
                 xscale = c(0, col_cm), yscale = c(0, height_cm))
  kids <- gList(textGrob(legend_name, x = nx(0.06), y = nx(height_cm * 0.975), hjust = 0, vjust = 1,
                         gp = gpar(fontsize = base_size * 1.05, fontfamily = "sans")))
  if (is.null(col_lim)) {
    wrap_w <- max(6, floor((col_cm - 0.55) / (base_size * 0.5 / 72 * 2.54)))
    vals <- vapply(vals, wrap_label, character(1), width = wrap_w)
    n <- length(vals)
    pos <- seq(height_cm * 0.90, height_cm * 0.03, length.out = n + 1)
    y <- (pos[-1] + pos[-(n + 1)]) / 2
    kids <- gList(kids,
                  rectGrob(x = nx(0.10), y = nx(y), width = nx(0.26), height = nx(0.26),
                           just = c("left", "centre"), gp = gpar(fill = cols, col = NA)),
                  textGrob(vals, x = nx(0.44), y = nx(y), hjust = 0, vjust = 0.5,
                           gp = gpar(fontsize = base_size * 0.95,
                                     fontfamily = "sans",
                                     lineheight = 0.95)))
  } else {
    kids <- gList(kids,
                  rasterGrob(as.raster(matrix(rev(cols), ncol = 1)), x = nx(0.10),
                             y = nx(height_cm * 0.5), width = nx(0.26), height = nx(height_cm * 0.62),
                             interpolate = TRUE),
                  textGrob(sprintf("%g", col_lim[2]), x = nx(0.44), y = nx(height_cm * 0.84),
                           hjust = 0, vjust = 0.5,
                           gp = gpar(fontsize = base_size * 0.95, fontfamily = "sans")),
                  textGrob(sprintf("%g", col_lim[1]), x = nx(0.44), y = nx(height_cm * 0.16),
                           hjust = 0, vjust = 0.5,
                           gp = gpar(fontsize = base_size * 0.95, fontfamily = "sans")))
  }
  gTree(children = kids, vp = vp,
        widthDetails = unit(col_cm, "cm"), heightDetails = unit(height_cm, "cm"))
}

per_slice_dir <- file.path(dirname(getwd()), "06plate", "per_slice")

save_slice <- function(dat, section, fig_name, label, width = 5.5, height = 5.5, ...) {
  dir.create(per_slice_dir, showWarnings = FALSE)
  p <- slice_plot(dat[dat$section == section, ], label = label, ...)
  ggsave(file.path(per_slice_dir, sprintf("%s_%s.pdf", fig_name, section)),
         p, width = width, height = height, units = "cm", dpi = 600)
  ggsave(file.path(per_slice_dir, sprintf("%s_%s.png", fig_name, section)),
         p, width = width, height = height, units = "cm", dpi = 600)
}

save_gtable <- function(g, fig_name, width, height) {
  pdf(sprintf("%s.pdf", fig_name), width = width / 2.54, height = height / 2.54)
  grid.newpage(); grid.draw(g); dev.off()
  png(sprintf("%s.png", fig_name), width = width, height = height, units = "cm", res = 600)
  grid.newpage(); grid.draw(g); dev.off()
}

slice_gtable <- function(panels, leg, panel_cm = 4.3, col_cm = 2.6) {
  gt <- gtable(widths = unit.c(rep(unit(panel_cm, "cm"), 3), unit(col_cm, "cm")),
               heights = unit.c(rep(unit(panel_cm, "cm"), 3), unit(1, "null")))
  k <- 1
  for (r in 1:3) for (c in 1:3) {
    gt <- gtable_add_grob(gt, ggplotGrob(panels[[k]]), t = r, l = c)
    k <- k + 1
  }
  gtable_add_grob(gt, leg, t = 1, b = 3, l = 4)
}

combine_slices <- function(dat, sections, labels, label, col_name, pal, legend_name,
                           col_lim = NULL, legend_ncol = 1, point_size = 0.4, base_size = 8,
                           panel_cm = 4.3, col_cm = 2.6) {
  panels <- lapply(sections, function(s) {
    slice_plot(dat[dat$section == s, ], col_name, pal, legend_name,
               title = labels[match(s, sections)], point_size = point_size, col_lim = col_lim,
               legend_ncol = legend_ncol, label = label, base_size = base_size) +
      theme(legend.position = "none")
  })
  vals <- if (is.null(col_lim)) levels(dat[[col_name]]) else NULL
  cols <- if (is.null(col_lim)) unname(pal[seq_along(vals)]) else pal
  leg <- legend_grob(vals, cols, legend_name, col_lim = col_lim, base_size = base_size,
                     col_cm = col_cm, height_cm = 3 * panel_cm)
  slice_gtable(panels, leg, panel_cm = panel_cm, col_cm = col_cm)
}
