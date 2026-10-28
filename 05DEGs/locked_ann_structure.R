rm(list = ls()); gc()
library(neuralnet)
library(ggplot2)
library(scales)

mod_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
proj_dir <- dirname(mod_dir)
out_stem <- Sys.getenv("ANN_OUT", unset = "ANN_Network_Structure")
n_hidden <- 5
seed <- 798
top_n <- 12
col_pos <- "#C0392B"
col_neg <- "#2166AC"
col_node <- "#4D4D4D"
fs_gene <- 7
fs_hid <- 5.5
fs_leg <- 7
pt <- 72.27 / 25.4
fig_w <- 4.3
fig_h <- 3.3

hub_genes <- read.csv(file.path(mod_dir, "intersection_gene.csv"), stringsAsFactors = FALSE)
key_genes <- hub_genes$symbol

train_data <- read.csv(file.path(proj_dir, "00bulkdata", "dat.GSE87466.csv"), row.names = 1, check.names = FALSE)
train_group <- read.csv(file.path(proj_dir, "00bulkdata", "group.GSE87466.csv"), stringsAsFactors = FALSE)
colnames(train_group) <- c("sample", "group")

train_df <- as.data.frame(t(train_data[key_genes, ]))
train_df$label <- ifelse(train_group$group == "UC", 1, 0)
formula_obj <- as.formula(paste("label ~", paste(key_genes, collapse = " + ")))

tr_mean <- colMeans(train_df[, key_genes])
tr_sd <- apply(train_df[, key_genes], 2, sd)
train_df_scaled <- train_df
train_df_scaled[, key_genes] <- scale(train_df[, key_genes], center = tr_mean, scale = tr_sd)

set.seed(seed)
ann_model <- neuralnet(formula_obj, data = train_df_scaled, hidden = n_hidden, linear.output = FALSE)

w1 <- ann_model$weights[[1]][[1]]
w2 <- ann_model$weights[[1]][[2]]

n_in <- length(key_genes)
in_y <- seq(4, 1, length.out = n_in)
hid_y <- seq(4.5, 1.5, length.out = n_hidden)
out_y <- mean(range(in_y))

nodes <- rbind(
  data.frame(layer = "Input", name = key_genes, x = 0, y = in_y, size = 3.0),
  data.frame(layer = "Hidden", name = paste0("H", seq_len(n_hidden)), x = 1, y = hid_y, size = 4.2),
  data.frame(layer = "Output", name = "Output", x = 2, y = out_y, size = 3.6))

e1 <- data.frame(x = 0, y = rep(in_y, times = n_hidden), xend = 1, yend = rep(hid_y, each = n_in),
                 w = as.vector(w1[-1, , drop = FALSE]))
e2 <- data.frame(x = 1, y = hid_y, xend = 2, yend = rep(out_y, times = n_hidden),
                 w = as.vector(w2[-1, , drop = FALSE]))
e1 <- head(e1[order(-abs(e1$w)), ], top_n)
e2 <- head(e2[order(-abs(e2$w)), ], top_n)
e1$sgn <- ifelse(e1$w > 0, "positive", "negative")
e2$sgn <- ifelse(e2$w > 0, "positive", "negative")

arrow_head <- arrow(length = unit(0.05, "inches"), type = "closed")
p <- ggplot() +
  geom_curve(data = e1, aes(x = x, y = y, xend = xend, yend = yend, linewidth = abs(w), colour = sgn),
             curvature = 0.12, arrow = arrow_head, alpha = 0.75) +
  geom_curve(data = e2, aes(x = x, y = y, xend = xend, yend = yend, linewidth = abs(w), colour = sgn),
             curvature = 0.12, arrow = arrow_head, alpha = 0.75) +
  geom_point(data = nodes, aes(x = x, y = y, size = size), shape = 21, fill = "white",
             colour = col_node, stroke = 0.4) +
  geom_text(data = subset(nodes, layer == "Hidden"), aes(x = x, y = y, label = name),
            size = fs_hid / pt, colour = "grey20", family = "sans") +
  geom_text(data = subset(nodes, layer == "Input"), aes(x = x - 0.13, y = y, label = name),
            hjust = 1, size = fs_gene / pt, family = "sans") +
  geom_text(data = subset(nodes, layer == "Output"), aes(x = x + 0.13, y = y, label = name),
            hjust = 0, size = fs_gene / pt, family = "sans") +
  geom_text(data = data.frame(x = c(0, 1, 2), y = 4.85, name = c("Input", "Hidden", "Output")),
            aes(x = x, y = y, label = name), size = fs_hid / pt, colour = "grey45", family = "sans") +
  scale_linewidth(range = c(0.2, 1.2), name = "|w|") +
  scale_colour_manual(values = c(positive = col_pos, negative = col_neg), name = "sign") +
  scale_size_identity() +
  coord_cartesian(xlim = c(-0.68, 2.55), ylim = c(0.6, 5.1)) +
  theme_void(base_size = 8, base_family = "sans") +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = fs_leg), legend.text = element_text(size = fs_leg),
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.3, "cm"),
        legend.spacing.x = unit(0.12, "cm"), legend.box.spacing = unit(3, "pt"),
        plot.margin = margin(4, 6, 2, 4)) +
  guides(linewidth = guide_legend(order = 1),
         colour = guide_legend(order = 2, override.aes = list(linewidth = 0.8, alpha = 1)))

ggsave(file.path(mod_dir, paste0(out_stem, ".pdf")), p, width = fig_w, height = fig_h)
ggsave(file.path(mod_dir, paste0(out_stem, ".png")), p, width = fig_w, height = fig_h, dpi = 600)
