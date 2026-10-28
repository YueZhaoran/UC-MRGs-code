library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())

key_cell <- "Stromal"
groups <- c("UC", "Normal")
top_n <- 20
top_path <- 15

cell_cols <- c("B/Plasma" = "#2166AC", "T cell" = "#C0392B", "Epithelial" = "#4D9221",
               "Myeloid" = "#7B3294", "Stromal" = "#F4A582", "NK/ILC" = "#00A087",
               "Unassigned" = "#BDBDBD")
celltype_lv <- names(cell_cols)

comm <- read.csv("01_comm_all.csv", stringsAsFactors = FALSE)
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))
sec_lv <- sheet$section

read_net <- function(f) {
  m <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  m[, celltype_lv, drop = FALSE][celltype_lv, , drop = FALSE]
}
net_count <- lapply(groups, function(g) read_net(paste0("01_net_count_", g, ".csv")))
net_weight <- lapply(groups, function(g) read_net(paste0("01_net_weight_", g, ".csv")))
names(net_count) <- names(net_weight) <- groups

net_long <- data.frame()
for (g in groups) {
  for (m in c("count", "weight")) {
    mat <- if (m == "count") net_count[[g]] else net_weight[[g]]
    idx <- which(mat != 0, arr.ind = TRUE)
    net_long <- rbind(net_long, data.frame(
      group = g, measure = m, source = rownames(mat)[idx[, 1]],
      target = colnames(mat)[idx[, 2]], value = mat[idx]))
  }
}
net_long$source <- factor(net_long$source, levels = celltype_lv)
net_long$target <- factor(net_long$target, levels = rev(celltype_lv))

sec_tab <- comm %>%
  filter(source == key_cell | target == key_cell) %>%
  group_by(section, group) %>%
  summarise(n_pair = n(), n_out = sum(source == key_cell), n_in = sum(target == key_cell),
            .groups = "drop") %>%
  mutate(section = factor(section, levels = sec_lv))

heat_plot <- function(m, title, measure, lim) {
  idx <- which(m != 0, arr.ind = TRUE)
  d <- data.frame(source = factor(rownames(m)[idx[, 1]], levels = celltype_lv),
                  target = factor(colnames(m)[idx[, 2]], levels = rev(celltype_lv)),
                  value = m[idx])
  ggplot(d, aes(source, target, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = if (measure == "count") value else sprintf("%.1e", value)),
              size = 1.9, colour = "black") +
    scale_fill_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B"),
                         limits = c(0, lim), name = ifelse(measure == "count",
                                                           "Interactions", "Strength")) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 8, base_family = "sans") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.8, "cm"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6))
}

lim_count <- max(sapply(net_count, max))
lim_weight <- max(sapply(net_weight, max))

p_full <- heat_plot(net_count[["UC"]], "UC (7 sections)", "count", lim_count) +
  heat_plot(net_count[["Normal"]], "Normal (2 sections)", "count", lim_count) +
  plot_layout(ncol = 2)
ggsave("01_net_count_heatmap.pdf", p_full, width = 15.5, height = 8, units = "cm", dpi = 600)
ggsave("01_net_count_heatmap.png", p_full, width = 15.5, height = 8, units = "cm", dpi = 600)

p_full_w <- heat_plot(net_weight[["UC"]], "UC (7 sections)", "weight", lim_weight) +
  heat_plot(net_weight[["Normal"]], "Normal (2 sections)", "weight", lim_weight) +
  plot_layout(ncol = 2)
ggsave("02_net_weight_heatmap.pdf", p_full_w, width = 15.5, height = 8, units = "cm", dpi = 600)
ggsave("02_net_weight_heatmap.png", p_full_w, width = 15.5, height = 8, units = "cm", dpi = 600)

p_key <- net_long %>%
  filter(measure == "weight", source == key_cell | target == key_cell) %>%
  mutate(pair = ifelse(source == key_cell,
                       paste0(key_cell, " -> ", target), paste0(source, " -> ", key_cell)),
         pair = factor(pair, levels = unique(pair[order(-value)]))) %>%
  ggplot(aes(pair, value, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_manual(values = c(UC = "#C0392B", Normal = "#2166AC"), name = "Group") +
  labs(x = NULL, y = "Communication strength") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        legend.position = "top", legend.key.size = unit(0.3, "cm"))
ggsave("03_keycell_pair.pdf", p_key, width = 12, height = 8, units = "cm", dpi = 600)
ggsave("03_keycell_pair.png", p_key, width = 12, height = 8, units = "cm", dpi = 600)

p_sec <- ggplot(sec_tab, aes(section, n_pair, fill = group)) +
  geom_col(width = 0.68) +
  scale_fill_manual(values = c(UC = "#C0392B", Normal = "#2166AC"), name = "Group") +
  labs(x = NULL, y = paste0("Interactions involving ", key_cell)) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "top", legend.key.size = unit(0.3, "cm"))
ggsave("04_keycell_section.pdf", p_sec, width = 9, height = 7, units = "cm", dpi = 600)
ggsave("04_keycell_section.png", p_sec, width = 9, height = 7, units = "cm", dpi = 600)

plate_10c <- (p_full / (p_key | p_sec)) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("05_plate_10c.pdf", plate_10c, width = 15.5, height = 18, units = "cm", dpi = 600)
ggsave("05_plate_10c.png", plate_10c, width = 15.5, height = 18, units = "cm", dpi = 600)

key_dir <- net_long %>% filter(measure == "weight", source == key_cell)
lr_sum <- comm %>%
  filter(source == key_cell) %>%
  group_by(interaction_name_2) %>%
  summarise(prob = sum(prob), .groups = "drop") %>%
  arrange(desc(prob)) %>%
  slice_head(n = top_n)
lr_use <- lr_sum$interaction_name_2

bubble_dat <- comm %>%
  filter(source == key_cell, interaction_name_2 %in% lr_use) %>%
  mutate(section = factor(section, levels = sec_lv),
         lr = factor(interaction_name_2, levels = rev(lr_use)))

bubble_dat <- complete(bubble_dat, section, lr, fill = list(prob = 0, pval = 1)) %>%
  mutate(section = factor(section, levels = sec_lv),
         lr = factor(lr, levels = rev(lr_use)))

p_bubble <- ggplot(bubble_dat, aes(section, lr)) +
  geom_point(aes(size = -log10(pval + 1e-4), colour = prob)) +
  scale_colour_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                         limits = c(0, max(bubble_dat$prob)),
                         name = "Communication\nprobability") +
  scale_size_continuous(range = c(0.6, 3.2), name = "-log10(P)") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"), legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))
ggsave("06_keycell_lr_bubble.pdf", p_bubble, width = 12, height = 10, units = "cm", dpi = 600)
ggsave("06_keycell_lr_bubble.png", p_bubble, width = 12, height = 10, units = "cm", dpi = 600)

path_dat <- comm %>%
  filter(source == key_cell) %>%
  group_by(pathway_name, group) %>%
  summarise(prob = sum(prob), n = n(), .groups = "drop") %>%
  group_by(pathway_name) %>%
  mutate(total = sum(prob)) %>%
  ungroup() %>%
  filter(pathway_name %in% (arrange(., desc(total)) %>% distinct(pathway_name) %>% head(top_path) %>% pull(pathway_name)))
path_dat$pathway_name <- factor(path_dat$pathway_name,
                                levels = rev(unique(path_dat$pathway_name[order(path_dat$total)])))

p_path <- ggplot(path_dat, aes(pathway_name, prob, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  coord_flip() +
  scale_fill_manual(values = c(UC = "#C0392B", Normal = "#2166AC"), name = "Group") +
  labs(x = NULL, y = "Total communication probability") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 6.5),
        legend.position = "top", legend.key.size = unit(0.3, "cm"))
ggsave("07_keycell_pathway.pdf", p_path, width = 10, height = 10, units = "cm", dpi = 600)
ggsave("07_keycell_pathway.png", p_path, width = 10, height = 10, units = "cm", dpi = 600)

plate_10d <- (p_bubble | p_path) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11, family = "sans"))
ggsave("08_plate_10d.pdf", plate_10d, width = 15.5, height = 12, units = "cm", dpi = 600)
ggsave("08_plate_10d.png", plate_10d, width = 15.5, height = 12, units = "cm", dpi = 600)

n_uc <- sum(sheet$group == "UC")
n_nc <- sum(sheet$group == "Normal")
n_spot <- 11444
n_int_nc <- sort(sapply(groups, function(g) sum(net_count[[g]])))

