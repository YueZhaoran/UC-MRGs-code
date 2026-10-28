library(Seurat)
library(CellChat)
library(ggplot2)
library(jsonlite)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- dirname(getwd())
set.seed(123)

n_workers <- 8
spot_size <- 65
db_search <- "Secreted Signaling"

cell_cols <- c("B/Plasma" = "#2166AC", "T cell" = "#C0392B", "Epithelial" = "#4D9221",
               "Myeloid" = "#7B3294", "Stromal" = "#F4A582", "NK/ILC" = "#00A087",
               "Unassigned" = "#BDBDBD")
celltype_lv <- names(cell_cols)

seu_list <- readRDS(file.path(project_dir, "02cluster", "02_seu_list_cluster.rds"))
prop <- readRDS(file.path(project_dir, "05deconv", "06_spotlight_prop.rds"))
sheet <- read.csv(file.path(project_dir, "_sample_sheet.csv"))

st_list <- list()
for (i in 1:length(seu_list)) {
  sec <- names(seu_list)[i]
  seu <- seu_list[[i]]
  pr <- prop[prop$sample == sec, ]
  pr <- pr[match(colnames(seu), pr$barcode), celltype_lv]
  seu$labels <- factor(celltype_lv[max.col(as.matrix(pr), ties.method = "first")],
                       levels = celltype_lv)
  st_list[[sec]] <- seu
}
saveRDS(st_list, "00_st_list.rds")
for (sec in names(st_list)) {
}

cellchat_spatial <- list()

for (sec in names(st_list)) {
  st <- st_list[[sec]]
  gsm <- sheet$gsm[match(sec, sheet$section)]

  data.input <- GetAssayData(st, assay = "Spatial", layer = "data")
  meta <- data.frame(labels = st$labels, row.names = colnames(data.input))
  coords <- as.matrix(GetTissueCoordinates(st)[colnames(data.input), 1:2])
  colnames(coords) <- c("x_cent", "y_cent")

  sf <- fromJSON(file.path(project_dir, "00raw", paste0(gsm, "_", sec),
                           "spatial", "scalefactors_json.json"))
  spatial.factors <- data.frame(ratio = spot_size / sf$spot_diameter_fullres,
                                tol = spot_size / 2)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = coords,
                             spatial.factors = spatial.factors)
  db_use <- subsetDB(CellChatDB.human, search = db_search)
  cellchat@DB <- db_use

  cellchat <- subsetData(cellchat)
  plan("multisession", workers = n_workers)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, contact.knn.k = 6, scale.distance = 2)
  plan("sequential")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  cellchat_spatial[[sec]] <- cellchat

  comm <- subsetCommunication(cellchat)
  comm$section <- sec
  comm$group <- sheet$group[match(sec, sheet$section)]
  write.csv(comm, paste0("01_comm_", sec, ".csv"), row.names = FALSE)
}

saveRDS(cellchat_spatial, "01_cellchat_spatial_list.rds")

comm_all <- data.frame()
for (sec in names(cellchat_spatial)) {
  tmp <- subsetCommunication(cellchat_spatial[[sec]])
  tmp$section <- sec
  tmp$group <- sheet$group[match(sec, sheet$section)]
  comm_all <- rbind(comm_all, tmp)
}
write.csv(comm_all, "01_comm_all.csv", row.names = FALSE)

for (sec in names(cellchat_spatial)) {
  st <- st_list[[sec]]
  cellchat <- cellchat_spatial[[sec]]

  coords <- as.matrix(GetTissueCoordinates(st)[colnames(st), 1:2])
  rownames(coords) <- colnames(st)
  colnames(coords) <- c("x_cent", "y_cent")
  meta_plot <- data.frame(labels = st$labels, samples = sec,
                          row.names = colnames(st))

  p_count <- netVisual_spatial(net = cellchat@net$count, coordinates = coords,
                               meta = meta_plot,
                               title.name = paste0("Spatial cell-cell communication (count) - ",
                                        sec),
                               edge.width.max = 6, vertex.weight = 20, vertex.label.cex = 1.5,
                               edge.curved = 0.2, alpha.edge = 0.7, weight.scale = TRUE)
  ggsave(paste0("01_spatial_count_", sec, ".pdf"), p_count, width = 15.5, height = 12,
         units = "cm", dpi = 600)
  ggsave(paste0("01_spatial_count_", sec, ".png"), p_count, width = 15.5, height = 12,
         units = "cm", dpi = 600)

  p_weight <- netVisual_spatial(net = cellchat@net$weight, coordinates = coords,
                                meta = meta_plot,
                                title.name = paste0("Spatial cell-cell communication (weight) - ",
                                        sec),
                                edge.width.max = 6, vertex.weight = 20, vertex.label.cex = 1.5,
                                edge.curved = 0.2, alpha.edge = 0.7, weight.scale = TRUE)
  ggsave(paste0("01_spatial_weight_", sec, ".pdf"), p_weight, width = 15.5, height = 12,
         units = "cm", dpi = 600)
  ggsave(paste0("01_spatial_weight_", sec, ".png"), p_weight, width = 15.5, height = 12,
         units = "cm", dpi = 600)
}

nodes <- data.frame()
edges <- data.frame()
net_sum <- data.frame()

for (sec in names(cellchat_spatial)) {
  st <- st_list[[sec]]
  coords <- as.matrix(GetTissueCoordinates(st)[colnames(st), 1:2])
  rownames(coords) <- colnames(st)
  labels <- as.character(st$labels)
  names(labels) <- colnames(st)

  count <- as.matrix(cellchat_spatial[[sec]]@net$count)
  weight <- as.matrix(cellchat_spatial[[sec]]@net$weight)

  for (ct in rownames(count)) {
    idx <- labels == ct
    nodes <- rbind(nodes, data.frame(
      section = sec, group = sheet$group[match(sec, sheet$section)], celltype = ct,
      n_spot = sum(idx), x_cent = median(coords[idx, 1]), y_cent = median(coords[idx, 2])))
  }

  nz <- which(count > 0 | weight > 0, arr.ind = TRUE)
  if (nrow(nz) > 0) {
    edges <- rbind(edges, data.frame(
      section = sec, group = sheet$group[match(sec, sheet$section)],
      source = rownames(count)[nz[, 1]], target = colnames(count)[nz[, 2]],
      count = count[nz], weight = weight[nz]))
  }

  net_sum <- rbind(net_sum, data.frame(
    section = sec, group = sheet$group[match(sec, sheet$section)], n_spot = ncol(st),
    n_celltype = nrow(count), total_count = sum(count), total_weight = sum(weight)))
}

write.csv(nodes, "01_network_nodes.csv", row.names = FALSE)
write.csv(edges, "01_network_edges.csv", row.names = FALSE)
write.csv(net_sum, "01_network_summary.csv", row.names = FALSE)

grp_sum <- data.frame()
for (g in unique(sheet$group)) {
  secs <- sheet$section[sheet$group == g]
  net_count <- Reduce("+", lapply(cellchat_spatial[secs], function(x) x@net$count))
  net_weight <- Reduce("+", lapply(cellchat_spatial[secs], function(x) x@net$weight))
  write.csv(as.data.frame(net_count), paste0("01_net_count_", g, ".csv"))
  write.csv(as.data.frame(net_weight), paste0("01_net_weight_", g, ".csv"))
  grp_sum <- rbind(grp_sum, data.frame(
    group = g, n_section = length(secs), n_pair = sum(net_count > 0),
    total_count = sum(net_count), total_weight = signif(sum(net_weight), 6)))
}
write.csv(grp_sum, "01_group_summary.csv", row.names = FALSE)
