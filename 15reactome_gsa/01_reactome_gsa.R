rm(list = ls()); gc()

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
project_dir <- file.path(getwd(), "..")

library(Seurat)
library(ReactomeGSA)
library(pheatmap)


target_celltype <- "Stromal"
top_n           <- 20

scData <- readRDS(file.path(project_dir, "01singleCell", "05_scData_annotated.rds"))

Idents(scData) <- "external_annotation"
gsva_result <- analyse_sc_clusters(scData, verbose = FALSE)

significant_pathways <- pathways(gsva_result, p = 0.05)

rownames(significant_pathways) <- significant_pathways$Name
significant_pathways <- significant_pathways[, -1, drop = FALSE]

pathway_matrix <- as.matrix(significant_pathways)
storage.mode(pathway_matrix) <- "numeric"
write.csv(pathway_matrix, "01_reactome_gsa_pathways.csv")

colnames(pathway_matrix) <- sub("\\.Seurat$", "", colnames(pathway_matrix))
colnames(pathway_matrix) <- gsub("_", " ", colnames(pathway_matrix))
colnames(pathway_matrix)[colnames(pathway_matrix) == "BPlasma"] <- "B/Plasma"
colnames(pathway_matrix)[colnames(pathway_matrix) == "NKILC"] <- "NK/ILC"

other <- setdiff(colnames(pathway_matrix), target_celltype)
diff_tab <- data.frame(
  pathway = rownames(pathway_matrix),
  mean_keycell = pathway_matrix[, target_celltype],
  mean_other = rowMeans(pathway_matrix[, other, drop = FALSE]))
diff_tab$delta <- diff_tab$mean_keycell - diff_tab$mean_other
diff_tab <- diff_tab[order(-diff_tab$delta), ]
write.csv(diff_tab, "01_keycell_vs_other_pathways.csv", row.names = FALSE)


top_per_cell <- unique(unlist(apply(pathway_matrix, 2, function(x) {
  names(sort(x, decreasing = TRUE)[1:top_n])
})))
pathway_matrix_top <- pathway_matrix[rownames(pathway_matrix) %in% top_per_cell, ]

pdf("02_reactome_gsa_heatmap.pdf", width = 12, height = 10)
pheatmap(pathway_matrix_top,
         main = paste0("Top ", top_n, " pathways per cell type"),
         fontsize_row = 6, fontsize_col = 10,
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "row")
dev.off()

png("02_reactome_gsa_heatmap.png", width = 12, height = 10, units = "in", res = 300)
pheatmap(pathway_matrix_top,
         main = paste0("Top ", top_n, " pathways per cell type"),
         fontsize_row = 6, fontsize_col = 10,
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "row")
dev.off()
