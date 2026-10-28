rm(list = ls()); gc()

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
project_dir <- file.path(getwd(), "..")
set.seed(123)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggpubr)
library(tidyr)
library(viridis)
library(stringr)
library(AUCell)
library(UCell)
library(singscore)
library(GSVA)
library(GSEABase)
library(RColorBrewer)
library(BiocParallel)
library(scales)

dir_csv   <- "01_CSV"
dir_plots <- "02_Plots"
dir_data  <- "03_Data"
block_size <- 1000

scRNA <- readRDS(file.path(project_dir, "01singleCell", "05_scData_annotated.rds"))
DefaultAssay(scRNA) <- "RNA"
if ("JoinLayers" %in% ls("package:Seurat")) scRNA <- JoinLayers(scRNA, assay = "RNA")

key_genes <- unique(read.csv(file.path(project_dir, "05DEGs", "hub_gene.csv"))$symbol)
genes_found <- key_genes[key_genes %in% rownames(scRNA)]
geneset <- list(KeyGenes = genes_found)

matrix_sparse <- GetAssayData(scRNA, assay = "RNA", layer = "data")
cells_rankings <- AUCell_buildRankings(matrix_sparse, nCores = 1, plotStats = FALSE)
AUC <- AUCell_calcAUC(geneset, cells_rankings, nCores = 1, aucMaxRank = nrow(cells_rankings) * 0.1)
auc <- as.numeric(getAUC(AUC)["KeyGenes", ])
rm(matrix_sparse, cells_rankings, AUC); gc()

scRNA <- AddModuleScore_UCell(scRNA, features = geneset, name = "KeyGenes_UCell_", slot = "data")
ucell_score <- scRNA@meta.data[[grep("KeyGenes_UCell_", colnames(scRNA@meta.data), value = TRUE)]]

matrix_sparse <- GetAssayData(scRNA, assay = "RNA", layer = "data")
n_cells <- ncol(matrix_sparse)
n_blocks <- ceiling(n_cells / block_size)
singscore_scores <- numeric(n_cells)
for (i in 1:n_blocks) {
  idx <- ((i - 1) * block_size + 1):min(i * block_size, n_cells)
  expr_block <- as.matrix(matrix_sparse[, colnames(matrix_sparse)[idx]])
  ranked_expr <- rankGenes(expr_block)
  singscore_scores[idx] <- simpleScore(rankData = ranked_expr, upSet = genes_found,
                                       centerScore = TRUE)$TotalScore
  rm(expr_block, ranked_expr); gc()
}
rm(matrix_sparse); gc()

scRNA <- AddModuleScore(scRNA, features = geneset, name = "KeyGenes_Add_")
add_score <- scRNA@meta.data[[grep("KeyGenes_Add_", colnames(scRNA@meta.data), value = TRUE)]]

matrix_sparse <- GetAssayData(scRNA, assay = "RNA", layer = "data")
ssgsea_scores <- numeric(n_cells)
for (i in 1:n_blocks) {
  idx <- ((i - 1) * block_size + 1):min(i * block_size, n_cells)
  expr_block <- as.matrix(matrix_sparse[, colnames(matrix_sparse)[idx]])
  ssgsea_par <- ssgseaParam(exprData = expr_block,
                            geneSets = geneset,
                            alpha = 0.25,
                            normalize = FALSE)
  block_scores <- as.numeric(gsva(ssgsea_par, BPPARAM = SerialParam(), verbose = FALSE)[1, ])
  block_scores[is.nan(block_scores)] <- 0
  ssgsea_scores[idx] <- block_scores
  rm(expr_block, ssgsea_par, block_scores); gc()
}
rm(matrix_sparse); gc()

methods <- c("AUCell", "UCell", "singscore", "ssGSEA", "AddModuleScore")
score_df <- data.frame(AUCell = auc, UCell = ucell_score, singscore = singscore_scores,
                       ssGSEA = ssgsea_scores, AddModuleScore = add_score)
scRNA@meta.data <- cbind(scRNA@meta.data, score_df)

score_norm <- as.data.frame(apply(scale(scRNA@meta.data[, methods]), 2, scales::rescale))
scRNA$KeyGenes_Combined <- rowMeans(score_norm, na.rm = TRUE)
saveRDS(scRNA, file = file.path(dir_data, "scRNA_with_KeyGenes_scores.rds"))

existing_types <- levels(scRNA$external_annotation)
if (is.null(existing_types)) existing_types <- unique(scRNA$external_annotation)
celltype_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(existing_types))
names(celltype_colors) <- existing_types
unique_groups <- unique(scRNA$group)
group_colors <- c("#377EB8", "#E41A1C")[1:length(unique_groups)]
names(group_colors) <- unique_groups

for (method in methods) {
  p <- VlnPlot(scRNA, features = method, group.by = "external_annotation", pt.size = 0,
               cols = celltype_colors) +
    ggtitle(paste(method, "Score of Key Genes")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "none")
  ggsave(file.path(dir_plots, paste0("01_Violin_", method, ".pdf")), p, width = 10, height = 6)
  ggsave(file.path(dir_plots, paste0("01_Violin_", method, ".png")),
         p,
         width = 10,
         height = 6,
         dpi = 300)
}

meta_data_long <- scRNA@meta.data %>%
  dplyr::select(external_annotation, all_of(methods)) %>%
  pivot_longer(cols = -external_annotation, names_to = "Method", values_to = "Score")
write.csv(meta_data_long,
          file.path(dir_csv, "02_Violin_Combined_Methods_Data.csv"),
          row.names = FALSE)

p_combined <- ggplot(meta_data_long,
                     aes(x = external_annotation, y = Score, fill = external_annotation)) +
  geom_violin(trim = FALSE, scale = "width", linewidth = 0.3) +
  facet_wrap(~ Method, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = celltype_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Key Gene Scores by Method and Cell Type", x = "", y = "Score")
ggsave(file.path(dir_plots, "02_Violin_Combined_Methods.pdf"), p_combined, width = 10, height = 15)
ggsave(file.path(dir_plots, "02_Violin_Combined_Methods.png"),
       p_combined,
       width = 10,
       height = 15,
       dpi = 300)

p_combined_violin <- VlnPlot(scRNA,
                             features = "KeyGenes_Combined",
                             group.by = "external_annotation",

                             pt.size = 0, cols = celltype_colors) +
  ggtitle("Key Gene Combined Activity Score by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(dir_plots, "03_Violin_Combined_FinalScore_by_CellType.pdf"),
       p_combined_violin,
       width = 10,
       height = 6)
ggsave(file.path(dir_plots, "03_Violin_Combined_FinalScore_by_CellType.png"),
       p_combined_violin,
       width = 10,
       height = 6,
       dpi = 300)

p_group_violin <- VlnPlot(scRNA, features = "KeyGenes_Combined", group.by = "external_annotation",
                          split.by = "group", pt.size = 0) +
  scale_fill_manual(values = group_colors) +
  ggtitle("Key Gene Combined Score: UC vs HC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(dir_plots, "04_Violin_Combined_FinalScore_by_Group.pdf"),
       p_group_violin,
       width = 12,
       height = 6)
ggsave(file.path(dir_plots, "04_Violin_Combined_FinalScore_by_Group.png"),
       p_group_violin,
       width = 12,
       height = 6,
       dpi = 300)

write.csv(scRNA@meta.data %>% dplyr::select(external_annotation, group, KeyGenes_Combined),
          file.path(dir_csv, "03_04_Violin_FinalScore_Data.csv"), row.names = FALSE)

summary_table_celltype <- scRNA@meta.data %>%
  group_by(external_annotation) %>%
  summarise(AUCell_Mean = mean(AUCell, na.rm = TRUE),
            UCell_Mean = mean(UCell, na.rm = TRUE),
            singscore_Mean = mean(singscore, na.rm = TRUE),
            ssGSEA_Mean = mean(ssGSEA, na.rm = TRUE),
            AddModuleScore_Mean = mean(AddModuleScore, na.rm = TRUE),
            KeyGenes_Combined_Mean = mean(KeyGenes_Combined, na.rm = TRUE))
write.csv(summary_table_celltype,
          file.path(dir_csv, "Summary_Mean_Scores_by_Celltype.csv"),
          row.names = FALSE)

summary_table_group <- scRNA@meta.data %>%
  group_by(external_annotation, group) %>%
  summarise(AUCell_Mean = mean(AUCell, na.rm = TRUE),
            UCell_Mean = mean(UCell, na.rm = TRUE),
            singscore_Mean = mean(singscore, na.rm = TRUE),
            ssGSEA_Mean = mean(ssGSEA, na.rm = TRUE),
            AddModuleScore_Mean = mean(AddModuleScore, na.rm = TRUE),
            KeyGenes_Combined_Mean = mean(KeyGenes_Combined, na.rm = TRUE),
            .groups = "drop")
write.csv(summary_table_group,
          file.path(dir_csv, "Summary_Mean_Scores_by_Celltype_and_Group.csv"),
          row.names = FALSE)
