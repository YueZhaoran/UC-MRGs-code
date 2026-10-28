rm(list = ls())
gc()

# ===== 路径（RStudio 与 Rscript 都能定位到本模块目录）=====
if (rstudioapi::isAvailable()) {
  module_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  .args <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .args[grep("^--file=", .args)])
  module_dir <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
}
setwd(module_dir)
project_dir <- dirname(module_dir)
cat("工作目录：", getwd(), "\n")

library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)
library(forcats)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(KEGG.db)


# 读取数据
scData <- readRDS(file.path("../", "01singleCell", "05_scData_annotated.rds"))

# 选择关键细胞类型
target_celltype <- "Stromal"
scData <- subset(scData, subset = external_annotation == target_celltype)

# 只保留 UC 组（疾病组）
scData <- subset(scData, subset = group == "UC")
n_cell_input <- ncol(scData)

# SingleR 精细注释把 Endothelial_cells、Neurons 误标进 Stromal（UC 侧实测 98 个内皮 + 2 个神经元），
# 二者非间充质谱系，会污染基因调控网络；与 2.12.4 Monocle3、2.12.3 CellOracle、2.15 scTenifoldNet 同一剔除口径
drop_subtypes <- c("Endothelial_cells", "Neurons")
scData <- subset(scData, subset = !labels %in% drop_subtypes)
cat(target_celltype, "（UC）细胞数：", n_cell_input, " → 剔除 ", paste(drop_subtypes, collapse = "、"),
    " 后 ", ncol(scData), "\n", sep = "")

# 处理 Seurat v5 的图层
if ("JoinLayers" %in% ls("package:Seurat")) {
  scData <- JoinLayers(scData, assay = "RNA")
}
DefaultAssay(scData) <- "RNA"

# 标准化、找高变基因
scData <- NormalizeData(scData, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
scData <- FindVariableFeatures(scData, assay = "RNA", selection.method = "vst", nfeatures = 3000, verbose = FALSE)
hvg <- VariableFeatures(scData)
hvg <- unique(hvg)

# 提取 counts 矩阵
countMatrix <- GetAssayData(scData, assay = "RNA", layer = "counts")

# 目标基因

genes_to_ko <- c("ADAMTS1", "FBN1", "SPARC","VCAM1")

# 网络构建参数（显式写出，与 scTenifoldKnk 包默认值一致；图注数字从这里派生）
nc_nNet   <- 10     # 每个方向构建的 PCR 网络数
nc_nCells <- 500    # 每次网络的细胞自助抽样数（有放回）
nc_nComp  <- 3      # 主成分数
seed_use  <- 123
nCores    <- min(parallel::detectCores(), 4)
kegg_src  <- "local KEGG.db (offline annotation)"
n_gene_network <- NA_integer_

# 过滤低表达基因：至少 5 个细胞表达，或为目标基因
keep_genes <- Matrix::rowSums(countMatrix > 0) >= 5 | rownames(countMatrix) %in% genes_to_ko
keep_cells <- Matrix::colSums(countMatrix > 0) >= 200
countMatrix <- countMatrix[keep_genes, keep_cells, drop = FALSE]

# 循环每个基因进行虚拟敲除
for (gene in genes_to_ko) {
  message("\n======================================")
  message("Running virtual knockout for: ", gene)
  
  if (!gene %in% rownames(countMatrix)) {
    warning("Gene ", gene, " not found in count matrix. Skip.")
    next
  }
  
  expr_cells <- sum(countMatrix[gene, ] > 0)
  cat(gene, " expressed in ", expr_cells, " cells\n")
  if (expr_cells < 10) {
    warning("Gene ", gene, " expressed in fewer than 10 cells. Skip.")
    next
  }
  
  # 构建用于 scTenifoldKnk 的基因集（高变基因 + 目标基因）
  hvg_use <- unique(c(hvg, genes_to_ko))
  hvg_use <- intersect(hvg_use, rownames(countMatrix))
  
  safe_min_cells <- max(50, round(ncol(countMatrix) * 0.03))
  hvg_use <- hvg_use[Matrix::rowSums(countMatrix[hvg_use, , drop = FALSE] > 0) >= safe_min_cells | hvg_use %in% genes_to_ko]
  if (!gene %in% hvg_use) hvg_use <- unique(c(hvg_use, gene))
  
  countMatrix_hvg <- countMatrix[hvg_use, , drop = FALSE]
  countMatrix_hvg <- as(countMatrix_hvg, "dgCMatrix")
  n_gene_network <- nrow(countMatrix_hvg)
  cat("参与网络构建的基因数：", n_gene_network, "\n")

  # 运行 scTenifoldKnk
  set.seed(seed_use)
  result <- tryCatch({
    scTenifoldKnk(countMatrix = countMatrix_hvg, gKO = gene, qc = FALSE,
                  nc_nNet = nc_nNet, nc_nCells = nc_nCells, nc_nComp = nc_nComp,
                  nCores = nCores)
  }, error = function(e) {
    message("scTenifoldKnk failed for ", gene, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(result)) next
  
  df_KO <- result[["diffRegulation"]]
  if (is.null(df_KO) || nrow(df_KO) == 0) {
    warning("No diffRegulation result for ", gene)
    next
  }
  
  if (!"gene" %in% colnames(df_KO)) df_KO$gene <- rownames(df_KO)
  df_KO$logFC <- ifelse(df_KO$FC > 0, log2(df_KO$FC), NA)
  
  # 保存结果
  write.csv(df_KO, file.path(paste0(gene, "_df_KO.result.csv")), row.names = FALSE)
  save(df_KO, file = file.path(paste0(gene, "_df_KO.rda")))
  
  # 绘制气泡图（火山图风格）
  df_plot <- df_KO %>%
    filter(!is.na(logFC), !is.na(p.value), p.value > 0) %>%
    mutate(neg_log10_pval = -log10(p.value),
           sig = ifelse(p.value < 0.05, "P < 0.05", "NS"))
  
  if (nrow(df_plot) == 0) next
  
  top_label <- df_plot %>% arrange(p.value) %>% head(10)
  
  p1 <- ggplot(df_plot, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(size = abs(logFC), color = logFC), alpha = 0.75) +
    scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0, name = "log2FC") +
    scale_size_continuous(range = c(1.5, 6), name = "|log2FC|") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(data = top_label, aes(label = gene), size = 3.5, max.overlaps = Inf) +
    labs(x = "log2 Fold Change", y = "-log10(P value)", title = paste0(gene, " virtual knockout")) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(paste0(gene, "_bubble_plot.pdf")), p1, width = 7, height = 5.5)
  
  # Top10 高响应柱状图（扰动幅度最大，非传统上调）
  df_top10_high <- df_KO %>%
    filter(!is.na(logFC), !is.na(p.value), p.value > 0) %>%
    arrange(desc(logFC)) %>%
    head(10) %>%
    mutate(neg_log10_pval = -log10(p.value), gene = fct_reorder(gene, logFC))

  if (nrow(df_top10_high) > 0) {
    p2 <- ggplot(df_top10_high, aes(x = logFC, y = gene)) +
      geom_col(aes(fill = neg_log10_pval), width = 0.7) +
      geom_text(aes(label = format(logFC, scientific = TRUE, digits = 3)), hjust = -0.15, size = 3.5) +
      scale_fill_gradient(low = "#C6DBEF", high = "#08519C", name = "-log10(P)") +
      labs(x = "log2 Fold Change", y = NULL, title = paste0(gene, " KO: Top 10 highest-response genes")) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      expand_limits(x = max(df_top10_high$logFC, na.rm = TRUE) * 1.3)
    ggsave(file.path(paste0(gene, "_TOP10_high_response_bar_plot.pdf")), p2, width = 7, height = 5)
  }

  # Top10 低响应柱状图（扰动幅度最小，非传统下调）
  df_top10_low <- df_KO %>%
    filter(!is.na(logFC), !is.na(p.value), p.value > 0) %>%
    arrange(logFC) %>%
    head(10) %>%
    mutate(neg_log10_pval = -log10(p.value), gene = fct_reorder(gene, logFC))

  if (nrow(df_top10_low) > 0) {
    p3 <- ggplot(df_top10_low, aes(x = logFC, y = gene)) +
      geom_col(aes(fill = neg_log10_pval), width = 0.7) +
      geom_text(aes(label = format(logFC, scientific = TRUE, digits = 3)), hjust = 1.1, size = 3.5, color = "white") +
      scale_fill_gradient(low = "#FEE0D2", high = "#A50F15", name = "-log10(P)") +
      labs(x = "log2 Fold Change", y = NULL, title = paste0(gene, " KO: Top 10 lowest-response genes")) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      expand_limits(x = min(df_top10_low$logFC, na.rm = TRUE) * 1.3)
    ggsave(file.path(paste0(gene, "_TOP10_low_response_bar_plot.pdf")), p3, width = 7, height = 5)
  }
  
  # 保存全部差异基因按P值排序
  df_top_sig <- df_KO %>% filter(!is.na(logFC), !is.na(p.value)) %>% arrange(p.value)
  write.csv(df_top_sig, file.path(paste0(gene, "_KO_all_genes_sorted_by_pvalue.csv")), row.names = FALSE)
  
  # --- 富集分析 ---
  message("Running GO and KEGG enrichment for ", gene)
  
  sig_genes <- df_KO %>% filter(p.value < 0.05, !is.na(logFC), abs(logFC) > 0.5) %>% pull(gene)
  if (length(sig_genes) < 10) {
    warning("Significant genes less than 10, skip enrichment.")
    next
  }
  
  gene_trans <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)
  
  if (is.null(gene_trans) || nrow(gene_trans) == 0) next
  
  # GO 富集
  go_res <- enrichGO(gene = gene_trans$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(go_res) && nrow(go_res) > 0) {
    write.csv(as.data.frame(go_res), file.path(paste0(gene, "_GO_ALL_results.csv")), row.names = FALSE)
    p_go <- dotplot(go_res, showCategory = 5, split = "ONTOLOGY", title = paste0("GO Enrichment - ", gene, " KO")) +
      facet_grid(ONTOLOGY ~ ., scale = "free") +
      theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(paste0(gene, "_GO_ALL_dotplot.pdf")), p_go, width = 8, height = 8)
  }
  
  # KEGG 富集（使用本地数据避免网络问题）
  kegg_res <- tryCatch({
    enrichKEGG(gene = gene_trans$ENTREZID, organism = 'hsa', pAdjustMethod = "BH",
               pvalueCutoff = 0.99, qvalueCutoff = 0.99, use_internal_data = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    kegg_df <- as.data.frame(kegg_res) %>% filter(pvalue < 0.05) %>% arrange(desc(Count))
    if (nrow(kegg_df) > 0) {
      # 走本地 KEGG.db 时 Description 会是通路编号（setReadable 无效），用 KEGGPATHID2NAME 补回名称
      name_map <- unlist(mget(kegg_df$ID, KEGGPATHID2NAME, ifnotfound = NA))
      kegg_df$Description <- ifelse(is.na(name_map[kegg_df$ID]), kegg_df$Description,
                                    as.character(name_map[kegg_df$ID]))
      write.csv(kegg_df, file.path(paste0(gene, "_KEGG_results.csv")), row.names = FALSE)
      top_kegg <- kegg_df %>% arrange(p.adjust) %>% head(10) %>%
        mutate(GeneRatio_Num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
               Description = fct_reorder(Description, Count))
      p_kegg <- ggplot(top_kegg, aes(x = GeneRatio_Num, y = Description)) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
        scale_size(range = c(3, 8), name = "Count") +
        theme_bw() +
        labs(title = paste0("KEGG Enrichment - ", gene, " KO"), x = "GeneRatio", y = NULL) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      ggsave(file.path(paste0(gene, "_KEGG_dotplot.pdf")), p_kegg, width = 8, height = 6)
    }
  }
  message("Finished KO for: ", gene)
}

# ===== 运行元信息（02_ko_figures.R 的图注数字全部从这里取，不手打）=====
writeLines(c(
  paste0("celltype=", target_celltype),
  paste0("group=UC"),
  paste0("n_cell_input=", n_cell_input),
  paste0("n_cell_used=", ncol(scData)),
  paste0("drop_subtypes=", paste(drop_subtypes, collapse = ",")),
  paste0("n_gene_network=", n_gene_network),
  paste0("nhvg=3000"),
  paste0("nNet=", nc_nNet),
  paste0("nCells=", nc_nCells),
  paste0("nComp=", nc_nComp),
  paste0("seed=", seed_use),
  paste0("drg_rule=p<0.05 & |log2FC|>0.5"),
  paste0("go_rule=enrichGO ont=ALL, p<0.05, BH"),
  paste0("kegg_source=", kegg_src),
  paste0("run_date=", format(Sys.Date()))), "03_run_meta.txt")
cat("\n已写出 03_run_meta.txt\n")
print(readLines("03_run_meta.txt"))

