rm(list = ls())
gc()

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir <- file.path(getwd(), "..")

library(dplyr)
library(Seurat)
library(scTenifoldNet)
library(ggplot2)
library(ggrepel)
library(forcats)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(KEGG.db)
library(AnnotationDbi)

outdir <- "."
target_celltype <- "Stromal"
genes_to_oe <- readLines("target_genes.txt")
nc_nNet   <- 10
nc_nCells <- 500
nc_nComp  <- 3
seed_use  <- 123
kegg_src  <- NA_character_
n_cell_input <- NA_integer_
n_gene_network <- NA_integer_

fb <- readRDS(file.path(project_dir, "01singleCell", "05_scData_annotated.rds"))
Idents(fb) <- fb$external_annotation

scData <- subset(fb, idents = target_celltype, subset = group == "UC")
n_cell_input <- ncol(scData)

drop_subtypes <- c("Endothelial_cells", "Neurons")
scData <- subset(scData, subset = !labels %in% drop_subtypes)

if ("JoinLayers" %in% ls("package:Seurat")) scData <- JoinLayers(scData, assay = "RNA")
DefaultAssay(scData) <- "RNA"

scData <- NormalizeData(scData, assay = "RNA", normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
scData <- FindVariableFeatures(scData, assay = "RNA", selection.method = "vst",
                               nfeatures = 2000, verbose = FALSE)
hvg <- unique(VariableFeatures(scData))

countMatrix <- GetAssayData(scData, assay = "RNA", layer = "counts")

keep_genes <- Matrix::rowSums(countMatrix > 0) >= 5 | rownames(countMatrix) %in% genes_to_oe
keep_cells <- Matrix::colSums(countMatrix > 0) >= 200
countMatrix <- countMatrix[keep_genes, keep_cells, drop = FALSE]

nCores <- min(parallel::detectCores(), 4)

for (gene in genes_to_oe) {

  if (!gene %in% rownames(countMatrix)) {
    warning("Gene ", gene, " not found in count matrix. Skip.")
    next
  }

  expr_cells <- sum(countMatrix[gene, ] > 0)
  if (expr_cells < 10) {
    warning("Gene ", gene, " expressed in fewer than 10 cells. Skip.")
    next
  }

  hvg_use <- unique(c(hvg, genes_to_oe))
  hvg_use <- intersect(hvg_use, rownames(countMatrix))
  safe_min_cells <- max(50, round(ncol(countMatrix) * 0.03))
  hvg_use <- hvg_use[Matrix::rowSums(countMatrix[hvg_use, , drop = FALSE] > 0) >= safe_min_cells |
                       hvg_use %in% genes_to_oe]
  if (!gene %in% hvg_use) hvg_use <- unique(c(hvg_use, gene))
  n_gene_network <- length(hvg_use)

  countMatrix_hvg <- countMatrix[hvg_use, , drop = FALSE]
  countMatrix_hvg <- as(countMatrix_hvg, "dgCMatrix")

  countMatrix_oe <- countMatrix_hvg
  countMatrix_oe[gene, ] <- countMatrix_oe[gene, ] * 2

  set.seed(seed_use)
  result <- tryCatch({
    scTenifoldNet(X = countMatrix_hvg, Y = countMatrix_oe, qc = FALSE,
                  nc_nNet = nc_nNet, nc_nCells = nc_nCells, nc_nComp = nc_nComp,
                  nCores = nCores)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(result)) next

  df_OE <- result[["diffRegulation"]]
  if (is.null(df_OE) || nrow(df_OE) == 0) {
    warning("No diffRegulation result for ", gene)
    next
  }
  if (!"gene" %in% colnames(df_OE)) df_OE$gene <- rownames(df_OE)

  df_OE$logFC <- ifelse(df_OE$FC > 0, log2(df_OE$FC), NA)
  df_OE$neg_log10_pval <- ifelse(df_OE$p.value > 0, -log10(df_OE$p.value), NA)
  df_OE$sig <- ifelse(df_OE$p.value < 0.05, "P < 0.05", "NS")

  write.csv(df_OE, file.path(outdir, paste0(gene, "_OE_diffRegulation.csv")), row.names = FALSE)
  save(df_OE, file = file.path(outdir, paste0(gene, "_OE_diffRegulation.rda")))

  df_plot <- df_OE %>% filter(!is.na(logFC), is.finite(logFC), !is.na(neg_log10_pval))
  if (nrow(df_plot) > 0) {
    top_label <- df_plot %>% arrange(p.value) %>% head(10)
    p1 <- ggplot(df_plot, aes(x = logFC, y = neg_log10_pval)) +
      geom_point(aes(size = abs(logFC), color = logFC), alpha = 0.75) +
      scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B",
                            midpoint = 0, name = "log2FC") +
      scale_size_continuous(range = c(1.5, 6), name = "|log2FC|") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      ggrepel::geom_text_repel(data = top_label, aes(label = gene), size = 3.5,
                               max.overlaps = Inf) +
      labs(x = "log2 Fold Change", y = "-log10(P value)",
           title = paste0(gene, " virtual overexpression in ", target_celltype)) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(outdir, "figures", paste0(gene, "_OE_volcano.pdf")),
           p1,
           width = 7,
           height = 5.5)
    ggsave(file.path(outdir, "figures", paste0(gene, "_OE_volcano.png")),
           p1,
           width = 7,
           height = 5.5,
           dpi = 300)
  }

  df_top10_up <- df_OE %>%
    filter(!is.na(logFC), is.finite(logFC), !is.na(p.value), p.value > 0) %>%
    arrange(desc(logFC)) %>%
    head(10) %>%
    mutate(neg_log10_pval = -log10(p.value), gene = fct_reorder(gene, logFC))
  if (nrow(df_top10_up) > 0) {
    p2 <- ggplot(df_top10_up, aes(x = logFC, y = gene)) +
      geom_col(aes(fill = neg_log10_pval), width = 0.7) +
      geom_text(aes(label = format(logFC, scientific = TRUE, digits = 3)),
                hjust = -0.15, size = 3.5) +
      scale_fill_gradient(low = "#C6DBEF", high = "#08519C", name = "-log10(P)") +
      labs(x = "log2 Fold Change", y = NULL,
           title = paste0(gene, " OE: Top 10 up-regulated genes")) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      expand_limits(x = max(df_top10_up$logFC, na.rm = TRUE) * 1.3)
    ggsave(file.path(outdir, "figures", paste0(gene, "_OE_top10_up.pdf")),
           p2,
           width = 7,
           height = 5)
    ggsave(file.path(outdir, "figures", paste0(gene, "_OE_top10_up.png")),
           p2,
           width = 7,
           height = 5,
           dpi = 300)
  }

  df_top10_down <- df_OE %>%
    filter(!is.na(logFC), is.finite(logFC), !is.na(p.value), p.value > 0) %>%
    arrange(logFC) %>%
    head(10) %>%
    mutate(neg_log10_pval = -log10(p.value), gene = fct_reorder(gene, logFC))
  if (nrow(df_top10_down) > 0) {
    p3 <- ggplot(df_top10_down, aes(x = logFC, y = gene)) +
      geom_col(aes(fill = neg_log10_pval), width = 0.7) +
      geom_text(aes(label = format(logFC, scientific = TRUE, digits = 3)),
                hjust = 1.1, size = 3.5, color = "white") +
      scale_fill_gradient(low = "#FEE0D2", high = "#A50F15", name = "-log10(P)") +
      labs(x = "log2 Fold Change", y = NULL,
           title = paste0(gene, " OE: Top 10 down-regulated genes")) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      expand_limits(x = min(df_top10_down$logFC, na.rm = TRUE) * 1.3)
    ggsave(file.path(outdir, "figures", paste0(gene, "_OE_top10_down.pdf")),
           p3,
           width = 7,
           height = 5)
    ggsave(file.path(outdir, "figures", paste0(gene, "_OE_top10_down.png")),
           p3,
           width = 7,
           height = 5,
           dpi = 300)
  }

  sig_genes <- df_OE %>%
    filter(p.value < 0.05, !is.na(logFC), is.finite(logFC)) %>%
    arrange(p.value) %>%
    head(100) %>%
    pull(gene)
  sig_genes <- setdiff(sig_genes, gene)
  if (length(sig_genes) < 10) {
    warning("Significant genes less than 10, skip enrichment.")
    next
  }

  gene_trans <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)
  if (is.null(gene_trans) || nrow(gene_trans) == 0) next

  kegg_res <- tryCatch({
    res <- enrichKEGG(gene = gene_trans$ENTREZID, organism = 'hsa', pAdjustMethod = "BH",
                      pvalueCutoff = 0.99, qvalueCutoff = 0.99, use_internal_data = FALSE)
    kegg_src <<- "online KEGG REST API (rest.kegg.jp)"
    res
  }, error = function(e) {
    tryCatch({
      res <- enrichKEGG(gene = gene_trans$ENTREZID, organism = 'hsa',
                        pAdjustMethod = "BH", pvalueCutoff = 0.99,
                        qvalueCutoff = 0.99, use_internal_data = TRUE)
      kegg_src <<- "local KEGG.db (offline annotation)"
      res
    }, error = function(e2) NULL)
  })

  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    kegg_df <- as.data.frame(kegg_res) %>% filter(pvalue < 0.05) %>% arrange(desc(Count))
    if (nrow(kegg_df) > 0) {
      name_map <- unlist(mget(kegg_df$ID, KEGGPATHID2NAME, ifnotfound = NA))
      kegg_df$Description <- ifelse(is.na(name_map[kegg_df$ID]), kegg_df$Description,
                                    as.character(name_map[kegg_df$ID]))
      write.csv(kegg_df, file.path(outdir, paste0(gene, "_OE_KEGG_results.csv")), row.names = FALSE)
      top_kegg <- kegg_df %>% arrange(p.adjust) %>% head(10) %>%
        mutate(GeneRatio_Num = sapply(strsplit(GeneRatio, "/"),
                                      function(x) as.numeric(x[1]) / as.numeric(x[2])),
               Description = fct_reorder(Description, Count))
      p_kegg <- ggplot(top_kegg, aes(x = GeneRatio_Num, y = Description)) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
        scale_size(range = c(3, 8), name = "Count") +
        theme_bw() +
        labs(title = paste0("KEGG Enrichment - ", gene, " OE in ", target_celltype),
             x = "GeneRatio", y = NULL) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      ggsave(file.path(outdir, "figures", paste0(gene, "_OE_KEGG_dotplot.pdf")),
             p_kegg,
             width = 8,
             height = 6)
      ggsave(file.path(outdir, "figures", paste0(gene, "_OE_KEGG_dotplot.png")),
             p_kegg,
             width = 8,
             height = 6,
             dpi = 300)
    } else {
    }
  } else {
  }
}

writeLines(c(
  paste0("celltype=", target_celltype),
  paste0("group=UC"),
  paste0("n_cell_input=", n_cell_input),
  paste0("n_cell_used=", ncol(scData)),
  paste0("drop_subtypes=", paste(drop_subtypes, collapse = ",")),
  paste0("n_gene_network=", n_gene_network),
  paste0("nNet=", nc_nNet),
  paste0("nCells=", nc_nCells),
  paste0("nComp=", nc_nComp),
  paste0("seed=", seed_use),
  paste0("kegg_source=", ifelse(is.na(kegg_src), "none", kegg_src)),
  paste0("run_date=", format(Sys.Date()))), "03_run_meta.txt")
