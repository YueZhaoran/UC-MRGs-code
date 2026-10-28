library(mclust)

if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

bs <- read.csv("09_spot_cluster.csv")
seu <- read.csv(file.path(dirname(getwd()), "02cluster", "02_spot_cluster.csv"))
merged <- merge(bs[, c("barcode", "section", "group", "cluster")],
                seu[, c("barcode", "section", "cluster")],
                by = c("barcode", "section"), suffixes = c("_bs", "_seurat"))

ari_tab <- data.frame()
for (sec in unique(merged$section)) {
  sub <- merged[merged$section == sec, ]
  tab <- table(sub$cluster_bs, sub$cluster_seurat)
  jac <- max(sapply(seq_len(nrow(tab)), function(i) {
    max(sapply(seq_len(ncol(tab)), function(j) {
      inter <- tab[i, j]
      if (inter == 0) return(0)
      inter / (sum(tab[i, ]) + sum(tab[, j]) - inter)
    }))
  }))
  ari_tab <- rbind(ari_tab, data.frame(
    section = sec, group = sub$group[1],
    n_spot = nrow(sub),
    n_cluster_bayesspace = length(unique(sub$cluster_bs)),
    n_cluster_seurat = length(unique(sub$cluster_seurat)),
    ARI = round(mclust::adjustedRandIndex(sub$cluster_bs, sub$cluster_seurat), 3),
    max_jaccard = round(jac, 3)))
}
ari_tab$ARI_all_sections <- round(mclust::adjustedRandIndex(merged$cluster_bs, merged$cluster_seurat),
                                  3)
write.csv(ari_tab, "09_cluster_ari.csv", row.names = FALSE)
