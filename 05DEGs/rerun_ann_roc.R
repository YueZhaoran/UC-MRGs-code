rm(list = ls()); gc()
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grep("^--file=", .args)])
  setwd(dirname(normalizePath(.file)))
}
library(neuralnet)
library(pROC)

# 读取关键基因
hub_genes <- read.csv("intersection_gene.csv", stringsAsFactors = FALSE)
key_genes <- hub_genes$symbol   # VCAM1, SPARC, FBN1, ADAMTS1

train_data <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
train_group <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(train_group) <- c("sample", "group")
val_data <- read.csv("../00bulkdata/dat.GSE75214.csv", row.names = 1, check.names = FALSE)
val_group <- read.csv("../00bulkdata/group.GSE75214.csv", stringsAsFactors = FALSE)
colnames(val_group) <- c("sample", "group")

train_df <- as.data.frame(t(train_data[key_genes, ])); train_df$label <- ifelse(train_group$group == "UC", 1, 0)
val_df   <- as.data.frame(t(val_data[key_genes, ]));   val_df$label   <- ifelse(val_group$group == "UC", 1, 0)

formula_obj <- as.formula(paste("label ~", paste(key_genes, collapse = " + ")))

# 正确的外部验证：仅用训练集计算 mean/SD 并锁定，验证集套用同一组参数
tr_mean <- colMeans(train_df[, key_genes])
tr_sd   <- apply(train_df[, key_genes], 2, sd)

train_df_scaled <- train_df
train_df_scaled[, key_genes] <- scale(train_df[, key_genes], center = tr_mean, scale = tr_sd)

set.seed(798)
ann_model <- neuralnet(formula_obj, data = train_df_scaled, hidden = 5, linear.output = FALSE)

# 训练集 ROC
train_pred <- compute(ann_model, train_df_scaled[, key_genes])
roc_train <- roc(train_df_scaled$label, train_pred$net.result[, 1], quiet = TRUE)

# 验证集：套用训练集锁定的 mean/SD
val_df_scaled <- val_df
val_df_scaled[, key_genes] <- scale(val_df[, key_genes], center = tr_mean, scale = tr_sd)
val_pred <- compute(ann_model, val_df_scaled[, key_genes])
roc_val <- roc(val_df_scaled$label, val_pred$net.result[, 1], quiet = TRUE)

# 出图，样式与原报告一致（标题、颜色、legend）
pdf("Rerun_Training_ANN_ROC.pdf", width = 6, height = 6)
plot(roc_train, col = "#D55E00", main = "Training Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_train), 3)), col = "#D55E00", lwd = 2)
dev.off()
png("Rerun_Training_ANN_ROC.png", width = 6, height = 6, units = "in", res = 600)
plot(roc_train, col = "#D55E00", main = "Training Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_train), 3)), col = "#D55E00", lwd = 2)
dev.off()

pdf("Rerun_Validation_ANN_ROC.pdf", width = 6, height = 6)
plot(roc_val, col = "#0072B2", main = "Validation Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_val), 3)), col = "#0072B2", lwd = 2)
dev.off()
png("Rerun_Validation_ANN_ROC.png", width = 6, height = 6, units = "in", res = 600)
plot(roc_val, col = "#0072B2", main = "Validation Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_val), 3)), col = "#0072B2", lwd = 2)
dev.off()

cat("训练集 AUC:", round(auc(roc_train), 3), "\n")
cat("验证集 AUC:", round(auc(roc_val), 3), "\n")