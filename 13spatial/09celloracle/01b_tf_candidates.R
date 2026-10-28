rm(list = ls()); gc()

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
project_dir <- file.path(getwd(), "..", "..")
setwd(getwd())

key_genes   <- c("ADAMTS1", "FBN1", "SPARC", "VCAM1")
min_targets <- 2

trrust_file <- file.path(project_dir, "08molecular_network", "01.TF_KeyGene_Network.csv")

out_all   <- "01_tf_candidates_all.csv"
out_least <- "01_tf_candidates_at_least2.csv"
out_sum   <- "01_tf_candidates_summary.csv"

net <- read.csv(trrust_file, stringsAsFactors = FALSE)

net <- data.frame(TF = toupper(trimws(net$TF)),
                  target = toupper(trimws(net$Target)),
                  direction = net$Direction,
                  stringsAsFactors = FALSE)

net <- net[net$target %in% key_genes, ]

exist_targets <- sort(unique(net$target))
missing <- setdiff(key_genes, exist_targets)
if (length(missing) > 0) cat("警告：以下关键基因在 TRRUST 表中无记录：", paste(missing, collapse = ", "), "\n")

tf_summary <- do.call(rbind, lapply(split(net, net$TF), function(d) {
  data.frame(TF = d$TF[1],
             target_count = length(unique(d$target)),
             target_list = paste(sort(unique(d$target)), collapse = ", "),
             direction = paste(sort(unique(d$direction)), collapse = "/"),
             stringsAsFactors = FALSE)
}))
tf_summary <- tf_summary[order(-tf_summary$target_count, tf_summary$TF), ]

for (g in key_genes) {
  tf_summary[[paste0("hits_", g)]] <- as.integer(sapply(tf_summary$TF, function(tf) {
    any(net$TF == tf & net$target == g)
  }))
}

write.csv(tf_summary, out_sum, row.names = FALSE)
write.csv(tf_summary[, c("TF", "target_count", "target_list", "direction")],
          out_all,
          row.names = FALSE)

cand <- tf_summary[tf_summary$target_count >= min_targets, ]
write.csv(cand, out_least, row.names = FALSE)

for (g in key_genes) {
  n <- sum(tf_summary[[paste0("hits_", g)]])
}
