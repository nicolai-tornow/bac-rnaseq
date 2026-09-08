#!/usr/bin/env Rscript
suppressMessages({library(sva)})
args <- commandArgs(trailingOnly = TRUE)
mat_f <- args[1]; meta_f <- args[2]; out_dir <- args[3]

mat <- as.matrix(read.delim(mat_f, row.names = 1, check.names = FALSE))
meta <- read.delim(meta_f, row.names = 1, check.names = FALSE)
mat <- mat[, rownames(meta), drop = FALSE]
batch <- factor(meta$batch)
mod <- model.matrix(~ factor(meta$condition))

corrected <- ComBat(dat = mat, batch = batch, mod = mod)
write.table(data.frame(Gene = rownames(corrected), corrected, check.names = FALSE),
            file.path(out_dir, "corrected_vst.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pca_plot <- function(m, path, ttl) {
  p <- prcomp(t(m), scale. = TRUE)
  pdf(path, width = 5, height = 4)
  plot(p$x[, 1], p$x[, 2], col = as.integer(batch), pch = 19,
       xlab = "PC1", ylab = "PC2", main = ttl)
  legend("topright", legend = levels(batch), col = seq_along(levels(batch)), pch = 19)
  dev.off()
}
pca_plot(mat, file.path(out_dir, "pca_before.pdf"), "PCA before ComBat")
pca_plot(corrected, file.path(out_dir, "pca_after.pdf"), "PCA after ComBat")
