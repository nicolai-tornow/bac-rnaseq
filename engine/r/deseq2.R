#!/usr/bin/env Rscript
suppressMessages({library(DESeq2)})
args <- commandArgs(trailingOnly = TRUE)
counts_f <- args[1]; coldata_f <- args[2]; contrasts_f <- args[3]
out_dir <- args[4]; use_batch <- as.integer(args[5])

dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

cnt <- read.delim(counts_f, row.names = 1, check.names = FALSE)
col <- read.delim(coldata_f, row.names = 1, check.names = FALSE)
cnt <- cnt[, rownames(col), drop = FALSE]        # align sample order
col$condition <- factor(col$condition)
design <- ~condition
if (use_batch == 1 && "batch" %in% colnames(col) && nlevels(factor(col$batch)) > 1) {
  col$batch <- factor(col$batch)
  design <- ~batch + condition
}

dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt)), colData = col, design = design)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

write.table(data.frame(Gene = rownames(counts(dds, normalized = TRUE)),
                       counts(dds, normalized = TRUE), check.names = FALSE),
            file.path(out_dir, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
vsd <- assay(vst(dds, blind = FALSE))
write.table(data.frame(Gene = rownames(vsd), vsd, check.names = FALSE),
            file.path(out_dir, "vst_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sf <- sizeFactors(dds)
write.table(data.frame(sample = names(sf), size_factor = sf),
            file.path(out_dir, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cons <- read.delim(contrasts_f, header = FALSE, col.names = c("name", "num", "den"))
for (i in seq_len(nrow(cons))) {
  con <- c("condition", as.character(cons$num[i]), as.character(cons$den[i]))
  res <- results(dds, contrast = con, alpha = 0.05)
  res <- lfcShrink(dds, contrast = con, res = res, type = "normal")
  res <- res[order(res$padj, res$pvalue), ]
  df <- data.frame(Gene = rownames(res), baseMean = res$baseMean,
                   log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                   stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  write.table(df, file.path(out_dir, "results", paste0(cons$name[i], ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
