suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
})

data("airway", package = "airway")
sample_data <- as.data.frame(colData(airway))
sample_data$cell <- factor(sample_data$cell)
sample_data$dex <- relevel(factor(sample_data$dex), ref = "untrt")

if (!all(c("untrt", "trt") %in% levels(sample_data$dex))) {
  stop("Expected dex treatment levels 'untrt' and 'trt'.", call. = FALSE)
}

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, quiet = TRUE)

result <- results(dds, contrast = c("dex", "trt", "untrt"))
result_table <- as.data.frame(result)
result_table$gene_id <- rownames(result_table)
result_table <- result_table[
  order(result_table$padj, result_table$pvalue, na.last = TRUE),
  c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
]
write.csv(result_table, "deseq2_results.csv", row.names = FALSE, na = "")

significant <- result_table[!is.na(result_table$padj) & result_table$padj < 0.05, ]
upregulated <- sum(significant$log2FoldChange > 0)
downregulated <- sum(significant$log2FoldChange < 0)
sample_count <- ncol(dds)
donor_count <- nlevels(sample_data$cell)

summary_lines <- c(
  "Dataset: Bioconductor airway RNA-seq count data (GEO GSE52778)",
  "Comparison: dexamethasone-treated (trt) vs untreated (untrt)",
  "Model: ~ cell + dex; cell line is included to account for paired donors",
  sprintf("Samples: %d; donor cell lines: %d", sample_count, donor_count),
  "Significance threshold: Benjamini-Hochberg adjusted p-value < 0.05",
  sprintf("Significant genes: %d (%d higher with dex, %d lower with dex)",
          nrow(significant), upregulated, downregulated),
  paste(
    "Conclusion: dexamethasone treatment is associated with differential",
    "gene expression in this airway smooth-muscle dataset after accounting",
    "for donor cell line. This result describes this dataset; it is not",
    "independent validation or evidence of a general treatment effect."
  )
)
writeLines(summary_lines, "analysis_summary.txt")

png("ma_plot.png", width = 1200, height = 900, res = 150)
plotMA(result, ylim = c(-5, 5), main = "Dexamethasone vs untreated")
dev.off()

adjusted_p <- pmax(result_table$padj, .Machine$double.xmin, na.rm = FALSE)
volcano_y <- -log10(adjusted_p)
volcano_y[is.na(result_table$padj)] <- NA_real_
png("volcano_plot.png", width = 1200, height = 900, res = 150)
plot(
  result_table$log2FoldChange,
  volcano_y,
  pch = 16,
  cex = 0.45,
  col = ifelse(!is.na(result_table$padj) & result_table$padj < 0.05, "firebrick", "grey40"),
  xlab = "Log2 fold change (treated / untreated)",
  ylab = "-Log10 adjusted p-value",
  main = "Airway RNA-seq differential expression"
)
abline(h = -log10(0.05), lty = 2)
dev.off()

writeLines(capture.output(sessionInfo()), "session_info.txt")
