args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !file.exists(args[[1]])) {
  stop("Provide the RNA-seq analysis script path.", call. = FALSE)
}

local({
analysis_script <- normalizePath(args[[1]])
original_wd <- getwd()
work_dir <- tempfile("seqcheckflow-rnaseq-smoke-")
dir.create(work_dir)
on.exit({
  setwd(original_wd)
  unlink(work_dir, recursive = TRUE)
}, add = TRUE)
setwd(work_dir)
set.seed(42)

sample_data <- expand.grid(
  cell = c("donor1", "donor2", "donor3", "donor4"),
  dex = c("untrt", "trt"),
  stringsAsFactors = FALSE
)
sample_data <- sample_data[order(sample_data$cell, sample_data$dex), ]
sample_data$sample <- paste(sample_data$cell, sample_data$dex, sep = "_")
genes <- sprintf("ENSG%06d", seq_len(120))
transcripts <- sprintf("ENST%06d", seq_len(120))
write.table(
  data.frame(transcripts, genes),
  "tx2gene.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

manifest <- sample_data[c("sample", "cell", "dex")]
manifest$quant_file <- character(nrow(sample_data))
for (sample_index in seq_len(nrow(sample_data))) {
  treated <- sample_data$dex[[sample_index]] == "trt"
  donor_factor <- match(sample_data$cell[[sample_index]], unique(sample_data$cell))
  effect <- rep(1, length(genes))
  effect[seq_len(20)] <- 4
  effect[21:40] <- 0.25
  counts <- round(rnbinom(
    length(genes),
    mu = 80 * donor_factor * ifelse(treated, effect, 1),
    size = 25
  ))
  quant_file <- paste0(sample_data$sample[[sample_index]], ".quant.sf")
  write.table(
    data.frame(
      Name = transcripts,
      Length = rep(1000, length(genes)),
      EffectiveLength = rep(900, length(genes)),
      TPM = counts / sum(counts) * 1e6,
      NumReads = counts
    ),
    quant_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  manifest$quant_file[[sample_index]] <- quant_file
}
write.table(
  manifest,
  "quant_manifest.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

output <- system2(
  "Rscript",
  c(shQuote(analysis_script), "quant_manifest.tsv", "tx2gene.tsv"),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(output, "status")
if (!is.null(status) && status != 0) {
  stop(paste(output, collapse = "\n"), call. = FALSE)
}

expected <- c(
  "deseq2_results.csv",
  "analysis_summary.txt",
  "ma_plot.png",
  "volcano_plot.png",
  "pca_plot.png",
  "session_info.txt"
)
if (!all(file.exists(expected)) || !all(file.info(expected)$size > 0)) {
  stop("The analysis did not produce all expected outputs.", call. = FALSE)
}
message("RNA-seq DESeq2 smoke test passed.")
})
