#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/02_5_export_cell_metadata.R",
      "<input_barebones_seurat_rds>",
      "<output_metadata_tsv>"
    )
  )
}

input_rds <- args[1]
output_tsv <- args[2]

dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
cells <- colnames(obj)

if (length(cells) == 0) {
  stop("No cells found in Seurat object.")
}

meta <- obj@meta.data
meta <- meta[cells, , drop = FALSE]

library_size <- rep(NA_real_, length(cells))
if ("library_size" %in% colnames(meta)) {
  library_size <- suppressWarnings(as.numeric(meta$library_size))
} else if ("nCount_peaks" %in% colnames(meta)) {
  library_size <- suppressWarnings(as.numeric(meta$nCount_peaks))
}

cell_type <- rep(NA_character_, length(cells))
if ("cell_type" %in% colnames(meta)) {
  cell_type <- as.character(meta$cell_type)
}

out <- data.frame(
  cell_id = cells,
  cell_type = cell_type,
  library_size = library_size,
  stringsAsFactors = FALSE
)

write.table(
  out,
  file = output_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("DONE\n")
