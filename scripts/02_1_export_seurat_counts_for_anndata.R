#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/02_1_export_seurat_counts_for_anndata.R",
      "<input_barebones_seurat_rds>",
      "<output_sparse_export_dir>",
      "[chunk_size_cells=5000]"
    )
  )
}

input_rds <- args[1]
output_dir <- args[2]
chunk_size <- if (length(args) >= 3) as.integer(args[3]) else 5000L

if (is.na(chunk_size) || chunk_size <= 0L) {
  stop("chunk_size_cells must be a positive integer")
}

cat("input_rds: ", input_rds, "\n", sep = "")
cat("output_sparse_export_dir: ", output_dir, "\n", sep = "")
cat("chunk_size_cells: ", chunk_size, "\n", sep = "")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "peaks"

counts <- GetAssayData(obj, assay = "peaks", layer = "counts")
cat("counts matrix dimensions, peaks x cells:\n")
print(dim(counts))

n_peaks <- nrow(counts)
n_cells <- ncol(counts)

if (n_peaks == 0 || n_cells == 0) {
  stop("Counts matrix is empty.")
}

cell_md <- obj@meta.data[colnames(counts), , drop = FALSE]
if (!identical(rownames(cell_md), colnames(counts))) {
  stop("Cell metadata rownames do not match count matrix column names.")
}

required_cols <- c("cell_type", "library_size")
missing_cols <- setdiff(required_cols, colnames(cell_md))
if (length(missing_cols) > 0) {
  stop(paste("Missing required metadata columns:", paste(missing_cols, collapse = ", ")))
}

cell_md$cell_name <- rownames(cell_md)
cell_md <- cell_md[, c("cell_name", setdiff(colnames(cell_md), "cell_name")), drop = FALSE]
write.table(
  cell_md,
  file = file.path(output_dir, "obs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

peak_names <- rownames(counts)
peak_gr <- StringToGRanges(peak_names, sep = c("-", "-"))
peak_md <- data.frame(
  peak_name = peak_names,
  chrom = as.character(seqnames(peak_gr)),
  start = start(peak_gr),
  end = end(peak_gr),
  stringsAsFactors = FALSE
)
write.table(
  peak_md,
  file = file.path(output_dir, "var.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write_int64_le <- function(x, con) {
  x <- as.numeric(x)

  if (length(x) == 0) {
    return(invisible(NULL))
  }
  if (any(is.na(x))) {
    stop("Cannot write NA values as int64")
  }
  if (any(x < 0)) {
    stop("Negative values are not allowed in write_int64_le")
  }
  if (any(x > 9007199254740991)) {
    stop("Values exceed exact integer range of double (2^53 - 1)")
  }

  x <- floor(x)
  bytes <- raw(length(x) * 8L)
  for (k in 0:7) {
    byte_vals <- floor(x / (256^k)) %% 256
    bytes[seq.int(k + 1L, length(bytes), by = 8L)] <- as.raw(byte_vals)
  }
  writeBin(bytes, con)
}

coerce_block_to_dgC <- function(x) {
  tryCatch(
    as(x, "dgCMatrix"),
    error = function(e1) {
      stop(paste("Failed to coerce chunk to dgCMatrix:", conditionMessage(e1)))
    }
  )
}

data_file <- file.path(output_dir, "data_float32.bin")
indices_file <- file.path(output_dir, "indices_int64.bin")
indptr_file <- file.path(output_dir, "indptr_int64.bin")
shape_file <- file.path(output_dir, "shape.tsv")

con_data <- file(data_file, open = "wb")
con_indices <- file(indices_file, open = "wb")
con_indptr <- file(indptr_file, open = "wb")

on.exit({
  try(close(con_data), silent = TRUE)
  try(close(con_indices), silent = TRUE)
  try(close(con_indptr), silent = TRUE)
}, add = TRUE)

# CSR is cell x peak, so indptr starts with one zero for the first cell row.
write_int64_le(0, con_indptr)

total_nnz <- 0
starts <- seq.int(1L, n_cells, by = chunk_size)

cat("starting chunked raw-count CSR export, cell x peak\n")

for (s in starts) {
  e <- min(s + chunk_size - 1L, n_cells)
  idx <- s:e
  n_chunk_cells <- length(idx)

  cat("processing cells ", s, "-", e, "\n", sep = "")

  # counts block is peak x cell; transpose to cell x peak for AnnData.X.
  block <- counts[, idx, drop = FALSE]
  block <- coerce_block_to_dgC(block)
  block_t <- as(t(block), "TsparseMatrix")

  chunk_nnz <- as.numeric(length(block_t@x))
  cat("chunk_nnz: ", chunk_nnz, "\n", sep = "")

  if (chunk_nnz > 0) {
    ord <- order(block_t@i, block_t@j)
    row_idx0 <- as.numeric(block_t@i[ord])
    col_idx0 <- as.numeric(block_t@j[ord])
    x_vals <- as.numeric(block_t@x[ord])

    writeBin(x_vals, con_data, size = 4, endian = "little")
    write_int64_le(col_idx0, con_indices)

    row_counts <- tabulate(as.integer(row_idx0 + 1), nbins = n_chunk_cells)
    p_chunk <- cumsum(c(0, as.numeric(row_counts)))
  } else {
    p_chunk <- rep(0, n_chunk_cells + 1L)
  }

  p_adj <- p_chunk[-1] + total_nnz
  write_int64_le(p_adj, con_indptr)
  total_nnz <- total_nnz + chunk_nnz

  rm(block, block_t, p_chunk, p_adj)
  gc()
}

close(con_data)
close(con_indices)
close(con_indptr)

shape_df <- data.frame(
  n_rows = n_cells,
  n_cols = n_peaks,
  nnz = total_nnz
)
write.table(
  shape_df,
  file = shape_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("total_nnz: ", total_nnz, "\n", sep = "")
cat("DONE\n")
