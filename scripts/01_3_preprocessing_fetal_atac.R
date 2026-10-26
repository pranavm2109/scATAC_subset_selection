#Step 1: libs
library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(BPCells)
library(R.utils)

print("libs done")

#Step 2: Files
raw_dir <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/raw/fetal_atac"
unzipped_dir <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/fetal_atac/unzipped_rds"
bp_dir <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/fetal_atac/bpcells_counts"

dir.create(unzipped_dir, recursive = TRUE, showWarnings = FALSE)

gz_files <- sort(list.files(
  raw_dir,
  pattern = "\\.seurat\\.RDS\\.gz$",
  full.names = TRUE
))

print(gz_files)
print(length(gz_files))

#Step 3: Gunzip all files to processed directory
unzipped_files <- character(length(gz_files))

for (i in seq_along(gz_files)) {
  gz_file <- gz_files[i]
  out_file <- file.path(unzipped_dir, sub("\\.gz$", "", basename(gz_file)))

  R.utils::gunzip(
    filename = gz_file,
    destname = out_file,
    overwrite = TRUE,
    remove = FALSE
  )

  unzipped_files[i] <- out_file
}

print("all files decompressed")
print(unzipped_files)

#Step 4: Read peaks assay from each object, write each to BPCells, reopen, and combine
part_dir <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/fetal_atac/bpcells_parts"
dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)

bp_list <- vector("list", length(unzipped_files))
meta_list <- vector("list", length(unzipped_files))

for (i in seq_along(unzipped_files)) {
  f <- unzipped_files[i]
  obj <- readRDS(f)

  counts_mat <- GetAssayData(obj, assay = "peaks", layer = "counts")
  md <- obj@meta.data[, "cell_type", drop = FALSE]

  md <- md[colnames(counts_mat), , drop = FALSE]
  stopifnot(identical(rownames(md), colnames(counts_mat)))

  this_bp_dir <- file.path(
    part_dir,
    sub("\\.seurat\\.RDS$", "", basename(f))
  )

  if (!dir.exists(this_bp_dir)) {
    write_matrix_dir(
      mat = counts_mat,
      dir = this_bp_dir
    )
  }

  bp_list[[i]] <- open_matrix_dir(dir = this_bp_dir)
  meta_list[[i]] <- md

  print(paste("done:", basename(f), "features =", nrow(counts_mat), "cells =", ncol(counts_mat)))

  rm(obj, counts_mat, md)
  gc()
}

combined_counts <- do.call(cbind, bp_list)
combined_md <- do.call(rbind, meta_list)

stopifnot(ncol(combined_counts) == nrow(combined_md))
stopifnot(identical(colnames(combined_counts), rownames(combined_md)))

print(dim(combined_counts))
print(dim(combined_md))

rm(bp_list, meta_list)
gc()

#Step 5: Write combined matrix to BPCells and reopen
if (!dir.exists(bp_dir)) {
  write_matrix_dir(
    mat = combined_counts,
    dir = bp_dir
  )
}

fetal_mat <- open_matrix_dir(dir = bp_dir)

print(fetal_mat)
print(dim(fetal_mat))

rm(combined_counts)
gc()

#Step 6: Create object
fetal_atac <- CreateSeuratObject(
  counts = fetal_mat,
  assay = "peaks",
  meta.data = combined_md
)

DefaultAssay(fetal_atac) <- "peaks"

print(fetal_atac)
print(dim(fetal_atac))
print(head(fetal_atac@meta.data))

#Step 7: Keep top accessible peaks with Signac
target_n_peaks <- min(100000, nrow(fetal_atac))
min_cutoff_q <- max(0, 100 * (1 - target_n_peaks / nrow(fetal_atac)))
min_cutoff <- paste0("q", signif(min_cutoff_q, 6))

fetal_atac <- FindTopFeatures(
  object = fetal_atac,
  assay = "peaks",
  min.cutoff = min_cutoff
)

selected_peaks <- VariableFeatures(fetal_atac)
if (length(selected_peaks) > target_n_peaks) {
  selected_peaks <- selected_peaks[seq_len(target_n_peaks)]
}

fetal_atac <- fetal_atac[selected_peaks, ]
message("Retained ", nrow(fetal_atac), " peaks using Signac::FindTopFeatures with min.cutoff=", min_cutoff)

dim(fetal_atac)

#Step 8: Keep raw peak counts only
peak_mat <- GetAssayData(fetal_atac, assay = "peaks", layer = "counts")
str(peak_mat)

print("raw counts retained in peaks counts layer")

#Step 9: Keep barebones metadata
lib_size <- Matrix::colSums(GetAssayData(fetal_atac, assay = "peaks", layer = "counts"))
fetal_atac@meta.data$library_size <- as.numeric(lib_size[colnames(fetal_atac)])
fetal_atac@meta.data <- fetal_atac@meta.data[, c("cell_type", "library_size"), drop = FALSE]
stopifnot(all(c("cell_type", "library_size") %in% colnames(fetal_atac@meta.data)))

#Step 10: Save barebones Seurat object
str(fetal_atac)
saveRDS(
  fetal_atac,
  file = "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/fetal_atac/fetal_atac_barebones_seurat.rds"
)

message("Saved barebones fetal ATAC Seurat object.")
