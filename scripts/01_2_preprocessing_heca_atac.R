#Step 1: libs
library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(AnnotationHub)
library(Matrix)
library(data.table)
library(zellkonverter)
library(SingleCellExperiment)
library(BPCells)

print("libs done")

#Step 2: Data

h5ad_path <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/raw/heca_atac/ATAC-all.h5ad"

bp_mat <- open_matrix_anndata_hdf5(h5ad_path)
bp_mat
dim(bp_mat)

#Step 3: Write BPCells matrix to disk if not already written

bp_dir <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/heca_atac/bpcells_counts"

if (!dir.exists(bp_dir)) {
  write_matrix_dir(
    mat = bp_mat,
    dir = bp_dir
  )
}

# always open after
heca_mat <- open_matrix_dir(dir = bp_dir)

heca_mat
dim(heca_mat)

#Step 4: Load metadata only
heca_sce <- readH5AD(
  h5ad_path,
  use_hdf5 = TRUE
)

md <- as.data.frame(colData(heca_sce))
rownames(md) <- colnames(bp_mat)

md <- md[colnames(heca_mat), , drop = FALSE]
md <- md[, "cell_type", drop = FALSE]

stopifnot(identical(rownames(md), colnames(heca_mat)))

head(md)
str(md)

#Step 5: Create object
heca_seurat <- CreateSeuratObject(counts = heca_mat, assay = "peaks", meta.data = md)
str(heca_seurat)

#Step 6: libs

DefaultAssay(heca_seurat) <- "peaks"

print("default assay set")

#Step 7: extract counts layer
peak_mat <- heca_seurat[["peaks"]]$counts
print(peak_mat)
str(peak_mat)

print("counts layer extracted")

#Step 8: Keep top accessible peaks with Signac
target_n_peaks <- min(100000, nrow(heca_seurat))
min_cutoff_q <- max(0, 100 * (1 - target_n_peaks / nrow(heca_seurat)))
min_cutoff <- paste0("q", signif(min_cutoff_q, 6))

heca_seurat <- FindTopFeatures(
  object = heca_seurat,
  assay = "peaks",
  min.cutoff = min_cutoff
)

selected_peaks <- VariableFeatures(heca_seurat)
if (length(selected_peaks) > target_n_peaks) {
  selected_peaks <- selected_peaks[seq_len(target_n_peaks)]
}

heca_seurat <- heca_seurat[selected_peaks, ]
message("Retained ", nrow(heca_seurat), " peaks using Signac::FindTopFeatures with min.cutoff=", min_cutoff)

dim(heca_seurat)

#Step 9: Keep raw peak counts only
peak_mat <- GetAssayData(heca_seurat, assay = "peaks", layer = "counts")
str(peak_mat)

print("raw counts retained in peaks counts layer")

#Step 10: keep barebones metadata
lib_size <- Matrix::colSums(GetAssayData(heca_seurat, assay = "peaks", layer = "counts"))
heca_seurat@meta.data$library_size <- as.numeric(lib_size[colnames(heca_seurat)])
heca_seurat@meta.data <- heca_seurat@meta.data[, c("cell_type", "library_size"), drop = FALSE]
stopifnot(all(c("cell_type", "library_size") %in% colnames(heca_seurat@meta.data)))

print("barebones metadata retained")

#Step 11: inspect object
str(heca_seurat)

#Step 12: save barebones Seurat object
saveRDS(
  heca_seurat,
  file = "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/heca_atac/heca_atac_barebones_seurat.rds"
)

message("Saved barebones HECA Seurat object.")

