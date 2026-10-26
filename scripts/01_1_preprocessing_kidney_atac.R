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
library(R.utils)

print("libs done")

#Step 2: Gunzip counts file to processed directory
raw_counts <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/raw/kidney_atac/GSE211785_EXPORT_snATACseq_counts.rds.gz"
processed_counts <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/kidney_atac/GSE211785_EXPORT_snATACseq_counts.rds"

dir.create(dirname(processed_counts), recursive = TRUE, showWarnings = FALSE)

R.utils::gunzip(
  filename = raw_counts,
  destname = processed_counts,
  overwrite = TRUE,
  remove = FALSE
)

print("counts file decompressed")

#Step 3: Data
atac_counts <- readRDS(processed_counts)
str(atac_counts)

bp_dir <- "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/kidney_atac/bpcells_counts"

if (!dir.exists(bp_dir)) {
  write_matrix_dir(
    mat = atac_counts,
    dir = bp_dir
  )
}

kidney_mat <- open_matrix_dir(dir = bp_dir)

kidney_mat
dim(kidney_mat)

atac_metadata <- fread("/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/raw/kidney_atac/GSE211785_snATAC-seq_metadata.txt.gz", header = TRUE)
str(atac_metadata)

md <- as.data.frame(atac_metadata)
rownames(md) <- md$V1

md <- md[colnames(kidney_mat), , drop = FALSE]
stopifnot(identical(rownames(md), colnames(kidney_mat)))

kidney_atac <- CreateSeuratObject(
  counts = kidney_mat,
  assay = "peaks",
  meta.data = md
)

kidney_atac

DefaultAssay(kidney_atac) <- "peaks"
head(kidney_atac@meta.data)
dim(kidney_atac)

peak_gr <- StringToGRanges(rownames(kidney_atac), sep = c("-", "-"))
peaks.keep <- seqnames(peak_gr) %in% standardChromosomes(peak_gr)

kidney_atac <- kidney_atac[as.vector(peaks.keep), ]

table(seqnames(peak_gr) %in% standardChromosomes(peak_gr))
dim(kidney_atac)

#Step 4: Keep top accessible peaks with Signac
target_n_peaks <- min(100000, nrow(kidney_atac))
min_cutoff_q <- max(0, 100 * (1 - target_n_peaks / nrow(kidney_atac)))
min_cutoff <- paste0("q", signif(min_cutoff_q, 6))

kidney_atac <- FindTopFeatures(
  object = kidney_atac,
  assay = "peaks",
  min.cutoff = min_cutoff
)

selected_peaks <- VariableFeatures(kidney_atac)
if (length(selected_peaks) > target_n_peaks) {
  selected_peaks <- selected_peaks[seq_len(target_n_peaks)]
}

kidney_atac <- kidney_atac[selected_peaks, ]
message("Retained ", nrow(kidney_atac), " peaks using Signac::FindTopFeatures with min.cutoff=", min_cutoff)

dim(kidney_atac)

#Step 5: Keep raw peak counts only
peak_mat <- GetAssayData(kidney_atac, assay = "peaks", layer = "counts")
str(peak_mat)

print("raw counts retained in peaks counts layer")

#Step 6: Keep barebones metadata
lib_size <- Matrix::colSums(GetAssayData(kidney_atac, assay = "peaks", layer = "counts"))
kidney_atac@meta.data$library_size <- as.numeric(lib_size[colnames(kidney_atac)])
kidney_atac@meta.data <- kidney_atac@meta.data[, c("cell_type", "library_size"), drop = FALSE]
stopifnot(all(c("cell_type", "library_size") %in% colnames(kidney_atac@meta.data)))

#Step 7: Save barebones Seurat object
str(kidney_atac)
saveRDS(
  kidney_atac,
  file = "/home/groups/CEDAR/mishrap/scATAC_subset_selection/data/processed/kidney_atac/kidney_atac_barebones_seurat.rds"
)

message("Saved barebones Seurat object.")
