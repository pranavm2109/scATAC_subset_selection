#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/08_1_factor_loading_all_cell_types_heatmap.R",
      "<input_seurat_rds>",
      "<input_peakvi_latent_tsv>",
      "<input_celltype_proportions_tsv>",
      "<atlas_name>",
      "<output_heatmap_png>",
      "<output_summary_tsv>"
    )
  )
}

input_rds <- args[1]
input_latent_tsv <- args[2]
input_prop_tsv <- args[3]
atlas_name <- args[4]
output_png <- args[5]
output_tsv <- args[6]

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
if (!"cell_type" %in% colnames(obj@meta.data)) stop("cell_type column missing in Seurat metadata.")
cells <- colnames(obj)
cell_type <- as.character(obj@meta.data[cells, "cell_type", drop = TRUE])

latent_df <- read.delim(input_latent_tsv, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(latent_df) < 2) stop("PeakVI latent TSV must contain cell_name plus latent dimensions.")
latent_cells <- as.character(latent_df[[1]])
if (anyDuplicated(latent_cells)) stop("Duplicate cell names found in PeakVI latent TSV.")
rownames(latent_df) <- latent_cells
missing_cells <- setdiff(cells, latent_cells)
if (length(missing_cells) > 0) stop(paste0("PeakVI latent TSV missing ", length(missing_cells), " Seurat cells."))

Z <- as.matrix(latent_df[cells, -1, drop = FALSE])
storage.mode(Z) <- "double"
if (any(!is.finite(Z))) stop("PeakVI latent matrix contains non-finite values.")
colnames(Z) <- paste0("PeakVI_", seq_len(ncol(Z)))

prop_df <- read.delim(input_prop_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ratio_cols <- grep("^ratio_.*_vs_atlas$", colnames(prop_df), value = TRUE)
if (length(ratio_cols) == 0) stop("No ratio_*_vs_atlas columns found in proportions TSV.")
prop_df$mean_ratio <- rowMeans(prop_df[, ratio_cols, drop = FALSE], na.rm = TRUE)
prop_df <- prop_df[order(prop_df$mean_ratio), , drop = FALSE]

ct_levels <- prop_df$cell_type
ct_levels <- ct_levels[ct_levels %in% unique(cell_type)]
if (length(ct_levels) == 0) stop("No overlapping cell types between proportions table and Seurat metadata.")

mean_mat <- matrix(NA_real_, nrow = length(ct_levels), ncol = ncol(Z))
rownames(mean_mat) <- ct_levels
colnames(mean_mat) <- colnames(Z)

for (i in seq_along(ct_levels)) {
  ct_i <- ct_levels[i]
  cells_i <- cells[cell_type == ct_i]
  if (length(cells_i) > 0) {
    mean_mat[i, ] <- colMeans(Z[cells_i, , drop = FALSE])
  }
}

z_mat <- t(scale(t(mean_mat)))
z_mat[is.na(z_mat)] <- 0

long_df <- data.frame(
  cell_type = rep(rownames(z_mat), times = ncol(z_mat)),
  latent_dimension = rep(colnames(z_mat), each = nrow(z_mat)),
  z_score = as.vector(z_mat),
  stringsAsFactors = FALSE
)

p <- ggplot(long_df, aes(x = latent_dimension, y = cell_type, fill = z_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = paste0(atlas_name, ": PeakVI latent heatmap (all cell types)"),
    subtitle = "Cell types ordered by mean subset enrichment ratio; rows are z-scored",
    x = "PeakVI latent dimension",
    y = "Cell type",
    fill = "Row z-score"
  )

plot_height <- max(6, 0.25 * length(ct_levels))
ggsave(output_png, p, width = 12, height = plot_height, dpi = 300, units = "in", limitsize = FALSE)

summary_df <- data.frame(
  atlas = atlas_name,
  cell_type = ct_levels,
  mean_enrichment_ratio = prop_df$mean_ratio[match(ct_levels, prop_df$cell_type)],
  stringsAsFactors = FALSE
)
write.table(summary_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("DONE\n")
