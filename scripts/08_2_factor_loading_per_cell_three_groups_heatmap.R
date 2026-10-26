#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/08_2_factor_loading_per_cell_three_groups_heatmap.R",
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
prop_df <- prop_df[prop_df$cell_type %in% unique(cell_type), , drop = FALSE]
if (nrow(prop_df) < 3) stop("Need at least 3 overlapping cell types.")

least_ct <- prop_df$cell_type[1]
mid_ct <- prop_df$cell_type[ceiling(nrow(prop_df) / 2)]
most_ct <- prop_df$cell_type[nrow(prop_df)]
target_ct <- c(least_ct, mid_ct, most_ct)

cell_meta <- data.frame(
  cell = cells,
  cell_type = cell_type,
  stringsAsFactors = FALSE
) %>%
  filter(cell_type %in% target_ct) %>%
  mutate(cell_type = factor(cell_type, levels = target_ct))

if (nrow(cell_meta) == 0) stop("No cells found for selected cell types.")

cell_meta$mean_abs_latent <- rowMeans(abs(Z[cell_meta$cell, , drop = FALSE]))
cell_meta <- cell_meta %>% arrange(cell_type, desc(mean_abs_latent))
cell_meta$cell_order <- seq_len(nrow(cell_meta))

Z_sub <- t(Z[cell_meta$cell, , drop = FALSE])
Z_z <- t(scale(t(Z_sub)))
Z_z[is.na(Z_z)] <- 0

long_df <- data.frame(
  latent_dimension = rep(rownames(Z_z), times = ncol(Z_z)),
  cell = rep(colnames(Z_z), each = nrow(Z_z)),
  z_score = as.vector(Z_z),
  stringsAsFactors = FALSE
) %>%
  left_join(cell_meta[, c("cell", "cell_order", "cell_type")], by = "cell")

boundary_df <- cell_meta %>%
  group_by(cell_type) %>%
  summarize(end_order = max(cell_order), .groups = "drop") %>%
  mutate(boundary_x = end_order + 0.5)

label_df <- cell_meta %>%
  group_by(cell_type) %>%
  summarize(
    start = min(cell_order),
    end = max(cell_order),
    mid = (start + end) / 2,
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0(as.character(cell_type), " (n=", n_cells, ")"))

p <- ggplot(long_df, aes(x = cell_order, y = latent_dimension, fill = z_score)) +
  geom_raster() +
  geom_vline(
    data = boundary_df[-nrow(boundary_df), , drop = FALSE],
    aes(xintercept = boundary_x),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.1
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Latent value\n(row z-score)"
  ) +
  scale_x_continuous(breaks = label_df$mid, labels = label_df$label) +
  labs(
    title = paste0(atlas_name, ": Per-cell PeakVI latent heatmap"),
    subtitle = paste0("Least selected: ", least_ct, " | Intermediate: ", mid_ct, " | Most selected: ", most_ct),
    x = "Cells grouped by selected cell type",
    y = "PeakVI latent dimension"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold")
  )

ggsave(output_png, p, width = 18, height = 8, dpi = 300, units = "in", limitsize = FALSE)

summary_df <- data.frame(
  atlas = atlas_name,
  group = c("least_selected", "intermediate_selected", "most_selected"),
  cell_type = c(least_ct, mid_ct, most_ct),
  n_cells = as.integer(table(cell_meta$cell_type)[target_ct]),
  mean_enrichment_ratio = prop_df$mean_ratio[match(target_ct, prop_df$cell_type)],
  selection_rule = "cell types sorted by mean enrichment ratio across subset sizes",
  stringsAsFactors = FALSE
)
write.table(summary_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("DONE\n")
