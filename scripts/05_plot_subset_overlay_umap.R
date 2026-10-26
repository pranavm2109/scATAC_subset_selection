#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

cat("libs loaded\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/05_plot_subset_overlay_umap.R",
      "<input_barebones_seurat_rds>",
      "<input_peakvi_umap_tsv>",
      "<subset_cells_txt>",
      "<atlas_name>",
      "<fraction_label>",
      "<output_figure_png>",
      "<color_by=cell_type|library_size>"
    )
  )
}

input_rds <- args[1]
input_umap_tsv <- args[2]
subset_cells_file <- args[3]
atlas_name <- args[4]
fraction_label <- args[5]
output_file <- args[6]
color_by <- args[7]
reduction_name <- "umap"

cat("input_rds: ", input_rds, "\n", sep = "")
cat("input_umap_tsv: ", input_umap_tsv, "\n", sep = "")
cat("subset_cells_file: ", subset_cells_file, "\n", sep = "")
cat("atlas_name: ", atlas_name, "\n", sep = "")
cat("fraction_label: ", fraction_label, "\n", sep = "")
cat("color_by: ", color_by, "\n", sep = "")
cat("output_file: ", output_file, "\n", sep = "")

if (!grepl("\\.png$", output_file, ignore.case = TRUE)) {
  stop("UMAP figures must be saved as PNG files.")
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
if (!"cell_type" %in% colnames(obj@meta.data)) {
  stop("cell_type column not found in object metadata.")
}

umap_df <- read.delim(
  input_umap_tsv,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_umap_cols <- c("cell_name", "UMAP_1", "UMAP_2")
missing_umap_cols <- setdiff(required_umap_cols, colnames(umap_df))
if (length(missing_umap_cols) > 0) {
  stop(
    "PeakVI UMAP TSV is missing columns: ",
    paste(missing_umap_cols, collapse = ", ")
  )
}
if (anyDuplicated(umap_df$cell_name)) {
  stop("PeakVI UMAP TSV contains duplicate cell names.")
}

all_cells <- colnames(obj)
missing_cells <- setdiff(all_cells, umap_df$cell_name)
extra_cells <- setdiff(umap_df$cell_name, all_cells)
if (length(missing_cells) > 0 || length(extra_cells) > 0) {
  stop(
    "Cell-name mismatch between Seurat object and PeakVI UMAP TSV: ",
    length(missing_cells), " missing from UMAP; ",
    length(extra_cells), " extra in UMAP."
  )
}

umap_df <- umap_df[match(all_cells, umap_df$cell_name), , drop = FALSE]
umap_matrix <- as.matrix(umap_df[, c("UMAP_1", "UMAP_2"), drop = FALSE])
storage.mode(umap_matrix) <- "double"
if (any(!is.finite(umap_matrix))) {
  stop("PeakVI UMAP TSV contains non-finite coordinates.")
}
rownames(umap_matrix) <- all_cells
obj[[reduction_name]] <- CreateDimReducObject(
  embeddings = umap_matrix,
  key = "UMAP_",
  assay = DefaultAssay(obj)
)

subset_cells <- readLines(subset_cells_file, warn = FALSE)
subset_cells <- unique(subset_cells[nzchar(subset_cells)])
unmatched_subset_cells <- setdiff(subset_cells, all_cells)
if (length(unmatched_subset_cells) > 0) {
  stop(
    length(unmatched_subset_cells),
    " selected subset cell IDs do not match the Seurat object."
  )
}
if (length(subset_cells) == 0) {
  stop("No subset cell IDs were supplied.")
}

overlay_df <- as.data.frame(umap_matrix[subset_cells, , drop = FALSE])
colnames(overlay_df) <- c("UMAP_1", "UMAP_2")

if (color_by == "cell_type") {
  base_plot <- DimPlot(
    object = obj,
    reduction = reduction_name,
    group.by = "cell_type",
    pt.size = 0.12,
    raster = FALSE,
    label = TRUE,
    repel = TRUE,
    label.size = 4.6,
    label.box = TRUE
  )
} else if (color_by == "library_size") {
  if (!"library_size" %in% colnames(obj@meta.data)) {
    stop("library_size column not found in object metadata.")
  }
  if (!is.numeric(obj@meta.data$library_size)) {
    stop("library_size metadata column exists but is not numeric.")
  }
  lib_vals <- obj@meta.data$library_size
  if (all(is.na(lib_vals))) {
    stop("library_size metadata column is all NA.")
  }
  lib_sd <- sd(lib_vals, na.rm = TRUE)
  if (is.na(lib_sd) || lib_sd == 0) {
    stop("library_size has zero variance; cannot color by library size.")
  }

  base_plot <- FeaturePlot(
    object = obj,
    reduction = reduction_name,
    features = "library_size",
    pt.size = 0.12,
    raster = FALSE,
    min.cutoff = "q01",
    max.cutoff = "q99",
    cols = c("#f0f0f0", "#74add1", "#4575b4", "#313695")
  )
} else {
  stop("color_by must be 'cell_type' or 'library_size'")
}

p <- base_plot +
  geom_point(
    data = overlay_df,
    mapping = aes(x = UMAP_1, y = UMAP_2, shape = "Selected subset cells"),
    inherit.aes = FALSE,
    color = "black",
    size = 0.95,
    alpha = 1
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Selected subset cells" = 16)
  ) +
  labs(
    title = paste0(atlas_name, " | PeakVI UMAP | subset ", fraction_label),
    subtitle = paste0("Highlighted subset cells: ", length(subset_cells)),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 13),
    plot.margin = margin(10, 20, 10, 10),
    legend.position = "right",
    legend.text = element_text(size = 14)
  ) +
  coord_cartesian(clip = "off")

if (color_by == "cell_type") {
  p <- p +
    guides(
      color = "none",
      shape = guide_legend(override.aes = list(size = 3.2, color = "black"))
    )
} else {
  p <- p +
    labs(color = "Library size") +
    guides(
      shape = guide_legend(override.aes = list(size = 3.2, color = "black")),
      colour = guide_colorbar(order = 1),
      color = guide_colorbar(order = 1)
    ) +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
}

ggsave(
  filename = output_file,
  plot = p,
  width = 20,
  height = 16,
  dpi = 400,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)

cat("figure saved\n")
cat("DONE\n")
