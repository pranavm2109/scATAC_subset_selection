#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/06_4_plot_cell_type_fraction_heatmap.R",
      "<input_tidy_fraction_tsv>",
      "<atlas_name>",
      "<output_png>"
    )
  )
}

input_tidy_tsv <- args[1]
atlas_name <- args[2]
output_png <- args[3]

subset_levels <- c("02%", "04%", "06%", "08%", "10%", "100%")

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)

df <- read.delim(input_tidy_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

required_cols <- c("atlas", "subset_size", "cell_type", "n_cells", "fraction")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns in tidy fractions TSV: ", paste(missing_cols, collapse = ", ")))
}

unexpected_subset_sizes <- setdiff(unique(df$subset_size), subset_levels)
if (length(unexpected_subset_sizes) > 0) {
  stop(paste0("Unexpected subset_size values: ", paste(unexpected_subset_sizes, collapse = ", ")))
}

df$subset_size <- factor(df$subset_size, levels = subset_levels)

fraction_sums <- aggregate(fraction ~ subset_size, data = df, sum)
bad_sums <- fraction_sums[abs(fraction_sums$fraction - 1) > 1e-6, , drop = FALSE]
if (nrow(bad_sums) > 0) {
  stop(
    paste0(
      "Fractions do not sum to approximately 1 for subset_size values: ",
      paste(
        paste0(as.character(bad_sums$subset_size), " (sum=", signif(bad_sums$fraction, 12), ")"),
        collapse = ", "
      )
    )
  )
}

full_df <- df[df$subset_size == "100%", c("cell_type", "fraction"), drop = FALSE]
full_df <- full_df[order(full_df$fraction, decreasing = TRUE), , drop = FALSE]
df$cell_type <- factor(df$cell_type, levels = full_df$cell_type)

n_cell_types <- length(levels(df$cell_type))
plot_width <- max(12, 0.18 * n_cell_types)

p <- ggplot(df, aes(x = cell_type, y = subset_size, fill = fraction)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08306B",
    name = "Fraction\nof cells",
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  labs(
    title = paste0(atlas_name, ": cell-type fractions across selected subsets"),
    subtitle = "Fill shows fraction of cells in each subset; 100% is the full atlas, not an enrichment ratio",
    x = "Cell type",
    y = "Subset size"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = output_png,
  plot = p,
  width = plot_width,
  height = 5.5,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

cat("DONE\n")
