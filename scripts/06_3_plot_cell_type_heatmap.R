#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/06_3_plot_cell_type_heatmap.R",
      "<input_cell_type_proportions_tsv>",
      "<atlas_name>",
      "<output_png>"
    )
  )
}

input_tsv <- args[1]
atlas_name <- args[2]
output_png <- args[3]

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)

df <- read.delim(input_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

required_cols <- c(
  "cell_type",
  "ratio_02pct_vs_atlas",
  "ratio_04pct_vs_atlas",
  "ratio_06pct_vs_atlas",
  "ratio_08pct_vs_atlas",
  "ratio_10pct_vs_atlas"
)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns in input table: ", paste(missing_cols, collapse = ", ")))
}

plot_df <- df %>%
  select(all_of(required_cols)) %>%
  pivot_longer(
    cols = starts_with("ratio_"),
    names_to = "subset",
    values_to = "enrichment_score"
  ) %>%
  mutate(
    subset = recode(
      subset,
      ratio_02pct_vs_atlas = "2%",
      ratio_04pct_vs_atlas = "4%",
      ratio_06pct_vs_atlas = "6%",
      ratio_08pct_vs_atlas = "8%",
      ratio_10pct_vs_atlas = "10%"
    ),
    subset = factor(subset, levels = c("2%", "4%", "6%", "8%", "10%"))
  )

p <- ggplot(plot_df, aes(x = cell_type, y = subset, fill = enrichment_score)) +
  geom_tile(color = "grey90", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 1,
    name = "Enrichment score\n(subset/atlas)"
  ) +
  labs(
    title = paste0(atlas_name, ": Cell-type enrichment heatmap"),
    subtitle = "Rows are subset sizes; values are subset-to-atlas enrichment ratios",
    x = "Cell type",
    y = "Subset size"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 11),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = output_png,
  plot = p,
  width = 16,
  height = 5,
  dpi = 300,
  units = "in",
  limitsize = FALSE
)

cat("DONE\n")
