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
      "Rscript scripts/06_2_plot_cell_type_alluvial.R",
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

if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  stop("Package 'ggalluvial' is required for this script. Please install it in your R env.")
}

df <- read.delim(input_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

required_cols <- c(
  "cell_type",
  "subset_02pct_count",
  "subset_04pct_count",
  "subset_06pct_count",
  "subset_08pct_count",
  "subset_10pct_count",
  "atlas_count"
)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns in input table: ", paste(missing_cols, collapse = ", ")))
}

plot_df <- df %>%
  select(
    cell_type,
    subset_02pct_count,
    subset_04pct_count,
    subset_06pct_count,
    subset_08pct_count,
    subset_10pct_count,
    atlas_count
  ) %>%
  pivot_longer(
    cols = -cell_type,
    names_to = "stage",
    values_to = "count"
  ) %>%
  mutate(
    stage = recode(
      stage,
      subset_02pct_count = "2%",
      subset_04pct_count = "4%",
      subset_06pct_count = "6%",
      subset_08pct_count = "8%",
      subset_10pct_count = "10%",
      atlas_count = "100%"
    ),
    stage = factor(stage, levels = c("2%", "4%", "6%", "8%", "10%", "100%"))
  )

p <- ggplot(
  plot_df,
  aes(
    x = stage,
    stratum = cell_type,
    alluvium = cell_type,
    y = count,
    fill = cell_type
  )
) +
  ggalluvial::geom_alluvium(alpha = 0.8, decreasing = FALSE) +
  ggalluvial::geom_stratum(width = 0.35, color = "grey25", size = 0.2) +
  labs(
    title = paste0(atlas_name, ": Cell-type composition flow across subset sizes"),
    subtitle = "2%, 4%, 6%, 8%, 10%, and full atlas (100%)",
    x = "Subset size",
    y = "Number of cells",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )

ggsave(
  filename = output_png,
  plot = p,
  width = 14,
  height = 8,
  dpi = 300,
  units = "in",
  limitsize = FALSE
)

cat("DONE\n")
