#!/usr/bin/env Rscript

# LEGACY DIAGNOSTIC: this script summarizes the former shift-min raw PeakVI
# feature representation. It is retained for historical reproducibility and is
# intentionally excluded from the active landmark-feature production workflow.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/07_peakvi_feature_coverage_vs_subset_enrichment.R",
      "<input_seurat_rds>",
      "<input_peakvi_latent_tsv>",
      "<input_proportion_tsv>",
      "<atlas_name>",
      "<output_table_tsv>",
      "<output_pearson_png>",
      "<output_spearman_png>",
      "[log_transform_ratio=true|false]"
    )
  )
}

input_rds <- args[1]
input_latent_tsv <- args[2]
input_prop_tsv <- args[3]
atlas_name <- args[4]
output_table_tsv <- args[5]
output_pearson_png <- args[6]
output_spearman_png <- args[7]
log_transform_ratio <- if (length(args) >= 8) tolower(args[8]) else "true"

if (!(log_transform_ratio %in% c("true", "false"))) {
  stop("log_transform_ratio must be 'true' or 'false'")
}
use_log_transform <- identical(log_transform_ratio, "true")

make_nonnegative_shift_min <- function(mat) {
  mins <- apply(mat, 2, min, na.rm = TRUE)
  sweep(mat, 2, mins, FUN = "-")
}

dir.create(dirname(output_table_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_pearson_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_spearman_png), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
if (!"cell_type" %in% colnames(obj@meta.data)) stop("cell_type not found in metadata.")

cells <- colnames(obj)
cell_type <- as.character(obj@meta.data[cells, "cell_type", drop = TRUE])

latent_df <- read.delim(input_latent_tsv, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(latent_df) < 2) stop("PeakVI latent TSV must contain cell_name plus latent dimensions.")
if (colnames(latent_df)[1] != "cell_name") {
  warning("First latent TSV column is not named cell_name; treating it as cell names.")
}
latent_cells <- as.character(latent_df[[1]])
if (anyDuplicated(latent_cells)) stop("Duplicate cell names found in PeakVI latent TSV.")
rownames(latent_df) <- latent_cells
missing_cells <- setdiff(cells, latent_cells)
if (length(missing_cells) > 0) stop(paste0("PeakVI latent TSV missing ", length(missing_cells), " Seurat cells."))

latent_mat <- as.matrix(latent_df[cells, -1, drop = FALSE])
storage.mode(latent_mat) <- "double"
if (any(!is.finite(latent_mat))) stop("PeakVI latent matrix contains non-finite values.")

# Match the subset-selection feature space: PeakVI latents are shifted to be nonnegative.
latent_nonnegative <- make_nonnegative_shift_min(latent_mat)

ct <- factor(cell_type)
G <- sparse.model.matrix(~ 0 + ct)
type_levels <- levels(ct)
colnames(G) <- type_levels
rownames(G) <- cells

sum_by_type <- t(G) %*% latent_nonnegative
n_by_type <- as.numeric(colSums(G))
names(n_by_type) <- type_levels

mean_by_type <- sweep(as.matrix(sum_by_type), 1, n_by_type, FUN = "/")
rownames(mean_by_type) <- type_levels
coverage_score <- rowSums(sqrt(pmax(mean_by_type, 0)))

coverage_df <- data.frame(
  cell_type = type_levels,
  peakvi_coverage_score = as.numeric(coverage_score),
  n_cells = n_by_type,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(peakvi_coverage_score))

prop_df <- read.delim(input_prop_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ratio_cols <- grep("^ratio_.*_vs_atlas$", colnames(prop_df), value = TRUE)
if (length(ratio_cols) == 0) stop("No ratio_*_vs_atlas columns found in proportion table.")

df <- prop_df %>% left_join(coverage_df, by = "cell_type")
write.table(df, file = output_table_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

df_long <- df %>%
  select(cell_type, peakvi_coverage_score, all_of(ratio_cols)) %>%
  pivot_longer(cols = all_of(ratio_cols), names_to = "subset_col", values_to = "enrichment_ratio") %>%
  mutate(frac = sub("^ratio_(.*)_vs_atlas$", "\\1", subset_col))

if (use_log_transform) {
  df_long <- df_long %>% mutate(enrichment_for_plot = log(enrichment_ratio + 1e-8))
  y_axis_label <- "log(subset proportion / atlas proportion)"
  title_mode <- "cell-type log-enrichment vs PeakVI latent coverage score"
} else {
  df_long <- df_long %>% mutate(enrichment_for_plot = enrichment_ratio)
  y_axis_label <- "subset proportion / atlas proportion"
  title_mode <- "cell-type enrichment ratio vs PeakVI latent coverage score"
}

pearson_df <- df_long %>%
  group_by(frac) %>%
  summarize(r = cor(peakvi_coverage_score, enrichment_for_plot, use = "complete.obs", method = "pearson"), .groups = "drop") %>%
  mutate(legend_label = sprintf("%s (r = %.2f)", frac, r))

spearman_df <- df_long %>%
  group_by(frac) %>%
  summarize(rho = cor(peakvi_coverage_score, enrichment_for_plot, use = "complete.obs", method = "spearman"), .groups = "drop") %>%
  mutate(legend_label = sprintf("%s (rho = %.2f)", frac, rho))

df_long_pearson <- df_long %>% left_join(pearson_df, by = "frac")
df_long_spearman <- df_long %>% left_join(spearman_df, by = "frac")

coverage_x_label <- "Latent coverage score (sum sqrt of cell-type mean shifted PeakVI latent dimensions)"
plot_subtitle_common <- "Raw cell-type scatter; positive correlations indicate monotone enrichment with latent coverage"

p_pearson <- ggplot(df_long_pearson, aes(x = peakvi_coverage_score, y = enrichment_for_plot, color = legend_label)) +
  geom_point(alpha = 0.8, size = 2.2) +
  labs(
    title = paste0(atlas_name, ": ", title_mode),
    subtitle = paste0(plot_subtitle_common, "; Pearson r shown in legend"),
    x = coverage_x_label,
    y = y_axis_label,
    color = "Subset"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

p_spearman <- ggplot(df_long_spearman, aes(x = peakvi_coverage_score, y = enrichment_for_plot, color = legend_label)) +
  geom_point(alpha = 0.8, size = 2.2) +
  labs(
    title = paste0(atlas_name, ": ", title_mode),
    subtitle = paste0(plot_subtitle_common, "; Spearman rho shown in legend"),
    x = coverage_x_label,
    y = y_axis_label,
    color = "Subset"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(output_pearson_png, p_pearson, width = 13.33, height = 7.5, units = "in", dpi = 300)
ggsave(output_spearman_png, p_spearman, width = 13.33, height = 7.5, units = "in", dpi = 300)

cat("DONE\n")
