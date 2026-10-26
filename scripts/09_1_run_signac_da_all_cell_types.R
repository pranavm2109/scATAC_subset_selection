#!/usr/bin/env Rscript

# Run annotation-blind Signac/BPCells differential accessibility for every
# annotated cell type versus the rest. Labels are used only for evaluation
# after subsets have already been selected.

suppressPackageStartupMessages({
  library(BPCells)
  library(Signac)
  library(Seurat)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8L || length(args) > 13L) {
  stop(
    paste(
      "Usage: Rscript scripts/09_1_run_signac_da_all_cell_types.R",
      "<input_seurat.rds> <selected_cells.txt|FULL> <atlas> <dataset_size>",
      "<output_da_peaks.tsv> <output_status.tsv>",
      "<output_reference_effects.tsv> <full_reference_da.tsv|SELF>",
      "[min_target_cells=10] [min_rest_cells=20] [fdr_threshold=0.05]",
      "[abs_log2fc_threshold=0.25] [min_pct=0.01]"
    )
  )
}

input_rds <- args[[1]]
selected_cells_path <- args[[2]]
atlas_name <- args[[3]]
dataset_size <- args[[4]]
output_da_tsv <- args[[5]]
output_status_tsv <- args[[6]]
output_reference_tsv <- args[[7]]
full_reference_path <- args[[8]]
min_target_cells <- if (length(args) >= 9L) as.integer(args[[9]]) else 10L
min_rest_cells <- if (length(args) >= 10L) as.integer(args[[10]]) else 20L
fdr_threshold <- if (length(args) >= 11L) as.numeric(args[[11]]) else 0.05
abs_log2fc_threshold <- if (length(args) >= 12L) as.numeric(args[[12]]) else 0.25
min_pct <- if (length(args) >= 13L) as.numeric(args[[13]]) else 0.01

allowed_sizes <- c("02pct", "04pct", "06pct", "08pct", "10pct", "100pct")
if (!dataset_size %in% allowed_sizes) {
  stop("dataset_size must be one of: ", paste(allowed_sizes, collapse = ", "))
}
if (min_target_cells < 2L || min_rest_cells < 2L) {
  stop("min_target_cells and min_rest_cells must both be at least 2")
}
if (!is.finite(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold > 1) {
  stop("fdr_threshold must be greater than 0 and no greater than 1")
}
if (!is.finite(abs_log2fc_threshold) || abs_log2fc_threshold < 0) {
  stop("abs_log2fc_threshold must be nonnegative")
}
if (!is.finite(min_pct) || min_pct < 0 || min_pct > 1) {
  stop("min_pct must be between 0 and 1")
}

empty_da_table <- function() {
  tibble(
    atlas = character(),
    dataset_size = character(),
    cell_type = character(),
    peak = character(),
    avg_log2FC = double(),
    pct_target = double(),
    pct_rest = double(),
    p_value = double(),
    adjusted_p_value = double(),
    direction = character(),
    n_target_cells = integer(),
    n_rest_cells = integer()
  )
}

empty_reference_table <- function() {
  tibble(
    atlas = character(),
    dataset_size = character(),
    cell_type = character(),
    peak = character(),
    full_avg_log2FC = double(),
    subset_avg_log2FC = double(),
    full_direction = character(),
    subset_direction = character(),
    subset_p_value = double(),
    subset_adjusted_p_value = double(),
    subset_pct_target = double(),
    subset_pct_rest = double(),
    subset_is_significant = logical(),
    effect_available = logical(),
    same_direction_effect = logical(),
    direction_switched_effect = logical(),
    strict_same_direction_capture = logical(),
    n_target_cells = integer(),
    n_rest_cells = integer(),
    subset_analysis_status = character()
  )
}

cat("Reading Seurat object:", input_rds, "\n")
atlas <- readRDS(input_rds)
if (!inherits(atlas, "Seurat")) {
  stop("Input RDS is not a Seurat object")
}
if (!"peaks" %in% Assays(atlas)) {
  stop("Input Seurat object does not contain a 'peaks' assay")
}

required_metadata <- c("cell_type", "library_size")
missing_metadata <- setdiff(required_metadata, colnames(atlas[[]]))
if (length(missing_metadata) > 0L) {
  stop(
    "Input Seurat object is missing required metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

metadata <- atlas[[]]
if (anyNA(metadata$cell_type) ||
    any(trimws(as.character(metadata$cell_type)) == "")) {
  stop("cell_type contains missing or empty values")
}
if (anyNA(metadata$library_size) ||
    any(!is.finite(as.numeric(metadata$library_size))) ||
    any(as.numeric(metadata$library_size) < 0)) {
  stop("library_size must contain finite, nonnegative values")
}

full_cell_types <- sort(unique(as.character(metadata$cell_type)))
if (length(full_cell_types) < 2L) {
  stop("At least two cell types are required for target-versus-rest DA")
}

is_full_atlas <- identical(toupper(selected_cells_path), "FULL")
if (dataset_size == "100pct" && !is_full_atlas) {
  stop("The 100pct dataset must use FULL as the selected-cells argument")
}
if (dataset_size != "100pct" && is_full_atlas) {
  stop("Subset datasets must provide a selected-cells file")
}
if (dataset_size == "100pct" &&
    !identical(toupper(full_reference_path), "SELF")) {
  stop("The 100pct dataset must use SELF as the reference argument")
}
if (dataset_size != "100pct" &&
    identical(toupper(full_reference_path), "SELF")) {
  stop("Subset datasets must provide the full-atlas DA reference table")
}

if (!is_full_atlas) {
  selected_cells <- scan(selected_cells_path, what = character(), quiet = TRUE)
  if (length(selected_cells) == 0L) {
    stop("Selected-cells file is empty: ", selected_cells_path)
  }
  if (anyDuplicated(selected_cells)) {
    stop("Selected-cells file contains duplicate cell names: ", selected_cells_path)
  }
  missing_cells <- setdiff(selected_cells, colnames(atlas))
  if (length(missing_cells) > 0L) {
    stop(
      length(missing_cells),
      " selected cells are absent from the Seurat object. First examples: ",
      paste(head(missing_cells, 10L), collapse = ", ")
    )
  }
  atlas <- subset(atlas, cells = selected_cells)
}

if (ncol(atlas) < min_target_cells + min_rest_cells) {
  stop("Dataset is too small for the requested target/rest minimums")
}

DefaultAssay(atlas) <- "peaks"
Idents(atlas) <- factor(as.character(atlas$cell_type), levels = full_cell_types)
cell_groups <- as.character(Idents(atlas))
cell_type_counts <- tibble(
  cell_type = full_cell_types,
  n_target_cells = as.integer(vapply(
    full_cell_types,
    function(cell_type_name) sum(cell_groups == cell_type_name),
    integer(1)
  ))
) %>%
  mutate(n_rest_cells = as.integer(ncol(atlas) - n_target_cells))

eligible_cell_types <- cell_type_counts %>%
  filter(
    n_target_cells >= min_target_cells,
    n_rest_cells >= min_rest_cells
  ) %>%
  pull(cell_type)

cat(
  "Running DA for", atlas_name, dataset_size, "with",
  format(ncol(atlas), big.mark = ","), "cells,",
  format(nrow(atlas), big.mark = ","), "peaks, and",
  length(eligible_cell_types), "eligible cell types\n"
)

marker_table_all <- tibble(
  cell_type = character(),
  peak = character(),
  p_value = double(),
  foreground_mean = double(),
  background_mean = double(),
  pct_target = double(),
  pct_rest = double(),
  n_target_cells = integer(),
  n_rest_cells = integer(),
  avg_log2FC = double(),
  adjusted_p_value = double()
)
if (length(eligible_cell_types) > 0L) {
  counts_matrix <- LayerData(atlas[["peaks"]], layer = "counts")
  if (!inherits(counts_matrix, "IterableMatrix")) {
    stop("The peaks count layer must be a BPCells IterableMatrix")
  }

  da_tmp_root <- Sys.getenv(
    "SIGNAC_DA_TMPDIR",
    unset = Sys.getenv("TMPDIR", unset = tempdir())
  )
  if (!nzchar(da_tmp_root)) {
    da_tmp_root <- tempdir()
  }
  dir.create(da_tmp_root, recursive = TRUE, showWarnings = FALSE)
  row_major_dir <- file.path(
    da_tmp_root,
    paste0("signac_da_row_major_", atlas_name, "_", dataset_size, "_", Sys.getpid())
  )
  if (dir.exists(row_major_dir)) {
    unlink(row_major_dir, recursive = TRUE, force = TRUE)
  }
  on.exit(unlink(row_major_dir, recursive = TRUE, force = TRUE), add = TRUE)

  cat("Creating one row-major BPCells matrix at:", row_major_dir, "\n")
  row_major_counts <- transpose_storage_order(
    counts_matrix,
    outdir = row_major_dir,
    tmpdir = da_tmp_root
  )

  cat("Running one-pass BPCells Wilcoxon tests and detection fractions\n")
  marker_raw <- marker_features(
    row_major_counts,
    groups = cell_groups,
    method = "wilcoxon"
  )
  marker_detected <- marker_features(
    row_major_counts > 0,
    groups = cell_groups,
    method = "wilcoxon"
  )

  marker_table_all <- marker_raw %>%
    transmute(
      cell_type = as.character(foreground),
      peak = as.character(feature),
      p_value = as.numeric(p_val_raw),
      foreground_mean = as.numeric(foreground_mean),
      background_mean = as.numeric(background_mean)
    ) %>%
    filter(cell_type %in% eligible_cell_types) %>%
    left_join(
      marker_detected %>%
        transmute(
          cell_type = as.character(foreground),
          peak = as.character(feature),
          pct_target = round(as.numeric(foreground_mean), 6),
          pct_rest = round(as.numeric(background_mean), 6)
        ) %>%
        filter(cell_type %in% eligible_cell_types),
      by = c("cell_type", "peak")
    ) %>%
    left_join(cell_type_counts, by = "cell_type") %>%
    group_by(cell_type) %>%
    mutate(
      avg_log2FC = log2(
        (foreground_mean * n_target_cells + 1) / n_target_cells
      ) - log2(
        (background_mean * n_rest_cells + 1) / n_rest_cells
      ),
      adjusted_p_value = p.adjust(p_value, method = "BH")
    ) %>%
    ungroup()

  rm(counts_matrix, row_major_counts, marker_raw, marker_detected)
  invisible(gc(verbose = FALSE))
  unlink(row_major_dir, recursive = TRUE, force = TRUE)
}

marker_table_tested <- marker_table_all %>%
  filter(pmax(pct_target, pct_rest) >= min_pct)

significant_table <- marker_table_tested %>%
  filter(
    is.finite(avg_log2FC),
    is.finite(adjusted_p_value),
    adjusted_p_value <= fdr_threshold,
    abs(avg_log2FC) >= abs_log2fc_threshold
  )

da_table <- significant_table %>%
  transmute(
    atlas = .env$atlas_name,
    dataset_size = .env$dataset_size,
    cell_type,
    peak,
    avg_log2FC,
    pct_target,
    pct_rest,
    p_value,
    adjusted_p_value,
    direction = if_else(avg_log2FC > 0, "open", "closed"),
    n_target_cells = as.integer(n_target_cells),
    n_rest_cells = as.integer(n_rest_cells)
  )

tested_summary <- marker_table_tested %>%
  count(cell_type, name = "n_tested_peaks")
significant_summary <- significant_table %>%
  group_by(cell_type) %>%
  summarise(
    n_significant_peaks = n(),
    n_open_peaks = sum(avg_log2FC > 0),
    n_closed_peaks = sum(avg_log2FC < 0),
    .groups = "drop"
  )

# Build from the complete count table so every full-atlas cell type has one
# status row, including types that are absent or underpowered in a subset.
status_table <- cell_type_counts %>%
  left_join(tested_summary, by = "cell_type") %>%
  left_join(significant_summary, by = "cell_type") %>%
  mutate(
    atlas = .env$atlas_name,
    dataset_size = .env$dataset_size,
    n_dataset_cells = as.integer(ncol(.env$atlas)),
    n_tested_peaks = coalesce(as.integer(n_tested_peaks), 0L),
    n_significant_peaks = coalesce(as.integer(n_significant_peaks), 0L),
    n_open_peaks = coalesce(as.integer(n_open_peaks), 0L),
    n_closed_peaks = coalesce(as.integer(n_closed_peaks), 0L),
    status = case_when(
      n_target_cells < min_target_cells ~ "insufficient_target_cells",
      n_rest_cells < min_rest_cells ~ "insufficient_rest_cells",
      TRUE ~ "success"
    ),
    error_message = case_when(
      n_target_cells < min_target_cells ~ paste0(
        "Requires at least ", min_target_cells, " target cells"
      ),
      n_rest_cells < min_rest_cells ~ paste0(
        "Requires at least ", min_rest_cells, " rest cells"
      ),
      TRUE ~ ""
    ),
    test = "wilcox",
    multiple_testing = "BH_within_cell_type",
    latent_variable = "none",
    fdr_threshold = .env$fdr_threshold,
    abs_log2fc_threshold = .env$abs_log2fc_threshold,
    min_pct = .env$min_pct
  ) %>%
  select(
    atlas, dataset_size, cell_type, n_dataset_cells,
    n_target_cells, n_rest_cells, n_tested_peaks,
    n_significant_peaks, n_open_peaks, n_closed_peaks,
    status, error_message, test, multiple_testing, latent_variable,
    fdr_threshold, abs_log2fc_threshold, min_pct
  )

if (nrow(status_table) != length(full_cell_types) ||
    !setequal(status_table$cell_type, full_cell_types)) {
  stop(
    "Internal error while constructing the per-cell-type status table: ",
    "expected ", length(full_cell_types), " unique cell types, produced ",
    nrow(status_table), " rows for ",
    n_distinct(status_table$cell_type), " unique cell types"
  )
}

if (dataset_size == "100pct") {
  reference_effects <- if (nrow(da_table) > 0L) {
    da_table %>%
      transmute(
        atlas,
        dataset_size,
        cell_type,
        peak,
        full_avg_log2FC = avg_log2FC,
        subset_avg_log2FC = avg_log2FC,
        full_direction = direction,
        subset_direction = direction,
        subset_p_value = p_value,
        subset_adjusted_p_value = adjusted_p_value,
        subset_pct_target = pct_target,
        subset_pct_rest = pct_rest,
        subset_is_significant = TRUE,
        effect_available = TRUE,
        same_direction_effect = TRUE,
        direction_switched_effect = FALSE,
        strict_same_direction_capture = TRUE,
        n_target_cells,
        n_rest_cells,
        subset_analysis_status = "success"
      )
  } else {
    empty_reference_table()
  }
} else {
  full_reference <- read_tsv(full_reference_path, show_col_types = FALSE)
  required_reference_columns <- c(
    "atlas", "dataset_size", "cell_type", "peak", "avg_log2FC", "direction"
  )
  missing_reference_columns <- setdiff(required_reference_columns, names(full_reference))
  if (length(missing_reference_columns) > 0L) {
    stop(
      "Full-atlas reference table is missing columns: ",
      paste(missing_reference_columns, collapse = ", ")
    )
  }
  if (nrow(full_reference) > 0L &&
      (any(full_reference$atlas != atlas_name) ||
       any(full_reference$dataset_size != "100pct"))) {
    stop("Full-atlas DA reference does not match the requested atlas")
  }
  if (anyDuplicated(full_reference[c("cell_type", "peak")])) {
    stop("Full-atlas DA reference contains duplicate cell_type/peak pairs")
  }

  subset_effects <- marker_table_all %>%
    transmute(
      cell_type,
      peak,
      subset_avg_log2FC = avg_log2FC,
      subset_p_value = p_value,
      subset_adjusted_p_value = adjusted_p_value,
      subset_pct_target = pct_target,
      subset_pct_rest = pct_rest
    )
  status_join <- status_table %>%
    select(
      cell_type,
      n_target_cells,
      n_rest_cells,
      subset_analysis_status = status
    )

  reference_effects <- full_reference %>%
    transmute(
      atlas = atlas_name,
      dataset_size = .env$dataset_size,
      cell_type = as.character(cell_type),
      peak = as.character(peak),
      full_avg_log2FC = as.numeric(avg_log2FC),
      full_direction = as.character(direction)
    ) %>%
    left_join(subset_effects, by = c("cell_type", "peak")) %>%
    left_join(status_join, by = "cell_type") %>%
    mutate(
      effect_available = (
        subset_analysis_status == "success" &
          is.finite(subset_avg_log2FC)
      ),
      subset_direction = case_when(
        !effect_available ~ NA_character_,
        subset_avg_log2FC > 0 ~ "open",
        subset_avg_log2FC < 0 ~ "closed",
        TRUE ~ "zero"
      ),
      subset_is_significant = (
        effect_available &
          pmax(subset_pct_target, subset_pct_rest) >= min_pct &
          is.finite(subset_adjusted_p_value) &
          subset_adjusted_p_value <= fdr_threshold &
          abs(subset_avg_log2FC) >= abs_log2fc_threshold
      ),
      same_direction_effect = (
        effect_available &
          subset_direction %in% c("open", "closed") &
          subset_direction == full_direction
      ),
      direction_switched_effect = (
        effect_available &
          subset_direction %in% c("open", "closed") &
          subset_direction != full_direction
      ),
      strict_same_direction_capture = (
        subset_is_significant & same_direction_effect
      )
    ) %>%
    select(
      atlas, dataset_size, cell_type, peak,
      full_avg_log2FC, subset_avg_log2FC,
      full_direction, subset_direction,
      subset_p_value, subset_adjusted_p_value,
      subset_pct_target, subset_pct_rest,
      subset_is_significant, effect_available,
      same_direction_effect, direction_switched_effect,
      strict_same_direction_capture,
      n_target_cells, n_rest_cells, subset_analysis_status
    )
}

dir.create(dirname(output_da_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_status_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_reference_tsv), recursive = TRUE, showWarnings = FALSE)
write_tsv(da_table, output_da_tsv)
write_tsv(status_table, output_status_tsv)
write_tsv(reference_effects, output_reference_tsv)

cat("Wrote significant DA peaks:", output_da_tsv, "\n")
cat("Wrote per-cell-type status:", output_status_tsv, "\n")
cat("Wrote atlas-reference effects:", output_reference_tsv, "\n")
cat(
  "Successful cell types:", sum(status_table$status == "success"),
  "| insufficient:", sum(grepl("^insufficient_", status_table$status)), "\n"
)
cat("Multiple-testing correction: BH FDR within each cell type\n")
cat("DONE\n")
