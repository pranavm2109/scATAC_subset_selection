#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
subset_sizes <- c("02pct", "04pct", "06pct", "08pct", "10pct")

if (length(args) != 20L) {
  stop(
    paste(
      "Usage: Rscript scripts/09_2_compute_signac_da_peak_capture.R",
      "<atlas> <full_da.tsv> <full_status.tsv> <output_summary.tsv>",
      "<output_details.tsv> <five_subset_da_tsvs>",
      "<five_subset_status_tsvs> <five_reference_effect_tsvs>"
    )
  )
}

atlas_name <- args[[1]]
full_da_path <- args[[2]]
full_status_path <- args[[3]]
output_summary <- args[[4]]
output_details <- args[[5]]
subset_da_paths <- args[6:10]
subset_status_paths <- args[11:15]
reference_effect_paths <- args[16:20]

required_da_columns <- c(
  "atlas", "dataset_size", "cell_type", "peak", "avg_log2FC",
  "adjusted_p_value", "direction", "n_target_cells", "n_rest_cells"
)
required_status_columns <- c(
  "atlas", "dataset_size", "cell_type", "n_target_cells",
  "n_rest_cells", "status"
)
required_reference_columns <- c(
  "atlas", "dataset_size", "cell_type", "peak",
  "full_avg_log2FC", "subset_avg_log2FC",
  "full_direction", "subset_direction",
  "subset_adjusted_p_value", "subset_is_significant",
  "effect_available", "same_direction_effect",
  "direction_switched_effect", "strict_same_direction_capture",
  "subset_analysis_status"
)

as_logical <- function(x) {
  if (is.logical(x)) {
    return(replace(x, is.na(x), FALSE))
  }
  if (is.numeric(x)) {
    return(replace(x, is.na(x), 0) != 0)
  }
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = method))
}

read_da <- function(path, expected_size) {
  x <- read_tsv(path, show_col_types = FALSE)
  missing <- setdiff(required_da_columns, names(x))
  if (length(missing) > 0L) {
    stop(path, " is missing DA columns: ", paste(missing, collapse = ", "))
  }
  x <- x %>%
    mutate(
      atlas = as.character(atlas),
      dataset_size = as.character(dataset_size),
      cell_type = as.character(cell_type),
      peak = as.character(peak),
      avg_log2FC = as.numeric(avg_log2FC),
      adjusted_p_value = as.numeric(adjusted_p_value),
      direction = as.character(direction)
    )
  if (nrow(x) > 0L) {
    if (any(x$atlas != atlas_name) || any(x$dataset_size != expected_size)) {
      stop("Atlas or dataset_size mismatch in ", path)
    }
    if (any(!x$direction %in% c("open", "closed"))) {
      stop("Invalid DA direction in ", path)
    }
    if (anyDuplicated(x[c("cell_type", "peak")])) {
      stop("Duplicate cell_type/peak rows in ", path)
    }
  }
  x
}

read_status <- function(path, expected_size) {
  x <- read_tsv(path, show_col_types = FALSE)
  missing <- setdiff(required_status_columns, names(x))
  if (length(missing) > 0L) {
    stop(path, " is missing status columns: ", paste(missing, collapse = ", "))
  }
  x <- x %>%
    mutate(
      atlas = as.character(atlas),
      dataset_size = as.character(dataset_size),
      cell_type = as.character(cell_type),
      status = as.character(status),
      n_target_cells = as.integer(n_target_cells),
      n_rest_cells = as.integer(n_rest_cells)
    )
  if (any(x$atlas != atlas_name) || any(x$dataset_size != expected_size)) {
    stop("Atlas or dataset_size mismatch in ", path)
  }
  if (anyDuplicated(x$cell_type)) {
    stop("Status table has duplicate cell types: ", path)
  }
  x
}

read_reference <- function(path, expected_size) {
  x <- read_tsv(path, show_col_types = FALSE)
  missing <- setdiff(required_reference_columns, names(x))
  if (length(missing) > 0L) {
    stop(path, " is missing reference-effect columns: ", paste(missing, collapse = ", "))
  }
  x <- x %>%
    mutate(
      atlas = as.character(atlas),
      dataset_size = as.character(dataset_size),
      cell_type = as.character(cell_type),
      peak = as.character(peak),
      full_avg_log2FC = as.numeric(full_avg_log2FC),
      subset_avg_log2FC = as.numeric(subset_avg_log2FC),
      full_direction = as.character(full_direction),
      subset_direction = as.character(subset_direction),
      subset_adjusted_p_value = as.numeric(subset_adjusted_p_value),
      subset_is_significant = as_logical(subset_is_significant),
      effect_available = as_logical(effect_available),
      same_direction_effect = as_logical(same_direction_effect),
      direction_switched_effect = as_logical(direction_switched_effect),
      strict_same_direction_capture = as_logical(strict_same_direction_capture),
      subset_analysis_status = as.character(subset_analysis_status)
    )
  if (nrow(x) > 0L) {
    if (any(x$atlas != atlas_name) || any(x$dataset_size != expected_size)) {
      stop("Atlas or dataset_size mismatch in ", path)
    }
    if (anyDuplicated(x[c("cell_type", "peak")])) {
      stop("Duplicate reference cell_type/peak rows in ", path)
    }
  }
  x
}

full_da <- read_da(full_da_path, "100pct")
full_status <- read_status(full_status_path, "100pct")
if (nrow(full_status) == 0L) {
  stop("Full-atlas status table is empty")
}

summary_results <- list()
detail_results <- list()

for (i in seq_along(subset_sizes)) {
  subset_size <- subset_sizes[[i]]
  subset_da <- read_da(subset_da_paths[[i]], subset_size)
  subset_status <- read_status(subset_status_paths[[i]], subset_size)
  reference <- read_reference(reference_effect_paths[[i]], subset_size)

  missing_status_types <- setdiff(full_status$cell_type, subset_status$cell_type)
  if (length(missing_status_types) > 0L) {
    stop(
      "Subset status table ", subset_size, " is missing cell types: ",
      paste(head(missing_status_types, 10L), collapse = ", ")
    )
  }

  details <- reference %>%
    left_join(
      full_da %>%
        select(
          cell_type,
          peak,
          full_adjusted_p_value = adjusted_p_value
        ),
      by = c("cell_type", "peak")
    ) %>%
    transmute(
      atlas,
      cell_type,
      subset_size = dataset_size,
      peak,
      full_direction,
      subset_direction,
      full_avg_log2FC,
      subset_avg_log2FC,
      full_adjusted_p_value,
      subset_adjusted_p_value,
      is_da_atlas = TRUE,
      is_da_subset = subset_is_significant,
      effect_available,
      same_direction_effect,
      same_direction_overlap = strict_same_direction_capture,
      direction_switched = direction_switched_effect,
      subset_analysis_status
    )
  detail_results[[subset_size]] <- details

  for (cell_type_name in full_status$cell_type) {
    full_status_row <- full_status %>% filter(cell_type == cell_type_name)
    subset_status_row <- subset_status %>% filter(cell_type == cell_type_name)
    full_cell_da <- full_da %>% filter(cell_type == cell_type_name)
    subset_cell_da <- subset_da %>% filter(cell_type == cell_type_name)
    reference_cell <- reference %>% filter(cell_type == cell_type_name)

    for (direction_name in c("open", "closed")) {
      atlas_peaks <- full_cell_da %>%
        filter(direction == direction_name) %>%
        pull(peak)
      subset_peaks <- subset_cell_da %>%
        filter(direction == direction_name) %>%
        pull(peak)
      direction_reference <- reference_cell %>%
        filter(full_direction == direction_name)
      available <- direction_reference %>% filter(effect_available)

      n_atlas <- length(atlas_peaks)
      n_subset <- length(subset_peaks)
      n_available <- nrow(available)
      n_same_effect <- sum(available$same_direction_effect)
      n_sign_switched <- sum(available$direction_switched_effect)
      n_strict <- sum(direction_reference$strict_same_direction_capture)
      full_ok <- identical(full_status_row$status[[1]], "success")
      subset_ok <- identical(subset_status_row$status[[1]], "success")

      if (!full_ok) {
        analysis_status <- paste0("full_atlas_", full_status_row$status[[1]])
      } else if (!subset_ok) {
        analysis_status <- paste0("subset_", subset_status_row$status[[1]])
      } else if (n_atlas == 0L) {
        analysis_status <- "no_atlas_da_peaks"
      } else if (n_available == 0L) {
        analysis_status <- "no_subset_effects_available"
      } else {
        analysis_status <- "success"
      }

      summary_results[[paste(subset_size, cell_type_name, direction_name)]] <- tibble(
        atlas = atlas_name,
        cell_type = cell_type_name,
        subset_size = subset_size,
        direction = direction_name,
        n_cells_full = full_status_row$n_target_cells[[1]],
        n_cells_subset = subset_status_row$n_target_cells[[1]],
        n_da_peaks_atlas = n_atlas,
        n_da_peaks_subset = n_subset,
        n_reference_effects_available = n_available,
        n_same_direction_effect = n_same_effect,
        n_same_direction_overlap = n_strict,
        n_sign_switched = n_sign_switched,
        percent_effect_direction_concordant = if (n_available > 0L) {
          100 * n_same_effect / n_available
        } else {
          NA_real_
        },
        percent_atlas_da_direction_concordant = if (n_atlas > 0L) {
          100 * n_same_effect / n_atlas
        } else {
          NA_real_
        },
        percent_atlas_da_captured = if (full_ok && subset_ok && n_atlas > 0L) {
          100 * n_strict / n_atlas
        } else {
          NA_real_
        },
        percent_atlas_da_sign_switched = if (n_available > 0L) {
          100 * n_sign_switched / n_available
        } else {
          NA_real_
        },
        effect_spearman_rho = safe_cor(
          available$full_avg_log2FC,
          available$subset_avg_log2FC,
          "spearman"
        ),
        effect_pearson_r = safe_cor(
          available$full_avg_log2FC,
          available$subset_avg_log2FC,
          "pearson"
        ),
        median_absolute_effect_error = if (n_available > 0L) {
          median(abs(available$subset_avg_log2FC - available$full_avg_log2FC))
        } else {
          NA_real_
        },
        precision = if (n_subset > 0L) n_strict / n_subset else NA_real_,
        jaccard = if ((n_atlas + n_subset - n_strict) > 0L) {
          n_strict / (n_atlas + n_subset - n_strict)
        } else {
          NA_real_
        },
        full_analysis_status = full_status_row$status[[1]],
        subset_analysis_status = subset_status_row$status[[1]],
        analysis_status = analysis_status
      )
    }
  }
}

summary_table <- bind_rows(summary_results) %>%
  mutate(
    subset_size = factor(subset_size, levels = subset_sizes, ordered = TRUE),
    direction = factor(direction, levels = c("open", "closed"), ordered = TRUE)
  ) %>%
  arrange(cell_type, subset_size, direction) %>%
  mutate(
    subset_size = as.character(subset_size),
    direction = as.character(direction)
  )

details_table <- bind_rows(detail_results) %>%
  mutate(subset_size = factor(subset_size, levels = subset_sizes, ordered = TRUE)) %>%
  arrange(cell_type, subset_size, peak) %>%
  mutate(subset_size = as.character(subset_size))

expected_rows <- nrow(full_status) * length(subset_sizes) * 2L
if (nrow(summary_table) != expected_rows) {
  stop("Expected ", expected_rows, " summary rows; produced ", nrow(summary_table))
}

dir.create(dirname(output_summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_details), recursive = TRUE, showWarnings = FALSE)
write_tsv(summary_table, output_summary)
write_tsv(details_table, output_details)

cat("Wrote recapitulation summary:", output_summary, "\n")
cat("Wrote atlas-reference peak details:", output_details, "\n")
cat("Sign switches are counted from available subset effects, regardless of significance.\n")
cat("Strict capture still requires subset significance and matching direction.\n")
cat("DONE\n")
