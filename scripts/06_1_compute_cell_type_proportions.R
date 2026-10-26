#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

cat("libs loaded\n")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/06_1_compute_cell_type_proportions.R",
      "<input_seurat_rds>",
      "<subset_dir>",
      "<output_wide_tsv>",
      "[output_tidy_tsv]",
      "<fractions_csv>",
      "[subset_file_prefix=selected_cells_peakvi]"
    )
  )
}

input_rds <- args[1]
subset_dir <- args[2]
output_tsv <- args[3]

# Backward compatible argument handling:
# old: Rscript 06_1 ... <output_wide_tsv> <fractions_csv> [prefix]
# new: Rscript 06_1 ... <output_wide_tsv> <output_tidy_tsv> <fractions_csv> [prefix]
looks_like_tsv_path <- grepl("\\.tsv$", args[4]) || grepl("/", args[4], fixed = TRUE)
if (length(args) >= 5 && looks_like_tsv_path) {
  output_tidy_tsv <- args[4]
  fractions_arg <- args[5]
  subset_file_prefix <- if (length(args) >= 6) args[6] else "selected_cells_peakvi"
} else {
  output_tidy_tsv <- NA_character_
  fractions_arg <- args[4]
  subset_file_prefix <- if (length(args) >= 5) args[5] else "selected_cells_peakvi"
}

fractions <- strsplit(fractions_arg, ",", fixed = TRUE)[[1]]
fractions <- trimws(fractions)

if (length(fractions) == 0) {
  stop("No fractions provided")
}

subset_size_label <- function(frac) {
  sub("pct$", "%", frac)
}

check_fraction_sum <- function(df, subset_label, tolerance = 1e-6) {
  fraction_sum <- sum(df$fraction, na.rm = FALSE)
  if (!is.finite(fraction_sum) || abs(fraction_sum - 1) > tolerance) {
    stop(
      paste0(
        "Fractions for subset_size=", subset_label,
        " sum to ", signif(fraction_sum, 12),
        ", expected approximately 1."
      )
    )
  }
}

cat("input_rds: ", input_rds, "\n", sep = "")
cat("subset_dir: ", subset_dir, "\n", sep = "")
cat("output_wide_tsv: ", output_tsv, "\n", sep = "")
cat("output_tidy_tsv: ", ifelse(is.na(output_tidy_tsv), "not requested", output_tidy_tsv), "\n", sep = "")
cat("fractions: ", paste(fractions, collapse = ","), "\n", sep = "")
cat("subset_file_prefix: ", subset_file_prefix, "\n", sep = "")

dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
if (!is.na(output_tidy_tsv)) {
  dir.create(dirname(output_tidy_tsv), recursive = TRUE, showWarnings = FALSE)
}

obj <- readRDS(input_rds)

if (!"cell_type" %in% colnames(obj@meta.data)) {
  stop("cell_type column not found in Seurat metadata.")
}

atlas_cells <- colnames(obj)
atlas_n <- length(atlas_cells)

atlas_cell_type <- as.character(obj@meta.data$cell_type)
names(atlas_cell_type) <- atlas_cells

if (atlas_n == 0) {
  stop("No cells in atlas object.")
}

atlas_counts <- table(atlas_cell_type)
cell_types <- sort(names(atlas_counts))

result <- data.frame(
  cell_type = cell_types,
  atlas_count = as.numeric(atlas_counts[cell_types]),
  stringsAsFactors = FALSE
)
result$atlas_proportion <- result$atlas_count / atlas_n

tidy_rows <- list()

for (frac in fractions) {
  subset_file <- file.path(subset_dir, paste0(subset_file_prefix, "_", frac, ".txt"))
  if (!file.exists(subset_file)) {
    stop(paste0("Missing subset file: ", subset_file))
  }

  subset_cells <- readLines(subset_file, warn = FALSE)
  subset_cells <- unique(subset_cells[nzchar(subset_cells)])
  subset_cells <- intersect(subset_cells, atlas_cells)

  subset_n <- length(subset_cells)
  if (subset_n == 0) {
    stop(paste0("No matching cells found for subset file: ", subset_file))
  }

  subset_types <- atlas_cell_type[subset_cells]
  subset_counts <- table(subset_types)

  subset_count_col <- paste0("subset_", frac, "_count")
  subset_prop_col <- paste0("subset_", frac, "_proportion")
  ratio_col <- paste0("ratio_", frac, "_vs_atlas")

  result[[subset_count_col]] <- as.numeric(subset_counts[result$cell_type])
  result[[subset_count_col]][is.na(result[[subset_count_col]])] <- 0

  result[[subset_prop_col]] <- result[[subset_count_col]] / subset_n

  result[[ratio_col]] <- ifelse(
    result$atlas_proportion > 0,
    result[[subset_prop_col]] / result$atlas_proportion,
    NA_real_
  )

  tidy_subset <- data.frame(
    atlas = basename(dirname(input_rds)),
    subset_size = subset_size_label(frac),
    cell_type = result$cell_type,
    n_cells = result[[subset_count_col]],
    fraction = result[[subset_prop_col]],
    stringsAsFactors = FALSE
  )
  check_fraction_sum(tidy_subset, subset_size_label(frac))
  tidy_rows[[length(tidy_rows) + 1]] <- tidy_subset
}

tidy_atlas <- data.frame(
  atlas = basename(dirname(input_rds)),
  subset_size = "100%",
  cell_type = result$cell_type,
  n_cells = result$atlas_count,
  fraction = result$atlas_proportion,
  stringsAsFactors = FALSE
)
check_fraction_sum(tidy_atlas, "100%")
tidy_rows[[length(tidy_rows) + 1]] <- tidy_atlas
tidy_result <- do.call(rbind, tidy_rows)

write.table(
  result,
  file = output_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("saved proportions table\n")

if (!is.na(output_tidy_tsv)) {
  write.table(
    tidy_result,
    file = output_tidy_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  cat("saved tidy cell-type fractions table\n")
}

cat("DONE\n")
