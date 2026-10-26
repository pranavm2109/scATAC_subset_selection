#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
subset_sizes <- c("02pct", "04pct", "06pct", "08pct", "10pct")
dataset_sizes <- c("100pct", subset_sizes)
dataset_labels <- c(
  "100pct" = "Full atlas",
  "02pct" = "2% subset",
  "04pct" = "4% subset",
  "06pct" = "6% subset",
  "08pct" = "8% subset",
  "10pct" = "10% subset"
)

if (length(args) != 19L) {
  stop(
    paste(
      "Usage: Rscript scripts/09_3_plot_signac_da_capture.R",
      "<atlas> <full_da.tsv> <full_status.tsv> <capture_summary.tsv>",
      "<histogram_dir> <capture_heatmap.png> <trajectories.png>",
      "<bubble_heatmap.png> <histogram_manifest.tsv>",
      "<five_subset_da_tsvs> <five_subset_status_tsvs>"
    )
  )
}

atlas_name <- args[[1]]
full_da_path <- args[[2]]
full_status_path <- args[[3]]
capture_summary_path <- args[[4]]
histogram_dir <- args[[5]]
output_heatmap <- args[[6]]
output_trajectories <- args[[7]]
output_bubble <- args[[8]]
output_manifest <- args[[9]]
subset_da_paths <- args[10:14]
subset_status_paths <- args[15:19]

required_da_columns <- c(
  "atlas", "dataset_size", "cell_type", "peak", "avg_log2FC",
  "direction"
)
required_status_columns <- c(
  "atlas", "dataset_size", "cell_type", "n_target_cells", "status"
)
required_capture_columns <- c(
  "atlas", "cell_type", "subset_size", "direction",
  "n_cells_full", "n_cells_subset", "n_da_peaks_atlas",
  "n_da_peaks_subset", "n_same_direction_overlap",
  "n_sign_switched", "percent_atlas_da_sign_switched",
  "percent_effect_direction_concordant",
  "percent_atlas_da_captured", "precision", "jaccard",
  "full_analysis_status", "subset_analysis_status", "analysis_status"
)

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
      direction = as.character(direction)
    )
  if (nrow(x) > 0L) {
    if (any(x$atlas != atlas_name) || any(x$dataset_size != expected_size)) {
      stop("Atlas or dataset_size mismatch in ", path)
    }
    if (any(!x$direction %in% c("open", "closed"))) {
      stop("Invalid direction in ", path)
    }
    if (any(
      (x$direction == "open" & x$avg_log2FC <= 0) |
        (x$direction == "closed" & x$avg_log2FC >= 0)
    )) {
      stop("DA direction and avg_log2FC sign disagree in ", path)
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
      n_target_cells = as.integer(n_target_cells),
      status = as.character(status)
    )
  if (any(x$atlas != atlas_name) || any(x$dataset_size != expected_size)) {
    stop("Atlas or dataset_size mismatch in ", path)
  }
  if (anyDuplicated(x$cell_type)) {
    stop("Status table has duplicate cell types: ", path)
  }
  x
}

full_da <- read_da(full_da_path, "100pct")
full_status <- read_status(full_status_path, "100pct")
da_tables <- list(full_da)
status_tables <- list(full_status)

for (i in seq_along(subset_sizes)) {
  da_tables[[i + 1L]] <- read_da(subset_da_paths[[i]], subset_sizes[[i]])
  status_tables[[i + 1L]] <- read_status(
    subset_status_paths[[i]],
    subset_sizes[[i]]
  )
}

all_da <- bind_rows(da_tables) %>%
  mutate(
    dataset_size = factor(dataset_size, levels = dataset_sizes, ordered = TRUE),
    direction = factor(direction, levels = c("open", "closed"), ordered = TRUE)
  )
all_status <- bind_rows(status_tables) %>%
  mutate(
    dataset_size = factor(dataset_size, levels = dataset_sizes, ordered = TRUE)
  )

cell_types <- full_status$cell_type
if (length(cell_types) == 0L) {
  stop("Full-atlas status table contains no cell types")
}
for (size_name in dataset_sizes) {
  observed <- all_status %>%
    filter(as.character(dataset_size) == size_name) %>%
    pull(cell_type)
  missing_types <- setdiff(cell_types, observed)
  if (length(missing_types) > 0L) {
    stop(
      "Status table for ", size_name, " is missing cell types: ",
      paste(head(missing_types, 10L), collapse = ", ")
    )
  }
}

capture <- read_tsv(capture_summary_path, show_col_types = FALSE)
missing_capture <- setdiff(required_capture_columns, names(capture))
if (length(missing_capture) > 0L) {
  stop(
    capture_summary_path, " is missing capture columns: ",
    paste(missing_capture, collapse = ", ")
  )
}
capture <- capture %>%
  mutate(
    atlas = as.character(atlas),
    cell_type = as.character(cell_type),
    subset_size = factor(
      as.character(subset_size),
      levels = subset_sizes,
      ordered = TRUE
    ),
    direction = factor(
      as.character(direction),
      levels = c("open", "closed"),
      ordered = TRUE
    ),
    percent_atlas_da_captured = as.numeric(percent_atlas_da_captured),
    percent_effect_direction_concordant = as.numeric(percent_effect_direction_concordant),
    n_da_peaks_atlas = as.numeric(n_da_peaks_atlas),
    estimable = (
      full_analysis_status == "success" &
        subset_analysis_status == "success" &
        is.finite(percent_effect_direction_concordant)
    )
  )
if (any(capture$atlas != atlas_name)) {
  stop("Atlas mismatch in capture summary")
}

safe_filename <- function(value) {
  result <- gsub("[^A-Za-z0-9._-]+", "_", value)
  result <- gsub("^_+|_+$", "", result)
  if (nchar(result) == 0L) {
    result <- "cell_type"
  }
  substr(result, 1L, 120L)
}

dir.create(histogram_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_heatmap), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_trajectories), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_bubble), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_manifest), recursive = TRUE, showWarnings = FALSE)

# Remove stale PNGs so the directory exactly reflects the current cell types.
old_histograms <- list.files(
  histogram_dir,
  pattern = "\\.png$",
  full.names = TRUE
)
if (length(old_histograms) > 0L) {
  unlink(old_histograms)
}

direction_colors <- c("open" = "#D95F02", "closed" = "#2C7FB8")
manifest_rows <- vector("list", length(cell_types))

for (cell_type_index in seq_along(cell_types)) {
  cell_type_name <- cell_types[[cell_type_index]]
  cell_da <- all_da %>% filter(cell_type == cell_type_name)
  cell_status <- all_status %>% filter(cell_type == cell_type_name)

  max_abs_effect <- suppressWarnings(max(abs(cell_da$avg_log2FC), na.rm = TRUE))
  if (!is.finite(max_abs_effect) || max_abs_effect <= 0) {
    max_abs_effect <- 1
  }
  max_abs_effect <- max_abs_effect * 1.05
  bin_width <- (2 * max_abs_effect) / 40

  annotations <- tibble(
    dataset_size = factor(dataset_sizes, levels = dataset_sizes, ordered = TRUE)
  ) %>%
    left_join(
      cell_da %>%
        group_by(dataset_size) %>%
        summarise(
          n_open = sum(direction == "open"),
          n_closed = sum(direction == "closed"),
          .groups = "drop"
        ),
      by = "dataset_size"
    ) %>%
    mutate(
      n_open = replace_na(n_open, 0L),
      n_closed = replace_na(n_closed, 0L)
    ) %>%
    left_join(
      cell_status %>% select(dataset_size, status, n_target_cells),
      by = "dataset_size"
    ) %>%
    mutate(
      label = paste0(
        "Open: ", n_open, "\nClosed: ", n_closed,
        "\nCells: ", comma(n_target_cells),
        if_else(status == "success", "", paste0("\n", status))
      )
    )

  histogram_plot <- ggplot(
    cell_da,
    aes(x = avg_log2FC, fill = direction)
  ) +
    geom_histogram(
      binwidth = bin_width,
      boundary = 0,
      color = "white",
      linewidth = 0.15,
      na.rm = TRUE
    ) +
    geom_vline(xintercept = 0, color = "#303030", linewidth = 0.55) +
    geom_text(
      data = annotations,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.1,
      size = 3.2,
      lineheight = 1.05
    ) +
    facet_wrap(
      ~dataset_size,
      ncol = 2,
      labeller = as_labeller(dataset_labels),
      drop = FALSE
    ) +
    scale_fill_manual(values = direction_colors, drop = FALSE) +
    scale_x_continuous(
      limits = c(-max_abs_effect, max_abs_effect),
      breaks = pretty_breaks(n = 5)
    ) +
    labs(
      title = paste0(atlas_name, ": ", cell_type_name),
      subtitle = paste(
        "Signed effect sizes for significant target-vs-rest Signac DA peaks;",
        "positive is open and negative is closed"
      ),
      x = "Signac DA effect size (average log2 fold change)",
      y = "Number of significant DA peaks",
      fill = "Direction"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      axis.text = element_text(color = "black"),
      panel.spacing = unit(1, "lines")
    )

  histogram_file <- file.path(
    histogram_dir,
    sprintf(
      "%03d_%s_signed_effect_histogram.png",
      cell_type_index,
      safe_filename(cell_type_name)
    )
  )
  ggsave(
    histogram_file,
    histogram_plot,
    width = 15,
    height = 9,
    units = "in",
    dpi = 320,
    bg = "white"
  )

  manifest_rows[[cell_type_index]] <- tibble(
    atlas = atlas_name,
    cell_type = cell_type_name,
    histogram_file = histogram_file
  )
}

manifest <- bind_rows(manifest_rows)
write_tsv(manifest, output_manifest)

cell_type_order <- full_status %>%
  arrange(desc(n_target_cells), cell_type) %>%
  pull(cell_type)
capture <- capture %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(cell_type_order)),
    display_value = if_else(estimable, percent_effect_direction_concordant, NA_real_)
  )

figure_height <- max(8, 0.24 * length(cell_types) + 3)

capture_heatmap <- ggplot(
  capture,
  aes(x = subset_size, y = cell_type, fill = display_value)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(
    ~direction,
    ncol = 2,
    labeller = as_labeller(c("open" = "Open peaks", "closed" = "Closed peaks"))
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 100),
    oob = squish,
    na.value = "#D9D9D9",
    name = "Atlas DA effects\ndirection-concordant (%)"
  ) +
  scale_x_discrete(labels = c(
    "02pct" = "2%", "04pct" = "4%", "06pct" = "6%", "08pct" = "8%", "10pct" = "10%"
  )) +
  labs(
    title = paste0(atlas_name, ": direction-specific DA-peak capture"),
    subtitle = paste(
      "Percentage of full-atlas DA effects retaining",
      "the same direction in each subset; gray means not estimable"
    ),
    x = "Selected subset size",
    y = "Cell type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black", size = 8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  output_heatmap,
  capture_heatmap,
  width = 14,
  height = figure_height,
  units = "in",
  dpi = 320,
  bg = "white",
  limitsize = FALSE
)

median_trajectories <- capture %>%
  filter(estimable) %>%
  group_by(direction, subset_size) %>%
  summarise(
    median_capture = median(percent_effect_direction_concordant, na.rm = TRUE),
    .groups = "drop"
  )

trajectory_plot <- ggplot(
  capture %>% filter(estimable),
  aes(
    x = subset_size,
    y = percent_effect_direction_concordant,
    group = cell_type
  )
) +
  geom_line(color = "#7A8793", alpha = 0.24, linewidth = 0.45) +
  geom_point(color = "#7A8793", alpha = 0.32, size = 1.1) +
  geom_line(
    data = median_trajectories,
    aes(x = subset_size, y = median_capture, group = 1),
    inherit.aes = FALSE,
    color = "#B2182B",
    linewidth = 1.4
  ) +
  geom_point(
    data = median_trajectories,
    aes(x = subset_size, y = median_capture),
    inherit.aes = FALSE,
    color = "#B2182B",
    size = 2.8
  ) +
  facet_wrap(
    ~direction,
    ncol = 2,
    labeller = as_labeller(c("open" = "Open peaks", "closed" = "Closed peaks"))
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_discrete(labels = c(
    "02pct" = "2%", "04pct" = "4%", "06pct" = "6%", "08pct" = "8%", "10pct" = "10%"
  )) +
  labs(
    title = paste0(atlas_name, ": DA-peak capture trajectories"),
    subtitle = paste(
      "Thin lines represent estimable cell types;",
      "the red line is the median across cell types"
    ),
    x = "Selected subset size",
    y = "Atlas DA effects with matching direction"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black")
  )

ggsave(
  output_trajectories,
  trajectory_plot,
  width = 13,
  height = 6.5,
  units = "in",
  dpi = 320,
  bg = "white"
)

estimable_bubbles <- capture %>% filter(estimable)
missing_bubbles <- capture %>% filter(!estimable)

bubble_plot <- ggplot() +
  geom_point(
    data = estimable_bubbles,
    aes(
      x = subset_size,
      y = cell_type,
      size = n_da_peaks_atlas,
      fill = percent_effect_direction_concordant
    ),
    shape = 21,
    color = "#333333",
    stroke = 0.25,
    alpha = 0.92
  ) +
  geom_point(
    data = missing_bubbles,
    aes(x = subset_size, y = cell_type),
    shape = 4,
    color = "#8C8C8C",
    size = 2.4,
    stroke = 0.7
  ) +
  facet_wrap(
    ~direction,
    ncol = 2,
    labeller = as_labeller(c("open" = "Open peaks", "closed" = "Closed peaks"))
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 100),
    oob = squish,
    name = "Atlas DA effects\ndirection-concordant (%)"
  ) +
  scale_size_area(
    max_size = 9,
    breaks = pretty_breaks(n = 4),
    name = "Full-atlas\nDA peaks"
  ) +
  scale_x_discrete(labels = c(
    "02pct" = "2%", "04pct" = "4%", "06pct" = "6%", "08pct" = "8%", "10pct" = "10%"
  )) +
  labs(
    title = paste0(atlas_name, ": DA-peak capture overview"),
    subtitle = paste(
      "Color shows effect-direction concordance; bubble size shows the number",
      "of full-atlas DA peaks; crosses indicate non-estimable comparisons"
    ),
    x = "Selected subset size",
    y = "Cell type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black", size = 8),
    panel.grid.major = element_line(color = "#ECECEC", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  output_bubble,
  bubble_plot,
  width = 15,
  height = figure_height,
  units = "in",
  dpi = 320,
  bg = "white",
  limitsize = FALSE
)

cat("Wrote", nrow(manifest), "cell-type effect histograms to:", histogram_dir, "\n")
cat("Wrote capture heatmap:", output_heatmap, "\n")
cat("Wrote capture trajectories:", output_trajectories, "\n")
cat("Wrote bubble heatmap:", output_bubble, "\n")
cat("Wrote histogram manifest:", output_manifest, "\n")
cat("DONE\n")
