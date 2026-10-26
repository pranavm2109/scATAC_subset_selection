#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

cat("libs loaded\n")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/02_4_extract_cell_names.R",
      "<input_barebones_seurat_rds>",
      "<output_cell_names_txt>"
    )
  )
}

input_rds <- args[1]
output_cells <- args[2]

cat("input_rds: ", input_rds, "\n", sep = "")
cat("output_cells: ", output_cells, "\n", sep = "")

dir.create(dirname(output_cells), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "peaks"

cell_names <- colnames(obj)

cat("extracted cell names\n")
print(length(cell_names))
print(head(cell_names))

write.table(
  cell_names,
  file = output_cells,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("saved cell names\n")
cat("DONE\n")