suppressPackageStartupMessages({library(BPCells);library(Signac);library(Seurat);library(readr);library(Matrix)})
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
base <- 'publication_figures/benchmarks/random_DA_stage_A_2001'
manifest <- read_tsv(file.path(base,'preflight/membership_manifest.tsv'),show_col_types=FALSE)
rows <- list()
for (a in unique(manifest$atlas)) {
  obj <- readRDS(file.path('data/processed',a,paste0(a,'_barebones_seurat.rds')))
  mat <- LayerData(obj[['peaks']],layer='counts')
  stopifnot(inherits(mat,'IterableMatrix'),!anyDuplicated(colnames(obj)),identical(colnames(mat),colnames(obj)))
  # Force a tiny disk-backed read to check backing-store access.
  stopifnot(all(is.finite(as.matrix(mat[1:2,1:2]))))
  for (f in c('02pct','10pct')) {
    r <- manifest[manifest$atlas==a & manifest$subset_size==f,]
    ids <- readLines(r$ids_file)
    stopifnot(length(ids)==r$n_selected,!anyDuplicated(ids),all(ids %in% colnames(obj)),all(ids %in% colnames(mat)))
    stopifnot(!anyNA(obj[[]][ids,'cell_type']),!anyNA(obj[[]][ids,'library_size']))
    rows[[length(rows)+1]] <- data.frame(atlas=a,subset_size=f,n_selected=length(ids),seurat_valid=TRUE,BPCells_valid=TRUE)
  }
  if (a=='kidney_atac') {
    # Small real-input smoke test: same production script, parameters and full reference definition.
    top <- names(sort(table(obj$cell_type),decreasing=TRUE))[1:2]
    cells <- unlist(lapply(top,function(t)head(colnames(obj)[obj$cell_type==t],40)),use.names=FALSE)
    ref <- read_tsv(file.path('results/signac_da',a,'100pct/signac_da_peaks.tsv'),show_col_types=FALSE)
    features <- head(unique(ref$peak[ref$cell_type %in% top]),100)
    stopifnot(length(features)>0,all(features %in% rownames(mat)))
    tiny <- Matrix(as.matrix(mat[features,cells]),sparse=TRUE)
    backing <- file.path(base,'preflight/smoke_counts')
    if (dir.exists(backing)) stop('Smoke output exists; inspect rather than overwrite')
    write_matrix_dir(mat=tiny,dir=backing)
    metadata <- obj[[]][cells,,drop=FALSE]
    small <- CreateSeuratObject(counts=open_matrix_dir(backing),assay='peaks',meta.data=metadata)
    saveRDS(small,file.path(base,'preflight/smoke_input.rds'))
    writeLines(cells,file.path(base,'preflight/smoke_ids.txt'))
    write_tsv(ref[ref$cell_type %in% top & ref$peak %in% features,],file.path(base,'preflight/smoke_reference.tsv'))
  }
  rm(obj,mat);gc()
}
write_tsv(do.call(rbind,rows),file.path(base,'preflight/seurat_validation.tsv'))
writeLines(capture.output(sessionInfo()),file.path(base,'preflight/R_sessionInfo.txt'))
cat('PASS: six memberships validated against Seurat and BPCells; small real input prepared\n')
