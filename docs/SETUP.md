# Setup guide (template)

Goal: standardize 4 heterogeneous scATAC atlases into a single downstream-compatible format, then run:
1) subset selection on raw peaks
2) subset selection on NMF k=100 factors
…compare runtime and evaluate subsets by differential binding.

## 1) Add atlases
For each atlas, create:
- `data/raw/<atlas>/` and put raw downloads there.

Fill out `config/atlases.tsv` with:
- atlas_id
- format (mtx/rds/h5ad/arrow/etc)
- paths to counts/fragments/peaks/metadata
- genome build (hg38/mm10/etc)
- notes on cell_type column name

## 2) Preprocessing stage (standardization)
Write one preprocessing script per atlas under `scripts/preprocessing/` OR one script that dispatches by atlas_id.

Output contract (important):
- Save: `data/processed/<atlas>/seurat_standardized.rds`
- The Seurat object must include:
  - ChromatinAssay in `assays$peaks` (or consistent name you choose)
  - `meta.data$cell_type` populated
  - cell barcodes are consistent everywhere
  - genome build recorded (in object misc, or a sidecar file)

## 3) NMF stage
Run NMF (k=100) and export:
- `results/nmf/<atlas>/nmf_k100_cell_by_factor.mtx`
- `results/nmf/<atlas>/nmf_k100_cells.txt`
- `results/nmf/<atlas>/nmf_k100_factors.txt`
Also store runtime:
- `results/nmf/<atlas>/nmf_runtime.tsv`

## 4) Subset selection stage
For each atlas and each fraction (2/4/6/8/10%):
- raw peaks pipeline outputs:
  - `results/subsets/raw_peaks/<atlas>/subset_<frac>.txt`
- nmf pipeline outputs:
  - `results/subsets/nmf_k100/<atlas>/subset_<frac>.txt`
Also store runtime:
- `results/subsets/.../<atlas>/runtime.tsv`

## 5) Evaluation stage (differential binding)
For each atlas + subset fraction + pipeline:
- run differential binding on subset vs full (or subset vs complement, depending on your design)
- output:
  - `results/evaluation/differential_binding/<atlas>/<pipeline>/db_results.rds`
  - `results/evaluation/differential_binding/<atlas>/<pipeline>/summary_metrics.tsv`

## 6) Benchmarking stage (runtime + plots)
Aggregate runtimes across:
- preprocessing
- NMF
- subset selection (raw vs nmf)
Produce:
- `results/benchmarking/runtime_summary.tsv`
- plots under `results/benchmarking/runtime_plots/`

## 7) Figures
Final figures only go under `figures/`.
Exploratory plots should stay in `results/` or `notebooks/`.

