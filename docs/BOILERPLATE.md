# Boilerplate writing guide (template)

This file is a checklist for what to write and where. Replace each placeholder with your actual content.

## docs/datasets.md
Write:
- dataset name, accession, link, citation
- organism, tissue(s), number of cells (raw + post-QC)
- file formats received (mtx/rds/h5ad/etc)
- any special preprocessing notes

## docs/methods_preprocessing.md
Write:
- QC thresholds
- feature selection (e.g., FindTopFeatures min.cutoff)
- normalization (TF-IDF)
- how peaks were filtered / restricted to standard chromosomes
- what “standardized Seurat object” guarantees

## docs/methods_nmf.md
Write:
- input matrix used for NMF (TF-IDF slot, etc.)
- NMF implementation (RcppML)
- k=100 choice rationale
- export formats (MatrixMarket + labels)

## docs/methods_subset_selection.md
Write:
- Apricot function used (FeatureBasedSelection)
- optimizer settings (two-stage)
- concave function (sqrt) and why
- fractions / subset sizes
- how you ran it for peaks vs factors

## docs/methods_evaluation.md
Write:
- differential binding method (which test, design matrix, covariates)
- what constitutes “subset performance”
- any additional preservation metrics

## docs/runtime_benchmarking.md
Write:
- what you timed (preprocess, NMF, selection)
- hardware/resources
- how you ensured fairness (same node type, same threads)
- how you summarized runtimes

## docs/reproducibility.md
Write:
- conda env YAMLs
- exact package versions
- git commit hashes
- random seeds

