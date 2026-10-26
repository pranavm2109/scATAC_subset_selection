# Installation (template)

This doc is intentionally a template. Fill in the exact commands that work on your machines.

## Platforms
- Local machine: ______________________
- HPC #1: _____________________________
- HPC #2: _____________________________

## Dependencies (high-level)
You will likely need:
- R (version: ________)
- Python (version: ________)
- Seurat + Signac (versions: ________)
- RcppML (for NMF) (version: ________)
- Python apricot-select (version: ________)
- Matrix IO tools (readMM/writeMM)
- Optional: zellkonverter / reticulate (for h5ad -> R)

## Recommended environment strategy
Use one conda env for R-heavy steps and one for Python-heavy steps.

- `env/r_env.yaml`  (R + Seurat/Signac + Bioconductor deps)
- `env/py_env.yaml` (Python + apricot-select + scientific stack)

Put the YAMLs in `env/` and document them here.

## Install on HPC (example)
1) Create envs
- `conda env create -f env/r_env.yaml`
- `conda env create -f env/py_env.yaml`

2) Verify key packages
- R: `R -q -e 'packageVersion("Seurat"); packageVersion("Signac")'`
- Python: `python -c "import apricot, numpy, scipy; print('ok')"`

3) Record versions
Write down exact versions used to run the paper under:
- `docs/reproducibility.md`

## Troubleshooting notes
Document any cluster-specific fixes here:
- module loads
- libc / gcc toolchain issues
- Bioconductor mismatch fixes
- hdf5 issues

