#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args():
    p = argparse.ArgumentParser(description="Build PeakVI AnnData from sparse CSR export.")
    p.add_argument("--input_dir", required=True, help="Sparse export dir from 02_1_export_seurat_counts_for_anndata.R")
    p.add_argument("--output_h5ad", required=True, help="Output PeakVI input h5ad")
    p.add_argument("--binarize", action="store_true", help="Store X as binary accessibility instead of raw counts")
    return p.parse_args()


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_h5ad = Path(args.output_h5ad)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    shape = pd.read_csv(input_dir / "shape.tsv", sep="\t")
    n_rows = int(shape.loc[0, "n_rows"])
    n_cols = int(shape.loc[0, "n_cols"])
    nnz = int(shape.loc[0, "nnz"])

    print(f"loading sparse CSR: {n_rows:,} cells x {n_cols:,} peaks; nnz={nnz:,}", flush=True)

    data = np.memmap(input_dir / "data_float32.bin", dtype=np.float32, mode="r", shape=(nnz,))
    indices = np.memmap(input_dir / "indices_int64.bin", dtype=np.int64, mode="r", shape=(nnz,))
    indptr = np.memmap(input_dir / "indptr_int64.bin", dtype=np.int64, mode="r", shape=(n_rows + 1,))

    if int(indptr[0]) != 0:
        raise ValueError("CSR indptr must start at 0")
    if int(indptr[-1]) != nnz:
        raise ValueError("CSR indptr last value must equal nnz")

    X = sparse.csr_matrix((data, indices, indptr), shape=(n_rows, n_cols), copy=False)
    if args.binarize:
        X = X.copy()
        X.data = (X.data > 0).astype(np.float32, copy=False)
        X.eliminate_zeros()

    obs = pd.read_csv(input_dir / "obs.tsv", sep="\t", dtype={"cell_name": str})
    var = pd.read_csv(input_dir / "var.tsv", sep="\t", dtype={"peak_name": str})

    if obs.shape[0] != n_rows:
        raise ValueError(f"obs rows ({obs.shape[0]}) do not match matrix rows ({n_rows})")
    if var.shape[0] != n_cols:
        raise ValueError(f"var rows ({var.shape[0]}) do not match matrix columns ({n_cols})")
    if "cell_name" not in obs.columns:
        raise ValueError("obs.tsv must contain cell_name")
    if "peak_name" not in var.columns:
        raise ValueError("var.tsv must contain peak_name")

    obs = obs.set_index("cell_name", drop=True)
    var = var.set_index("peak_name", drop=True)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Keep a named counts layer too. This is sparse and shares the same values.
    adata.layers["counts"] = adata.X.copy()

    print(f"writing h5ad: {output_h5ad}", flush=True)
    adata.write_h5ad(output_h5ad)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
