#!/usr/bin/env python3
"""Compute UMAP directly from the PeakVI latent representation."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Scanpy neighbors and UMAP on PeakVI latent dimensions."
    )
    parser.add_argument("--latent_tsv", required=True)
    parser.add_argument("--output_umap_tsv", required=True)
    parser.add_argument("--output_done", required=True)
    parser.add_argument("--n_neighbors", type=int, default=15)
    parser.add_argument("--min_dist", type=float, default=0.2)
    parser.add_argument("--random_seed", type=int, default=123)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    latent_tsv = Path(args.latent_tsv)
    output_umap_tsv = Path(args.output_umap_tsv)
    output_done = Path(args.output_done)

    if args.n_neighbors < 2:
        raise ValueError("--n_neighbors must be at least 2")
    if not 0 <= args.min_dist <= 1:
        raise ValueError("--min_dist must be between 0 and 1")

    print(f"Reading PeakVI latent representation: {latent_tsv}", flush=True)
    latent_df = pd.read_csv(latent_tsv, sep="\t", index_col=0)
    if latent_df.empty:
        raise ValueError("PeakVI latent TSV is empty")
    if not latent_df.index.is_unique:
        raise ValueError("Cell names in the PeakVI latent TSV are not unique")

    latent = latent_df.to_numpy(dtype=np.float32, copy=True)
    if not np.isfinite(latent).all():
        raise ValueError("PeakVI latent representation contains non-finite values")

    obs = pd.DataFrame(index=latent_df.index.astype(str))
    adata = ad.AnnData(obs=obs)
    adata.obsm["X_peakvi"] = latent

    print(
        f"Running neighbors and UMAP on {adata.n_obs:,} cells x "
        f"{latent.shape[1]} PeakVI dimensions",
        flush=True,
    )
    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        use_rep="X_peakvi",
        random_state=args.random_seed,
    )
    sc.tl.umap(
        adata,
        min_dist=args.min_dist,
        random_state=args.random_seed,
    )

    umap = np.asarray(adata.obsm["X_umap"], dtype=np.float32)
    umap_df = pd.DataFrame(
        {
            "cell_name": adata.obs_names,
            "UMAP_1": umap[:, 0],
            "UMAP_2": umap[:, 1],
        }
    )

    output_umap_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_done.parent.mkdir(parents=True, exist_ok=True)
    umap_df.to_csv(output_umap_tsv, sep="\t", index=False)
    output_done.write_text(
        "\n".join(
            [
                "done",
                f"n_cells={adata.n_obs}",
                f"n_peakvi_latent={latent.shape[1]}",
                f"n_neighbors={args.n_neighbors}",
                f"min_dist={args.min_dist}",
                f"random_seed={args.random_seed}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote UMAP coordinates: {output_umap_tsv}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
