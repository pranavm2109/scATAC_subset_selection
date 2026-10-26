#!/usr/bin/env python3
"""Compute annotation-blind kNN sparsity weights in PeakVI latent space."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from pynndescent import NNDescent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--latent_npy", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--output_npz", required=True)
    parser.add_argument("--output_summary", required=True)
    parser.add_argument("--knn_k", type=int, default=30)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--random_seed", type=int, default=123)
    return parser.parse_args()


def read_cell_names(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any():
        raise ValueError(f"Cell-name file is empty or contains missing values: {path}")
    cell_names = values.to_numpy(dtype=str)
    if len(np.unique(cell_names)) != len(cell_names):
        raise ValueError(f"Cell names are not unique: {path}")
    return cell_names


def cell_name_checksum(cell_names: np.ndarray) -> str:
    digest = hashlib.sha256()
    for cell_name in cell_names:
        digest.update(cell_name.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def compute_knn_distances(
    latent: np.ndarray,
    k: int,
    threads: int,
    random_seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    if k < 1 or k >= len(latent):
        raise ValueError("--knn_k must be between 1 and n_cells - 1")

    index = NNDescent(
        latent,
        n_neighbors=k + 1,
        metric="euclidean",
        random_state=random_seed,
        n_jobs=threads,
        low_memory=True,
        verbose=True,
    )
    neighbors, distances = index.neighbor_graph
    mean_distance = np.empty(len(latent), dtype=np.float32)
    radius = np.empty(len(latent), dtype=np.float32)

    for start in range(0, len(latent), 100_000):
        stop = min(start + 100_000, len(latent))
        neighbor_block = neighbors[start:stop]
        distance_block = distances[start:stop]
        rows = np.arange(start, stop, dtype=neighbor_block.dtype)[:, None]
        non_self = neighbor_block != rows
        use = non_self & (np.cumsum(non_self, axis=1) <= k)
        if np.any(use.sum(axis=1) != k):
            raise RuntimeError("NN-descent returned fewer than k non-self neighbors")
        used_distances = np.where(use, distance_block, np.nan)
        mean_distance[start:stop] = np.nanmean(used_distances, axis=1)
        radius[start:stop] = np.nanmax(used_distances, axis=1)

    return mean_distance, radius


def percentile_rank(values: np.ndarray) -> np.ndarray:
    return (
        pd.Series(values).rank(method="average", pct=True).to_numpy() * 100
    ).astype(np.float32)


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be positive")

    latent = np.array(
        np.load(args.latent_npy, mmap_mode="r"),
        dtype=np.float32,
        order="C",
    )
    cell_names = read_cell_names(Path(args.cell_names))
    if latent.ndim != 2 or latent.shape[0] != len(cell_names):
        raise ValueError("PeakVI latent rows do not match the cell-name file")
    if not np.isfinite(latent).all():
        raise ValueError("PeakVI latent matrix contains non-finite values")

    print(
        f"Computing annotation-blind {args.knn_k}-NN sparsity for "
        f"{args.atlas}: {latent.shape[0]:,} cells x {latent.shape[1]} dimensions",
        flush=True,
    )
    mean_distance, radius = compute_knn_distances(
        latent,
        args.knn_k,
        args.threads,
        args.random_seed,
    )
    mean_percentile = percentile_rank(mean_distance)
    radius_percentile = percentile_rank(radius)

    output_npz = Path(args.output_npz)
    output_summary = Path(args.output_summary)
    output_npz.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    checksum = cell_name_checksum(cell_names)

    np.savez_compressed(
        output_npz,
        mean_knn_distance=mean_distance,
        knn_radius=radius,
        mean_knn_distance_percentile=mean_percentile,
        knn_radius_percentile=radius_percentile,
        cell_name_sha256=np.array(checksum),
        knn_k=np.array(args.knn_k, dtype=np.int32),
    )
    pd.DataFrame(
        {
            "atlas": [args.atlas],
            "n_cells": [latent.shape[0]],
            "n_peakvi_dimensions": [latent.shape[1]],
            "knn_k": [args.knn_k],
            "distance_metric": ["euclidean"],
            "neighbor_method": ["pynndescent.NNDescent"],
            "random_seed": [args.random_seed],
            "cell_name_sha256": [checksum],
            "mean_knn_distance_min": [float(mean_distance.min())],
            "mean_knn_distance_median": [float(np.median(mean_distance))],
            "mean_knn_distance_max": [float(mean_distance.max())],
        }
    ).to_csv(output_summary, sep="\t", index=False)
    print(f"Wrote density metrics: {output_npz}", flush=True)
    print(f"Wrote density summary: {output_summary}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
