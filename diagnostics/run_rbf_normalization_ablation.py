#!/usr/bin/env python3
"""Run one annotation-blind RBF row-normalization selection at 2%.

This diagnostic changes only cell-wise scaling of the existing unweighted RBF
landmark matrix. It does not load cell-type annotations. The production subset
selection and its outputs are not modified.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from apricot import FeatureBasedSelection


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--method", choices=["l1", "singleton_capped"], required=True)
    parser.add_argument("--feature_matrix", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--landmark_metadata", required=True)
    parser.add_argument("--baseline_selected_cells", required=True)
    parser.add_argument("--output_selected_cells", required=True)
    parser.add_argument("--output_summary_tsv", required=True)
    parser.add_argument("--fraction", type=float, default=0.02)
    parser.add_argument("--min_gain_ratio", type=float, default=0.5)
    parser.add_argument("--max_gain_ratio", type=float, default=2.0)
    parser.add_argument("--chunk_size", type=int, default=25000)
    parser.add_argument("--random_seed", type=int, default=123)
    return parser.parse_args()


def read_ids(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any() or values.duplicated().any():
        raise ValueError(f"Cell IDs must be nonempty, unique, and complete: {path}")
    return values.to_numpy(dtype=str)


def quantiles(values: np.ndarray, prefix: str) -> dict[str, float]:
    probs = [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1]
    names = ["min", "q01", "q10", "median", "q90", "q99", "max"]
    return {
        f"{prefix}_{name}": float(value)
        for name, value in zip(names, np.quantile(values, probs))
    }


def main() -> None:
    args = parse_args()
    if not 0 < args.fraction <= 1:
        raise ValueError("--fraction must be > 0 and <= 1")
    if not 0 < args.min_gain_ratio <= args.max_gain_ratio:
        raise ValueError("Require 0 < min_gain_ratio <= max_gain_ratio")

    cell_ids = read_ids(Path(args.cell_names))
    baseline_ids = read_ids(Path(args.baseline_selected_cells))
    features = np.load(args.feature_matrix, mmap_mode="r")
    metadata = json.loads(Path(args.landmark_metadata).read_text(encoding="utf-8"))

    if features.ndim != 2 or features.shape[0] != len(cell_ids):
        raise ValueError("RBF feature rows do not match cell names")
    if str(metadata.get("atlas")) != args.atlas:
        raise ValueError("Atlas does not match landmark metadata")
    if tuple(metadata.get("feature_matrix_shape", [])) != features.shape:
        raise ValueError("Feature matrix shape does not match landmark metadata")
    if metadata.get("peakvi_latent_distribution") != "normal":
        raise ValueError("This ablation expects Gaussian PeakVI latent dimensions")
    if int(metadata.get("n_landmarks", -1)) != features.shape[1]:
        raise ValueError("Landmark count does not match feature matrix")

    subset_size = max(int(round(len(cell_ids) * args.fraction)), 1)
    if len(baseline_ids) != subset_size:
        raise ValueError(
            f"Baseline contains {len(baseline_ids):,} cells; expected {subset_size:,}"
        )
    if np.any(pd.Index(cell_ids).get_indexer(baseline_ids) < 0):
        raise ValueError("Baseline selection contains unknown cell IDs")

    row_sum = np.empty(len(cell_ids), dtype=np.float64)
    singleton_gain = np.empty(len(cell_ids), dtype=np.float64)
    for start in range(0, len(cell_ids), args.chunk_size):
        stop = min(start + args.chunk_size, len(cell_ids))
        block = np.asarray(features[start:stop], dtype=np.float32)
        if not np.isfinite(block).all() or np.any(block < 0):
            raise ValueError("RBF matrix contains negative or non-finite values")
        row_sum[start:stop] = block.sum(axis=1, dtype=np.float64)
        singleton_gain[start:stop] = np.sqrt(block).sum(axis=1, dtype=np.float64)

    if np.any(row_sum <= 0) or np.any(singleton_gain <= 0):
        raise ValueError("RBF normalization requires every row to have positive mass")

    if args.method == "l1":
        target = float(np.median(row_sum))
        scale = target / row_sum
        normalization_target = "median_rbf_row_sum"
    else:
        target = float(np.median(singleton_gain))
        raw_gain_ratio = target / singleton_gain
        gain_ratio = np.clip(
            raw_gain_ratio,
            args.min_gain_ratio,
            args.max_gain_ratio,
        )
        scale = np.square(gain_ratio)
        normalization_target = "median_singleton_sqrt_gain"

    normalized = np.empty(features.shape, dtype=np.float32)
    for start in range(0, len(cell_ids), args.chunk_size):
        stop = min(start + args.chunk_size, len(cell_ids))
        np.multiply(
            features[start:stop],
            scale[start:stop, None],
            out=normalized[start:stop],
            casting="unsafe",
        )

    print(
        f"{args.atlas}: selecting {subset_size:,}/{len(cell_ids):,} cells; "
        f"method={args.method}; features={features.shape[1]:,}",
        flush=True,
    )
    print(
        f"Scale min/median/max={scale.min():.6g}/"
        f"{np.median(scale):.6g}/{scale.max():.6g}",
        flush=True,
    )

    selector = FeatureBasedSelection(
        subset_size,
        "sqrt",
        optimizer="two-stage",
        random_state=args.random_seed,
        verbose=True,
    )
    started = time.time()
    selector.fit(normalized)
    elapsed = time.time() - started
    ranking = np.asarray(selector.ranking, dtype=np.int64)
    gains = np.asarray(selector.gains, dtype=np.float64)
    del normalized

    if len(ranking) != subset_size or len(np.unique(ranking)) != subset_size:
        raise RuntimeError("Apricot returned an invalid selection")

    output_selected = Path(args.output_selected_cells)
    output_summary = Path(args.output_summary_tsv)
    output_selected.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    pd.Series(cell_ids[ranking]).to_csv(output_selected, index=False, header=False)

    summary = {
        "atlas": args.atlas,
        "method": args.method,
        "fraction": args.fraction,
        "n_cells": len(cell_ids),
        "n_selected": subset_size,
        "n_features": features.shape[1],
        "normalization_target": normalization_target,
        "target_value": target,
        "min_gain_ratio": args.min_gain_ratio if args.method == "singleton_capped" else np.nan,
        "max_gain_ratio": args.max_gain_ratio if args.method == "singleton_capped" else np.nan,
        "elapsed_seconds": elapsed,
        "first_gain": float(gains[0]),
        "last_gain": float(gains[-1]),
        "sum_selection_gains": float(gains.sum()),
        **quantiles(row_sum, "original_row_sum"),
        **quantiles(singleton_gain, "original_singleton_gain"),
        **quantiles(scale, "feature_scale"),
    }
    pd.DataFrame([summary]).to_csv(output_summary, sep="\t", index=False)
    print(f"Wrote selected cells: {output_selected}", flush=True)
    print(f"Wrote run summary: {output_summary}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
