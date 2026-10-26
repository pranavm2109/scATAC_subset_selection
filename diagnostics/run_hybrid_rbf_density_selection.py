#!/usr/bin/env python3
"""Run one annotation-blind hybrid RBF+density-RBF Apricot selection.

The first feature block is the original RBF landmark matrix X. The second is
the exact density-weighted copy wX, where the quadratic kNN-sparsity weight is
capped at 2.0. The blocks are concatenated without block normalization or an
additional tuning coefficient. Cell-type metadata is not accepted.
"""

from __future__ import annotations

import argparse
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd
from apricot import FeatureBasedSelection


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature_matrix", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--sparsity_npz", required=True)
    parser.add_argument("--baseline_selected_cells", required=True)
    parser.add_argument("--fraction", type=float, required=True)
    parser.add_argument("--fraction_label", required=True)
    parser.add_argument("--output_selected_cells", required=True)
    parser.add_argument("--output_selected_metrics", required=True)
    parser.add_argument("--output_done", required=True)
    parser.add_argument("--max_weight", type=float, default=2.0)
    parser.add_argument("--random_seed", type=int, default=123)
    return parser.parse_args()


def read_ids(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any():
        raise ValueError(f"Empty or missing cell IDs: {path}")
    result = values.to_numpy(dtype=str)
    if len(np.unique(result)) != len(result):
        raise ValueError(f"Duplicate cell IDs: {path}")
    return result


def cell_name_checksum(cell_ids: np.ndarray) -> str:
    digest = hashlib.sha256()
    for cell_id in cell_ids:
        digest.update(str(cell_id).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    if not 0 < args.fraction <= 1:
        raise ValueError("--fraction must be > 0 and <= 1")
    if not np.isclose(args.max_weight, 2.0):
        raise ValueError("This frozen ablation requires --max_weight 2.0")

    cell_ids = read_ids(Path(args.cell_names))
    features = np.load(args.feature_matrix, mmap_mode="r")
    if features.ndim != 2 or features.shape[0] != len(cell_ids):
        raise ValueError("RBF feature rows do not match the cell-name file")
    if float(features.min()) < 0 or not np.isfinite(features).all():
        raise ValueError("RBF feature matrix contains negative or non-finite values")

    with np.load(args.sparsity_npz, allow_pickle=False) as sparsity:
        required = {
            "mean_knn_distance",
            "mean_knn_distance_percentile",
            "cell_name_sha256",
            "knn_k",
        }
        missing = required - set(sparsity.files)
        if missing:
            raise ValueError(f"Sparsity NPZ is missing arrays: {sorted(missing)}")
        mean_distance = np.asarray(sparsity["mean_knn_distance"], dtype=np.float32)
        percentile = np.asarray(
            sparsity["mean_knn_distance_percentile"], dtype=np.float32
        )
        expected_checksum = str(sparsity["cell_name_sha256"].item())
        knn_k = int(sparsity["knn_k"].item())

    if len(percentile) != len(cell_ids) or len(mean_distance) != len(cell_ids):
        raise ValueError("Sparsity vectors and cell IDs have different lengths")
    if expected_checksum != cell_name_checksum(cell_ids):
        raise ValueError("Sparsity metrics and RBF features use different cell orders")
    if not np.isfinite(percentile).all() or np.any((percentile < 0) | (percentile > 100)):
        raise ValueError("Invalid kNN sparsity percentiles")

    baseline_ids = read_ids(Path(args.baseline_selected_cells))
    subset_size = max(int(round(len(cell_ids) * args.fraction)), 1)
    if len(baseline_ids) != subset_size:
        raise ValueError(
            f"Baseline has {len(baseline_ids):,} cells; expected {subset_size:,}"
        )
    if np.any(pd.Index(cell_ids).get_indexer(baseline_ids) < 0):
        raise ValueError("Baseline selection contains unknown cell IDs")

    scaled = percentile / np.float32(100.0)
    weights = np.float32(1.0) + np.square(scaled)
    weights = np.minimum(weights, np.float32(2.0))
    print(
        f"Selecting {subset_size:,}/{len(cell_ids):,} cells from "
        f"{2 * features.shape[1]:,} hybrid features "
        f"({features.shape[1]:,} X + {features.shape[1]:,} wX); "
        f"k={knn_k}; max_weight=2.0",
        flush=True,
    )
    print(
        f"Weight min/median/max: {weights.min():.4f}/"
        f"{np.median(weights):.4f}/{weights.max():.4f}",
        flush=True,
    )

    # The weighting is cell-specific and annotation-blind. No metadata has been
    # loaded, and this script has no metadata argument.
    hybrid_features = np.empty(
        (features.shape[0], 2 * features.shape[1]), dtype=np.float32
    )
    hybrid_features[:, : features.shape[1]] = features
    np.multiply(
        features,
        weights[:, None],
        out=hybrid_features[:, features.shape[1] :],
    )
    selector = FeatureBasedSelection(
        subset_size,
        "sqrt",
        optimizer="two-stage",
        random_state=args.random_seed,
        verbose=True,
    )
    start = time.time()
    selector.fit(hybrid_features)
    elapsed = time.time() - start
    ranking = np.asarray(selector.ranking, dtype=np.int64)
    gains = np.asarray(selector.gains, dtype=np.float64)
    del hybrid_features

    if len(ranking) != subset_size or len(np.unique(ranking)) != subset_size:
        raise RuntimeError("Apricot returned an invalid hybrid selection")

    output_selected = Path(args.output_selected_cells)
    output_metrics = Path(args.output_selected_metrics)
    output_done = Path(args.output_done)
    for output in (output_selected, output_metrics, output_done):
        output.parent.mkdir(parents=True, exist_ok=True)

    pd.Series(cell_ids[ranking]).to_csv(output_selected, index=False, header=False)
    pd.DataFrame(
        {
            "selected_rank": np.arange(1, subset_size + 1),
            "cell_id": cell_ids[ranking],
            "hybrid_selection_gain": gains,
            "sparsity_percentile": percentile[ranking],
            "density_weight": weights[ranking],
            "mean_knn_distance": mean_distance[ranking],
        }
    ).to_csv(output_metrics, sep="\t", index=False)
    output_done.write_text("DONE\n", encoding="utf-8")
    print(f"Selection completed in {elapsed:.2f} seconds", flush=True)
    print(f"Wrote selected cells: {output_selected}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
