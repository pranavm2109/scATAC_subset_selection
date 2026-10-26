#!/usr/bin/env python3
"""Diagnose why density weighting trades Mes for Macula_Densa in kidney.

All latent-space, sparsity, landmark, RBF, singleton, and completed-set gain
quantities are computed before cell-type annotations are loaded. Labels are
joined only for the final retrospective Mes-versus-Macula_Densa summaries.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


FRACTIONS = ("02pct", "04pct", "06pct", "08pct", "10pct")
TARGETS = ("Mes", "Macula_Densa")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project_dir", default=".")
    parser.add_argument("--output_tsv", required=True)
    parser.add_argument("--output_done", required=True)
    parser.add_argument("--chunk_size", type=int, default=5000)
    parser.add_argument("--substantial_activation_threshold", type=float, default=0.5)
    parser.add_argument("--strongest_feature_count", type=int, default=10)
    return parser.parse_args()


def read_ids(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any():
        raise ValueError(f"Empty or missing cell IDs: {path}")
    result = values.to_numpy(dtype=str)
    if len(np.unique(result)) != len(result):
        raise ValueError(f"Duplicate cell IDs: {path}")
    return result


def checksum(cell_ids: np.ndarray) -> str:
    digest = hashlib.sha256()
    for cell_id in cell_ids:
        digest.update(str(cell_id).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def selected_positions(root: Path, cell_index: pd.Index, method: str) -> dict[str, np.ndarray]:
    output = {}
    for fraction in FRACTIONS:
        if method == "unweighted":
            path = root / (
                "results/subset_selection/peakvi/kidney_atac/"
                f"selected_cells_peakvi_{fraction}.txt"
            )
        else:
            path = root / (
                "diagnostics/density_weighted_generalization/selections/"
                f"kidney_atac/{fraction}/"
                f"selected_cells_density_max_2p0_{fraction}.txt"
            )
        ids = read_ids(path)
        positions = cell_index.get_indexer(ids)
        if np.any(positions < 0):
            raise ValueError(f"Unknown selected cells: {path}")
        output[fraction] = positions.astype(np.int64)
    return output


def intrinsic_rbf_metrics(
    features: np.ndarray,
    weights: np.ndarray,
    sigma: float,
    threshold: float,
    chunk_size: int,
) -> dict[str, np.ndarray]:
    n_cells = features.shape[0]
    output = {
        "max_activation": np.empty(n_cells, dtype=np.float32),
        "nearest_landmark_distance": np.empty(n_cells, dtype=np.float32),
        "n_substantial_activations": np.empty(n_cells, dtype=np.int32),
        "unweighted_singleton_score": np.empty(n_cells, dtype=np.float32),
        "density_weighted_singleton_score": np.empty(n_cells, dtype=np.float32),
    }
    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        current = np.asarray(features[start:stop], dtype=np.float32)
        maximum = current.max(axis=1)
        output["max_activation"][start:stop] = maximum
        output["nearest_landmark_distance"][start:stop] = (
            sigma
            * np.sqrt(-2.0 * np.log(np.clip(maximum, 1e-30, 1.0)))
        )
        output["n_substantial_activations"][start:stop] = (
            current >= threshold
        ).sum(axis=1)
        singleton = np.sqrt(current).sum(axis=1)
        output["unweighted_singleton_score"][start:stop] = singleton
        output["density_weighted_singleton_score"][start:stop] = (
            singleton * np.sqrt(weights[start:stop])
        )
    output["nearest_landmark_distance_percentile"] = (
        pd.Series(output["nearest_landmark_distance"])
        .rank(method="average", pct=True)
        .to_numpy(dtype=np.float32)
        * 100
    )
    return output


def completed_gain_percentiles(
    features: np.ndarray,
    weights: np.ndarray | None,
    selections: dict[str, np.ndarray],
    chunk_size: int,
) -> dict[str, np.ndarray]:
    output = {}
    n_cells, n_features = features.shape
    for fraction, positions in selections.items():
        coverage = np.zeros(n_features, dtype=np.float64)
        for start in range(0, len(positions), chunk_size):
            current_positions = positions[start : start + chunk_size]
            current = np.asarray(features[current_positions], dtype=np.float64)
            if weights is not None:
                current *= weights[current_positions, None]
            coverage += current.sum(axis=0)
        baseline = np.sqrt(coverage).sum()
        gains = np.empty(n_cells, dtype=np.float64)
        print(f"Computing {fraction} completed-set residual gains", flush=True)
        for start in range(0, n_cells, chunk_size):
            stop = min(start + chunk_size, n_cells)
            current = np.asarray(features[start:stop], dtype=np.float64)
            if weights is not None:
                current *= weights[start:stop, None]
            gains[start:stop] = np.sqrt(current + coverage).sum(axis=1) - baseline
        gains[positions] = np.nan
        unselected = np.isfinite(gains)
        percentile = np.full(n_cells, np.nan, dtype=np.float32)
        percentile[unselected] = (
            pd.Series(gains[unselected])
            .rank(method="average", pct=True)
            .to_numpy(dtype=np.float32)
            * 100
        )
        output[fraction] = percentile
    return output


def quantiles(values: np.ndarray, prefix: str) -> dict[str, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            f"{prefix}_mean": np.nan,
            f"{prefix}_q10": np.nan,
            f"{prefix}_q25": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_q75": np.nan,
            f"{prefix}_q90": np.nan,
        }
    q10, q25, median, q75, q90 = np.quantile(values, (0.10, 0.25, 0.50, 0.75, 0.90))
    return {
        f"{prefix}_mean": float(np.mean(values)),
        f"{prefix}_q10": float(q10),
        f"{prefix}_q25": float(q25),
        f"{prefix}_median": float(median),
        f"{prefix}_q75": float(q75),
        f"{prefix}_q90": float(q90),
    }


def main() -> None:
    args = parse_args()
    if args.chunk_size < 1 or args.strongest_feature_count < 1:
        raise ValueError("Chunk and strongest-feature counts must be positive")
    if not 0 < args.substantial_activation_threshold <= 1:
        raise ValueError("Activation threshold must be in (0, 1]")
    root = Path(args.project_dir).resolve()
    cell_ids = read_ids(root / "results/peakvi/kidney_atac/cell_names.txt")
    cell_index = pd.Index(cell_ids)

    features = np.load(
        root / (
            "results/peakvi_landmark_features/kidney_atac/"
            "peakvi_rbf_landmark_features.npy"
        ),
        mmap_mode="r",
    )
    latent = np.load(
        root / "results/peakvi/kidney_atac/peakvi_latent.npy", mmap_mode="r"
    )
    centroids = np.load(
        root / (
            "results/peakvi_landmark_features/kidney_atac/"
            "peakvi_landmark_centroids.npy"
        ),
        mmap_mode="r",
    )
    if features.shape[0] != len(cell_ids) or latent.shape[0] != len(cell_ids):
        raise ValueError("Kidney cell order does not match latent/RBF rows")
    if latent.shape[1] != centroids.shape[1] or features.shape[1] != len(centroids):
        raise ValueError("PeakVI, landmark, and RBF dimensions do not align")
    with open(
        root / (
            "results/peakvi_landmark_features/kidney_atac/"
            "peakvi_landmark_metadata.json"
        ),
        encoding="utf-8",
    ) as handle:
        landmark_metadata = json.load(handle)
    sigma = float(landmark_metadata["rbf_sigma"])

    with np.load(
        root / (
            "diagnostics/peakvi_celltype_separability/kidney_atac/"
            "per_cell_knn_sparsity.npz"
        ),
        allow_pickle=False,
    ) as sparsity:
        sparsity_percentile = np.asarray(
            sparsity["mean_knn_distance_percentile"], dtype=np.float32
        )
        expected_checksum = str(sparsity["cell_name_sha256"].item())
    if len(sparsity_percentile) != len(cell_ids) or expected_checksum != checksum(cell_ids):
        raise ValueError("Kidney sparsity metrics have the wrong cell order")
    weights = np.float32(1.0) + np.square(
        sparsity_percentile / np.float32(100.0)
    )
    weights = np.minimum(weights, np.float32(2.0))

    unweighted_selections = selected_positions(root, cell_index, "unweighted")
    density_selections = selected_positions(root, cell_index, "density")
    intrinsic = intrinsic_rbf_metrics(
        features,
        weights,
        sigma,
        args.substantial_activation_threshold,
        args.chunk_size,
    )
    unweighted_residual = completed_gain_percentiles(
        features, None, unweighted_selections, args.chunk_size
    )
    density_residual = completed_gain_percentiles(
        features, weights, density_selections, args.chunk_size
    )

    # All annotation-blind quantities above are complete before this join.
    metadata_path = root / (
        "results/cell_metadata/kidney_atac/kidney_atac_cell_metadata.tsv"
    )
    metadata = pd.read_csv(metadata_path, sep="\t", dtype={"cell_id": str})
    labels = metadata.set_index("cell_id").reindex(cell_ids)["cell_type"]
    if labels.isna().any():
        raise ValueError("Some kidney cells lack cell_type metadata")
    labels = labels.astype(str).to_numpy()
    separability = pd.read_csv(
        root / (
            "diagnostics/peakvi_celltype_separability/kidney_atac/"
            "cell_type_separability.tsv"
        ),
        sep="\t",
    ).set_index("cell_type")

    rows = []
    for target in TARGETS:
        target_positions = np.flatnonzero(labels == target)
        if not len(target_positions):
            raise ValueError(f"Required kidney cell type is absent: {target}")
        if target not in separability.index:
            raise ValueError(f"Missing kNN purity for: {target}")
        target_features = np.asarray(features[target_positions], dtype=np.float64)
        feature_means = target_features.mean(axis=0)
        strongest = np.argsort(-feature_means)[: args.strongest_feature_count]
        row = {
            "cell_type": target,
            "full_atlas_n": len(target_positions),
            "gaussian_peakvi_knn_purity": float(
                separability.loc[target, "knn_purity"]
            ),
            "max_weight": 2.0,
            "substantial_rbf_activation_threshold": (
                args.substantial_activation_threshold
            ),
            "strongest_rbf_landmark_features": ";".join(
                f"RBF_landmark_{index + 1}" for index in strongest
            ),
            "strongest_rbf_landmark_mean_activations": ";".join(
                f"{feature_means[index]:.8g}" for index in strongest
            ),
        }
        row.update(
            quantiles(
                sparsity_percentile[target_positions], "sparsity_percentile"
            )
        )
        row.update(quantiles(weights[target_positions], "density_weight"))
        row.update(
            quantiles(
                intrinsic["nearest_landmark_distance"][target_positions],
                "nearest_landmark_distance",
            )
        )
        row.update(
            quantiles(
                intrinsic["nearest_landmark_distance_percentile"][target_positions],
                "nearest_landmark_distance_percentile",
            )
        )
        row.update(
            quantiles(
                intrinsic["max_activation"][target_positions],
                "maximum_rbf_activation",
            )
        )
        row.update(
            quantiles(
                intrinsic["n_substantial_activations"][target_positions].astype(float),
                "n_landmarks_with_substantial_activation",
            )
        )
        row.update(
            quantiles(
                intrinsic["unweighted_singleton_score"][target_positions],
                "unweighted_singleton_score",
            )
        )
        row.update(
            quantiles(
                intrinsic["density_weighted_singleton_score"][target_positions],
                "density_weighted_singleton_score",
            )
        )
        for fraction in FRACTIONS:
            row[f"unweighted_selected_n_{fraction}"] = int(
                np.isin(unweighted_selections[fraction], target_positions).sum()
            )
            row[f"density_max_2p0_selected_n_{fraction}"] = int(
                np.isin(density_selections[fraction], target_positions).sum()
            )
            row.update(
                quantiles(
                    unweighted_residual[fraction][target_positions],
                    f"unweighted_residual_gain_percentile_{fraction}",
                )
            )
            row.update(
                quantiles(
                    density_residual[fraction][target_positions],
                    f"density_max_2p0_residual_gain_percentile_{fraction}",
                )
            )
        rows.append(row)

    output = Path(args.output_tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    done = Path(args.output_done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("DONE\n", encoding="utf-8")
    print(f"Wrote focused Mes/Macula_Densa diagnostic: {output}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
