#!/usr/bin/env python3
"""Summarize the frozen max-weight-2 density ablation retrospectively.

Selection files are already complete before this script loads cell-type
annotations. It compares unweighted and density-weighted selections for fetal
and kidney at 2%, 4%, 6%, 8%, and 10%, including the original unweighted RBF
objective and annotation-blind selected-cell sparsity summaries.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ATLASES = ("fetal_atac", "kidney_atac")
FRACTIONS = ("02pct", "04pct", "06pct", "08pct", "10pct")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project_dir", default=".")
    parser.add_argument("--output_summary", required=True)
    parser.add_argument("--output_cell_types", required=True)
    parser.add_argument("--output_changes", required=True)
    parser.add_argument("--output_done", required=True)
    parser.add_argument("--objective_chunk_size", type=int, default=10000)
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


def density_path(root: Path, atlas: str, fraction: str) -> Path:
    if atlas == "fetal_atac" and fraction == "02pct":
        return root / (
            "diagnostics/density_weighted_ablation/fetal_atac/02pct/"
            "selected_cells_density_max_2p0_02pct.txt"
        )
    return root / (
        f"diagnostics/density_weighted_generalization/selections/{atlas}/{fraction}/"
        f"selected_cells_density_max_2p0_{fraction}.txt"
    )


def status(count: int) -> str:
    if count == 0:
        return "absent"
    if count < 10:
        return "underpowered"
    return "evaluable"


def change_class(unweighted_n: int, density_n: int) -> str:
    if unweighted_n == 0 and density_n > 0:
        return "rescued"
    if unweighted_n > 0 and density_n == 0:
        return "lost"
    if 1 <= unweighted_n <= 9 and density_n >= 10:
        return "improved_underpowered_to_evaluable"
    if unweighted_n >= 10 and 1 <= density_n <= 9:
        return "worsened_evaluable_to_underpowered"
    return "unchanged"


def safe_correlation(full: np.ndarray, selected: np.ndarray, method: str) -> float:
    if len(full) < 2 or np.std(full) == 0 or np.std(selected) == 0:
        return float("nan")
    return float(pd.Series(full).corr(pd.Series(selected), method=method))


def objective(features: np.ndarray, positions: np.ndarray, chunk_size: int) -> float:
    coverage = np.zeros(features.shape[1], dtype=np.float64)
    for start in range(0, len(positions), chunk_size):
        current = positions[start : start + chunk_size]
        coverage += np.asarray(features[current], dtype=np.float64).sum(axis=0)
    return float(np.sqrt(coverage).sum())


def selected_summary(
    method: str,
    selected_n: np.ndarray,
    full_fraction: np.ndarray,
    selected_fraction: np.ndarray,
    selected_positions: np.ndarray,
    sparsity_percentile: np.ndarray,
    unweighted_objective: float,
    baseline_objective: float,
) -> dict[str, object]:
    statuses = np.asarray([status(int(value)) for value in selected_n])
    selected_sparsity = sparsity_percentile[selected_positions]
    retention = unweighted_objective / baseline_objective
    return {
        "method": method,
        "n_total_cell_types": len(selected_n),
        "n_retained": int(np.sum(selected_n > 0)),
        "n_absent": int(np.sum(selected_n == 0)),
        "n_underpowered": int(np.sum(statuses == "underpowered")),
        "n_evaluable": int(np.sum(statuses == "evaluable")),
        "pearson_proportion_correlation": safe_correlation(
            full_fraction, selected_fraction, "pearson"
        ),
        "spearman_proportion_correlation": safe_correlation(
            full_fraction, selected_fraction, "spearman"
        ),
        "original_unweighted_objective": unweighted_objective,
        "objective_retention_vs_unweighted": retention,
        "objective_retention_vs_unweighted_percent": 100.0 * retention,
        "median_selected_sparsity_percentile": float(np.median(selected_sparsity)),
        "mean_selected_sparsity_percentile": float(np.mean(selected_sparsity)),
        "proportion_selected_above_75th_sparsity_percentile": float(
            np.mean(selected_sparsity > 75)
        ),
        "proportion_selected_above_90th_sparsity_percentile": float(
            np.mean(selected_sparsity > 90)
        ),
    }


def main() -> None:
    args = parse_args()
    if args.objective_chunk_size < 1:
        raise ValueError("--objective_chunk_size must be >= 1")
    root = Path(args.project_dir).resolve()
    summary_rows: list[dict[str, object]] = []
    cell_type_rows: list[dict[str, object]] = []
    change_rows: list[dict[str, object]] = []

    # All selections are file inputs and are complete before annotations are
    # loaded here. Labels are used only for retrospective evaluation.
    for atlas in ATLASES:
        cell_ids = read_ids(root / f"results/peakvi/{atlas}/cell_names.txt")
        index = pd.Index(cell_ids)
        metadata_path = root / f"results/cell_metadata/{atlas}/{atlas}_cell_metadata.tsv"
        metadata = pd.read_csv(metadata_path, sep="\t", dtype={"cell_id": str})
        required = {"cell_id", "cell_type"}
        missing = required - set(metadata.columns)
        if missing:
            raise ValueError(f"{metadata_path} is missing columns: {sorted(missing)}")
        if metadata["cell_id"].duplicated().any():
            raise ValueError(f"Duplicate metadata cell IDs: {metadata_path}")
        labels = metadata.set_index("cell_id").reindex(cell_ids)["cell_type"]
        if labels.isna().any():
            raise ValueError(f"Some {atlas} cells lack cell_type metadata")
        labels = labels.astype(str).to_numpy()
        cell_types, codes = np.unique(labels, return_inverse=True)
        full_n = np.bincount(codes, minlength=len(cell_types))
        full_fraction = full_n / len(cell_ids)

        features = np.load(
            root / (
                f"results/peakvi_landmark_features/{atlas}/"
                "peakvi_rbf_landmark_features.npy"
            ),
            mmap_mode="r",
        )
        if features.ndim != 2 or features.shape[0] != len(cell_ids):
            raise ValueError(f"Feature matrix and cell IDs differ for {atlas}")
        sparsity_path = root / (
            f"diagnostics/peakvi_celltype_separability/{atlas}/"
            "per_cell_knn_sparsity.npz"
        )
        with np.load(sparsity_path, allow_pickle=False) as sparsity:
            percentile = np.asarray(
                sparsity["mean_knn_distance_percentile"], dtype=np.float32
            )
            expected_checksum = str(sparsity["cell_name_sha256"].item())
        if len(percentile) != len(cell_ids) or expected_checksum != checksum(cell_ids):
            raise ValueError(f"Sparsity metrics have the wrong cell order for {atlas}")

        for fraction in FRACTIONS:
            paths = {
                "unweighted": root / (
                    f"results/subset_selection/peakvi/{atlas}/"
                    f"selected_cells_peakvi_{fraction}.txt"
                ),
                "density_max_2p0": density_path(root, atlas, fraction),
            }
            positions: dict[str, np.ndarray] = {}
            counts: dict[str, np.ndarray] = {}
            proportions: dict[str, np.ndarray] = {}
            objectives: dict[str, float] = {}
            selected_sizes: dict[str, int] = {}
            for method, path in paths.items():
                selected_ids = read_ids(path)
                selected_positions = index.get_indexer(selected_ids)
                if np.any(selected_positions < 0):
                    raise ValueError(f"Unknown cells in {path}")
                positions[method] = selected_positions.astype(np.int64)
                selected_sizes[method] = len(selected_positions)
                counts[method] = np.bincount(
                    codes[selected_positions], minlength=len(cell_types)
                )
                proportions[method] = counts[method] / len(selected_positions)
                objectives[method] = objective(
                    features, positions[method], args.objective_chunk_size
                )

            if selected_sizes["unweighted"] != selected_sizes["density_max_2p0"]:
                raise ValueError(f"Selection sizes differ for {atlas} {fraction}")
            baseline_objective = objectives["unweighted"]
            for method in ("unweighted", "density_max_2p0"):
                row = selected_summary(
                    method,
                    counts[method],
                    full_fraction,
                    proportions[method],
                    positions[method],
                    percentile,
                    objectives[method],
                    baseline_objective,
                )
                row.update(
                    {
                        "atlas": atlas,
                        "subset_size": fraction,
                        "n_selected": selected_sizes[method],
                    }
                )
                summary_rows.append(row)

            classes = []
            for cell_type_index, cell_type in enumerate(cell_types):
                unweighted_n = int(counts["unweighted"][cell_type_index])
                density_n = int(counts["density_max_2p0"][cell_type_index])
                current_class = change_class(unweighted_n, density_n)
                classes.append(current_class)
                cell_type_rows.append(
                    {
                        "atlas": atlas,
                        "subset_size": fraction,
                        "cell_type": cell_type,
                        "full_n": int(full_n[cell_type_index]),
                        "unweighted_selected_n": unweighted_n,
                        "density_selected_n": density_n,
                        "unweighted_status": status(unweighted_n),
                        "density_status": status(density_n),
                        "change_class": current_class,
                    }
                )

            classes = np.asarray(classes)
            change_rows.append(
                {
                    "atlas": atlas,
                    "subset_size": fraction,
                    "populations_rescued": ";".join(
                        cell_types[classes == "rescued"]
                    ),
                    "populations_lost": ";".join(cell_types[classes == "lost"]),
                    "populations_promoted_underpowered_to_evaluable": ";".join(
                        cell_types[classes == "improved_underpowered_to_evaluable"]
                    ),
                    "populations_demoted_evaluable_to_underpowered": ";".join(
                        cell_types[classes == "worsened_evaluable_to_underpowered"]
                    ),
                }
            )

    outputs = {
        Path(args.output_summary): pd.DataFrame(summary_rows),
        Path(args.output_cell_types): pd.DataFrame(cell_type_rows),
        Path(args.output_changes): pd.DataFrame(change_rows),
    }
    for path, frame in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, sep="\t", index=False)
        print(f"Wrote: {path}", flush=True)
    done = Path(args.output_done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("DONE\n", encoding="utf-8")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
