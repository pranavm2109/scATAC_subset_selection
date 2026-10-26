#!/usr/bin/env python3
"""Retrospectively evaluate the kidney-only density-cap ablation.

All selected-cell files are completed before this script loads cell-type
annotations. The original unweighted RBF objective is evaluated on every set.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


FRACTIONS = ("02pct", "04pct", "06pct", "08pct", "10pct")
CAPS = ("1p5", "1p6", "1p7", "1p8", "1p9", "2p0")
CAP_VALUES = {label: float(label.replace("p", ".")) for label in CAPS}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project_dir", default=".")
    parser.add_argument("--output_summary", required=True)
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


def cap_selection_path(root: Path, cap: str, fraction: str) -> Path:
    if cap == "2p0":
        return root / (
            "diagnostics/density_weighted_generalization/selections/"
            f"kidney_atac/{fraction}/selected_cells_density_max_2p0_{fraction}.txt"
        )
    return root / (
        "diagnostics/kidney_density_cap_ablation/selections/"
        f"{cap}/{fraction}/selected_cells_density_max_{cap}_{fraction}.txt"
    )


def status(counts: np.ndarray) -> np.ndarray:
    return np.where(
        counts == 0,
        "absent",
        np.where(counts < 10, "underpowered", "evaluable"),
    )


def objective(features: np.ndarray, positions: np.ndarray, chunk_size: int) -> float:
    coverage = np.zeros(features.shape[1], dtype=np.float64)
    for start in range(0, len(positions), chunk_size):
        current = positions[start : start + chunk_size]
        coverage += np.asarray(features[current], dtype=np.float64).sum(axis=0)
    return float(np.sqrt(coverage).sum())


def safe_correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(pd.Series(x).corr(pd.Series(y), method=method))


def main() -> None:
    args = parse_args()
    if args.objective_chunk_size < 1:
        raise ValueError("--objective_chunk_size must be >= 1")
    root = Path(args.project_dir).resolve()
    cell_ids = read_ids(root / "results/peakvi/kidney_atac/cell_names.txt")
    cell_index = pd.Index(cell_ids)

    # Selection is complete before labels are loaded below.
    metadata_path = root / (
        "results/cell_metadata/kidney_atac/kidney_atac_cell_metadata.tsv"
    )
    metadata = pd.read_csv(metadata_path, sep="\t", dtype={"cell_id": str})
    required = {"cell_id", "cell_type"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    labels = metadata.set_index("cell_id").reindex(cell_ids)["cell_type"]
    if labels.isna().any():
        raise ValueError("Some kidney cells lack cell_type metadata")
    labels = labels.astype(str).to_numpy()
    cell_types, codes = np.unique(labels, return_inverse=True)
    full_n = np.bincount(codes, minlength=len(cell_types))
    full_fraction = full_n / len(cell_ids)
    type_index = {cell_type: index for index, cell_type in enumerate(cell_types)}
    for target in ("Mes", "Macula_Densa"):
        if target not in type_index:
            raise ValueError(f"Required kidney cell type is absent: {target}")

    features = np.load(
        root / (
            "results/peakvi_landmark_features/kidney_atac/"
            "peakvi_rbf_landmark_features.npy"
        ),
        mmap_mode="r",
    )
    if features.ndim != 2 or features.shape[0] != len(cell_ids):
        raise ValueError("Kidney RBF features and cell names do not align")

    rows: list[dict[str, object]] = []
    for fraction in FRACTIONS:
        baseline_path = root / (
            "results/subset_selection/peakvi/kidney_atac/"
            f"selected_cells_peakvi_{fraction}.txt"
        )
        baseline_ids = read_ids(baseline_path)
        baseline_positions = cell_index.get_indexer(baseline_ids)
        if np.any(baseline_positions < 0):
            raise ValueError(f"Unknown baseline cells: {baseline_path}")
        baseline_counts = np.bincount(
            codes[baseline_positions], minlength=len(cell_types)
        )
        baseline_status = status(baseline_counts)
        baseline_objective = objective(
            features, baseline_positions, args.objective_chunk_size
        )

        for cap in CAPS:
            selected_path = cap_selection_path(root, cap, fraction)
            selected_ids = read_ids(selected_path)
            selected_positions = cell_index.get_indexer(selected_ids)
            if np.any(selected_positions < 0):
                raise ValueError(f"Unknown selected cells: {selected_path}")
            if len(selected_positions) != len(baseline_positions):
                raise ValueError(f"Subset size mismatch: {selected_path}")
            selected_counts = np.bincount(
                codes[selected_positions], minlength=len(cell_types)
            )
            selected_status = status(selected_counts)
            selected_fraction = selected_counts / len(selected_positions)
            newly_lost = cell_types[(baseline_counts > 0) & (selected_counts == 0)]
            newly_rescued = cell_types[(baseline_counts == 0) & (selected_counts > 0)]
            demoted = cell_types[
                (baseline_status == "evaluable")
                & (selected_status == "underpowered")
            ]
            current_objective = objective(
                features, selected_positions, args.objective_chunk_size
            )
            rows.append(
                {
                    "max_weight": CAP_VALUES[cap],
                    "subset_size": fraction,
                    "n_retained": int(np.sum(selected_counts > 0)),
                    "n_absent": int(np.sum(selected_counts == 0)),
                    "n_underpowered": int(np.sum(selected_status == "underpowered")),
                    "n_evaluable": int(np.sum(selected_status == "evaluable")),
                    "Mes_selected_n": int(selected_counts[type_index["Mes"]]),
                    "Macula_Densa_selected_n": int(
                        selected_counts[type_index["Macula_Densa"]]
                    ),
                    "absent_cell_types": ";".join(
                        cell_types[selected_counts == 0]
                    ),
                    "newly_lost_cell_types": ";".join(newly_lost),
                    "newly_rescued_cell_types": ";".join(newly_rescued),
                    "evaluable_to_underpowered_cell_types": ";".join(demoted),
                    "unweighted_objective_retention": (
                        current_objective / baseline_objective
                    ),
                    "pearson_proportion_correlation": safe_correlation(
                        full_fraction, selected_fraction, "pearson"
                    ),
                    "spearman_proportion_correlation": safe_correlation(
                        full_fraction, selected_fraction, "spearman"
                    ),
                    "perfect_retention": bool(np.all(selected_counts > 0)),
                    "retains_mes_and_macula_densa": bool(
                        selected_counts[type_index["Mes"]] > 0
                        and selected_counts[type_index["Macula_Densa"]] > 0
                    ),
                    "no_newly_lost_cell_types": len(newly_lost) == 0,
                }
            )

    summary = pd.DataFrame(rows).sort_values(
        ["max_weight", "subset_size"]
    ).reset_index(drop=True)
    success_by_cap = summary.groupby("max_weight", sort=True).apply(
        lambda frame: bool(
            frame["perfect_retention"].all()
            and frame["retains_mes_and_macula_densa"].all()
            and frame["no_newly_lost_cell_types"].all()
        ),
    )
    successful_caps = success_by_cap[success_by_cap].index.to_numpy(dtype=float)
    smallest = float(successful_caps.min()) if len(successful_caps) else np.nan
    summary["cap_success_all_subset_sizes"] = summary["max_weight"].map(
        success_by_cap
    )
    summary["smallest_successful_cap_all_subset_sizes"] = smallest

    output = Path(args.output_summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, sep="\t", index=False)
    done = Path(args.output_done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("DONE\n", encoding="utf-8")
    print(f"Wrote kidney density-cap summary: {output}", flush=True)
    print(
        "Smallest successful cap across all subset sizes: "
        + ("none" if not np.isfinite(smallest) else f"{smallest:g}"),
        flush=True,
    )
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
