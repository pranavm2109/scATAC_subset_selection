#!/usr/bin/env python3
"""Evaluate unweighted, density-only, and hybrid selections retrospectively."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ATLASES = ("kidney_atac", "fetal_atac")
FRACTIONS = ("02pct", "04pct", "06pct", "08pct", "10pct")
METHODS = ("unweighted", "density_only", "hybrid")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project_dir", default=".")
    parser.add_argument("--output_summary", required=True)
    parser.add_argument("--output_cell_types", required=True)
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


def selection_path(root: Path, atlas: str, fraction: str, method: str) -> Path:
    if method == "unweighted":
        return root / (
            f"results/subset_selection/peakvi/{atlas}/"
            f"selected_cells_peakvi_{fraction}.txt"
        )
    if method == "density_only":
        if atlas == "fetal_atac" and fraction == "02pct":
            return root / (
                "diagnostics/density_weighted_ablation/fetal_atac/02pct/"
                "selected_cells_density_max_2p0_02pct.txt"
            )
        return root / (
            "diagnostics/density_weighted_generalization/selections/"
            f"{atlas}/{fraction}/"
            f"selected_cells_density_max_2p0_{fraction}.txt"
        )
    if method == "hybrid":
        return root / (
            "diagnostics/hybrid_rbf_density_ablation/selections/"
            f"{atlas}/{fraction}/selected_cells_hybrid_{fraction}.txt"
        )
    raise ValueError(f"Unknown method: {method}")


def status(counts: np.ndarray) -> np.ndarray:
    return np.where(
        counts == 0,
        "absent",
        np.where(counts < 10, "underpowered", "evaluable"),
    )


def objective(
    features: np.ndarray,
    positions: np.ndarray,
    weights: np.ndarray | None,
    chunk_size: int,
) -> float:
    coverage = np.zeros(features.shape[1], dtype=np.float64)
    for start in range(0, len(positions), chunk_size):
        current_positions = positions[start : start + chunk_size]
        current = np.asarray(features[current_positions], dtype=np.float64)
        if weights is not None:
            current *= weights[current_positions, None]
        coverage += current.sum(axis=0)
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
    summary_rows: list[dict[str, object]] = []
    cell_type_rows: list[dict[str, object]] = []

    # Snakemake requires every selection before this evaluator starts. Labels
    # are therefore used only after all annotation-blind selections exist.
    for atlas in ATLASES:
        cell_ids = read_ids(root / f"results/peakvi/{atlas}/cell_names.txt")
        cell_index = pd.Index(cell_ids)
        metadata_path = root / f"results/cell_metadata/{atlas}/{atlas}_cell_metadata.tsv"
        metadata = pd.read_csv(metadata_path, sep="\t", dtype={"cell_id": str})
        required = {"cell_id", "cell_type"}
        missing = required - set(metadata.columns)
        if missing:
            raise ValueError(f"{metadata_path} is missing: {sorted(missing)}")
        labels = metadata.set_index("cell_id").reindex(cell_ids)["cell_type"]
        if labels.isna().any():
            raise ValueError(f"Some {atlas} cells lack cell_type metadata")
        labels = labels.astype(str).to_numpy()
        cell_types, codes = np.unique(labels, return_inverse=True)
        full_n = np.bincount(codes, minlength=len(cell_types))
        full_fraction = full_n / len(cell_ids)
        type_index = {cell_type: index for index, cell_type in enumerate(cell_types)}

        features = np.load(
            root / (
                f"results/peakvi_landmark_features/{atlas}/"
                "peakvi_rbf_landmark_features.npy"
            ),
            mmap_mode="r",
        )
        if features.ndim != 2 or features.shape[0] != len(cell_ids):
            raise ValueError(f"RBF features and cell IDs differ for {atlas}")
        with np.load(
            root / (
                f"diagnostics/peakvi_celltype_separability/{atlas}/"
                "per_cell_knn_sparsity.npz"
            ),
            allow_pickle=False,
        ) as sparsity:
            percentile = np.asarray(
                sparsity["mean_knn_distance_percentile"], dtype=np.float32
            )
            expected_checksum = str(sparsity["cell_name_sha256"].item())
        if len(percentile) != len(cell_ids) or expected_checksum != checksum(cell_ids):
            raise ValueError(f"Sparsity metrics have the wrong cell order for {atlas}")
        weights = np.float32(1.0) + np.square(percentile / np.float32(100.0))
        weights = np.minimum(weights, np.float32(2.0))

        for fraction in FRACTIONS:
            positions: dict[str, np.ndarray] = {}
            counts: dict[str, np.ndarray] = {}
            proportions: dict[str, np.ndarray] = {}
            original_objectives: dict[str, float] = {}
            density_objectives: dict[str, float] = {}
            for method in METHODS:
                path = selection_path(root, atlas, fraction, method)
                selected_ids = read_ids(path)
                current_positions = cell_index.get_indexer(selected_ids)
                if np.any(current_positions < 0):
                    raise ValueError(f"Unknown selected cells: {path}")
                positions[method] = current_positions.astype(np.int64)
                counts[method] = np.bincount(
                    codes[current_positions], minlength=len(cell_types)
                )
                proportions[method] = counts[method] / len(current_positions)
                original_objectives[method] = objective(
                    features,
                    positions[method],
                    None,
                    args.objective_chunk_size,
                )
                density_objectives[method] = objective(
                    features,
                    positions[method],
                    weights,
                    args.objective_chunk_size,
                )

            selected_sizes = {method: len(positions[method]) for method in METHODS}
            if len(set(selected_sizes.values())) != 1:
                raise ValueError(f"Selection sizes differ for {atlas} {fraction}")
            original_baseline = original_objectives["unweighted"]
            density_baseline = density_objectives["density_only"]
            unweighted_retained = counts["unweighted"] > 0
            density_retained = counts["density_only"] > 0
            unweighted_absent = counts["unweighted"] == 0

            for method in METHODS:
                current_status = status(counts[method])
                lost_vs_unweighted = cell_types[
                    unweighted_retained & (counts[method] == 0)
                ]
                lost_vs_density = cell_types[
                    density_retained & (counts[method] == 0)
                ]
                rescued_from_unweighted = cell_types[
                    unweighted_absent & (counts[method] > 0)
                ]
                summary_rows.append(
                    {
                        "atlas": atlas,
                        "subset_size": fraction,
                        "method": method,
                        "n_total_cell_types": len(cell_types),
                        "n_retained": int(np.sum(counts[method] > 0)),
                        "n_absent": int(np.sum(counts[method] == 0)),
                        "n_underpowered": int(
                            np.sum(current_status == "underpowered")
                        ),
                        "n_evaluable": int(np.sum(current_status == "evaluable")),
                        "pearson_proportion_correlation": safe_correlation(
                            full_fraction, proportions[method], "pearson"
                        ),
                        "spearman_proportion_correlation": safe_correlation(
                            full_fraction, proportions[method], "spearman"
                        ),
                        "original_unweighted_objective_retention": (
                            original_objectives[method] / original_baseline
                        ),
                        "density_weighted_objective_retention": (
                            density_objectives[method] / density_baseline
                        ),
                        "original_unweighted_objective": original_objectives[method],
                        "density_weighted_objective": density_objectives[method],
                        "Mes_selected_n": (
                            int(counts[method][type_index["Mes"]])
                            if atlas == "kidney_atac" and "Mes" in type_index
                            else np.nan
                        ),
                        "Macula_Densa_selected_n": (
                            int(counts[method][type_index["Macula_Densa"]])
                            if atlas == "kidney_atac" and "Macula_Densa" in type_index
                            else np.nan
                        ),
                        "cell_types_lost_vs_unweighted": ";".join(
                            lost_vs_unweighted
                        ),
                        "cell_types_lost_vs_density_only": ";".join(
                            lost_vs_density
                        ),
                        "unweighted_absent_cell_types_rescued": ";".join(
                            rescued_from_unweighted
                        ),
                    }
                )

            for index, cell_type in enumerate(cell_types):
                cell_type_rows.append(
                    {
                        "atlas": atlas,
                        "subset_size": fraction,
                        "cell_type": cell_type,
                        "full_n": int(full_n[index]),
                        "unweighted_selected_n": int(counts["unweighted"][index]),
                        "density_selected_n": int(counts["density_only"][index]),
                        "hybrid_selected_n": int(counts["hybrid"][index]),
                        "unweighted_status": status(counts["unweighted"])[index],
                        "density_status": status(counts["density_only"])[index],
                        "hybrid_status": status(counts["hybrid"])[index],
                    }
                )

    summary = pd.DataFrame(summary_rows)
    cell_type_table = pd.DataFrame(cell_type_rows)
    output_summary = Path(args.output_summary)
    output_cell_types = Path(args.output_cell_types)
    output_done = Path(args.output_done)
    for output in (output_summary, output_cell_types, output_done):
        output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_summary, sep="\t", index=False)
    cell_type_table.to_csv(output_cell_types, sep="\t", index=False)
    output_done.write_text("DONE\n", encoding="utf-8")
    print(f"Wrote hybrid method summary: {output_summary}", flush=True)
    print(f"Wrote hybrid cell-type details: {output_cell_types}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
