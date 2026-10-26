#!/usr/bin/env python3
"""Run a controlled, annotation-blind density-weighted selection ablation.

The existing unweighted PeakVI-RBF-Apricot subset is the baseline. New variants
multiply each cell's nonnegative RBF feature row by a bounded weight derived
only from its percentile-ranked mean kNN distance in Gaussian PeakVI space:

    weight = 1 + (max_weight - 1) * (sparsity_percentile / 100)^2

Cell-type annotations are loaded only after every weighted selection is done.
They are used solely for retrospective evaluation, never for weighting or
selection. This is a temporary diagnostic and does not alter production files.
"""

from __future__ import annotations

import argparse
import hashlib
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from apricot import FeatureBasedSelection


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare baseline selection with capped density weighting."
    )
    parser.add_argument("--feature_matrix", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--sparsity_npz", required=True)
    parser.add_argument("--baseline_selected_cells", required=True)
    parser.add_argument("--metadata_tsv", required=True)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--fraction", type=float, required=True)
    parser.add_argument("--fraction_label", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--max_weights", default="1.25,1.5,2.0")
    parser.add_argument("--rare_threshold", type=float, default=0.01)
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


def parse_max_weights(value: str) -> list[float]:
    weights = [float(item) for item in value.split(",")]
    if not weights or any(not np.isfinite(item) or item <= 1 for item in weights):
        raise ValueError("--max_weights must contain comma-separated values > 1")
    if len(set(weights)) != len(weights):
        raise ValueError("--max_weights contains duplicate values")
    return sorted(weights)


def method_name(max_weight: float) -> str:
    return "density_max_" + str(max_weight).replace(".", "p")


def load_inputs(args: argparse.Namespace):
    cell_ids = read_ids(Path(args.cell_names))
    features = np.load(args.feature_matrix, mmap_mode="r")
    if features.ndim != 2 or features.shape[0] != len(cell_ids):
        raise ValueError(
            f"Feature matrix shape {features.shape} does not match "
            f"{len(cell_ids):,} cell IDs"
        )
    if not np.issubdtype(features.dtype, np.number):
        raise ValueError("RBF feature matrix must be numeric")
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
    positions = pd.Index(cell_ids).get_indexer(baseline_ids)
    if np.any(positions < 0):
        raise ValueError("Baseline selected cells are missing from the feature matrix")

    expected_subset_size = max(int(round(len(cell_ids) * args.fraction)), 1)
    if len(baseline_ids) != expected_subset_size:
        raise ValueError(
            f"Baseline contains {len(baseline_ids):,} cells; "
            f"expected {expected_subset_size:,}"
        )
    return (
        cell_ids,
        features,
        mean_distance,
        percentile,
        knn_k,
        positions.astype(np.int64),
        expected_subset_size,
    )


def selection_objective(features: np.ndarray, positions: np.ndarray) -> float:
    coverage = np.asarray(features[positions], dtype=np.float64).sum(axis=0)
    return float(np.sqrt(coverage).sum())


def run_weighted_selection(
    features: np.ndarray,
    sparsity_percentile: np.ndarray,
    subset_size: int,
    max_weight: float,
    random_seed: int,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    # Quadratic scaling concentrates the intervention in the sparse tail while
    # retaining a strict, interpretable maximum influence per cell.
    scaled = sparsity_percentile.astype(np.float32) / np.float32(100.0)
    weights = np.float32(1.0) + np.float32(max_weight - 1.0) * np.square(scaled)
    weights = np.minimum(weights, np.float32(max_weight))

    print(
        f"Running density variant max_weight={max_weight:g}; "
        f"weight min/median/max={weights.min():.4f}/"
        f"{np.median(weights):.4f}/{weights.max():.4f}",
        flush=True,
    )
    start = time.time()
    weighted_features = np.empty(features.shape, dtype=np.float32)
    np.multiply(features, weights[:, None], out=weighted_features)
    selector = FeatureBasedSelection(
        subset_size,
        "sqrt",
        optimizer="two-stage",
        random_state=random_seed,
        verbose=True,
    )
    selector.fit(weighted_features)
    elapsed = time.time() - start
    ranking = np.asarray(selector.ranking, dtype=np.int64)
    gains = np.asarray(selector.gains, dtype=np.float64)
    del weighted_features

    if len(ranking) != subset_size or len(np.unique(ranking)) != subset_size:
        raise RuntimeError("Apricot returned an invalid weighted selection")
    return ranking, gains, elapsed, weights


def load_labels_after_selection(path: Path, cell_ids: np.ndarray) -> np.ndarray:
    metadata = pd.read_csv(path, sep="\t", dtype={"cell_id": str})
    required = {"cell_id", "cell_type"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    if metadata["cell_id"].duplicated().any():
        raise ValueError("Metadata cell IDs are not unique")
    metadata = metadata.set_index("cell_id").reindex(cell_ids)
    if metadata["cell_type"].isna().any():
        raise ValueError("Some feature-matrix cells lack cell_type metadata")
    return metadata["cell_type"].astype(str).to_numpy()


def evaluate(
    atlas: str,
    fraction_label: str,
    labels: np.ndarray,
    selections: dict[str, np.ndarray],
    sparsity_percentile: np.ndarray,
    objectives: dict[str, float],
    elapsed: dict[str, float],
    rare_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cell_types, codes = np.unique(labels, return_inverse=True)
    full_n = np.bincount(codes, minlength=len(cell_types))
    full_fraction = full_n / len(labels)
    rare = full_fraction < rare_threshold
    baseline_retained = None
    detail_rows = []
    method_counts = {}

    for method, positions in selections.items():
        selected_n = np.bincount(codes[positions], minlength=len(cell_types))
        method_counts[method] = selected_n
        retained = selected_n > 0
        if method == "baseline":
            baseline_retained = retained
        for code, cell_type in enumerate(cell_types):
            detail_rows.append(
                {
                    "atlas": atlas,
                    "subset_size": fraction_label,
                    "method": method,
                    "cell_type": cell_type,
                    "full_n": int(full_n[code]),
                    "full_fraction": full_fraction[code],
                    "is_rare": bool(rare[code]),
                    "selected_n": int(selected_n[code]),
                    "selected_fraction": selected_n[code] / len(positions),
                    "retained": bool(retained[code]),
                    "enrichment": (
                        (selected_n[code] / len(positions)) / full_fraction[code]
                    ),
                }
            )

    if baseline_retained is None:
        raise RuntimeError("Baseline evaluation is missing")

    summary_rows = []
    baseline_objective = objectives["baseline"]
    for method, positions in selections.items():
        selected_n = method_counts[method]
        retained = selected_n > 0
        rescued = cell_types[retained & ~baseline_retained]
        lost = cell_types[~retained & baseline_retained]
        selected_sparsity = sparsity_percentile[positions]
        summary_rows.append(
            {
                "atlas": atlas,
                "subset_size": fraction_label,
                "method": method,
                "n_selected": len(positions),
                "n_cell_types_retained": int(retained.sum()),
                "n_cell_types_total": len(cell_types),
                "n_rare_cell_types_retained": int((retained & rare).sum()),
                "n_rare_cell_types_total": int(rare.sum()),
                "n_cell_types_rescued_vs_baseline": len(rescued),
                "n_cell_types_lost_vs_baseline": len(lost),
                "rescued_cell_types": ";".join(rescued),
                "lost_cell_types": ";".join(lost),
                "unweighted_rbf_objective": objectives[method],
                "objective_relative_to_baseline": objectives[method] / baseline_objective,
                "selected_sparsity_percentile_median": float(
                    np.median(selected_sparsity)
                ),
                "selected_sparsity_percentile_q90": float(
                    np.quantile(selected_sparsity, 0.90)
                ),
                "selection_time_seconds": elapsed[method],
            }
        )
    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


def plot_results(
    summary: pd.DataFrame,
    selected_sparsity: dict[str, np.ndarray],
    output_png: Path,
) -> None:
    methods = summary["method"].tolist()
    labels = [
        "Baseline" if method == "baseline" else method.replace("density_max_", "Max ").replace("p", ".")
        for method in methods
    ]
    x = np.arange(len(methods))
    figure, axes = plt.subplots(1, 3, figsize=(16, 5.2), dpi=220)

    axes[0].bar(x - 0.18, summary["n_cell_types_retained"], width=0.36, label="All types")
    axes[0].bar(
        x + 0.18,
        summary["n_rare_cell_types_retained"],
        width=0.36,
        label="Rare types",
    )
    axes[0].set_ylabel("Cell types retained")
    axes[0].legend(frameon=False)

    axes[1].bar(x, summary["objective_relative_to_baseline"], color="#4C78A8")
    axes[1].axhline(1, color="black", linestyle="--", linewidth=1)
    axes[1].set_ylabel("Unweighted RBF objective / baseline")

    axes[2].boxplot(
        [selected_sparsity[method] for method in methods],
        labels=labels,
        showfliers=False,
    )
    axes[2].set_ylabel("Selected-cell kNN sparsity percentile")

    for axis in axes[:2]:
        axis.set_xticks(x, labels, rotation=25, ha="right")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#E4E4E4", linewidth=0.6)
    figure.suptitle(
        f"{summary.loc[0, 'atlas']} {summary.loc[0, 'subset_size']} "
        "annotation-blind density-weighting ablation",
        fontsize=16,
        fontweight="bold",
    )
    figure.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if not 0 < args.fraction <= 1:
        raise ValueError("--fraction must be > 0 and <= 1")
    if not 0 < args.rare_threshold < 1:
        raise ValueError("--rare_threshold must be > 0 and < 1")
    max_weights = parse_max_weights(args.max_weights)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    (
        cell_ids,
        features,
        mean_distance,
        sparsity_percentile,
        knn_k,
        baseline_positions,
        subset_size,
    ) = load_inputs(args)
    print(
        f"Loaded {args.atlas}: {features.shape[0]:,} cells x "
        f"{features.shape[1]:,} RBF features; subset={subset_size:,}; k={knn_k}",
        flush=True,
    )

    selections = {"baseline": baseline_positions}
    objectives = {"baseline": selection_objective(features, baseline_positions)}
    elapsed = {"baseline": 0.0}
    selected_sparsity = {"baseline": sparsity_percentile[baseline_positions]}

    # No metadata or labels have been loaded at this point.
    for max_weight in max_weights:
        method = method_name(max_weight)
        ranking, gains, seconds, weights = run_weighted_selection(
            features,
            sparsity_percentile,
            subset_size,
            max_weight,
            args.random_seed,
        )
        selections[method] = ranking
        objectives[method] = selection_objective(features, ranking)
        elapsed[method] = seconds
        selected_sparsity[method] = sparsity_percentile[ranking]
        pd.Series(cell_ids[ranking]).to_csv(
            output_dir / f"selected_cells_{method}_{args.fraction_label}.txt",
            index=False,
            header=False,
        )
        pd.DataFrame(
            {
                "selected_rank": np.arange(1, len(ranking) + 1),
                "cell_id": cell_ids[ranking],
                "selection_gain": gains,
                "sparsity_percentile": sparsity_percentile[ranking],
                "density_weight": weights[ranking],
                "mean_knn_distance": mean_distance[ranking],
            }
        ).to_csv(
            output_dir / f"selected_cell_metrics_{method}_{args.fraction_label}.tsv",
            sep="\t",
            index=False,
        )

    print("All annotation-blind selections complete; loading labels for evaluation", flush=True)
    labels = load_labels_after_selection(Path(args.metadata_tsv), cell_ids)
    summary, details = evaluate(
        args.atlas,
        args.fraction_label,
        labels,
        selections,
        sparsity_percentile,
        objectives,
        elapsed,
        args.rare_threshold,
    )
    summary.to_csv(output_dir / "density_ablation_summary.tsv", sep="\t", index=False)
    details.to_csv(
        output_dir / "density_ablation_cell_type_details.tsv", sep="\t", index=False
    )
    plot_results(summary, selected_sparsity, output_dir / "density_ablation_overview.png")
    (output_dir / "density_ablation.done").write_text("DONE\n", encoding="utf-8")
    print(f"Wrote density ablation outputs: {output_dir}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
