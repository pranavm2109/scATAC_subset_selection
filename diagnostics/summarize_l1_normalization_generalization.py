#!/usr/bin/env python3
"""Evaluate baseline versus L1-normalized RBF selection across 2-10%."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.neighbors import NearestNeighbors


FRACTIONS = {
    "02pct": 0.02,
    "04pct": 0.04,
    "06pct": 0.06,
    "08pct": 0.08,
    "10pct": 0.10,
}
METHODS = ["baseline", "l1"]
COLORS = {"baseline": "#607487", "l1": "#E07A3F"}
LABELS = {"baseline": "Baseline unnormalized", "l1": "L1 normalized"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--cell_metadata", required=True)
    parser.add_argument("--umap_tsv", required=True)
    parser.add_argument("--baseline_selected", nargs=5, required=True)
    parser.add_argument("--l1_selected", nargs=5, required=True)
    parser.add_argument("--output_method_summary", required=True)
    parser.add_argument("--output_cell_type_details", required=True)
    parser.add_argument("--output_retention_changes", required=True)
    parser.add_argument("--output_retention_png", required=True)
    parser.add_argument("--output_composition_png", required=True)
    parser.add_argument("--output_umap_coverage_png", required=True)
    parser.add_argument("--cell_type_column", default="cell_type")
    parser.add_argument("--rare_threshold", type=float, default=0.01)
    parser.add_argument("--chunk_size", type=int, default=50000)
    parser.add_argument("--dpi", type=int, default=320)
    return parser.parse_args()


def read_ids(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any() or values.duplicated().any():
        raise ValueError(f"Invalid cell IDs: {path}")
    return values.to_numpy(dtype=str)


def align_table(path: Path, cell_ids: np.ndarray, required: list[str]) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    missing = set(required) - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    id_columns = [column for column in ("cell_name", "cell_id") if column in table.columns]
    if not id_columns:
        raise ValueError(f"{path} needs a cell_name or cell_id column")
    table = table.set_index(id_columns[0])
    table.index = table.index.astype(str)
    if table.index.has_duplicates:
        raise ValueError(f"Duplicate cell identifiers in {path}")
    missing_ids = pd.Index(cell_ids).difference(table.index)
    if len(missing_ids):
        raise ValueError(f"{path} lacks {len(missing_ids):,} expected cells")
    return table.loc[cell_ids]


def safe_correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    result = pearsonr(x, y) if method == "pearson" else spearmanr(x, y)
    return float(result.statistic)


def nearest_selected_summary(
    coordinates: np.ndarray,
    selected_indices: np.ndarray,
    chunk_size: int,
) -> dict[str, float]:
    model = NearestNeighbors(n_neighbors=1, algorithm="kd_tree", n_jobs=1)
    model.fit(coordinates[selected_indices])
    distances = np.empty(len(coordinates), dtype=np.float64)
    for start in range(0, len(coordinates), chunk_size):
        stop = min(start + chunk_size, len(coordinates))
        distances[start:stop] = model.kneighbors(
            coordinates[start:stop], return_distance=True
        )[0][:, 0]
    return {
        "umap_nearest_selected_mean": float(np.mean(distances)),
        "umap_nearest_selected_median": float(np.median(distances)),
        "umap_nearest_selected_q90": float(np.quantile(distances, 0.90)),
        "umap_nearest_selected_q99": float(np.quantile(distances, 0.99)),
    }


def main() -> None:
    args = parse_args()
    if not 0 < args.rare_threshold < 1:
        raise ValueError("--rare_threshold must be between 0 and 1")

    cell_ids = read_ids(Path(args.cell_names))
    cell_index = pd.Index(cell_ids)
    metadata = align_table(
        Path(args.cell_metadata), cell_ids, [args.cell_type_column]
    )
    labels = metadata[args.cell_type_column].fillna("NA").astype(str)
    umap = align_table(Path(args.umap_tsv), cell_ids, ["UMAP_1", "UMAP_2"])
    coordinates = umap[["UMAP_1", "UMAP_2"]].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=float)
    if not np.isfinite(coordinates).all():
        raise ValueError("UMAP coordinates contain non-finite values")

    paths = {
        "baseline": dict(zip(FRACTIONS, map(Path, args.baseline_selected))),
        "l1": dict(zip(FRACTIONS, map(Path, args.l1_selected))),
    }
    full_counts = labels.value_counts().sort_index()
    full_proportions = full_counts / len(labels)
    rare = full_proportions < args.rare_threshold
    method_rows: list[dict] = []
    cell_rows: list[dict] = []
    retention_by_key: dict[tuple[str, str], pd.Series] = {}

    for fraction_label, fraction in FRACTIONS.items():
        expected_size = max(int(round(len(cell_ids) * fraction)), 1)
        for method in METHODS:
            selected_ids = read_ids(paths[method][fraction_label])
            if len(selected_ids) != expected_size:
                raise ValueError(
                    f"{method} {fraction_label} has {len(selected_ids):,} cells; "
                    f"expected {expected_size:,}"
                )
            selected_indices = cell_index.get_indexer(selected_ids)
            if np.any(selected_indices < 0):
                raise ValueError(f"{method} {fraction_label} contains unknown cells")

            subset_counts = (
                labels.iloc[selected_indices]
                .value_counts()
                .reindex(full_counts.index, fill_value=0)
            )
            subset_proportions = subset_counts / expected_size
            retained = subset_counts > 0
            retention_by_key[(method, fraction_label)] = retained
            coverage = nearest_selected_summary(
                coordinates, selected_indices, args.chunk_size
            )
            method_rows.append(
                {
                    "atlas": args.atlas,
                    "fraction_label": fraction_label,
                    "fraction": fraction,
                    "method": method,
                    "n_selected": expected_size,
                    "n_cell_types_full": len(full_counts),
                    "n_cell_types_retained": int(retained.sum()),
                    "n_cell_types_missed": int((~retained).sum()),
                    "all_cell_types_retained": bool(retained.all()),
                    "n_rare_types_full": int(rare.sum()),
                    "n_rare_types_retained": int((retained & rare).sum()),
                    "all_rare_types_retained": bool((retained | ~rare).all()),
                    "composition_pearson_r": safe_correlation(
                        full_proportions.to_numpy(),
                        subset_proportions.to_numpy(),
                        "pearson",
                    ),
                    "composition_spearman_rho": safe_correlation(
                        full_proportions.to_numpy(),
                        subset_proportions.to_numpy(),
                        "spearman",
                    ),
                    "composition_total_variation": float(
                        0.5 * np.abs(subset_proportions - full_proportions).sum()
                    ),
                    **coverage,
                }
            )
            for cell_type in full_counts.index:
                cell_rows.append(
                    {
                        "atlas": args.atlas,
                        "fraction_label": fraction_label,
                        "fraction": fraction,
                        "method": method,
                        "cell_type": cell_type,
                        "full_count": int(full_counts[cell_type]),
                        "full_proportion": float(full_proportions[cell_type]),
                        "is_rare": bool(rare[cell_type]),
                        "subset_count": int(subset_counts[cell_type]),
                        "subset_proportion": float(subset_proportions[cell_type]),
                        "enrichment_ratio": float(
                            subset_proportions[cell_type] / full_proportions[cell_type]
                        ),
                        "retained": bool(retained[cell_type]),
                    }
                )

    summary = pd.DataFrame(method_rows)
    details = pd.DataFrame(cell_rows)
    change_rows = []
    for fraction_label in FRACTIONS:
        baseline = retention_by_key[("baseline", fraction_label)]
        l1 = retention_by_key[("l1", fraction_label)]
        for cell_type in full_counts.index:
            status = "unchanged"
            if l1[cell_type] and not baseline[cell_type]:
                status = "gained"
            elif baseline[cell_type] and not l1[cell_type]:
                status = "lost"
            change_rows.append(
                {
                    "atlas": args.atlas,
                    "fraction_label": fraction_label,
                    "cell_type": cell_type,
                    "full_count": int(full_counts[cell_type]),
                    "is_rare": bool(rare[cell_type]),
                    "baseline_retained": bool(baseline[cell_type]),
                    "l1_retained": bool(l1[cell_type]),
                    "retention_change": status,
                }
            )
    changes = pd.DataFrame(change_rows)
    change_counts = (
        changes.groupby(["fraction_label", "retention_change"])
        .size()
        .unstack(fill_value=0)
    )
    for status in ("gained", "lost"):
        summary[f"n_types_{status}_vs_baseline"] = summary.apply(
            lambda row: int(change_counts.get(status, pd.Series()).get(
                row["fraction_label"], 0
            )) if row["method"] == "l1" else 0,
            axis=1,
        )

    outputs = [
        Path(args.output_method_summary),
        Path(args.output_cell_type_details),
        Path(args.output_retention_changes),
        Path(args.output_retention_png),
        Path(args.output_composition_png),
        Path(args.output_umap_coverage_png),
    ]
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outputs[0], sep="\t", index=False)
    details.to_csv(outputs[1], sep="\t", index=False)
    changes.to_csv(outputs[2], sep="\t", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_values = np.array([100 * value for value in FRACTIONS.values()])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=args.dpi)
    for method in METHODS:
        current = summary[summary["method"] == method]
        axes[0].plot(
            x_values, current["n_cell_types_retained"], marker="o",
            lw=2, color=COLORS[method], label=LABELS[method],
        )
        axes[1].plot(
            x_values, current["n_rare_types_retained"], marker="o",
            lw=2, color=COLORS[method], label=LABELS[method],
        )
    axes[0].axhline(len(full_counts), color="black", ls="--", lw=1)
    axes[1].axhline(int(rare.sum()), color="black", ls="--", lw=1)
    axes[0].set_title("All annotated cell types")
    axes[1].set_title(f"Rare cell types (<{100 * args.rare_threshold:g}%)")
    axes[0].set_ylabel("Number retained")
    for ax in axes:
        ax.set_xlabel("Subset size (%)")
        ax.set_xticks(x_values)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].legend(frameon=False)
    fig.suptitle(f"{args.atlas}: cell-type retention across subset sizes", fontsize=17, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outputs[3], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    metrics = [
        ("composition_pearson_r", "Pearson r", True),
        ("composition_spearman_rho", "Spearman rho", True),
        ("composition_total_variation", "Total variation (lower is better)", False),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), dpi=args.dpi)
    for ax, (column, title, _) in zip(axes, metrics):
        for method in METHODS:
            current = summary[summary["method"] == method]
            ax.plot(x_values, current[column], marker="o", lw=2, color=COLORS[method], label=LABELS[method])
        ax.set_title(title)
        ax.set_xlabel("Subset size (%)")
        ax.set_xticks(x_values)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False)
    fig.suptitle(f"{args.atlas}: composition fidelity", fontsize=17, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outputs[4], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=args.dpi)
    for method in METHODS:
        current = summary[summary["method"] == method]
        axes[0].plot(x_values, current["umap_nearest_selected_q90"], marker="o", lw=2, color=COLORS[method], label=LABELS[method])
        axes[1].plot(x_values, current["umap_nearest_selected_q99"], marker="o", lw=2, color=COLORS[method], label=LABELS[method])
    axes[0].set_title("90th percentile")
    axes[1].set_title("99th percentile")
    axes[0].set_ylabel("Nearest-selected UMAP distance")
    for ax in axes:
        ax.set_xlabel("Subset size (%)")
        ax.set_xticks(x_values)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].legend(frameon=False)
    fig.suptitle(f"{args.atlas}: PeakVI UMAP coverage (lower is better)", fontsize=17, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outputs[5], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(summary.to_string(index=False), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
