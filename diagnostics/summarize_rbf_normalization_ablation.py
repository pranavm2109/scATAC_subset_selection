#!/usr/bin/env python3
"""Compare baseline, L1, and capped singleton-normalized 2% selections."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.neighbors import NearestNeighbors


METHODS = ["baseline", "l1", "singleton_capped"]
METHOD_LABELS = {
    "baseline": "Baseline\nunnormalized",
    "l1": "L1 row\nnormalized",
    "singleton_capped": "Capped singleton-gain\nnormalized",
}
COLORS = {
    "baseline": "#607487",
    "l1": "#E07A3F",
    "singleton_capped": "#2A9D8F",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--feature_matrix", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--cell_metadata", required=True)
    parser.add_argument("--umap_tsv", required=True)
    parser.add_argument("--baseline_selected", required=True)
    parser.add_argument("--l1_selected", required=True)
    parser.add_argument("--singleton_selected", required=True)
    parser.add_argument("--l1_run_summary", required=True)
    parser.add_argument("--singleton_run_summary", required=True)
    parser.add_argument("--output_method_summary", required=True)
    parser.add_argument("--output_cell_type_details", required=True)
    parser.add_argument("--output_retention_changes", required=True)
    parser.add_argument("--output_overview_png", required=True)
    parser.add_argument("--output_composition_png", required=True)
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
    table = pd.read_csv(path, sep="\t", dtype={"cell_name": str})
    missing = set(required) - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    id_columns = [column for column in ("cell_name", "cell_id") if column in table.columns]
    if id_columns:
        table = table.set_index(id_columns[0])
    elif table.index.name != "cell_name" and len(table) == len(cell_ids):
        table.index = cell_ids
    missing_ids = pd.Index(cell_ids).difference(table.index.astype(str))
    if len(missing_ids):
        raise ValueError(f"{path} lacks {len(missing_ids):,} expected cells")
    table.index = table.index.astype(str)
    return table.loc[cell_ids]


def safe_correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    result = pearsonr(x, y) if method == "pearson" else spearmanr(x, y)
    return float(result.statistic)


def nearest_umap_summary(coords: np.ndarray, selected: np.ndarray, chunk: int) -> dict[str, float]:
    model = NearestNeighbors(n_neighbors=1, algorithm="kd_tree", n_jobs=1)
    model.fit(coords[selected])
    distances = np.empty(len(coords), dtype=np.float64)
    for start in range(0, len(coords), chunk):
        stop = min(start + chunk, len(coords))
        distances[start:stop] = model.kneighbors(
            coords[start:stop], return_distance=True
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
    metadata = align_table(
        Path(args.cell_metadata), cell_ids, [args.cell_type_column]
    )
    labels = metadata[args.cell_type_column].fillna("NA").astype(str)
    umap = align_table(Path(args.umap_tsv), cell_ids, ["UMAP_1", "UMAP_2"])
    coords = umap[["UMAP_1", "UMAP_2"]].apply(pd.to_numeric).to_numpy(float)
    if not np.isfinite(coords).all():
        raise ValueError("UMAP contains non-finite coordinates")

    features = np.load(args.feature_matrix, mmap_mode="r")
    if features.shape[0] != len(cell_ids):
        raise ValueError("RBF rows and cell IDs differ")

    selected_paths = {
        "baseline": Path(args.baseline_selected),
        "l1": Path(args.l1_selected),
        "singleton_capped": Path(args.singleton_selected),
    }
    selections: dict[str, np.ndarray] = {}
    selected_masks: dict[str, np.ndarray] = {}
    index = pd.Index(cell_ids)
    for method, path in selected_paths.items():
        selected_ids = read_ids(path)
        selected_indices = index.get_indexer(selected_ids)
        if np.any(selected_indices < 0):
            raise ValueError(f"{method} selection contains unknown cells")
        selections[method] = selected_indices
        mask = np.zeros(len(cell_ids), dtype=bool)
        mask[selected_indices] = True
        selected_masks[method] = mask

    selected_sizes = {len(value) for value in selections.values()}
    if len(selected_sizes) != 1:
        raise ValueError("All methods must select the same number of cells")

    full_counts = labels.value_counts().sort_index()
    full_prop = full_counts / len(labels)
    rare = full_prop < args.rare_threshold
    cell_rows = []
    method_rows = []
    baseline_retained = None
    baseline_objective = None

    for method in METHODS:
        selected_labels = labels.iloc[selections[method]]
        counts = selected_labels.value_counts().reindex(full_counts.index, fill_value=0)
        prop = counts / len(selected_labels)
        retained = counts > 0
        if method == "baseline":
            baseline_retained = retained.copy()

        coverage = np.zeros(features.shape[1], dtype=np.float64)
        for start in range(0, len(selections[method]), args.chunk_size):
            chosen = selections[method][start : start + args.chunk_size]
            coverage += np.asarray(features[chosen], dtype=np.float64).sum(axis=0)
        objective = float(np.sqrt(coverage).sum())
        if method == "baseline":
            baseline_objective = objective

        pearson = safe_correlation(full_prop.to_numpy(), prop.to_numpy(), "pearson")
        spearman = safe_correlation(full_prop.to_numpy(), prop.to_numpy(), "spearman")
        umap_stats = nearest_umap_summary(coords, selections[method], args.chunk_size)
        method_rows.append(
            {
                "atlas": args.atlas,
                "method": method,
                "n_selected": len(selections[method]),
                "n_cell_types_full": len(full_counts),
                "n_cell_types_retained": int(retained.sum()),
                "n_cell_types_missed": int((~retained).sum()),
                "n_rare_types_full": int(rare.sum()),
                "n_rare_types_retained": int((retained & rare).sum()),
                "n_types_gained_vs_baseline": int((retained & ~baseline_retained).sum()),
                "n_types_lost_vs_baseline": int((~retained & baseline_retained).sum()),
                "composition_pearson_r": pearson,
                "composition_spearman_rho": spearman,
                "composition_total_variation": float(0.5 * np.abs(prop - full_prop).sum()),
                "original_rbf_objective": objective,
                "original_rbf_objective_ratio_to_baseline": objective / baseline_objective,
                **umap_stats,
            }
        )
        for cell_type in full_counts.index:
            cell_rows.append(
                {
                    "atlas": args.atlas,
                    "method": method,
                    "cell_type": cell_type,
                    "full_count": int(full_counts[cell_type]),
                    "full_proportion": float(full_prop[cell_type]),
                    "is_rare": bool(rare[cell_type]),
                    "subset_count": int(counts[cell_type]),
                    "subset_proportion": float(prop[cell_type]),
                    "enrichment_ratio": float(prop[cell_type] / full_prop[cell_type]),
                    "retained": bool(retained[cell_type]),
                }
            )

    method_summary = pd.DataFrame(method_rows)
    details = pd.DataFrame(cell_rows)
    baseline_status = details.query("method == 'baseline'").set_index("cell_type")["retained"]
    changes = details.query("method != 'baseline'").copy()
    changes["baseline_retained"] = changes["cell_type"].map(baseline_status)
    changes["retention_change"] = np.select(
        [changes["retained"] & ~changes["baseline_retained"],
         ~changes["retained"] & changes["baseline_retained"]],
        ["gained", "lost"],
        default="unchanged",
    )

    outputs = [
        Path(args.output_method_summary), Path(args.output_cell_type_details),
        Path(args.output_retention_changes), Path(args.output_overview_png),
        Path(args.output_composition_png),
    ]
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
    method_summary.to_csv(outputs[0], sep="\t", index=False)
    details.to_csv(outputs[1], sep="\t", index=False)
    changes.to_csv(outputs[2], sep="\t", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), dpi=args.dpi)
    x = np.arange(len(METHODS))
    labels_short = [METHOD_LABELS[m] for m in METHODS]
    axes[0, 0].bar(x, method_summary["n_cell_types_retained"], color=[COLORS[m] for m in METHODS])
    axes[0, 0].axhline(len(full_counts), color="black", ls="--", lw=1)
    axes[0, 0].set_ylabel("Cell types retained")
    axes[0, 0].set_title("Cell-type retention at 2%")
    for i, value in enumerate(method_summary["n_cell_types_retained"]):
        axes[0, 0].text(i, value, f"{value}/{len(full_counts)}", ha="center", va="bottom")

    axes[0, 1].bar(x, method_summary["n_rare_types_retained"], color=[COLORS[m] for m in METHODS])
    axes[0, 1].axhline(int(rare.sum()), color="black", ls="--", lw=1)
    axes[0, 1].set_ylabel("Rare cell types retained")
    axes[0, 1].set_title(f"Rare retention (<{100*args.rare_threshold:g}% of atlas)")
    for i, value in enumerate(method_summary["n_rare_types_retained"]):
        axes[0, 1].text(i, value, f"{value}/{int(rare.sum())}", ha="center", va="bottom")

    axes[1, 0].bar(x, method_summary["original_rbf_objective_ratio_to_baseline"], color=[COLORS[m] for m in METHODS])
    axes[1, 0].axhline(1, color="black", ls="--", lw=1)
    axes[1, 0].set_ylabel("Objective / baseline objective")
    axes[1, 0].set_title("Coverage under original unnormalized RBF objective")

    axes[1, 1].bar(x, method_summary["umap_nearest_selected_q90"], color=[COLORS[m] for m in METHODS])
    axes[1, 1].set_ylabel("90th percentile nearest-selected distance")
    axes[1, 1].set_title("PeakVI UMAP coverage (lower is better)")
    for ax in axes.flat:
        ax.set_xticks(x, labels_short)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{args.atlas}: RBF normalization ablation at 2%", fontsize=18, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outputs[3], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4), dpi=args.dpi, sharex=True, sharey=True)
    max_value = float(details[["full_proportion", "subset_proportion"]].max().max()) * 1.05
    for ax, method in zip(axes, METHODS):
        current = details[details["method"] == method]
        colors = np.where(current["is_rare"], "#E76F51", "#607487")
        ax.scatter(current["full_proportion"], current["subset_proportion"], c=colors, s=20, alpha=0.75)
        ax.plot([0, max_value], [0, max_value], "k--", lw=1)
        row = method_summary[method_summary["method"] == method].iloc[0]
        ax.set_title(
            f"{METHOD_LABELS[method].replace(chr(10), ' ')}\n"
            f"Pearson={row['composition_pearson_r']:.3f}; Spearman={row['composition_spearman_rho']:.3f}\n"
            f"missed={int(row['n_cell_types_missed'])}"
        )
        ax.set_xlabel("Full-atlas proportion")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("2% subset proportion")
    fig.suptitle(f"{args.atlas}: composition fidelity", fontsize=17, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outputs[4], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(method_summary.to_string(index=False), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
