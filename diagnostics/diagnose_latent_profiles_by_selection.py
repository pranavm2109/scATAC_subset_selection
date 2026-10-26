#!/usr/bin/env python3
"""Compare PeakVI latent profiles of retained and missed cell types.

For each annotated cell type, this diagnostic computes its mean value in every
PeakVI latent dimension. Values are standardized against the atlas-wide mean
and standard deviation of the corresponding latent dimension, making rows
comparable without forcing every cell type to have the same within-row range.

The figure orders cell types from highest to lowest mean enrichment across the
five selected subsets and adds a retained/missed tile for every subset size.
Cell-type labels are used only for retrospective evaluation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SUBSET_ORDER = ["02pct", "04pct", "06pct", "08pct", "10pct"]
SUBSET_LABELS = {
    "02pct": "2%",
    "04pct": "4%",
    "06pct": "6%",
    "08pct": "8%",
    "10pct": "10%",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot cell-type mean PeakVI latent profiles ordered by subset "
            "enrichment, with retained/missed status at every subset size."
        )
    )
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--latent_npy", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--cell_metadata", required=True)
    parser.add_argument("--proportions_tsv", required=True)
    parser.add_argument("--output_png", required=True)
    parser.add_argument("--output_summary_tsv", required=True)
    parser.add_argument("--cell_type_column", default="cell_type")
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=100000,
        help="PeakVI latent rows processed per chunk [default: 100000].",
    )
    parser.add_argument(
        "--z_clip",
        type=float,
        default=3.0,
        help="Symmetric color limit for atlas-standardized means [default: 3].",
    )
    parser.add_argument("--dpi", type=int, default=320)
    return parser.parse_args()


def read_cell_names(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if len(values) > 1 and values.iloc[0] in {"cell_name", "cell_id"}:
        values = values.iloc[1:]
    cell_names = values.astype(str).to_numpy()
    if len(cell_names) == 0:
        raise ValueError(f"Cell-name file is empty: {path}")
    if len(np.unique(cell_names)) != len(cell_names):
        raise ValueError(f"Cell names are not unique: {path}")
    return cell_names


def read_aligned_cell_types(
    path: Path,
    cell_names: np.ndarray,
    cell_type_column: str,
) -> pd.Series:
    metadata = pd.read_csv(path, sep="\t", dtype=str)
    if cell_type_column not in metadata.columns:
        raise ValueError(
            f"Cell metadata lacks '{cell_type_column}': {path}. "
            f"Available columns: {list(metadata.columns)}"
        )
    id_candidates = [
        column
        for column in ("cell_id", "cell_name", metadata.columns[0])
        if column in metadata.columns
    ]
    cell_id_column = id_candidates[0]
    if metadata[cell_id_column].duplicated().any():
        raise ValueError(
            f"Cell identifiers are not unique in {path}: {cell_id_column}"
        )

    annotations = metadata.set_index(cell_id_column)[cell_type_column]
    aligned = annotations.reindex(cell_names)
    if aligned.isna().any():
        missing = cell_names[aligned.isna().to_numpy()][:10]
        raise ValueError(
            f"Metadata is missing {aligned.isna().sum():,} PeakVI cells. "
            f"Examples: {missing.tolist()}"
        )
    return aligned.astype(str)


def summarize_latent_by_cell_type(
    latent: np.ndarray,
    cell_types: pd.Series,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    categories = np.array(sorted(cell_types.unique()), dtype=object)
    category_index = {value: index for index, value in enumerate(categories)}
    codes = cell_types.map(category_index).to_numpy(dtype=np.int32)
    n_types = len(categories)
    n_dimensions = latent.shape[1]

    group_sums = np.zeros((n_types, n_dimensions), dtype=np.float64)
    group_counts = np.bincount(codes, minlength=n_types).astype(np.int64)
    atlas_sum = np.zeros(n_dimensions, dtype=np.float64)
    atlas_sum_squares = np.zeros(n_dimensions, dtype=np.float64)

    for start in range(0, latent.shape[0], chunk_size):
        stop = min(start + chunk_size, latent.shape[0])
        block = np.asarray(latent[start:stop], dtype=np.float64)
        if not np.isfinite(block).all():
            raise ValueError(
                f"PeakVI latent matrix contains non-finite values in rows "
                f"{start:,}-{stop:,}"
            )
        block_codes = codes[start:stop]
        atlas_sum += block.sum(axis=0)
        atlas_sum_squares += np.square(block).sum(axis=0)
        for code in np.unique(block_codes):
            group_sums[code] += block[block_codes == code].sum(axis=0)
        print(
            f"Summarized latent rows {start:,}-{stop:,} of {latent.shape[0]:,}",
            flush=True,
        )

    if np.any(group_counts == 0):
        raise ValueError("At least one annotated cell type has zero aligned cells")
    group_means = group_sums / group_counts[:, None]
    atlas_mean = atlas_sum / latent.shape[0]
    atlas_variance = (
        atlas_sum_squares / latent.shape[0] - np.square(atlas_mean)
    )
    atlas_variance = np.maximum(atlas_variance, 0.0)
    atlas_sd = np.sqrt(atlas_variance)
    if np.any(atlas_sd <= np.finfo(np.float64).eps):
        constant = np.flatnonzero(atlas_sd <= np.finfo(np.float64).eps) + 1
        raise ValueError(
            "PeakVI contains constant latent dimensions: "
            + ", ".join(map(str, constant.tolist()))
        )
    standardized_means = (group_means - atlas_mean) / atlas_sd
    return categories, group_counts, group_means, standardized_means, atlas_sd


def read_selection_metrics(
    path: Path,
    categories: np.ndarray,
) -> pd.DataFrame:
    proportions = pd.read_csv(path, sep="\t")
    required = {"cell_type", "atlas_count", "atlas_proportion"}
    for subset in SUBSET_ORDER:
        required.add(f"subset_{subset}_count")
        required.add(f"ratio_{subset}_vs_atlas")
    missing_columns = required - set(proportions.columns)
    if missing_columns:
        raise ValueError(
            f"Proportions TSV is missing columns: {sorted(missing_columns)}"
        )
    if proportions["cell_type"].duplicated().any():
        raise ValueError("Proportions TSV contains duplicate cell types")

    proportions["cell_type"] = proportions["cell_type"].astype(str)
    proportions = proportions.set_index("cell_type").reindex(categories)
    if proportions["atlas_count"].isna().any():
        missing = proportions.index[proportions["atlas_count"].isna()].tolist()
        raise ValueError(
            "Proportions TSV is missing annotated cell types: "
            + "; ".join(missing)
        )

    ratio_columns = [f"ratio_{subset}_vs_atlas" for subset in SUBSET_ORDER]
    ratios = proportions[ratio_columns].apply(pd.to_numeric, errors="coerce")
    if ratios.isna().all(axis=1).any():
        bad = ratios.index[ratios.isna().all(axis=1)].tolist()
        raise ValueError(
            "No finite enrichment ratios were found for: " + "; ".join(bad)
        )
    proportions["mean_enrichment_ratio"] = ratios.mean(axis=1, skipna=True)

    for subset in SUBSET_ORDER:
        count_column = f"subset_{subset}_count"
        counts = pd.to_numeric(proportions[count_column], errors="raise")
        proportions[f"retained_{subset}"] = counts > 0
    retained_columns = [f"retained_{subset}" for subset in SUBSET_ORDER]
    proportions["n_subsets_retained"] = proportions[retained_columns].sum(axis=1)
    proportions["retention_class"] = np.select(
        [
            proportions["n_subsets_retained"] == 0,
            proportions["n_subsets_retained"] == len(SUBSET_ORDER),
        ],
        ["Never retained", "Retained in all subsets"],
        default="Retained in some subsets",
    )
    return proportions.reset_index()


def make_summary_table(
    atlas: str,
    categories: np.ndarray,
    group_counts: np.ndarray,
    group_means: np.ndarray,
    standardized_means: np.ndarray,
    metrics: pd.DataFrame,
) -> pd.DataFrame:
    metric_lookup = metrics.set_index("cell_type")
    summary = pd.DataFrame(
        {
            "atlas": atlas,
            "cell_type": categories,
            "n_cells": group_counts,
        }
    )
    metric_columns = [
        "atlas_proportion",
        "mean_enrichment_ratio",
        *[f"retained_{subset}" for subset in SUBSET_ORDER],
        "n_subsets_retained",
        "retention_class",
    ]
    for column in metric_columns:
        summary[column] = metric_lookup.loc[categories, column].to_numpy()
    for index in range(group_means.shape[1]):
        summary[f"PeakVI_{index + 1}_mean"] = group_means[:, index]
        summary[f"PeakVI_{index + 1}_atlas_z"] = standardized_means[:, index]
    return summary.sort_values(
        ["mean_enrichment_ratio", "n_cells", "cell_type"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def plot_diagnostic(
    atlas: str,
    summary: pd.DataFrame,
    output_path: Path,
    z_clip: float,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    z_columns = [
        column
        for column in summary.columns
        if column.startswith("PeakVI_") and column.endswith("_atlas_z")
    ]
    retained_columns = [f"retained_{subset}" for subset in SUBSET_ORDER]
    latent_values = summary[z_columns].to_numpy(dtype=float)
    retained = summary[retained_columns].to_numpy(dtype=bool)

    n_types = len(summary)
    figure_height = max(8.5, 0.24 * n_types + 2.8)
    figure = plt.figure(figsize=(18, figure_height), dpi=dpi)
    grid = figure.add_gridspec(
        1,
        2,
        width_ratios=[12, 1.35],
        left=0.24,
        right=0.96,
        bottom=0.08,
        top=0.91,
        wspace=0.035,
    )
    latent_axis = figure.add_subplot(grid[0, 0])
    retention_axis = figure.add_subplot(grid[0, 1])

    latent_image = latent_axis.imshow(
        latent_values,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-z_clip,
        vmax=z_clip,
        rasterized=True,
    )
    cell_labels = [
        f"{row.cell_type}  ({row.mean_enrichment_ratio:.2f}x)"
        for row in summary.itertuples(index=False)
    ]
    latent_axis.set_yticks(np.arange(n_types))
    latent_axis.set_yticklabels(cell_labels, fontsize=7.5)
    latent_axis.set_xticks(np.arange(len(z_columns)))
    latent_axis.set_xticklabels(
        [f"{index + 1}" for index in range(len(z_columns))],
        rotation=90,
        fontsize=7,
    )
    latent_axis.set_xlabel("PeakVI latent dimension", fontsize=11)
    latent_axis.set_ylabel(
        "Cell type (mean enrichment across subsets)",
        fontsize=11,
    )
    latent_axis.tick_params(length=0)
    for spine in latent_axis.spines.values():
        spine.set_visible(False)

    retention_axis.imshow(
        retained.astype(int),
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(["#D95F4C", "#167D7F"]),
        vmin=0,
        vmax=1,
    )
    retention_axis.set_xticks(np.arange(len(SUBSET_ORDER)))
    retention_axis.set_xticklabels(
        [SUBSET_LABELS[subset] for subset in SUBSET_ORDER],
        rotation=90,
        fontsize=8,
    )
    retention_axis.set_yticks([])
    retention_axis.set_xlabel("Subset size", fontsize=10)
    retention_axis.tick_params(length=0)
    for spine in retention_axis.spines.values():
        spine.set_visible(False)

    colorbar = figure.colorbar(
        latent_image,
        ax=[latent_axis, retention_axis],
        orientation="horizontal",
        fraction=0.025,
        pad=0.045,
        aspect=45,
    )
    colorbar.set_label(
        "Cell-type mean relative to atlas mean (atlas SD units)",
        fontsize=10,
    )
    figure.legend(
        handles=[
            Patch(facecolor="#167D7F", label="Retained (>=1 cell)"),
            Patch(facecolor="#D95F4C", label="Missed (0 cells)"),
        ],
        loc="upper right",
        bbox_to_anchor=(0.965, 0.975),
        frameon=False,
        fontsize=9,
    )
    figure.suptitle(
        f"{atlas}: PeakVI latent profiles versus subset retention",
        fontsize=20,
        fontweight="bold",
        y=0.975,
    )
    figure.text(
        0.5,
        0.94,
        (
            "Cell types are ordered from highest to lowest mean enrichment "
            "across the 2-10% subsets; labels were used only for evaluation."
        ),
        ha="center",
        va="center",
        fontsize=11,
        color="#333333",
    )
    figure.text(
        0.01,
        0.012,
        (
            "Values are cell-type mean latent coordinates standardized using "
            "the atlas-wide mean and SD of each dimension; colors are clipped "
            f"at +/-{z_clip:g} SD for display."
        ),
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#4D5960",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if args.chunk_size < 1:
        raise ValueError("--chunk_size must be positive")
    if not np.isfinite(args.z_clip) or args.z_clip <= 0:
        raise ValueError("--z_clip must be a positive finite number")

    latent_path = Path(args.latent_npy)
    cell_names_path = Path(args.cell_names)
    metadata_path = Path(args.cell_metadata)
    proportions_path = Path(args.proportions_tsv)
    output_png = Path(args.output_png)
    output_summary = Path(args.output_summary_tsv)

    cell_names = read_cell_names(cell_names_path)
    latent = np.load(latent_path, mmap_mode="r")
    if latent.ndim != 2:
        raise ValueError("PeakVI latent array must be two-dimensional")
    if latent.shape[0] != len(cell_names):
        raise ValueError(
            f"Latent rows ({latent.shape[0]:,}) do not match cell names "
            f"({len(cell_names):,})"
        )
    if latent.shape[1] != 60:
        print(
            f"WARNING: expected 60 PeakVI dimensions, found {latent.shape[1]}",
            flush=True,
        )

    cell_types = read_aligned_cell_types(
        metadata_path,
        cell_names,
        args.cell_type_column,
    )
    categories, counts, means, standardized, atlas_sd = (
        summarize_latent_by_cell_type(latent, cell_types, args.chunk_size)
    )
    metrics = read_selection_metrics(proportions_path, categories)
    summary = make_summary_table(
        args.atlas,
        categories,
        counts,
        means,
        standardized,
        metrics,
    )

    output_summary.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_summary, sep="\t", index=False)
    plot_diagnostic(
        args.atlas,
        summary,
        output_png,
        args.z_clip,
        args.dpi,
    )

    never_retained = summary.loc[
        summary["retention_class"] == "Never retained", "cell_type"
    ].tolist()
    print(f"Wrote latent-profile heatmap: {output_png}", flush=True)
    print(f"Wrote latent-profile summary: {output_summary}", flush=True)
    print(
        f"Cell types retained in all subsets: "
        f"{(summary['n_subsets_retained'] == len(SUBSET_ORDER)).sum()}/"
        f"{len(summary)}",
        flush=True,
    )
    if never_retained:
        print("Never-retained cell types: " + "; ".join(never_retained), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
