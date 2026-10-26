#!/usr/bin/env python3
"""Diagnose how PeakVI features interact with Apricot FeatureBasedSelection.

This is a read-only retrospective diagnostic. It uses the exact nonnegative
RBF landmark feature matrix consumed by Apricot, computes
annotation-blind per-cell utility and completed-subset marginal gains, and only
then joins cell-type annotations for interpretation.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SUBSETS = ("02pct", "04pct", "06pct", "08pct", "10pct")
SUBSET_LABELS = dict(zip(SUBSETS, ("2%", "4%", "6%", "8%", "10%")))
STATUSES = ("absent", "underpowered", "evaluable")
STATUS_COLORS = {
    "absent": "#C83E4D",
    "underpowered": "#F2A541",
    "evaluable": "#4C78A8",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Diagnose the current PeakVI-to-Apricot feature objective."
    )
    parser.add_argument("--latent_tsv", required=True)
    parser.add_argument("--metadata_tsv", required=True)
    parser.add_argument("--separability_tsv", required=True)
    parser.add_argument("--selected_cells", nargs="+", required=True)
    parser.add_argument("--selection_summaries", nargs="+", required=True)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--output_table", required=True)
    parser.add_argument("--output_latent_summary", required=True)
    parser.add_argument("--output_feature_profiles", required=True)
    parser.add_argument("--output_figures_dir", required=True)
    parser.add_argument("--chunk_size", type=int, default=100000)
    parser.add_argument("--high_purity_threshold", type=float, default=0.75)
    parser.add_argument("--max_labels", type=int, default=15)
    parser.add_argument("--simplex_tolerance", type=float, default=1e-5)
    return parser.parse_args()


def subset_from_path(path: Path) -> str:
    match = re.search(r"_(02pct|04pct|06pct|08pct|10pct)(?:\.txt|\.tsv)$", path.name)
    if match is None:
        raise ValueError(f"Cannot infer subset size from filename: {path}")
    return match.group(1)


def read_ids(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any():
        raise ValueError(f"Empty or missing cell IDs in {path}")
    values = values.to_numpy(dtype=str)
    if len(np.unique(values)) != len(values):
        raise ValueError(f"Duplicate cell IDs in {path}")
    return values


def validate_selection_summaries(paths: list[str], atlas: str) -> tuple[Path, int]:
    observed = {}
    for value in paths:
        path = Path(value)
        subset = subset_from_path(path)
        frame = pd.read_csv(path, sep="\t")
        if len(frame) != 1:
            raise ValueError(f"Selection summary must have one row: {path}")
        required = {
            "atlas",
            "fraction_label",
            "nonnegative_transform",
            "selection_function",
            "n_latent",
            "feature_construction_method",
            "n_landmark_features",
            "rbf_feature_matrix",
        }
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        row = frame.iloc[0]
        if str(row["atlas"]) != atlas or str(row["fraction_label"]) != subset:
            raise ValueError(f"Atlas or subset mismatch in {path}")
        observed[subset] = (
            str(row["nonnegative_transform"]),
            str(row["selection_function"]),
            int(row["n_latent"]),
            str(row["feature_construction_method"]),
            int(row["n_landmark_features"]),
            str(row["rbf_feature_matrix"]),
        )

    if set(observed) != set(SUBSETS):
        raise ValueError("Exactly one selection summary is required per subset")
    settings = set(observed.values())
    if len(settings) != 1:
        raise ValueError(f"Selection settings differ among subsets: {observed}")

    transform, function, _, feature_method, _, feature_matrix = next(iter(settings))
    if (
        transform != "none"
        or function != "sqrt"
        or feature_method != "peakvi_rbf_landmarks"
    ):
        raise ValueError(
            "This diagnostic implements the current none + RBF landmarks + "
            f"sqrt objective, but summaries report {transform} + "
            f"{feature_method} + {function}"
        )
    return Path(feature_matrix), next(iter(settings))[4]


def load_selection_features(
    path: Path, expected_cells: int, expected_features: int
) -> np.ndarray:
    print(f"Memory-mapping the exact RBF feature matrix used by Apricot: {path}", flush=True)
    features = np.load(path, mmap_mode="r")
    if features.ndim != 2 or features.shape[0] != expected_cells:
        raise ValueError(
            "RBF feature matrix shape does not match the latent cell order: "
            f"features={features.shape}, expected_cells={expected_cells:,}"
        )
    if features.shape[1] != expected_features:
        raise ValueError(
            f"RBF matrix has {features.shape[1]} features; expected {expected_features}"
        )
    if not np.issubdtype(features.dtype, np.number):
        raise ValueError(f"RBF feature matrix is not numeric: {features.dtype}")
    # RBF similarities are intrinsically nonnegative, so no min shift is used.
    if float(features.min()) < 0 or not np.isfinite(features).all():
        raise ValueError("RBF feature matrix contains negative or non-finite values")
    return features


def load_latent(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    print(f"Reading the exact latent TSV used by Apricot: {path}", flush=True)
    frame = pd.read_csv(path, sep="\t")
    if frame.shape[1] < 2:
        raise ValueError("Latent TSV needs a cell-name column and latent dimensions")
    cell_ids = frame.iloc[:, 0].astype(str).to_numpy()
    if len(np.unique(cell_ids)) != len(cell_ids):
        raise ValueError("Latent TSV cell IDs are not unique")
    dimension_names = frame.columns[1:].astype(str).tolist()
    raw = frame.iloc[:, 1:].to_numpy(dtype=np.float64, copy=True)
    del frame
    if not np.isfinite(raw).all():
        raise ValueError("PeakVI latent matrix contains non-finite values")
    return raw, cell_ids, dimension_names


def latent_distribution_summary(
    atlas: str,
    raw: np.ndarray,
    dimension_names: list[str],
    tolerance: float,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    dimension_min = raw.min(axis=0)
    dimension_max = raw.max(axis=0)
    row_sum = raw.sum(axis=1)
    fraction_negative = float(np.mean(raw < 0))
    fraction_zero = float(np.mean(raw == 0))
    fraction_rows_near_one = float(np.mean(np.abs(row_sum - 1) <= tolerance))
    nonnegative = bool(fraction_negative == 0)
    approximately_simplex = bool(nonnegative and fraction_rows_near_one >= 0.99)

    common = {
        "atlas": atlas,
        "n_cells": raw.shape[0],
        "n_dimensions": raw.shape[1],
        "nonnegative_transform_reproduced": "none",
        "selection_function_reproduced": "sqrt",
        "simplex_tolerance": tolerance,
    }
    rows = [
        {
            **common,
            "record_type": "global",
            "dimension": "ALL",
            "raw_min": float(raw.min()),
            "raw_max": float(raw.max()),
            "transformed_min": np.nan,
            "transformed_max": np.nan,
            "fraction_values_negative": fraction_negative,
            "fraction_values_zero": fraction_zero,
            "row_sum_min": float(row_sum.min()),
            "row_sum_median": float(np.median(row_sum)),
            "row_sum_mean": float(row_sum.mean()),
            "row_sum_max": float(row_sum.max()),
            "fraction_row_sums_near_one": fraction_rows_near_one,
            "all_values_nonnegative": nonnegative,
            "approximately_simplex_valued": approximately_simplex,
            "shift_min_needed_for_nonnegativity": np.nan,
        }
    ]
    for index, name in enumerate(dimension_names):
        rows.append(
            {
                **common,
                "record_type": "dimension",
                "dimension": name,
                "raw_min": float(dimension_min[index]),
                "raw_max": float(dimension_max[index]),
                "transformed_min": np.nan,
                "transformed_max": np.nan,
            }
        )
    return pd.DataFrame(rows), dimension_min, row_sum


def per_cell_intrinsic(
    raw: np.ndarray, transformed: np.ndarray, chunk_size: int
) -> dict[str, np.ndarray]:
    n_cells = len(raw)
    output = {
        "singleton_score": np.empty(n_cells),
        "transformed_l1": np.empty(n_cells),
        "transformed_l2": np.empty(n_cells),
        "transformed_max": np.empty(n_cells),
        "transformed_n_positive": np.empty(n_cells, dtype=np.int16),
        "transformed_entropy": np.empty(n_cells),
        "raw_max": np.empty(n_cells),
        "raw_l2": np.empty(n_cells),
    }
    print("Calculating annotation-blind singleton and latent metrics", flush=True)
    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        x = transformed[start:stop]
        z = raw[start:stop]
        l1 = x.sum(axis=1)
        x_log_x = np.zeros_like(x)
        np.log(x, out=x_log_x, where=x > 0)
        x_log_x *= x
        entropy = np.zeros(len(x))
        nonzero_sum = l1 > 0
        entropy[nonzero_sum] = (
            np.log(l1[nonzero_sum])
            - x_log_x.sum(axis=1)[nonzero_sum] / l1[nonzero_sum]
        )
        output["singleton_score"][start:stop] = np.sqrt(x).sum(axis=1)
        output["transformed_l1"][start:stop] = l1
        output["transformed_l2"][start:stop] = np.sqrt(np.square(x).sum(axis=1))
        output["transformed_max"][start:stop] = x.max(axis=1)
        output["transformed_n_positive"][start:stop] = (x > 0).sum(axis=1)
        output["transformed_entropy"][start:stop] = entropy
        output["raw_max"][start:stop] = z.max(axis=1)
        output["raw_l2"][start:stop] = np.sqrt(np.square(z).sum(axis=1))
    return output


def load_selected_positions(
    paths: list[str], cell_ids: np.ndarray
) -> dict[str, np.ndarray]:
    cell_index = pd.Index(cell_ids)
    selected_positions = {}
    for value in paths:
        path = Path(value)
        subset = subset_from_path(path)
        selected_ids = read_ids(path)
        positions = cell_index.get_indexer(selected_ids)
        if np.any(positions < 0):
            examples = selected_ids[positions < 0][:5]
            raise ValueError(
                f"{subset} selected cells absent from latent TSV: "
                + ", ".join(examples)
            )
        selected_positions[subset] = positions
    if set(selected_positions) != set(SUBSETS):
        raise ValueError("Exactly one selected-cell file is required per subset")
    return selected_positions


def completed_subset_gains(
    transformed: np.ndarray,
    selected_positions: dict[str, np.ndarray],
    chunk_size: int,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    dict[str, np.ndarray],
]:
    gains = {}
    gain_percentiles = {}
    accumulated = {}
    n_cells, n_dimensions = transformed.shape
    for subset in SUBSETS:
        positions = selected_positions[subset]
        coverage = np.zeros(n_dimensions)
        for start in range(0, len(positions), chunk_size):
            coverage += transformed[positions[start : start + chunk_size]].sum(axis=0)
        baseline = np.sqrt(coverage).sum()
        delta = np.empty(n_cells)
        print(
            f"Calculating completed-set marginal gains for {subset} "
            f"({len(positions):,} selected cells)",
            flush=True,
        )
        for start in range(0, n_cells, chunk_size):
            stop = min(start + chunk_size, n_cells)
            delta[start:stop] = (
                np.sqrt(transformed[start:stop] + coverage).sum(axis=1) - baseline
            )
        delta[positions] = np.nan
        unselected = np.isfinite(delta)
        percentile = np.full(n_cells, np.nan)
        percentile[unselected] = (
            pd.Series(delta[unselected]).rank(method="average", pct=True).to_numpy()
            * 100
        )
        gains[subset] = delta
        gain_percentiles[subset] = percentile
        accumulated[subset] = coverage
    return gains, gain_percentiles, accumulated


def load_retrospective_labels(path: Path, cell_ids: np.ndarray) -> np.ndarray:
    metadata = pd.read_csv(path, sep="\t", dtype={"cell_id": str})
    missing = {"cell_id", "cell_type"} - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    if metadata["cell_id"].duplicated().any():
        raise ValueError("Metadata cell IDs are not unique")
    metadata = metadata.set_index("cell_id").reindex(cell_ids)
    if metadata["cell_type"].isna().any():
        examples = metadata.index[metadata["cell_type"].isna()][:5]
        raise ValueError("Missing metadata for cells: " + ", ".join(examples))
    return metadata["cell_type"].astype(str).to_numpy()


def group_statistics(
    values: np.ndarray,
    order: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
    quantiles: bool = False,
) -> dict[str, np.ndarray]:
    result = {
        "mean": np.full(len(starts), np.nan),
        "median": np.full(len(starts), np.nan),
    }
    if quantiles:
        result["q90"] = np.full(len(starts), np.nan)
        result["q99"] = np.full(len(starts), np.nan)
        result["n_nonmissing"] = np.zeros(len(starts), dtype=int)
    ordered = values[order]
    for code, (start, stop) in enumerate(zip(starts, stops)):
        current = ordered[start:stop]
        current = current[np.isfinite(current)]
        if not len(current):
            continue
        result["mean"][code] = current.mean()
        result["median"][code] = np.median(current)
        if quantiles:
            result["q90"][code], result["q99"][code] = np.quantile(
                current, (0.90, 0.99)
            )
            result["n_nonmissing"][code] = len(current)
    return result


def build_cell_type_table(
    atlas: str,
    labels: np.ndarray,
    intrinsic: dict[str, np.ndarray],
    gains: dict[str, np.ndarray],
    gain_percentiles: dict[str, np.ndarray],
    selected_positions: dict[str, np.ndarray],
    separability_path: Path,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    cell_types, codes = np.unique(labels, return_inverse=True)
    full_n = np.bincount(codes, minlength=len(cell_types))
    starts = np.r_[0, np.cumsum(full_n)[:-1]]
    stops = np.cumsum(full_n)
    order = np.argsort(codes, kind="stable")
    table = pd.DataFrame(
        {
            "atlas": atlas,
            "cell_type": cell_types,
            "full_n": full_n,
            "full_fraction": full_n / len(labels),
        }
    )

    separability = pd.read_csv(separability_path, sep="\t")
    separability_metrics = [
        "knn_purity",
        "mean_knn_distance",
        "median_knn_distance",
        "median_knn_distance_percentile",
        "q90_knn_distance_percentile",
        "mean_knn_radius",
        "median_knn_radius",
        "median_knn_radius_percentile",
        "q90_knn_radius_percentile",
    ]
    required = {"cell_type", "full_n", *separability_metrics}
    missing = required - set(separability.columns)
    if missing:
        raise ValueError(f"Separability TSV is missing columns: {sorted(missing)}")
    indexed_separability = separability.set_index("cell_type")
    prior_n = table["cell_type"].map(indexed_separability["full_n"])
    if not np.array_equal(prior_n.to_numpy(dtype=int), full_n):
        raise ValueError("Separability output does not match current metadata")
    for metric in separability_metrics:
        table[metric] = table["cell_type"].map(indexed_separability[metric])
        if table[metric].isna().any():
            raise ValueError(f"Separability metric contains missing values: {metric}")

    for metric, values in intrinsic.items():
        stats = group_statistics(values, order, starts, stops)
        table[f"mean_{metric}"] = stats["mean"]
        table[f"median_{metric}"] = stats["median"]

    for subset in SUBSETS:
        selected_n = np.bincount(
            codes[selected_positions[subset]], minlength=len(cell_types)
        )
        residual = group_statistics(
            gains[subset], order, starts, stops, quantiles=True
        )
        residual_percentile = group_statistics(
            gain_percentiles[subset], order, starts, stops, quantiles=True
        )
        table[f"selected_n_{subset}"] = selected_n
        table[f"status_{subset}"] = np.where(
            selected_n == 0,
            "absent",
            np.where(selected_n < 10, "underpowered", "evaluable"),
        )
        table[f"unselected_n_with_gain_{subset}"] = residual["n_nonmissing"]
        table[f"mean_residual_gain_{subset}"] = residual["mean"]
        table[f"median_residual_gain_{subset}"] = residual["median"]
        table[f"q90_residual_gain_{subset}"] = residual["q90"]
        table[f"q99_residual_gain_{subset}"] = residual["q99"]
        table[f"median_residual_gain_percentile_{subset}"] = (
            residual_percentile["median"]
        )
        table[f"q90_residual_gain_percentile_{subset}"] = (
            residual_percentile["q90"]
        )
    return table.sort_values(["full_n", "cell_type"]).reset_index(drop=True), codes, full_n


def type_feature_means(
    transformed: np.ndarray, codes: np.ndarray, full_n: np.ndarray
) -> np.ndarray:
    means = np.empty((len(full_n), transformed.shape[1]))
    for dimension in range(transformed.shape[1]):
        means[:, dimension] = (
            np.bincount(codes, weights=transformed[:, dimension], minlength=len(full_n))
            / full_n
        )
    return means


def build_feature_profiles(
    atlas: str,
    table: pd.DataFrame,
    selection_features: np.ndarray,
    feature_names: list[str],
    codes: np.ndarray,
    full_n: np.ndarray,
    accumulated: dict[str, np.ndarray],
    purity_threshold: float,
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    type_means = type_feature_means(selection_features, codes, full_n)
    atlas_feature_mean = selection_features.mean(axis=0)
    all_accumulated = selection_features.sum(axis=0)

    missed = np.zeros(len(table), dtype=bool)
    for subset in SUBSETS:
        missed |= table[f"status_{subset}"].isin(("absent", "underpowered")).to_numpy()
    qualifying_table = table[(table["knn_purity"] >= purity_threshold) & missed].copy()
    qualifying_table = qualifying_table.sort_values("knn_purity", ascending=False)
    cell_type_to_code = {
        cell_type: code for code, cell_type in enumerate(np.unique(table["cell_type"]))
    }
    qualifying_codes = np.array(
        [cell_type_to_code[value] for value in qualifying_table["cell_type"]], dtype=int
    )

    rows = []
    for _, type_row in qualifying_table.iterrows():
        cell_type = type_row["cell_type"]
        code = cell_type_to_code[cell_type]
        feature_mean = type_means[code]
        relative = np.divide(
            feature_mean,
            atlas_feature_mean,
            out=np.full_like(feature_mean, np.nan),
            where=atlas_feature_mean > 0,
        )
        top_features = np.argsort(-np.nan_to_num(relative, nan=-np.inf))[:5]
        top_rank = {feature: rank + 1 for rank, feature in enumerate(top_features)}
        for subset in SUBSETS:
            coverage = accumulated[subset]
            for feature, name in enumerate(feature_names):
                rows.append(
                    {
                        "atlas": atlas,
                        "cell_type": cell_type,
                        "knn_purity": type_row["knn_purity"],
                        "subset_size": subset,
                        "retention_status": type_row[f"status_{subset}"],
                        "feature": name,
                        "feature_index": feature + 1,
                        "mean_selection_feature_type": feature_mean[feature],
                        "atlas_mean_selection_feature": atlas_feature_mean[feature],
                        "selection_feature_elevation_difference": (
                            feature_mean[feature] - atlas_feature_mean[feature]
                        ),
                        "selection_feature_relative_elevation": relative[feature],
                        "top5_relative_elevation_rank": top_rank.get(feature, np.nan),
                        "selected_accumulated_feature": coverage[feature],
                        "selected_concave_coverage": np.sqrt(coverage[feature]),
                        "selected_fraction_of_atlas_feature_mass": (
                            coverage[feature] / all_accumulated[feature]
                            if all_accumulated[feature] > 0
                            else np.nan
                        ),
                    }
                )
    columns = [
        "atlas",
        "cell_type",
        "knn_purity",
        "subset_size",
        "retention_status",
        "feature",
        "feature_index",
        "mean_selection_feature_type",
        "atlas_mean_selection_feature",
        "selection_feature_elevation_difference",
        "selection_feature_relative_elevation",
        "top5_relative_elevation_rank",
        "selected_accumulated_feature",
        "selected_concave_coverage",
        "selected_fraction_of_atlas_feature_mass",
    ]
    relative_heatmap = np.divide(
        type_means[qualifying_codes],
        atlas_feature_mean,
        out=np.full_like(type_means[qualifying_codes], np.nan),
        where=atlas_feature_mean > 0,
    )
    return (
        pd.DataFrame(rows, columns=columns),
        relative_heatmap,
        qualifying_table["cell_type"].tolist(),
    )


def point_sizes(values: pd.Series) -> np.ndarray:
    values = np.sqrt(values.to_numpy(dtype=float))
    if values.max() == values.min():
        return np.full(len(values), 80.0)
    return 30 + 230 * (values - values.min()) / (values.max() - values.min())


def annotate_missed(axis, table, x_column, y_column, status_column, threshold, maximum):
    labels = table[
        (table["knn_purity"] >= threshold)
        & table[status_column].isin(("absent", "underpowered"))
        & table[y_column].notna()
    ].sort_values("knn_purity", ascending=False).head(maximum)
    offsets = ((6, 7), (6, -11), (-6, 7), (-6, -11))
    for number, (_, row) in enumerate(labels.iterrows()):
        offset = offsets[number % len(offsets)]
        axis.annotate(
            row["cell_type"],
            (row[x_column], row[y_column]),
            xytext=offset,
            textcoords="offset points",
            ha="left" if offset[0] > 0 else "right",
            va="bottom" if offset[1] > 0 else "top",
            fontsize=7.5,
            arrowprops={"arrowstyle": "-", "lw": 0.5, "color": "#999999"},
        )


def diagnostic_scatter(
    atlas: str,
    table: pd.DataFrame,
    y_column: str,
    y_label: str,
    title: str,
    status_column: str,
    output_path: Path,
    threshold: float,
    max_labels: int,
) -> None:
    figure, axis = plt.subplots(figsize=(11, 7.2))
    sizes = point_sizes(table["full_n"])
    for current in STATUSES:
        mask = table[status_column] == current
        axis.scatter(
            table.loc[mask, "knn_purity"],
            table.loc[mask, y_column],
            s=sizes[mask.to_numpy()],
            color=STATUS_COLORS[current],
            alpha=0.76,
            edgecolor="white",
            linewidth=0.45,
            label=f"{current} (n={int(mask.sum())})",
        )
    axis.axvline(threshold, color="#777777", linestyle=":", linewidth=1)
    annotate_missed(
        axis, table, "knn_purity", y_column, status_column, threshold, max_labels
    )
    axis.set_xlim(-0.02, 1.02)
    axis.set_ylim(bottom=0)
    axis.set_xlabel("Mean 30-nearest-neighbor purity in PeakVI latent space")
    axis.set_ylabel(y_label)
    axis.set_title(f"{atlas}: {title}", fontsize=15, fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(color="#E5E5E5", linewidth=0.6, alpha=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def feature_heatmap(
    atlas: str,
    matrix: np.ndarray,
    cell_types: list[str],
    dimension_names: list[str],
    output_path: Path,
) -> None:
    height = max(4.5, min(24, 2.8 + 0.32 * max(len(cell_types), 1)))
    figure, axis = plt.subplots(figsize=(17, height))
    if not len(cell_types):
        axis.axis("off")
        axis.text(
            0.5,
            0.5,
            "No cell type met both criteria:\nkNN purity >= 0.75 and absent/underpowered",
            ha="center",
            va="center",
            transform=axis.transAxes,
        )
    else:
        from matplotlib.colors import TwoSlopeNorm

        finite_max = float(np.nanmax(matrix)) if np.isfinite(matrix).any() else 1.0
        norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=max(1.01, finite_max))
        image = axis.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap="RdBu_r",
            norm=norm,
        )
        axis.set_yticks(np.arange(len(cell_types)), labels=cell_types, fontsize=8)
        axis.set_xticks(
            np.arange(len(dimension_names)), labels=dimension_names, rotation=90, fontsize=6
        )
        axis.set_xlabel("PeakVI dimension")
        axis.set_ylabel("High-purity missed or underpowered cell type")
        colorbar = figure.colorbar(image, ax=axis, fraction=0.025, pad=0.015)
        colorbar.set_label(
            "Cell-type mean / atlas-wide mean RBF selection feature (1 = atlas mean)"
        )
    axis.set_title(
        f"{atlas}: relative feature elevation of high-purity missed populations",
        fontsize=14,
        fontweight="bold",
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if args.chunk_size < 1:
        raise ValueError("--chunk_size must be positive")
    if not 0 <= args.high_purity_threshold <= 1:
        raise ValueError("--high_purity_threshold must be between 0 and 1")
    feature_matrix_path, n_landmarks = validate_selection_summaries(
        args.selection_summaries, args.atlas
    )

    raw, cell_ids, latent_dimension_names = load_latent(Path(args.latent_tsv))
    latent_summary, _, _ = latent_distribution_summary(
        args.atlas, raw, latent_dimension_names, args.simplex_tolerance
    )
    transformed = load_selection_features(
        feature_matrix_path,
        expected_cells=len(cell_ids),
        expected_features=n_landmarks,
    )
    selection_feature_names = [
        f"RBF_landmark_{index}" for index in range(1, n_landmarks + 1)
    ]
    intrinsic = per_cell_intrinsic(raw, transformed, args.chunk_size)
    selected_positions = load_selected_positions(args.selected_cells, cell_ids)
    gains, gain_percentiles, accumulated = completed_subset_gains(
        transformed, selected_positions, args.chunk_size
    )

    # Cell-type annotations are deliberately loaded only after every per-cell
    # utility and marginal-gain quantity above has been calculated.
    labels = load_retrospective_labels(Path(args.metadata_tsv), cell_ids)
    table, codes, full_n = build_cell_type_table(
        args.atlas,
        labels,
        intrinsic,
        gains,
        gain_percentiles,
        selected_positions,
        Path(args.separability_tsv),
    )
    feature_profiles, heatmap_matrix, heatmap_types = build_feature_profiles(
        args.atlas,
        table,
        transformed,
        selection_feature_names,
        codes,
        full_n,
        accumulated,
        args.high_purity_threshold,
    )

    output_table = Path(args.output_table)
    output_summary = Path(args.output_latent_summary)
    output_profiles = Path(args.output_feature_profiles)
    figures_dir = Path(args.output_figures_dir)
    for path in (output_table, output_summary, output_profiles):
        path.parent.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_table, sep="\t", index=False)
    latent_summary.to_csv(output_summary, sep="\t", index=False)
    feature_profiles.to_csv(output_profiles, sep="\t", index=False)

    diagnostic_scatter(
        args.atlas,
        table,
        "median_singleton_score",
        "Median singleton FeatureBasedSelection score",
        "PeakVI purity versus intrinsic singleton utility",
        "status_02pct",
        figures_dir / "knn_purity_vs_singleton_score.png",
        args.high_purity_threshold,
        args.max_labels,
    )
    for subset in SUBSETS:
        diagnostic_scatter(
            args.atlas,
            table,
            f"median_residual_gain_{subset}",
            "Median counterfactual marginal gain among unselected cells",
            f"PeakVI purity versus residual gain ({SUBSET_LABELS[subset]} subset)",
            f"status_{subset}",
            figures_dir / f"knn_purity_vs_residual_gain_{subset}.png",
            args.high_purity_threshold,
            args.max_labels,
        )
    feature_heatmap(
        args.atlas,
        heatmap_matrix,
        heatmap_types,
        selection_feature_names,
        figures_dir / "high_purity_missed_feature_profiles.png",
    )

    print(f"Wrote cell-type table: {output_table}", flush=True)
    print(f"Wrote latent summary: {output_summary}", flush=True)
    print(f"Wrote feature profiles: {output_profiles}", flush=True)
    print(f"Wrote diagnostic PNGs: {figures_dir}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
