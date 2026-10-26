#!/usr/bin/env python3
"""Diagnose RBF activation breadth and Apricot marginal-gain competition.

This retrospective diagnostic uses the exact nonnegative RBF matrix consumed by
Apricot. It computes annotation-blind per-cell activation metrics and exact
counterfactual marginal gains after each completed subset, then joins cell-type
labels for aggregation and abundance-matched comparisons. It does not rerun or
change subset selection.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SUBSETS = ("02pct", "04pct", "06pct", "08pct", "10pct")
SUBSET_LABELS = ("2%", "4%", "6%", "8%", "10%")
THRESHOLDS = (0.10, 0.25, 0.50, 0.75)
INTRINSIC = (
    "total_rbf_mass", "effective_landmark_support", "max_rbf_similarity",
    "n_active_0p10", "n_active_0p25", "n_active_0p50", "n_active_0p75",
    "initial_marginal_gain",
)


def parse_args():
    p = argparse.ArgumentParser()
    for name in (
        "atlas", "feature_matrix", "landmark_metadata", "cell_names",
        "cell_metadata", "landmark_counts", "output_cell_type_tsv",
        "output_matches_tsv", "output_associations_tsv",
        "output_per_cell_npz", "output_summary_json", "output_figures_dir",
    ):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--selected_cells", nargs="+", required=True)
    p.add_argument("--cell_type_column", default="cell_type")
    p.add_argument("--chunk_size", type=int, default=10000)
    p.add_argument("--n_controls", type=int, default=3)
    p.add_argument("--random_seed", type=int, default=123)
    p.add_argument("--dpi", type=int, default=320)
    return p.parse_args()


def subset_from_path(path):
    match = re.search(r"_(02pct|04pct|06pct|08pct|10pct)\.txt$", Path(path).name)
    if not match:
        raise ValueError(f"Cannot infer subset from {path}")
    return match.group(1)


def read_names(path):
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if len(values) > 1 and values.iloc[0] in {"cell_name", "cell_id"}:
        values = values.iloc[1:]
    names = values.astype(str).to_numpy()
    if not len(names) or len(np.unique(names)) != len(names):
        raise ValueError(f"Cell names are empty or non-unique: {path}")
    return names


def read_labels(path, names, column):
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if column not in frame:
        raise ValueError(f"Metadata lacks {column}: {path}")
    id_column = next(
        candidate for candidate in ("cell_id", "cell_name", frame.columns[0])
        if candidate in frame
    )
    if frame[id_column].duplicated().any():
        raise ValueError(f"Metadata cell IDs are duplicated: {path}")
    labels = frame.set_index(id_column)[column].reindex(names)
    if labels.isna().any():
        raise ValueError(f"Metadata lacks {labels.isna().sum():,} PeakVI cells")
    return labels.astype(str).to_numpy()


def selected_positions(paths, names):
    index = pd.Index(names)
    result = {}
    for value in paths:
        subset = subset_from_path(value)
        ids = pd.read_csv(value, header=None, dtype=str).iloc[:, 0].to_numpy()
        positions = index.get_indexer(ids)
        if np.any(positions < 0):
            raise ValueError(f"{subset} contains unknown selected cells")
        if len(np.unique(positions)) != len(positions):
            raise ValueError(f"{subset} contains duplicate selected cells")
        result[subset] = positions
    if set(result) != set(SUBSETS):
        raise ValueError("Exactly one selected-cell file is required per subset")
    return result


def rank_percentiles(values):
    finite = np.isfinite(values)
    output = np.full(len(values), np.nan, dtype=np.float32)
    output[finite] = (
        pd.Series(values[finite]).rank(method="average", pct=True).to_numpy()
        * 100
    )
    return output


def group_summary(labels, metrics, selected):
    cell_types, codes = np.unique(labels, return_inverse=True)
    n_cells = np.bincount(codes, minlength=len(cell_types))
    rows = []
    for code, cell_type in enumerate(cell_types):
        mask = codes == code
        row = {
            "atlas": None,
            "cell_type": cell_type,
            "n_cells": int(n_cells[code]),
            "atlas_fraction": float(n_cells[code] / len(labels)),
        }
        for metric in INTRINSIC:
            values = metrics[metric][mask]
            row[f"median_{metric}"] = float(np.nanmedian(values))
            row[f"q90_{metric}"] = float(np.nanquantile(values, 0.90))
        enrichments = []
        retained = 0
        for subset in SUBSETS:
            selected_n = int(np.sum(mask[selected[subset]]))
            subset_n = len(selected[subset])
            enrichment = (
                (selected_n / subset_n) / (n_cells[code] / len(labels))
            )
            retained += int(selected_n > 0)
            enrichments.append(enrichment)
            row[f"selected_n_{subset}"] = selected_n
            row[f"enrichment_{subset}"] = enrichment
            row[f"status_{subset}"] = (
                "absent" if selected_n == 0 else
                "underpowered" if selected_n < 10 else "evaluable"
            )
            for stem in ("residual_gain", "residual_gain_percentile"):
                values = metrics[f"{stem}_{subset}"][mask]
                values = values[np.isfinite(values)]
                row[f"median_{stem}_{subset}"] = (
                    float(np.median(values)) if len(values) else np.nan
                )
                row[f"q90_{stem}_{subset}"] = (
                    float(np.quantile(values, 0.90)) if len(values) else np.nan
                )
        row["mean_enrichment_ratio"] = float(np.mean(enrichments))
        row["n_subsets_retained"] = retained
        row["retention_class"] = (
            "Never retained" if retained == 0 else
            "Retained in all subsets" if retained == len(SUBSETS) else
            "Retained in some subsets"
        )
        rows.append(row)
    return pd.DataFrame(rows), codes


def partial_spearman(x, y, abundance):
    rx, ry, ra = stats.rankdata(x), stats.rankdata(y), stats.rankdata(abundance)
    design = np.column_stack([np.ones(len(x)), ra])
    ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def association_table(atlas, table):
    rows = []
    tests = [
        (f"median_{metric}", "mean_enrichment_ratio")
        for metric in INTRINSIC
    ]
    tests.extend(
        (f"median_residual_gain_percentile_{subset}", f"enrichment_{subset}")
        for subset in SUBSETS
    )
    for metric, outcome in tests:
        current = table[[metric, outcome, "n_cells"]].dropna()
        x = current[metric].to_numpy()
        y = current[outcome].to_numpy()
        if len(current) < 3 or np.std(x) == 0 or np.std(y) == 0:
            rho = p_value = adjusted = np.nan
        else:
            result = stats.spearmanr(x, y)
            rho, p_value = float(result.statistic), float(result.pvalue)
            adjusted = partial_spearman(
                x, y, np.log10(current.n_cells.to_numpy())
            )
        rows.append({
            "atlas": atlas, "metric": metric, "outcome": outcome,
            "spearman_rho": rho, "spearman_p_value": p_value,
            "partial_spearman_controlling_log10_n_cells": adjusted,
            "n_cell_types": len(current),
        })
    output = pd.DataFrame(rows)
    output["spearman_fdr"] = np.nan
    finite = output["spearman_p_value"].notna()
    if finite.any():
        p_values = output.loc[finite, "spearman_p_value"].to_numpy()
        order = np.argsort(p_values)
        sorted_p = p_values[order]
        sorted_q = sorted_p * len(sorted_p) / np.arange(1, len(sorted_p) + 1)
        sorted_q = np.minimum.accumulate(sorted_q[::-1])[::-1]
        q_values = np.empty_like(sorted_q)
        q_values[order] = np.minimum(sorted_q, 1.0)
        output.loc[finite, "spearman_fdr"] = q_values
    return output


def metric_percentiles(table, columns):
    result = table.copy()
    for column in columns:
        result[f"{column}_percentile"] = (
            result[column].rank(method="average", pct=True) * 100
        )
    return result


def matched_controls(atlas, table, n_controls):
    comparison_metrics = [
        "median_total_rbf_mass", "median_effective_landmark_support",
        "median_max_rbf_similarity", "median_n_active_0p50",
        "median_initial_marginal_gain",
    ]
    ranked = metric_percentiles(table, comparison_metrics)
    rows = []
    for subset in SUBSETS:
        targets = ranked[ranked[f"selected_n_{subset}"] == 0]
        controls = ranked[ranked[f"selected_n_{subset}"] > 0]
        for _, target in targets.iterrows():
            candidates = controls.copy()
            candidates["abundance_log_distance"] = np.abs(
                np.log10(candidates.n_cells) - np.log10(target.n_cells)
            )
            candidates = candidates.sort_values(
                ["abundance_log_distance", "cell_type"]
            ).head(n_controls)
            for rank, (_, control) in enumerate(candidates.iterrows(), 1):
                row = {
                    "atlas": atlas, "subset_size": subset,
                    "missed_cell_type": target.cell_type,
                    "missed_n_cells": int(target.n_cells),
                    "control_rank": rank,
                    "control_cell_type": control.cell_type,
                    "control_n_cells": int(control.n_cells),
                    "abundance_log10_difference": float(
                        np.log10(target.n_cells) - np.log10(control.n_cells)
                    ),
                    "missed_residual_gain_percentile_median": target[
                        f"median_residual_gain_percentile_{subset}"
                    ],
                    "control_residual_gain_percentile_median": control[
                        f"median_residual_gain_percentile_{subset}"
                    ],
                }
                for metric in comparison_metrics:
                    row[f"missed_{metric}"] = target[metric]
                    row[f"control_{metric}"] = control[metric]
                    row[f"difference_{metric}_percentile"] = (
                        target[f"{metric}_percentile"]
                        - control[f"{metric}_percentile"]
                    )
                row["difference_residual_gain_percentile_median"] = (
                    row["missed_residual_gain_percentile_median"]
                    - row["control_residual_gain_percentile_median"]
                )
                rows.append(row)
    return pd.DataFrame(rows)


def plot_activation(atlas, table, path, dpi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [
        ("median_total_rbf_mass", "Median total RBF mass"),
        ("median_effective_landmark_support", "Median effective landmark support"),
        ("median_max_rbf_similarity", "Median maximum RBF similarity"),
        ("median_n_active_0p50", "Median landmarks with similarity >= 0.50"),
    ]
    colors = {
        "Retained in all subsets": "#287C8E",
        "Retained in some subsets": "#E69F00",
        "Never retained": "#C43C39",
    }
    figure, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=dpi)
    for axis, (metric, label) in zip(axes.flat, panels):
        for group, current in table.groupby("retention_class"):
            axis.scatter(
                current.mean_enrichment_ratio, current[metric],
                s=25 + 18 * np.log10(current.n_cells),
                color=colors[group], alpha=0.72, edgecolor="white",
                linewidth=0.3, label=group,
            )
        for row in table[
            table.retention_class != "Retained in all subsets"
        ].itertuples():
            axis.annotate(
                row.cell_type,
                (row.mean_enrichment_ratio, getattr(row, metric)),
                xytext=(3, 3), textcoords="offset points", fontsize=6.5,
            )
        axis.axvline(1, color="black", linestyle="--", linewidth=0.8)
        axis.set_xlabel("Mean subset enrichment ratio")
        axis.set_ylabel(label)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(alpha=0.2)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    figure.suptitle(
        f"{atlas}: RBF activation breadth versus cell-type enrichment",
        fontsize=18, fontweight="bold",
    )
    figure.tight_layout(rect=[0, 0.06, 1, 0.95])
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def plot_residual_heatmap(atlas, table, path, dpi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ordered = table.sort_values(
        ["mean_enrichment_ratio", "cell_type"], ascending=[False, True]
    )
    matrix = ordered[
        [f"median_residual_gain_percentile_{subset}" for subset in SUBSETS]
    ].to_numpy()
    height = max(8, 0.22 * len(ordered) + 2.5)
    figure, axis = plt.subplots(figsize=(9, height), dpi=dpi)
    image = axis.imshow(
        matrix, aspect="auto", cmap="viridis", vmin=0, vmax=100,
        interpolation="nearest",
    )
    axis.set_xticks(range(len(SUBSETS)), SUBSET_LABELS)
    axis.set_yticks(
        range(len(ordered)),
        [f"{row.cell_type} ({row.mean_enrichment_ratio:.2f}x)"
         for row in ordered.itertuples()],
        fontsize=7,
    )
    axis.set_xlabel("Completed subset")
    axis.set_ylabel("Cell type, ordered by mean enrichment")
    axis.set_title(
        f"{atlas}: residual marginal-gain percentile among unselected cells",
        fontweight="bold",
    )
    colorbar = figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02)
    colorbar.set_label("Median percentile (100 = highest residual gain)")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def plot_matches(atlas, matches, path, dpi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    columns = [
        "difference_median_total_rbf_mass_percentile",
        "difference_median_effective_landmark_support_percentile",
        "difference_median_max_rbf_similarity_percentile",
        "difference_median_n_active_0p50_percentile",
        "difference_median_initial_marginal_gain_percentile",
        "difference_residual_gain_percentile_median",
    ]
    labels = [
        "Total mass", "Effective support", "Max similarity",
        "Active >=0.5", "Initial gain", "Residual gain",
    ]
    if matches.empty:
        matrix = np.zeros((1, len(columns)))
        row_labels = ["No missed cell types"]
    else:
        grouped = matches.groupby(
            ["subset_size", "missed_cell_type"], sort=False
        )[columns].mean().reset_index()
        grouped["_order"] = grouped.subset_size.map(
            {value: index for index, value in enumerate(SUBSETS)}
        )
        grouped = grouped.sort_values(["_order", "missed_cell_type"])
        matrix = grouped[columns].to_numpy()
        row_labels = [
            f"{row.missed_cell_type} ({row.subset_size})"
            for row in grouped.itertuples()
        ]
    limit = max(10, float(np.nanmax(np.abs(matrix))))
    height = max(5, 0.28 * len(row_labels) + 2.5)
    figure, axis = plt.subplots(figsize=(11, height), dpi=dpi)
    image = axis.imshow(
        matrix, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit,
        interpolation="nearest",
    )
    axis.set_xticks(range(len(labels)), labels, rotation=30, ha="right")
    axis.set_yticks(range(len(row_labels)), row_labels, fontsize=7)
    axis.set_title(
        f"{atlas}: missed types versus abundance-matched retained controls",
        fontweight="bold",
    )
    colorbar = figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02)
    colorbar.set_label("Missed minus matched-control percentile points")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    if args.chunk_size < 1 or args.n_controls < 1:
        raise ValueError("Chunk size and number of controls must be positive")

    names = read_names(args.cell_names)
    labels = read_labels(
        args.cell_metadata, names, args.cell_type_column
    )
    selected = selected_positions(args.selected_cells, names)
    features = np.load(args.feature_matrix, mmap_mode="r")
    if features.ndim != 2 or features.shape[0] != len(names):
        raise ValueError("RBF feature matrix and cell names are incompatible")
    with open(args.landmark_metadata, encoding="utf-8") as handle:
        metadata = json.load(handle)
    if metadata.get("atlas") != args.atlas:
        raise ValueError("Atlas mismatch in landmark metadata")
    if metadata.get("feature_construction_method") != "peakvi_rbf_landmarks":
        raise ValueError("Metadata does not describe RBF landmarks")
    if metadata.get("selection_function") != "sqrt":
        raise ValueError("Diagnostic requires the production sqrt objective")
    if metadata.get("min_shift_applied") is not False:
        raise ValueError("Diagnostic requires no min-shift")
    if tuple(metadata["feature_matrix_shape"]) != features.shape:
        raise ValueError("RBF matrix shape disagrees with metadata")

    coverage = {}
    for subset in SUBSETS:
        total = np.zeros(features.shape[1], dtype=np.float64)
        positions = selected[subset]
        for start in range(0, len(positions), args.chunk_size):
            total += np.asarray(
                features[positions[start:start + args.chunk_size]],
                dtype=np.float32,
            ).sum(0, dtype=np.float64)
        coverage[subset] = total.astype(np.float32)
        print(f"Built completed {subset} coverage", flush=True)

    n_cells = len(names)
    metrics = {
        name: np.empty(n_cells, dtype=np.float32)
        for name in INTRINSIC
    }
    for subset in SUBSETS:
        metrics[f"residual_gain_{subset}"] = np.empty(
            n_cells, dtype=np.float32
        )

    for start in range(0, n_cells, args.chunk_size):
        stop = min(start + args.chunk_size, n_cells)
        block = np.asarray(features[start:stop], dtype=np.float32)
        mass = block.sum(1, dtype=np.float32)
        xlogx = np.zeros_like(block)
        np.log(block, out=xlogx, where=block > 0)
        xlogx *= block
        entropy = np.log(mass) - xlogx.sum(1) / mass
        metrics["total_rbf_mass"][start:stop] = mass
        metrics["effective_landmark_support"][start:stop] = np.exp(entropy)
        metrics["max_rbf_similarity"][start:stop] = block.max(1)
        for threshold in THRESHOLDS:
            threshold_label = f"{threshold:.2f}".replace(".", "p")
            key = f"n_active_{threshold_label}"
            metrics[key][start:stop] = (block >= threshold).sum(1)
        metrics["initial_marginal_gain"][start:stop] = np.sqrt(block).sum(1)
        for subset in SUBSETS:
            current_coverage = coverage[subset]
            baseline = np.sqrt(current_coverage, dtype=np.float32).sum(
                dtype=np.float64
            )
            metrics[f"residual_gain_{subset}"][start:stop] = (
                np.sqrt(block + current_coverage).sum(1) - baseline
            )
        print(f"Processed RBF rows {start:,}-{stop:,} of {n_cells:,}", flush=True)

    metrics["initial_marginal_gain_percentile"] = rank_percentiles(
        metrics["initial_marginal_gain"]
    )
    for subset in SUBSETS:
        key = f"residual_gain_{subset}"
        metrics[key][selected[subset]] = np.nan
        metrics[f"residual_gain_percentile_{subset}"] = rank_percentiles(
            metrics[key]
        )

    table, _ = group_summary(labels, metrics, selected)
    table["atlas"] = args.atlas
    landmark_counts = pd.read_csv(args.landmark_counts, sep="\t")
    required = {"cell_type", "n_landmarks"}
    if required - set(landmark_counts):
        raise ValueError("Landmark-count table lacks required columns")
    table = table.merge(
        landmark_counts[["cell_type", "n_landmarks"]],
        on="cell_type", how="left", validate="one_to_one",
    )
    if table.n_landmarks.isna().any():
        raise ValueError("Landmark-count table lacks annotated cell types")
    table = table.sort_values(
        ["mean_enrichment_ratio", "cell_type"], ascending=[False, True]
    )
    associations = association_table(args.atlas, table)
    matches = matched_controls(args.atlas, table, args.n_controls)

    for value in (
        args.output_cell_type_tsv, args.output_matches_tsv,
        args.output_associations_tsv, args.output_per_cell_npz,
        args.output_summary_json,
    ):
        Path(value).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_cell_type_tsv, sep="\t", index=False)
    matches.to_csv(args.output_matches_tsv, sep="\t", index=False)
    associations.to_csv(args.output_associations_tsv, sep="\t", index=False)
    np.savez(
        args.output_per_cell_npz,
        cell_name_order_file=np.array([args.cell_names]),
        **{key: value for key, value in metrics.items()},
    )

    figures = Path(args.output_figures_dir)
    plot_activation(
        args.atlas, table, figures / "rbf_activation_vs_enrichment.png",
        args.dpi,
    )
    plot_residual_heatmap(
        args.atlas, table,
        figures / "residual_gain_percentile_heatmap.png", args.dpi,
    )
    plot_matches(
        args.atlas, matches,
        figures / "missed_vs_abundance_matched_controls.png", args.dpi,
    )
    summary = {
        "atlas": args.atlas,
        "n_cells": n_cells,
        "n_rbf_landmarks": features.shape[1],
        "selection_objective": "sum_j sqrt(sum_selected_i x_ij)",
        "activation_thresholds": list(THRESHOLDS),
        "n_abundance_matched_controls": args.n_controls,
        "per_cell_metrics_are_annotation_blind": True,
        "labels_joined_after_per_cell_calculation": True,
        "n_cell_types": len(table),
        "n_never_retained": int(
            (table.retention_class == "Never retained").sum()
        ),
        "n_zero_landmark_types": int((table.n_landmarks == 0).sum()),
        "figures_directory": str(figures),
    }
    with open(args.output_summary_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(f"Wrote cell-type summary: {args.output_cell_type_tsv}", flush=True)
    print(f"Wrote matched controls: {args.output_matches_tsv}", flush=True)
    print(f"Wrote associations: {args.output_associations_tsv}", flush=True)
    print(f"Wrote figures: {figures}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
