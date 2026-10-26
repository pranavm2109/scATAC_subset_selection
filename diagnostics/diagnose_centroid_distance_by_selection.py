#!/usr/bin/env python3
"""Test whether own-landmark distances explain cell-type underselection.

Centroids are assigned by their nearest observed cell's annotation, matching
diagnose_landmark_celltype_coverage.py. Statistics use every assigned centroid;
three reproducibly sampled centroids are shown per cell type. Distances use the
unmodified Gaussian PeakVI space and are normalized by the production RBF sigma.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SUBSETS = ["02pct", "04pct", "06pct", "08pct", "10pct"]


def parse_args():
    p = argparse.ArgumentParser()
    for name in (
        "atlas", "latent_npy", "centroids", "landmark_metadata",
        "landmark_assignments", "cell_names", "cell_metadata",
        "proportions_tsv", "output_cell_type_tsv", "output_centroid_tsv",
        "output_statistics_tsv", "output_summary_png", "output_violin_dir",
        "output_manifest_tsv",
    ):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--cell_type_column", default="cell_type")
    p.add_argument("--max_cells_per_type", type=int, default=20000)
    p.add_argument("--distance_chunk_size", type=int, default=5000)
    p.add_argument("--n_permutations", type=int, default=10000)
    p.add_argument("--n_bootstrap", type=int, default=5000)
    p.add_argument("--random_seed", type=int, default=123)
    p.add_argument("--dpi", type=int, default=320)
    return p.parse_args()


def read_names(path):
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if len(values) > 1 and values.iloc[0] in {"cell_name", "cell_id"}:
        values = values.iloc[1:]
    names = values.astype(str).to_numpy()
    if not len(names) or len(np.unique(names)) != len(names):
        raise ValueError(f"Cell names are empty or non-unique: {path}")
    return names


def read_types(path, names, column):
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if column not in frame:
        raise ValueError(f"Metadata lacks {column}: {path}")
    id_column = next(
        candidate for candidate in ("cell_id", "cell_name", frame.columns[0])
        if candidate in frame
    )
    if frame[id_column].duplicated().any():
        raise ValueError(f"Metadata cell identifiers are duplicated: {path}")
    aligned = frame.set_index(id_column)[column].reindex(names)
    if aligned.isna().any():
        raise ValueError(f"Metadata lacks {aligned.isna().sum():,} PeakVI cells")
    return aligned.astype(str).to_numpy()


def rng_for(seed, atlas, cell_type, purpose):
    token = f"{seed}|{atlas}|{cell_type}|{purpose}".encode()
    derived = int.from_bytes(hashlib.sha256(token).digest()[:8], "little")
    return np.random.default_rng(derived)


def euclidean_distances(points, centers):
    squared = (
        np.square(points).sum(1, keepdims=True)
        + np.square(centers).sum(1)[None, :]
        - 2 * points @ centers.T
    )
    np.maximum(squared, 0, out=squared)
    np.sqrt(squared, out=squared)
    return squared


def describe(values):
    if not len(values):
        return {"mean": np.nan, "p10": np.nan, "median": np.nan, "p90": np.nan}
    p10, median, p90 = np.quantile(values, [0.1, 0.5, 0.9])
    return {
        "mean": float(values.mean()), "p10": float(p10),
        "median": float(median), "p90": float(p90),
    }


def safe_slug(value):
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_.") or "cell_type"
    digest = hashlib.sha1(value.encode()).hexdigest()[:8]
    return f"{stem[:90]}_{digest}"


def cliffs_delta(x, y):
    return float(np.sign(x[:, None] - y[None, :]).mean())


def partial_spearman(x, y, covariates):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    ranked_covariates = np.column_stack([
        stats.rankdata(covariates[:, index])
        for index in range(covariates.shape[1])
    ])
    design = np.column_stack([np.ones(len(x)), ranked_covariates])
    residual_x = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    residual_y = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    return float(np.corrcoef(residual_x, residual_y)[0, 1])


def inference(x, y, n_permutations, n_bootstrap, rng):
    observed = float(stats.spearmanr(x, y).statistic)
    exceedances = 0
    for _ in range(n_permutations):
        permuted = stats.spearmanr(x, rng.permutation(y)).statistic
        exceedances += int(abs(permuted) >= abs(observed))
    bootstraps = []
    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(x), len(x))
        if len(np.unique(x[indices])) > 1 and len(np.unique(y[indices])) > 1:
            bootstraps.append(stats.spearmanr(x[indices], y[indices]).statistic)
    low, high = (
        np.quantile(bootstraps, [0.025, 0.975])
        if bootstraps else (np.nan, np.nan)
    )
    return observed, (exceedances + 1) / (n_permutations + 1), low, high


def plot_violins(atlas, cell_type, ids, values, metadata, output, dpi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(12.6, 5.4), dpi=dpi)
    present = np.concatenate(values) if values else np.array([])
    upper = (
        max(float(np.quantile(present, 0.995) * 1.08), 0.25)
        if len(present) else 1
    )
    for index, axis in enumerate(axes):
        if index < len(values):
            violin = axis.violinplot(
                values[index], positions=[1], widths=0.72, showextrema=False
            )
            for body in violin["bodies"]:
                body.set_facecolor("#287C8E")
                body.set_edgecolor("#163F49")
                body.set_alpha(0.82)
            median = np.median(values[index])
            axis.scatter(1, median, s=32, color="#D95F02", zorder=3)
            axis.set_title(f"Centroid {ids[index]}\nmedian = {median:.3f}")
        else:
            axis.text(
                0.5, 0.5, "No additional\nassigned centroid",
                ha="center", va="center", transform=axis.transAxes,
            )
            axis.set_title(f"Centroid {index + 1}: unavailable")
        axis.set_xlim(0.5, 1.5)
        axis.set_ylim(0, upper)
        axis.set_xticks([])
        axis.grid(axis="y", alpha=0.2)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    axes[0].set_ylabel("Euclidean distance / atlas RBF sigma")
    figure.suptitle(
        f"{atlas}: {cell_type}\nCell-to-own-landmark distances",
        fontsize=15, fontweight="bold",
    )
    figure.text(
        0.5, 0.89,
        f"assigned centroids={metadata['n_centroids']}; "
        f"cells={metadata['n_cells']:,}; sampled={metadata['n_sampled']:,}; "
        f"mean enrichment={metadata['enrichment']:.3f}x",
        ha="center", fontsize=9,
    )
    if not values:
        axes[0].text(
            0.5, 0.2, "No direct landmark anchor", ha="center",
            transform=axes[0].transAxes, color="#A33A2B", fontweight="bold",
        )
    figure.tight_layout(rect=[0, 0.02, 1, 0.84])
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def plot_overview(atlas, table, tests, output, dpi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eligible = table.dropna(subset=["median_normalized_distance"])
    colors = {
        "Retained in all subsets": "#287C8E",
        "Retained in some subsets": "#E69F00",
        "Never retained": "#C43C39",
    }
    figure, axes = plt.subplots(1, 2, figsize=(15, 6), dpi=dpi)
    for label, current in eligible.groupby("retention_class"):
        axes[0].scatter(
            current["median_normalized_distance"],
            current["mean_enrichment_ratio"],
            s=28 + 22 * np.log10(current["n_cells"]),
            color=colors[label], alpha=0.72, edgecolor="white",
            linewidth=0.35, label=label,
        )
    axes[0].axhline(1, color="black", linestyle="--", linewidth=0.9)
    for row in eligible[
        eligible["retention_class"] != "Retained in all subsets"
    ].itertuples():
        axes[0].annotate(
            row.cell_type,
            (row.median_normalized_distance, row.mean_enrichment_ratio),
            xytext=(3, 3), textcoords="offset points", fontsize=7,
        )
    axes[0].set_xlabel(
        "Median cell-to-assigned-centroid distance / RBF sigma"
    )
    axes[0].set_ylabel("Mean subset enrichment ratio")
    axes[0].set_title("Distance versus subset enrichment", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)

    under = eligible.loc[
        eligible.mean_enrichment_ratio < 1, "median_normalized_distance"
    ].to_numpy()
    over = eligible.loc[
        eligible.mean_enrichment_ratio > 1, "median_normalized_distance"
    ].to_numpy()
    for position, values, color in (
        (1, under, "#D95F02"), (2, over, "#287C8E")
    ):
        if len(values):
            box = axes[1].boxplot(
                [values], positions=[position], widths=0.55,
                patch_artist=True, showfliers=False,
            )
            box["boxes"][0].set_facecolor(color)
            box["boxes"][0].set_alpha(0.65)
            jitter = np.random.default_rng(991 + position).normal(
                position, 0.045, len(values)
            )
            axes[1].scatter(
                jitter, values, s=16, color="#30383D", alpha=0.55
            )
    axes[1].set_xticks(
        [1, 2],
        [f"Underselected\n(n={len(under)})", f"Overselected\n(n={len(over)})"],
    )
    axes[1].set_ylabel("Cell-type median normalized distance")
    axes[1].set_title("Cell-type-level comparison", fontweight="bold")
    result = tests[tests.test == "spearman_distance_vs_enrichment"].iloc[0]
    axes[1].text(
        0.04, 0.96,
        f"Spearman rho={result.estimate:.3f}\n"
        f"permutation p={result.p_value:.4g}\n"
        f"95% bootstrap CI [{result.ci_low:.3f}, {result.ci_high:.3f}]",
        transform=axes[1].transAxes, va="top", fontsize=9,
        bbox=dict(facecolor="white", edgecolor="#CCCCCC", alpha=0.9),
    )
    for axis in axes:
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle(
        f"{atlas}: do landmark distances explain underselection?",
        fontsize=18, fontweight="bold",
    )
    figure.tight_layout(rect=[0, 0, 1, 0.94])
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    if min(args.max_cells_per_type, args.distance_chunk_size) < 1:
        raise ValueError("Sampling and chunk sizes must be positive")

    latent = np.load(args.latent_npy, mmap_mode="r")
    centers = np.load(args.centroids, mmap_mode="r")
    names = read_names(args.cell_names)
    cell_types = read_types(args.cell_metadata, names, args.cell_type_column)
    if latent.shape[0] != len(names) or latent.shape[1] != centers.shape[1]:
        raise ValueError("Latent, cell-name, and centroid shapes are incompatible")

    with open(args.landmark_metadata, encoding="utf-8") as handle:
        landmark_metadata = json.load(handle)
    if landmark_metadata.get("atlas") != args.atlas:
        raise ValueError("Atlas mismatch in landmark metadata")
    if landmark_metadata.get("peakvi_latent_distribution") != "normal":
        raise ValueError("Expected Gaussian PeakVI landmarks")
    if landmark_metadata.get("min_shift_applied") is not False:
        raise ValueError("Landmark metadata indicates a min-shift")
    sigma = float(landmark_metadata["rbf_sigma"])
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError(f"Invalid RBF sigma: {sigma}")

    assignments = pd.read_csv(args.landmark_assignments, sep="\t")
    required = {"atlas", "landmark_id", "nearest_cell_type"}
    if required - set(assignments):
        raise ValueError("Landmark assignments lack required columns")
    if set(assignments["atlas"].astype(str)) != {args.atlas}:
        raise ValueError("Atlas mismatch in landmark assignments")
    assignments["landmark_id"] = pd.to_numeric(
        assignments.landmark_id, errors="raise"
    ).astype(int)
    if assignments.landmark_id.duplicated().any():
        raise ValueError("Landmark IDs are duplicated")
    if assignments.landmark_id.min() < 1 or assignments.landmark_id.max() > len(centers):
        raise ValueError("Landmark IDs are outside the centroid matrix")

    proportions = pd.read_csv(args.proportions_tsv, sep="\t")
    ratio_columns = [f"ratio_{subset}_vs_atlas" for subset in SUBSETS]
    count_columns = [f"subset_{subset}_count" for subset in SUBSETS]
    if {"cell_type", *ratio_columns, *count_columns} - set(proportions):
        raise ValueError("Proportions table lacks required subset columns")
    proportions.cell_type = proportions.cell_type.astype(str)
    if proportions.cell_type.duplicated().any():
        raise ValueError("Proportions table contains duplicate cell types")
    proportions = proportions.set_index("cell_type")

    violin_dir = Path(args.output_violin_dir)
    violin_dir.mkdir(parents=True, exist_ok=True)
    type_rows, centroid_rows, manifest_rows = [], [], []
    for cell_type in sorted(np.unique(cell_types)):
        if cell_type not in proportions.index:
            raise ValueError(f"Cell type missing from proportions: {cell_type}")
        indices = np.flatnonzero(cell_types == cell_type)
        cell_rng = rng_for(args.random_seed, args.atlas, cell_type, "cells")
        sampled = (
            np.sort(cell_rng.choice(
                indices, args.max_cells_per_type, replace=False
            ))
            if len(indices) > args.max_cells_per_type else indices
        )
        points = np.asarray(latent[sampled], dtype=np.float32)
        type_mean = points.mean(0, dtype=np.float64)
        radius = float(np.sqrt(np.mean(np.square(points - type_mean).sum(1))))
        assigned = np.sort(assignments.loc[
            assignments.nearest_cell_type.astype(str) == cell_type,
            "landmark_id",
        ].to_numpy(dtype=int))
        enrichment = float(np.nanmean(pd.to_numeric(
            proportions.loc[cell_type, ratio_columns], errors="coerce"
        )))
        subset_counts = pd.to_numeric(
            proportions.loc[cell_type, count_columns], errors="raise"
        ).to_numpy()
        retained = int((subset_counts > 0).sum())
        retention_class = (
            "Never retained" if retained == 0 else
            "Retained in all subsets" if retained == len(SUBSETS) else
            "Retained in some subsets"
        )

        chosen = np.array([], dtype=int)
        violin_values = []
        if len(assigned):
            chosen = np.sort(rng_for(
                args.random_seed, args.atlas, cell_type, "centroids"
            ).choice(assigned, min(3, len(assigned)), replace=False))
            blocks = []
            own_centers = np.asarray(centers[assigned - 1], dtype=np.float32)
            for start in range(0, len(points), args.distance_chunk_size):
                blocks.append(
                    euclidean_distances(
                        points[start:start + args.distance_chunk_size],
                        own_centers,
                    ) / sigma
                )
            matrix = np.vstack(blocks)
            for column, landmark_id in enumerate(assigned):
                normalized = matrix[:, column]
                summary = describe(normalized)
                offset = float(np.linalg.norm(
                    np.asarray(centers[landmark_id - 1]) - type_mean
                ))
                centroid_rows.append({
                    "atlas": args.atlas,
                    "cell_type": cell_type,
                    "landmark_id": landmark_id,
                    "selected_for_violin": landmark_id in chosen,
                    "n_cells": len(indices),
                    "n_cells_sampled": len(sampled),
                    "mean_normalized_distance": summary["mean"],
                    "median_normalized_distance": summary["median"],
                    "p90_normalized_distance": summary["p90"],
                    "median_raw_distance": summary["median"] * sigma,
                    "median_rbf_similarity": float(np.median(
                        np.exp(-0.5 * np.square(normalized))
                    )),
                    "centroid_to_type_mean_distance": offset,
                    "centroid_to_type_mean_normalized": offset / sigma,
                    "within_type_rms_radius": radius,
                    "within_type_rms_radius_normalized": radius / sigma,
                })
            flat = matrix.ravel()
            combined = describe(flat)
            median_similarity = float(np.median(
                np.exp(-0.5 * np.square(flat))
            ))
            violin_values = [
                matrix[:, np.where(assigned == landmark_id)[0][0]]
                for landmark_id in chosen
            ]
        else:
            combined = describe(np.array([]))
            median_similarity = np.nan

        type_rows.append({
            "atlas": args.atlas,
            "cell_type": cell_type,
            "n_cells": len(indices),
            "n_cells_sampled": len(sampled),
            "n_assigned_centroids": len(assigned),
            "centroid_coverage_status": (
                "at_least_three" if len(assigned) >= 3 else
                "one_or_two" if len(assigned) else "zero"
            ),
            "mean_enrichment_ratio": enrichment,
            "selection_group": (
                "underselected" if enrichment < 1 else
                "overselected" if enrichment > 1 else "balanced"
            ),
            "n_subsets_retained": retained,
            "retention_class": retention_class,
            "mean_normalized_distance": combined["mean"],
            "median_normalized_distance": combined["median"],
            "p10_normalized_distance": combined["p10"],
            "p90_normalized_distance": combined["p90"],
            "median_raw_distance": combined["median"] * sigma,
            "median_rbf_similarity": median_similarity,
            "within_type_rms_radius": radius,
            "within_type_rms_radius_normalized": radius / sigma,
        })
        output = violin_dir / f"{safe_slug(cell_type)}.png"
        plot_violins(
            args.atlas, cell_type, chosen.tolist(), violin_values,
            {
                "n_centroids": len(assigned), "n_cells": len(indices),
                "n_sampled": len(sampled), "enrichment": enrichment,
            },
            output, args.dpi,
        )
        manifest_rows.append({
            "atlas": args.atlas,
            "cell_type": cell_type,
            "mean_enrichment_ratio": enrichment,
            "n_assigned_centroids": len(assigned),
            "violin_landmark_ids": ";".join(map(str, chosen)),
            "violin_png": str(output),
        })
        print(
            f"{cell_type}: cells={len(indices):,}, centroids={len(assigned)}, "
            f"enrichment={enrichment:.3f}", flush=True,
        )

    type_table = pd.DataFrame(type_rows).sort_values(
        ["mean_enrichment_ratio", "cell_type"], ascending=[False, True]
    )
    centroid_table = pd.DataFrame(centroid_rows)
    manifest = pd.DataFrame(manifest_rows).sort_values(
        ["mean_enrichment_ratio", "cell_type"], ascending=[False, True]
    )
    eligible = type_table.dropna(
        subset=["median_normalized_distance", "mean_enrichment_ratio"]
    )
    x = eligible.median_normalized_distance.to_numpy()
    y = eligible.mean_enrichment_ratio.to_numpy()
    if len(eligible) < 3:
        raise ValueError("Fewer than three cell types have assigned landmarks")
    rng = np.random.default_rng(args.random_seed)
    rho, permutation_p, low, high = inference(
        x, y, args.n_permutations, args.n_bootstrap, rng
    )
    partial = partial_spearman(
        x, y, np.column_stack([
            np.log10(eligible.n_cells.to_numpy()),
            eligible.n_assigned_centroids.to_numpy(),
        ])
    )
    under = eligible.loc[
        eligible.mean_enrichment_ratio < 1, "median_normalized_distance"
    ].to_numpy()
    over = eligible.loc[
        eligible.mean_enrichment_ratio > 1, "median_normalized_distance"
    ].to_numpy()
    if not len(under) or not len(over):
        raise ValueError("Need both underselected and overselected cell types")
    mann_whitney = stats.mannwhitneyu(under, over, alternative="two-sided")
    tests = pd.DataFrame([
        {
            "atlas": args.atlas,
            "test": "spearman_distance_vs_enrichment",
            "estimate": rho, "p_value": permutation_p,
            "ci_low": low, "ci_high": high,
            "n_cell_types": len(eligible),
            "details": "cell-type bootstrap CI and label-permutation p-value",
        },
        {
            "atlas": args.atlas,
            "test": "partial_spearman_controlling_abundance_and_landmarks",
            "estimate": partial, "p_value": np.nan,
            "ci_low": np.nan, "ci_high": np.nan,
            "n_cell_types": len(eligible),
            "details": "controls rank log10(n_cells) and n_assigned_centroids",
        },
        {
            "atlas": args.atlas,
            "test": "mann_whitney_underselected_vs_overselected",
            "estimate": float(mann_whitney.statistic),
            "p_value": float(mann_whitney.pvalue),
            "ci_low": np.nan, "ci_high": np.nan,
            "n_cell_types": len(under) + len(over),
            "details": f"underselected n={len(under)}; overselected n={len(over)}",
        },
        {
            "atlas": args.atlas,
            "test": "cliffs_delta_underselected_minus_overselected",
            "estimate": cliffs_delta(under, over), "p_value": np.nan,
            "ci_low": np.nan, "ci_high": np.nan,
            "n_cell_types": len(under) + len(over),
            "details": "positive means larger distances in underselected types",
        },
    ])

    for value in (
        args.output_cell_type_tsv, args.output_centroid_tsv,
        args.output_statistics_tsv, args.output_manifest_tsv,
    ):
        Path(value).parent.mkdir(parents=True, exist_ok=True)
    type_table.to_csv(args.output_cell_type_tsv, sep="\t", index=False)
    centroid_table.to_csv(args.output_centroid_tsv, sep="\t", index=False)
    tests.to_csv(args.output_statistics_tsv, sep="\t", index=False)
    manifest.to_csv(args.output_manifest_tsv, sep="\t", index=False)
    plot_overview(
        args.atlas, type_table, tests, Path(args.output_summary_png), args.dpi
    )
    print(f"Wrote overview: {args.output_summary_png}", flush=True)
    print(f"Wrote violin panels: {args.output_violin_dir}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
