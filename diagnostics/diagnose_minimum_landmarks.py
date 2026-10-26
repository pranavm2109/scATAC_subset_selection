#!/usr/bin/env python3
"""Estimate stable landmark requirements for annotated cell populations.

This diagnostic refits annotation-blind MiniBatchKMeans landmarks across a grid
of landmark counts and random seeds. Cell-type labels are introduced only after
each fit to evaluate which populations the landmarks represent.

Three retrospective capture definitions are reported:

1. Anchor capture: at least one centroid's nearest observed cell has the type.
2. Dominant capture: the type is the largest population in at least one
   centroid's Voronoi cluster.
3. Enriched capture: at least one cluster contains a minimum number of cells of
   the type and exceeds a specified enrichment over its atlas-wide frequency.

The primary estimate is the minimum stable anchor count: the first tested
landmark count that meets the requested seed reproducibility threshold and
continues to meet it at every larger tested count. This avoids choosing a
one-off favorable K-means realization.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import time

import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans


def parse_integer_list(value: str, name: str) -> list[int]:
    try:
        parsed = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"{name} must be a comma-separated list of integers"
        ) from error
    if not parsed or any(item < 1 for item in parsed):
        raise argparse.ArgumentTypeError(f"{name} values must all be positive")
    if len(set(parsed)) != len(parsed):
        raise argparse.ArgumentTypeError(f"{name} contains duplicate values")
    return sorted(parsed)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate the minimum stable number of PeakVI K-means landmarks "
            "needed to cover each annotated cell type."
        )
    )
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--latent_npy", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--cell_metadata", required=True)
    parser.add_argument("--output_run_details_tsv", required=True)
    parser.add_argument("--output_cell_type_tsv", required=True)
    parser.add_argument("--output_atlas_tsv", required=True)
    parser.add_argument("--output_heatmap_png", required=True)
    parser.add_argument("--output_curve_png", required=True)
    parser.add_argument("--cell_type_column", default="cell_type")
    parser.add_argument(
        "--landmark_counts",
        default="256,512,768,1024,1536,2048,3072,4096",
        help="Comma-separated landmark counts to test.",
    )
    parser.add_argument(
        "--random_seeds",
        default="123,456,789",
        help="Comma-separated K-means random seeds.",
    )
    parser.add_argument(
        "--required_seed_fraction",
        type=float,
        default=1.0,
        help=(
            "Required capture fraction across seeds for stable coverage "
            "[default: 1.0]."
        ),
    )
    parser.add_argument(
        "--min_type_cells_enriched",
        type=int,
        default=5,
        help="Minimum cells of a type in an enriched landmark cluster.",
    )
    parser.add_argument(
        "--min_cluster_enrichment",
        type=float,
        default=5.0,
        help="Minimum cell-type enrichment for enriched capture.",
    )
    parser.add_argument("--kmeans_batch_size", type=int, default=8192)
    parser.add_argument("--kmeans_max_iter", type=int, default=100)
    parser.add_argument("--kmeans_n_init", type=int, default=3)
    parser.add_argument("--kmeans_max_no_improvement", type=int, default=20)
    parser.add_argument("--kmeans_reassignment_ratio", type=float, default=0.01)
    parser.add_argument("--distance_chunk_size", type=int, default=100000)
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
        examples = cell_names[aligned.isna().to_numpy()][:10]
        raise ValueError(
            f"Metadata is missing {aligned.isna().sum():,} PeakVI cells. "
            f"Examples: {examples.tolist()}"
        )
    return aligned.astype(str)


def nearest_anchor_cells(
    latent: np.ndarray,
    labels: np.ndarray,
    centroids: np.ndarray,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_landmarks = centroids.shape[0]
    best_indices = np.full(n_landmarks, -1, dtype=np.int64)
    best_squared_distance = np.full(n_landmarks, np.inf, dtype=np.float64)

    for start in range(0, latent.shape[0], chunk_size):
        stop = min(start + chunk_size, latent.shape[0])
        block = np.asarray(latent[start:stop], dtype=np.float32)
        block_labels = labels[start:stop]
        differences = block - centroids[block_labels]
        squared_distance = np.einsum(
            "ij,ij->i", differences, differences, optimize=True
        )
        for landmark in np.unique(block_labels):
            positions = np.flatnonzero(block_labels == landmark)
            local_position = positions[np.argmin(squared_distance[positions])]
            local_distance = float(squared_distance[local_position])
            if local_distance < best_squared_distance[landmark]:
                best_squared_distance[landmark] = local_distance
                best_indices[landmark] = start + local_position

    # A MiniBatchKMeans center can rarely receive no final labels. Find its
    # nearest cell directly so every fitted centroid remains auditable.
    missing = np.flatnonzero(best_indices < 0)
    if len(missing) > 0:
        missing_centroids = centroids[missing]
        missing_norms = np.square(missing_centroids).sum(axis=1)[None, :]
        for start in range(0, latent.shape[0], chunk_size):
            stop = min(start + chunk_size, latent.shape[0])
            block = np.asarray(latent[start:stop], dtype=np.float32)
            block_norms = np.square(block).sum(axis=1, keepdims=True)
            distances = block @ missing_centroids.T
            distances *= -2.0
            distances += block_norms
            distances += missing_norms
            np.maximum(distances, 0.0, out=distances)
            local_indices = np.argmin(distances, axis=0)
            columns = np.arange(len(missing))
            local_distances = distances[local_indices, columns]
            improved = local_distances < best_squared_distance[missing]
            improved_landmarks = missing[improved]
            best_squared_distance[improved_landmarks] = local_distances[improved]
            best_indices[improved_landmarks] = start + local_indices[improved]

    if np.any(best_indices < 0) or not np.isfinite(best_squared_distance).all():
        raise ValueError("Failed to identify a nearest cell for every landmark")
    return best_indices, np.sqrt(best_squared_distance)


def evaluate_fit(
    atlas: str,
    latent: np.ndarray,
    cell_type_codes: np.ndarray,
    cell_type_names: np.ndarray,
    atlas_counts: np.ndarray,
    n_landmarks: int,
    random_seed: int,
    args: argparse.Namespace,
) -> pd.DataFrame:
    init_size = max(3 * args.kmeans_batch_size, 3 * n_landmarks)
    model = MiniBatchKMeans(
        n_clusters=n_landmarks,
        init="k-means++",
        max_iter=args.kmeans_max_iter,
        batch_size=args.kmeans_batch_size,
        verbose=0,
        compute_labels=True,
        random_state=random_seed,
        tol=0.0,
        max_no_improvement=args.kmeans_max_no_improvement,
        init_size=init_size,
        n_init=args.kmeans_n_init,
        reassignment_ratio=args.kmeans_reassignment_ratio,
    )
    print(
        f"Fitting {atlas}: landmarks={n_landmarks:,}, seed={random_seed}",
        flush=True,
    )
    start_time = time.time()
    model.fit(latent)
    fit_seconds = time.time() - start_time
    labels = np.asarray(model.labels_, dtype=np.int32)
    centroids = np.asarray(model.cluster_centers_, dtype=np.float32)
    if len(labels) != latent.shape[0]:
        raise ValueError("MiniBatchKMeans did not return one label per cell")

    anchor_indices, anchor_distances = nearest_anchor_cells(
        latent,
        labels,
        centroids,
        args.distance_chunk_size,
    )
    anchor_codes = cell_type_codes[anchor_indices]
    n_types = len(cell_type_names)
    cluster_type_counts = np.bincount(
        labels.astype(np.int64) * n_types + cell_type_codes,
        minlength=n_landmarks * n_types,
    ).reshape(n_landmarks, n_types)
    cluster_sizes = cluster_type_counts.sum(axis=1)
    nonempty_clusters = cluster_sizes > 0
    dominant_codes = np.argmax(
        cluster_type_counts[nonempty_clusters],
        axis=1,
    )
    cluster_purity = np.zeros_like(cluster_type_counts, dtype=np.float64)
    cluster_purity[nonempty_clusters] = (
        cluster_type_counts[nonempty_clusters]
        / cluster_sizes[nonempty_clusters, None]
    )
    atlas_fraction = atlas_counts / atlas_counts.sum()
    cluster_enrichment = cluster_purity / atlas_fraction[None, :]
    enough_cells = cluster_type_counts >= args.min_type_cells_enriched
    enriched = enough_cells & (
        cluster_enrichment >= args.min_cluster_enrichment
    )

    n_anchor = np.bincount(anchor_codes, minlength=n_types)
    n_dominant = np.bincount(dominant_codes, minlength=n_types)
    n_enriched = enriched.sum(axis=0)
    max_purity = cluster_purity.max(axis=0)
    max_enrichment = cluster_enrichment.max(axis=0)
    max_type_cells = cluster_type_counts.max(axis=0)

    return pd.DataFrame(
        {
            "atlas": atlas,
            "landmark_count": n_landmarks,
            "random_seed": random_seed,
            "cell_type": cell_type_names,
            "atlas_cell_count": atlas_counts,
            "atlas_cell_fraction": atlas_fraction,
            "n_anchor_landmarks": n_anchor,
            "n_dominant_landmarks": n_dominant,
            "n_enriched_landmarks": n_enriched,
            "anchor_captured": n_anchor > 0,
            "dominant_captured": n_dominant > 0,
            "enriched_captured": n_enriched > 0,
            "maximum_cluster_purity": max_purity,
            "maximum_cluster_enrichment": max_enrichment,
            "maximum_type_cells_in_cluster": max_type_cells,
            "median_anchor_distance": float(np.median(anchor_distances)),
            "maximum_anchor_distance": float(np.max(anchor_distances)),
            "kmeans_inertia": float(model.inertia_),
            "kmeans_iterations": int(getattr(model, "n_iter_", -1)),
            "kmeans_steps": int(getattr(model, "n_steps_", -1)),
            "fit_seconds": fit_seconds,
        }
    )


def minimum_stable_count(
    landmark_counts: list[int],
    capture_rates: np.ndarray,
    threshold: float,
) -> float:
    meets_threshold = capture_rates >= threshold
    for index, landmark_count in enumerate(landmark_counts):
        if np.all(meets_threshold[index:]):
            return float(landmark_count)
    return np.nan


def aggregate_cell_types(
    run_details: pd.DataFrame,
    landmark_counts: list[int],
    required_seed_fraction: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = (
        run_details.groupby(["cell_type", "landmark_count"], as_index=False)
        .agg(
            atlas_cell_count=("atlas_cell_count", "first"),
            atlas_cell_fraction=("atlas_cell_fraction", "first"),
            anchor_capture_rate=("anchor_captured", "mean"),
            dominant_capture_rate=("dominant_captured", "mean"),
            enriched_capture_rate=("enriched_captured", "mean"),
            median_anchor_landmarks=("n_anchor_landmarks", "median"),
            median_dominant_landmarks=("n_dominant_landmarks", "median"),
            median_enriched_landmarks=("n_enriched_landmarks", "median"),
            median_maximum_cluster_purity=("maximum_cluster_purity", "median"),
            median_maximum_cluster_enrichment=(
                "maximum_cluster_enrichment",
                "median",
            ),
        )
    )

    rows = []
    for cell_type, cell_data in grouped.groupby("cell_type", sort=False):
        aligned = cell_data.set_index("landmark_count").reindex(landmark_counts)
        if aligned.isna().any().any():
            raise ValueError(f"Incomplete landmark grid for cell type: {cell_type}")
        row = {
            "cell_type": cell_type,
            "atlas_cell_count": int(aligned["atlas_cell_count"].iloc[0]),
            "atlas_cell_fraction": float(aligned["atlas_cell_fraction"].iloc[0]),
        }
        for criterion in ("anchor", "dominant", "enriched"):
            rates = aligned[f"{criterion}_capture_rate"].to_numpy(dtype=float)
            row[f"minimum_stable_{criterion}_landmarks"] = minimum_stable_count(
                landmark_counts,
                rates,
                required_seed_fraction,
            )
            for landmark_count, rate in zip(landmark_counts, rates):
                row[f"{criterion}_capture_rate_{landmark_count}"] = rate
        rows.append(row)
    cell_summary = pd.DataFrame(rows)
    cell_summary = cell_summary.sort_values(
        ["minimum_stable_anchor_landmarks", "atlas_cell_count", "cell_type"],
        ascending=[False, True, True],
        na_position="first",
    ).reset_index(drop=True)
    return cell_summary, grouped


def build_atlas_summary(
    atlas: str,
    cell_summary: pd.DataFrame,
    run_details: pd.DataFrame,
    landmark_counts: list[int],
    seeds: list[int],
    args: argparse.Namespace,
) -> pd.DataFrame:
    row = {
        "atlas": atlas,
        "n_cells": int(
            run_details.groupby("cell_type")["atlas_cell_count"].first().sum()
        ),
        "n_cell_types": len(cell_summary),
        "landmark_counts_tested": ",".join(map(str, landmark_counts)),
        "random_seeds": ",".join(map(str, seeds)),
        "required_seed_fraction": args.required_seed_fraction,
        "min_type_cells_enriched": args.min_type_cells_enriched,
        "min_cluster_enrichment": args.min_cluster_enrichment,
        "annotation_guided_evaluation": True,
    }
    for criterion in ("anchor", "dominant", "enriched"):
        column = f"minimum_stable_{criterion}_landmarks"
        missing = cell_summary.loc[cell_summary[column].isna(), "cell_type"].tolist()
        row[f"all_types_stably_{criterion}_captured"] = len(missing) == 0
        row[f"cell_types_without_stable_{criterion}_capture"] = ";".join(missing)
        row[f"atlas_minimum_stable_{criterion}_landmarks"] = (
            int(cell_summary[column].max()) if not missing else np.nan
        )
    return pd.DataFrame([row])


def plot_heatmap(
    atlas: str,
    cell_summary: pd.DataFrame,
    landmark_counts: list[int],
    output_path: Path,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rate_columns = [f"anchor_capture_rate_{count}" for count in landmark_counts]
    values = cell_summary[rate_columns].to_numpy(dtype=float)
    labels = [
        f"{row.cell_type} (n={row.atlas_cell_count:,})"
        for row in cell_summary.itertuples(index=False)
    ]
    height = max(8.5, 0.23 * len(cell_summary) + 2.5)
    figure, axis = plt.subplots(figsize=(10.5, height), dpi=dpi)
    image = axis.imshow(
        values,
        aspect="auto",
        interpolation="nearest",
        cmap="YlGnBu",
        vmin=0,
        vmax=1,
        rasterized=True,
    )
    axis.set_xticks(np.arange(len(landmark_counts)))
    axis.set_xticklabels([f"{count:,}" for count in landmark_counts], rotation=45)
    axis.set_yticks(np.arange(len(labels)))
    axis.set_yticklabels(labels, fontsize=7.5)
    axis.set_xlabel("Number of MiniBatchKMeans landmarks")
    axis.set_ylabel("Cell type")
    axis.set_title(
        f"{atlas}: reproducible landmark anchor coverage",
        fontsize=17,
        fontweight="bold",
        pad=14,
    )
    axis.tick_params(length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    colorbar = figure.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
    colorbar.set_label("Fraction of random seeds with >=1 anchor landmark")
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def plot_curves(
    atlas: str,
    run_details: pd.DataFrame,
    landmark_counts: list[int],
    output_path: Path,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {
        "anchor": "#126782",
        "dominant": "#E07A5F",
        "enriched": "#5B8E7D",
    }
    labels = {
        "anchor": "Nearest-cell anchor",
        "dominant": "Dominant cluster type",
        "enriched": "Enriched cluster",
    }
    per_seed = (
        run_details.groupby(["landmark_count", "random_seed"], as_index=False)
        .agg(
            anchor=("anchor_captured", "mean"),
            dominant=("dominant_captured", "mean"),
            enriched=("enriched_captured", "mean"),
        )
    )

    figure, axis = plt.subplots(figsize=(10, 6.5), dpi=dpi)
    for criterion in ("anchor", "dominant", "enriched"):
        aggregate = per_seed.groupby("landmark_count")[criterion].agg(
            ["mean", "min", "max"]
        ).reindex(landmark_counts)
        x = np.asarray(landmark_counts)
        axis.plot(
            x,
            100 * aggregate["mean"],
            marker="o",
            linewidth=2.2,
            color=colors[criterion],
            label=labels[criterion],
        )
        axis.fill_between(
            x,
            100 * aggregate["min"],
            100 * aggregate["max"],
            color=colors[criterion],
            alpha=0.14,
            linewidth=0,
        )
    axis.set_xscale("log", base=2)
    axis.set_xticks(landmark_counts)
    axis.set_xticklabels([f"{count:,}" for count in landmark_counts], rotation=35)
    axis.set_ylim(0, 102)
    axis.set_xlabel("Number of MiniBatchKMeans landmarks")
    axis.set_ylabel("Annotated cell types captured (%)")
    axis.set_title(
        f"{atlas}: cell-type landmark coverage versus landmark count",
        fontsize=16,
        fontweight="bold",
    )
    axis.grid(alpha=0.25, linewidth=0.6)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.legend(frameon=False, loc="lower right")
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    landmark_counts = parse_integer_list(args.landmark_counts, "landmark_counts")
    seeds = parse_integer_list(args.random_seeds, "random_seeds")
    if not 0 < args.required_seed_fraction <= 1:
        raise ValueError("--required_seed_fraction must be in (0, 1]")
    if args.min_type_cells_enriched < 1:
        raise ValueError("--min_type_cells_enriched must be positive")
    if args.min_cluster_enrichment <= 1:
        raise ValueError("--min_cluster_enrichment must be greater than 1")
    if args.distance_chunk_size < 1:
        raise ValueError("--distance_chunk_size must be positive")
    if args.kmeans_batch_size < max(landmark_counts):
        raise ValueError(
            "--kmeans_batch_size must be at least the largest landmark count"
        )

    cell_names = read_cell_names(Path(args.cell_names))
    latent = np.load(Path(args.latent_npy), mmap_mode="r")
    if latent.ndim != 2:
        raise ValueError("PeakVI latent array must be two-dimensional")
    if latent.shape[0] != len(cell_names):
        raise ValueError(
            f"Latent rows ({latent.shape[0]:,}) do not match cell names "
            f"({len(cell_names):,})"
        )
    if latent.shape[0] < max(landmark_counts):
        raise ValueError("Atlas has fewer cells than the largest landmark count")
    if latent.dtype != np.float32:
        print(
            f"WARNING: latent dtype is {latent.dtype}; K-means expects float32",
            flush=True,
        )

    cell_types = read_aligned_cell_types(
        Path(args.cell_metadata),
        cell_names,
        args.cell_type_column,
    )
    cell_type_names = np.array(sorted(cell_types.unique()), dtype=object)
    type_lookup = {name: index for index, name in enumerate(cell_type_names)}
    cell_type_codes = cell_types.map(type_lookup).to_numpy(dtype=np.int32)
    atlas_counts = np.bincount(
        cell_type_codes,
        minlength=len(cell_type_names),
    ).astype(np.int64)

    run_frames = []
    for landmark_count in landmark_counts:
        for seed in seeds:
            run_frames.append(
                evaluate_fit(
                    args.atlas,
                    latent,
                    cell_type_codes,
                    cell_type_names,
                    atlas_counts,
                    landmark_count,
                    seed,
                    args,
                )
            )
    run_details = pd.concat(run_frames, ignore_index=True)
    cell_summary, grouped = aggregate_cell_types(
        run_details,
        landmark_counts,
        args.required_seed_fraction,
    )
    atlas_summary = build_atlas_summary(
        args.atlas,
        cell_summary,
        run_details,
        landmark_counts,
        seeds,
        args,
    )

    output_run_details = Path(args.output_run_details_tsv)
    output_cell_types = Path(args.output_cell_type_tsv)
    output_atlas = Path(args.output_atlas_tsv)
    for path in (output_run_details, output_cell_types, output_atlas):
        path.parent.mkdir(parents=True, exist_ok=True)
    run_details.to_csv(output_run_details, sep="\t", index=False)
    cell_summary.to_csv(output_cell_types, sep="\t", index=False)
    atlas_summary.to_csv(output_atlas, sep="\t", index=False)
    plot_heatmap(
        args.atlas,
        cell_summary,
        landmark_counts,
        Path(args.output_heatmap_png),
        args.dpi,
    )
    plot_curves(
        args.atlas,
        run_details,
        landmark_counts,
        Path(args.output_curve_png),
        args.dpi,
    )

    print(f"Wrote run details: {output_run_details}", flush=True)
    print(f"Wrote cell-type minimums: {output_cell_types}", flush=True)
    print(f"Wrote atlas recommendation: {output_atlas}", flush=True)
    print(f"Wrote capture heatmap: {args.output_heatmap_png}", flush=True)
    print(f"Wrote capture curves: {args.output_curve_png}", flush=True)
    print(atlas_summary.to_string(index=False), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
