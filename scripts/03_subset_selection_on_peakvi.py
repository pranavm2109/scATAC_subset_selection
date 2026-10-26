#!/usr/bin/env python3
"""Build PeakVI landmark features and select one Apricot subset.

Production divergence:
    PeakVI 60D Gaussian latent
    -> 1024 annotation-blind MiniBatchKMeans landmarks
    -> RBF similarities
    -> row-wise L1 normalization to the atlas median RBF mass
    -> Apricot FeatureBasedSelection(sqrt)

No cell-type annotations are loaded or used by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from apricot import FeatureBasedSelection
from sklearn.cluster import MiniBatchKMeans


FEATURE_METHOD = "peakvi_rbf_landmarks"
SELECTION_FUNCTION = "sqrt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build PeakVI RBF landmark features or select one subset."
    )
    parser.add_argument("--mode", choices=["build_features", "select"], required=True)
    parser.add_argument("--atlas", required=True)

    parser.add_argument("--latent_npy")
    parser.add_argument("--cell_names")
    parser.add_argument("--peakvi_run_summary")
    parser.add_argument("--feature_matrix", required=True)
    parser.add_argument("--centroids", required=True)
    parser.add_argument("--landmark_metadata", required=True)
    parser.add_argument("--sigma_tsv", required=True)

    parser.add_argument("--fraction", type=float)
    parser.add_argument("--outdir")

    parser.add_argument(
        "--latent_distribution",
        choices=["ln", "normal"],
        default="normal",
    )
    parser.add_argument("--n_landmarks", type=int, default=1024)
    parser.add_argument("--random_seed", type=int, default=123)
    parser.add_argument("--kmeans_batch_size", type=int, default=8192)
    parser.add_argument("--kmeans_max_iter", type=int, default=100)
    parser.add_argument("--kmeans_n_init", type=int, default=3)
    parser.add_argument("--kmeans_max_no_improvement", type=int, default=20)
    parser.add_argument("--kmeans_reassignment_ratio", type=float, default=0.01)
    parser.add_argument("--distance_chunk_size", type=int, default=50000)
    parser.add_argument("--nonnegative_tolerance", type=float, default=1e-7)
    parser.add_argument("--row_sum_tolerance", type=float, default=1e-4)
    parser.add_argument(
        "--selection_function",
        choices=[SELECTION_FUNCTION],
        default=SELECTION_FUNCTION,
    )
    return parser.parse_args()


def require_arguments(args: argparse.Namespace, names: tuple[str, ...]) -> None:
    missing = [name for name in names if getattr(args, name) in (None, "")]
    if missing:
        flags = ", ".join(f"--{name}" for name in missing)
        raise ValueError(f"{args.mode} requires: {flags}")


def read_cell_names(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any():
        raise ValueError(f"Cell-name file is empty or contains missing values: {path}")
    cell_names = values.to_numpy(dtype=str)
    if len(np.unique(cell_names)) != len(cell_names):
        raise ValueError(f"Cell names are not unique: {path}")
    return cell_names


def cell_name_checksum(cell_names: np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in cell_names:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def fraction_label(fraction: float) -> str:
    if fraction <= 0 or fraction > 1:
        raise ValueError("--fraction must be greater than 0 and no more than 1")
    percent = fraction * 100
    if not np.isclose(percent, round(percent)):
        raise ValueError("--fraction must correspond to a whole percentage")
    return f"{int(round(percent)):02d}pct"


def validate_latent(
    latent: np.ndarray,
    latent_distribution: str,
    nonnegative_tolerance: float,
    row_sum_tolerance: float,
    chunk_size: int,
) -> dict[str, float | int | bool | str]:
    if latent.ndim != 2:
        raise ValueError(f"PeakVI latent array must be 2D; found {latent.shape}")
    if not np.issubdtype(latent.dtype, np.number):
        raise ValueError(f"PeakVI latent array is not numeric: {latent.dtype}")

    n_cells, n_latent = latent.shape
    global_min = math.inf
    global_max = -math.inf
    negative_count = 0
    total_values = 0
    row_sums = np.empty(n_cells, dtype=np.float64)

    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        block = np.asarray(latent[start:stop])
        if not np.isfinite(block).all():
            raise ValueError("PeakVI latent array contains non-finite values")
        global_min = min(global_min, float(block.min()))
        global_max = max(global_max, float(block.max()))
        negative_count += int((block < 0).sum())
        total_values += block.size
        row_sums[start:stop] = block.sum(axis=1, dtype=np.float64)

    fraction_negative = negative_count / total_values
    row_sum_close = np.isclose(
        row_sums,
        1.0,
        atol=row_sum_tolerance,
        rtol=row_sum_tolerance,
    )

    if latent_distribution == "ln":
        if global_min < -nonnegative_tolerance:
            raise ValueError(
                "Logistic-normal PeakVI latent values are not nonnegative "
                f"within tolerance: minimum={global_min:.9g}"
            )
        if not row_sum_close.all():
            raise ValueError(
                "Logistic-normal PeakVI latent rows are not approximately "
                f"simplex-valued: {int((~row_sum_close).sum()):,}/{n_cells:,} "
                "rows fall outside tolerance"
            )

    summary = {
        "n_cells": n_cells,
        "n_latent": n_latent,
        "latent_distribution": latent_distribution,
        "raw_latent_min": global_min,
        "raw_latent_max": global_max,
        "fraction_negative": fraction_negative,
        "row_sum_min": float(row_sums.min()),
        "row_sum_median": float(np.median(row_sums)),
        "row_sum_mean": float(row_sums.mean()),
        "row_sum_max": float(row_sums.max()),
        "fraction_row_sums_near_one": float(row_sum_close.mean()),
        "approximately_simplex_valued": bool(row_sum_close.all()),
    }
    print(f"Raw latent shape: {latent.shape}", flush=True)
    print(f"PeakVI latent distribution: {latent_distribution}", flush=True)
    print(f"Raw latent min/max: {global_min:.9g} / {global_max:.9g}", flush=True)
    print(f"Fraction of raw latent values < 0: {fraction_negative:.9g}", flush=True)
    print(
        "Raw latent row sums, min/median/mean/max: "
        f"{summary['row_sum_min']:.9g} / {summary['row_sum_median']:.9g} / "
        f"{summary['row_sum_mean']:.9g} / {summary['row_sum_max']:.9g}",
        flush=True,
    )
    print("No min-shift or other transformation will be applied.", flush=True)
    return summary

def squared_distances(
    points: np.ndarray,
    centroids: np.ndarray,
    centroid_norms: np.ndarray,
) -> np.ndarray:
    point_norms = np.square(points).sum(axis=1, keepdims=True)
    distances = points @ centroids.T
    distances *= -2.0
    distances += point_norms
    distances += centroid_norms
    np.maximum(distances, 0.0, out=distances)
    return distances


def build_landmark_features(args: argparse.Namespace) -> None:
    require_arguments(args, ("latent_npy", "cell_names", "peakvi_run_summary"))
    if args.n_landmarks < 2:
        raise ValueError("--n_landmarks must be at least 2")
    if args.kmeans_batch_size < args.n_landmarks:
        raise ValueError("--kmeans_batch_size must be at least --n_landmarks")
    if args.distance_chunk_size < 1:
        raise ValueError("--distance_chunk_size must be positive")

    latent_path = Path(args.latent_npy)
    cell_names_path = Path(args.cell_names)
    peakvi_summary_path = Path(args.peakvi_run_summary)
    feature_path = Path(args.feature_matrix)
    centroids_path = Path(args.centroids)
    metadata_path = Path(args.landmark_metadata)
    sigma_path = Path(args.sigma_tsv)
    for path in (feature_path, centroids_path, metadata_path, sigma_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    with open(peakvi_summary_path, encoding="utf-8") as handle:
        peakvi_summary = json.load(handle)
    saved_distribution = peakvi_summary.get("latent_distribution")
    if saved_distribution != args.latent_distribution:
        raise ValueError(
            "PeakVI distribution mismatch: "
            f"expected {args.latent_distribution}, summary reports "
            f"{saved_distribution}"
        )

    latent = np.load(latent_path, mmap_mode="r")
    cell_names = read_cell_names(cell_names_path)
    if len(cell_names) != latent.shape[0]:
        raise ValueError(
            f"Cell names ({len(cell_names):,}) and latent rows "
            f"({latent.shape[0]:,}) differ"
        )
    if latent.shape[0] < args.n_landmarks:
        raise ValueError("The atlas has fewer cells than requested landmarks")

    latent_summary = validate_latent(
        latent,
        args.latent_distribution,
        args.nonnegative_tolerance,
        args.row_sum_tolerance,
        args.distance_chunk_size,
    )
    if latent.dtype != np.float32:
        print(
            f"WARNING: canonical latent dtype is {latent.dtype}; "
            "MiniBatchKMeans will operate in float32.",
            flush=True,
        )
        latent_for_fit = np.asarray(latent, dtype=np.float32)
    else:
        latent_for_fit = latent

    init_size = max(3 * args.kmeans_batch_size, 3 * args.n_landmarks)
    print(
        "Fitting annotation-blind MiniBatchKMeans landmarks: "
        f"n_clusters={args.n_landmarks}, batch_size={args.kmeans_batch_size}, "
        f"max_iter={args.kmeans_max_iter}, n_init={args.kmeans_n_init}, "
        f"random_state={args.random_seed}",
        flush=True,
    )
    fit_start = time.time()
    model = MiniBatchKMeans(
        n_clusters=args.n_landmarks,
        init="k-means++",
        max_iter=args.kmeans_max_iter,
        batch_size=args.kmeans_batch_size,
        verbose=0,
        compute_labels=True,
        random_state=args.random_seed,
        tol=0.0,
        max_no_improvement=args.kmeans_max_no_improvement,
        init_size=init_size,
        n_init=args.kmeans_n_init,
        reassignment_ratio=args.kmeans_reassignment_ratio,
    )
    model.fit(latent_for_fit)
    fit_seconds = time.time() - fit_start
    centroids = np.asarray(model.cluster_centers_, dtype=np.float32, order="C")
    del latent_for_fit
    np.save(centroids_path, centroids)
    print(
        "MiniBatchKMeans completed: "
        f"n_iter={getattr(model, 'n_iter_', 'NA')}, "
        f"n_steps={getattr(model, 'n_steps_', 'NA')}, "
        f"inertia={float(model.inertia_):.9g}, seconds={fit_seconds:.2f}",
        flush=True,
    )

    centroid_norms = np.square(centroids).sum(axis=1)[None, :]
    nearest_distances = np.empty(latent.shape[0], dtype=np.float32)
    print("Pass 1/2: calculating nearest-landmark distances", flush=True)
    for start in range(0, latent.shape[0], args.distance_chunk_size):
        stop = min(start + args.distance_chunk_size, latent.shape[0])
        points = np.asarray(latent[start:stop], dtype=np.float32)
        distances = squared_distances(points, centroids, centroid_norms)
        nearest_distances[start:stop] = np.sqrt(distances.min(axis=1))
        del points, distances

    sigma = float(np.median(nearest_distances))
    if not np.isfinite(sigma) or sigma <= np.finfo(np.float32).eps:
        raise ValueError(
            "The median nearest-landmark distance produced an invalid RBF "
            f"bandwidth: sigma={sigma}"
        )
    distance_quantiles = np.quantile(nearest_distances, [0.0, 0.5, 0.9, 0.99, 1.0])
    print(
        "Nearest-landmark distance min/median/p90/p99/max: "
        + " / ".join(f"{value:.9g}" for value in distance_quantiles),
        flush=True,
    )
    print(f"RBF sigma (median nearest-landmark distance): {sigma:.9g}", flush=True)

    shape = (latent.shape[0], args.n_landmarks)
    feature_bytes = int(np.prod(shape, dtype=np.int64) * np.dtype(np.float32).itemsize)
    print(
        f"RBF feature matrix shape: {shape}; dtype=float32; "
        f"size={feature_bytes / 1024**3:.3f} GiB",
        flush=True,
    )
    features = np.lib.format.open_memmap(
        feature_path,
        mode="w+",
        dtype=np.float32,
        shape=shape,
    )
    scale = np.float32(-1.0 / (2.0 * sigma * sigma))
    feature_min = math.inf
    feature_max = -math.inf
    print("Pass 2/2: writing RBF similarities directly to final memmap", flush=True)
    for start in range(0, latent.shape[0], args.distance_chunk_size):
        stop = min(start + args.distance_chunk_size, latent.shape[0])
        points = np.asarray(latent[start:stop], dtype=np.float32)
        distances = squared_distances(points, centroids, centroid_norms)
        distances *= scale
        np.exp(distances, out=distances)
        feature_min = min(feature_min, float(distances.min()))
        feature_max = max(feature_max, float(distances.max()))
        features[start:stop] = distances
        del points, distances
    features.flush()
    del features

    if not np.isfinite(feature_min) or not np.isfinite(feature_max):
        raise ValueError("RBF feature matrix contains non-finite values")
    if feature_min < 0 or feature_max > 1.0 + 1e-6:
        raise ValueError(
            f"RBF feature range is invalid: {feature_min} to {feature_max}"
        )
    print(f"RBF feature min/max: {feature_min:.9g} / {feature_max:.9g}", flush=True)

    checksum = cell_name_checksum(cell_names)
    metadata = {
        "atlas": args.atlas,
        "feature_construction_method": FEATURE_METHOD,
        "canonical_peakvi_latent": str(latent_path),
        "peakvi_run_summary": str(peakvi_summary_path),
        "peakvi_latent_distribution": args.latent_distribution,
        "cell_names_file": str(cell_names_path),
        "cell_name_sha256": checksum,
        "n_cells": int(latent.shape[0]),
        "n_peakvi_dimensions": int(latent.shape[1]),
        "n_landmarks": int(args.n_landmarks),
        "landmark_centroids": str(centroids_path),
        "rbf_feature_matrix": str(feature_path),
        "feature_matrix_dtype": "float32",
        "feature_matrix_shape": list(shape),
        "feature_matrix_bytes": feature_bytes,
        "feature_matrix_gib": feature_bytes / 1024**3,
        "rbf_sigma": sigma,
        "sigma_definition": "median Euclidean distance to nearest landmark",
        "min_shift_applied": False,
        "nonnegative_transform": "none",
        "selection_class": "apricot.FeatureBasedSelection",
        "selection_function": SELECTION_FUNCTION,
        "random_seed": int(args.random_seed),
        "minibatch_kmeans": {
            "n_clusters": int(args.n_landmarks),
            "init": "k-means++",
            "batch_size": int(args.kmeans_batch_size),
            "max_iter": int(args.kmeans_max_iter),
            "n_init": int(args.kmeans_n_init),
            "max_no_improvement": int(args.kmeans_max_no_improvement),
            "reassignment_ratio": float(args.kmeans_reassignment_ratio),
            "tol": 0.0,
            "init_size": int(init_size),
            "random_state": int(args.random_seed),
            "n_iter": int(getattr(model, "n_iter_", -1)),
            "n_steps": int(getattr(model, "n_steps_", -1)),
            "inertia": float(model.inertia_),
            "fit_seconds": fit_seconds,
        },
        "nearest_landmark_distance": {
            "minimum": float(distance_quantiles[0]),
            "median": float(distance_quantiles[1]),
            "mean": float(nearest_distances.mean()),
            "p90": float(distance_quantiles[2]),
            "p99": float(distance_quantiles[3]),
            "maximum": float(distance_quantiles[4]),
        },
        "raw_latent_validation": latent_summary,
        "rbf_feature_min": feature_min,
        "rbf_feature_max": feature_max,
        "distance_chunk_size": int(args.distance_chunk_size),
        "sklearn_version": sklearn.__version__,
    }
    with open(metadata_path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    pd.DataFrame(
        {
            "atlas": [args.atlas],
            "feature_construction_method": [FEATURE_METHOD],
            "sigma": [sigma],
            "sigma_definition": [metadata["sigma_definition"]],
            "n_landmarks": [args.n_landmarks],
            "n_peakvi_dimensions": [latent.shape[1]],
            "peakvi_latent_distribution": [args.latent_distribution],
            "random_seed": [args.random_seed],
            "min_shift_applied": [False],
        }
    ).to_csv(sigma_path, sep="\t", index=False)
    print(f"Wrote landmark centroids: {centroids_path}", flush=True)
    print(f"Wrote landmark metadata: {metadata_path}", flush=True)
    print(f"Wrote sigma summary: {sigma_path}", flush=True)
    print("DONE", flush=True)


def select_subset(args: argparse.Namespace) -> None:
    require_arguments(args, ("cell_names", "fraction", "outdir"))
    frac_label = fraction_label(args.fraction)
    feature_path = Path(args.feature_matrix)
    centroids_path = Path(args.centroids)
    metadata_path = Path(args.landmark_metadata)
    sigma_path = Path(args.sigma_tsv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(metadata_path, encoding="utf-8") as handle:
        metadata = json.load(handle)
    if metadata.get("feature_construction_method") != FEATURE_METHOD:
        raise ValueError("Landmark metadata does not describe PeakVI RBF features")
    if metadata.get("min_shift_applied") is not False:
        raise ValueError("Landmark metadata unexpectedly reports a min-shift")
    if metadata.get("selection_function") != SELECTION_FUNCTION:
        raise ValueError("Landmark metadata does not specify sqrt selection")
    if str(metadata.get("atlas")) != args.atlas:
        raise ValueError("Atlas mismatch in landmark metadata")

    cell_names = read_cell_names(Path(args.cell_names))
    if cell_name_checksum(cell_names) != metadata.get("cell_name_sha256"):
        raise ValueError("Cell-name checksum does not match landmark metadata")

    features = np.load(feature_path, mmap_mode="r")
    centroids = np.load(centroids_path, mmap_mode="r")
    sigma_frame = pd.read_csv(sigma_path, sep="\t")
    sigma = float(sigma_frame.loc[0, "sigma"])
    expected_shape = tuple(metadata["feature_matrix_shape"])
    if features.shape != expected_shape:
        raise ValueError(
            f"RBF feature shape {features.shape} != metadata {expected_shape}"
        )
    if len(cell_names) != features.shape[0]:
        raise ValueError("Cell names and RBF feature rows differ")
    if centroids.shape != (
        int(metadata["n_landmarks"]),
        int(metadata["n_peakvi_dimensions"]),
    ):
        raise ValueError("Landmark centroid shape does not match metadata")
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError(f"Invalid saved sigma: {sigma}")

    subset_size = max(int(round(features.shape[0] * args.fraction)), 1)
    print(f"Atlas: {args.atlas}", flush=True)
    print(f"Fraction: {args.fraction} ({frac_label})", flush=True)
    print(
        f"Feature construction: PeakVI {metadata['n_peakvi_dimensions']}D "
        f"{metadata['peakvi_latent_distribution']} -> "
        f"{metadata['n_landmarks']} RBF landmark features",
        flush=True,
    )
    print("Min-shift applied: false", flush=True)
    print(f"RBF sigma: {sigma:.9g}", flush=True)
    print(
        f"RBF feature matrix: shape={features.shape}, dtype={features.dtype}, "
        f"min/max={metadata['rbf_feature_min']:.9g}/"
        f"{metadata['rbf_feature_max']:.9g}, "
        f"disk size={metadata['feature_matrix_gib']:.3f} GiB",
        flush=True,
    )
    print(f"Expected selected cells: {subset_size:,}", flush=True)
    row_sums = np.empty(features.shape[0], dtype=np.float64)
    for start in range(0, features.shape[0], args.distance_chunk_size):
        stop = min(start + args.distance_chunk_size, features.shape[0])
        block = np.asarray(features[start:stop], dtype=np.float32)
        if not np.isfinite(block).all() or np.any(block < 0):
            raise ValueError("RBF matrix contains negative or non-finite values")
        row_sums[start:stop] = block.sum(axis=1, dtype=np.float64)
    if np.any(row_sums <= 0):
        raise ValueError("L1 normalization requires positive RBF mass in every row")

    row_sum_target = float(np.median(row_sums))
    row_scales = row_sum_target / row_sums
    normalized_features = np.empty(features.shape, dtype=np.float32)
    for start in range(0, features.shape[0], args.distance_chunk_size):
        stop = min(start + args.distance_chunk_size, features.shape[0])
        np.multiply(
            features[start:stop],
            row_scales[start:stop, None],
            out=normalized_features[start:stop],
            casting="unsafe",
        )

    print(
        "Starting apricot.FeatureBasedSelection on L1-normalized, unweighted "
        "RBF features with sqrt concave function",
        flush=True,
    )
    print(
        f"RBF row-sum target={row_sum_target:.9g}; "
        f"scale min/median/max={row_scales.min():.6g}/"
        f"{np.median(row_scales):.6g}/{row_scales.max():.6g}",
        flush=True,
    )

    start_time = time.time()
    selector = FeatureBasedSelection(
        subset_size,
        SELECTION_FUNCTION,
        random_state=args.random_seed,
        verbose=True,
    )
    selector.fit(normalized_features)
    elapsed = time.time() - start_time
    del normalized_features

    ranking = np.asarray(selector.ranking, dtype=np.int64)
    if len(ranking) != subset_size:
        raise RuntimeError(
            f"Expected {subset_size:,} selected indices, got {len(ranking):,}"
        )
    if np.any(ranking < 0) or np.any(ranking >= len(cell_names)):
        raise RuntimeError("Apricot returned an out-of-range cell index")
    if len(np.unique(ranking)) != subset_size:
        raise RuntimeError("Apricot returned duplicate selected indices")
    selected_cells = cell_names[ranking]
    if len(np.unique(selected_cells)) != subset_size:
        raise RuntimeError("Selected cell names are not unique")

    selected_cells_file = outdir / f"selected_cells_peakvi_{frac_label}.txt"
    summary_file = outdir / f"subset_selection_summary_peakvi_{frac_label}.tsv"
    pd.Series(selected_cells).to_csv(
        selected_cells_file,
        index=False,
        header=False,
    )

    gains = np.asarray(selector.gains, dtype=float)
    summary = pd.DataFrame(
        {
            "atlas": [args.atlas],
            "fraction": [args.fraction],
            "fraction_label": [frac_label],
            "subset_size": [subset_size],
            "n_selected_unique": [len(np.unique(selected_cells))],
            "time_seconds": [elapsed],
            "n_cells_total": [features.shape[0]],
            "n_latent": [metadata["n_peakvi_dimensions"]],
            "peakvi_latent_distribution": [
                metadata["peakvi_latent_distribution"]
            ],
            "feature_construction_method": [FEATURE_METHOD],
            "selection_feature_method": ["rbf_landmarks_l1_normalized"],
            "rbf_row_normalization": ["l1_to_atlas_median_row_sum"],
            "rbf_row_sum_target": [row_sum_target],
            "rbf_row_sum_min": [float(row_sums.min())],
            "rbf_row_sum_median": [float(np.median(row_sums))],
            "rbf_row_sum_max": [float(row_sums.max())],
            "rbf_row_scale_min": [float(row_scales.min())],
            "rbf_row_scale_median": [float(np.median(row_scales))],
            "rbf_row_scale_max": [float(row_scales.max())],
            "n_landmark_features": [metadata["n_landmarks"]],
            "n_selection_features": [int(metadata["n_landmarks"])],
            "landmark_centroids": [str(centroids_path)],
            "landmark_metadata": [str(metadata_path)],
            "rbf_feature_matrix": [str(feature_path)],
            "rbf_sigma": [sigma],
            "sigma_definition": [metadata["sigma_definition"]],
            "nonnegative_transform": ["none"],
            "min_shift_applied": [False],
            "selection_class": ["apricot.FeatureBasedSelection"],
            "selection_function": [SELECTION_FUNCTION],
            "selection_optimizer": [selector.optimizer],
            "density_weighting": [False],
            "density_knn_k": [np.nan],
            "density_weight_formula": ["none"],
            "density_max_weight": [1.0],
            "density_sparsity_file": [""],
            "selected_density_weight_median": [np.nan],
            "selected_density_weight_q90": [np.nan],
            "selected_mean_knn_distance_median": [np.nan],
            "random_seed_landmarks": [metadata["random_seed"]],
            "random_seed_selection": [args.random_seed],
            "selected_gain_min": [float(gains.min())],
            "selected_gain_median": [float(np.median(gains))],
            "selected_gain_mean": [float(gains.mean())],
            "selected_gain_max": [float(gains.max())],
            "selected_cells_file": [str(selected_cells_file)],
        }
    )
    summary.to_csv(summary_file, sep="\t", index=False)

    print(
        f"Completed {frac_label}: selected {len(selected_cells):,} unique cells "
        f"in {elapsed:.2f} seconds",
        flush=True,
    )
    print(f"Selected cells: {selected_cells_file}", flush=True)
    print(f"Selection summary: {summary_file}", flush=True)
    print("DONE", flush=True)


def main() -> None:
    args = parse_args()
    if args.mode == "build_features":
        build_landmark_features(args)
    else:
        select_subset(args)


if __name__ == "__main__":
    main()
