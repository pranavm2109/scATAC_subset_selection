#!/usr/bin/env python3
"""Diagnose whether PeakVI RBF landmarks cover annotated cell populations.

The RBF landmarks are MiniBatchKMeans centroids in PeakVI latent space, not
observed cells. The saved UMAP contains coordinates only for observed cells.
For an auditable visualization, each centroid is therefore represented by its
nearest observed cell. Because RBF similarity decreases monotonically with
distance, the nearest cell is the row with maximum similarity for that
landmark in the saved RBF feature matrix.

Cell-type labels are used only after landmark construction for this diagnostic.
They did not influence PeakVI training, landmark placement, or subset selection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Map PeakVI RBF centroids to their nearest observed cells, plot "
            "the resulting landmark anchors on the existing UMAP, and count "
            "landmarks by annotated cell type."
        )
    )
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--feature_matrix", required=True)
    parser.add_argument("--centroids", required=True)
    parser.add_argument("--landmark_metadata", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--cell_metadata", required=True)
    parser.add_argument("--umap_tsv", required=True)
    parser.add_argument("--output_umap_png", required=True)
    parser.add_argument("--output_cell_type_tsv", required=True)
    parser.add_argument("--output_landmark_tsv", required=True)
    parser.add_argument("--cell_type_column", default="cell_type")
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=25000,
        help="Rows read from the RBF feature matrix per pass [default: 25000].",
    )
    parser.add_argument(
        "--max_background_cells",
        type=int,
        default=500000,
        help=(
            "Maximum number of atlas cells drawn behind landmarks; the table "
            "always uses every cell [default: 500000]."
        ),
    )
    parser.add_argument("--random_seed", type=int, default=123)
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


def cell_name_checksum(cell_names: np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in cell_names:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_aligned_annotations(
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


def read_aligned_umap(path: Path, cell_names: np.ndarray) -> pd.DataFrame:
    # Numeric-looking identifiers must remain strings for canonical alignment.
    umap = pd.read_csv(path, sep="\t", dtype={"cell_name": str})
    required = {"cell_name", "UMAP_1", "UMAP_2"}
    missing_columns = required - set(umap.columns)
    if missing_columns:
        raise ValueError(
            f"UMAP TSV is missing columns {sorted(missing_columns)}: {path}"
        )
    if umap["cell_name"].duplicated().any():
        raise ValueError(f"UMAP cell names are not unique: {path}")

    aligned = umap.set_index("cell_name").reindex(cell_names)
    if aligned[["UMAP_1", "UMAP_2"]].isna().any().any():
        missing = aligned.index[
            aligned[["UMAP_1", "UMAP_2"]].isna().any(axis=1)
        ][:10]
        raise ValueError(
            f"UMAP is missing coordinates for PeakVI cells. "
            f"Examples: {missing.tolist()}"
        )
    aligned = aligned[["UMAP_1", "UMAP_2"]].astype(np.float32)
    if not np.isfinite(aligned.to_numpy()).all():
        raise ValueError("UMAP coordinates contain non-finite values")
    return aligned


def nearest_cells_from_rbf(
    features: np.ndarray,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the cell index and maximum RBF similarity for each landmark."""
    n_cells, n_landmarks = features.shape
    best_indices = np.full(n_landmarks, -1, dtype=np.int64)
    best_similarity = np.full(n_landmarks, -np.inf, dtype=np.float32)

    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        block = np.asarray(features[start:stop], dtype=np.float32)
        local_indices = np.argmax(block, axis=0)
        columns = np.arange(n_landmarks)
        local_similarity = block[local_indices, columns]
        improved = local_similarity > best_similarity
        best_similarity[improved] = local_similarity[improved]
        best_indices[improved] = start + local_indices[improved]
        print(
            f"Scanned RBF rows {start:,}-{stop:,} of {n_cells:,}",
            flush=True,
        )

    if np.any(best_indices < 0) or not np.isfinite(best_similarity).all():
        raise ValueError("Failed to identify a nearest cell for every landmark")
    if np.any(best_similarity <= 0) or np.any(best_similarity > 1 + 1e-6):
        raise ValueError(
            "Nearest-cell RBF similarities fall outside the expected (0, 1] range"
        )
    return best_indices, best_similarity


def build_cell_type_summary(
    atlas: str,
    cell_types: pd.Series,
    landmark_table: pd.DataFrame,
) -> pd.DataFrame:
    atlas_counts = cell_types.value_counts(sort=False).rename("n_cells")
    landmark_counts = (
        landmark_table["nearest_cell_type"]
        .value_counts(sort=False)
        .rename("n_landmarks")
    )
    unique_anchor_counts = (
        landmark_table.groupby("nearest_cell_type")["nearest_cell_id"]
        .nunique()
        .rename("n_unique_anchor_cells")
    )

    summary = pd.concat(
        [atlas_counts, landmark_counts, unique_anchor_counts], axis=1
    ).fillna(0)
    summary[["n_cells", "n_landmarks", "n_unique_anchor_cells"]] = summary[
        ["n_cells", "n_landmarks", "n_unique_anchor_cells"]
    ].astype(int)
    summary.index.name = "cell_type"
    summary = summary.reset_index()
    summary.insert(0, "atlas", atlas)

    n_cells = int(summary["n_cells"].sum())
    n_landmarks = int(summary["n_landmarks"].sum())
    summary["atlas_cell_fraction"] = summary["n_cells"] / n_cells
    summary["landmark_fraction"] = summary["n_landmarks"] / n_landmarks
    summary["expected_landmarks_by_cell_fraction"] = (
        summary["atlas_cell_fraction"] * n_landmarks
    )
    summary["landmark_enrichment_vs_abundance"] = np.where(
        summary["expected_landmarks_by_cell_fraction"] > 0,
        summary["n_landmarks"]
        / summary["expected_landmarks_by_cell_fraction"],
        np.nan,
    )
    summary["has_landmark"] = summary["n_landmarks"] > 0
    return summary.sort_values(
        ["has_landmark", "n_landmarks", "n_cells", "cell_type"],
        ascending=[True, True, True, True],
    ).reset_index(drop=True)


def plot_landmarks(
    atlas: str,
    umap: pd.DataFrame,
    nearest_indices: np.ndarray,
    summary: pd.DataFrame,
    output_path: Path,
    max_background_cells: int,
    random_seed: int,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    coordinates = umap.to_numpy(dtype=np.float32)
    rng = np.random.default_rng(random_seed)
    if len(coordinates) > max_background_cells:
        background_indices = np.sort(
            rng.choice(
                len(coordinates),
                size=max_background_cells,
                replace=False,
            )
        )
        background_note = (
            f"Background: deterministic sample of {max_background_cells:,} "
            f"of {len(coordinates):,} cells"
        )
    else:
        background_indices = np.arange(len(coordinates))
        background_note = f"Background: all {len(coordinates):,} cells"

    anchor_coordinates = coordinates[nearest_indices]
    n_landmarks = len(nearest_indices)
    n_unique_anchors = len(np.unique(nearest_indices))
    n_types_with_landmark = int(summary["has_landmark"].sum())
    n_types = len(summary)

    figure, axis = plt.subplots(figsize=(11, 9), dpi=dpi)
    axis.scatter(
        coordinates[background_indices, 0],
        coordinates[background_indices, 1],
        s=0.35,
        c="#C9D0D5",
        alpha=0.42,
        linewidths=0,
        rasterized=True,
    )
    axis.scatter(
        anchor_coordinates[:, 0],
        anchor_coordinates[:, 1],
        s=17,
        c="black",
        alpha=0.88,
        edgecolors="white",
        linewidths=0.28,
        rasterized=True,
        zorder=3,
    )

    axis.set_title(
        f"{atlas}: PeakVI RBF landmark coverage",
        fontsize=19,
        fontweight="bold",
        pad=24,
    )
    axis.text(
        0.5,
        1.01,
        (
            f"{n_landmarks:,} centroids mapped to {n_unique_anchors:,} nearest "
            f"observed cells; {n_types_with_landmark}/{n_types} cell types "
            "receive at least one landmark"
        ),
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=11,
        color="#30383D",
    )
    axis.set_xlabel("PeakVI UMAP 1", fontsize=12)
    axis.set_ylabel("PeakVI UMAP 2", fontsize=12)
    axis.set_aspect("equal", adjustable="datalim")
    axis.tick_params(labelsize=9, colors="#48545B")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#C9D0D5",
            markeredgecolor="none",
            markersize=7,
            label="Atlas cells",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="black",
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=8,
            label="Nearest observed cell to an RBF centroid",
        ),
    ]
    axis.legend(handles=legend, loc="best", frameon=False, fontsize=10)
    figure.text(
        0.01,
        0.01,
        (
            f"{background_note}. Black points are cell-coordinate proxies for "
            "centroids; centroids themselves were not embedded by the saved UMAP."
        ),
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#4D5960",
    )
    figure.tight_layout(rect=[0, 0.035, 1, 1])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, facecolor="white", bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if args.chunk_size < 1:
        raise ValueError("--chunk_size must be positive")
    if args.max_background_cells < 1:
        raise ValueError("--max_background_cells must be positive")

    feature_path = Path(args.feature_matrix)
    centroids_path = Path(args.centroids)
    metadata_path = Path(args.landmark_metadata)
    cell_names_path = Path(args.cell_names)
    cell_metadata_path = Path(args.cell_metadata)
    umap_path = Path(args.umap_tsv)
    output_umap = Path(args.output_umap_png)
    output_cell_types = Path(args.output_cell_type_tsv)
    output_landmarks = Path(args.output_landmark_tsv)

    with metadata_path.open(encoding="utf-8") as handle:
        landmark_metadata = json.load(handle)
    if landmark_metadata.get("atlas") != args.atlas:
        raise ValueError(
            f"Atlas mismatch: requested {args.atlas}, landmark metadata reports "
            f"{landmark_metadata.get('atlas')}"
        )

    cell_names = read_cell_names(cell_names_path)
    expected_checksum = landmark_metadata.get("cell_name_sha256")
    observed_checksum = cell_name_checksum(cell_names)
    if expected_checksum and observed_checksum != expected_checksum:
        raise ValueError(
            "Cell-name checksum does not match the landmark feature metadata"
        )

    features = np.load(feature_path, mmap_mode="r")
    centroids = np.load(centroids_path, mmap_mode="r")
    if features.ndim != 2 or centroids.ndim != 2:
        raise ValueError("Feature matrix and centroids must both be two-dimensional")
    if features.shape[0] != len(cell_names):
        raise ValueError(
            f"Feature rows ({features.shape[0]:,}) do not match cell names "
            f"({len(cell_names):,})"
        )
    if features.shape[1] != centroids.shape[0]:
        raise ValueError(
            f"Feature columns ({features.shape[1]:,}) do not match centroids "
            f"({centroids.shape[0]:,})"
        )

    cell_types = read_aligned_annotations(
        cell_metadata_path,
        cell_names,
        args.cell_type_column,
    )
    umap = read_aligned_umap(umap_path, cell_names)

    print(
        f"Finding nearest observed cells for {features.shape[1]:,} landmarks "
        f"across {features.shape[0]:,} atlas cells",
        flush=True,
    )
    nearest_indices, nearest_similarity = nearest_cells_from_rbf(
        features,
        args.chunk_size,
    )

    sigma = float(landmark_metadata["rbf_sigma"])
    nearest_distance = np.sqrt(
        np.maximum(
            -2.0 * sigma * sigma * np.log(np.clip(nearest_similarity, 1e-30, 1.0)),
            0.0,
        )
    )
    landmark_table = pd.DataFrame(
        {
            "atlas": args.atlas,
            "landmark_id": np.arange(1, features.shape[1] + 1),
            "nearest_cell_index": nearest_indices,
            "nearest_cell_id": cell_names[nearest_indices],
            "nearest_cell_type": cell_types.iloc[nearest_indices].to_numpy(),
            "nearest_rbf_similarity": nearest_similarity,
            "nearest_latent_distance": nearest_distance,
            "UMAP_1": umap.iloc[nearest_indices]["UMAP_1"].to_numpy(),
            "UMAP_2": umap.iloc[nearest_indices]["UMAP_2"].to_numpy(),
        }
    )
    summary = build_cell_type_summary(args.atlas, cell_types, landmark_table)

    output_cell_types.parent.mkdir(parents=True, exist_ok=True)
    output_landmarks.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_cell_types, sep="\t", index=False)
    landmark_table.to_csv(output_landmarks, sep="\t", index=False)
    plot_landmarks(
        atlas=args.atlas,
        umap=umap,
        nearest_indices=nearest_indices,
        summary=summary,
        output_path=output_umap,
        max_background_cells=args.max_background_cells,
        random_seed=args.random_seed,
        dpi=args.dpi,
    )

    missing_types = summary.loc[~summary["has_landmark"], "cell_type"].tolist()
    print(f"Wrote landmark UMAP: {output_umap}", flush=True)
    print(f"Wrote cell-type landmark counts: {output_cell_types}", flush=True)
    print(f"Wrote landmark-level assignments: {output_landmarks}", flush=True)
    print(
        f"Cell types with >=1 landmark: {summary['has_landmark'].sum()}/"
        f"{len(summary)}",
        flush=True,
    )
    if missing_types:
        print("Cell types with zero landmarks: " + "; ".join(missing_types), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
