#!/usr/bin/env python3
"""Retrospectively diagnose cell-type separability in PeakVI latent space.

Annotations are used only after PeakVI, its neighbor graph, and the existing
Apricot selections have been produced. PyNNDescent provides scalable
approximate neighbor search; silhouette values use a bounded stratified sample
to avoid an atlas-wide pairwise distance matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pynndescent import NNDescent
from sklearn.metrics import silhouette_samples


SUBSETS = ("02pct", "04pct", "06pct", "08pct", "10pct")
SUBSET_LABELS = dict(zip(SUBSETS, ("2%", "4%", "6%", "8%", "10%")))
STATUSES = ("absent", "underpowered", "evaluable")
COLORS = {
    "absent": "#C83E4D",
    "underpowered": "#F2A541",
    "evaluable": "#4C78A8",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare PeakVI separability with existing Apricot retention."
    )
    parser.add_argument("--latent_npy", required=True)
    parser.add_argument("--cell_ids", required=True)
    parser.add_argument("--metadata_tsv", required=True)
    parser.add_argument("--selected_cells", nargs="+", required=True)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--output_table", required=True)
    parser.add_argument("--output_summary", required=True)
    parser.add_argument("--output_cell_sparsity_npz", required=True)
    parser.add_argument("--output_figures_dir", required=True)
    parser.add_argument("--knn_k", type=int, default=30)
    parser.add_argument("--silhouette_max_cells", type=int, default=10000)
    parser.add_argument("--high_purity_threshold", type=float, default=0.75)
    parser.add_argument("--max_labels", type=int, default=12)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--random_seed", type=int, default=123)
    return parser.parse_args()


def cell_name_checksum(cell_ids: np.ndarray) -> str:
    digest = hashlib.sha256()
    for cell_id in cell_ids:
        digest.update(str(cell_id).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_ids(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any():
        raise ValueError(f"Empty or missing cell IDs in {path}")
    return values.to_numpy(dtype=str)


def subset_from_filename(path: Path) -> str:
    match = re.search(r"_(02pct|04pct|06pct|08pct|10pct)\.txt$", path.name)
    if match is None:
        raise ValueError(f"Cannot infer subset size from {path}")
    return match.group(1)


def load_inputs(args: argparse.Namespace):
    # PyNNDescent's Numba projection-tree kernels require writable C-contiguous
    # input. At 60 dimensions, even HECA occupies only about 350 MB as float32.
    latent = np.array(np.load(args.latent_npy, mmap_mode="r"), dtype=np.float32, order="C")
    if latent.ndim != 2 or not np.issubdtype(latent.dtype, np.number):
        raise ValueError(f"Invalid PeakVI latent array: shape={latent.shape}")

    cell_ids = read_ids(Path(args.cell_ids))
    if latent.shape[0] != len(cell_ids):
        raise ValueError(
            f"Latent rows ({latent.shape[0]:,}) != cell IDs ({len(cell_ids):,})"
        )
    if len(np.unique(cell_ids)) != len(cell_ids):
        raise ValueError("PeakVI cell IDs are not unique")

    metadata = pd.read_csv(args.metadata_tsv, sep="\t", dtype={"cell_id": str})
    missing = {"cell_id", "cell_type"} - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    if metadata["cell_id"].duplicated().any():
        raise ValueError("Metadata cell IDs are not unique")
    metadata = metadata.set_index("cell_id").reindex(cell_ids)
    if metadata["cell_type"].isna().any():
        examples = metadata.index[metadata["cell_type"].isna()][:5]
        raise ValueError(
            "PeakVI cells missing cell_type metadata: " + ", ".join(examples)
        )
    labels = metadata["cell_type"].astype(str).to_numpy()

    selected = {}
    for value in args.selected_cells:
        path = Path(value)
        subset = subset_from_filename(path)
        if subset in selected:
            raise ValueError(f"Duplicate selected-cell input for {subset}")
        selected[subset] = read_ids(path)
        if len(np.unique(selected[subset])) != len(selected[subset]):
            raise ValueError(f"Duplicate selected cells in {path}")
    missing_subsets = set(SUBSETS) - set(selected)
    if missing_subsets:
        raise ValueError(f"Missing subsets: {sorted(missing_subsets)}")
    return latent, cell_ids, labels, selected


def knn_metrics(
    latent: np.ndarray,
    label_codes: np.ndarray,
    k: int,
    threads: int,
    random_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if k < 1 or k >= len(latent):
        raise ValueError("--knn_k must be between 1 and n_cells - 1")
    print(
        f"Building NN-descent graph: {len(latent):,} cells x "
        f"{latent.shape[1]} dimensions, k={k}",
        flush=True,
    )
    index = NNDescent(
        latent,
        n_neighbors=k + 1,
        metric="euclidean",
        random_state=random_seed,
        n_jobs=threads,
        low_memory=True,
        verbose=True,
    )
    neighbors, distances = index.neighbor_graph
    purity = np.empty(len(latent), dtype=np.float32)
    mean_distance = np.empty(len(latent), dtype=np.float32)
    radius = np.empty(len(latent), dtype=np.float32)
    self_rows = 0

    for start in range(0, len(latent), 100_000):
        stop = min(start + 100_000, len(latent))
        block = neighbors[start:stop]
        distance_block = distances[start:stop]
        rows = np.arange(start, stop, dtype=block.dtype)[:, None]
        non_self = block != rows
        self_rows += int(np.any(~non_self, axis=1).sum())

        # Retain exactly the first k non-self neighbors.
        use = non_self & (np.cumsum(non_self, axis=1) <= k)
        if np.any(use.sum(axis=1) != k):
            raise RuntimeError("NN-descent returned fewer than k non-self neighbors")
        same_type = label_codes[block] == label_codes[start:stop, None]
        purity[start:stop] = (same_type & use).sum(axis=1) / k

        used_distances = np.where(use, distance_block, np.nan)
        mean_distance[start:stop] = np.nanmean(used_distances, axis=1)
        radius[start:stop] = np.nanmax(used_distances, axis=1)

    print(
        f"Neighbor rows containing self: {self_rows / len(latent):.4%}",
        flush=True,
    )
    return purity, mean_distance, radius


def stratified_indices(
    label_codes: np.ndarray, max_cells: int, random_seed: int
) -> np.ndarray:
    if max_cells >= len(label_codes):
        return np.arange(len(label_codes))
    n_types = int(label_codes.max()) + 1
    if max_cells < n_types:
        raise ValueError("--silhouette_max_cells is smaller than n_cell_types")

    counts = np.bincount(label_codes, minlength=n_types)
    minimum = np.minimum(counts, 2)
    if minimum.sum() > max_cells:
        minimum = np.minimum(counts, 1)
    remaining = max_cells - int(minimum.sum())
    capacity = counts - minimum
    exact = remaining * capacity / capacity.sum() if capacity.sum() else capacity
    extra = np.minimum(np.floor(exact).astype(int), capacity)
    leftover = remaining - int(extra.sum())
    for code in np.argsort(-(exact - extra)):
        if leftover == 0:
            break
        if extra[code] < capacity[code]:
            extra[code] += 1
            leftover -= 1

    rng = np.random.default_rng(random_seed)
    sampled = []
    for code, quota in enumerate(minimum + extra):
        if quota:
            candidates = np.flatnonzero(label_codes == code)
            sampled.append(rng.choice(candidates, int(quota), replace=False))
    return np.sort(np.concatenate(sampled))


def silhouette_by_type(
    latent: np.ndarray,
    label_codes: np.ndarray,
    n_types: int,
    max_cells: int,
    random_seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    indices = stratified_indices(label_codes, max_cells, random_seed)
    sampled_labels = label_codes[indices]
    sampled_n = np.bincount(sampled_labels, minlength=n_types)
    means = np.full(n_types, np.nan)
    print(
        f"Computing silhouette on a stratified sample of {len(indices):,} cells",
        flush=True,
    )
    try:
        values = silhouette_samples(
            np.asarray(latent[indices], dtype=np.float32),
            sampled_labels,
            metric="euclidean",
        )
        for code in range(n_types):
            current = values[sampled_labels == code]
            if len(current):
                means[code] = current.mean()
    except Exception as error:
        # kNN purity is primary, so silhouette failure is non-fatal.
        print(f"WARNING: silhouette failed; writing NA: {error}", flush=True)
    return means, sampled_n


def status(values: np.ndarray) -> np.ndarray:
    return np.where(
        values == 0,
        "absent",
        np.where(values < 10, "underpowered", "evaluable"),
    )


def build_table(
    cell_ids: np.ndarray,
    labels: np.ndarray,
    purity: np.ndarray,
    mean_distance: np.ndarray,
    radius: np.ndarray,
    silhouette: np.ndarray,
    silhouette_n: np.ndarray,
    selected: dict[str, np.ndarray],
) -> pd.DataFrame:
    cell_types, codes = np.unique(labels, return_inverse=True)
    full_n = np.bincount(codes, minlength=len(cell_types))
    purity_mean = np.bincount(codes, weights=purity) / full_n

    # Higher percentiles indicate cells in sparser latent-space neighborhoods.
    mean_distance_percentile = (
        pd.Series(mean_distance).rank(method="average", pct=True).to_numpy() * 100
    )
    radius_percentile = (
        pd.Series(radius).rank(method="average", pct=True).to_numpy() * 100
    )

    def grouped(values: np.ndarray, reducer) -> np.ndarray:
        return np.array(
            [reducer(values[codes == code]) for code in range(len(cell_types))],
            dtype=float,
        )

    table = pd.DataFrame(
        {
            "cell_type": cell_types,
            "full_n": full_n,
            "full_fraction": full_n / len(cell_ids),
            "knn_purity": purity_mean,
            "mean_knn_distance": grouped(mean_distance, np.mean),
            "median_knn_distance": grouped(mean_distance, np.median),
            "median_knn_distance_percentile": grouped(
                mean_distance_percentile, np.median
            ),
            "q90_knn_distance_percentile": grouped(
                mean_distance_percentile, lambda x: np.quantile(x, 0.90)
            ),
            "mean_knn_radius": grouped(radius, np.mean),
            "median_knn_radius": grouped(radius, np.median),
            "median_knn_radius_percentile": grouped(radius_percentile, np.median),
            "q90_knn_radius_percentile": grouped(
                radius_percentile, lambda x: np.quantile(x, 0.90)
            ),
            "silhouette": silhouette,
            "silhouette_sample_n": silhouette_n,
        }
    )
    cell_index = pd.Index(cell_ids)
    for subset in SUBSETS:
        positions = cell_index.get_indexer(selected[subset])
        if np.any(positions < 0):
            examples = selected[subset][positions < 0][:5]
            raise ValueError(
                f"{subset} selected cells absent from PeakVI: " + ", ".join(examples)
            )
        selected_n = np.bincount(codes[positions], minlength=len(cell_types))
        subset_fraction = selected_n / len(positions)
        full_fraction = full_n / len(cell_ids)
        table[f"selected_n_{subset}"] = selected_n
        table[f"selected_fraction_of_type_{subset}"] = selected_n / full_n
        table[f"subset_fraction_{subset}"] = subset_fraction
        table[f"enrichment_{subset}"] = subset_fraction / full_fraction
        table[f"status_{subset}"] = status(selected_n)
    return table.sort_values(["full_n", "cell_type"]).reset_index(drop=True)


def build_summary(atlas: str, table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subset in SUBSETS:
        status_column = f"status_{subset}"
        row = {"atlas": atlas, "subset_size": subset}
        for current in STATUSES:
            mask = table[status_column] == current
            row[f"n_{current}_cell_types"] = int(mask.sum())
            row[f"median_knn_purity_{current}"] = (
                float(table.loc[mask, "knn_purity"].median()) if mask.any() else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def scaled_sizes(values: pd.Series) -> np.ndarray:
    values = np.sqrt(values.to_numpy(dtype=float))
    if values.max() == values.min():
        return np.full(len(values), 80.0)
    return 28 + 240 * (values - values.min()) / (values.max() - values.min())


def labels_to_show(
    table: pd.DataFrame,
    subset: str,
    threshold: float,
    max_labels: int,
) -> pd.DataFrame:
    status_column = f"status_{subset}"
    candidates = table[
        table[status_column].isin(("absent", "underpowered"))
    ].copy()
    candidates["above_threshold"] = candidates["knn_purity"] >= threshold
    candidates["status_priority"] = candidates[status_column].map(
        {"absent": 0, "underpowered": 1}
    )
    return candidates.sort_values(
        ["above_threshold", "knn_purity", "status_priority", "full_n"],
        ascending=[False, False, True, True],
    ).head(max_labels)


def plot_subset(
    atlas: str,
    table: pd.DataFrame,
    subset: str,
    output_path: Path,
    threshold: float,
    max_labels: int,
) -> None:
    status_column = f"status_{subset}"
    enrichment_column = f"enrichment_{subset}"
    sizes = scaled_sizes(table["full_n"])
    figure, axis = plt.subplots(figsize=(11.5, 7.5), dpi=220)
    for current in STATUSES:
        mask = table[status_column] == current
        axis.scatter(
            table.loc[mask, "knn_purity"],
            table.loc[mask, enrichment_column],
            s=sizes[mask.to_numpy()],
            color=COLORS[current],
            alpha=0.76,
            edgecolor="white",
            linewidth=0.45,
            label=f"{current} (n={int(mask.sum())})",
        )
    axis.axhline(1, color="#333333", linestyle="--", linewidth=1)
    axis.axvline(threshold, color="#777777", linestyle=":", linewidth=1)

    offsets = ((6, 7), (6, -12), (-6, 7), (-6, -12))
    for number, (_, row) in enumerate(
        labels_to_show(table, subset, threshold, max_labels).iterrows()
    ):
        offset = offsets[number % len(offsets)]
        axis.annotate(
            row["cell_type"],
            (row["knn_purity"], row[enrichment_column]),
            xytext=offset,
            textcoords="offset points",
            ha="left" if offset[0] > 0 else "right",
            va="bottom" if offset[1] > 0 else "top",
            fontsize=8,
            arrowprops={"arrowstyle": "-", "color": "#999999", "lw": 0.55},
        )
    axis.set_xlim(-0.02, 1.02)
    axis.set_ylim(bottom=-0.03)
    axis.set_xlabel("Mean 30-nearest-neighbor purity in PeakVI latent space")
    axis.set_ylabel("Apricot subset enrichment ratio")
    axis.set_title(
        f"{atlas}: PeakVI separability vs Apricot retention ({SUBSET_LABELS[subset]})",
        fontsize=15,
        fontweight="bold",
    )
    axis.text(
        0.01,
        1.01,
        "Labels prioritize absent or underpowered types with high purity; point size reflects full-atlas cell count.",
        transform=axis.transAxes,
        fontsize=9.5,
        color="#444444",
    )
    axis.legend(frameon=False, loc="upper right")
    axis.grid(color="#E5E5E5", linewidth=0.6, alpha=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if not 0 <= args.high_purity_threshold <= 1:
        raise ValueError("--high_purity_threshold must be between 0 and 1")
    if args.max_labels < 0 or args.threads < 1:
        raise ValueError("--max_labels must be nonnegative and --threads positive")

    latent, cell_ids, labels, selected = load_inputs(args)
    cell_types, codes = np.unique(labels, return_inverse=True)
    print(
        f"Loaded {args.atlas}: {latent.shape[0]:,} cells x {latent.shape[1]} "
        f"dimensions; {len(cell_types)} retrospective cell types",
        flush=True,
    )
    purity, mean_distance, radius = knn_metrics(
        latent, codes, args.knn_k, args.threads, args.random_seed
    )
    sparsity_output = Path(args.output_cell_sparsity_npz)
    sparsity_output.parent.mkdir(parents=True, exist_ok=True)
    mean_distance_percentile = (
        pd.Series(mean_distance).rank(method="average", pct=True).to_numpy() * 100
    ).astype(np.float32)
    radius_percentile = (
        pd.Series(radius).rank(method="average", pct=True).to_numpy() * 100
    ).astype(np.float32)
    np.savez_compressed(
        sparsity_output,
        mean_knn_distance=mean_distance.astype(np.float32),
        knn_radius=radius.astype(np.float32),
        mean_knn_distance_percentile=mean_distance_percentile,
        knn_radius_percentile=radius_percentile,
        cell_name_sha256=np.array(cell_name_checksum(cell_ids)),
        knn_k=np.array(args.knn_k, dtype=np.int32),
    )
    print(f"Wrote label-blind per-cell sparsity metrics: {sparsity_output}", flush=True)

    silhouette, silhouette_n = silhouette_by_type(
        latent,
        codes,
        len(cell_types),
        args.silhouette_max_cells,
        args.random_seed,
    )
    table = build_table(
        cell_ids,
        labels,
        purity,
        mean_distance,
        radius,
        silhouette,
        silhouette_n,
        selected,
    )
    summary = build_summary(args.atlas, table)

    output_table = Path(args.output_table)
    output_summary = Path(args.output_summary)
    figures_dir = Path(args.output_figures_dir)
    output_table.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_table, sep="\t", index=False)
    summary.to_csv(output_summary, sep="\t", index=False)
    for subset in SUBSETS:
        plot_subset(
            args.atlas,
            table,
            subset,
            figures_dir / f"separability_vs_enrichment_{subset}.png",
            args.high_purity_threshold,
            args.max_labels,
        )

    print(f"Wrote: {output_table}", flush=True)
    print(f"Wrote: {output_summary}", flush=True)
    print(f"Wrote five PNGs under: {figures_dir}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
