#!/usr/bin/env python3
"""Diagnose greedy competition along an existing hybrid 2% selection path.

This script is deliberately read-only with respect to selection. It replays the
already selected cell order and reconstructs the exact FeatureBasedSelection
square-root objective state. Cell-type annotations are loaded only after the
selection trajectory has been read and validated.

Best gains within evaluated cell types are exact. Global ranks are sample-based
estimates with stale-gain upper bounds: gains for a fixed annotation-blind random probe
are refreshed at regular checkpoints, and stale gains between checkpoints are
valid upper bounds because feature-based marginal gains can only decrease.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
from pathlib import Path

import numpy as np
import pandas as pd


KIDNEY_FOCAL_TYPES = ("Macula_Densa", "Mes")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature_matrix", required=True)
    parser.add_argument("--cell_names", required=True)
    parser.add_argument("--sparsity_npz", required=True)
    parser.add_argument("--metadata_tsv", required=True)
    parser.add_argument("--selected_cells", required=True)
    parser.add_argument("--selected_metrics", required=True)
    parser.add_argument("--atlas", required=True, choices=["kidney_atac", "fetal_atac"])
    parser.add_argument("--output_trajectory", required=True)
    parser.add_argument("--output_summary", required=True)
    parser.add_argument("--output_png", required=True)
    parser.add_argument("--output_done", required=True)
    parser.add_argument("--cell_type_column", default="cell_type")
    parser.add_argument("--max_weight", type=float, default=2.0)
    parser.add_argument("--comparison_types", type=int, default=4)
    parser.add_argument("--rank_probe_size", type=int, default=4096)
    parser.add_argument("--rank_checkpoints", type=int, default=128)
    parser.add_argument("--gain_chunk_size", type=int, default=1024)
    parser.add_argument("--random_seed", type=int, default=123)
    return parser.parse_args()


def read_ids(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any():
        raise ValueError(f"Empty or missing cell IDs: {path}")
    result = values.to_numpy(dtype=str)
    if len(np.unique(result)) != len(result):
        raise ValueError(f"Duplicate cell IDs: {path}")
    return result


def checksum(cell_ids: np.ndarray) -> str:
    digest = hashlib.sha256()
    for cell_id in cell_ids:
        digest.update(str(cell_id).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


class HybridObjective:
    """Evaluate concat(X, wX) square-root marginal gains without materializing it."""

    def __init__(self, features: np.ndarray, weights: np.ndarray, chunk_size: int):
        self.features = features
        self.weights = weights.astype(np.float32, copy=False)
        self.chunk_size = chunk_size
        self.coverage_x = np.zeros(features.shape[1], dtype=np.float64)
        self.coverage_wx = np.zeros(features.shape[1], dtype=np.float64)

    def gains(self, positions: np.ndarray) -> np.ndarray:
        positions = np.asarray(positions, dtype=np.int64)
        output = np.empty(len(positions), dtype=np.float64)
        sqrt_x = np.sqrt(self.coverage_x)
        sqrt_wx = np.sqrt(self.coverage_wx)
        for start in range(0, len(positions), self.chunk_size):
            stop = min(start + self.chunk_size, len(positions))
            current_positions = positions[start:stop]
            block32 = np.asarray(self.features[current_positions], dtype=np.float32)
            block = block32.astype(np.float64)
            first = np.sqrt(block + self.coverage_x).sum(axis=1) - sqrt_x.sum()
            # The original hybrid run materialized wX in float32.
            weighted = np.multiply(
                block32,
                self.weights[current_positions, None],
            ).astype(np.float64)
            second = np.sqrt(weighted + self.coverage_wx).sum(axis=1) - sqrt_wx.sum()
            output[start:stop] = first + second
        return output

    def add(self, position: int) -> None:
        row32 = np.asarray(self.features[position], dtype=np.float32)
        row = row32.astype(np.float64)
        self.coverage_x += row
        self.coverage_wx += np.multiply(row32, self.weights[position]).astype(np.float64)


class ExactTypeBestTracker:
    """Find an exact within-type maximum using diminishing-return upper bounds."""

    def __init__(self, positions: np.ndarray, objective: HybridObjective):
        self.positions = np.asarray(positions, dtype=np.int64)
        initial = objective.gains(self.positions)
        self.heap = [(-float(gain), int(position), 0) for gain, position in zip(initial, self.positions)]
        heapq.heapify(self.heap)

    def best(
        self,
        iteration: int,
        available: np.ndarray,
        objective: HybridObjective,
    ) -> tuple[int | None, float]:
        while self.heap:
            negative_bound, position, evaluated_iteration = heapq.heappop(self.heap)
            if not available[position]:
                continue
            if evaluated_iteration == iteration:
                heapq.heappush(
                    self.heap,
                    (negative_bound, position, evaluated_iteration),
                )
                return position, -negative_bound
            current_gain = float(objective.gains(np.asarray([position]))[0])
            heapq.heappush(self.heap, (-current_gain, position, iteration))
        return None, float("nan")


def choose_evaluation_types(
    atlas: str,
    labels: np.ndarray,
    selected_positions: np.ndarray,
    comparison_count: int,
) -> pd.DataFrame:
    cell_types, full_counts = np.unique(labels, return_counts=True)
    selected_counts = pd.Series(labels[selected_positions]).value_counts()
    table = pd.DataFrame({"cell_type": cell_types, "full_atlas_n": full_counts})
    table["final_2pct_selected_n"] = (
        table["cell_type"].map(selected_counts).fillna(0).astype(int)
    )
    table["final_status"] = np.where(
        table["final_2pct_selected_n"] == 0, "absent", "retained"
    )

    missed = table.loc[table["final_status"] == "absent", "cell_type"].tolist()
    if atlas == "kidney_atac":
        missing_required = set(KIDNEY_FOCAL_TYPES) - set(cell_types)
        if missing_required:
            raise ValueError(f"Kidney annotations are missing: {sorted(missing_required)}")
        focal = list(dict.fromkeys([*KIDNEY_FOCAL_TYPES, *missed]))
    else:
        focal = missed
    if not focal:
        raise ValueError(f"No missed 2% cell types were found for {atlas}")

    focal_counts = table.set_index("cell_type").loc[focal, "full_atlas_n"].to_numpy()
    retained = table[
        (table["final_status"] == "retained") & (~table["cell_type"].isin(focal))
    ].copy()
    retained["match_distance"] = retained["full_atlas_n"].map(
        lambda count: np.min(np.abs(np.log1p(count) - np.log1p(focal_counts)))
    )
    comparison = (
        retained.sort_values(["match_distance", "full_atlas_n", "cell_type"])
        .head(comparison_count)["cell_type"]
        .tolist()
    )

    roles = []
    for cell_type in focal:
        if cell_type in missed:
            role = "missed"
        else:
            role = "kidney_explicit_retained"
        roles.append({"cell_type": cell_type, "evaluation_group": role})
    roles.extend(
        {"cell_type": cell_type, "evaluation_group": "matched_rare_retained"}
        for cell_type in comparison
    )
    return pd.DataFrame(roles).merge(table, on="cell_type", how="left", validate="one_to_one")


def checkpoint_iterations(n_iterations: int, n_checkpoints: int) -> set[int]:
    count = min(max(n_checkpoints, 1), n_iterations)
    points = np.linspace(1, n_iterations, num=count, dtype=int)
    return set(points.tolist()) | {1, n_iterations}


def main() -> None:
    args = parse_args()
    if not np.isclose(args.max_weight, 2.0):
        raise ValueError("This frozen diagnostic requires --max_weight 2.0")
    if args.comparison_types < 1 or args.rank_probe_size < 100 or args.rank_checkpoints < 2:
        raise ValueError("Invalid comparison/rank diagnostic settings")
    if args.gain_chunk_size < 1:
        raise ValueError("--gain_chunk_size must be >= 1")

    cell_ids = read_ids(Path(args.cell_names))
    cell_index = pd.Index(cell_ids)
    features = np.load(args.feature_matrix, mmap_mode="r")
    if features.ndim != 2 or features.shape[0] != len(cell_ids):
        raise ValueError("RBF feature rows do not match cell names")
    if float(features.min()) < 0 or not np.isfinite(features).all():
        raise ValueError("RBF feature matrix contains invalid values")

    with np.load(args.sparsity_npz, allow_pickle=False) as sparsity:
        percentile = np.asarray(
            sparsity["mean_knn_distance_percentile"], dtype=np.float32
        )
        expected_checksum = str(sparsity["cell_name_sha256"].item())
    if len(percentile) != len(cell_ids) or expected_checksum != checksum(cell_ids):
        raise ValueError("Sparsity metrics and RBF features use different cell orders")
    weights = np.minimum(
        np.float32(1.0) + np.square(percentile / np.float32(100.0)),
        np.float32(args.max_weight),
    )

    # The completed selection is loaded before labels. It is never modified.
    selected_ids = read_ids(Path(args.selected_cells))
    selected_positions = cell_index.get_indexer(selected_ids)
    if np.any(selected_positions < 0):
        raise ValueError("Hybrid selection contains unknown cell IDs")
    expected_size = max(int(round(len(cell_ids) * 0.02)), 1)
    if len(selected_positions) != expected_size:
        raise ValueError(
            f"Hybrid selection has {len(selected_positions):,} cells; expected {expected_size:,}"
        )

    metrics = pd.read_csv(args.selected_metrics, sep="\t", dtype={"cell_id": str})
    required_metrics = {"selected_rank", "cell_id", "hybrid_selection_gain"}
    if required_metrics - set(metrics.columns):
        raise ValueError("Selected metrics table lacks required trajectory columns")
    metrics = metrics.sort_values("selected_rank").reset_index(drop=True)
    if not np.array_equal(metrics["cell_id"].to_numpy(dtype=str), selected_ids):
        raise ValueError("Selected-cell and selected-metrics trajectories differ")

    # Labels enter only now, after the fixed trajectory has been established.
    metadata = pd.read_csv(args.metadata_tsv, sep="\t", dtype={"cell_id": str})
    required_metadata = {"cell_id", args.cell_type_column}
    missing_metadata = required_metadata - set(metadata.columns)
    if missing_metadata:
        raise ValueError(f"Metadata is missing: {sorted(missing_metadata)}")
    labels = (
        metadata.set_index("cell_id")
        .reindex(cell_ids)[args.cell_type_column]
    )
    if labels.isna().any():
        raise ValueError("Some cells lack cell-type annotations")
    labels = labels.astype(str).to_numpy()

    evaluation = choose_evaluation_types(
        args.atlas, labels, selected_positions, args.comparison_types
    )
    print("Retrospectively evaluated cell types:", flush=True)
    print(evaluation.to_string(index=False), flush=True)

    objective = HybridObjective(features, weights, args.gain_chunk_size)
    trackers = {
        row.cell_type: ExactTypeBestTracker(
            np.flatnonzero(labels == row.cell_type), objective
        )
        for row in evaluation.itertuples(index=False)
    }
    available = np.ones(len(cell_ids), dtype=bool)
    type_remaining = pd.Series(labels).value_counts().astype(int).to_dict()

    rng = np.random.default_rng(args.random_seed)
    probe_size = min(args.rank_probe_size, len(cell_ids))
    probe_positions = np.sort(rng.choice(len(cell_ids), size=probe_size, replace=False))
    probe_gains = np.full(probe_size, np.nan, dtype=np.float64)
    checkpoints = checkpoint_iterations(len(selected_positions), args.rank_checkpoints)
    last_checkpoint = 0
    trajectory_rows: list[dict[str, object]] = []

    for iteration, selected_position in enumerate(selected_positions, start=1):
        selected_gain = float(objective.gains(np.asarray([selected_position]))[0])
        logged_gain = float(metrics.loc[iteration - 1, "hybrid_selection_gain"])
        if not np.isclose(selected_gain, logged_gain, rtol=2e-5, atol=2e-6):
            raise RuntimeError(
                f"Replay gain mismatch at iteration {iteration}: "
                f"computed={selected_gain}, logged={logged_gain}"
            )

        if iteration in checkpoints:
            probe_available = available[probe_positions]
            probe_gains[:] = np.nan
            probe_gains[probe_available] = objective.gains(
                probe_positions[probe_available]
            )
            last_checkpoint = iteration

        n_remaining = len(cell_ids) - iteration + 1
        probe_available = available[probe_positions]
        n_probe_remaining = int(probe_available.sum())
        for row in evaluation.itertuples(index=False):
            best_position, best_gain = trackers[row.cell_type].best(
                iteration, available, objective
            )
            n_type_remaining = int(type_remaining[row.cell_type])
            if best_position is None:
                rank_estimate = np.nan
                gain_ratio = np.nan
                gain_difference = np.nan
                best_id = ""
            else:
                n_probe_higher = int(
                    np.sum(probe_available & (probe_gains > best_gain))
                )
                rank_estimate = 1.0 + (
                    n_probe_higher / max(n_probe_remaining, 1)
                ) * max(n_remaining - 1, 0)
                gain_ratio = best_gain / selected_gain if selected_gain > 0 else np.nan
                gain_difference = selected_gain - best_gain
                best_id = cell_ids[best_position]

            trajectory_rows.append(
                {
                    "atlas": args.atlas,
                    "iteration": iteration,
                    "budget_fraction_consumed": iteration / len(selected_positions),
                    "evaluated_cell_type": row.cell_type,
                    "evaluation_group": row.evaluation_group,
                    "n_remaining_cells_from_type": n_type_remaining,
                    "best_candidate_cell_id": best_id,
                    "best_candidate_gain": best_gain,
                    "selected_cell_id": cell_ids[selected_position],
                    "selected_gain": selected_gain,
                    "gain_ratio": gain_ratio,
                    "gain_difference": gain_difference,
                    "best_candidate_rank": rank_estimate,
                    "rank_is_exact": False,
                    "rank_method": "sampled_stale_upper_bound_rank_estimate",
                    "rank_probe_size_remaining": n_probe_remaining,
                    "rank_gain_checkpoint_iteration": last_checkpoint,
                }
            )

        available[selected_position] = False
        selected_type = labels[selected_position]
        type_remaining[selected_type] -= 1
        objective.add(int(selected_position))
        if iteration == 1 or iteration % 500 == 0 or iteration == len(selected_positions):
            print(f"Replayed {iteration:,}/{len(selected_positions):,} iterations", flush=True)

    trajectory = pd.DataFrame(trajectory_rows)
    summary_rows = []
    for row in evaluation.itertuples(index=False):
        current = trajectory[trajectory["evaluated_cell_type"] == row.cell_type]
        finite = current[np.isfinite(current["gain_ratio"])]
        if finite.empty:
            closest_iteration = np.nan
            minimum_rank = np.nan
            median_rank = np.nan
            best_ratio = np.nan
            median_ratio = np.nan
        else:
            closest = finite.loc[finite["gain_ratio"].idxmax()]
            closest_iteration = int(closest["iteration"])
            minimum_rank = float(finite["best_candidate_rank"].min())
            median_rank = float(finite["best_candidate_rank"].median())
            best_ratio = float(finite["gain_ratio"].max())
            median_ratio = float(finite["gain_ratio"].median())
        summary_rows.append(
            {
                "atlas": args.atlas,
                "cell_type": row.cell_type,
                "evaluation_group": row.evaluation_group,
                "full_atlas_n": int(row.full_atlas_n),
                "final_2pct_selected_n": int(row.final_2pct_selected_n),
                "final_status": row.final_status,
                "minimum_best_candidate_rank": minimum_rank,
                "median_best_candidate_rank": median_rank,
                "best_gain_ratio": best_ratio,
                "median_gain_ratio": median_ratio,
                "iteration_of_closest_competition": closest_iteration,
                "rank_is_exact": False,
                "rank_method": "sampled_stale_upper_bound_rank_estimate",
                "rank_probe_size": probe_size,
                "rank_checkpoints": len(checkpoints),
            }
        )
    summary = pd.DataFrame(summary_rows)

    output_trajectory = Path(args.output_trajectory)
    output_summary = Path(args.output_summary)
    output_png = Path(args.output_png)
    output_done = Path(args.output_done)
    for output in (output_trajectory, output_summary, output_png, output_done):
        output.parent.mkdir(parents=True, exist_ok=True)
    trajectory.to_csv(output_trajectory, sep="\t", index=False)
    summary.to_csv(output_summary, sep="\t", index=False)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(13, 7.5), dpi=220)
    for row in evaluation.itertuples(index=False):
        current = trajectory[trajectory["evaluated_cell_type"] == row.cell_type]
        linewidth = 2.0 if row.evaluation_group != "matched_rare_retained" else 1.1
        alpha = 0.95 if row.evaluation_group != "matched_rare_retained" else 0.65
        linestyle = "-" if row.final_status == "absent" else "--"
        axis.plot(
            current["budget_fraction_consumed"] * 100,
            current["gain_ratio"],
            label=f"{row.cell_type} ({row.final_status})",
            linewidth=linewidth,
            alpha=alpha,
            linestyle=linestyle,
        )
    axis.axhline(1.0, color="black", linestyle=":", linewidth=1.2)
    axis.set_xlabel("Percentage of the 2% selection budget consumed")
    axis.set_ylabel("Best cell-type candidate gain / selected-cell gain")
    axis.set_title(f"{args.atlas}: greedy competition during hybrid 2% selection")
    axis.text(
        0.01,
        0.02,
        "Values near 1 indicate close competition with the selected cell.\n"
        "Solid lines are missed populations; dashed lines are retained comparisons.",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
    )
    axis.grid(alpha=0.22, linewidth=0.6)
    axis.legend(loc="upper right", fontsize=8, frameon=False, ncol=2)
    figure.tight_layout()
    figure.savefig(output_png, dpi=320, facecolor="white")
    plt.close(figure)

    output_done.write_text("DONE\n", encoding="utf-8")
    print(f"Wrote trajectory: {output_trajectory}", flush=True)
    print(f"Wrote summary: {output_summary}", flush=True)
    print(f"Wrote figure: {output_png}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
