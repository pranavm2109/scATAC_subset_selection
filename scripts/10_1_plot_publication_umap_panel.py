#!/usr/bin/env python3
"""Create a publication-ready 3x5 PeakVI UMAP subset-coverage panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ATLAS_LABELS = {
    "kidney_atac": "Kidney atlas",
    "heca_atac": "HECA atlas",
    "fetal_atac": "Fetal atlas",
}
FRACTION_LABELS = {
    "02pct": "2% subset",
    "04pct": "4% subset",
    "06pct": "6% subset",
    "08pct": "8% subset",
    "10pct": "10% subset",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlases", nargs="+", required=True)
    parser.add_argument("--fractions", nargs="+", required=True)
    parser.add_argument("--umap_tsvs", nargs="+", required=True)
    parser.add_argument("--metadata_tsvs", nargs="+", required=True)
    parser.add_argument(
        "--selected_cells",
        nargs="+",
        required=True,
        help="Selected-cell files in atlas-major, fraction-minor order.",
    )
    parser.add_argument("--output_png", required=True)
    parser.add_argument("--output_pdf", required=True)
    parser.add_argument("--output_summary_tsv", required=True)
    parser.add_argument("--max_background_cells", type=int, default=250000)
    parser.add_argument("--random_seed", type=int, default=123)
    parser.add_argument("--dpi", type=int, default=400)
    return parser.parse_args()


def read_selected(path: Path) -> np.ndarray:
    values = pd.read_csv(path, header=None, dtype=str).iloc[:, 0]
    if values.empty or values.isna().any() or values.duplicated().any():
        raise ValueError(f"Selected-cell IDs must be nonempty and unique: {path}")
    return values.to_numpy(dtype=str)


def load_atlas(umap_path: Path, metadata_path: Path) -> pd.DataFrame:
    umap = pd.read_csv(umap_path, sep="\t", dtype={"cell_name": str})
    metadata = pd.read_csv(metadata_path, sep="\t", dtype={"cell_id": str})
    required_umap = {"cell_name", "UMAP_1", "UMAP_2"}
    required_metadata = {"cell_id", "cell_type"}
    if missing := required_umap.difference(umap.columns):
        raise ValueError(f"{umap_path} is missing columns: {sorted(missing)}")
    if missing := required_metadata.difference(metadata.columns):
        raise ValueError(f"{metadata_path} is missing columns: {sorted(missing)}")
    if umap["cell_name"].duplicated().any() or metadata["cell_id"].duplicated().any():
        raise ValueError("UMAP and metadata cell IDs must be unique")
    if set(umap["cell_name"]) != set(metadata["cell_id"]):
        raise ValueError(
            f"Cell IDs differ between {umap_path} and {metadata_path}"
        )
    frame = umap.merge(
        metadata[["cell_id", "cell_type"]],
        left_on="cell_name",
        right_on="cell_id",
        how="left",
        validate="one_to_one",
    )
    frame["UMAP_1"] = pd.to_numeric(frame["UMAP_1"], errors="coerce")
    frame["UMAP_2"] = pd.to_numeric(frame["UMAP_2"], errors="coerce")
    if (
        not np.isfinite(frame[["UMAP_1", "UMAP_2"]].to_numpy()).all()
        or frame["cell_type"].isna().any()
    ):
        raise ValueError("UMAP coordinates and cell-type metadata must be complete")
    return frame.drop(columns="cell_id")


def main() -> None:
    args = parse_args()
    n_atlases = len(args.atlases)
    n_fractions = len(args.fractions)
    if (n_atlases, n_fractions) != (3, 5):
        raise ValueError("The publication panel requires exactly 3 atlases and 5 fractions")
    if len(args.umap_tsvs) != n_atlases or len(args.metadata_tsvs) != n_atlases:
        raise ValueError("Supply one UMAP and metadata TSV per atlas")
    if len(args.selected_cells) != n_atlases * n_fractions:
        raise ValueError("Supply 9 selected-cell files in atlas-major order")
    if set(args.fractions) != set(FRACTION_LABELS):
        raise ValueError("Fractions must be 02pct, 04pct, 06pct, 08pct, and 10pct")
    if args.max_background_cells < 1000:
        raise ValueError("--max_background_cells must be at least 1000")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(args.random_seed)
    frames = [
        load_atlas(Path(umap_path), Path(metadata_path))
        for umap_path, metadata_path in zip(args.umap_tsvs, args.metadata_tsvs)
    ]

    figure, axes = plt.subplots(
        nrows=3,
        ncols=5,
        figsize=(25.5, 14.5),
        constrained_layout=False,
    )
    summary_rows: list[dict[str, object]] = []

    for row_index, (atlas, frame) in enumerate(zip(args.atlases, frames)):
        all_ids = frame["cell_name"].to_numpy(dtype=str)
        all_types = frame["cell_type"].astype(str)
        n_types_full = all_types.nunique()
        if len(frame) > args.max_background_cells:
            background_indices = np.sort(
                rng.choice(
                    len(frame),
                    size=args.max_background_cells,
                    replace=False,
                )
            )
            background = frame.iloc[background_indices]
        else:
            background = frame

        x_min, x_max = frame["UMAP_1"].min(), frame["UMAP_1"].max()
        y_min, y_max = frame["UMAP_2"].min(), frame["UMAP_2"].max()
        x_pad = max((x_max - x_min) * 0.03, 1e-6)
        y_pad = max((y_max - y_min) * 0.03, 1e-6)

        for column_index, fraction in enumerate(args.fractions):
            ax = axes[row_index, column_index]
            selected_path = Path(
                args.selected_cells[row_index * n_fractions + column_index]
            )
            selected_ids = read_selected(selected_path)
            missing = np.setdiff1d(selected_ids, all_ids)
            if len(missing):
                raise ValueError(
                    f"{selected_path} contains {len(missing)} IDs absent from {atlas}"
                )
            selected = frame.set_index("cell_name").loc[selected_ids].reset_index()
            n_types_selected = selected["cell_type"].astype(str).nunique()
            fraction_observed = len(selected) / len(frame)

            ax.scatter(
                background["UMAP_1"],
                background["UMAP_2"],
                s=0.45,
                c="#D7DDE1",
                alpha=0.55,
                linewidths=0,
                rasterized=True,
                zorder=1,
            )
            ax.scatter(
                selected["UMAP_1"],
                selected["UMAP_2"],
                s=1.6,
                c="#D95F0E",
                alpha=0.82,
                linewidths=0,
                rasterized=True,
                zorder=2,
            )
            ax.set_xlim(x_min - x_pad, x_max + x_pad)
            ax.set_ylim(y_min - y_pad, y_max + y_pad)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#B8C0C5")
                spine.set_linewidth(0.65)

            if row_index == 0:
                ax.set_title(
                    FRACTION_LABELS.get(fraction, fraction),
                    fontsize=15,
                    fontweight="bold",
                    pad=10,
                )
            if column_index == 0:
                ax.set_ylabel(
                    ATLAS_LABELS.get(atlas, atlas),
                    fontsize=14,
                    fontweight="bold",
                    labelpad=13,
                )

            annotation = (
                f"{len(selected):,} cells\n"
                f"{n_types_selected}/{n_types_full} cell types"
            )
            ax.text(
                0.025,
                0.975,
                annotation,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.5,
                color="#17242D",
                bbox={
                    "boxstyle": "round,pad=0.32",
                    "facecolor": "white",
                    "edgecolor": "#CAD1D5",
                    "alpha": 0.92,
                },
                zorder=3,
            )

            summary_rows.append(
                {
                    "atlas": atlas,
                    "fraction": fraction,
                    "n_atlas_cells": len(frame),
                    "n_selected_cells": len(selected),
                    "observed_fraction": fraction_observed,
                    "n_cell_types_atlas": n_types_full,
                    "n_cell_types_selected": n_types_selected,
                    "n_cell_types_missing": n_types_full - n_types_selected,
                    "selected_cells_file": str(selected_path),
                }
            )

    figure.suptitle(
        "Annotation-blind subsets preserve broad PeakVI manifold coverage",
        fontsize=21,
        fontweight="bold",
        x=0.5,
        y=0.985,
    )
    figure.text(
        0.5,
        0.952,
        (
            "Gaussian PeakVI latent space; 1,024 RBF landmarks; "
            "L1-normalized FeatureBasedSelection"
        ),
        ha="center",
        va="top",
        fontsize=11.5,
        color="#40515B",
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#D7DDE1",
            markeredgecolor="none",
            markersize=7,
            label="Full atlas",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#D95F0E",
            markeredgecolor="none",
            markersize=7,
            label="Selected subset",
        ),
    ]
    figure.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=11,
        bbox_to_anchor=(0.5, 0.018),
    )
    figure.text(
        0.01,
        0.006,
        (
            "Cell-type annotations are used only to report retention after selection; "
            "they are not used by the subset-selection algorithm."
        ),
        ha="left",
        va="bottom",
        fontsize=9,
        color="#4A5961",
    )
    figure.subplots_adjust(
        left=0.065,
        right=0.99,
        top=0.925,
        bottom=0.065,
        hspace=0.075,
        wspace=0.045,
    )

    output_png = Path(args.output_png)
    output_pdf = Path(args.output_pdf)
    output_summary = Path(args.output_summary_tsv)
    for output in (output_png, output_pdf, output_summary):
        output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=args.dpi, facecolor="white", bbox_inches="tight")
    figure.savefig(output_pdf, facecolor="white", bbox_inches="tight")
    plt.close(figure)
    pd.DataFrame(summary_rows).to_csv(output_summary, sep="\t", index=False)

    print(f"Wrote UMAP panel: {output_png}", flush=True)
    print(f"Wrote vector UMAP panel: {output_pdf}", flush=True)
    print(f"Wrote panel summary: {output_summary}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
