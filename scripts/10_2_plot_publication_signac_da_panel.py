#!/usr/bin/env python3
"""Create a publication-ready 3x5 Signac DA effect-recapitulation panel."""

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
    parser.add_argument("--details_tsvs", nargs="+", required=True)
    parser.add_argument("--output_png", required=True)
    parser.add_argument("--output_pdf", required=True)
    parser.add_argument("--output_summary_tsv", required=True)
    parser.add_argument("--display_quantile", type=float, default=0.9975)
    parser.add_argument("--dpi", type=int, default=400)
    return parser.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float).ne(0)
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "t", "1", "yes", "y"})
    )


def safe_correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    if method == "pearson":
        return float(np.corrcoef(x, y)[0, 1])
    return float(pd.Series(x).corr(pd.Series(y), method="spearman"))


def format_metric(value: float, digits: int = 2) -> str:
    return "NA" if not np.isfinite(value) else f"{value:.{digits}f}"


def load_details(path: Path, atlas: str, fractions: list[str]) -> pd.DataFrame:
    required = {
        "atlas",
        "cell_type",
        "subset_size",
        "peak",
        "full_avg_log2FC",
        "subset_avg_log2FC",
        "effect_available",
        "same_direction_effect",
        "same_direction_overlap",
        "direction_switched",
        "subset_analysis_status",
    }
    frame = pd.read_csv(
        path, sep="\t", usecols=lambda column: column in required, low_memory=False
    )
    if missing := required.difference(frame.columns):
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    frame["atlas"] = frame["atlas"].astype(str)
    frame["subset_size"] = frame["subset_size"].astype(str)
    if len(frame) and set(frame["atlas"]) != {atlas}:
        raise ValueError(f"Atlas mismatch in {path}")
    unexpected = set(frame["subset_size"]).difference(fractions)
    if unexpected:
        raise ValueError(f"Unexpected subset sizes in {path}: {sorted(unexpected)}")
    frame["full_avg_log2FC"] = pd.to_numeric(
        frame["full_avg_log2FC"], errors="coerce"
    )
    frame["subset_avg_log2FC"] = pd.to_numeric(
        frame["subset_avg_log2FC"], errors="coerce"
    )
    for column in (
        "effect_available",
        "same_direction_effect",
        "same_direction_overlap",
        "direction_switched",
    ):
        frame[column] = as_bool(frame[column])
    if frame.duplicated(["cell_type", "subset_size", "peak"]).any():
        raise ValueError(f"Duplicate atlas-reference peak rows in {path}")
    return frame


def main() -> None:
    args = parse_args()
    if len(args.atlases) != 3 or len(args.fractions) != 5:
        raise ValueError("The publication panel requires exactly 3 atlases and 5 fractions")
    if len(args.details_tsvs) != 3:
        raise ValueError("Supply one DA overlap-details TSV per atlas")
    if set(args.fractions) != set(FRACTION_LABELS):
        raise ValueError("Fractions must be 02pct, 04pct, 06pct, 08pct, and 10pct")
    if not 0.95 <= args.display_quantile < 1:
        raise ValueError("--display_quantile must be at least 0.95 and less than 1")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(
        nrows=3,
        ncols=5,
        figsize=(25.8, 14.8),
        constrained_layout=False,
    )
    summary_rows: list[dict[str, object]] = []

    for row_index, (atlas, path) in enumerate(zip(args.atlases, args.details_tsvs)):
        print(f"Loading DA details for {atlas}", flush=True)
        frame = load_details(Path(path), atlas, args.fractions)
        available_all = frame.loc[frame["effect_available"]].copy()
        absolute_effects = np.abs(
            available_all[["full_avg_log2FC", "subset_avg_log2FC"]]
            .to_numpy(dtype=float)
            .ravel()
        )
        absolute_effects = absolute_effects[np.isfinite(absolute_effects)]
        if len(absolute_effects):
            axis_limit = float(np.quantile(absolute_effects, args.display_quantile))
            axis_limit = max(axis_limit, 0.25)
        else:
            axis_limit = 1.0
        axis_limit *= 1.04

        for column_index, fraction in enumerate(args.fractions):
            ax = axes[row_index, column_index]
            panel_all = frame.loc[frame["subset_size"] == fraction].copy()
            panel = panel_all.loc[panel_all["effect_available"]].copy()
            x = panel["full_avg_log2FC"].to_numpy(dtype=float)
            y = panel["subset_avg_log2FC"].to_numpy(dtype=float)
            finite = np.isfinite(x) & np.isfinite(y)
            x = x[finite]
            y = y[finite]
            n_reference = len(panel_all)
            n_available = len(x)
            n_same = int(panel.loc[finite, "same_direction_effect"].sum())
            n_switched = int(panel.loc[finite, "direction_switched"].sum())
            n_strict = int(panel_all["same_direction_overlap"].sum())
            direction_concordance = (
                n_same / n_available if n_available else float("nan")
            )
            sign_switch_fraction = (
                n_switched / n_available if n_available else float("nan")
            )
            strict_recovery = (
                n_strict / n_reference if n_reference else float("nan")
            )
            spearman = safe_correlation(x, y, "spearman")
            pearson = safe_correlation(x, y, "pearson")
            median_error = (
                float(np.median(np.abs(y - x))) if n_available else float("nan")
            )
            n_cell_types_reference = panel_all["cell_type"].nunique()
            n_cell_types_evaluable = panel.loc[
                panel["subset_analysis_status"] == "success", "cell_type"
            ].nunique()

            if n_available:
                x_display = np.clip(x, -axis_limit, axis_limit)
                y_display = np.clip(y, -axis_limit, axis_limit)
                ax.hexbin(
                    x_display,
                    y_display,
                    gridsize=52,
                    extent=(-axis_limit, axis_limit, -axis_limit, axis_limit),
                    mincnt=1,
                    bins="log",
                    cmap="cividis",
                    linewidths=0,
                    rasterized=True,
                )
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No evaluable\natlas-reference effects",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=12,
                    color="#59666D",
                )

            ax.plot(
                [-axis_limit, axis_limit],
                [-axis_limit, axis_limit],
                linestyle="--",
                linewidth=1.0,
                color="#D95F0E",
                alpha=0.95,
                zorder=3,
            )
            ax.axhline(0, color="#AAB2B7", linewidth=0.55, zorder=0)
            ax.axvline(0, color="#AAB2B7", linewidth=0.55, zorder=0)
            ax.set_xlim(-axis_limit, axis_limit)
            ax.set_ylim(-axis_limit, axis_limit)
            ax.set_aspect("equal", adjustable="box")
            ax.tick_params(axis="both", labelsize=8.5, colors="#26343B")
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            for spine in ("left", "bottom"):
                ax.spines[spine].set_color("#8B989F")
                ax.spines[spine].set_linewidth(0.7)

            if row_index == 0:
                ax.set_title(
                    FRACTION_LABELS.get(fraction, fraction),
                    fontsize=15,
                    fontweight="bold",
                    pad=10,
                )
            if column_index == 0:
                ax.set_ylabel(
                    (
                        f"{ATLAS_LABELS.get(atlas, atlas)}\n"
                        "Subset average log2 fold change"
                    ),
                    fontsize=11.5,
                    fontweight="bold",
                    labelpad=9,
                )
            else:
                ax.set_ylabel("Subset average log2 fold change", fontsize=10)
            if row_index == 2:
                ax.set_xlabel("Full-atlas average log2 fold change", fontsize=10.5)
            else:
                ax.set_xlabel("Full-atlas average log2 fold change", fontsize=10)

            annotation = "\n".join(
                [
                    f"Spearman rho = {format_metric(spearman)}",
                    f"Direction match = {100 * direction_concordance:.1f}%"
                    if np.isfinite(direction_concordance)
                    else "Direction match = NA",
                    f"Sign switched = {n_switched:,}/{n_available:,}",
                    f"Strict recovery = {100 * strict_recovery:.1f}%"
                    if np.isfinite(strict_recovery)
                    else "Strict recovery = NA",
                ]
            )
            ax.text(
                0.035,
                0.965,
                annotation,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.1,
                linespacing=1.22,
                bbox={
                    "boxstyle": "round,pad=0.35",
                    "facecolor": "white",
                    "edgecolor": "#C9D0D4",
                    "alpha": 0.92,
                },
                zorder=4,
            )

            summary_rows.append(
                {
                    "atlas": atlas,
                    "subset_size": fraction,
                    "n_atlas_reference_peak_cell_type_pairs": n_reference,
                    "n_effects_available": n_available,
                    "fraction_effects_available": (
                        n_available / n_reference if n_reference else np.nan
                    ),
                    "n_reference_cell_types": n_cell_types_reference,
                    "n_evaluable_cell_types": n_cell_types_evaluable,
                    "spearman_rho": spearman,
                    "pearson_r": pearson,
                    "direction_concordance": direction_concordance,
                    "n_sign_switched": n_switched,
                    "sign_switch_fraction": sign_switch_fraction,
                    "n_strict_same_direction_recovered": n_strict,
                    "strict_recovery_fraction": strict_recovery,
                    "median_absolute_effect_error": median_error,
                    "display_axis_limit": axis_limit,
                }
            )

        del frame, available_all, absolute_effects, panel_all, panel, x, y, finite

    figure.suptitle(
        "Small annotation-blind subsets recapitulate atlas-level DA effect structure",
        fontsize=20.5,
        fontweight="bold",
        y=0.985,
    )
    figure.text(
        0.5,
        0.952,
        (
            "Each hexagon summarizes full-atlas significant peak-cell-type pairs; "
            "the dashed line marks exact effect agreement"
        ),
        ha="center",
        va="top",
        fontsize=11.2,
        color="#40515B",
    )
    figure.text(
        0.01,
        0.008,
        (
            "Primary metrics use subset effect estimates at atlas-defined DA peaks "
            "(BH FDR < 0.05 and |average log2FC| >= 0.25 in the full atlas). "
            "Strict recovery additionally requires subset significance and matching direction. "
            "Metrics use unmodified values; axes are clipped only for display."
        ),
        ha="left",
        va="bottom",
        fontsize=8.7,
        color="#4A5961",
        wrap=True,
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.99,
        top=0.922,
        bottom=0.085,
        hspace=0.28,
        wspace=0.21,
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

    print(f"Wrote DA panel: {output_png}", flush=True)
    print(f"Wrote vector DA panel: {output_pdf}", flush=True)
    print(f"Wrote panel summary: {output_summary}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
