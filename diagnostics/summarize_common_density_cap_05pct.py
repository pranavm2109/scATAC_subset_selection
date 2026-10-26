#!/usr/bin/env python3
"""Summarize a shared density-cap sweep across all atlases at 5%."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summaries", nargs="+", required=True)
    parser.add_argument("--details", nargs="+", required=True)
    parser.add_argument("--output_summary", required=True)
    parser.add_argument("--output_recommendation", required=True)
    parser.add_argument("--output_done", required=True)
    return parser.parse_args()


def cap_from_method(method: str) -> float:
    match = re.fullmatch(r"density_max_([0-9]+)p([0-9]+)", str(method))
    if match is None:
        raise ValueError(f"Cannot parse density cap from method: {method}")
    return float(f"{match.group(1)}.{match.group(2)}")


def main() -> None:
    args = parse_args()
    summaries = pd.concat(
        [pd.read_csv(path, sep="\t") for path in args.summaries],
        ignore_index=True,
    )
    details = pd.concat(
        [pd.read_csv(path, sep="\t") for path in args.details],
        ignore_index=True,
    )

    summaries = summaries[summaries["method"].str.startswith("density_max_")].copy()
    details = details[details["method"].str.startswith("density_max_")].copy()
    if summaries.empty or details.empty:
        raise ValueError("No density-cap results were found")
    if set(summaries["subset_size"].astype(str)) != {"05pct"}:
        raise ValueError("The common-cap diagnostic must contain only 05pct results")

    summaries["density_max_weight"] = summaries["method"].map(cap_from_method)
    details["density_max_weight"] = details["method"].map(cap_from_method)
    atlases = sorted(summaries["atlas"].astype(str).unique())

    rows = []
    for cap, cap_summary in summaries.groupby("density_max_weight", sort=True):
        cap_details = details[details["density_max_weight"] == cap]
        absent = cap_details[~cap_details["retained"].astype(bool)]
        absent_by_atlas = []
        for atlas in atlases:
            names = sorted(
                absent.loc[absent["atlas"] == atlas, "cell_type"].astype(str)
            )
            absent_by_atlas.append(f"{atlas}=" + (";".join(names) if names else "NONE"))

        complete_atlases = int((cap_summary["n_cell_types_retained"] == cap_summary["n_cell_types_total"]).sum())
        objective = pd.to_numeric(
            cap_summary["objective_relative_to_baseline"], errors="coerce"
        )
        rows.append(
            {
                "density_max_weight": cap,
                "n_atlases": len(atlases),
                "n_atlases_with_all_cell_types": complete_atlases,
                "all_atlases_retain_all_cell_types": complete_atlases == len(atlases),
                "total_absent_cell_types": int((~cap_details["retained"].astype(bool)).sum()),
                "worst_atlas_absent_cell_types": int(
                    cap_details.assign(absent=~cap_details["retained"].astype(bool))
                    .groupby("atlas")["absent"]
                    .sum()
                    .max()
                ),
                "minimum_objective_retention": float(objective.min()),
                "mean_objective_retention": float(objective.mean()),
                "absent_cell_types_by_atlas": " | ".join(absent_by_atlas),
            }
        )

    result = pd.DataFrame(rows).sort_values("density_max_weight").reset_index(drop=True)
    eligible = result[result["all_atlases_retain_all_cell_types"]].copy()
    if eligible.empty:
        recommendation = pd.DataFrame(
            {
                "status": ["NO_COMMON_CAP_FOUND"],
                "recommended_density_max_weight": [np.nan],
                "selection_rule": [
                    "No tested cap retained every annotated cell type in all atlases at 5%."
                ],
                "warning": [
                    "Do not claim complete retention; expand or redesign the diagnostic."
                ],
            }
        )
    else:
        eligible = eligible.sort_values(
            ["minimum_objective_retention", "mean_objective_retention", "density_max_weight"],
            ascending=[False, False, True],
        )
        best = eligible.iloc[0]
        recommendation = pd.DataFrame(
            {
                "status": ["COMMON_CAP_FOUND"],
                "recommended_density_max_weight": [best["density_max_weight"]],
                "selection_rule": [
                    "Among caps retaining all types in all atlases, maximize the worst-atlas objective retention."
                ],
                "warning": [
                    "Cell-type labels were used retrospectively to tune this hyperparameter; validate on held-out data."
                ],
            }
        )

    output_summary = Path(args.output_summary)
    output_recommendation = Path(args.output_recommendation)
    output_done = Path(args.output_done)
    for path in (output_summary, output_recommendation, output_done):
        path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_summary, sep="\t", index=False)
    recommendation.to_csv(output_recommendation, sep="\t", index=False)
    output_done.write_text("DONE\n", encoding="utf-8")
    print(f"Wrote common-cap summary: {output_summary}", flush=True)
    print(f"Wrote recommendation: {output_recommendation}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
