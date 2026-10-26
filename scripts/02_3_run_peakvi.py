#!/usr/bin/env python3
"""
Train PeakVI on a PeakVI-ready AnnData object.

Expected input:
  - AnnData .X contains raw peak accessibility counts, cells x peaks.
  - AnnData .obs contains cell metadata such as cell_type and library_size.

Main outputs:
  - Saved PeakVI model directory.
  - Cell-by-latent-factor TSV for downstream subset selection.
  - Optional AnnData with adata.obsm['X_peakvi'] added.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def str_to_bool(value: str) -> bool:
    value = str(value).strip().lower()
    if value in {"true", "t", "1", "yes", "y"}:
        return True
    if value in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got: {value}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Train scvi-tools PeakVI on raw scATAC peak counts.")
    p.add_argument("--input_h5ad", required=True, help="Input AnnData file from 02_2 AnnData builder.")
    p.add_argument("--atlas", required=True, help="Atlas name for output metadata.")
    p.add_argument("--output_dir", required=True, help="Directory for model and auxiliary outputs.")
    p.add_argument("--latent_tsv", required=True, help="Output TSV containing cells x PeakVI latent dimensions.")
    p.add_argument("--output_h5ad", default=None, help="Optional output h5ad with obsm['X_peakvi'] added.")
    p.add_argument("--history_tsv", default=None, help="Optional TSV of PeakVI training history/loss curves.")
    p.add_argument("--convergence_png", default=None, help="Optional PNG convergence plot of training history/loss curves.")
    p.add_argument("--n_latent", type=int, default=60, help="PeakVI latent dimensionality [default: 60].")
    p.add_argument("--n_hidden", type=int, default=256, help="Hidden units per layer [default: 256].")
    p.add_argument("--n_layers_encoder", type=int, default=2, help="Encoder layers [default: 2].")
    p.add_argument("--n_layers_decoder", type=int, default=2, help="Decoder layers [default: 2].")
    p.add_argument("--dropout_rate", type=float, default=0.1, help="Dropout rate [default: 0.1].")
    p.add_argument("--latent_distribution", choices=["ln", "normal"], default="normal", help="Latent distribution for PEAKVI [default: normal].")
    p.add_argument("--max_epochs", default="auto", help="Training epochs, or auto to use scvi-tools n_obs heuristic [default: auto].")
    p.add_argument("--epoch_multiplier", type=int, default=2, help="Train for this many times the resolved heuristic/explicit max_epochs [default: 1].")
    p.add_argument("--epochs_cap", type=int, default=400, help="Epoch cap used when --max_epochs auto [default: 400].")
    p.add_argument("--decay_at_n_obs", type=int, default=20000, help="Observation count where auto epoch decay begins [default: 20000].")
    p.add_argument("--batch_size", type=int, default=1024, help="Training batch size [default: 1024].")
    p.add_argument("--batch_key", default=None, help="Optional obs column for batch correction.")
    p.add_argument("--binarize", type=str_to_bool, default=True, help="Convert counts > 0 to 1 before training [default: true].")
    p.add_argument("--early_stopping", type=str_to_bool, default=True, help="Use early stopping [default: true].")
    p.add_argument("--early_stopping_patience", type=int, default=10, help="Epochs without improvement before early stopping [default: 10].")
    p.add_argument("--dataloader_workers", type=int, default=0, help="Number of dataloader workers when supported by scvi-tools [default: 0]. Use 0 for large shared-memory-limited jobs.")
    p.add_argument("--accelerator", choices=["auto", "cpu", "gpu"], default="auto", help="Training accelerator [default: auto].")
    p.add_argument("--devices", default="auto", help="Lightning devices argument [default: auto].")
    p.add_argument("--random_seed", type=int, default=123, help="Random seed [default: 123].")
    return p.parse_args()


def import_dependencies():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp
        import torch
        import scvi
    except ImportError as e:
        missing = getattr(e, "name", str(e))
        raise SystemExit(
            "Missing required Python dependency for PeakVI: "
            f"{missing}\nInstall/activate an environment containing scvi-tools, anndata, "
            "torch, numpy, pandas, and scipy."
        ) from e
    return ad, np, pd, sp, torch, scvi


def binarize_matrix_in_place(adata, np, sp):
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr(copy=False)
        adata.X.data = (adata.X.data > 0).astype(np.float32, copy=False)
        adata.X.eliminate_zeros()
    else:
        adata.X = (adata.X > 0).astype(np.float32, copy=False)


def choose_accelerator(requested: str, torch) -> str:
    if requested == "auto":
        return "gpu" if torch.cuda.is_available() else "cpu"
    if requested == "gpu" and not torch.cuda.is_available():
        raise SystemExit("--accelerator gpu was requested, but torch.cuda.is_available() is false.")
    return requested


def resolve_max_epochs(max_epochs_arg: str, n_obs: int, epochs_cap: int, decay_at_n_obs: int) -> int:
    if str(max_epochs_arg).strip().lower() != "auto":
        try:
            max_epochs = int(max_epochs_arg)
        except ValueError as e:
            raise SystemExit(f"--max_epochs must be an integer or auto, got: {max_epochs_arg}") from e
        if max_epochs < 1:
            raise SystemExit("--max_epochs must be >= 1")
        return max_epochs

    try:
        from scvi.model._utils import get_max_epochs_heuristic

        return int(
            get_max_epochs_heuristic(
                n_obs=n_obs,
                epochs_cap=epochs_cap,
                decay_at_n_obs=decay_at_n_obs,
            )
        )
    except Exception as e:
        print(
            "Could not import scvi.model._utils.get_max_epochs_heuristic; "
            f"using local equivalent formula. Reason: {e}",
            flush=True,
        )
        return max(1, min(round((decay_at_n_obs / n_obs) * epochs_cap), epochs_cap))


def train_peakvi_with_compat(
    model,
    args: argparse.Namespace,
    accelerator: str,
    precomputed_max_epochs: int,
):
    """Handle small scvi-tools / Lightning API differences across versions."""
    train_kwargs = {
        "max_epochs": precomputed_max_epochs,
        "batch_size": args.batch_size,
        "early_stopping": args.early_stopping,
    }

    if args.early_stopping:
        train_kwargs["early_stopping_patience"] = args.early_stopping_patience

    optional_speed_kwargs = {
        "datasplitter_kwargs": {"num_workers": args.dataloader_workers},
    }

    # Modern scvi-tools versions use Lightning's accelerator/devices args.
    try:
        model.train(
            **train_kwargs,
            **optional_speed_kwargs,
            accelerator=accelerator,
            devices=args.devices,
        )
        return
    except TypeError:
        pass

    try:
        model.train(
            **train_kwargs,
            accelerator=accelerator,
            devices=args.devices,
        )
        return
    except TypeError:
        pass

    # Older scvi-tools versions used use_gpu.
    use_gpu = accelerator == "gpu"
    try:
        model.train(**train_kwargs, **optional_speed_kwargs, use_gpu=use_gpu)
        return
    except TypeError:
        pass

    try:
        model.train(**train_kwargs, use_gpu=use_gpu)
        return
    except TypeError:
        # Last fallback: do not pass backend-specific args.
        model.train(**train_kwargs)


def training_history_to_dataframe(model, pd):
    """Convert scvi-tools training history into a long data frame when available."""
    history = getattr(model, "history", None)
    rows = []

    if history is None:
        return pd.DataFrame(columns=["metric", "epoch", "value"])

    def add_values(metric_name, values):
        if hasattr(values, "to_frame"):
            values = values.to_frame()

        if hasattr(values, "reset_index") and hasattr(values, "columns"):
            df = values.reset_index()
            epoch_col = df.columns[0]
            value_cols = [col for col in df.columns if col != epoch_col]
            for _, row in df.iterrows():
                for value_col in value_cols:
                    value = row[value_col]
                    if pd.notna(value):
                        rows.append(
                            {
                                "metric": str(value_col) if str(value_col) != "0" else str(metric_name),
                                "epoch": row[epoch_col],
                                "value": value,
                            }
                        )
            return

        try:
            for epoch, value in enumerate(values, start=1):
                if pd.notna(value):
                    rows.append({"metric": str(metric_name), "epoch": epoch, "value": value})
        except TypeError:
            if pd.notna(values):
                rows.append({"metric": str(metric_name), "epoch": 1, "value": values})

    if isinstance(history, dict):
        for metric_name, values in history.items():
            add_values(metric_name, values)
    else:
        add_values("training_history", history)

    out = pd.DataFrame(rows, columns=["metric", "epoch", "value"])
    if out.empty:
        return out

    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    out["epoch"] = pd.to_numeric(out["epoch"], errors="coerce")
    out = out.dropna(subset=["epoch", "value"])
    out["epoch"] = out["epoch"].astype(int) + (0 if out["epoch"].min() == 1 else 1)
    return out.sort_values(["metric", "epoch"]).reset_index(drop=True)


def save_convergence_outputs(model, pd, history_tsv: Path | None, convergence_png: Path | None, heuristic_epochs: int, train_epochs: int) -> None:
    history_df = training_history_to_dataframe(model, pd)

    if history_tsv is not None:
        history_tsv.parent.mkdir(parents=True, exist_ok=True)
        history_df.to_csv(history_tsv, sep="\t", index=False)
        print(f"Wrote training history TSV: {history_tsv}", flush=True)

    if convergence_png is None:
        return

    convergence_png.parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise SystemExit(f"matplotlib is required to write convergence PNG: {e}") from e

    if history_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6), dpi=180)
        ax.text(
            0.5,
            0.5,
            "No PeakVI training history was found in model.history",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.axvline(heuristic_epochs, color="black", linestyle="--", linewidth=1.2, label=f"heuristic t = {heuristic_epochs}")
        ax.axvline(train_epochs, color="gray", linestyle=":", linewidth=1.1, label=f"training ceiling = {train_epochs}")
        ax.set_xlabel("Epoch")
        ax.set_title("PeakVI convergence check")
        ax.legend(loc="best", fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(convergence_png)
        plt.close(fig)
        print(f"Wrote placeholder convergence PNG: {convergence_png}", flush=True)
        return

    metric_order = list(dict.fromkeys(history_df["metric"].astype(str)))
    loss_like = [
        metric
        for metric in metric_order
        if any(token in metric.lower() for token in ("loss", "elbo", "reconstruction", "kl", "validation"))
    ]
    plot_metrics = loss_like if loss_like else metric_order

    fig, ax = plt.subplots(figsize=(10, 6), dpi=180)
    for metric in plot_metrics:
        df = history_df[history_df["metric"].astype(str) == metric]
        ax.plot(df["epoch"], df["value"], linewidth=1.6, label=metric)

    ax.axvline(heuristic_epochs, color="black", linestyle="--", linewidth=1.2, label=f"heuristic t = {heuristic_epochs}")
    ax.axvline(train_epochs, color="gray", linestyle=":", linewidth=1.1, label=f"training ceiling = {train_epochs}")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Training history value")
    ax.set_title("PeakVI convergence check")
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.legend(loc="best", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(convergence_png)
    plt.close(fig)
    print(f"Wrote convergence PNG: {convergence_png}", flush=True)


def main() -> None:
    args = parse_args()
    ad, np, pd, sp, torch, scvi = import_dependencies()

    input_h5ad = Path(args.input_h5ad)
    output_dir = Path(args.output_dir)
    latent_tsv = Path(args.latent_tsv)
    output_h5ad = Path(args.output_h5ad) if args.output_h5ad else None
    history_tsv = Path(args.history_tsv) if args.history_tsv else None
    convergence_png = Path(args.convergence_png) if args.convergence_png else None
    model_dir = output_dir / "model"

    output_dir.mkdir(parents=True, exist_ok=True)
    latent_tsv.parent.mkdir(parents=True, exist_ok=True)
    if output_h5ad is not None:
        output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    print(f"scvi-tools version: {scvi.__version__}", flush=True)
    print(f"torch version: {torch.__version__}", flush=True)
    print(f"CUDA available: {torch.cuda.is_available()}", flush=True)
    if torch.cuda.is_available():
        print(f"CUDA device count: {torch.cuda.device_count()}", flush=True)
        print(f"CUDA device 0: {torch.cuda.get_device_name(0)}", flush=True)

    scvi.settings.seed = args.random_seed
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("high")
    accelerator = choose_accelerator(args.accelerator, torch)

    print(f"Reading AnnData: {input_h5ad}", flush=True)
    adata = ad.read_h5ad(input_h5ad)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if args.batch_key is not None and args.batch_key not in adata.obs.columns:
        raise SystemExit(f"batch_key '{args.batch_key}' was not found in adata.obs")

    if args.binarize:
        print("Binarizing adata.X in memory for PeakVI: counts > 0 become 1", flush=True)
        binarize_matrix_in_place(adata, np, sp)

    heuristic_max_epochs = resolve_max_epochs(
        args.max_epochs,
        n_obs=adata.n_obs,
        epochs_cap=args.epochs_cap,
        decay_at_n_obs=args.decay_at_n_obs,
    )
    if args.epoch_multiplier < 1:
        raise SystemExit("--epoch_multiplier must be >= 1")
    training_max_epochs = int(heuristic_max_epochs * args.epoch_multiplier)

    print(f"AnnData shape, cells x peaks: {adata.n_obs:,} x {adata.n_vars:,}", flush=True)
    print(f"Training accelerator: {accelerator}", flush=True)
    print(f"Latent dimensions: {args.n_latent}", flush=True)
    print(f"Requested max epochs: {args.max_epochs}", flush=True)
    print(f"Heuristic/explicit epoch count t: {heuristic_max_epochs}", flush=True)
    print(f"Epoch multiplier: {args.epoch_multiplier}", flush=True)
    print(f"Precomputed max epochs passed to train(): {training_max_epochs}", flush=True)
    print(f"Auto epoch settings: epochs_cap={args.epochs_cap}; decay_at_n_obs={args.decay_at_n_obs}", flush=True)
    print(f"Batch size: {args.batch_size}", flush=True)
    print(f"Early stopping: {args.early_stopping}; patience: {args.early_stopping_patience}", flush=True)
    print(f"Dataloader workers: {args.dataloader_workers}", flush=True)

    scvi.model.PEAKVI.setup_anndata(adata, batch_key=args.batch_key)
    model = scvi.model.PEAKVI(
        adata,
        n_latent=args.n_latent,
        n_hidden=args.n_hidden,
        n_layers_encoder=args.n_layers_encoder,
        n_layers_decoder=args.n_layers_decoder,
        dropout_rate=args.dropout_rate,
        latent_distribution=args.latent_distribution,
    )

    train_peakvi_with_compat(model, args, accelerator, training_max_epochs)
    save_convergence_outputs(model, pd, history_tsv, convergence_png, heuristic_max_epochs, training_max_epochs)

    print("Extracting latent representation", flush=True)
    latent = model.get_latent_representation()
    latent = np.asarray(latent, dtype=np.float32)
    latent_cols = [f"PEAKVI_{i}" for i in range(1, latent.shape[1] + 1)]
    latent_df = pd.DataFrame(latent, index=adata.obs_names, columns=latent_cols)
    latent_df.index.name = "cell_name"
    latent_df.to_csv(latent_tsv, sep="\t")

    np.save(output_dir / "peakvi_latent.npy", latent)
    pd.Series(adata.obs_names, name="cell_name").to_csv(output_dir / "cell_names.txt", index=False, header=False)

    print(f"Saving model to: {model_dir}", flush=True)
    model.save(model_dir, overwrite=True, save_anndata=False)

    adata.obsm["X_peakvi"] = latent
    if output_h5ad is not None:
        print(f"Writing AnnData with X_peakvi: {output_h5ad}", flush=True)
        adata.write_h5ad(output_h5ad)

    summary = {
        "atlas": args.atlas,
        "input_h5ad": str(input_h5ad),
        "n_cells": int(adata.n_obs),
        "n_peaks": int(adata.n_vars),
        "n_latent": int(args.n_latent),
        "n_hidden": int(args.n_hidden),
        "n_layers_encoder": int(args.n_layers_encoder),
        "n_layers_decoder": int(args.n_layers_decoder),
        "dropout_rate": float(args.dropout_rate),
        "latent_distribution": args.latent_distribution,
        "requested_max_epochs": str(args.max_epochs),
        "heuristic_epochs": int(heuristic_max_epochs),
        "epoch_multiplier": int(args.epoch_multiplier),
        "max_epochs": int(training_max_epochs),
        "epochs_cap": int(args.epochs_cap),
        "decay_at_n_obs": int(args.decay_at_n_obs),
        "batch_size": int(args.batch_size),
        "batch_key": args.batch_key,
        "binarize": bool(args.binarize),
        "early_stopping": bool(args.early_stopping),
        "early_stopping_patience": int(args.early_stopping_patience),
        "dataloader_workers": int(args.dataloader_workers),
        "accelerator": accelerator,
        "devices": str(args.devices),
        "random_seed": int(args.random_seed),
        "latent_tsv": str(latent_tsv),
        "history_tsv": str(history_tsv) if history_tsv is not None else None,
        "convergence_png": str(convergence_png) if convergence_png is not None else None,
        "model_dir": str(model_dir),
        "output_h5ad": str(output_h5ad) if output_h5ad is not None else None,
        "scvi_version": scvi.__version__,
        "torch_version": torch.__version__,
    }

    with open(output_dir / "peakvi_run_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("DONE", flush=True)


if __name__ == "__main__":
    main()
