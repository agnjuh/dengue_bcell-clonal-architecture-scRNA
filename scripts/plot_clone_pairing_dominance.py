#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REQUIRED_COLS = [
    "sample",
    "clone_id",
    "n_barcodes_with_pairs",
    "dominant_pair_frac",
    "pairing_entropy_bits",
    "clone_size_sequences",
]


def load_sample_tsv(tsv: Path, sample: str, min_clone_barcodes: int) -> pd.DataFrame:
    if not tsv.exists():
        raise SystemExit(f"Missing file: {tsv}")

    df = pd.read_csv(tsv, sep="\t")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise SystemExit(f"{tsv}: missing columns: {missing}")

    df = df.copy()
    df["sample"] = sample

    # numeric
    for c in ["n_barcodes_with_pairs", "dominant_pair_frac", "pairing_entropy_bits", "clone_size_sequences"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["n_barcodes_with_pairs", "dominant_pair_frac", "pairing_entropy_bits", "clone_size_sequences"])
    df = df[df["n_barcodes_with_pairs"] >= int(min_clone_barcodes)].copy()

    # log10 size (avoid -inf)
    df = df[df["clone_size_sequences"] > 0].copy()
    df["log10_clone_size_sequences"] = np.log10(df["clone_size_sequences"].astype(float))

    # clip dominance to [0,1] just in case
    df["dominant_pair_frac"] = df["dominant_pair_frac"].clip(lower=0.0, upper=1.0)

    return df


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    # Spearman via ranking then Pearson on ranks (no scipy dependency)
    rx = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    ry = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    if rx.size < 3:
        return np.nan
    return pearson_r(rx, ry)


def pearson_r(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 3:
        return np.nan
    x = x.astype(float)
    y = y.astype(float)
    x = x - np.mean(x)
    y = y - np.mean(y)
    denom = (np.sqrt(np.sum(x * x)) * np.sqrt(np.sum(y * y)))
    if denom == 0:
        return np.nan
    return float(np.sum(x * y) / denom)


def compute_correlations(df: pd.DataFrame, sample: str) -> dict:
    x = df["log10_clone_size_sequences"].to_numpy(dtype=float)
    dom = df["dominant_pair_frac"].to_numpy(dtype=float)
    ent = df["pairing_entropy_bits"].to_numpy(dtype=float)

    out = {
        "sample": sample,
        "n_clones": int(df.shape[0]),
        "min_clone_barcodes": int(df["n_barcodes_with_pairs"].min()) if df.shape[0] else np.nan,
        "median_clone_size_sequences": float(df["clone_size_sequences"].median()) if df.shape[0] else np.nan,
        "median_entropy_bits": float(df["pairing_entropy_bits"].median()) if df.shape[0] else np.nan,
        "median_dominant_pair_frac": float(df["dominant_pair_frac"].median()) if df.shape[0] else np.nan,
        # size vs dominance
        "spearman_rho_log10size_vs_dominance": spearman_rho(x, dom),
        "pearson_r_log10size_vs_dominance": pearson_r(x, dom),
        # size vs entropy
        "spearman_rho_log10size_vs_entropy": spearman_rho(x, ent),
        "pearson_r_log10size_vs_entropy": pearson_r(x, ent),
        # dominance vs entropy
        "spearman_rho_dominance_vs_entropy": spearman_rho(dom, ent),
        "pearson_r_dominance_vs_entropy": pearson_r(dom, ent),
    }
    return out


def _grid(n: int):
    if n <= 1:
        return 1, 1
    if n == 2:
        return 1, 2
    return 2, 2


def _outer_labels(ax, r, c, nrows, ncols, xlabel, ylabel):
    # only bottom row gets x tick labels + xlabel
    if r == nrows - 1:
        ax.set_xlabel(xlabel)
        ax.tick_params(labelbottom=True)
    else:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)

    # only left col gets y tick labels + ylabel
    if c == 0:
        ax.set_ylabel(ylabel)
        ax.tick_params(labelleft=True)
    else:
        ax.set_ylabel("")
        ax.tick_params(labelleft=False)


def plot_2x2(per_sample: dict, samples: list, out_png: Path, kind: str, dpi: int):
    n = len(samples)
    nrows, ncols = _grid(n)

    # reasonable sizing (avoids giant whitespace)
    fig_w = 11 if (nrows, ncols) == (2, 2) else 10
    fig_h = 8.5 if (nrows, ncols) == (2, 2) else 4.8
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h))
    if hasattr(axes, "flatten"):
        axes = axes.flatten()
    else:
        axes = [axes]

    # determine global limits for comparability
    all_x = np.concatenate([per_sample[s]["log10_clone_size_sequences"].to_numpy(float) for s in samples])
    x_min, x_max = float(np.nanmin(all_x)), float(np.nanmax(all_x))
    x_pad = 0.05 * (x_max - x_min) if x_max > x_min else 0.2
    xlim = (x_min - x_pad, x_max + x_pad)

    if kind == "dominance_vs_logsize":
        all_y = np.concatenate([per_sample[s]["dominant_pair_frac"].to_numpy(float) for s in samples])
        ylim = (-0.02, 1.02)
        xlabel = "log10(clone size; sequences)"
        ylabel = "dominant VH–VL pair fraction"
        title = "Clone pairing dominance vs clonotype expansion"
    elif kind == "entropy_vs_logsize":
        all_y = np.concatenate([per_sample[s]["pairing_entropy_bits"].to_numpy(float) for s in samples])
        y_min, y_max = float(np.nanmin(all_y)), float(np.nanmax(all_y))
        y_pad = 0.08 * (y_max - y_min) if y_max > y_min else 0.2
        ylim = (max(0.0, y_min - y_pad), y_max + y_pad)
        xlabel = "log10(clone size; sequences)"
        ylabel = "pairing entropy (bits)"
        title = "Pairing entropy vs clonotype expansion"
    elif kind == "dominance_vs_entropy":
        all_y = np.concatenate([per_sample[s]["dominant_pair_frac"].to_numpy(float) for s in samples])
        ylim = (-0.02, 1.02)
        xlabel = "pairing entropy (bits)"
        ylabel = "dominant VH–VL pair fraction"
        title = "Dominance vs entropy"
    else:
        raise ValueError(f"Unknown kind: {kind}")

    for i, s in enumerate(samples):
        ax = axes[i]
        df = per_sample[s]

        r = i // ncols
        c = i % ncols

        if kind == "dominance_vs_logsize":
            ax.scatter(df["log10_clone_size_sequences"], df["dominant_pair_frac"], s=16, alpha=0.65, linewidths=0)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            _outer_labels(ax, r, c, nrows, ncols, xlabel, ylabel)

        elif kind == "entropy_vs_logsize":
            ax.scatter(df["log10_clone_size_sequences"], df["pairing_entropy_bits"], s=16, alpha=0.65, linewidths=0)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            _outer_labels(ax, r, c, nrows, ncols, xlabel, ylabel)

        elif kind == "dominance_vs_entropy":
            ax.scatter(df["pairing_entropy_bits"], df["dominant_pair_frac"], s=16, alpha=0.65, linewidths=0)
            ax.set_ylim(*ylim)
            # xlim from data
            x2 = df["pairing_entropy_bits"].to_numpy(float)
            if x2.size:
                xmin, xmax = float(np.nanmin(x2)), float(np.nanmax(x2))
                pad = 0.08 * (xmax - xmin) if xmax > xmin else 0.2
                ax.set_xlim(xmin - pad, xmax + pad)
            _outer_labels(ax, r, c, nrows, ncols, xlabel, ylabel)

        ax.set_title(f"{s} (n={df.shape[0]})", fontsize=12)
        ax.grid(False)

    # turn off unused panels
    for j in range(n, nrows * ncols):
        axes[j].axis("off")

    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    metrics_dir = Path(args.metrics_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # load
    per_sample = {}
    corr_rows = []

    for s in args.samples:
        tsv = metrics_dir / f"{s}.clone_pairing_dominance.tsv"
        df = load_sample_tsv(tsv, s, args.min_clone_barcodes)
        if df.empty:
            raise SystemExit(f"No rows after filtering for sample {s} (min_clone_barcodes={args.min_clone_barcodes})")
        per_sample[s] = df
        corr_rows.append(compute_correlations(df, s))

    corr_df = pd.DataFrame(corr_rows)
    corr_out = out_dir / "clone_pairing_dominance.correlations.tsv"
    corr_df.to_csv(corr_out, sep="\t", index=False)

    # plots
    plot_2x2(
        per_sample,
        args.samples,
        out_dir / "dominance_vs_log10_clone_size.png",
        kind="dominance_vs_logsize",
        dpi=args.dpi,
    )
    plot_2x2(
        per_sample,
        args.samples,
        out_dir / "entropy_vs_log10_clone_size.png",
        kind="entropy_vs_logsize",
        dpi=args.dpi,
    )
    plot_2x2(
        per_sample,
        args.samples,
        out_dir / "dominance_vs_entropy.png",
        kind="dominance_vs_entropy",
        dpi=args.dpi,
    )

    print(f"Wrote: {corr_out}")
    print(f"Wrote: {out_dir / 'dominance_vs_log10_clone_size.png'}")
    print(f"Wrote: {out_dir / 'entropy_vs_log10_clone_size.png'}")
    print(f"Wrote: {out_dir / 'dominance_vs_entropy.png'}")


if __name__ == "__main__":
    main()
