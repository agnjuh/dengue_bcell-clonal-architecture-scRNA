#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_summary(summary_path: Path) -> dict:
    if not summary_path.exists():
        return {}
    df = pd.read_csv(summary_path, sep="\t")
    if df.empty:
        return {}
    row = df.iloc[0]
    out = {}
    for c in df.columns:
        try:
            out[c] = float(row[c])
        except Exception:
            out[c] = row[c]
    return out


def load_per_clone(tsv: Path, sample: str, min_clone_barcodes: int):
    df = pd.read_csv(tsv, sep="\t")
    if df.empty:
        return df

    df["pairing_entropy_bits"] = pd.to_numeric(df["pairing_entropy_bits"], errors="coerce")
    df["clone_size_barcodes"] = pd.to_numeric(df["clone_size_barcodes"], errors="coerce")

    df = df.dropna(subset=["pairing_entropy_bits", "clone_size_barcodes"])
    df = df[df["clone_size_barcodes"] >= min_clone_barcodes].copy()
    df["sample"] = sample
    return df


def panel(ax, df, sample, stats, point_size, alpha):
    if df.empty:
        ax.set_title(sample)
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    x = df["clone_size_barcodes"].to_numpy()
    y = df["pairing_entropy_bits"].to_numpy()

    ax.scatter(x, y, s=point_size, alpha=alpha)
    ax.set_xscale("log")
    ax.set_title(sample)

    s1 = stats.get("spearman_rho_size_vs_entropy", np.nan)
    s2 = stats.get("spearman_rho_log10size_vs_entropy", np.nan)
    p1 = stats.get("pearson_r_log10size_vs_entropy", np.nan)

    def fmt(v):
        if isinstance(v, float) and np.isnan(v):
            return "NA"
        return f"{v:.3f}"

    text = (
        f"n={len(df)}\n"
        f"Spearman(size): {fmt(s1)}\n"
        f"Spearman(log10): {fmt(s2)}\n"
        f"Pearson(log10): {fmt(p1)}"
    )

    ax.text(
        0.02, 0.98, text,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.75)
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2)
    ap.add_argument("--point_size", type=float, default=10.0)
    ap.add_argument("--alpha", type=float, default=0.45)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    metrics_dir = Path(args.metrics_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_sample = {}
    stats_sample = {}

    for s in args.samples:
        per_clone_tsv = metrics_dir / f"{s}.pairing_entropy_per_clone.tsv"
        summary_tsv = metrics_dir / f"{s}.pairing_entropy_vs_clone_size.summary.tsv"

        per_sample[s] = load_per_clone(per_clone_tsv, s, args.min_clone_barcodes)
        stats_sample[s] = read_summary(summary_tsv)

    n = len(args.samples)
    if n == 1:
        nrows, ncols = 1, 1
    elif n == 2:
        nrows, ncols = 1, 2
    else:
        nrows, ncols = 2, 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))
    axes = np.array(axes).reshape(-1)

    for i in range(nrows * ncols):
        ax = axes[i]
        if i >= n:
            ax.axis("off")
            continue
        s = args.samples[i]
        panel(ax, per_sample[s], s, stats_sample[s], args.point_size, args.alpha)

    fig.supxlabel("Clone size (barcodes, log scale)")
    fig.supylabel("Pairing entropy (bits)")
    fig.suptitle("Pairing entropy vs clonotype expansion", y=0.98)

    fig.savefig(out_dir / "ALL.pairing_entropy_vs_clone_size.2x2.png",
                dpi=args.dpi, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    # Combined colored plot
    fig, ax = plt.subplots(figsize=(7, 5))

    for s in args.samples:
        df = per_sample[s]
        if df.empty:
            continue
        ax.scatter(
            df["clone_size_barcodes"],
            df["pairing_entropy_bits"],
            s=args.point_size,
            alpha=args.alpha,
            label=s
        )

    ax.set_xscale("log")
    ax.legend(frameon=False)
    fig.supxlabel("Clone size (barcodes, log scale)")
    fig.supylabel("Pairing entropy (bits)")
    fig.suptitle("Pairing entropy vs clonotype expansion (ALL)", y=0.98)

    fig.savefig(out_dir / "ALL.pairing_entropy_vs_clone_size.scatter.png",
                dpi=args.dpi, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


if __name__ == "__main__":
    main()
