#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


def load_data(sample, metrics_dir, min_clone_barcodes):
    tsv = Path(metrics_dir) / f"{sample}.pairing_entropy_per_clone.tsv"
    df = pd.read_csv(tsv, sep="\t")

    df = df[df["clone_size_barcodes"] >= min_clone_barcodes].copy()
    df["log10_size"] = np.log10(df["clone_size_barcodes"])

    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2)
    ap.add_argument("--dpi", type=int, default=300)

    args = ap.parse_args()

    n = len(args.samples)
    nrows = 2
    ncols = 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 9))
    axes = axes.flatten()

    for i, sample in enumerate(args.samples):
        df = load_data(sample, args.metrics_dir, args.min_clone_barcodes)
        ax = axes[i]

        rho, p = spearmanr(
            df["clone_size_barcodes"],
            df["n_distinct_vhvl_pairs"]
        )

        ax.scatter(
            df["clone_size_barcodes"],
            df["n_distinct_vhvl_pairs"],
            alpha=0.6
        )

        ax.set_xscale("log")
        ax.set_title(sample)
        ax.set_xlabel("Clone size (barcodes, log scale)")
        ax.set_ylabel("n distinct VH–VL pairs")

        ax.text(
            0.05,
            0.95,
            f"n={len(df)}\nSpearman: {rho:.3f}",
            transform=ax.transAxes,
            verticalalignment="top"
        )

    fig.suptitle("Clone size vs distinct VH–VL pair count", fontsize=14)

    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=args.dpi)
    plt.close(fig)


if __name__ == "__main__":
    main()
