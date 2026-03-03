#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_data(sample, metrics_dir, min_clone_barcodes):
    tsv = Path(metrics_dir) / f"{sample}.pairing_entropy_per_clone.tsv"
    df = pd.read_csv(tsv, sep="\t")

    df = df[df["clone_size_barcodes"] >= min_clone_barcodes].copy()
    df["nonzero_entropy"] = df["pairing_entropy_bits"] > 0.0

    return df


def compute_probability_curve(df):
    grouped = (
        df.groupby("clone_size_barcodes")
        .agg(
            n_total=("nonzero_entropy", "size"),
            n_nonzero=("nonzero_entropy", "sum")
        )
        .reset_index()
    )

    grouped["prob_nonzero"] = grouped["n_nonzero"] / grouped["n_total"]

    return grouped.sort_values("clone_size_barcodes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2)
    ap.add_argument("--dpi", type=int, default=300)

    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(8, 6))

    for sample in args.samples:
        df = load_data(sample, args.metrics_dir, args.min_clone_barcodes)
        curve = compute_probability_curve(df)

        ax.plot(
            curve["clone_size_barcodes"],
            curve["prob_nonzero"],
            marker="o",
            label=sample
        )

    ax.set_xscale("log")
    ax.set_xlabel("Clone size (barcodes, log scale)")
    ax.set_ylabel("P(entropy > 0)")
    ax.set_ylim(0, 1)
    ax.set_title("Probability of non-zero pairing entropy vs clone size")

    ax.legend()

    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=args.dpi)
    plt.close(fig)


if __name__ == "__main__":
    main()
