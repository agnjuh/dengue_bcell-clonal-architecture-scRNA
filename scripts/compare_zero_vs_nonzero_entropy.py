#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def rank_biserial(u_stat, n1, n2):
    return 1 - (2 * u_stat) / (n1 * n2)


def analyze_sample(sample, metrics_dir, min_clone_barcodes):
    tsv = Path(metrics_dir) / f"{sample}.pairing_entropy_per_clone.tsv"
    df = pd.read_csv(tsv, sep="\t")

    df = df[df["clone_size_barcodes"] >= min_clone_barcodes].copy()

    if df.empty:
        return None

    df["log10_size"] = np.log10(df["clone_size_barcodes"])

    zero = df[df["pairing_entropy_bits"] == 0.0]
    nonzero = df[df["pairing_entropy_bits"] > 0.0]

    if len(zero) < 3 or len(nonzero) < 3:
        return None

    u, p = mannwhitneyu(
        zero["log10_size"],
        nonzero["log10_size"],
        alternative="two-sided"
    )

    rbc = rank_biserial(u, len(zero), len(nonzero))

    return {
        "sample": sample,
        "n_zero": len(zero),
        "n_nonzero": len(nonzero),
        "median_size_zero": zero["clone_size_barcodes"].median(),
        "median_size_nonzero": nonzero["clone_size_barcodes"].median(),
        "mannwhitney_p": p,
        "rank_biserial_effect_size": rbc
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2)

    args = ap.parse_args()

    rows = []
    for s in args.samples:
        res = analyze_sample(s, args.metrics_dir, args.min_clone_barcodes)
        if res:
            rows.append(res)

    if not rows:
        print("No valid samples.")
        return

    out = pd.DataFrame(rows)
    out.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
