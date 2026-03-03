#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_data(sample: str, metrics_dir: str, min_clone_barcodes: int) -> pd.DataFrame:
    tsv = Path(metrics_dir) / f"{sample}.pairing_entropy_per_clone.tsv"
    if not tsv.exists():
        raise SystemExit(f"Missing input: {tsv}")

    df = pd.read_csv(tsv, sep="\t")
    req = {"clone_size_barcodes", "pairing_entropy_bits"}
    missing = req - set(df.columns)
    if missing:
        raise SystemExit(f"{tsv} missing columns: {sorted(missing)}")

    df = df[df["clone_size_barcodes"] >= min_clone_barcodes].copy()
    df["nonzero_entropy"] = df["pairing_entropy_bits"].astype(float) > 0.0
    df["clone_size_barcodes"] = df["clone_size_barcodes"].astype(int)

    return df


def assign_quantile_bins(df: pd.DataFrame, n_bins: int = 3) -> pd.Series:
    """
    Robust quantile binning for discrete heavy-tied sizes.

    Strategy:
      1) Try qcut on size with duplicates='drop'
      2) If bins collapse (<2 bins), do rank-based qcut (breaks ties deterministically)
    Returns a categorical Series (ordered) with labels like Q1/Q2/Q3 (or fewer if needed).
    """
    s = df["clone_size_barcodes"].astype(float)

    # Attempt 1: qcut directly (may fail on tied edges)
    try:
        b = pd.qcut(s, q=n_bins, duplicates="drop")
        # if everything collapsed into a single bin, fallback
        if b.cat.categories.size < 2:
            raise ValueError("qcut collapsed bins")
        return b
    except Exception:
        pass

    r = s.rank(method="average")
    jitter = (np.arange(len(r)) + 1) * 1e-9
    rj = r.to_numpy() + jitter

    b = pd.qcut(pd.Series(rj, index=df.index), q=n_bins, duplicates="drop")
    if b.cat.categories.size < 2:
        # absolute fallback: single bin
        b = pd.Series(["All"] * len(df), index=df.index, dtype="category")
    return b


def summarize_by_bin(df: pd.DataFrame, bins: pd.Series) -> pd.DataFrame:
    tmp = df.copy()
    tmp["size_bin"] = bins

    grouped = (
        tmp.groupby("size_bin", observed=True)
        .agg(
            n_total=("nonzero_entropy", "size"),
            n_nonzero=("nonzero_entropy", "sum"),
            median_size=("clone_size_barcodes", "median"),
            min_size=("clone_size_barcodes", "min"),
            max_size=("clone_size_barcodes", "max"),
        )
        .reset_index()
    )
    grouped["prob_nonzero"] = grouped["n_nonzero"] / grouped["n_total"]
    return grouped


def nice_bin_labels(grouped: pd.DataFrame) -> list[str]:
    # Label bins as e.g. "Q1 (2–2; med=2)"
    out = []
    for i, r in grouped.iterrows():
        bname = str(r["size_bin"])
        # If bin names are intervals, give Q1/Q2... instead
        q = f"Q{i+1}"
        lo = int(r["min_size"])
        hi = int(r["max_size"])
        med = float(r["median_size"])
        out.append(f"{q} ({lo}–{hi}; med={med:g})")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2)
    ap.add_argument("--n_bins", type=int, default=3)
    ap.add_argument("--dpi", type=int, default=300)

    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(9, 6))

   
    per_sample_summary = {}
    max_bins = 0

    for sample in args.samples:
        df = load_data(sample, args.metrics_dir, args.min_clone_barcodes)
        bins = assign_quantile_bins(df, n_bins=args.n_bins)
        summ = summarize_by_bin(df, bins)
        per_sample_summary[sample] = summ
        max_bins = max(max_bins, len(summ))

    # Plot: x positions = 0..(k-1) per sample
    for sample in args.samples:
        summ = per_sample_summary[sample].copy()
        x = np.arange(len(summ))
        ax.plot(x, summ["prob_nonzero"], marker="o", label=sample)

    ax.set_ylim(0, 1)
    ax.set_ylabel("P(entropy > 0)")
    ax.set_title("Probability of non-zero pairing entropy\n(quantile-binned clone size)")

    # Use the first sample’s bin labels if all samples have same number of bins,
    # otherwise keep generic Q1/Q2/Q3...
    same_k = len({len(per_sample_summary[s]) for s in args.samples}) == 1
    if same_k and args.samples:
        labels = nice_bin_labels(per_sample_summary[args.samples[0]])
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=0)
    else:
        ax.set_xticks(np.arange(max_bins))
        ax.set_xticklabels([f"Q{i+1}" for i in range(max_bins)])

    ax.legend()

    out = Path(args.out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)


if __name__ == "__main__":
    main()
