#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def corr_safe(x, y, method: str) -> float:
    s = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if s.shape[0] < 3:
        return float("nan")
    if s["x"].nunique() < 2 or s["y"].nunique() < 2:
        return float("nan")
    return float(s["x"].corr(s["y"], method=method))


def main():
    ap = argparse.ArgumentParser(description="Master dominance/entropy correlations across samples")
    ap.add_argument("--samples", nargs="+", required=True, help="Sample names")
    ap.add_argument("--metrics_dir", required=True, help="Directory containing per-sample dominance TSVs")
    ap.add_argument("--min_clone_barcodes", type=int, default=2, help="Filter clones by n_barcodes_with_pairs")
    ap.add_argument("--out_tsv", required=True, help="Output TSV path")
    args = ap.parse_args()

    metrics_dir = Path(args.metrics_dir)
    rows = []

    for sample in args.samples:
        tsv = metrics_dir / f"{sample}.clone_pairing_dominance.tsv"
        df = pd.read_csv(tsv, sep="\t")

        if "n_barcodes_with_pairs" not in df.columns:
            raise SystemExit(f"Missing n_barcodes_with_pairs in {tsv}")
        df = df[df["n_barcodes_with_pairs"] >= args.min_clone_barcodes].copy()

        n_clones = int(df.shape[0])
        med_size_seq = float(np.median(df["clone_size_sequences"])) if n_clones else float("nan")
        med_entropy = float(np.median(df["pairing_entropy_bits"])) if n_clones else float("nan")
        med_dom = float(np.median(df["dominant_pair_frac"])) if n_clones else float("nan")

        if n_clones:
            log10size = np.log10(df["clone_size_sequences"].astype(float).replace(0, np.nan))
        else:
            log10size = pd.Series([], dtype=float)

        rho_log10_dom = corr_safe(log10size, df.get("dominant_pair_frac", pd.Series([], dtype=float)).astype(float), "spearman")
        r_log10_dom   = corr_safe(log10size, df.get("dominant_pair_frac", pd.Series([], dtype=float)).astype(float), "pearson")

        rho_log10_ent = corr_safe(log10size, df.get("pairing_entropy_bits", pd.Series([], dtype=float)).astype(float), "spearman")
        r_log10_ent   = corr_safe(log10size, df.get("pairing_entropy_bits", pd.Series([], dtype=float)).astype(float), "pearson")

        rho_dom_ent   = corr_safe(df.get("dominant_pair_frac", pd.Series([], dtype=float)).astype(float),
                                 df.get("pairing_entropy_bits", pd.Series([], dtype=float)).astype(float),
                                 "spearman")
        r_dom_ent     = corr_safe(df.get("dominant_pair_frac", pd.Series([], dtype=float)).astype(float),
                                 df.get("pairing_entropy_bits", pd.Series([], dtype=float)).astype(float),
                                 "pearson")

        rows.append({
            "sample": sample,
            "n_clones": n_clones,
            "min_clone_barcodes": int(args.min_clone_barcodes),
            "median_clone_size_sequences": med_size_seq,
            "median_entropy_bits": med_entropy,
            "median_dominant_pair_frac": med_dom,
            "spearman_rho_log10size_vs_dominance": rho_log10_dom,
            "pearson_r_log10size_vs_dominance": r_log10_dom,
            "spearman_rho_log10size_vs_entropy": rho_log10_ent,
            "pearson_r_log10size_vs_entropy": r_log10_ent,
            "spearman_rho_dominance_vs_entropy": rho_dom_ent,
            "pearson_r_dominance_vs_entropy": r_dom_ent
        })

    out = pd.DataFrame(rows)
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    print(f"Wrote: {out.shape[0]} rows -> {out_path}")


if __name__ == "__main__":
    main()
