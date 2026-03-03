
#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


BARCODE_RE = re.compile(r"([ACGT]{16}-\d+)")


def extract_barcode(s: str) -> str:
    if s is None:
        return ""
    m = BARCODE_RE.search(str(s))
    return m.group(1) if m else ""


def shannon_entropy_from_counts(counts: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return np.nan
    p = counts / total
    p = p[p > 0]
    H = float(-(p * np.log2(p)).sum())
    # remove negative zero artefact
    if abs(H) < 1e-12:
        H = 0.0
    return H


def spearman_corr(x: pd.Series, y: pd.Series) -> float:
    xr = x.rank(method="average")
    yr = y.rank(method="average")
    if xr.nunique() < 2 or yr.nunique() < 2:
        return np.nan
    return float(xr.corr(yr))


def pearson_corr(x: pd.Series, y: pd.Series) -> float:
    if x.nunique() < 2 or y.nunique() < 2:
        return np.nan
    return float(x.corr(y))


def pick_first_existing(df: pd.DataFrame, candidates) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--pairs_tsv", required=True)
    ap.add_argument("--clones_tsv", required=True)
    ap.add_argument("--out_per_clone_tsv", required=True)
    ap.add_argument("--out_summary_tsv", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs_tsv, sep="\t", dtype=str).fillna("")
    clones = pd.read_csv(args.clones_tsv, sep="\t", dtype=str).fillna("")

    # barcode in pairs
    pairs_barcode_col = pick_first_existing(pairs, ["barcode", "cell_barcode"])
    if not pairs_barcode_col:
        for c in pairs.columns:
            bc = pairs[c].map(extract_barcode)
            if (bc != "").any():
                pairs["barcode"] = bc
                pairs_barcode_col = "barcode"
                break

    if not pairs_barcode_col:
        raise SystemExit("No barcode column found in pairs file.")

    pairs = pairs.rename(columns={pairs_barcode_col: "barcode"}).copy()

    # clone_id + sequence_id in clones
    clone_id_col = pick_first_existing(clones, ["clone_id", "CLONE_ID"])
    seq_id_col = pick_first_existing(clones, ["sequence_id", "SEQUENCE_ID"])

    if not clone_id_col or not seq_id_col:
        raise SystemExit("Missing clone_id or sequence_id in clones file.")

    clones = clones.rename(columns={clone_id_col: "clone_id", seq_id_col: "sequence_id"}).copy()
    clones["barcode"] = clones["sequence_id"].map(extract_barcode)

    # join
    joined = pairs.merge(
        clones[["barcode", "clone_id"]],
        on="barcode",
        how="inner"
    )

    if joined.empty:
        raise SystemExit("Join failed: no overlapping barcodes.")

    #  clone sizes
    clone_sizes = (
        clones.groupby("clone_id")["barcode"]
        .nunique()
        .rename("clone_size_barcodes")
        .reset_index()
    )

    # entropy per clone
    joined["pair"] = joined["v_gene_H"] + "|" + joined["v_gene_L"]

    rows = []
    for cid, g in joined.groupby("clone_id"):
        vc = g["pair"].value_counts()
        H = shannon_entropy_from_counts(vc.values)
        rows.append({
            "sample": args.sample,
            "clone_id": cid,
            "pairing_entropy_bits": H,
            "n_distinct_vhvl_pairs": int(vc.shape[0]),
        })

    per_clone = pd.DataFrame(rows)
    per_clone = per_clone.merge(clone_sizes, on="clone_id", how="left")

    # write per-clone
    Path(args.out_per_clone_tsv).parent.mkdir(parents=True, exist_ok=True)
    per_clone.to_csv(args.out_per_clone_tsv, sep="\t", index=False)

    # correlation (filter expanded clones)
    df_corr = per_clone[
        per_clone["clone_size_barcodes"] >= args.min_clone_barcodes
    ].copy()

    if df_corr.shape[0] >= 3:
        size = df_corr["clone_size_barcodes"].astype(float)
        entropy = df_corr["pairing_entropy_bits"].astype(float)

        rho = spearman_corr(size, entropy)

        log_size = np.log10(size)
        rho_log = spearman_corr(log_size, entropy)
        r_log = pearson_corr(log_size, entropy)
    else:
        rho = rho_log = r_log = np.nan

    summary = pd.DataFrame([{
        "sample": args.sample,
        "n_clones_for_corr": int(df_corr.shape[0]),
        "spearman_rho_size_vs_entropy": rho,
        "spearman_rho_log10size_vs_entropy": rho_log,
        "pearson_r_log10size_vs_entropy": r_log,
        "min_clone_barcodes": args.min_clone_barcodes
    }])

    Path(args.out_summary_tsv).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out_summary_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()