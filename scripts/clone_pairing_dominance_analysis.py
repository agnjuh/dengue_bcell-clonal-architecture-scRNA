#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


BARCODE_RE = re.compile(r'([ACGT]+-\d+)')

def extract_barcode(s: str) -> str:
    """Try to extract 10x barcode like AACTG...-1 from a sequence_id header."""
    if not isinstance(s, str):
        return ""
    m = BARCODE_RE.search(s)
    return m.group(1) if m else ""

def extract_contig_id(s: str) -> str:
    """
    Try to extract a contig identifier from a sequence_id header.

    Supports common patterns:
      - contig_id=XYZ
      - ..._contig_1
      - ...|contig=XYZ
    """
    if not isinstance(s, str):
        return ""

    # contig_id=...
    m = re.search(r'(?:^|[|;\s])contig_id=([^\s|;]+)', s)
    if m:
        return m.group(1)

    # contig=...
    m = re.search(r'(?:^|[|;\s])contig=([^\s|;]+)', s)
    if m:
        return m.group(1)

    # barcode_contig_#
    m = re.search(r'([ACGT]+-\d+_contig_\d+)', s)
    if m:
        return m.group(1)

    # generic: something containing "contig" without spaces
    m = re.search(r'([^\s|;]*contig[^\s|;]*)', s)
    if m:
        return m.group(1)

    return ""

def shannon_entropy_bits(counts: np.ndarray) -> float:
    counts = counts.astype(float)
    tot = counts.sum()
    if tot <= 0:
        return 0.0
    p = counts / tot
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    h = -np.sum(p * np.log2(p))
    # avoid "-0.0"
    return float(0.0 if abs(h) < 1e-12 else h)

def main():
    ap = argparse.ArgumentParser(description="Clone-level VH–VL pairing dominance (requires join pairs->clones).")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--pairs_dir", required=True)
    ap.add_argument("--clones_dir", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--min_clone_barcodes", type=int, default=2, help="minimum paired barcodes per clone to report")

    args = ap.parse_args()

    pairs_path = Path(args.pairs_dir) / f"{args.sample}.vhvl_pairs.tsv"
    clones_path = Path(args.clones_dir) / f"{args.sample}.clones.tsv"
    out_path = Path(args.out_tsv)

    if not pairs_path.exists():
        raise SystemExit(f"Missing pairs file: {pairs_path}")
    if not clones_path.exists():
        raise SystemExit(f"Missing clones file: {clones_path}")

    pairs = pd.read_csv(pairs_path, sep="\t", dtype=str).fillna("")
    clones = pd.read_csv(clones_path, sep="\t", dtype=str).fillna("")

    # Validate expected cols
    for col in ["barcode", "contig_id_H", "v_gene_H", "v_gene_L"]:
        if col not in pairs.columns:
            raise SystemExit(f"Pairs file missing required column '{col}': {pairs_path}")

    for col in ["sequence_id", "clone_id"]:
        if col not in clones.columns:
            raise SystemExit(f"Clones file missing required column '{col}': {clones_path}")

    # Prefer only successfully paired rows if that column exists
    if "pairing_status" in pairs.columns:
        # keep only clearly paired if present; otherwise keep all
        paired_mask = pairs["pairing_status"].astype(str).str.lower().isin(["paired", "true", "t", "1", "yes"])
        if paired_mask.any():
            pairs = pairs.loc[paired_mask].copy()

    # Build mapping from clone sequences -> contig_id / barcode
    clmap = clones[["sequence_id", "clone_id"]].copy()
    clmap["barcode_from_seqid"] = clmap["sequence_id"].map(extract_barcode)
    clmap["contig_id_from_seqid"] = clmap["sequence_id"].map(extract_contig_id)

    # 1) Try contig-level join
    merged = None
    if (clmap["contig_id_from_seqid"] != "").any():
        m1 = pairs.merge(
            clmap,
            left_on="contig_id_H",
            right_on="contig_id_from_seqid",
            how="inner",
            suffixes=("", "_cl")
        )
        if len(m1) > 0:
            merged = m1

    # 2) Fallback: barcode-level join
    if merged is None:
        if (clmap["barcode_from_seqid"] == "").all():
            raise SystemExit(
                "Could not extract contig_id or barcode from clones.sequence_id, so cannot join pairs to clones.\n"
                f"Example sequence_id: {clones['sequence_id'].iloc[0] if len(clones) else 'NA'}"
            )
        m2 = pairs.merge(
            clmap,
            left_on="barcode",
            right_on="barcode_from_seqid",
            how="inner",
            suffixes=("", "_cl")
        )
        if len(m2) == 0:
            raise SystemExit(
                "Join produced 0 rows (tried contig_id_H→sequence_id-derived contig_id, then barcode→sequence_id-derived barcode)."
            )
        merged = m2

    # Clone size in sequences (from clones table)
    clone_size_sequences = clones.groupby("clone_id")["sequence_id"].nunique().rename("clone_size_sequences").reset_index()

    # Define VH–VL pair key (V genes)
    merged["vhvl_pair"] = merged["v_gene_H"].astype(str) + "||" + merged["v_gene_L"].astype(str)

    rows = []
    for clone_id, g in merged.groupby("clone_id", sort=False):
        n_rows = int(len(g))
        n_barcodes = int(g["barcode"].nunique())
        if n_barcodes < args.min_clone_barcodes:
            continue

        # For dominance we want counts per VH–VL pair among barcodes.
        # If multiple rows per barcode exist, collapse to one per barcode+pair.
        gb = (
            g[["barcode", "vhvl_pair"]]
            .drop_duplicates()
            .groupby("vhvl_pair")["barcode"]
            .nunique()
            .sort_values(ascending=False)
        )

        counts = gb.values.astype(float)
        n_pairs = int(gb.shape[0])
        dom_pair = str(gb.index[0]) if n_pairs > 0 else ""
        dom_frac = float((counts[0] / counts.sum()) if counts.sum() > 0 else 0.0)
        ent = shannon_entropy_bits(counts)

        rows.append({
            "sample": args.sample,
            "clone_id": clone_id,
            "n_joined_rows": n_rows,
            "n_barcodes_with_pairs": n_barcodes,
            "n_distinct_vhvl_pairs": n_pairs,
            "dominant_vhvl_pair": dom_pair,
            "dominant_pair_frac": dom_frac,
            "pairing_entropy_bits": ent,
        })

    out = pd.DataFrame(rows)
    out = out.merge(clone_size_sequences, on="clone_id", how="left")
    out["clone_size_sequences"] = out["clone_size_sequences"].fillna(0).astype(int)

    # sort: largest paired barcode clones first
    if not out.empty:
        out = out.sort_values(["n_barcodes_with_pairs", "dominant_pair_frac"], ascending=[False, False])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    print(f"Wrote: {out_path}")
    print(f"Rows: {len(out)}")
    if merged is not None:
        print(f"Joined rows total: {len(merged)}")


if __name__ == "__main__":
    main()
