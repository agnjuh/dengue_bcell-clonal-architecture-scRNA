#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def pick_seq_col(df):
    for c in ["sequence_id", "SEQUENCE_ID", "sequenceid", "seq_id"]:
        if c in df.columns:
            return c
    raise ValueError(
        "Could not find sequence id column in clones file. "
        f"Columns: {list(df.columns)}"
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clones_tsv", required=True)
    ap.add_argument("--vhvl_pairs_tsv", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_summary", required=True)
    ap.add_argument(
        "--include_ambiguous",
        action="store_true",
        help="If set, keep ambiguous pairings as well (pairing_status preserved).",
    )
    args = ap.parse_args()

    clones = pd.read_csv(args.clones_tsv, sep="\t", dtype=str).fillna("")
    pairs  = pd.read_csv(args.vhvl_pairs_tsv, sep="\t", dtype=str).fillna("")

    seq_col = pick_seq_col(clones)

    if "contig_id_H" not in pairs.columns:
        raise ValueError(f"vhvl_pairs file missing contig_id_H: {args.vhvl_pairs_tsv}")

    if "pairing_status" not in pairs.columns:
        raise ValueError(f"vhvl_pairs file missing pairing_status: {args.vhvl_pairs_tsv}")

    # default: only unique pairing rows
    if not args.include_ambiguous:
        pairs = pairs[pairs["pairing_status"] == "unique"].copy()

    merged = clones.merge(
        pairs,
        left_on=seq_col,
        right_on="contig_id_H",
        how="left",
        suffixes=("", "_vhvl"),
    )

    total = len(clones)
    matched = int((merged["pairing_status"] != "").sum())
    unique = int((merged["pairing_status"] == "unique").sum())
    ambiguous = int((merged["pairing_status"] == "ambiguous").sum())

    light_igk = int((merged.get("chain_L", "") == "IGK").sum())
    light_igl = int((merged.get("chain_L", "") == "IGL").sum())

    out_summary = pd.DataFrame(
        [
            ["clones_rows", str(total)],
            ["clones_rows_with_pairing", str(matched)],
            ["paired_unique_rows", str(unique)],
            ["paired_ambiguous_rows", str(ambiguous)],
            ["light_chain_rows_IGK", str(light_igk)],
            ["light_chain_rows_IGL", str(light_igl)],
        ],
        columns=["key", "value"],
    )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_summary).parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(args.out_tsv, sep="\t", index=False)
    out_summary.to_csv(args.out_summary, sep="\t", index=False)

if __name__ == "__main__":
    main()
