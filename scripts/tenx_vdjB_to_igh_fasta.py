#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd


def open_maybe_gzip(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def build_sequence_from_regions(row: pd.Series) -> str:
    parts = []
    for col in ["fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt"]:
        val = row.get(col, "")
        if pd.isna(val):
            val = ""
        parts.append(str(val))
    seq = "".join(parts).strip().lower()
    return seq


def sanitize_nt(seq: str) -> str:
    # allow IUPAC ambiguity codes commonly seen in germline/VDJ contexts
    allowed = set("acgtnryswkmbdhv")
    seq = seq.strip().lower().replace(" ", "").replace("\n", "")
    bad = {c for c in seq if c not in allowed}
    if bad:
        # do not hard-fail: warn and keep sequence (downstream tools might still accept it)
        sys.stderr.write(f"WARNING: sequence contains non-IUPAC chars: {sorted(bad)}\n")
    return seq


def main():
    ap = argparse.ArgumentParser(description="Export 10x VDJ-B contig annotations -> FASTA (robust to missing `sequence` column).")
    ap.add_argument("-i", "--input", required=True, help="10x filtered_contig_annotations.csv(.gz)")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA")
    ap.add_argument("--chain", default="IGH", help="Chain to filter (default: IGH)")
    ap.add_argument("--only-productive", action="store_true", help="Keep only productive contigs")
    ap.add_argument("--min-umis", type=int, default=0, help="Optional: filter by UMIs >= N")
    args = ap.parse_args()

    inp = args.input
    outp = args.output

    df = pd.read_csv(inp, compression="infer", low_memory=False)

    # Basic required columns
    required_base = ["barcode", "contig_id", "chain"]
    missing_base = [c for c in required_base if c not in df.columns]
    if missing_base:
        raise SystemExit(f"Missing required columns in {inp}: {missing_base}\nAvailable columns: {list(df.columns)}")

    # Filter chain
    df = df[df["chain"].astype(str) == str(args.chain)].copy()

    # Productive filter (if requested and column exists)
    if args.only_productive:
        if "productive" not in df.columns:
            raise SystemExit(f"--only-productive requested but `productive` column is missing in {inp}")
        df = df[df["productive"].astype(str).str.lower().isin(["true", "t", "1", "yes"])].copy()

    # UMI filter (if requested and column exists)
    if args.min_umis > 0:
        if "umis" not in df.columns:
            sys.stderr.write("WARNING: --min-umis requested but `umis` column missing; skipping UMI filter.\n")
        else:
            df = df[pd.to_numeric(df["umis"], errors="coerce").fillna(0).astype(int) >= args.min_umis].copy()

    if df.empty:
        raise SystemExit(f"No contigs left after filtering (chain={args.chain}, only_productive={args.only_productive}).")

    # Determine how to obtain nucleotide sequence
    has_sequence_col = "sequence" in df.columns
    has_region_cols = all(c in df.columns for c in ["fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt"])

    if not has_sequence_col and not has_region_cols:
        raise SystemExit(
            f"Cannot export FASTA: no `sequence` column and missing region columns.\n"
            f"Available columns: {list(df.columns)}"
        )

    # Prepare output
    Path(outp).parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with open(outp, "w") as fout:
        for _, row in df.iterrows():
            if has_sequence_col:
                seq = row.get("sequence", "")
                if pd.isna(seq):
                    seq = ""
                seq = str(seq)
            else:
                seq = build_sequence_from_regions(row)

            seq = sanitize_nt(seq)
            if not seq:
                continue

            barcode = str(row.get("barcode", "NA"))
            contig_id = str(row.get("contig_id", "NA"))

            # Optional metadata (safe even if missing)
            v_gene = str(row.get("v_gene", "NA"))
            d_gene = str(row.get("d_gene", "NA"))
            j_gene = str(row.get("j_gene", "NA"))
            c_gene = str(row.get("c_gene", "NA"))
            productive = str(row.get("productive", "NA"))
            umis = str(row.get("umis", "NA"))
            reads = str(row.get("reads", "NA"))

            header = f">{barcode}|{contig_id}|chain={args.chain}|v={v_gene}|d={d_gene}|j={j_gene}|c={c_gene}|productive={productive}|umis={umis}|reads={reads}"
            fout.write(header + "\n")

            # wrap 80
            for i in range(0, len(seq), 80):
                fout.write(seq[i:i+80] + "\n")

            n_written += 1

    if n_written == 0:
        raise SystemExit("No sequences written (all sequences empty after parsing).")

    sys.stderr.write(f"Wrote {n_written} FASTA records to {outp}\n")


if __name__ == "__main__":
    main()
