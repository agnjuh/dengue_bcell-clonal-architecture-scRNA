import sys
import pandas as pd
from pathlib import Path

inp = Path(sys.argv[1])
out_fa = Path(sys.argv[2])
out_tsv = Path(sys.argv[3])
donor = sys.argv[4]

df = pd.read_csv(inp, compression="gzip")

# keep high-confidence, productive, real cells
df = df[(df["is_cell"] == True) & (df["high_confidence"] == True) & (df["productive"] == True)].copy()

# reconstruct variable-region nucleotide sequence
parts = ["fwr1_nt","cdr1_nt","fwr2_nt","cdr2_nt","fwr3_nt","cdr3_nt","fwr4_nt"]
for c in parts:
    if c not in df.columns:
        raise SystemExit(f"Missing column: {c}")

df["vregion_nt"] = df[parts].fillna("").astype(str).agg("".join, axis=1)
df = df[df["vregion_nt"].str.len() > 100].copy()  # sanity filter

# sequence_id links to cell+contig
df["sequence_id"] = donor + "|" + df["barcode"].astype(str) + "|" + df["contig_id"].astype(str)

# write fasta
with out_fa.open("w") as f:
    for sid, seq in zip(df["sequence_id"], df["vregion_nt"]):
        f.write(f">{sid}\n{seq}\n")

# write minimal table for joining later
keep = ["sequence_id","barcode","contig_id","chain","v_gene","d_gene","j_gene","c_gene","cdr3","cdr3_nt","reads","umis"]
for c in keep:
    if c not in df.columns:
        df[c] = ""
df["donor"] = donor
df[["sequence_id","donor"] + keep[1:]].to_csv(out_tsv, sep="\t", index=False)

print(f"Wrote FASTA: {out_fa}  n={len(df)}")
print(f"Wrote TSV:  {out_tsv}")
