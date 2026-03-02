#!/usr/bin/env python3
import argparse
import re
import pandas as pd

def clean_junction(s: str) -> str:
    s = s.upper()
    s = re.sub(r"[^ACGTN]", "N", s)
    s = re.sub(r"N{2,}", "N", s)
    return s

def parse_igblast_outfmt7(path):
    junction = {}
    current = None
    want_next = False

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("# Query:"):
                current = line.split(":", 1)[1].strip()
                want_next = False

            elif line.startswith("# V-(D)-J junction details"):
                want_next = True

            elif want_next:
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) >= 5 and current:
                    s = "".join(parts[:5])
                    s = re.sub(r"\s+", "", s)
                    s = s.replace("(", "").replace(")", "").replace("-", "")
                    junction[current] = s
                want_next = False

    return junction

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--airr", required=True)
    ap.add_argument("--igblast7", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    m = pd.read_csv(args.idmap, sep="\t", header=None, names=["short_id", "orig_id"], dtype=str)
    short2orig = dict(zip(m["short_id"], m["orig_id"]))

    short_junc = parse_igblast_outfmt7(args.igblast7)
    orig_junc = { short2orig[sid]: j for sid, j in short_junc.items() if sid in short2orig }

    df = pd.read_csv(args.airr, sep="\t", dtype=str).fillna("")
    prod = df.get("productive", "").replace({"True":"T","False":"F","TRUE":"T","FALSE":"F"})
    prod = prod.where(prod.isin(["T","F"]), "F")

    junc = df["sequence_id"].map(orig_junc).fillna("").map(clean_junction)

    out = pd.DataFrame({
        "sequence_id": df.get("sequence_id",""),
        "sequence_input": df.get("sequence",""),
        "v_call": df.get("v_call",""),
        "d_call": df.get("d_call",""),
        "j_call": df.get("j_call",""),
        "junction": junc,
        "junction_length": junc.str.len().astype(str).replace("0",""),
        "functional": prod
    })

    out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
