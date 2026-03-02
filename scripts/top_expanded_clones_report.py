#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def pick_first_nonempty(series: pd.Series) -> str:
    series = series.dropna().astype(str)
    series = series[series.str.len() > 0]
    return series.iloc[0] if len(series) else ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shm_tsv", required=True, help="results/metrics/{sample}.shm.tsv")
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_summary", required=True)
    ap.add_argument("--top_n", type=int, default=30)
    args = ap.parse_args()

    df = pd.read_csv(args.shm_tsv, sep="\t", dtype=str).fillna("")
    needed = {"sequence_id", "clone_id", "is_expanded", "shm_rate_v", "v_call", "j_call"}
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing columns in {args.shm_tsv}: {missing}")

    exp = df[df["is_expanded"].astype(str).str.lower().eq("true")].copy()

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)

    if exp.empty:
        pd.DataFrame(columns=[
            "clone_id","clone_size","n_sequences_with_shm","shm_median","shm_mean","v_call","j_call"
        ]).to_csv(args.out_tsv, sep="\t", index=False)
        Path(args.out_summary).write_text("No expanded clones found.\n")
        return

    exp["shm_rate_v_num"] = pd.to_numeric(exp["shm_rate_v"], errors="coerce")

    rows = []
    for clone_id, g in exp.groupby("clone_id"):
        shm = g["shm_rate_v_num"].dropna()
        rows.append({
            "clone_id": clone_id,
            "clone_size": int(len(g)),
            "n_sequences_with_shm": int(shm.shape[0]),
            "shm_median": float(shm.median()) if shm.shape[0] else float("nan"),
            "shm_mean": float(shm.mean()) if shm.shape[0] else float("nan"),
            "v_call": pick_first_nonempty(g["v_call"]),
            "j_call": pick_first_nonempty(g["j_call"]),
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(["clone_size", "shm_median"], ascending=[False, False])

    top_n = max(1, int(args.top_n))
    out_top = out.head(top_n).copy()
    out_top.to_csv(args.out_tsv, sep="\t", index=False)

    txt = []
    txt.append(f"Input: {args.shm_tsv}")
    txt.append(f"Expanded clones (unique): {out.shape[0]}")
    txt.append(f"Expanded sequences (rows): {exp.shape[0]}")
    txt.append(f"TopN: {top_n}")
    txt.append("")
    txt.append("Top 10 (clone_id, size, shm_median):")
    for _, r in out.head(10).iterrows():
        txt.append(f"  {r['clone_id']}\tsize={int(r['clone_size'])}\tshm_median={r['shm_median']:.6f}")
    txt.append("")
    Path(args.out_summary).write_text("\n".join(txt) + "\n")

if __name__ == "__main__":
    main()
