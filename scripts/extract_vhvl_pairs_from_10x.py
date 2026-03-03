#!/usr/bin/env python3

import argparse
import pandas as pd


def pick_best_per_chain(df_chain: pd.DataFrame) -> pd.DataFrame:
    # Prefer highest UMIs, then reads
    for col in ["umis", "reads"]:
        if col in df_chain.columns:
            df_chain[col] = pd.to_numeric(df_chain[col], errors="coerce").fillna(0).astype(int)
        else:
            df_chain[col] = 0
    df_chain = df_chain.sort_values(["barcode", "umis", "reads"], ascending=[True, False, False])
    return df_chain.groupby("barcode", as_index=False).head(1).copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--vdj_csv", required=True, help="10x filtered_contig_annotations.csv or .csv.gz")
    ap.add_argument("--out_paired", required=True)
    ap.add_argument("--out_summary", required=True)
    ap.add_argument("--out_top_pairs", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.vdj_csv, dtype=str).fillna("")

    # Basic QC filters (10x uses 'true'/'false')
    df = df[(df["is_cell"] == "true") & (df["high_confidence"] == "true") & (df["productive"] == "true")].copy()

    # Keep immunoglobulin chains
    df = df[df["chain"].isin(["IGH", "IGK", "IGL"])].copy()

    # Barcode-level chain counts for pairing classification
    counts = df.groupby(["barcode", "chain"]).size().unstack(fill_value=0)
    has_igh = counts.get("IGH", pd.Series(dtype=int)) > 0
    has_l = (counts.get("IGK", pd.Series(dtype=int)) + counts.get("IGL", pd.Series(dtype=int))) > 0
    candidates = counts.index[has_igh & has_l]

    sub = counts.loc[candidates].copy()
    sub["light_total"] = sub.get("IGK", 0) + sub.get("IGL", 0)

    paired_unique_by_counts = int(((sub.get("IGH", 0) == 1) & (sub["light_total"] == 1)).sum())
    paired_ambiguous_by_counts = int(((sub.get("IGH", 0) > 1) | (sub["light_total"] > 1)).sum())

    df_cand = df[df["barcode"].isin(candidates)].copy()

    heavy = pick_best_per_chain(df_cand[df_cand["chain"] == "IGH"].copy())
    light = pick_best_per_chain(df_cand[df_cand["chain"].isin(["IGK", "IGL"])].copy())

    paired = heavy.merge(light, on="barcode", suffixes=("_H", "_L"), how="inner")

    # attach barcode-level counts for status
    sub2 = sub[["light_total"]].assign(IGH_n=sub.get("IGH", 0)).reset_index().rename(columns={"index": "barcode"})
    paired = paired.merge(sub2, on="barcode", how="left")

    paired["pairing_status"] = "ambiguous"
    paired.loc[(paired["IGH_n"] == 1) & (paired["light_total"] == 1), "pairing_status"] = "unique"

    # Output paired table (light is best-picked IGK/IGL)
    out_cols = [
        "barcode", "pairing_status",
        "contig_id_H", "chain_H", "v_gene_H", "d_gene_H", "j_gene_H", "c_gene_H", "cdr3_H", "reads_H", "umis_H",
        "contig_id_L", "chain_L", "v_gene_L", "j_gene_L", "c_gene_L", "cdr3_L", "reads_L", "umis_L",
    ]
    out_cols = [c for c in out_cols if c in paired.columns]
    paired[out_cols].to_csv(args.out_paired, sep="\t", index=False)

    # Summary TSV as key/value (easy to parse + stable)
    light_counts = paired["chain_L"].value_counts().to_dict() if "chain_L" in paired.columns else {}
    n_igk = int(light_counts.get("IGK", 0))
    n_igl = int(light_counts.get("IGL", 0))

    summary_rows = [
        ("sample", args.sample),
        ("barcodes_with_IGH_and_light", str(int(len(candidates)))),
        ("paired_unique_by_counts", str(paired_unique_by_counts)),
        ("paired_ambiguous_by_counts", str(paired_ambiguous_by_counts)),
        ("paired_rows_written_unique", str(int((paired["pairing_status"] == "unique").sum()))),
        ("paired_rows_written_ambiguous", str(int((paired["pairing_status"] == "ambiguous").sum()))),
        ("light_chain_best_pick_IGK", str(n_igk)),
        ("light_chain_best_pick_IGL", str(n_igl)),
    ]
    pd.DataFrame(summary_rows, columns=["key", "value"]).to_csv(args.out_summary, sep="\t", index=False)

    # Top VH–VL gene pairs (V gene × V gene)
    if "v_gene_H" in paired.columns and "v_gene_L" in paired.columns:
        tp = (
            paired.groupby(["v_gene_H", "v_gene_L"])
            .size()
            .reset_index(name="n")
            .sort_values("n", ascending=False)
            .head(50)
        )
        tp.to_csv(args.out_top_pairs, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["v_gene_H", "v_gene_L", "n"]).to_csv(args.out_top_pairs, sep="\t", index=False)


if __name__ == "__main__":
    main()
