#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_tsv", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--top_n", type=int, default=50)
    ap.add_argument("--max_genes", type=int, default=20)
    ap.add_argument("--log1p", action="store_true")
    ap.add_argument("--mask_zeros", action="store_true")
    ap.add_argument("--normalize", action="store_true",
                    help="normalize counts to relative frequencies per sample")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t")

    if df.empty:
        raise SystemExit(f"Empty input: {args.in_tsv}")

    df = df.sort_values("n", ascending=False).head(args.top_n)

    total_pairs = df["n"].sum()

    if args.normalize:
        df["value"] = df["n"] / total_pairs
    else:
        df["value"] = df["n"]

    top_H = df.groupby("v_gene_H")["value"].sum().sort_values(ascending=False).head(args.max_genes).index
    top_L = df.groupby("v_gene_L")["value"].sum().sort_values(ascending=False).head(args.max_genes).index

    df = df[df["v_gene_H"].isin(top_H) & df["v_gene_L"].isin(top_L)]

    if df.empty:
        raise SystemExit("No data left after filtering.")

    mat = df.pivot_table(index="v_gene_H",
                         columns="v_gene_L",
                         values="value",
                         aggfunc="sum",
                         fill_value=0)

    mat = mat.loc[
        mat.sum(axis=1).sort_values(ascending=False).index,
        mat.sum(axis=0).sort_values(ascending=False).index
    ]

    data = mat.values.astype(float)

    if args.log1p and not args.normalize:
        data = np.log1p(data)

    if args.mask_zeros:
        data = np.ma.masked_where(mat.values == 0, data)

    fig_w = max(6, 0.35 * mat.shape[1] + 2)
    fig_h = max(5, 0.35 * mat.shape[0] + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(data, aspect="auto", interpolation="nearest")

    title = f"{args.sample}: VH–VL pairing"
    if args.normalize:
        title += " (relative frequency)"
    ax.set_title(title)

    ax.set_xlabel("V gene (light)")
    ax.set_ylabel("V gene (heavy)")

    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=90)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index)

    cbar = fig.colorbar(im, ax=ax)
    if args.normalize:
        cbar.set_label("relative frequency")
    elif args.log1p:
        cbar.set_label("log1p(n)")
    else:
        cbar.set_label("n (paired barcodes)")

    fig.tight_layout()
    fig.savefig(args.out_png, dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    main()