#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec


def read_paired_unique(summary_tsv: Path) -> int:
    s = pd.read_csv(summary_tsv, sep="\t")
    if "key" not in s.columns or "value" not in s.columns:
        raise SystemExit(f"Unexpected summary format: {summary_tsv}")
    d = dict(zip(s["key"].astype(str), s["value"].astype(str)))
    if "paired_rows_written_unique" in d:
        return int(float(d["paired_rows_written_unique"]))
    if "paired_unique_by_counts" in d:
        return int(float(d["paired_unique_by_counts"]))
    raise SystemExit(f"paired unique not found in: {summary_tsv}")


def load_pairs_top(tsv_path: Path, paired_unique: int, top_n: int) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", dtype={"v_gene_H": str, "v_gene_L": str, "n": int})
    if df.empty:
        return pd.DataFrame(columns=["v_gene_H", "v_gene_L", "pct"])
    df = df.sort_values("n", ascending=False).head(top_n).copy()
    df["pct"] = 100.0 * df["n"] / float(paired_unique) if paired_unique > 0 else 0.0
    return df[["v_gene_H", "v_gene_L", "pct"]]


def build_matrix(df: pd.DataFrame, H_genes, L_genes) -> pd.DataFrame:
    mat = pd.DataFrame(0.0, index=H_genes, columns=L_genes)
    for _, r in df.iterrows():
        h, l, v = r["v_gene_H"], r["v_gene_L"], float(r["pct"])
        if h in mat.index and l in mat.columns:
            mat.loc[h, l] += v
    return mat


def compute_clip_vmax(values: np.ndarray, q: float, exclude_zeros: bool) -> float:
    x = values[np.isfinite(values)]
    if exclude_zeros:
        x = x[x > 0]
    if x.size == 0:
        return 0.0
    return float(np.quantile(x, q))


def _safe_supxlabel(fig, text: str, fontsize: float, y: float = 0.03):
    if hasattr(fig, "supxlabel"):
        fig.supxlabel(text, fontsize=fontsize, y=y)
    else:
        fig.text(0.5, y, text, ha="center", va="center", fontsize=fontsize)


def _safe_supylabel(fig, text: str, fontsize: float, x: float = 0.05):
    if hasattr(fig, "supylabel"):
        fig.supylabel(text, fontsize=fontsize, x=x)
    else:
        fig.text(x, 0.5, text, ha="center", va="center", rotation=90, fontsize=fontsize)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--top_n_pairs", type=int, default=50)
    ap.add_argument("--max_genes", type=int, default=15, help="max distinct V genes per axis (global)")

    # scaling
    ap.add_argument(
        "--scale_mode",
        choices=["linear", "log", "clip"],
        default="clip",
        help="linear: raw %, log: LogNorm over nonzero, clip: linear with vmax=quantile",
    )
    ap.add_argument("--clip_q", type=float, default=0.99, help="quantile for clip mode (e.g. 0.95-0.995)")
    ap.add_argument("--clip_exclude_zero", action="store_true", help="exclude zeros when computing clip vmax")
    ap.add_argument("--dpi", type=int, default=300)

    # label density + font
    ap.add_argument("--x_tick_step", type=int, default=1, help="show every Nth x tick label")
    ap.add_argument("--y_tick_step", type=int, default=1, help="show every Nth y tick label")
    ap.add_argument("--tick_font", type=float, default=10.0)
    ap.add_argument("--panel_title_font", type=float, default=12.0)
    ap.add_argument("--suptitle_font", type=float, default=13.0)

    # layout knobs
    ap.add_argument("--wspace", type=float, default=0.05)
    ap.add_argument("--hspace", type=float, default=0.14)
    ap.add_argument("--cbar_width", type=float, default=0.045, help="relative width of colorbar column")
    ap.add_argument("--top_margin", type=float, default=0.975, help="suptitle y position (0-1)")

    args = ap.parse_args()
    metrics_dir = Path(args.metrics_dir)

    # load per-sample dfs
    dfs = {}
    for s in args.samples:
        top_tsv = metrics_dir / f"{s}.vhvl_vgene_pairs.top{args.top_n_pairs}.tsv"
        sum_tsv = metrics_dir / f"{s}.vhvl_pairing.summary.tsv"
        paired_unique = read_paired_unique(sum_tsv)
        dfs[s] = load_pairs_top(top_tsv, paired_unique, args.top_n_pairs)

    all_df = pd.concat([dfs[s].assign(sample=s) for s in args.samples], ignore_index=True)
    if all_df.empty:
        raise SystemExit("No data loaded.")

    # choose global top genes (across samples)
    top_H = (
        all_df.groupby("v_gene_H")["pct"].sum()
        .sort_values(ascending=False)
        .head(args.max_genes)
        .index.tolist()
    )
    top_L = (
        all_df.groupby("v_gene_L")["pct"].sum()
        .sort_values(ascending=False)
        .head(args.max_genes)
        .index.tolist()
    )

    # build matrices
    mats = {}
    global_max = 0.0
    global_nonzero_min = np.inf
    for s in args.samples:
        df = dfs[s]
        df = df[df["v_gene_H"].isin(top_H) & df["v_gene_L"].isin(top_L)]
        mat = build_matrix(df, top_H, top_L)
        mats[s] = mat

        arr = mat.values.astype(float)
        if arr.size:
            global_max = max(global_max, float(np.nanmax(arr)))
            nz = arr[arr > 0]
            if nz.size:
                global_nonzero_min = min(global_nonzero_min, float(np.min(nz)))

    # decide grid
    n = len(args.samples)
    if n == 1:
        nrows, ncols = 1, 1
    elif n == 2:
        nrows, ncols = 1, 2
    else:
        nrows, ncols = 2, 2

    # figure sizing
    panel_w = max(4.8, 0.24 * len(top_L) + 2.2)
    panel_h = max(4.2, 0.22 * len(top_H) + 2.0)
    fig_w = panel_w * ncols + 1.1
    fig_h = panel_h * nrows + 0.6

    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = GridSpec(
        nrows=nrows,
        ncols=ncols + 1,
        figure=fig,
        width_ratios=[1] * ncols + [args.cbar_width],
        wspace=args.wspace,
        hspace=args.hspace,
    )

    fig.subplots_adjust(left=0.18, bottom=0.14, top=0.93)

    #  normalization
    norm = None
    vmin = 0.0
    vmax = global_max

    all_vals = np.concatenate([mats[s].values.ravel() for s in args.samples]).astype(float)

    if args.scale_mode == "clip":
        vmax = compute_clip_vmax(all_vals, args.clip_q, exclude_zeros=args.clip_exclude_zero)
        vmin = 0.0
        norm = None
    elif args.scale_mode == "log":
        if not np.isfinite(global_nonzero_min) or global_nonzero_min <= 0:
            raise SystemExit("Log scale requested but no non-zero values found.")
        vmin = global_nonzero_min
        vmax = max(global_max, vmin * 10.0)
        norm = LogNorm(vmin=vmin, vmax=vmax)

    # consistent axis limits
    xlim = (-0.5, len(top_L) - 0.5)
    ylim = (len(top_H) - 0.5, -0.5)

    # plot panels
    im = None
    for i in range(nrows * ncols):
        r = i // ncols
        c = i % ncols

        ax = fig.add_subplot(gs[r, c])

        if i >= n:
            ax.axis("off")
            continue

        sname = args.samples[i]
        arr = mats[sname].values.astype(float)

        if args.scale_mode == "log":
            im = ax.imshow(arr, aspect="auto", interpolation="nearest", norm=norm)
        else:
            im = ax.imshow(arr, aspect="auto", interpolation="nearest", vmin=vmin, vmax=vmax)

        ax.set_title(sname, fontsize=args.panel_title_font, pad=6)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        ax.set_xticks(range(len(top_L)))
        ax.set_yticks(range(len(top_H)))

        # X labels only on bottom row
        if r == nrows - 1:
            xlabels = [g if (j % args.x_tick_step == 0) else "" for j, g in enumerate(top_L)]
            ax.set_xticklabels(xlabels, rotation=90, fontsize=args.tick_font)
            ax.tick_params(labelbottom=True)
        else:
            ax.set_xticklabels([])
            ax.tick_params(labelbottom=False)

        # Y labels only on leftmost col
        if c == 0:
            ylabels = [g if (j % args.y_tick_step == 0) else "" for j, g in enumerate(top_H)]
            ax.set_yticklabels(ylabels, fontsize=args.tick_font)
            ax.tick_params(labelleft=True)
        else:
            ax.set_yticklabels([])
            ax.tick_params(labelleft=False)

        ax.tick_params(axis="both", which="both", length=2)

    # external colorbar
    cax = fig.add_subplot(gs[:, -1])
    cbar = fig.colorbar(im, cax=cax)

    if args.scale_mode == "log":
        cbar.set_label("Pair frequency (% unique paired barcodes, log)")
    elif args.scale_mode == "clip":
        cbar.set_label("Pair frequency (% unique paired barcodes)")
    else:
        cbar.set_label("Pair frequency (% unique paired barcodes)")

    fig.suptitle(
        "VH–VL V-gene pairing architecture",
        fontsize=args.suptitle_font,
        y=args.top_margin,
    )

    _safe_supxlabel(fig, "V gene (light)", fontsize=args.tick_font + 2.0, y=0.03)
    _safe_supylabel(fig, "V gene (heavy)", fontsize=args.tick_font + 2.0, x=0.05)

    out = Path(args.out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)


if __name__ == "__main__":
    main()