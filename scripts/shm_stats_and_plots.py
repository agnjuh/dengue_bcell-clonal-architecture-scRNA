import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def cliffs_delta(x, y):
    """
    Cliff's delta in [-1, 1]. Positive => x tends to be larger than y.
    x,y: 1D numpy arrays
    """
    # O(n*m)
    gt = 0
    lt = 0
    for xi in x:
        gt += np.sum(xi > y)
        lt += np.sum(xi < y)
    return (gt - lt) / (len(x) * len(y))

def mann_whitney(x, y):
    try:
        from scipy.stats import mannwhitneyu
    except Exception as e:
        raise SystemExit("scipy is required for Mann–Whitney U. Install: pip/conda install scipy") from e

    if len(x) < 5 or len(y) < 5:
        return np.nan
    # two-sided
    stat, p = mannwhitneyu(x, y, alternative="two-sided")
    return float(p)

def load_one(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    # normalize types
    # is_expanded is written as True/False strings by the script
    df["is_expanded"] = df["is_expanded"].astype(str).str.lower().map({"true": True, "false": False})
    df["shm_rate_v"] = pd.to_numeric(df["shm_rate_v"], errors="coerce")
    df = df.dropna(subset=["is_expanded", "shm_rate_v"])
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True, help="e.g. H1 H4 S1 S4")
    ap.add_argument("--metrics_dir", default="results/metrics")
    ap.add_argument("--out_dir", default="results/figures")
    ap.add_argument("--out_tsv", default="results/metrics/shm_stats.summary.tsv")
    args = ap.parse_args()

    metrics_dir = Path(args.metrics_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    long_rows = []

    for s in args.samples:
        f = metrics_dir / f"{s}.shm.tsv"
        if not f.exists():
            raise SystemExit(f"Missing file: {f}")
        df = load_one(f)
        df["sample"] = s
        long_rows.append(df)

        exp = df.loc[df["is_expanded"], "shm_rate_v"].to_numpy()
        non = df.loc[~df["is_expanded"], "shm_rate_v"].to_numpy()

        row = {
            "sample": s,
            "n_expanded": int(len(exp)),
            "n_nonexpanded": int(len(non)),
            "median_expanded": float(np.median(exp)) if len(exp) else np.nan,
            "median_nonexpanded": float(np.median(non)) if len(non) else np.nan,
            "mean_expanded": float(np.mean(exp)) if len(exp) else np.nan,
            "mean_nonexpanded": float(np.mean(non)) if len(non) else np.nan,
            "p_mannwhitney": mann_whitney(exp, non),
            "cliffs_delta_exp_vs_non": cliffs_delta(exp, non) if (len(exp) and len(non)) else np.nan,
        }
        all_rows.append(row)

    summary = pd.DataFrame(all_rows).sort_values("sample")
    summary.to_csv(args.out_tsv, sep="\t", index=False)

    # plotting
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        raise SystemExit(
            "matplotlib is required for plots. "
            
        ) from e

    long_df = pd.concat(long_rows, ignore_index=True)
    # label for plots
    long_df["group"] = np.where(long_df["is_expanded"], "Expanded", "Non-expanded")

    # 1) boxplot per sample (two boxes/sample)
    # Arrange data in consistent order
    samples = list(summary["sample"])
    groups = ["Non-expanded", "Expanded"]

    data = []
    positions = []
    labels = []
    pos = 1
    gap = 1.0
    for s in samples:
        for g in groups:
            v = long_df.loc[(long_df["sample"] == s) & (long_df["group"] == g), "shm_rate_v"].to_numpy()
            data.append(v)
            positions.append(pos)
            labels.append(f"{s}\n{g}")
            pos += 1
        pos += gap  # extra gap between samples

    fig = plt.figure(figsize=(14, 5))
    ax = fig.add_subplot(111)
    ax.boxplot(
        data,
        positions=positions,
        widths=0.6,
        showfliers=False,
        whis=[5, 95],
    )
    ax.set_title("SHM rate (1 - V identity): Expanded vs Non-expanded (per sample)")
    ax.set_ylabel("SHM rate")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=0)
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(out_dir / "SHM_box_by_sample.png", dpi=200)
    plt.close(fig)

    # 2) ECDF plot (one panel per sample; two curves)
    fig = plt.figure(figsize=(14, 8))
    nrows = 2
    ncols = int(np.ceil(len(samples) / nrows))
    for i, s in enumerate(samples, start=1):
        ax = fig.add_subplot(nrows, ncols, i)
        for g in groups:
            v = long_df.loc[(long_df["sample"] == s) & (long_df["group"] == g), "shm_rate_v"].to_numpy()
            v = v[~np.isnan(v)]
            if len(v) == 0:
                continue
            xs = np.sort(v)
            ys = np.arange(1, len(xs) + 1) / len(xs)
            ax.plot(xs, ys, label=g)
        ax.set_title(s)
        ax.set_xlabel("SHM rate")
        ax.set_ylabel("ECDF")
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend()
    fig.suptitle("SHM ECDF: Expanded vs Non-expanded", y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "SHM_ecdf_by_sample.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print("Wrote summary:", args.out_tsv)
    print("Wrote figures:", str(out_dir / "SHM_box_by_sample.png"), "and", str(out_dir / "SHM_ecdf_by_sample.png"))

if __name__ == "__main__":
    main()
