import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", default="results/metrics")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    for s in args.samples:
        f = Path(args.metrics_dir) / f"{s}.shm.tsv"
        df = pd.read_csv(f, sep="\t", dtype=str).fillna("")
        df["sample"] = s
        df["shm_rate_v"] = pd.to_numeric(df["shm_rate_v"], errors="coerce")
        df["is_expanded"] = df["is_expanded"].astype(str)
        rows.append(df[["sample","is_expanded","shm_rate_v"]])
    dat = pd.concat(rows, ignore_index=True).dropna(subset=["shm_rate_v"])

    # boxplot data
    samples = args.samples
    data_non = [dat[(dat.sample==s) & (dat.is_expanded=="False")]["shm_rate_v"].values for s in samples]
    data_exp = [dat[(dat.sample==s) & (dat.is_expanded=="True")]["shm_rate_v"].values for s in samples]

    fig = plt.figure(figsize=(10,4.5))
    ax = fig.add_subplot(111)

    # positions
    pos_non = [i*3 + 1 for i in range(len(samples))]
    pos_exp = [i*3 + 2 for i in range(len(samples))]

    b1 = ax.boxplot(data_non, positions=pos_non, widths=0.7, patch_artist=True, showfliers=False)
    b2 = ax.boxplot(data_exp, positions=pos_exp, widths=0.7, patch_artist=True, showfliers=False)

    ax.set_xticks([i*3 + 1.5 for i in range(len(samples))])
    ax.set_xticklabels(samples)
    ax.set_ylabel("SHM rate (1 - V identity)")
    ax.set_title("SHM by sample: non-expanded vs expanded")

    # legend proxy
    ax.plot([], [], label="Non-expanded")
    ax.plot([], [], label="Expanded")
    ax.legend(frameon=False)

    ax.set_xlim(0, len(samples)*3)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.7)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    main()
