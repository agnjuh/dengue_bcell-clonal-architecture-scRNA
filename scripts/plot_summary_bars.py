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

    rows=[]
    for s in args.samples:
        f = Path(args.metrics_dir) / f"{s}.shm.tsv"
        df = pd.read_csv(f, sep="\t", dtype=str).fillna("")
        df["shm_rate_v"] = pd.to_numeric(df["shm_rate_v"], errors="coerce")
        exp = (df["is_expanded"].astype(str)=="True")
        expanded_frac = float(exp.mean())
        median_shm = float(df["shm_rate_v"].median())
        rows.append({"sample": s, "expanded_frac": expanded_frac, "median_shm": median_shm})

    summ = pd.DataFrame(rows)

    fig = plt.figure(figsize=(10,4.5))
    ax1 = fig.add_subplot(111)

    x = range(len(summ))
    ax1.bar(x, summ["expanded_frac"])
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(summ["sample"])
    ax1.set_ylabel("Expanded fraction")
    ax1.set_title("Expanded fraction per sample (bars); SHM median (points)")
    ax1.set_ylim(0, max(0.35, float(summ["expanded_frac"].max())*1.2))

    ax2 = ax1.twinx()
    ax2.plot(list(x), summ["median_shm"], marker="o", linestyle="-")
    ax2.set_ylabel("Median SHM rate")

    ax1.grid(True, axis="y", linestyle=":", linewidth=0.7)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    main()
