import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--per_clone_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--max_points", type=int, default=40000, help="downsample if huge")
    args = ap.parse_args()

    per_dir = Path(args.per_clone_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # combined figure (one panel per sample)
    for sample in args.samples:
        inp = per_dir / f"{sample}.clone_size_vs_shm.per_clone.tsv"
        if not inp.exists():
            raise SystemExit(f"Missing: {inp}")

        df = pd.read_csv(inp, sep="\t")
        # safety
        for c in ["clone_size","shm_median","is_expanded"]:
            if c not in df.columns:
                raise SystemExit(f"{inp} missing column: {c}")

        # optional downsample (mostly for huge S1)
        if len(df) > args.max_points:
            df = df.sample(args.max_points, random_state=1)

        # log10(clone_size)
        x = np.log10(df["clone_size"].astype(float))
        y = df["shm_median"].astype(float)
        exp = df["is_expanded"].astype(bool)

        plt.figure(figsize=(6,5))
        # singletons first, then expanded on top
        plt.scatter(x[~exp], y[~exp], s=8, alpha=0.25, label="singleton (size=1)")
        plt.scatter(x[exp],  y[exp],  s=14, alpha=0.6,  label="expanded (size>1)")
        plt.xlabel("log10(clone size)")
        plt.ylabel("median SHM rate (1 - v_identity/100)")
        plt.title(f"{sample}: clone size vs SHM (per-clone)")
        plt.legend(frameon=False)
        plt.tight_layout()

        out_png = out_dir / f"{sample}.clone_size_vs_shm.scatter.png"
        plt.savefig(out_png, dpi=200)
        plt.close()

    # overview: all samples in one plot
    plt.figure(figsize=(7,5))
    for sample in args.samples:
        inp = per_dir / f"{sample}.clone_size_vs_shm.per_clone.tsv"
        df = pd.read_csv(inp, sep="\t")
        if len(df) > args.max_points:
            df = df.sample(args.max_points, random_state=1)

        x = np.log10(df["clone_size"].astype(float))
        y = df["shm_median"].astype(float)
        plt.scatter(x, y, s=8, alpha=0.35, label=sample)

    plt.xlabel("log10(clone size)")
    plt.ylabel("median SHM rate (1 - v_identity/100)")
    plt.title("Clone size vs SHM (per-clone), all samples")
    plt.legend(frameon=False)
    plt.tight_layout()
    out_png = out_dir / "ALL.clone_size_vs_shm.scatter.png"
    plt.savefig(out_png, dpi=200)
    plt.close()

    print("Wrote figures to:", out_dir)

if __name__ == "__main__":
    main()
