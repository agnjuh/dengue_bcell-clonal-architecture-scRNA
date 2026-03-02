import argparse
from pathlib import Path
import pandas as pd

def spearman_corr(x, y):
    # scipy nélkül: pandas Spearman
    s = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(s) < 3:
        return float("nan")
    return float(s["x"].corr(s["y"], method="spearman"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--clones_dir", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_per_sample_dir", required=True)
    args = ap.parse_args()

    metrics_dir = Path(args.metrics_dir)
    clones_dir  = Path(args.clones_dir)
    out_tsv     = Path(args.out_tsv)
    per_dir     = Path(args.out_per_sample_dir)
    per_dir.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    for sample in args.samples:
        shm_path   = metrics_dir / f"{sample}.shm.tsv"
        clones_path= clones_dir  / f"{sample}.clones.tsv"

        if not shm_path.exists():
            raise SystemExit(f"[{sample}] Missing SHM file: {shm_path}")
        if not clones_path.exists():
            raise SystemExit(f"[{sample}] Missing clones file: {clones_path}")

        shm = pd.read_csv(shm_path, sep="\t", dtype=str).fillna("")
        clones = pd.read_csv(clones_path, sep="\t", dtype=str).fillna("")

        # normalize column names
        shm.columns = [c.strip().lower() for c in shm.columns]
        clones.columns = [c.strip().lower() for c in clones.columns]

        # required columns
        need_shm = {"clone_id", "shm_rate_v"}
        need_clo = {"clone_id"}
        if not need_shm.issubset(set(shm.columns)):
            raise SystemExit(f"[{sample}] SHM columns missing. Have: {shm.columns.tolist()}")
        if not need_clo.issubset(set(clones.columns)):
            raise SystemExit(f"[{sample}] clones columns missing. Have: {clones.columns.tolist()}")

        # clone sizes
        clone_sizes = clones["clone_id"].value_counts().rename("clone_size").reset_index()
        clone_sizes = clone_sizes.rename(columns={"index": "clone_id"})

        # numeric shm
        shm["shm_rate_v"] = pd.to_numeric(shm["shm_rate_v"], errors="coerce")

        # merge per sequence -> then aggregate per clone
        m = shm.merge(clone_sizes, on="clone_id", how="left")

        # per-clone summary (median SHM)
        per_clone = (
            m.groupby("clone_id", as_index=False)
             .agg(clone_size=("clone_size", "first"),
                  shm_median=("shm_rate_v", "median"),
                  shm_mean=("shm_rate_v", "mean"),
                  n_seq=("shm_rate_v", "count"))
        )
        per_clone["is_expanded"] = per_clone["clone_size"] > 1

        # correlations
        rho_all = spearman_corr(per_clone["clone_size"], per_clone["shm_median"])
        rho_exp = spearman_corr(per_clone.loc[per_clone["is_expanded"], "clone_size"],
                                per_clone.loc[per_clone["is_expanded"], "shm_median"])

        # write per-sample table
        per_out = per_dir / f"{sample}.clone_size_vs_shm.per_clone.tsv"
        per_clone.sort_values(["clone_size", "shm_median"], ascending=[False, False]).to_csv(
            per_out, sep="\t", index=False
        )

        rows.append({
            "sample": sample,
            "n_clones": int(per_clone.shape[0]),
            "n_expanded_clones": int(per_clone["is_expanded"].sum()),
            "rho_spearman_all": rho_all,
            "rho_spearman_expanded_only": rho_exp,
            "per_clone_tsv": str(per_out)
        })

    out = pd.DataFrame(rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    print("Wrote:", out_tsv)
    print(out.to_string(index=False))

if __name__ == "__main__":
    main()
