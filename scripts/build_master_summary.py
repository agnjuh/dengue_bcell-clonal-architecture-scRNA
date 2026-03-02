import argparse
import pandas as pd
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--metrics_dir", required=True)
    ap.add_argument("--clones_dir", required=True)
    ap.add_argument("--spearman_tsv", required=True)
    ap.add_argument("--out_tsv", required=True)
    args = ap.parse_args()

    rows = []

    # Spearman summary (clone size vs SHM)
    spearman = pd.read_csv(args.spearman_tsv, sep="\t", dtype=str).fillna("")

    for s in args.samples:
        metrics_tsv = Path(args.metrics_dir) / f"{s}.shm.tsv"
        clones_tsv  = Path(args.clones_dir) / f"{s}.clones.tsv"

        shm = pd.read_csv(metrics_tsv, sep="\t", dtype=str).fillna("")
        clones = pd.read_csv(clones_tsv, sep="\t", dtype=str).fillna("")

        # total sequences
        n_total = len(shm)

        # clone stats
        sizes = clones.groupby("clone_id").size()
        n_clones = len(sizes)
        n_expanded_clones = int((sizes > 1).sum())
        expanded_fraction = float(n_expanded_clones / n_clones) if n_clones else 0.0
        max_clone_size = int(sizes.max()) if n_clones else 0

        # SHM stats
        shm["shm_rate_v"] = pd.to_numeric(shm["shm_rate_v"], errors="coerce")
        exp_mask = shm["is_expanded"] == "True"

        med_exp = float(shm.loc[exp_mask, "shm_rate_v"].median())
        med_non = float(shm.loc[~exp_mask, "shm_rate_v"].median())
        mean_exp = float(shm.loc[exp_mask, "shm_rate_v"].mean())
        mean_non = float(shm.loc[~exp_mask, "shm_rate_v"].mean())

        # Spearman
        sp = spearman[spearman["sample"] == s]
        rho_all = float(sp["rho_spearman_all"].iloc[0])
        rho_exp = float(sp["rho_spearman_expanded_only"].iloc[0])

        rows.append({
            "sample": s,
            "n_total_sequences": n_total,
            "n_clones": n_clones,
            "n_expanded_clones": n_expanded_clones,
            "expanded_clone_fraction": expanded_fraction,
            "max_clone_size": max_clone_size,
            "median_shm_expanded": med_exp,
            "median_shm_nonexpanded": med_non,
            "mean_shm_expanded": mean_exp,
            "mean_shm_nonexpanded": mean_non,
            "rho_spearman_all": rho_all,
            "rho_spearman_expanded_only": rho_exp
        })

    out = pd.DataFrame(rows).sort_values("sample")
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print("Wrote:", args.out_tsv)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
