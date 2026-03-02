import argparse
import pandas as pd
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--airr", required=True)
    ap.add_argument("--clones", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_summary", required=True)
    args = ap.parse_args()

    airr = pd.read_csv(args.airr, sep="\t", dtype=str).fillna("")
    clones = pd.read_csv(args.clones, sep="\t", dtype=str).fillna("")

    if "sequence_id" not in airr.columns:
        raise SystemExit("AIRR missing sequence_id")
    if "v_identity" not in airr.columns:
        raise SystemExit("AIRR missing v_identity")
    if "sequence_id" not in clones.columns or "clone_id" not in clones.columns:
        raise SystemExit("clones file missing sequence_id and/or clone_id")

    # v_identity (%) -> SHM rate
    v_id = pd.to_numeric(airr["v_identity"], errors="coerce")
    airr["shm_rate_v"] = 1.0 - (v_id / 100.0)

    # expanded clones
    sizes = clones.groupby("clone_id").size()
    expanded_ids = set(sizes[sizes > 1].index.tolist())
    clones["is_expanded"] = clones["clone_id"].isin(expanded_ids)

    # merge SHM onto clones table
    m = clones[["sequence_id", "clone_id", "is_expanded"]].merge(
        airr[["sequence_id", "v_call", "j_call", "v_identity", "shm_rate_v"]],
        on="sequence_id",
        how="left"
    )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(args.out_tsv, sep="\t", index=False)

    # summary
    def summarize(x):
        x = pd.to_numeric(x, errors="coerce").dropna()
        return (int(x.shape[0]),
                float(x.mean()) if x.shape[0] else float("nan"),
                float(x.median()) if x.shape[0] else float("nan"))

    n_e, mean_e, med_e = summarize(m.loc[m["is_expanded"], "shm_rate_v"])
    n_n, mean_n, med_n = summarize(m.loc[~m["is_expanded"], "shm_rate_v"])

    txt = (
        f"AIRR: {args.airr}\n"
        f"CLONES: {args.clones}\n"
        f"OUT: {args.out_tsv}\n\n"
        f"Expanded:   n={n_e} mean={mean_e:.6f} median={med_e:.6f}\n"
        f"Non-exp:    n={n_n} mean={mean_n:.6f} median={med_n:.6f}\n"
    )
    Path(args.out_summary).write_text(txt)

if __name__ == "__main__":
    main()
