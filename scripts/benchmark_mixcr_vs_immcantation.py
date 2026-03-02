import argparse
import pandas as pd
from pathlib import Path
from math import log

def entropy(counts):
    total = sum(counts)
    if total == 0:
        return 0.0
    e = 0.0
    for c in counts:
        if c <= 0:
            continue
        p = c / total
        e -= p * log(p)
    return e

def nmi(labels_a, labels_b):
    # Normalized Mutual Information (symmetric)
    # labels: list-like same length
    import pandas as pd
    df = pd.DataFrame({"a": labels_a, "b": labels_b})
    ct = pd.crosstab(df["a"], df["b"])
    n = ct.values.sum()
    if n == 0:
        return float("nan")

    # MI
    pa = ct.sum(axis=1) / n
    pb = ct.sum(axis=0) / n
    mi = 0.0
    for i in ct.index:
        for j in ct.columns:
            pij = ct.loc[i, j] / n
            if pij <= 0:
                continue
            mi += pij * log(pij / (pa.loc[i] * pb.loc[j]))
    ha = entropy((pa * n).tolist())
    hb = entropy((pb * n).tolist())
    if ha == 0 or hb == 0:
        return float("nan")
    return mi / ((ha * hb) ** 0.5)

def ari(labels_true, labels_pred):
    # Adjusted Rand Index without sklearn (avoid extra deps)
    import pandas as pd
    from math import comb

    df = pd.DataFrame({"t": labels_true, "p": labels_pred})
    ct = pd.crosstab(df["t"], df["p"]).values
    n = ct.sum()
    if n < 2:
        return float("nan")

    sum_comb_c = sum(comb(int(x), 2) for x in ct.flatten())
    sum_comb_rows = sum(comb(int(x), 2) for x in ct.sum(axis=1))
    sum_comb_cols = sum(comb(int(x), 2) for x in ct.sum(axis=0))
    total_pairs = comb(int(n), 2)

    expected = (sum_comb_rows * sum_comb_cols) / total_pairs if total_pairs else 0.0
    max_index = 0.5 * (sum_comb_rows + sum_comb_cols)
    denom = max_index - expected
    if denom == 0:
        return float("nan")
    return (sum_comb_c - expected) / denom

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--immcantation_clones_tsv", required=True)
    ap.add_argument("--mixcr_map_tsv", required=True)     # readId cloneId
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_summary", required=True)
    ap.add_argument("--min_expanded_size", type=int, default=2)
    args = ap.parse_args()

    imm = pd.read_csv(args.immcantation_clones_tsv, sep="\t", dtype=str).fillna("")
    if "sequence_id" not in imm.columns or "clone_id" not in imm.columns:
        raise SystemExit("Immcantation clones TSV must contain sequence_id and clone_id")

    mix = pd.read_csv(args.mixcr_map_tsv, sep="\t", dtype=str).fillna("")
    # MiXCR exportAlignments header typically: readId cloneId (tab)
    # Normalize column names
    mix.columns = [c.strip() for c in mix.columns]
    if "readId" not in mix.columns or "cloneId" not in mix.columns:
        # fallback: assume first two columns
        if mix.shape[1] >= 2:
            mix = mix.iloc[:, :2]
            mix.columns = ["readId", "cloneId"]
        else:
            raise SystemExit("MiXCR map TSV must contain readId and cloneId columns")

    m = imm[["sequence_id", "clone_id"]].merge(
        mix[["readId", "cloneId"]],
        left_on="sequence_id",
        right_on="readId",
        how="inner"
    ).drop(columns=["readId"])

    # clone sizes
    imm_sizes = imm.groupby("clone_id").size().rename("imm_clone_size")
    mix_sizes = m.groupby("cloneId").size().rename("mixcr_clone_size")

    m = m.join(imm_sizes, on="clone_id")
    m = m.join(mix_sizes, on="cloneId")

    # expanded flags (by size)
    m["imm_expanded"] = m["imm_clone_size"] >= args.min_expanded_size
    m["mixcr_expanded"] = m["mixcr_clone_size"] >= args.min_expanded_size

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(args.out_tsv, sep="\t", index=False)

    # metrics on all matched sequences
    ari_all = ari(m["clone_id"].tolist(), m["cloneId"].tolist())
    nmi_all = nmi(m["clone_id"].tolist(), m["cloneId"].tolist())

    # metrics on union of expanded (either method)
    exp = m[m["imm_expanded"] | m["mixcr_expanded"]].copy()
    ari_exp = ari(exp["clone_id"].tolist(), exp["cloneId"].tolist()) if len(exp) else float("nan")
    nmi_exp = nmi(exp["clone_id"].tolist(), exp["cloneId"].tolist()) if len(exp) else float("nan")

    txt = (
        f"sample\t{args.sample}\n"
        f"matched_sequences\t{len(m)}\n"
        f"ari_all\t{ari_all}\n"
        f"nmi_all\t{nmi_all}\n"
        f"expanded_union_n\t{len(exp)}\n"
        f"ari_expanded_union\t{ari_exp}\n"
        f"nmi_expanded_union\t{nmi_exp}\n"
    )
    Path(args.out_summary).write_text(txt)

if __name__ == "__main__":
    main()
