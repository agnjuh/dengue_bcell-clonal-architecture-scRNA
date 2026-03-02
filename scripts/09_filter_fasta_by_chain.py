import sys

meta_tsv = sys.argv[1]
fasta_in = sys.argv[2]
fasta_out = sys.argv[3]
chain_keep = sys.argv[4]  # e.g. IGH

keep = set()
with open(meta_tsv, "r") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx_sid = header.index("sequence_id")
    idx_chain = header.index("chain")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= max(idx_sid, idx_chain):
            continue
        if parts[idx_chain] == chain_keep:
            keep.add(parts[idx_sid])

n_in = 0
n_out = 0
with open(fasta_in, "r") as fin, open(fasta_out, "w") as fout:
    write = False
    for line in fin:
        if line.startswith(">"):
            n_in += 1
            sid = line[1:].strip()
            write = sid in keep
            if write:
                n_out += 1
                fout.write(line)
        else:
            if write:
                fout.write(line)

print(f"FASTA records in: {n_in}")
print(f"FASTA records out ({chain_keep}): {n_out}")
