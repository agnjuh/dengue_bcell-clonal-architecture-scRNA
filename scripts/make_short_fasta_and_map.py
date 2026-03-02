#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--out_fasta", required=True)
    ap.add_argument("-m", "--out_map", required=True)
    args = ap.parse_args()

    inp = Path(args.input)
    out_fa = Path(args.out_fasta)
    out_map = Path(args.out_map)

    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_map.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    with inp.open() as fh, out_fa.open("w") as fo, out_map.open("w") as mo:
        header = None
        seq_chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    n += 1
                    sid = f"{inp.stem.split('.')[0]}_{n:06d}"
                    seq = "".join(seq_chunks).lower()
                    mo.write(f"{sid}\t{header}\n")
                    fo.write(f">{sid}\n{seq}\n")
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            n += 1
            sid = f"{inp.stem.split('.')[0]}_{n:06d}"
            seq = "".join(seq_chunks).lower()
            mo.write(f"{sid}\t{header}\n")
            fo.write(f">{sid}\n{seq}\n")

if __name__ == "__main__":
    main()
