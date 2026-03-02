#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd


def run(cmd, log_path: Path, stdout_path: Path | None = None) -> None:
    """
    Run command, write combined stdout+stderr to log_path.
    If stdout_path is given, command stdout is redirected to that file and NOT duplicated in the log.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)

    with open(log_path, "w", encoding="utf-8") as log_fh:
        if stdout_path is None:
            p = subprocess.run(cmd, text=True, capture_output=True)
            log_fh.write(p.stdout or "")
            log_fh.write(p.stderr or "")
        else:
            with open(stdout_path, "w", encoding="utf-8") as out_fh:
                p = subprocess.run(cmd, text=True, stdout=out_fh, stderr=subprocess.PIPE)
            log_fh.write(p.stderr or "")

    if p.returncode != 0:
        raise RuntimeError(f"Command failed (see {log_path}):\n" + " ".join(map(str, cmd)))


def parse_args():
    ap = argparse.ArgumentParser(description="Run MiXCR v4 from FASTA and export clones + readId->cloneId map.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--species", default="hsa")
    ap.add_argument("--preset", default="generic-amplicon")
    ap.add_argument("--rna", action="store_true", help="Pass --rna to mixcr align")
    ap.add_argument("--left_boundary", choices=["floating", "rigid"], default="floating")
    ap.add_argument("--right_boundary_mode", choices=["floating", "rigid"], default="floating")
    ap.add_argument("--right_boundary_anchor", choices=["J", "C"], default="J")
    ap.add_argument("--assemble_by", default="CDR3", help="Gene feature(s) for clonotype assembly, e.g. CDR3 or FR3+CDR3")
    ap.add_argument("--force", action="store_true", help="Force overwrite MiXCR intermediate files")
    ap.add_argument("--out_map_tsv", required=True)
    ap.add_argument("--out_clones_tsv", required=True)
    return ap.parse_args()


def main():
    args = parse_args()

    fasta = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_map_tsv = Path(args.out_map_tsv)
    out_clones_tsv = Path(args.out_clones_tsv)
    out_map_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_clones_tsv.parent.mkdir(parents=True, exist_ok=True)

    # intermediates
    align_vdjca = outdir / "align.vdjca"
    assemble_clna = outdir / "assemble.clna"
    tmp_alignments_tsv = outdir / "exportAlignments.readIds_cloneId.tsv"

    # logs
    log_align = outdir / "mixcr.align.log"
    log_assemble = outdir / "mixcr.assemble.log"
    log_export_clones = outdir / "mixcr.exportClones.log"
    log_export_alignments = outdir / "mixcr.exportAlignments.log"

    # 1) align
    cmd_align = [
        "mixcr", "align",
        "--preset", args.preset,
        "-s", args.species,
        "-t", str(args.threads),
        "--assemble-clonotypes-by", args.assemble_by,
    ]
    if args.rna:
        cmd_align.append("--rna")

    if args.left_boundary == "floating":
        cmd_align.append("--floating-left-alignment-boundary")

    if args.right_boundary_mode == "floating":
        cmd_align.extend(["--floating-right-alignment-boundary", args.right_boundary_anchor])

    if args.force:
        cmd_align.append("--force-overwrite")

    cmd_align.extend([str(fasta), str(align_vdjca)])
    run(cmd_align, log_align)

    # 2) assemble
    cmd_assemble = [
        "mixcr", "assemble",
        "--write-alignments",
        "--assemble-clonotypes-by", args.assemble_by,
    ]
    if args.force:
        cmd_assemble.append("--force-overwrite")
    cmd_assemble.extend([str(align_vdjca), str(assemble_clna)])
    run(cmd_assemble, log_assemble)

    # 3) export clones TSV
    # MiXCR v4: use -nFeature / -aaFeature instead of -nSeqCDR3 / -aaSeqCDR3
    cmd_export_clones = [
        "mixcr", "exportClones",
        "-cloneId",
        "-readCount",
        "-readFraction",
        "-topChains",
        "-nFeature", "CDR3",
        "-aaFeature", "CDR3",
    ]
    if args.force:
        cmd_export_clones.append("--force-overwrite")
    cmd_export_clones.append(str(assemble_clna))
    run(cmd_export_clones, log_export_clones, stdout_path=out_clones_tsv)

    if not out_clones_tsv.exists() or out_clones_tsv.stat().st_size == 0:
        raise RuntimeError("Empty clones TSV produced.\n"
                           f"See: {log_export_clones}")

    # 4) export alignments from .clna to get readIds + cloneId
    # This is the reliable way to build readId->cloneId mapping in MiXCR v4.
    cmd_export_alignments = [
        "mixcr", "exportAlignments",
        "-cloneId",
        "-readIds",
        str(assemble_clna),
    ]
    if args.force:
        cmd_export_alignments.append("--force-overwrite")
    run(cmd_export_alignments, log_export_alignments, stdout_path=tmp_alignments_tsv)

    if not tmp_alignments_tsv.exists() or tmp_alignments_tsv.stat().st_size == 0:
        raise RuntimeError("Empty exportAlignments TSV produced.\n"
                           f"See: {log_export_alignments}")

    df = pd.read_csv(tmp_alignments_tsv, sep="\t", dtype=str).fillna("")
    cols_lower = {c.lower(): c for c in df.columns}

    if "cloneid" not in cols_lower:
        raise RuntimeError(f"Could not find cloneId column in exportAlignments output. Columns: {list(df.columns)}\nSee: {tmp_alignments_tsv}")

    # MiXCR may name it readIds or readId depending on command/version; handle both.
    read_col = None
    if "readids" in cols_lower:
        read_col = cols_lower["readids"]
    elif "readid" in cols_lower:
        read_col = cols_lower["readid"]
    else:
        raise RuntimeError(f"Could not find readId/readIds column in exportAlignments output. Columns: {list(df.columns)}\nSee: {tmp_alignments_tsv}")

    clone_col = cols_lower["cloneid"]

    out = df[[read_col, clone_col]].copy()
    out[read_col] = out[read_col].astype(str).str.strip()
    out = out[out[read_col] != ""]

    # If readIds is comma-separated list, explode; if single id, split yields 1 element.
    out[read_col] = out[read_col].str.split(r"\s*,\s*")
    out = out.explode(read_col)

    out = out.rename(columns={read_col: "readId", clone_col: "cloneId"})[["readId", "cloneId"]]
    out.to_csv(out_map_tsv, sep="\t", index=False)

    if not out_map_tsv.exists() or out_map_tsv.stat().st_size == 0:
        raise RuntimeError("Empty readId_to_cloneId TSV produced.")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(str(e), file=sys.stderr)
        raise