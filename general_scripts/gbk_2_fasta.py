#!/usr/bin/env python3

import argparse
import gzip
import os
import sys
from glob import glob
from multiprocessing import Pool
from Bio import SeqIO


def open_gbk(path):
    """Open GBK / GB / compressed variants transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def process_file(args):
    gbk_path, outdir = args

    base = os.path.basename(gbk_path)

    # Strip extensions
    for ext in (".gbk.gz", ".gb.gz", ".gbk", ".gb"):
        if base.endswith(ext):
            name = base[:-len(ext)]
            break
    else:
        return

    out_path = os.path.join(outdir, f"{name}.fa")

    records_out = []

    with open_gbk(gbk_path) as handle:
        contig_idx = 1
        for record in SeqIO.parse(handle, "genbank"):
            record.id = f"{name}_{contig_idx}"
            record.name = ""
            record.description = ""
            records_out.append(record)
            contig_idx += 1

    if not records_out:
        return f"[WARN] {gbk_path}: no records found"

    with open(out_path, "w") as out_fh:
        SeqIO.write(records_out, out_fh, "fasta")

    return f"[OK] {gbk_path} → {out_path}"


def collect_inputs(patterns):
    files = []
    for p in patterns:
        matches = glob(p)
        if matches:
            files.extend(matches)
        else:
            files.append(p)

    return [
        f for f in files
        if f.endswith((".gbk", ".gb", ".gbk.gz", ".gb.gz"))
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Convert GBK / GB files to FASTA (.fa) using Biopython (parallel)"
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input files or glob patterns (e.g. *.gbk qwe/asd/*)"
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of worker processes (default: 1)"
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    gbk_files = collect_inputs(args.inputs)

    if not gbk_files:
        print("No GenBank files found.", file=sys.stderr)
        sys.exit(1)

    jobs = [(f, args.outdir) for f in gbk_files]

    if args.threads == 1:
        for job in jobs:
            msg = process_file(job)
            if msg:
                print(msg)
    else:
        with Pool(processes=args.threads) as pool:
            for msg in pool.imap_unordered(process_file, jobs):
                if msg:
                    print(msg)


if __name__ == "__main__":
    main()

