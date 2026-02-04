#!/usr/bin/env python3

import argparse
import gzip
import os
import sys
from glob import glob
from multiprocessing import Pool
from Bio import SeqIO


def open_embl(path):
    """Open EMBL or EMBL.GZ transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def process_file(args):
    embl_path, outdir = args

    base = os.path.basename(embl_path)

    # Strip extensions
    if base.endswith(".embl.gz"):
        name = base[:-8]
    elif base.endswith(".embl"):
        name = base[:-5]
    else:
        return

    out_path = os.path.join(outdir, f"{name}.gbk")

    records = []

    with open_embl(embl_path) as handle:
        contig_idx = 1
        for record in SeqIO.parse(handle, "embl"):
            record.id = f"{name}_{contig_idx}"
            record.name = record.id
            record.description = ""
            records.append(record)
            contig_idx += 1

    if not records:
        return f"[WARN] {embl_path}: no records found"

    with open(out_path, "w") as out_fh:
        SeqIO.write(records, out_fh, "genbank")

    return f"[OK] {embl_path} → {out_path}"


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
        if f.endswith(".embl") or f.endswith(".embl.gz")
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Convert EMBL / EMBL.GZ files to GenBank (GBK) using Biopython (parallel)"
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input files or glob patterns (e.g. *.embl qwe/asd/*)"
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

    embl_files = collect_inputs(args.inputs)

    if not embl_files:
        print("No EMBL files found.", file=sys.stderr)
        sys.exit(1)

    jobs = [(f, args.outdir) for f in embl_files]

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

