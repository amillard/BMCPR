#!/usr/bin/env python3

import argparse
import gzip
import os
import sys
from glob import glob
from multiprocessing import Pool
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def open_gbk(path):
    """Open GBK / GB / compressed transparently."""
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

    out_path = os.path.join(outdir, f"{name}.faa")

    proteins = []

    with open_gbk(gbk_path) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                quals = feature.qualifiers
                # Only CDS with translation + locus_tag
                if "translation" not in quals or "locus_tag" not in quals:
                    continue
                seq = Seq(quals["translation"][0])
                locus_tag = quals["locus_tag"][0]
                proteins.append(
                    SeqRecord(
                        seq,
                        id=locus_tag,
                        name="",
                        description=""
                    )
                )

    if not proteins:
        return f"[WARN] {gbk_path}: no CDS with locus_tag + translation found"

    with open(out_path, "w") as out_fh:
        SeqIO.write(proteins, out_fh, "fasta")

    return f"[OK] {gbk_path} → {out_path} ({len(proteins)} proteins)"


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
        description="Convert GBK / GB files to FAA using locus_tag (parallel)"
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

