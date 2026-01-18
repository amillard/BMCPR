#!/usr/bin/env python3

"""
For Rachael and Matt for legacy software linked to MS  

Writes a  lookup table at the same time to allow backwards compatability 

Used for the high-throughput proteomics 

Extract protein sequences from a GenBank (.gbk/.gb/.gbff) file and write a
protein FASTA (.faa) with NCBI-like headers, using FAKE protein accessions.

Header format (one line), inspired by NCBI's *.faa export:
>lcl|{contig}_prot_{FAKEACC}_{idx} [locus_tag=XYZ] [protein=desc]
[protein_id={FAKEACC}] [location=complement(123..456)] [gbkey=CDS]

Notes
-----
- FAKE protein accessions are generated as: {prefix}{number:06d}.{version}
  e.g., WP_FAKE000001.1 then WP_FAKE000002.1, etc.
- If a CDS lacks a 'translation' qualifier and --translate is not given,
  it will be skipped (common for GBK files that do not include translations).
- With --translate, sequences are translated from the genomic DNA using the
  specified transl_table (default 11), honoring codon_start if present.
- Multi-exon CDS locations are written as join(a..b,c..d,...) and wrapped in
  complement(join(...)) for minus-strand features.
- Input may be plain or gzipped.

Usage
-----
python gbk_to_faa_like_ncbi.py input.gbk \
  -o out.faa \
  --prefix WP_FAKE \
  --start-index 1 \
  --version 1 \
  [--translate] \
  [--min-len-aa 0]

Outputs
-------
1) out.faa
   Protein FASTA. Sequences have no trailing '*' (stop) character.
2) out.lookup.tsv
   Tab-separated table mapping fake_accession -> annotations:
   fake_accession  original_protein_id  locus_tag  product  location  contig
"""

import argparse
import gzip
from pathlib import Path
from typing import Optional, List

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _int_pos(pos) -> int:
    """
    Convert a Biopython position (ExactPosition, BeforePosition, etc.) to
    1-based int. Biopython positions are 0-based, GBK strings are 1-based.
    """
    try:
        return int(pos) + 1
    except Exception:
        return int(getattr(pos, "position", pos)) + 1


def format_location(feature) -> str:
    """
    Return a GenBank-like location string:
      123..456
      complement(593..3922)
      join(123..456,789..900)
      complement(join(123..456,789..900))
    """
    loc = feature.location
    strand = loc.strand  # +1, -1, or None

    def part_to_str(fl: FeatureLocation) -> str:
        start = _int_pos(fl.start)
        end = int(fl.end)  # Biopython end is exclusive
        return f"{start}..{end}"

    if isinstance(loc, CompoundLocation):
        parts = [part_to_str(p) for p in loc.parts]
        joined = f"join({','.join(parts)})"
        if strand == -1:
            return f"complement({joined})"
        else:
            return joined
    else:
        s = part_to_str(loc)
        if strand == -1:
            return f"complement({s})"
        else:
            return s


def get_first_qualifier(quals: dict, keys: List[str], default: Optional[str] = None) -> Optional[str]:
    for k in keys:
        if k in quals and quals[k]:
            v = quals[k][0]
            if v:
                return v
    return default


def iter_cds_records(gbk_path: Path):
    with open_maybe_gzip(gbk_path, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for i, feat in enumerate(record.features, start=1):
                if feat.type != "CDS":
                    continue
                yield record, i, feat


def main():
    ap = argparse.ArgumentParser(
        description="Extract protein sequences from GenBank and write NCBI-like .faa with FAKE accessions, plus a lookup TSV."
    )
    ap.add_argument("gbk", type=Path, help="Input GenBank file (.gbk/.gb/.gbff, optionally .gz)")
    ap.add_argument("-o", "--out", type=Path, help="Output .faa (default: input basename + .faa)")
    ap.add_argument("--prefix", default="WP_FAKE", help="Fake protein accession prefix (default: %(default)s)")
    ap.add_argument("--start-index", type=int, default=1, help="Starting index for fake accessions (default: %(default)d)")
    ap.add_argument("--version", type=int, default=1, help="Version number appended to fake accessions (default: %(default)d)")
    ap.add_argument("--translate", action="store_true",
                    help="If a CDS lacks 'translation', translate from DNA sequence (standard table, honors codon_start).")
    ap.add_argument("--min-len-aa", type=int, default=0, help="Drop proteins shorter than this length (aa). Default 0 = keep all.")
    args = ap.parse_args()

    gbk_path: Path = args.gbk
    out_path: Path = args.out or gbk_path.with_suffix("").with_suffix(".faa")
    idx = args.start_index
    version = args.version

    n_written = 0
    n_skipped = 0

    # Write FASTA
    with open(out_path, "wt") as fout:
        for record, cds_idx, feat in iter_cds_records(gbk_path):
            quals = feat.qualifiers
            product = get_first_qualifier(quals, ["product", "gene", "note"], default="hypothetical protein")
            locus_tag = get_first_qualifier(quals, ["locus_tag", "gene"], default=f"CDS_{cds_idx}")
            gbkey = "CDS"
            contig = record.id

            # Build fake accession like WP_FAKE000001.1
            fake_acc = f"{args.prefix}{idx:06d}.{version}"

            # Protein sequence (prefer provided translation)
            aa = get_first_qualifier(quals, ["translation"], default=None)
            if aa:
                aa = aa.rstrip('*')
            if aa is None and args.translate:
                try:
                    codon_start = int(get_first_qualifier(quals, ["codon_start"], default="1"))
                except Exception:
                    codon_start = 1
                transl_table = int(get_first_qualifier(quals, ["transl_table"], default="11"))
                seq = feat.extract(record.seq)
                if codon_start > 1:
                    seq = seq[codon_start-1:]
                aa = str(seq.translate(table=transl_table, to_stop=True)).rstrip('*')

            if aa is None:
                n_skipped += 1
                idx += 1
                continue

            if args.min_len_aa > 0 and len(aa) < args.min_len_aa:
                n_skipped += 1
                idx += 1
                continue

            # Location
            loc_str = format_location(feat)

            # Compose header
            header_core = f"lcl|{contig}_prot_{fake_acc}_{cds_idx}"
            header_ann = f"[locus_tag={locus_tag}] [protein={product}] [protein_id={fake_acc}] [location={loc_str}] [gbkey={gbkey}]"
            header = f">{header_core} {header_ann}"

            # Write FASTA (wrap at 60 chars)
            fout.write(header + "\n")
            for i in range(0, len(aa), 60):
                fout.write(aa[i:i+60] + "\n")

            n_written += 1
            idx += 1

    # Write lookup TSV (parallel iteration to keep indices consistent even for skipped CDS)
    tsv_path = out_path.with_suffix(".lookup.tsv")
    with open(tsv_path, "wt") as tsv:
        tsv.write("fake_accession\toriginal_protein_id\tlocus_tag\tproduct\tlocation\tcontig\n")
        idx = args.start_index
        for record, cds_idx, feat in iter_cds_records(gbk_path):
            quals = feat.qualifiers
            product = get_first_qualifier(quals, ["product", "gene", "note"], default="hypothetical protein")
            locus_tag = get_first_qualifier(quals, ["locus_tag", "gene"], default=f"CDS_{cds_idx}")
            contig = record.id
            fake_acc = f"{args.prefix}{idx:06d}.{version}"
            loc_str = format_location(feat)
            aa = get_first_qualifier(quals, ["translation"], default=None)
            if aa is None and args.translate:
                try:
                    codon_start = int(get_first_qualifier(quals, ["codon_start"], default="1"))
                except Exception:
                    codon_start = 1
                transl_table = int(get_first_qualifier(quals, ["transl_table"], default="11"))
                seq = feat.extract(record.seq)
                if codon_start > 1:
                    seq = seq[codon_start-1:]
                aa = str(seq.translate(table=transl_table, to_stop=True)).rstrip('*')
            if aa:
                orig_prot_id = get_first_qualifier(quals, ["protein_id"], default="-")
                tsv.write(f"{fake_acc}\t{orig_prot_id}\t{locus_tag}\t{product}\t{loc_str}\t{contig}\n")
            idx += 1

    print(f"Wrote {n_written} protein sequences to {out_path}")
    if n_skipped:
        print(f"Skipped {n_skipped} CDS without translation (use --translate to derive AA).")
    print(f"Lookup table written to {tsv_path}")


if __name__ == "__main__":
    main()
