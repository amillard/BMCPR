# BMCPR — general_scripts (auto-generated)

---

## `general_scripts/1_create_fake_headers_for_MS_data.py`

**CLI (`--help`):**
```text
usage: 1_create_fake_headers_for_MS_data.py [-h] [-o OUT] [--prefix PREFIX]
                                            [--start-index START_INDEX]
                                            [--version VERSION] [--translate]
                                            [--min-len-aa MIN_LEN_AA]
                                            gbk

Extract protein sequences from GenBank and write NCBI-like .faa with FAKE
accessions, plus a lookup TSV.

positional arguments:
  gbk                   Input GenBank file (.gbk/.gb/.gbff, optionally .gz)

options:
  -h, --help            show this help message and exit
  -o OUT, --out OUT     Output .faa (default: input basename + .faa)
  --prefix PREFIX       Fake protein accession prefix (default: WP_FAKE)
  --start-index START_INDEX
                        Starting index for fake accessions (default: 1)
  --version VERSION     Version number appended to fake accessions (default:
                        1)
  --translate           If a CDS lacks 'translation', translate from DNA
                        sequence (standard table, honors codon_start).
  --min-len-aa MIN_LEN_AA
                        Drop proteins shorter than this length (aa). Default 0
                        = keep all.
```

**Example:**
```bash
# python general_scripts/1_create_fake_headers_for_MS_data.py --help
```

## `general_scripts/auto_readme.py`


**CLI (`--help`):** *(not detected / no argparse found)*

**Example:**
```bash
# python general_scripts/auto_readme.py --help
```

## `general_scripts/embl_2_faa.py`


**CLI (`--help`):**
```text
usage: embl_2_faa.py [-h] [-o OUTDIR] [-t THREADS] inputs [inputs ...]

Convert EMBL / EMBL.GZ files to FAA using locus_tag (parallel)

positional arguments:
  inputs                Input files or glob patterns (e.g. *.embl qwe/asd/*)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -t THREADS, --threads THREADS
                        Number of worker processes (default: 1)
```

**Example:**
```bash
# python general_scripts/embl_2_faa.py --help
```

## `general_scripts/embl_2_fasta.py`



**CLI (`--help`):**
```text
usage: embl_2_fasta.py [-h] [-o OUTDIR] [-t THREADS] inputs [inputs ...]

Convert EMBL / EMBL.GZ files to FASTA (.fa) using Biopython (parallel)

positional arguments:
  inputs                Input files or glob patterns (e.g. *.embl
                        qwe/asd/*.embl)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -t THREADS, --threads THREADS
                        Number of worker processes (default: 1)
```

**Example:**
```bash
# python general_scripts/embl_2_fasta.py --help
```

## `general_scripts/embl_2_gbk.py`

**CLI (`--help`):**
```text
usage: embl_2_gbk.py [-h] [-o OUTDIR] [-t THREADS] inputs [inputs ...]

Convert EMBL / EMBL.GZ files to GenBank (GBK) using Biopython (parallel)

positional arguments:
  inputs                Input files or glob patterns (e.g. *.embl qwe/asd/*)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -t THREADS, --threads THREADS
                        Number of worker processes (default: 1)
```

**Example:**
```bash
# python general_scripts/embl_2_gbk.py --help
```

## `general_scripts/embl_2_gff3.py`

**CLI (`--help`):**
```text
usage: embl_2_gff3.py [-h] [-o OUTDIR] [-t THREADS] inputs [inputs ...]

Convert EMBL / EMBL.GZ files to GFF3 using Biopython (parallel)

positional arguments:
  inputs                Input files or glob patterns (e.g. *.embl qwe/asd/*)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -t THREADS, --threads THREADS
                        Number of worker processes (default: 1)
```

**Example:**
```bash
# python general_scripts/embl_2_gff3.py --help
```

## `general_scripts/gbk_2_faa.py`

**CLI (`--help`):**
```text
usage: gbk_2_faa.py [-h] [-o OUTDIR] [-t THREADS] inputs [inputs ...]

Convert GBK / GB files to FAA using locus_tag (parallel)

positional arguments:
  inputs                Input files or glob patterns (e.g. *.gbk qwe/asd/*)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -t THREADS, --threads THREADS
                        Number of worker processes (default: 1)
```

**Example:**
```bash
# python general_scripts/gbk_2_faa.py --help
```

## `general_scripts/gbk_2_fasta.py`

**CLI (`--help`):**
```text
usage: gbk_2_fasta.py [-h] [-o OUTDIR] [-t THREADS] inputs [inputs ...]

Convert GBK / GB files to FASTA (.fa) using Biopython (parallel)

positional arguments:
  inputs                Input files or glob patterns (e.g. *.gbk qwe/asd/*)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -t THREADS, --threads THREADS
                        Number of worker processes (default: 1)
```

**Example:**
```bash
# python general_scripts/gbk_2_fasta.py --help
```

## `general_scripts/gbk_2_gff.py`

**CLI (`--help`):**
```text
usage: gbk_2_gff.py [-h] [-o OUTDIR] [-t THREADS] inputs [inputs ...]

Convert GBK / GB files to GFF (.gff) using Biopython/BCBio (parallel)

positional arguments:
  inputs                Input files or glob patterns (e.g. *.gbk qwe/asd/*)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -t THREADS, --threads THREADS
                        Number of worker processes (default: 1)
```

**Example:**
```bash
# python general_scripts/gbk_2_gff.py --help
```
