# TRUST4 Output File Reference

This document lists every file produced by `run-trust4` / `trust4` / `annotator`, the format of each row, and the meaning of every column. Files are ordered along the pipeline (extraction → assembly → annotation → reporting).

The `<prefix>` placeholder is the output prefix from `-o` (default: inferred from input file name).
The `[ ]` symbols mark fields that are optional / conditionally present.

Stages (controlled by `--stage`):

| Stage | Step | Tool |
|---|---|---|
| 0 | Candidate read extraction | `bam-extractor` / `fastq-extractor` |
| 1 | De novo assembly | `trust4` (main.cpp) |
| 2 | Annotation against IMGT | `annotator` |
| 3 | Reporting | `trust-simplerep.pl`, `trust-barcoderep.pl`, `trust-airr.pl`, `trust-augment-cdr3.pl` |

---

## 1. Stage-0 outputs — extracted candidate reads

These are intermediate FASTQ/FASTA files written by the extractor. They are deleted by `--clean 1`.

### `<prefix>_toassemble_1.fq` / `<prefix>_toassemble_2.fq` (paired-end) or `<prefix>_toassemble.fq` (single-end)

Standard 4-line FASTQ. One record per candidate read that hit a V/D/J/C k-mer.

| Line | Field | Description |
|---|---|---|
| 1 | `@<read_id>` | Original read ID from input BAM/FASTQ |
| 2 | sequence | Read sequence (may be reverse-complemented) |
| 3 | `+` | Separator |
| 4 | quality | Phred quality string |

### `<prefix>_toassemble_bc.fa` (only if `--barcode` was used)

Two-line FASTA, one record per read in `_toassemble`, in the same order.

| Line | Field | Description |
|---|---|---|
| 1 | `><read_id>` | Read ID (matches FASTQ) |
| 2 | barcode sequence | Cell barcode after extraction (and translation if `--barcodeTranslate` was set) |

### `<prefix>_toassemble_umi.fa` (only if `--UMI` was used)

Same layout as `_toassemble_bc.fa`, with the UMI sequence on line 2.

---

## 2. Stage-1 outputs — assembly

### `<prefix>_raw.out` and `<prefix>_final.out`

Per-contig assembly with per-position read coverage. `_raw.out` is the pre-mate-extension version; `_final.out` is the version that goes into annotation. Format is identical.

Each contig spans **6 lines** (or 2 lines if no posWeight is recorded):

| Line | Format | Description |
|---|---|---|
| 1 | `>assemble<i> <name>` *(bulk)* or `>BARCODE_<i> <name>` *(barcode mode)* | Consensus header. `<i>` is the internal contig index. `<name>` is the assembly name set during construction (often the seed gene). |
| 2 | nucleotide string | Consensus sequence, length L |
| 3 | `c1 c2 ... cL` | A counts at each position (space-separated, length L) |
| 4 | `c1 c2 ... cL` | C counts at each position |
| 5 | `c1 c2 ... cL` | G counts at each position |
| 6 | `c1 c2 ... cL` | T counts at each position |

Per-position read coverage = sum of the four count rows at that position. `min_internal_cov` / `mean_cov` / `max_cov` (in `_cdr3.out`) are computed from these four rows by `trust-augment-cdr3.pl`.

### `<prefix>_assembled_reads.fa`

Per-read FASTA showing which reads were incorporated into which assembly. Two lines per read:

| Line | Format | Description |
|---|---|---|
| 1 | `><read_id> <strand> <minCnt> <medianCnt>[ barcode:<bc>][ umi:<umi>]` | `<strand>` is 0/1; `<minCnt>` and `<medianCnt>` are k-mer support counts used during overlap; `barcode:` and `umi:` only present if those features were on |
| 2 | sequence | Read sequence (orientation already adjusted) |

---

## 3. Stage-2 outputs — annotation

### `<prefix>_annot.fa`

FASTA-like file. Header line carries all annotation; sequence line is the consensus.

```
>consensus_id length avg_cov V_annot D_annot J_annot C_annot CDR1 CDR2 CDR3
<consensus_seq>
```

Header fields (whitespace-separated):

| # | Field | Description |
|---|---|---|
| 1 | `consensus_id` | `assemble<i>` in bulk mode, `<barcode>_<i>` in barcode mode |
| 2 | `length` | Consensus length in bp |
| 3 | `avg_cov` | Average per-base read coverage (float) |
| 4 | `V_annot` | V gene call(s); see *gene annot* format below. Up to 3 candidates separated by `,`. `*` = no hit |
| 5 | `D_annot` | D gene call(s) |
| 6 | `J_annot` | J gene call(s) |
| 7 | `C_annot` | C gene call(s) |
| 8 | `CDR1` | CDR1 annotation; see *CDR annot* below. `*` = none |
| 9 | `CDR2` | CDR2 annotation |
| 10 | `CDR3` | CDR3 annotation |

Gene annot format (one per gene candidate):

```
gene_name(ref_gene_length):(consensus_start-consensus_end):(ref_start-ref_end):similarity
```

CDR annot format:

```
CDRx(consensus_start-consensus_end):score=sequence
```

For CDR1/2 the `score` is similarity to germline. For CDR3 it has special meaning: `0.00` = partial CDR3, `1.00` = CDR3 with imputed nucleotides, other = motif signal strength (100.00 strongest). All coordinates are 0-based half-open.

### `<prefix>_cdr3.out`

TSV without header. One row per CDR3 detected (a contig may emit multiple rows if it carries more than one CDR3 region). Columns 14–18 are appended by `trust-augment-cdr3.pl`.

| # | Column | Type | Description |
|---|---|---|---|
| 1 | `consensus_id` | string | `assemble<i>` (bulk) or `<barcode>_<i>` (barcode) |
| 2 | `index_within_consensus` | int | Index of this CDR3 within the contig (0-based) |
| 3 | `V_gene` | string | V gene call, `*` if missing |
| 4 | `D_gene` | string | D gene call |
| 5 | `J_gene` | string | J gene call |
| 6 | `C_gene` | string | C gene call |
| 7 | `CDR1` | string | CDR1 nucleotide sequence, `*` if missing |
| 8 | `CDR2` | string | CDR2 nucleotide sequence |
| 9 | `CDR3` | string | CDR3 nucleotide sequence |
| 10 | `CDR3_score` | float | `0.00` partial, `1.00` imputed, otherwise motif strength divided by 100 (max 1.00) |
| 11 | `read_fragment_count` | float | EM-distributed CDR3-level abundance. May be < 1 because one read can contribute fractional weight to multiple CDR3s |
| 12 | `CDR3_germline_similarity` | float | % similarity of CDR3 nt to germline |
| 13 | `complete_vdj_assembly` | int | 1 if the contig spans a complete V(D)J, 0 otherwise |
| 14 | `length` | int | Contig length in bp (from `_final.out`) |
| 15 | `min_internal_cov` | int | Minimum per-base coverage across the whole contig — the most useful confidence metric. Real contigs stay above some N at every position; spurious merges typically drop to 1× somewhere |
| 16 | `mean_cov` | float | Mean per-base coverage |
| 17 | `max_cov` | int | Maximum per-base coverage |
| 18 | `read_count` | int | **Only present when TRUST4 was run with `--outputReadAssignment`.** Raw integer count of reads assigned to this contig in `_assign.out`. Differs from column 11: this is per-contig, no EM, no fractional weights |

> ⚠ Cell-calling note: for "how many reads support this barcode/contig" use **column 18 (`read_count`)** when available; it is the closest analog to CellRanger's filtered UMI count. Column 11 (`read_fragment_count`) is fine for clonal abundance ranking but is fractional and CDR3-level, not contig-level.

### `<prefix>_assign.out` (only with `--outputReadAssignment`)

Two-column TSV without header. One row per assembled read.

| # | Column | Description |
|---|---|---|
| 1 | `read_id` | Read ID (matches `_toassemble`) |
| 2 | `contig_id` | The assembly this read was assigned to (`assemble<i>` or `<barcode>_<i>`) |

Reads that failed to be assigned are not written. `read_count` in `_cdr3.out` column 18 is `count(read_id) GROUP BY contig_id` over this file.

### `<prefix>_airr.tsv`

[AIRR Rearrangement format](https://docs.airr-community.org/en/latest/datarep/rearrangements.html). Tab-separated **with header**. Per-CDR3 rows (one row per `(contig, CDR3 index)`).

| # | Column | Description |
|---|---|---|
| 1 | `sequence_id` | `<contig>_<idx>` — combines columns 1+2 of `_cdr3.out` |
| 2 | `sequence` | Full contig consensus |
| 3 | `rev_comp` | `T`/`F`, whether sequence was reverse-complemented |
| 4 | `productive` | `T`/`F`, productive rearrangement |
| 5 | `locus` | IGH / IGK / IGL / TRA / TRB / TRD / TRG |
| 6 | `v_call` | V gene call |
| 7 | `d_call` | D gene call |
| 8 | `j_call` | J gene call |
| 9 | `c_call` | C gene call |
| 10 | `sequence_alignment` | Aligned sequence portion (V→J) |
| 11 | `germline_alignment` | Inferred germline alignment |
| 12 | `cdr1` | CDR1 nucleotide |
| 13 | `cdr2` | CDR2 nucleotide |
| 14 | `junction` | Junction (CDR3 + flanking codons) nucleotide |
| 15 | `junction_aa` | Junction amino acid |
| 16 | `v_cigar` | V alignment CIGAR |
| 17 | `d_cigar` | D alignment CIGAR |
| 18 | `j_cigar` | J alignment CIGAR |
| 19 | `c_cigar` | C alignment CIGAR |
| 20 | `v_identity` | V % identity |
| 21 | `j_identity` | J % identity |
| 22 | `cell_id` | Cell barcode (only in barcode mode; empty otherwise) |
| 23 | `complete_vdj` | `T`/`F`, complete VDJ assembly |
| 24 | `consensus_count` | Number of reads supporting the consensus (≈ `read_fragment_count` rounded) |

---

## 4. Stage-3 outputs — reports

### `<prefix>_report.tsv`

Simple bulk-style repertoire report. Tab-separated **with header** (header line starts with `#`).

| # | Column | Description |
|---|---|---|
| 1 | `count` | Read count (or barcode count in barcode mode) supporting this CDR3 |
| 2 | `frequency` | Proportion of `count` within the same chain class (IG and TR are normalized separately) |
| 3 | `CDR3nt` | CDR3 nucleotide sequence |
| 4 | `CDR3aa` | CDR3 amino acid sequence; `_` = stop codon, `?` = ambiguous (N in codon), `out_of_frame` if length not divisible by 3 |
| 5 | `V` | V gene |
| 6 | `D` | D gene |
| 7 | `J` | J gene |
| 8 | `C` | C gene |
| 9 | `cid` | Source consensus_id |
| 10 | `cid_full_length` | 1 if the source contig is full-length VDJ, 0 otherwise |

### `<prefix>_barcode_report.tsv` (barcode mode only)

Per-cell summary. Tab-separated **with header** (line starts with `#`). One row per barcode.

| # | Column | Description |
|---|---|---|
| 1 | `barcode` | Cell barcode (translated if `--barcodeTranslate` was used) |
| 2 | `cell_type` | Inferred cell type (`B`, `abT`, `gdT`, `multi`, etc.) based on chain combination |
| 3 | `chain1` | Most-abundant heavy / β / δ chain (IGH / TRB / TRD), as CSV (see below). `*` if absent |
| 4 | `chain2` | Most-abundant light / α / γ chain (IGK / IGL / TRA / TRG), as CSV |
| 5 | `secondary_chain1` | Next most-abundant heavy/β/δ chain. `*` if absent |
| 6 | `secondary_chain2` | Next most-abundant light/α/γ chain |

Each chain field is a 10-column CSV (no spaces):

```
V_gene,D_gene,J_gene,C_gene,cdr3_nt,cdr3_aa,read_cnt,consensus_id,CDR3_germline_similarity,consensus_complete_vdj
```

| # | CSV field | Description |
|---|---|---|
| 1 | `V_gene` | V gene call |
| 2 | `D_gene` | D gene call |
| 3 | `J_gene` | J gene call |
| 4 | `C_gene` | C gene call |
| 5 | `cdr3_nt` | CDR3 nucleotide |
| 6 | `cdr3_aa` | CDR3 amino acid |
| 7 | `read_cnt` | Read/UMI/barcode count for this chain in this cell |
| 8 | `consensus_id` | Source contig |
| 9 | `CDR3_germline_similarity` | % similarity |
| 10 | `consensus_complete_vdj` | 1 if full-length |

If a barcode has more than two chains in the same class, additional ones are appended into `chain2`/`secondary_chain*` separated by `;`.

### `<prefix>_barcode_airr.tsv` (barcode mode only)

Same column layout as `<prefix>_airr.tsv`, but rows are deduplicated to the chains listed in `_barcode_report.tsv` (i.e. the representative chain per barcode), and `cell_id` is always populated.

---

## 5. Files removed by `--clean`

| `--clean` value | Behavior |
|---|---|
| `0` (default) | Keep everything |
| `1` | Delete `_toassemble_*`, `_assembled_reads.fa`, `_final.out`, `_raw.out`, `_airr_align.tsv` |
| `2` | Same as `1`, plus delete `_barcode_report.tsv`, keep only AIRR files |

> Note: `_cdr3.out` augmentation by `trust-augment-cdr3.pl` reads `_final.out` and `_assign.out`, so if you plan to re-run augmentation later, run it **before** `--clean 1`.

---

## 6. Quick file-purpose cheat sheet

| File | Stage | Per-row unit | Use it for |
|---|---|---|---|
| `_toassemble_*.fq` | 0 | candidate read | re-running stage 1+ from cached extraction |
| `_raw.out` | 1 | contig + 4 cov rows | inspecting pre-mate-extension assemblies |
| `_assembled_reads.fa` | 1 | read | which reads went into which contig |
| `_final.out` | 1 | contig + 4 cov rows | per-base coverage source for QC |
| `_annot.fa` | 2 | contig | one-line annotation summary per contig |
| `_cdr3.out` | 2 | CDR3 | primary CDR3-level table; cell-calling input |
| `_assign.out` | 2 | read | exact read→contig mapping (only with `--outputReadAssignment`) |
| `_airr.tsv` | 2 | CDR3 | AIRR-compliant export |
| `_report.tsv` | 3 | CDR3 (collapsed) | bulk repertoire report |
| `_barcode_report.tsv` | 3 | barcode | per-cell chain summary |
| `_barcode_airr.tsv` | 3 | barcode chain | AIRR-compliant per-cell export |
