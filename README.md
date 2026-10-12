# megahit-rs

[![CI](https://github.com/2shakilrafi/megahit-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/2shakilrafi/megahit-rs/actions/workflows/ci.yml)

`megahit-rs` is an experimental, deterministic metagenome assembler written in
Rust. Version 0.1 implements a bounded-memory input pipeline and a packed,
fixed-*k* compacted de Bruijn graph. It is a working foundation for a modern
MEGAHIT-compatible assembler; it is **not yet a drop-in scientific replacement
for MEGAHIT**.

## What works in 0.1

- FASTA and FASTQ input, including gzip and bzip2 compression.
- MEGAHIT-style paired (`-1`/`-2`), interleaved (`--12`), and singleton/orphan
  (`-r`) input selection.
- Validation of paired record counts and interleaved record parity.
- Packed two-bit k-mers and `(k+1)`-mer edges for `k=3..63`.
- Reverse-complement-normalized graph evidence.
- Parallel, bounded-batch edge counting.
- Solid-edge filtering by coverage.
- Deterministic maximal non-branching path extraction.
- Reverse-complement deduplication.
- FASTA contigs, GFA 1.0 segments/links, and a JSON run report.
- Linux, macOS, x86-64, and ARM64-compatible Rust code.

## Important limitations

- Only one fixed *k* is used. MEGAHIT's iterative multi-*k* algorithm is not
  implemented yet.
- Paired files are validated, but pairing constraints are not yet used to
  resolve branches.
- The graph uses packed integers plus hash maps; it is compact but not yet a
  BOSS/succinct de Bruijn graph.
- Graph cleaning currently consists of abundance filtering and unitig
  compaction. Tip clipping, bubble classification, mercy edges, local assembly,
  and repeat resolution remain roadmap work.
- Assembly quality has not yet been validated on CAMI or production datasets.

Do not substitute version 0.1 for MEGAHIT in consequential scientific analyses
without independent validation.

## Build

```bash
cargo build --release
```

## Run

Paired-end reads:

```bash
./target/release/megahit-rs \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o assembly \
  -k 31 \
  --min-count 2 \
  --min-contig-len 200
```

Singleton or orphan reads can be supplied alongside paired reads:

```bash
./target/release/megahit-rs \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -r orphan.fastq.gz \
  -o assembly
```

Outputs:

- `final.contigs.fa`
- `final.graph.gfa`
- `assembly.stats.json`

## Development checks

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
```

## Roadmap

1. Streaming external-memory k-mer partitioning.
2. Succinct/BOSS graph storage with rank/select.
3. Tip clipping and topology-aware bubble handling.
4. Iterative multi-*k* graph construction.
5. Paired-read path support and local assembly.
6. Checkpoint/resume and MEGAHIT preset compatibility.
7. CAMI, mock-community, and isolate quality benchmarks against MEGAHIT,
   metaSPAdes, and Minia.

See [`ARCHITECTURE.md`](ARCHITECTURE.md) for the current design and the parity
gates required before a stable release.

## License

GPL-3.0-or-later.

