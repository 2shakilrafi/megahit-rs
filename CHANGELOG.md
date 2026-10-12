# Changelog

All notable changes to this project will be documented here. The project uses
[Semantic Versioning](https://semver.org/) after the initial experimental
series.

## 0.1.0 — 2026-10-11

### Added

- Streaming FASTA/FASTQ ingestion with gzip and bzip2 support.
- MEGAHIT-style paired, interleaved, and singleton input flags.
- Packed two-bit fixed-*k* de Bruijn edges for `k=3..63`.
- Parallel bounded-batch edge counting.
- Coverage filtering and deterministic unitig extraction.
- Linear and circular reverse-complement deduplication.
- FASTA, GFA 1.0, and JSON run outputs.
- Unit, integration, compressed-input, determinism, and overlap-assembly tests.

### Known limitations

- No iterative multi-*k* construction.
- No paired-read branch resolution or local assembly.
- No succinct/BOSS graph representation yet.
- No production metagenome validation yet.

