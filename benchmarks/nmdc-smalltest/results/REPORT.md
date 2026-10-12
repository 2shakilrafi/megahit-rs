# NMDC small-test benchmark

This report compares the experimental Rust assembler with MEGAHIT on the same interleaved paired-end read set. It is a development benchmark, not evidence that `megahit-rs` is ready for scientific use.

## Input and protocol

- Input: [NMDC/JGI `smalltest.int.fastq.gz`](https://portal.nersc.gov/project/m3408/test_data/smalltest.int.fastq.gz) (125,000 read pairs; 37,109,226 bases)
- SHA-256: `9cc9f1e2e09aed826efb31123efbff541c655e08aa0dc40621026f9bdcec8c42`
- Replicates: 3 per configuration; reported time and memory are medians
- Threads: 8; minimum k-mer count: 2; minimum contig: 200 bp
- Configurations: fixed k=31 for both assemblers, plus MEGAHIT's stock multi-k defaults
- Input was decompressed once before timed runs; run order was rotated between replicates
- Evidence: [GitHub Actions run 38193820760](https://github.com/2shakilrafi/megahit-rs/actions/runs/38193820760) at commit `59089d2ae7b48ef57ba3625c815135404cc5e363`

## Results

| Assembler | Median wall time | Median peak RSS | Contigs | Total bp | N50 | Largest | ≥1 kbp | Primary reads mapped | Properly paired |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| megahit-rs (k=31) | 28.20 s | 5,580.8 MiB | 19 | 4,247 | 225 | 256 | 0 | 0.04% | 0.01% |
| MEGAHIT (k=31) | 2.79 s | 281.2 MiB | 192 | 48,683 | 246 | 579 | 0 | 0.42% | 0.16% |
| MEGAHIT (default multi-k) | 10.52 s | 317.0 MiB | 216 | 84,754 | 372 | 1,879 | 2 | 0.55% | 0.21% |

## Direct comparisons

- Against single-k MEGAHIT, `megahit-rs` used 10.11× the wall time and 19.85× the peak memory.
- Against default multi-k MEGAHIT, `megahit-rs` used 2.68× the wall time and 17.60× the peak memory.
- Byte-identical FASTA output across replicates: megahit-rs (k=31)=yes, MEGAHIT (k=31)=no, MEGAHIT (default multi-k)=no.

## Environment

- Platform: Linux-6.17.0-1022-azure-x86_64-with-glibc2.39 (x86_64)
- CPU: AMD EPYC 9V45 96-Core Processor (4 logical CPUs available)
- Runner memory: 15,990 MiB
- Tools: megahit-rs 0.1.0; MEGAHIT v1.2.9; minimap2 2.26-r1175; samtools 1.19.2

## What this does and does not show

The dataset is a small official pipeline test input, so it is useful for regression, performance, and smoke testing. Read mapping back to an assembly is only a retention signal; it does not detect chimeras or prove correctness. Because this dataset has no gold-standard source genomes, N50 and mapping percentages must not be interpreted as biological accuracy. A CAMI or defined mock-community benchmark with reference-based misassembly, genome-fraction, and duplication metrics is still required before any scientific parity claim.
