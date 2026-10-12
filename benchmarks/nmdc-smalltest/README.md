# NMDC small-test benchmark

This benchmark compares `megahit-rs` with MEGAHIT 1.2.9 using the 22.6 MB
interleaved Illumina test input published by the NMDC metagenome assembly
workflow. It runs three configurations on the same machine:

1. `megahit-rs` at fixed k=31.
2. MEGAHIT constrained to k=31.
3. MEGAHIT with its default iterative multi-k settings.

Each configuration uses eight threads, minimum k-mer count 2, and minimum
contig length 200 bp. Three replicates are run in rotated order after warming
the compressed input. The runner records wall time, peak resident memory,
contig statistics, SHA-256 reproducibility, and paired read mapping with
Minimap2/Samtools.

Use the **NMDC small-test benchmark** workflow from the repository's Actions
tab. The workflow verifies both downloaded files against
`checksums.sha256` and uploads the JSON results, report, time files, and logs.

This is a fast development benchmark. The NMDC input does not include a
gold-standard reference, so the results cannot establish assembly correctness.
A truth-backed CAMI or defined mock-community benchmark remains a separate
release gate.
