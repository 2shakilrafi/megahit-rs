# Architecture and parity gates

## Current pipeline

1. `needletail` streams FASTA/FASTQ records and transparently decodes gzip or
   bzip2 input.
2. Reads enter bounded batches, preventing the parser from retaining the full
   dataset in memory.
3. Rayon workers count packed directed `(k+1)`-mer edges. Each observed edge is
   accompanied by its reverse complement so graph construction is independent
   of sequenced strand.
4. Batch-local maps are reduced and merged into the run-level edge table.
5. Edges below `--min-count` are removed.
6. The graph is compacted into maximal non-branching paths. Remaining isolated
   cycles are emitted as circular segments.
7. Reverse-complement-equivalent unitigs are deduplicated and sorted
   deterministically by length and sequence.
8. FASTA, GFA, and JSON outputs are written after successful graph construction.

## Why the graph is not yet "succinct"

DNA words are packed into `u128`, avoiding heap-allocated strings for graph
nodes. Adjacency and coverage still live in standard hash maps, however. A true
MEGAHIT-class implementation needs an external-memory builder and a succinct
edge representation supporting rank/select operations without materializing
all adjacency lists.

## Correctness invariants

- Ambiguous bases terminate a k-mer window; no edge crosses an `N` or other
  unsupported symbol.
- A read and its reverse complement contribute the same set of canonical
  biological edges.
- Palindromic `(k+1)`-mers are not double counted.
- Every retained directed edge is visited exactly once during compaction.
- Output ordering is independent of hash-map iteration and worker scheduling.
- Existing output files are not replaced unless `--overwrite` is explicit.

## Stable-release parity gates

The project must not claim MEGAHIT replacement status until it passes all of
the following:

1. **Functional:** paired, interleaved, singleton, gzip, bzip2, checkpoint,
   preset, GFA, and multi-*k* behavior are covered by integration tests.
2. **Assembly:** isolate and CAMI mock-community results are assessed with
   QUAST/metaQUAST, including misassemblies, genome fraction, duplication,
   NGA50, and mismatches—not N50 alone.
3. **Metagenomic:** taxonomic abundance and MAG recovery are tested across low,
   medium, and high-complexity communities.
4. **Performance:** wall time and peak RSS are recorded on fixed, checksummed
   corpora for Apple Silicon, Linux x86-64, and Linux ARM64.
5. **Reproducibility:** deterministic mode produces byte-identical FASTA, GFA,
   and statistics across repeated runs and thread counts.
6. **Robustness:** malformed records, truncated compression streams, disk-full
   conditions, and memory limits fail with actionable diagnostics.

