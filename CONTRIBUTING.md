# Contributing

Contributions are welcome, particularly when they add an independently
verifiable correctness or performance improvement.

## Ground rules

- Do not copy source code from MEGAHIT or another assembler into this project.
  Cite papers and describe independently implemented algorithms in the pull
  request.
- Treat assembly quality as the primary outcome. A speedup that increases
  misassemblies is not an improvement unless exposed as an explicit tradeoff.
- Include a regression test for every bug fix.
- Keep output deterministic unless nondeterminism is explicitly documented and
  justified.
- Never use N50 alone as evidence of assembly quality.

## Required checks

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo build --release
```

Performance changes should include the corpus checksum, exact command, CPU,
memory, operating system, wall time, and peak resident memory.

