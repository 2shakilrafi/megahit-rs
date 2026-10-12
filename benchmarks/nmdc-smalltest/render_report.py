#!/usr/bin/env python3
"""Render benchmark.json as a compact, reviewable Markdown report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


DISPLAY_NAMES = {
    "megahit-rs-k31": "megahit-rs (k=31)",
    "megahit-k31": "MEGAHIT (k=31)",
    "megahit-default": "MEGAHIT (default multi-k)",
}


def mib(value: int) -> float:
    return value / (1024**2)


def pct_ratio(numerator: float, denominator: float) -> str:
    if denominator == 0:
        return "n/a"
    return f"{numerator / denominator:.2f}×"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark_json", type=Path)
    parser.add_argument("output_markdown", type=Path)
    args = parser.parse_args()

    data = json.loads(args.benchmark_json.read_text())
    benchmark = data["benchmark"]
    summary = data["summary"]
    assembly = data["assembly_metrics"]
    mapping = data["read_mapping"]
    system = data["system"]
    versions = data["versions"]
    run_url = None
    if system.get("github_run_id"):
        run_url = (
            "https://github.com/2shakilrafi/megahit-rs/actions/runs/"
            + system["github_run_id"]
        )

    lines = [
        "# NMDC small-test benchmark",
        "",
        "This report compares the experimental Rust assembler with MEGAHIT on the "
        "same interleaved paired-end read set. It is a development benchmark, not "
        "evidence that `megahit-rs` is ready for scientific use.",
        "",
        "## Input and protocol",
        "",
        f"- Input: [NMDC/JGI `smalltest.int.fastq.gz`](https://portal.nersc.gov/project/m3408/test_data/smalltest.int.fastq.gz) "
        f"({benchmark['read_pairs']:,} read pairs; {benchmark['read_bases']:,} bases)",
        f"- SHA-256: `{benchmark['input_sha256']}`",
        f"- Replicates: {benchmark['replicates']} per configuration; reported time and memory are medians",
        f"- Threads: {benchmark['threads']}; minimum k-mer count: {benchmark['minimum_kmer_count']}; "
        f"minimum contig: {benchmark['minimum_contig_length']} bp",
        "- Configurations: fixed k=31 for both assemblers, plus MEGAHIT's stock multi-k defaults",
        "- Input was decompressed once before timed runs; run order was rotated between replicates",
        f"- Evidence: [GitHub Actions run {system['github_run_id']}]({run_url}) at commit `{system['github_sha']}`"
        if run_url
        else "- Evidence: local run",
        "",
        "## Results",
        "",
        "| Assembler | Median wall time | Median peak RSS | Contigs | Total bp | N50 | Largest | ≥1 kbp | Primary reads mapped | Properly paired |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name in ("megahit-rs-k31", "megahit-k31", "megahit-default"):
        timing = summary[name]
        metrics = assembly[name]
        mapped = mapping[name]
        lines.append(
            f"| {DISPLAY_NAMES[name]} | {timing['median_wall_seconds']:.2f} s | "
            f"{mib(timing['median_max_rss_bytes']):,.1f} MiB | {metrics['contigs']:,} | "
            f"{metrics['total_bp']:,} | {metrics.get('n50_bp', 0):,} | "
            f"{metrics.get('largest_bp', 0):,} | {metrics.get('contigs_ge_1000', 0):,} | "
            f"{mapped['primary_mapped_percent']:.2f}% | {mapped['properly_paired_percent']:.2f}% |"
        )

    rust = summary["megahit-rs-k31"]
    fixed = summary["megahit-k31"]
    default = summary["megahit-default"]
    lines.extend(
        [
            "",
            "## Direct comparisons",
            "",
            f"- Against single-k MEGAHIT, `megahit-rs` used "
            f"{pct_ratio(rust['median_wall_seconds'], fixed['median_wall_seconds'])} the wall time and "
            f"{pct_ratio(rust['median_max_rss_bytes'], fixed['median_max_rss_bytes'])} the peak memory.",
            f"- Against default multi-k MEGAHIT, `megahit-rs` used "
            f"{pct_ratio(rust['median_wall_seconds'], default['median_wall_seconds'])} the wall time and "
            f"{pct_ratio(rust['median_max_rss_bytes'], default['median_max_rss_bytes'])} the peak memory.",
            "- Byte-identical FASTA output across replicates: "
            + ", ".join(
                f"{DISPLAY_NAMES[name]}={'yes' if value else 'no'}"
                for name, value in data["byte_deterministic_across_replicates"].items()
            )
            + ".",
            "",
            "## Environment",
            "",
            f"- Platform: {system['platform']} ({system['machine']})",
            f"- CPU: {system['cpu_model']} ({system['logical_cpus']} logical CPUs available)",
            f"- Runner memory: {mib(system['memory_kib'] * 1024):,.0f} MiB" if system.get("memory_kib") else "- Runner memory: unavailable",
            f"- Tools: {versions['megahit_rs']}; {versions['megahit']}; "
            f"minimap2 {versions['minimap2']}; {versions['samtools']}",
            "",
            "## What this does and does not show",
            "",
            "The dataset is a small official pipeline test input, so it is useful for regression, "
            "performance, and smoke testing. Read mapping back to an assembly is only a retention "
            "signal; it does not detect chimeras or prove correctness. Because this dataset has no "
            "gold-standard source genomes, N50 and mapping percentages must not be interpreted as "
            "biological accuracy. A CAMI or defined mock-community benchmark with reference-based "
            "misassembly, genome-fraction, and duplication metrics is still required before any "
            "scientific parity claim.",
            "",
        ]
    )
    args.output_markdown.write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
