#!/usr/bin/env python3
"""Run a reproducible small metagenome benchmark for megahit-rs and MEGAHIT."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path


def split_interleaved(source: Path, read1: Path, read2: Path) -> None:
    """Split an interleaved FASTQ while validating record and mate structure."""
    if read1.exists() and read2.exists():
        return
    read1.parent.mkdir(parents=True, exist_ok=True)
    records = 0
    with gzip.open(source, "rt") as incoming, gzip.open(read1, "wt") as first, gzip.open(
        read2, "wt"
    ) as second:
        while True:
            pair = []
            for _ in range(2):
                lines = [incoming.readline() for _ in range(4)]
                if not lines[0]:
                    if pair:
                        raise ValueError("interleaved FASTQ ended after only one mate")
                    return
                if any(line == "" for line in lines[1:]):
                    raise ValueError("truncated FASTQ record")
                if not lines[0].startswith("@") or not lines[2].startswith("+"):
                    raise ValueError("invalid FASTQ record")
                if len(lines[1].rstrip()) != len(lines[3].rstrip()):
                    raise ValueError("sequence and quality lengths differ")
                pair.append(lines)
            first.writelines(pair[0])
            second.writelines(pair[1])
            records += 2


def warm_input(path: Path) -> None:
    with gzip.open(path, "rb") as handle:
        while handle.read(1024 * 1024):
            pass


def fastq_metrics(path: Path) -> dict[str, int | float]:
    lengths = []
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip()
            plus = handle.readline()
            quality = handle.readline().rstrip()
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("invalid FASTQ record")
            if len(sequence) != len(quality):
                raise ValueError("sequence and quality lengths differ")
            lengths.append(len(sequence))
    if len(lengths) % 2:
        raise ValueError("interleaved FASTQ has an odd record count")
    return {
        "read_pairs": len(lengths) // 2,
        "read_records": len(lengths),
        "read_bases": sum(lengths),
        "read_length_min": min(lengths),
        "read_length_max": max(lengths),
        "read_length_median": statistics.median(lengths),
    }


def parse_elapsed(value: str) -> float:
    fields = value.strip().split(":")
    if len(fields) == 2:
        minutes, seconds = fields
        return int(minutes) * 60 + float(seconds)
    if len(fields) == 3:
        hours, minutes, seconds = fields
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    return float(value)


def parse_time_file(path: Path) -> dict[str, float | int]:
    text = path.read_text()
    if platform.system() == "Darwin":
        match = re.search(r"^\s*([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys", text, re.M)
        rss = re.search(r"^\s*(\d+)\s+maximum resident set size", text, re.M)
        if not match or not rss:
            raise ValueError(f"could not parse macOS time output in {path}")
        return {
            "wall_seconds": float(match.group(1)),
            "user_seconds": float(match.group(2)),
            "sys_seconds": float(match.group(3)),
            "max_rss_bytes": int(rss.group(1)),
        }
    elapsed = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", text)
    user = re.search(r"User time \(seconds\):\s*([0-9.]+)", text)
    system = re.search(r"System time \(seconds\):\s*([0-9.]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    if not elapsed or not user or not system or not rss:
        raise ValueError(f"could not parse GNU time output in {path}")
    return {
        "wall_seconds": parse_elapsed(elapsed.group(1)),
        "user_seconds": float(user.group(1)),
        "sys_seconds": float(system.group(1)),
        "max_rss_bytes": int(rss.group(1)) * 1024,
    }


def time_prefix(path: Path) -> list[str]:
    if platform.system() == "Darwin":
        return ["/usr/bin/time", "-l", "-o", str(path)]
    return ["/usr/bin/time", "-v", "-o", str(path)]


def run_timed(name: str, run: int, command: list[str], output_dir: Path, env: dict[str, str]) -> dict:
    run_dir = output_dir / f"{name}-run{run}"
    run_dir.mkdir(parents=True, exist_ok=False)
    assembly_dir = run_dir / "assembly"
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    time_path = run_dir / "time.txt"
    full_command = command + ["-o", str(assembly_dir)]
    started = time.time()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        result = subprocess.run(
            [*time_prefix(time_path), *full_command],
            stdout=stdout,
            stderr=stderr,
            env=env,
            check=False,
        )
    if result.returncode != 0:
        tail = stderr_path.read_text(errors="replace")[-4000:]
        raise RuntimeError(f"{name} run {run} failed ({result.returncode}):\n{tail}")
    measurement = parse_time_file(time_path)
    measurement.update(
        {
            "assembler": name,
            "run": run,
            "started_unix": started,
            "command": full_command,
            "output_dir": str(assembly_dir),
            "fasta": str(assembly_dir / "final.contigs.fa"),
        }
    )
    return measurement


def fasta_sequences(path: Path):
    sequence: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    yield "".join(sequence).upper()
                    sequence = []
            elif line:
                sequence.append(line)
    if sequence:
        yield "".join(sequence).upper()


def nx(lengths: list[int], fraction: float) -> tuple[int, int]:
    threshold = sum(lengths) * fraction
    cumulative = 0
    for index, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if cumulative >= threshold:
            return length, index
    return 0, 0


def assembly_metrics(path: Path) -> dict:
    sequences = list(fasta_sequences(path))
    lengths = [len(sequence) for sequence in sequences]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if not lengths:
        return {"sha256": digest, "contigs": 0, "total_bp": 0}
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    gc = sum(sequence.count("G") + sequence.count("C") for sequence in sequences)
    known = sum(sum(sequence.count(base) for base in "ACGT") for sequence in sequences)
    return {
        "sha256": digest,
        "contigs": len(lengths),
        "total_bp": sum(lengths),
        "largest_bp": max(lengths),
        "mean_bp": statistics.fmean(lengths),
        "n50_bp": n50,
        "l50": l50,
        "n90_bp": n90,
        "l90": l90,
        "contigs_ge_500": sum(length >= 500 for length in lengths),
        "contigs_ge_1000": sum(length >= 1000 for length in lengths),
        "contigs_ge_5000": sum(length >= 5000 for length in lengths),
        "contigs_ge_10000": sum(length >= 10000 for length in lengths),
        "gc_percent": 100.0 * gc / known if known else 0.0,
        "ambiguous_bases": sum(len(sequence) for sequence in sequences) - known,
    }


def map_reads(
    minimap2: Path,
    samtools: Path,
    assembly: Path,
    read1: Path,
    read2: Path,
    threads: int,
    destination: Path,
    env: dict[str, str],
) -> dict:
    stderr_path = destination.with_suffix(".minimap2.log")
    with stderr_path.open("wb") as stderr:
        mapper = subprocess.Popen(
            [str(minimap2), "-ax", "sr", "-t", str(threads), str(assembly), str(read1), str(read2)],
            stdout=subprocess.PIPE,
            stderr=stderr,
            env=env,
        )
        assert mapper.stdout is not None
        flagstat = subprocess.run(
            [str(samtools), "flagstat", "-@", str(threads), "-O", "json", "-"],
            stdin=mapper.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            check=False,
        )
        mapper.stdout.close()
        mapper_status = mapper.wait()
    if mapper_status != 0 or flagstat.returncode != 0:
        raise RuntimeError(
            f"mapping failed (minimap2={mapper_status}, samtools={flagstat.returncode}): "
            + flagstat.stderr.decode(errors="replace")
        )
    destination.write_bytes(flagstat.stdout)
    report = json.loads(flagstat.stdout)
    passed = report["QC-passed reads"]
    primary = passed["primary"]
    primary_mapped = passed["primary mapped"]
    paired = passed["paired in sequencing"]
    properly_paired = passed["properly paired"]
    return {
        "primary_reads": primary,
        "primary_mapped": primary_mapped,
        "primary_mapped_percent": 100.0 * primary_mapped / primary if primary else 0.0,
        "paired_reads": paired,
        "properly_paired": properly_paired,
        "properly_paired_percent": 100.0 * properly_paired / paired if paired else 0.0,
    }


def median_summary(rows: list[dict]) -> dict:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        grouped.setdefault(row["assembler"], []).append(row)
    summary = {}
    for name, values in grouped.items():
        summary[name] = {
            "runs": len(values),
            "median_wall_seconds": statistics.median(row["wall_seconds"] for row in values),
            "wall_seconds": [row["wall_seconds"] for row in values],
            "median_max_rss_bytes": int(statistics.median(row["max_rss_bytes"] for row in values)),
            "max_rss_bytes": [row["max_rss_bytes"] for row in values],
        }
    return summary


def system_metadata() -> dict[str, str | int | None]:
    cpu_model = platform.processor() or None
    memory_kib = None
    cpuinfo = Path("/proc/cpuinfo")
    meminfo = Path("/proc/meminfo")
    if cpuinfo.exists():
        match = re.search(r"^model name\s*:\s*(.+)$", cpuinfo.read_text(), re.M)
        if match:
            cpu_model = match.group(1)
    if meminfo.exists():
        match = re.search(r"^MemTotal:\s*(\d+)\s+kB$", meminfo.read_text(), re.M)
        if match:
            memory_kib = int(match.group(1))
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_model": cpu_model,
        "logical_cpus": os.cpu_count(),
        "memory_kib": memory_kib,
        "python": platform.python_version(),
        "github_run_id": os.environ.get("GITHUB_RUN_ID"),
        "github_sha": os.environ.get("GITHUB_SHA"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--rust-bin", type=Path, required=True)
    parser.add_argument("--megahit-bin", type=Path, required=True)
    parser.add_argument("--minimap2-bin", type=Path, required=True)
    parser.add_argument("--samtools-bin", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--runs", type=int, default=3)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise SystemExit(f"output directory is not empty: {args.output}")

    read1 = args.output / "reads" / "smalltest_R1.fastq.gz"
    read2 = args.output / "reads" / "smalltest_R2.fastq.gz"
    split_interleaved(args.input.resolve(), read1, read2)
    warm_input(args.input)

    rust_bin = args.rust_bin.resolve()
    megahit = args.megahit_bin.resolve()
    minimap2 = args.minimap2_bin.resolve()
    samtools = args.samtools_bin.resolve()
    env = os.environ.copy()
    tool_dirs = {str(tool.parent) for tool in (megahit, minimap2, samtools)}
    env["PATH"] = os.pathsep.join(sorted(tool_dirs)) + os.pathsep + env.get("PATH", "")

    common = ["--12", str(args.input.resolve()), "--min-count", "2", "--min-contig-len", "200", "-t", str(args.threads)]
    commands = {
        "megahit-rs-k31": [str(rust_bin), *common, "-k", "31"],
        "megahit-k31": [str(megahit), *common, "--k-list", "31"],
        "megahit-default": [str(megahit), *common],
    }
    orders = [
        ["megahit-rs-k31", "megahit-k31", "megahit-default"],
        ["megahit-default", "megahit-rs-k31", "megahit-k31"],
        ["megahit-k31", "megahit-default", "megahit-rs-k31"],
    ]

    measurements = []
    for run in range(1, args.runs + 1):
        order = orders[(run - 1) % len(orders)]
        for name in order:
            print(f"running {name}, replicate {run}", flush=True)
            measurement = run_timed(name, run, commands[name], args.output, env)
            measurements.append(measurement)
            print(
                f"finished {name}, replicate {run}: {measurement['wall_seconds']:.2f}s, "
                f"{measurement['max_rss_bytes'] / (1024**2):.1f} MiB",
                flush=True,
            )

    metrics = {}
    mappings = {}
    deterministic = {}
    for name in commands:
        matching = sorted((row for row in measurements if row["assembler"] == name), key=lambda row: row["run"])
        per_run_metrics = [assembly_metrics(Path(row["fasta"])) for row in matching]
        hashes = [row["sha256"] for row in per_run_metrics]
        deterministic[name] = len(set(hashes)) == 1
        metrics[name] = per_run_metrics[0]
        destination = args.output / f"{name}-flagstat.json"
        mappings[name] = map_reads(
            minimap2,
            samtools,
            Path(matching[0]["fasta"]),
            read1,
            read2,
            args.threads,
            destination,
            env,
        )

    result = {
        "benchmark": {
            "input": str(args.input.resolve()),
            "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
            **fastq_metrics(args.input),
            "threads": args.threads,
            "replicates": args.runs,
            "minimum_contig_length": 200,
            "minimum_kmer_count": 2,
        },
        "versions": {
            "megahit_rs": subprocess.check_output([str(rust_bin), "--version"], text=True).strip(),
            "megahit": subprocess.check_output([str(megahit), "--version"], text=True, env=env).strip(),
            "minimap2": subprocess.check_output([str(minimap2), "--version"], text=True, env=env).strip(),
            "samtools": subprocess.check_output(
                [str(samtools), "--version"], encoding="utf-8", errors="replace", env=env
            ).splitlines()[0],
        },
        "system": system_metadata(),
        "timings": measurements,
        "summary": median_summary(measurements),
        "assembly_metrics": metrics,
        "read_mapping": mappings,
        "byte_deterministic_across_replicates": deterministic,
    }
    (args.output / "benchmark.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"summary": result["summary"], "assembly_metrics": metrics, "read_mapping": mappings}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
