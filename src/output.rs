use crate::Result;
use crate::pipeline::{PipelineConfig, PipelineResult};
use serde::Serialize;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub struct OutputPaths {
    pub fasta: PathBuf,
    pub gfa: PathBuf,
    pub stats: PathBuf,
}

#[derive(Serialize)]
struct RunReport<'a> {
    program: &'static str,
    version: &'static str,
    status: &'static str,
    algorithm: &'static str,
    paired_constraints_used: bool,
    k: usize,
    min_count: u32,
    min_contig_len: usize,
    threads: usize,
    elapsed_seconds: f64,
    reads: &'a crate::pipeline::ReadMetrics,
    assembly: &'a crate::graph::AssemblyMetrics,
}

pub fn write_outputs(
    output_dir: &Path,
    overwrite: bool,
    config: &PipelineConfig,
    result: &PipelineResult,
) -> Result<OutputPaths> {
    fs::create_dir_all(output_dir)?;
    let paths = OutputPaths {
        fasta: output_dir.join("final.contigs.fa"),
        gfa: output_dir.join("final.graph.gfa"),
        stats: output_dir.join("assembly.stats.json"),
    };

    if !overwrite {
        for path in [&paths.fasta, &paths.gfa, &paths.stats] {
            if path.exists() {
                return Err(format!(
                    "refusing to overwrite {}; pass --overwrite to replace generated files",
                    path.display()
                )
                .into());
            }
        }
    }

    write_fasta(&paths.fasta, result)?;
    write_gfa(&paths.gfa, config, result)?;
    write_stats(&paths.stats, config, result)?;
    Ok(paths)
}

fn write_fasta(path: &Path, result: &PipelineResult) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for segment in &result.assembly.segments {
        writeln!(
            writer,
            ">{} len={} mean_cov={:.2} min_cov={} circular={}",
            segment.id,
            segment.length,
            segment.mean_coverage,
            segment.min_coverage,
            if segment.circular { "yes" } else { "no" }
        )?;
        for line in segment.sequence.chunks(80) {
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_gfa(path: &Path, config: &PipelineConfig, result: &PipelineResult) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "H\tVN:Z:1.0\tKL:i:{}", config.k)?;
    for segment in &result.assembly.segments {
        writeln!(
            writer,
            "S\t{}\t{}\tLN:i:{}\tDP:f:{:.3}\tMC:i:{}\tTP:Z:{}",
            segment.id,
            String::from_utf8_lossy(&segment.sequence),
            segment.length,
            segment.mean_coverage,
            segment.min_coverage,
            if segment.circular {
                "circular"
            } else {
                "linear"
            }
        )?;
    }
    for link in &result.assembly.links {
        writeln!(
            writer,
            "L\t{}\t{}\t{}\t{}\t{}M",
            link.from, link.from_orientation, link.to, link.to_orientation, link.overlap
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_stats(path: &Path, config: &PipelineConfig, result: &PipelineResult) -> Result<()> {
    let report = RunReport {
        program: "megahit-rs",
        version: env!("CARGO_PKG_VERSION"),
        status: "experimental-fixed-k",
        algorithm: "packed fixed-k compacted de Bruijn graph",
        paired_constraints_used: false,
        k: config.k,
        min_count: config.min_count,
        min_contig_len: config.min_contig_len,
        threads: config.threads,
        elapsed_seconds: result.elapsed_seconds,
        reads: &result.reads,
        assembly: &result.assembly.metrics,
    };
    let mut writer = BufWriter::new(File::create(path)?);
    serde_json::to_writer_pretty(&mut writer, &report)?;
    writer.write_all(b"\n")?;
    writer.flush()?;
    Ok(())
}
