use clap::Parser;
use megahit_rs::output::write_outputs;
use megahit_rs::{InputPaths, PipelineConfig, Result, run_pipeline};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "megahit-rs",
    version,
    about = "Experimental fixed-k metagenome assembler in Rust",
    long_about = "A deterministic fixed-k compacted de Bruijn graph assembler. Version 0.1 accepts MEGAHIT-style read inputs but does not yet implement iterative multi-k assembly or paired-read graph resolution."
)]
struct Cli {
    /// Comma-separated read-1 FASTA/FASTQ files. May be repeated.
    #[arg(short = '1', value_delimiter = ',', value_name = "FILES")]
    paired_1: Vec<PathBuf>,

    /// Comma-separated read-2 FASTA/FASTQ files. May be repeated.
    #[arg(short = '2', value_delimiter = ',', value_name = "FILES")]
    paired_2: Vec<PathBuf>,

    /// Interleaved paired-end FASTA/FASTQ files.
    #[arg(long = "12", value_delimiter = ',', value_name = "FILES")]
    interleaved: Vec<PathBuf>,

    /// Singleton/orphan or other single-end FASTA/FASTQ files.
    #[arg(
        short = 'r',
        long = "read",
        value_delimiter = ',',
        value_name = "FILES"
    )]
    single: Vec<PathBuf>,

    /// Output directory.
    #[arg(short = 'o', long = "out-dir", default_value = "megahit-rs-out")]
    output_dir: PathBuf,

    /// Fixed de Bruijn graph k-mer length (version 0.1 supports 3..=63).
    #[arg(short = 'k', long, default_value_t = 31)]
    k: usize,

    /// Minimum observed coverage for a graph edge.
    #[arg(short = 'm', long, default_value_t = 2)]
    min_count: u32,

    /// Do not emit contigs shorter than this many bases.
    #[arg(long, default_value_t = 200)]
    min_contig_len: usize,

    /// Worker threads used for edge counting.
    #[arg(short = 't', long = "num-cpu-threads", default_value_t = default_threads())]
    threads: usize,

    /// Number of reads retained per parallel counting batch.
    #[arg(long, default_value_t = 20_000, hide = true)]
    batch_size: usize,

    /// Replace generated files if they already exist; unrelated files remain untouched.
    #[arg(long)]
    overwrite: bool,
}

fn default_threads() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    let inputs = InputPaths {
        paired_1: cli.paired_1,
        paired_2: cli.paired_2,
        interleaved: cli.interleaved,
        single: cli.single,
    };
    let config = PipelineConfig {
        k: cli.k,
        min_count: cli.min_count,
        min_contig_len: cli.min_contig_len,
        threads: cli.threads,
        batch_size: cli.batch_size,
    };

    eprintln!(
        "[megahit-rs] v{} fixed-k prototype: k={}, min-count={}, threads={}",
        env!("CARGO_PKG_VERSION"),
        config.k,
        config.min_count,
        config.threads
    );
    let result = run_pipeline(&inputs, &config)?;
    let paths = write_outputs(&cli.output_dir, cli.overwrite, &config, &result)?;
    eprintln!(
        "[megahit-rs] assembled {} contigs ({} bp, N50 {}) in {:.3}s",
        result.assembly.metrics.emitted_contigs,
        result.assembly.metrics.total_assembled_bases,
        result.assembly.metrics.n50,
        result.elapsed_seconds
    );
    eprintln!("[megahit-rs] FASTA: {}", paths.fasta.display());
    eprintln!("[megahit-rs] GFA:   {}", paths.gfa.display());
    eprintln!("[megahit-rs] stats: {}", paths.stats.display());
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("megahit-rs: error: {error}");
        std::process::exit(1);
    }
}
