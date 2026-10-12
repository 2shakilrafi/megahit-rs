use crate::Result;
use crate::dna::{EdgeKey, base_bits, count_sequence_edges, merge_counts};
use crate::graph::{Assembly, assemble_graph};
use needletail::parse_fastx_file;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Debug, Default)]
pub struct InputPaths {
    pub paired_1: Vec<PathBuf>,
    pub paired_2: Vec<PathBuf>,
    pub interleaved: Vec<PathBuf>,
    pub single: Vec<PathBuf>,
}

#[derive(Clone, Debug)]
pub struct PipelineConfig {
    pub k: usize,
    pub min_count: u32,
    pub min_contig_len: usize,
    pub threads: usize,
    pub batch_size: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ReadMetrics {
    pub input_files: usize,
    pub reads: u64,
    pub bases: u64,
    pub ambiguous_bases: u64,
    pub reads_shorter_than_k_plus_one: u64,
}

#[derive(Clone, Debug)]
pub struct PipelineResult {
    pub assembly: Assembly,
    pub reads: ReadMetrics,
    pub elapsed_seconds: f64,
}

pub fn run_pipeline(inputs: &InputPaths, config: &PipelineConfig) -> Result<PipelineResult> {
    validate(inputs, config)?;
    let started = Instant::now();
    let pool = ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .thread_name(|index| format!("megahit-rs-{index}"))
        .build()?;

    let mut counts: HashMap<EdgeKey, u32> = HashMap::new();
    let mut metrics = ReadMetrics::default();

    for (read_1, read_2) in inputs.paired_1.iter().zip(&inputs.paired_2) {
        let count_1 = process_file(read_1, config, &pool, &mut counts, &mut metrics)?;
        let count_2 = process_file(read_2, config, &pool, &mut counts, &mut metrics)?;
        if count_1 != count_2 {
            return Err(format!(
                "paired files contain different record counts: {} has {count_1}, {} has {count_2}",
                read_1.display(),
                read_2.display()
            )
            .into());
        }
    }

    for path in &inputs.interleaved {
        let records = process_file(path, config, &pool, &mut counts, &mut metrics)?;
        if records % 2 != 0 {
            return Err(format!(
                "interleaved input {} contains an odd number of records ({records})",
                path.display()
            )
            .into());
        }
    }

    for path in &inputs.single {
        process_file(path, config, &pool, &mut counts, &mut metrics)?;
    }

    if counts.is_empty() {
        return Err(format!(
            "no (k+1)-mers were observed at k={}; use a smaller k or verify the input",
            config.k
        )
        .into());
    }

    let assembly = assemble_graph(&counts, config.k, config.min_count, config.min_contig_len);
    if assembly.metrics.solid_edges == 0 {
        return Err(format!(
            "no edges survived --min-count {}; lower the threshold or increase sequencing depth",
            config.min_count
        )
        .into());
    }

    Ok(PipelineResult {
        assembly,
        reads: metrics,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    })
}

fn validate(inputs: &InputPaths, config: &PipelineConfig) -> Result<()> {
    if !(3..=63).contains(&config.k) {
        return Err("k must be between 3 and 63 in version 0.1".into());
    }
    if config.min_count == 0 {
        return Err("--min-count must be at least 1".into());
    }
    if config.threads == 0 {
        return Err("--num-cpu-threads must be at least 1".into());
    }
    if config.batch_size == 0 {
        return Err("--batch-size must be at least 1".into());
    }
    if inputs.paired_1.len() != inputs.paired_2.len() {
        return Err(format!(
            "-1 and -2 must contain the same number of files ({} versus {})",
            inputs.paired_1.len(),
            inputs.paired_2.len()
        )
        .into());
    }
    if inputs.paired_1.is_empty() && inputs.interleaved.is_empty() && inputs.single.is_empty() {
        return Err("provide paired reads with -1/-2, interleaved reads with --12, or singleton reads with -r".into());
    }
    Ok(())
}

fn process_file(
    path: &Path,
    config: &PipelineConfig,
    pool: &rayon::ThreadPool,
    global_counts: &mut HashMap<EdgeKey, u32>,
    metrics: &mut ReadMetrics,
) -> Result<u64> {
    eprintln!("[megahit-rs] reading {}", path.display());
    let mut reader = parse_fastx_file(path)
        .map_err(|error| format!("cannot open or parse {}: {error}", path.display()))?;
    let mut batch: Vec<Vec<u8>> = Vec::with_capacity(config.batch_size);
    let mut file_records = 0_u64;
    metrics.input_files += 1;

    while let Some(record) = reader.next() {
        let record = record.map_err(|error| {
            format!("invalid FASTA/FASTQ record in {}: {error}", path.display())
        })?;
        let sequence = record.seq().into_owned();
        file_records += 1;
        metrics.reads += 1;
        metrics.bases += sequence.len() as u64;
        metrics.ambiguous_bases += sequence
            .iter()
            .filter(|base| base_bits(**base).is_none())
            .count() as u64;
        if sequence.len() < config.k + 1 {
            metrics.reads_shorter_than_k_plus_one += 1;
        }
        batch.push(sequence);

        if batch.len() >= config.batch_size {
            merge_batch(global_counts, &batch, config.k, pool);
            batch.clear();
        }
    }

    if !batch.is_empty() {
        merge_batch(global_counts, &batch, config.k, pool);
    }
    eprintln!(
        "[megahit-rs] read {file_records} records from {}",
        path.display()
    );
    Ok(file_records)
}

fn merge_batch(
    global_counts: &mut HashMap<EdgeKey, u32>,
    batch: &[Vec<u8>],
    k: usize,
    pool: &rayon::ThreadPool,
) {
    let batch_counts = pool.install(|| {
        batch
            .par_iter()
            .fold(HashMap::new, |mut local, sequence| {
                count_sequence_edges(sequence, k, &mut local);
                local
            })
            .reduce(HashMap::new, |mut left, right| {
                merge_counts(&mut left, right);
                left
            })
    });
    merge_counts(global_counts, batch_counts);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validation_rejects_unpaired_file_lists() {
        let inputs = InputPaths {
            paired_1: vec![PathBuf::from("r1.fastq")],
            ..InputPaths::default()
        };
        let config = PipelineConfig {
            k: 31,
            min_count: 2,
            min_contig_len: 200,
            threads: 1,
            batch_size: 100,
        };
        assert!(validate(&inputs, &config).is_err());
    }
}
