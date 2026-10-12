use bzip2::Compression as BzCompression;
use bzip2::write::BzEncoder;
use flate2::Compression as GzCompression;
use flate2::write::GzEncoder;
use megahit_rs::output::write_outputs;
use megahit_rs::{InputPaths, PipelineConfig, run_pipeline};
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

static NEXT_TEMP_DIRECTORY: AtomicU64 = AtomicU64::new(0);

fn temporary_directory() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let sequence = NEXT_TEMP_DIRECTORY.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!(
        "megahit-rs-test-{}-{nonce}-{sequence}",
        std::process::id()
    ));
    fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn assembles_fastq_and_writes_all_outputs() {
    let temp = temporary_directory();
    let reads = temp.join("reads.fastq");
    let mut writer = fs::File::create(&reads).unwrap();
    for index in 0..4 {
        writeln!(writer, "@read{index}").unwrap();
        writeln!(writer, "AACCGTTAGCTACGATCGGATC").unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "IIIIIIIIIIIIIIIIIIIIII").unwrap();
    }
    drop(writer);

    let config = PipelineConfig {
        k: 7,
        min_count: 2,
        min_contig_len: 1,
        threads: 2,
        batch_size: 2,
    };
    let result = run_pipeline(
        &InputPaths {
            single: vec![reads],
            ..InputPaths::default()
        },
        &config,
    )
    .unwrap();
    assert_eq!(result.assembly.segments.len(), 1);
    assert_eq!(result.assembly.segments[0].length, 22);

    let output = temp.join("out");
    let paths = write_outputs(&output, false, &config, &result).unwrap();
    assert!(paths.fasta.exists());
    assert!(paths.gfa.exists());
    assert!(paths.stats.exists());
    assert!(
        fs::read_to_string(paths.fasta)
            .unwrap()
            .contains("AACCGTTAGCTACGATCGGATC")
    );
    assert!(
        fs::read_to_string(paths.gfa)
            .unwrap()
            .starts_with("H\tVN:Z:1.0")
    );

    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn compressed_inputs_are_deterministic_across_thread_counts() {
    let temp = temporary_directory();
    let fastq = repeated_fastq();
    let gzip_path = temp.join("reads.fastq.gz");
    let bzip_path = temp.join("reads.fastq.bz2");

    let mut gzip = GzEncoder::new(
        fs::File::create(&gzip_path).unwrap(),
        GzCompression::default(),
    );
    gzip.write_all(fastq.as_bytes()).unwrap();
    gzip.finish().unwrap();

    let mut bzip = BzEncoder::new(
        fs::File::create(&bzip_path).unwrap(),
        BzCompression::default(),
    );
    bzip.write_all(fastq.as_bytes()).unwrap();
    bzip.finish().unwrap();

    let single_thread = run_pipeline(
        &InputPaths {
            single: vec![gzip_path],
            ..InputPaths::default()
        },
        &PipelineConfig {
            k: 7,
            min_count: 2,
            min_contig_len: 1,
            threads: 1,
            batch_size: 2,
        },
    )
    .unwrap();
    let four_threads = run_pipeline(
        &InputPaths {
            single: vec![bzip_path],
            ..InputPaths::default()
        },
        &PipelineConfig {
            k: 7,
            min_count: 2,
            min_contig_len: 1,
            threads: 4,
            batch_size: 3,
        },
    )
    .unwrap();

    let first_sequences: Vec<&[u8]> = single_thread
        .assembly
        .segments
        .iter()
        .map(|segment| segment.sequence.as_slice())
        .collect();
    let second_sequences: Vec<&[u8]> = four_threads
        .assembly
        .segments
        .iter()
        .map(|segment| segment.sequence.as_slice())
        .collect();
    assert_eq!(first_sequences, second_sequences);
    assert_eq!(single_thread.assembly.links, four_threads.assembly.links);
    assert_eq!(
        single_thread.assembly.metrics.total_assembled_bases,
        four_threads.assembly.metrics.total_assembled_bases
    );

    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn rejects_mismatched_paired_record_counts() {
    let temp = temporary_directory();
    let read_1 = temp.join("R1.fastq");
    let read_2 = temp.join("R2.fastq");
    fs::write(&read_1, repeated_fastq()).unwrap();
    fs::write(
        &read_2,
        "@only_one\nAACCGTTAGCTACGATCGGATC\n+\nIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();

    let error = run_pipeline(
        &InputPaths {
            paired_1: vec![read_1],
            paired_2: vec![read_2],
            ..InputPaths::default()
        },
        &PipelineConfig {
            k: 7,
            min_count: 1,
            min_contig_len: 1,
            threads: 2,
            batch_size: 2,
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("different record counts"));
    fs::remove_dir_all(temp).unwrap();
}

fn repeated_fastq() -> String {
    (0..4)
        .map(|index| format!("@read{index}\nAACCGTTAGCTACGATCGGATC\n+\nIIIIIIIIIIIIIIIIIIIIII\n"))
        .collect()
}
