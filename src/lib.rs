//! Core library for the experimental `megahit-rs` assembler.
//!
//! Version 0.1 deliberately implements a small, testable slice of a modern
//! metagenome assembler: fixed-k edge counting, solid-edge filtering,
//! compacted de Bruijn graph traversal, and deterministic FASTA/GFA output.

pub mod dna;
pub mod graph;
pub mod output;
pub mod pipeline;

pub use graph::{Assembly, AssemblyMetrics, Link, Segment, assemble_graph};
pub use pipeline::{InputPaths, PipelineConfig, PipelineResult, ReadMetrics, run_pipeline};

use std::error::Error;

pub type Result<T> = std::result::Result<T, Box<dyn Error + Send + Sync>>;
