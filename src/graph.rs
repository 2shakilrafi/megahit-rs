use crate::dna::{EdgeKey, canonical_circular_sequence, canonical_sequence, decode_kmer};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Clone, Debug)]
struct Edge {
    key: EdgeKey,
    coverage: u32,
}

#[derive(Clone, Debug)]
struct RawUnitig {
    sequence: Vec<u8>,
    start_node: u128,
    end_node: u128,
    mean_coverage: f64,
    min_coverage: u32,
    circular: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct Segment {
    pub id: String,
    #[serde(skip)]
    pub sequence: Vec<u8>,
    pub length: usize,
    pub mean_coverage: f64,
    pub min_coverage: u32,
    pub circular: bool,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Link {
    pub from: String,
    pub from_orientation: char,
    pub to: String,
    pub to_orientation: char,
    pub overlap: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct AssemblyMetrics {
    pub observed_edges: usize,
    pub solid_edges: usize,
    pub graph_nodes: usize,
    pub raw_unitigs: usize,
    pub emitted_contigs: usize,
    pub total_assembled_bases: usize,
    pub longest_contig: usize,
    pub n50: usize,
}

#[derive(Clone, Debug)]
pub struct Assembly {
    pub segments: Vec<Segment>,
    pub links: Vec<Link>,
    pub metrics: AssemblyMetrics,
}

#[derive(Clone, Debug)]
struct SegmentDraft {
    canonical: Vec<u8>,
    mean_coverage: f64,
    min_coverage: u32,
    circular: bool,
}

#[derive(Clone, Debug)]
struct RawMapping {
    canonical: Vec<u8>,
    forward: bool,
}

/// Build and compact a directed fixed-k de Bruijn graph from counted edges.
/// Edges below `min_count` are discarded before topology is computed.
pub fn assemble_graph(
    counts: &HashMap<EdgeKey, u32>,
    k: usize,
    min_count: u32,
    min_contig_len: usize,
) -> Assembly {
    let mut edges: Vec<Edge> = counts
        .iter()
        .filter(|(_, coverage)| **coverage >= min_count)
        .map(|(key, coverage)| Edge {
            key: *key,
            coverage: *coverage,
        })
        .collect();
    edges.sort_unstable_by_key(|edge| edge.key);

    let mut outgoing: HashMap<u128, Vec<usize>> = HashMap::new();
    let mut incoming: HashMap<u128, Vec<usize>> = HashMap::new();
    for (index, edge) in edges.iter().enumerate() {
        outgoing.entry(edge.key.from).or_default().push(index);
        incoming.entry(edge.key.to).or_default().push(index);
    }

    let mut all_nodes = BTreeSet::new();
    for edge in &edges {
        all_nodes.insert(edge.key.from);
        all_nodes.insert(edge.key.to);
    }

    let mut visited = vec![false; edges.len()];
    let mut raw = Vec::new();

    // Paths beginning at a topology boundary (source, sink, or branch).
    for node in &all_nodes {
        let indegree = incoming.get(node).map_or(0, Vec::len);
        let outdegree = outgoing.get(node).map_or(0, Vec::len);
        if indegree == 1 && outdegree == 1 {
            continue;
        }
        if let Some(next_edges) = outgoing.get(node) {
            for edge_index in next_edges {
                if !visited[*edge_index] {
                    raw.push(walk_unitig(
                        *edge_index,
                        &edges,
                        &incoming,
                        &outgoing,
                        &mut visited,
                        k,
                        false,
                    ));
                }
            }
        }
    }

    // Any remaining edge belongs to an isolated cycle.
    for edge_index in 0..edges.len() {
        if !visited[edge_index] {
            raw.push(walk_unitig(
                edge_index,
                &edges,
                &incoming,
                &outgoing,
                &mut visited,
                k,
                true,
            ));
        }
    }

    let mut drafts: BTreeMap<Vec<u8>, SegmentDraft> = BTreeMap::new();
    let mut raw_mappings = Vec::with_capacity(raw.len());
    for unitig in &raw {
        let (canonical, forward) = if unitig.circular {
            canonical_circular_sequence(&unitig.sequence)
        } else {
            canonical_sequence(&unitig.sequence)
        };
        raw_mappings.push(RawMapping {
            canonical: canonical.clone(),
            forward,
        });

        if canonical.len() < min_contig_len {
            continue;
        }
        drafts
            .entry(canonical.clone())
            .and_modify(|draft| {
                // Reverse-complement twins should have identical support. Max
                // makes the result stable even in asymmetric edge cases.
                draft.mean_coverage = draft.mean_coverage.max(unitig.mean_coverage);
                draft.min_coverage = draft.min_coverage.max(unitig.min_coverage);
                draft.circular |= unitig.circular;
            })
            .or_insert(SegmentDraft {
                canonical,
                mean_coverage: unitig.mean_coverage,
                min_coverage: unitig.min_coverage,
                circular: unitig.circular,
            });
    }

    let mut ordered_drafts: Vec<SegmentDraft> = drafts.into_values().collect();
    ordered_drafts.sort_by(|left, right| {
        right
            .canonical
            .len()
            .cmp(&left.canonical.len())
            .then_with(|| left.canonical.cmp(&right.canonical))
    });

    let mut key_to_index: HashMap<Vec<u8>, usize> = HashMap::new();
    let segments: Vec<Segment> = ordered_drafts
        .into_iter()
        .enumerate()
        .map(|(index, draft)| {
            key_to_index.insert(draft.canonical.clone(), index);
            Segment {
                id: format!("ctg{:06}", index + 1),
                length: draft.canonical.len(),
                sequence: draft.canonical,
                mean_coverage: draft.mean_coverage,
                min_coverage: draft.min_coverage,
                circular: draft.circular,
            }
        })
        .collect();

    let mut starts: HashMap<u128, Vec<usize>> = HashMap::new();
    for (index, unitig) in raw.iter().enumerate() {
        if !unitig.circular && key_to_index.contains_key(&raw_mappings[index].canonical) {
            starts.entry(unitig.start_node).or_default().push(index);
        }
    }

    let mut link_set = BTreeSet::new();
    for (from_raw_index, from_unitig) in raw.iter().enumerate() {
        if from_unitig.circular {
            continue;
        }
        let Some(&from_segment_index) = key_to_index.get(&raw_mappings[from_raw_index].canonical)
        else {
            continue;
        };
        let Some(targets) = starts.get(&from_unitig.end_node) else {
            continue;
        };

        for to_raw_index in targets {
            let Some(&to_segment_index) = key_to_index.get(&raw_mappings[*to_raw_index].canonical)
            else {
                continue;
            };
            link_set.insert(Link {
                from: segments[from_segment_index].id.clone(),
                from_orientation: if raw_mappings[from_raw_index].forward {
                    '+'
                } else {
                    '-'
                },
                to: segments[to_segment_index].id.clone(),
                to_orientation: if raw_mappings[*to_raw_index].forward {
                    '+'
                } else {
                    '-'
                },
                overlap: k,
            });
        }
    }

    let lengths: Vec<usize> = segments.iter().map(|segment| segment.length).collect();
    let metrics = AssemblyMetrics {
        observed_edges: counts.len(),
        solid_edges: edges.len(),
        graph_nodes: all_nodes.len(),
        raw_unitigs: raw.len(),
        emitted_contigs: segments.len(),
        total_assembled_bases: lengths.iter().sum(),
        longest_contig: lengths.first().copied().unwrap_or(0),
        n50: n50(&lengths),
    };

    Assembly {
        segments,
        links: link_set.into_iter().collect(),
        metrics,
    }
}

#[allow(clippy::too_many_arguments)]
fn walk_unitig(
    first_edge: usize,
    edges: &[Edge],
    incoming: &HashMap<u128, Vec<usize>>,
    outgoing: &HashMap<u128, Vec<usize>>,
    visited: &mut [bool],
    k: usize,
    expect_cycle: bool,
) -> RawUnitig {
    let start_node = edges[first_edge].key.from;
    let mut sequence = decode_kmer(start_node, k);
    let mut coverage_sum = 0_u64;
    let mut min_coverage = u32::MAX;
    let mut edge_count = 0_usize;
    let mut edge_index = first_edge;
    let mut end_node = start_node;
    let mut circular = false;

    loop {
        if visited[edge_index] {
            break;
        }
        visited[edge_index] = true;
        let edge = &edges[edge_index];
        sequence.push(crate::dna::bits_base((edge.key.to & 0b11) as u8));
        coverage_sum += u64::from(edge.coverage);
        min_coverage = min_coverage.min(edge.coverage);
        edge_count += 1;
        end_node = edge.key.to;

        if expect_cycle && end_node == start_node {
            circular = true;
            break;
        }

        let indegree = incoming.get(&end_node).map_or(0, Vec::len);
        let outdegree = outgoing.get(&end_node).map_or(0, Vec::len);
        if indegree != 1 || outdegree != 1 {
            break;
        }
        let next_edge = outgoing[&end_node][0];
        if visited[next_edge] {
            circular = expect_cycle;
            break;
        }
        edge_index = next_edge;
    }

    if circular && sequence.len() > k {
        // The final k bases repeat the starting node when linearizing a cycle.
        sequence.truncate(sequence.len() - k);
    }

    RawUnitig {
        sequence,
        start_node,
        end_node,
        mean_coverage: coverage_sum as f64 / edge_count.max(1) as f64,
        min_coverage: if min_coverage == u32::MAX {
            0
        } else {
            min_coverage
        },
        circular,
    }
}

fn n50(lengths: &[usize]) -> usize {
    let total: usize = lengths.iter().sum();
    if total == 0 {
        return 0;
    }
    let threshold = total.div_ceil(2);
    let mut accumulated = 0;
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable_by(|left, right| right.cmp(left));
    for length in sorted {
        accumulated += length;
        if accumulated >= threshold {
            return length;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::count_sequence_edges;

    #[test]
    fn compacts_a_linear_graph_and_deduplicates_strands() {
        let mut counts = HashMap::new();
        for _ in 0..3 {
            count_sequence_edges(b"AACCGTTAGCTACGATCGGATC", 7, &mut counts);
        }
        let assembly = assemble_graph(&counts, 7, 2, 1);
        assert_eq!(assembly.segments.len(), 1);
        assert_eq!(assembly.segments[0].length, 22);
        assert_eq!(assembly.segments[0].sequence, b"AACCGTTAGCTACGATCGGATC");
        assert!(assembly.segments[0].mean_coverage >= 3.0);
    }

    #[test]
    fn filters_low_support_edges() {
        let mut counts = HashMap::new();
        count_sequence_edges(b"AACCGTTAGCTACGATCGGATC", 7, &mut counts);
        let assembly = assemble_graph(&counts, 7, 2, 1);
        assert!(assembly.segments.is_empty());
        assert_eq!(assembly.metrics.solid_edges, 0);
    }

    #[test]
    fn n50_handles_empty_and_regular_inputs() {
        assert_eq!(n50(&[]), 0);
        assert_eq!(n50(&[10, 6, 4]), 10);
        assert_eq!(n50(&[8, 7, 6, 5]), 7);
    }

    #[test]
    fn reconstructs_a_sequence_from_overlapping_reads() {
        let genome = pseudo_random_dna(500);
        let mut counts = HashMap::new();
        let read_length = 90;
        let mut starts: Vec<usize> = (0..=genome.len() - read_length).step_by(7).collect();
        if starts.last().copied() != Some(genome.len() - read_length) {
            starts.push(genome.len() - read_length);
        }
        for start in starts {
            count_sequence_edges(&genome[start..start + read_length], 31, &mut counts);
        }

        let assembly = assemble_graph(&counts, 31, 1, 1);
        let (canonical, _) = crate::dna::canonical_sequence(&genome);
        assert!(
            assembly
                .segments
                .iter()
                .any(|segment| segment.sequence == canonical),
            "the full genome was not reconstructed; lengths={:?}",
            assembly
                .segments
                .iter()
                .map(|segment| segment.length)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn extracts_and_deduplicates_a_circular_component() {
        let genome = pseudo_random_dna(101);
        let mut counts = HashMap::new();
        for start in 0..genome.len() {
            let read: Vec<u8> = (0..45)
                .map(|offset| genome[(start + offset) % genome.len()])
                .collect();
            count_sequence_edges(&read, 21, &mut counts);
        }

        let assembly = assemble_graph(&counts, 21, 2, 1);
        assert_eq!(assembly.segments.len(), 1);
        assert_eq!(assembly.segments[0].length, genome.len());
        assert!(assembly.segments[0].circular);
    }

    fn pseudo_random_dna(length: usize) -> Vec<u8> {
        let mut state = 0x6a09_e667_f3bc_c909_u64;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }
}
