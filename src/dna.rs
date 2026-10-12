use std::collections::HashMap;

/// A directed edge between two k-mers, represented by the corresponding
/// `(k + 1)`-mer prefix and suffix.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct EdgeKey {
    pub from: u128,
    pub to: u128,
}

#[inline]
pub fn base_bits(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

#[inline]
pub fn bits_base(bits: u8) -> u8 {
    match bits & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

#[inline]
fn mask_for_bases(len: usize) -> u128 {
    let bits = len * 2;
    if bits == 128 {
        u128::MAX
    } else {
        (1_u128 << bits) - 1
    }
}

/// Reverse-complement a packed DNA word containing exactly `len` bases.
pub fn reverse_complement_packed(mut word: u128, len: usize) -> u128 {
    let mut result = 0_u128;
    for _ in 0..len {
        let base = (word & 0b11) as u8;
        result = (result << 2) | u128::from(3 - base);
        word >>= 2;
    }
    result
}

pub fn decode_kmer(mut word: u128, k: usize) -> Vec<u8> {
    let mut sequence = vec![b'A'; k];
    for index in (0..k).rev() {
        sequence[index] = bits_base((word & 0b11) as u8);
        word >>= 2;
    }
    sequence
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base_bits(*base) {
            Some(bits) => bits_base(3 - bits),
            None => b'N',
        })
        .collect()
}

/// Return the lexicographically canonical strand and whether the input is in
/// the same (`+`) orientation as that canonical strand.
pub fn canonical_sequence(sequence: &[u8]) -> (Vec<u8>, bool) {
    let reverse = reverse_complement(sequence);
    if sequence <= reverse.as_slice() {
        (sequence.to_vec(), true)
    } else {
        (reverse, false)
    }
}

/// Canonicalize a circular sequence across both strand and rotation. The
/// returned boolean reports whether the selected representation came from the
/// input strand rather than its reverse complement.
pub fn canonical_circular_sequence(sequence: &[u8]) -> (Vec<u8>, bool) {
    if sequence.is_empty() {
        return (Vec::new(), true);
    }
    let forward = minimal_rotation(sequence);
    let reverse_complemented = reverse_complement(sequence);
    let reverse = minimal_rotation(&reverse_complemented);
    if forward <= reverse {
        (forward, true)
    } else {
        (reverse, false)
    }
}

fn minimal_rotation(sequence: &[u8]) -> Vec<u8> {
    let length = sequence.len();
    if length <= 1 {
        return sequence.to_vec();
    }

    // Booth's algorithm: linear-time lexicographically minimal rotation.
    let mut left = 0_usize;
    let mut right = 1_usize;
    let mut offset = 0_usize;
    while left < length && right < length && offset < length {
        let left_base = sequence[(left + offset) % length];
        let right_base = sequence[(right + offset) % length];
        if left_base == right_base {
            offset += 1;
            continue;
        }
        if left_base > right_base {
            left += offset + 1;
            if left <= right {
                left = right + 1;
            }
        } else {
            right += offset + 1;
            if right <= left {
                right = left + 1;
            }
        }
        offset = 0;
    }
    let start = left.min(right);
    sequence[start..]
        .iter()
        .chain(&sequence[..start])
        .copied()
        .collect()
}

/// Count directed de Bruijn edges in one sequence. Both the observed strand
/// and its reverse complement are represented, making the graph invariant to
/// which DNA strand happened to be sequenced. Palindromic edges are counted
/// once rather than twice.
pub fn count_sequence_edges(sequence: &[u8], k: usize, counts: &mut HashMap<EdgeKey, u32>) {
    debug_assert!((1..=63).contains(&k));

    let window_len = k + 1;
    let window_mask = mask_for_bases(window_len);
    let node_mask = mask_for_bases(k);
    let mut packed = 0_u128;
    let mut valid = 0_usize;

    for base in sequence {
        let Some(bits) = base_bits(*base) else {
            packed = 0;
            valid = 0;
            continue;
        };

        packed = ((packed << 2) | u128::from(bits)) & window_mask;
        valid += 1;
        if valid < window_len {
            continue;
        }

        let forward = EdgeKey {
            from: packed >> 2,
            to: packed & node_mask,
        };
        increment(counts, forward);

        let reverse_word = reverse_complement_packed(packed, window_len);
        let reverse = EdgeKey {
            from: reverse_word >> 2,
            to: reverse_word & node_mask,
        };
        if reverse != forward {
            increment(counts, reverse);
        }
    }
}

#[inline]
fn increment(counts: &mut HashMap<EdgeKey, u32>, edge: EdgeKey) {
    let count = counts.entry(edge).or_insert(0);
    *count = count.saturating_add(1);
}

pub fn merge_counts(target: &mut HashMap<EdgeKey, u32>, source: HashMap<EdgeKey, u32>) {
    for (edge, incoming) in source {
        let count = target.entry(edge).or_insert(0);
        *count = count.saturating_add(incoming);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_reverse_complement_round_trip() {
        let sequence = b"ACGTTGCA";
        let mut packed = 0_u128;
        for base in sequence {
            packed = (packed << 2) | u128::from(base_bits(*base).unwrap());
        }
        let reverse = reverse_complement_packed(packed, sequence.len());
        assert_eq!(decode_kmer(reverse, sequence.len()), b"TGCAACGT");
        assert_eq!(reverse_complement_packed(reverse, sequence.len()), packed);
    }

    #[test]
    fn ambiguous_bases_break_windows() {
        let mut counts = HashMap::new();
        count_sequence_edges(b"AACGTNACGT", 3, &mut counts);
        // ACGT occurs once on each side of N, so each strand-specific edge has
        // coverage two after reverse-complement normalization.
        let edge_word = [b'A', b'C', b'G', b'T']
            .into_iter()
            .fold(0_u128, |word, base| {
                (word << 2) | u128::from(base_bits(base).unwrap())
            });
        let mask = mask_for_bases(3);
        let edge = EdgeKey {
            from: edge_word >> 2,
            to: edge_word & mask,
        };
        assert_eq!(counts.get(&edge), Some(&2));
    }

    #[test]
    fn canonical_sequence_is_strand_invariant() {
        let (forward, _) = canonical_sequence(b"GTCGAA");
        let (reverse, _) = canonical_sequence(b"TTCGAC");
        assert_eq!(forward, reverse);
    }

    #[test]
    fn circular_canonicalization_is_rotation_and_strand_invariant() {
        let (first, _) = canonical_circular_sequence(b"GTTACCA");
        let (rotated, _) = canonical_circular_sequence(b"ACCAGTT");
        let reverse = reverse_complement(b"GTTACCA");
        let (reverse, _) = canonical_circular_sequence(&reverse);
        assert_eq!(first, b"AACTGGT");
        assert_eq!(first, rotated);
        assert_eq!(first, reverse);
    }
}
