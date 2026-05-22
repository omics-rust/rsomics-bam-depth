#![allow(clippy::cast_precision_loss)]

use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Clone)]
pub struct DepthOpts {
    pub min_mapq: u8,
    /// Cap reported depth at this value; 0 = no cap (matches `samtools depth -d 0`).
    pub max_depth: u32,
    pub skip_flags: u16,
}

impl Default for DepthOpts {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            max_depth: 0,
            skip_flags: 0x704,
        }
    }
}

fn tid(r: &bam::Record) -> Option<usize> {
    r.reference_sequence_id().transpose().ok().flatten()
}

fn start(r: &bam::Record) -> Option<usize> {
    r.alignment_start()
        .transpose()
        .ok()
        .flatten()
        .map(|p| p.get())
}

/// Per-base coverage like `samtools depth`: each read contributes +1 to every
/// reference position its aligned bases (CIGAR M/=/X) cover. Deletions and
/// reference skips advance the reference cursor without contributing; inserts
/// and clips do not touch the reference. Coverage is accumulated as
/// (start,+1)/(end,-1) interval events per contiguous aligned run, then swept
/// once per chromosome — O(runs·log runs) instead of O(covered bases).
pub fn compute_depth(
    input: &Path,
    output: &mut dyn Write,
    opts: &DepthOpts,
    workers: NonZero<usize>,
) -> Result<u64> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();

    let mut events: HashMap<usize, Vec<(usize, i64)>> = HashMap::new();

    for result in reader.records() {
        let record = result.map_err(RsomicsError::Io)?;
        let flags = record.flags();

        if (flags.bits() & opts.skip_flags) != 0 {
            continue;
        }

        let mq = record.mapping_quality().map_or(0, |q| q.get());
        if mq < opts.min_mapq {
            continue;
        }

        let Some(t) = tid(&record) else { continue };
        let Some(s) = start(&record) else { continue };

        let chrom_events = events.entry(t).or_default();
        let mut ref_pos = s;
        for op in record.cigar().iter() {
            let op = op.map_err(RsomicsError::Io)?;
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    chrom_events.push((ref_pos, 1));
                    chrom_events.push((ref_pos + len, -1));
                    ref_pos += len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len;
                }
                _ => {}
            }
        }
    }

    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut lines: u64 = 0;

    let mut tids: Vec<usize> = events.keys().copied().collect();
    tids.sort_unstable();

    for t in tids {
        let name = ref_names.get(t).map_or("*", |n| n.as_str());
        let evs = events.get_mut(&t).unwrap();
        evs.sort_unstable();

        let mut depth: i64 = 0;
        let mut i = 0;
        while i < evs.len() {
            let p = evs[i].0;
            while i < evs.len() && evs[i].0 == p {
                depth += evs[i].1;
                i += 1;
            }
            if depth > 0 {
                let until = if i < evs.len() { evs[i].0 } else { p };
                let reported = if opts.max_depth > 0 {
                    depth.min(i64::from(opts.max_depth))
                } else {
                    depth
                };
                for pos in p..until {
                    writeln!(out, "{name}\t{pos}\t{reported}").map_err(RsomicsError::Io)?;
                    lines += 1;
                }
            }
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(lines)
}
