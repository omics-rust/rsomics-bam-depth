//! `samtools depth` port over the `rsomics-bamio` raw-record read path.
//!
//! depth only needs each read's refID, start, and reference span (plus MAPQ for
//! the `-Q` filter and the FLAG for the default read filter). It never needs the
//! decoded seq/qual/cigar as noodles types, so iterating `bam::Record` — which
//! fully decodes and clones every record — is pure waste on the hot loop. Each
//! record is read as a raw payload via [`rsomics_bamio::raw`] and its fields are
//! read at fixed offsets; the packed CIGAR is summed in place for the reference
//! span. No decode/clone per record.
//!
//! Per-base coverage matches `samtools depth` (bam_plcmd.c): each read
//! contributes +1 to every reference position its aligned bases (CIGAR M/=/X)
//! cover; deletions and reference skips (D/N) advance the reference cursor
//! without contributing; inserts and clips do not touch the reference. Coverage
//! is accumulated as (start,+1)/(end,-1) interval events per contiguous aligned
//! run, then swept once per chromosome — O(runs·log runs) rather than
//! O(covered bases).

#![allow(clippy::cast_precision_loss)]

use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

// CIGAR op codes (BAM packed encoding, low nibble): M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8.
const CIGAR_MATCH: u8 = 0;
const CIGAR_DELETION: u8 = 2;
const CIGAR_SKIP: u8 = 3;
const CIGAR_SEQ_MATCH: u8 = 7;
const CIGAR_SEQ_MISMATCH: u8 = 8;

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
    let mut record = RawRecord::default();

    while raw::read_record(reader.get_mut(), &mut record)? != 0 {
        if (record.flags() & opts.skip_flags) != 0 {
            continue;
        }
        if record.mapping_quality() < opts.min_mapq {
            continue;
        }

        let tid = record.reference_sequence_id();
        let pos0 = record.alignment_start();
        if tid < 0 || pos0 < 0 {
            continue;
        }
        let t = tid as usize;
        // noodles' `alignment_start` is 1-based; the raw field is 0-based. depth's
        // output column is 1-based, so the run starts at the 0-based pos + 1.
        let start = pos0 as usize + 1;

        let chrom_events = events.entry(t).or_default();
        let mut ref_pos = start;
        for (kind, len) in record.cigar_ops() {
            let len = len as usize;
            match kind {
                CIGAR_MATCH | CIGAR_SEQ_MATCH | CIGAR_SEQ_MISMATCH => {
                    chrom_events.push((ref_pos, 1));
                    chrom_events.push((ref_pos + len, -1));
                    ref_pos += len;
                }
                CIGAR_DELETION | CIGAR_SKIP => {
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

    let mut line: Vec<u8> = Vec::with_capacity(64);
    let mut ib = itoa::Buffer::new();

    for t in tids {
        let name = ref_names.get(t).map_or("*", |n| n.as_str());
        let name_b = name.as_bytes();
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
                let rep = itoa::Buffer::new().format(reported).as_bytes().to_vec();
                for pos in p..until {
                    line.clear();
                    line.extend_from_slice(name_b);
                    line.push(b'\t');
                    line.extend_from_slice(ib.format(pos).as_bytes());
                    line.push(b'\t');
                    line.extend_from_slice(&rep);
                    line.push(b'\n');
                    out.write_all(&line).map_err(RsomicsError::Io)?;
                    lines += 1;
                }
            }
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(lines)
}
