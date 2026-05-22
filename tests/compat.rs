use std::path::Path;
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam-depth"))
}

fn fixture() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/small.bam"))
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn ours(args: &[&str]) -> String {
    let out = bin().args(args).arg(fixture()).output().unwrap();
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn samtools_depth(args: &[&str]) -> String {
    let out = Command::new("samtools")
        .arg("depth")
        .args(args)
        .arg(fixture())
        .output()
        .unwrap();
    assert!(out.status.success());
    String::from_utf8(out.stdout).unwrap()
}

// Per-base depth across each read's CIGAR span must match `samtools depth`
// exactly on coordinate-sorted input (the golden fixture is sorted).
#[test]
fn depth_matches_samtools() {
    if !samtools_available() {
        eprintln!("skipping: samtools not found");
        return;
    }
    assert_eq!(ours(&[]), samtools_depth(&[]));
}

#[test]
fn min_mapq_matches_samtools() {
    if !samtools_available() {
        eprintln!("skipping: samtools not found");
        return;
    }
    // samtools depth -Q filters by mapping quality; ours --min-mapq.
    for q in ["1", "30", "60"] {
        assert_eq!(
            ours(&["--min-mapq", q]),
            samtools_depth(&["-Q", q]),
            "min-mapq {q}"
        );
    }
}
