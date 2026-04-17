//! Integration test for `--quantMode TranscriptomeSAM`.
//!
//! Builds a tiny genome + 2-transcript GTF + a small FASTQ, runs ruSTAR
//! with `--quantMode TranscriptomeSAM`, and asserts that
//! `Aligned.toTranscriptome.out.bam` is produced, is a valid BAM file,
//! and contains at least one record.  Acts as a smoke test for the
//! end-to-end pipeline.

use assert_cmd::Command;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use tempfile::TempDir;

fn generate_genome_seq(seed: u32, length: usize) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    let mut state = seed;
    let mut seq = String::with_capacity(length);
    for _ in 0..length {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(bases[((state >> 16) & 3) as usize]);
    }
    seq
}

fn create_fasta(dir: &TempDir) -> (PathBuf, String) {
    let fasta_path = dir.path().join("genome.fa");
    let mut file = fs::File::create(&fasta_path).unwrap();
    // Single chromosome, 2000 bp pseudo-random content.
    let chr1_seq = generate_genome_seq(98765, 2000);
    writeln!(file, ">chr1").unwrap();
    writeln!(file, "{}", chr1_seq).unwrap();
    (fasta_path, chr1_seq)
}

fn create_gtf(dir: &TempDir) -> PathBuf {
    let gtf_path = dir.path().join("annotations.gtf");
    let mut file = fs::File::create(&gtf_path).unwrap();

    // Transcript T1: one exon [101, 400] forward (1-based inclusive)
    writeln!(
        file,
        "chr1\ttest\texon\t101\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";"
    )
    .unwrap();
    // Transcript T2: one exon [601, 900] forward
    writeln!(
        file,
        "chr1\ttest\texon\t601\t900\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\";"
    )
    .unwrap();

    gtf_path
}

fn create_fastq(dir: &TempDir, n_reads: usize, chr1_seq: &str) -> PathBuf {
    let fastq_path = dir.path().join("reads.fq");
    let mut file = fs::File::create(&fastq_path).unwrap();

    // Alternate reads from T1 region and T2 region.
    for i in 0..n_reads {
        writeln!(file, "@read{}", i + 1).unwrap();
        let start = if i % 2 == 0 { 120 } else { 620 };
        let off = (i * 3) % 180;
        let s = start + off;
        writeln!(file, "{}", &chr1_seq[s..s + 30]).unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII").unwrap();
    }

    fastq_path
}

#[test]
fn transcriptome_sam_end_to_end_smoke_test() {
    let tmpdir = TempDir::new().unwrap();

    let (fasta_path, chr1_seq) = create_fasta(&tmpdir);
    let gtf_path = create_gtf(&tmpdir);
    let fastq_path = create_fastq(&tmpdir, 20, &chr1_seq);

    let genome_dir = tmpdir.path().join("genome");
    let output_dir = tmpdir.path().join("output");
    fs::create_dir_all(&output_dir).unwrap();
    // STAR's prefix semantics: if prefix ends with `/`, it is treated as a
    // directory and files are placed inside it.  ruSTAR's `PathBuf::join`
    // also handles the trailing slash convention.
    let output_prefix = format!("{}/", output_dir.display());

    // Build genome index.
    Command::cargo_bin("ruSTAR")
        .unwrap()
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta_path.to_str().unwrap(),
            "--genomeSAindexNbases",
            "5",
        ])
        .assert()
        .success();

    // Run alignment with --quantMode TranscriptomeSAM.
    Command::cargo_bin("ruSTAR")
        .unwrap()
        .args([
            "--runMode",
            "alignReads",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--readFilesIn",
            fastq_path.to_str().unwrap(),
            "--runThreadN",
            "1",
            "--outFileNamePrefix",
            output_prefix.as_str(),
            "--sjdbGTFfile",
            gtf_path.to_str().unwrap(),
            "--quantMode",
            "TranscriptomeSAM",
            // Permissive mismatch filter so our tiny read set gets through.
            "--outFilterMismatchNmax",
            "20",
            "--outFilterScoreMinOverLread",
            "0.3",
            "--outFilterMatchNminOverLread",
            "0.3",
        ])
        .assert()
        .success();

    // The transcriptome BAM must exist.
    let tr_bam = output_dir.join("Aligned.toTranscriptome.out.bam");
    assert!(
        tr_bam.exists(),
        "Aligned.toTranscriptome.out.bam was not created at {:?}",
        tr_bam
    );

    // File size sanity: BAM header + at least BGZF EOF marker > 100 bytes.
    let meta = fs::metadata(&tr_bam).unwrap();
    assert!(
        meta.len() > 100,
        "transcriptome BAM is suspiciously small: {} bytes",
        meta.len()
    );

    // Validate it's readable as a BAM and check the header has our
    // two @SQ lines (T1 + T2).
    use noodles::bam;
    use std::fs::File;
    let mut reader = bam::io::Reader::new(File::open(&tr_bam).unwrap());
    let header = reader.read_header().expect("valid BAM header");
    let refs = header.reference_sequences();
    assert_eq!(
        refs.len(),
        2,
        "expected 2 @SQ entries (T1, T2), got {}",
        refs.len()
    );
    let keys: Vec<String> = refs
        .keys()
        .map(|k| String::from_utf8_lossy(k).to_string())
        .collect();
    assert!(keys.iter().any(|k| k == "T1"), "T1 missing from @SQ");
    assert!(keys.iter().any(|k| k == "T2"), "T2 missing from @SQ");
}
