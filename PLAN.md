# CGLTK — Project Plan

## htsio: SAM/BAM/CRAM reader

### Phase 1: samtools-based reader (current, branch: `htsio`)

Shell out to `samtools view` to read records.

- [x] `SamRecord` struct with mandatory SAM fields and flag helpers
- [x] `SamReader` — stream whole file via `samtools view`
- [x] `SamRegionReader` — stream a `chrom:start-end` region
- [x] Unit tests for SAM line parsing and flags
- [ ] Integration test with a small test BAM file
- [ ] Header parsing (read `@SQ`, `@RG`, etc.)
- [ ] Tag parsing (typed access to optional fields like `NM:i:`, `MD:Z:`)

### Phase 2: native bgzip/BAM reader

Replace the samtools dependency with a native Go implementation.

- [ ] bgzip block reader
- [ ] BAM binary record parsing
- [ ] BAI index reading for region queries
- [ ] CRAM support (stretch goal)

### Design notes

- Keep `SamRecord` and the reader interface stable across phases so the native reader is a drop-in replacement.
- The samtools reader starts lazily on first `Next()` call.
- Scanner buffer supports up to 10MB lines for long-read data.

## seqio

- [ ] FAIDX indexed FASTA reader (noted as TODO in `seqio/fasta.go`)

## ont-primers

- No known open items.

## align

- No known open items.
