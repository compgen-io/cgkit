# Iterative Consensus MSA for ONT UMI Reads

## Context

We need a multiple sequence alignment algorithm for building consensus sequences from ONT reads sharing the same UMI. Typical input: 1-100 highly similar RNA sequences, up to 10kb each. The dominant error mode is homopolymer indels. No genome reference is used (RNA splicing makes reference coordinates unusable). Future extensions include HP-compressed alignment and consensus calling.

## Approach: Incremental Consensus Alignment

All alignment steps use sequence-to-sequence pairwise alignment (reusing the existing aligner). No profile-profile alignment, no guide trees.

### Step 1 — All-pairs pairwise distances
Full Smith-Waterman for all N*(N-1)/2 pairs (parallelizable).

### Step 2 — Seed selection
Pick the pair with the **highest alignment score**. Tiebreak: longest aligned length. Still tied: first pair found. Align the seed pair → initial 2-sequence MSA.

### Step 3 — Incremental incorporation
Loop until all sequences are incorporated:
1. Compute consensus of current MSA (majority vote per column, skip gaps)
2. Align all remaining sequences to the consensus
3. Pick the best-scoring one, add it to the MSA
4. Repeat

Each addition improves the consensus, so each subsequent sequence gets aligned to an increasingly accurate target. Early additions are the most similar sequences anyway, so their alignments are robust even with a weaker initial consensus.

## New Files

| File | Purpose |
|------|---------|
| `align/msa.go` | Profile type, consensus, incremental alignment, top-level `MSA()` |
| `align/msa_test.go` | Tests for profile, consensus, and end-to-end MSA |
| `internal/cmd/seqcmd/msa.go` | CLI command `seq-msa` |

## Modified Files

| File | Change |
|------|--------|
| `internal/cmd/seqcmd/seqcmd.go` | Register `msaCmd` in `InitCmd()` |

## Implementation Steps

### Step 1: Profile Data Structure (`align/msa.go`)

Column-oriented storage — each column holds one base per sequence (or gap `-`).

```go
type ProfileColumn struct {
    Bases    []byte  // len == number of sequences
    HPRunLen []int   // nil for now; future HP-compressed mode
}

type Profile struct {
    Names   []string
    Columns []ProfileColumn
    NumSeqs int
}
```

- `NewProfileFromSeq(sq SeqQual) *Profile` — single-sequence profile
- `GappedSequences() []string` — extract gapped strings for output
- `Consensus() string` — majority vote per column (skip gaps, break ties alphabetically)

**Why column-oriented**: Makes consensus computation natural (iterate columns), gap insertion straightforward (insert an all-gap column), and is ready for future HP run-length metadata per column.

### Step 2: Seed Pair Selection (`align/msa.go`)

```go
func selectSeedPair(alignments [][]* PairwiseAlignment) (int, int)
```

From the all-pairs alignment results, pick the pair with highest score. Tiebreak by longest aligned length (expanded CIGAR length excluding soft clips). Still tied: first pair found.

### Step 3: Pairwise-to-Profile Conversion (`align/msa.go`)

```go
func profileFromAlignment(aln *PairwiseAlignment) *Profile
```

Convert a `PairwiseAlignment` (with its CIGAR) into a 2-sequence `Profile`. Walk the expanded CIGAR:
- `M`: both bases go in the same column
- `I`: read base in column, gap for the other sequence
- `D`: gap for the read, other sequence's base in column

### Step 4: Add Sequence to Profile (`align/msa.go`)

```go
func (p *Profile) AddSequence(seq seqio.SeqQual, aln *PairwiseAlignment) *Profile
```

Given an alignment of a new sequence against the profile's consensus, merge it into the profile. Walk the CIGAR to determine gap placement:
- `M`: new base appended to existing column
- `I`: insert a new all-gap column for existing sequences + the new base
- `D`: gap for the new sequence at existing column

### Step 5: Distance Matrix & Consensus (`align/msa.go`)

```go
func computeDistances(seqs []seqio.SeqQual, aligner PairwiseAligner, sem *utils.Semaphore) ([][]float32, [][]*PairwiseAlignment)
```

Store both distances and full alignment results (we need the alignments for the seed pair conversion and potentially for incorporation ordering).

**Distance metric**: `1.0 - (matches / alignedLength)` from `PairwiseAlignment.Matches()`.

### Step 6: Top-Level MSA Function (`align/msa.go`)

```go
type MSAOptions struct {
    AlignmentOpts *alignmentOptions
    MaxWorkers    int
}

func MSA(seqs []seqio.SeqQual, opts *MSAOptions) *Profile
```

Flow:
- 0 seqs → nil
- 1 seq → trivial single-sequence profile
- 2 seqs → single pairwise alignment converted to profile
- N seqs:
  1. Compute all-pairs pairwise alignments (parallel)
  2. Select seed pair (highest score, tiebreak longest)
  3. Build initial 2-seq profile from seed pair alignment
  4. Loop: consensus → align remaining → pick best → add → repeat

### Step 7: CLI Command (`internal/cmd/seqcmd/msa.go`)

```
cgkit seq-msa input.fasta [flags]
```

**Flags**:
- `--ont` / `--illumina` — preset selection (default: ONT)
- `--threads` / `-t` — max parallel workers (default: 1)
- `--output` / `-o` — output file (default: stdout)
- `--consensus` — output a single consensus sequence instead of the full MSA
- `--verbose` / `-v` — debug output

**Input**: FASTA or FASTQ file (auto-detected), all sequences loaded into memory. If FASTQ, assumes all reads in the file belong to the same alignment group.

**Output**: Either multi-sequence gapped FASTA (the MSA itself) or a consensus-called FASTA/FASTQ sequence (single record). Controlled by `--consensus`.

### Reused Existing Code

- `align.PairwiseAligner` / `NewLocalAligner()` / `NewGlobalAligner()` — all alignment steps
- `align.OntAlignmentDefaults()` — scoring presets for ONT reads
- `align.alignmentOptions` builder pattern — configuring the aligner
- `align.PairwiseAlignment.Matches()` — for distance computation
- `align.CigarExpand()` — for CIGAR walking during profile construction
- `utils.NewSemaphore()` — parallelism control
- `seqio.NewFastaFile()` / `seqio.NewFastqFile()` — reading input

### Future HP-Compressed Extension Points

No structural changes needed — `ProfileColumn.HPRunLen` is already there:

1. `NewProfileFromSeqHP()` — compress seq, store run lengths per column
2. `Consensus()` — majority-vote base + median/mode HP run length, expand back
3. HP-aware distance metric — compare compressed sequences for distance matrix

### Verification

**Unit tests** (`align/msa_test.go`):
- Profile construction and `GappedSequences()` round-trip
- `Consensus()` — majority vote, tie-breaking, gap handling
- Seed pair selection — verify highest score wins, tiebreak by length
- `profileFromAlignment()` — verify CIGAR→profile conversion
- `AddSequence()` — verify gap propagation when adding to existing profile
- End-to-end MSA with 4-5 sequences with known mutations
- Homopolymer-heavy input with ONT defaults

**Manual testing**:
- Run `seq-msa` on a small FASTA of similar ONT reads
- Compare MSA output with and without `--consensus`
- Verify incorporation order with `--verbose`
