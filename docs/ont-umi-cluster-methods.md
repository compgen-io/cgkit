# UMI Clustering for Oxford Nanopore Sequencing: Methods for Error-Tolerant Molecular Deduplication

## Abstract

Unique Molecular Identifiers (UMIs) enable accurate quantification of original molecules in sequencing libraries by tagging each molecule before amplification. However, Oxford Nanopore Technology (ONT) sequencing introduces systematic errors, particularly in homopolymer regions, that complicate UMI-based deduplication. We present a suite of clustering methods for collapsing UMI sequences in ONT data, including a novel homopolymer-aware edit distance metric, an adaptive statistical threshold for large components, and splice junction-aware read grouping for RNA applications. These methods are implemented in the `ont-umi-cluster` command of the CGLTK toolkit.

## Introduction

### The UMI deduplication problem

PCR amplification during library preparation creates duplicate copies of original molecules. UMIs --- short random sequences ligated to molecules before amplification --- allow computational identification of duplicates: reads sharing the same UMI and mapping to the same genomic locus are presumed to originate from the same pre-amplification molecule. Accurate UMI deduplication is critical for quantitative applications such as gene expression analysis and variant calling at low allele frequencies.

### Challenges with ONT UMIs

ONT sequencing introduces several complications for UMI-based deduplication:

1. **Elevated error rates.** ONT per-base error rates of 5--15% (depending on chemistry and basecaller) mean that identical UMIs may be read differently across copies of the same molecule.

2. **Homopolymer-length errors.** The dominant ONT error mode is miscounting bases in homopolymer runs (e.g., `AAAA` read as `AAA` or `AAAAA`). Standard Levenshtein distance penalizes these errors equally with substitutions, inflating distances between copies of the same molecule.

3. **Reduced alphabet.** ONT UMI designs typically use a 3-letter alphabet (A, C, G --- no T) to avoid homopolymer runs at synthesis boundaries. The smaller sequence space ($3^L$ vs. $4^L$) increases the probability of random collisions between independent UMIs.

4. **Long reads spanning splice junctions.** ONT reads in RNA applications can span entire transcripts, creating the possibility that reads from different transcript isoforms overlap positionally but represent distinct molecules.

### Overview

We describe a three-phase pipeline: (1) spatial grouping of reads by genomic position, (2) pairwise edit distance computation with optional homopolymer awareness, and (3) graph-based clustering with multiple available algorithms. We also present an adaptive statistical threshold that adjusts the effective edit distance per component based on the expected rate of random collisions, and splice junction-aware grouping for RNA applications.

## Methods

### Read grouping by genomic position

Reads from a coordinate-sorted BAM file are grouped into components using a union-find data structure. Two reads are candidates for the same component if they satisfy positional overlap criteria on the same strand.

#### Overlap detection

We define two modes of positional overlap, controlled by the `--overlap` parameter $g$ (default 50 bp):

**Default mode (match-both-ends).** Two reads are grouped if their 3' ends differ by at most $g$ bp. The 5' proximity constraint is enforced implicitly: because the BAM is coordinate-sorted, the streaming algorithm ejects reads whose start position falls more than $g$ behind the current position. Any two reads still in the active buffer therefore have 5' starts within $g$ of each other.

**Match-one-end mode.** Two reads are grouped if *either* their 5' starts or their 3' ends differ by at most $g$ bp. This mode accommodates 5' degradation, where reads from the same molecule may not share a common 5' boundary.

#### Bin-indexed detection

Naive pairwise comparison of buffered reads would require $O(B)$ comparisons per incoming read, where $B$ is the buffer size. Instead, reads are indexed in two hash maps keyed by $\lfloor \text{position} / g \rfloor$:

- `endIndex`: bins reads by 3' end position
- `startIndex`: bins reads by 5' start position

Each incoming read queries at most 3 bins (target $\pm 1$) per index. Combined with a root-skip optimization (reads whose union-find root matches the current read's root are skipped after a single `find()` call), detection is $O(\text{matches})$ per read rather than $O(B)$.

#### Splice junction matching

For RNA applications, positional overlap alone is insufficient to distinguish reads from different transcript isoforms. When enabled, splice junction matching adds an additional requirement: two reads must have compatible sets of splice junctions, extracted from CIGAR `N` operations.

**Junction extraction.** Each `N` operation in the CIGAR string produces a splice junction defined by its donor position (where the intron begins in reference coordinates) and acceptor position (where the next exon begins).

**Pre-merging.** Adjacent junctions separated by a gap of at most $w$ bp (default 20) are merged into a single spanning junction. This handles cases where a small exon is present in one read's alignment but missed in another --- a common occurrence with ONT alignments due to elevated indel rates near exon boundaries.

**Compatibility rules.** Let $J_a$ and $J_b$ denote the (pre-merged) junction sets of two reads. Two junctions $j_a$ and $j_b$ are considered matching if both $|j_a.\text{donor} - j_b.\text{donor}| \leq w$ and $|j_a.\text{acceptor} - j_b.\text{acceptor}| \leq w$.

- If both reads have no junctions: compatible.
- If one read has junctions and the other does not: incompatible.
- **Default mode:** $|J_a| = |J_b|$ and each junction matches pairwise (exact set match within tolerance).
- **Match-one-end mode:** one read's junction set may be a contiguous sub-sequence of the other's, anchored at the matching end. When 3' ends match, the shorter set must be a suffix of the longer set; when 5' starts match, the shorter set must be a prefix. This accommodates 5' truncation where the shorter molecule is missing junctions from the truncated end.

### Edit distance computation

#### Standard Levenshtein distance

UMI similarity is measured by Levenshtein edit distance on the normalized (slash-separated) UMI string. The standard recurrence is:

$$d[i][j] = \begin{cases} j & \text{if } i = 0 \\ i & \text{if } j = 0 \\ d[i-1][j-1] & \text{if } a_i = b_j \\ 1 + \min(d[i-1][j],\ d[i][j-1],\ d[i-1][j-1]) & \text{otherwise} \end{cases}$$

where $a_i$ and $b_j$ are the $i$-th and $j$-th characters of strings $a$ and $b$ respectively.

#### Ukkonen's cutoff

For bounded computation (threshold $T$), we apply Ukkonen's cutoff: after filling each row $i$, if $\min_j d[i][j] > T$, the function returns $T+1$ immediately. This is correct because the row minimum is monotonically non-decreasing:

$$\min_j d[i+1][j] \geq \min_j d[i][j]$$

For a threshold of 3 on 19-character UMIs, the vast majority of dissimilar pairs exit after 3--4 DP rows instead of the full $19 \times 19$ matrix, providing a 5--10$\times$ speedup.

#### Homopolymer-aware edit distance

ONT's dominant error mode is homopolymer-length variation (e.g., `AAAA` read as `AAA`). We introduce an HP-aware edit distance that tolerates one such error per UMI segment while preserving the cost of substitutions and long HP collapses.

**Augmented state.** The DP is extended with a binary state variable $f \in \{0, 1\}$ tracking whether the single free HP indel for the current UMI segment has been consumed:

$$d[i][j][f]$$

where $f = 0$ indicates the free indel is available, and $f = 1$ indicates it has been used.

**HP context.** A character at position $i$ in string $a$ is an HP continuation if $a_i = a_{i-1}$ and $a_i \neq$ `/`. The `/` separator delimits UMI segments (typically 4-mer blocks).

**Transitions.** For each cell $(i, j)$, three operations are considered:

*Match/substitution* (from $d[i-1][j-1][f]$):

$$\text{cost} = \begin{cases} 0 & \text{if } a_i = b_j \\ 1 & \text{otherwise} \end{cases}$$

State $f$ resets to 0 if either $a_i$ or $b_j$ is a separator; otherwise $f$ carries forward.

*Deletion of $a_i$* (from $d[i-1][j][f]$):

$$\text{cost} = \begin{cases} 1,\ f' = 0 & \text{if } a_i = \text{`/'} \quad \text{(separator)} \\ 0,\ f' = 1 & \text{if HP context and } f = 0 \quad \text{(free)} \\ 1,\ f' = 1 & \text{if HP context and } f = 1 \quad \text{(paid)} \\ 1,\ f' = f & \text{otherwise} \end{cases}$$

*Insertion of $b_j$* (from $d[i][j-1][f]$): symmetric to deletion, checking HP context in $b$.

**Segment-level free indel.** The free HP indel is shared between insertions and deletions within a segment --- only one total HP indel per segment is discounted. The state $f$ resets when either string crosses a `/` boundary, giving each UMI segment its own independent free indel.

**Final answer:** $\min(d[m][n][0],\ d[m][n][1])$.

**Properties.** Representative examples with the HP-aware metric:

| $a$ | $b$ | Standard | HP-aware | Explanation |
|-----|-----|----------|----------|-------------|
| AAGA | AAAA | 1 | 1 | Substitution preserved |
| AAACG | AACG | 1 | 0 | $\pm 1$ HP error absorbed |
| AAAA | AA | 2 | 1 | One free + one paid |
| AAAA | A | 3 | 2 | One free + two paid |
| GGGG | G | 3 | 2 | Long HP collapse penalized |
| CCGA | CGAA | 2 | 1 | Shared free across both strings |
| CCGA/CCCC | CGAA/CGGC | 4 | 3 | Segment 1: 1 (free HP del + paid HP ins); segment 2: 2 (two subs) |

This design avoids two failure modes of simpler approaches. Full HP compression (reducing all runs to length 1) distorts substitutions: `AAGA` $\to$ `AGA` and `AAAA` $\to$ `A`, inflating distance from 1 to 2. Making all HP indels free (cost 0) is too lenient: `AAAA/GGGG/CCGA/CCCC` and `CCCC/AAAA/CCGA/GGGG` (standard distance 12) appear at HP-distance 3, causing false clustering.

### All-pairs edge finding

For each component, all pairs of unique UMIs within the edit distance threshold $T$ are identified. The computation is parallelized by distributing rows round-robin across worker threads, each with its own DP buffer to avoid allocation contention.

The result is a set of edges $E = \{(i, j, d) : d(u_i, u_j) \leq T\}$, where $u_i$ and $u_j$ are normalized UMI strings.

### Adaptive threshold

#### Motivation

The edit distance threshold $T$ determines which UMI pairs are considered potential duplicates. A fixed threshold that works for small components may produce excessive false edges in large components, where the number of random collisions grows quadratically with component size.

#### Collision probability model

For $N$ independent random UMIs of length $L$ over a $k$-letter alphabet, the probability that two UMIs are at exactly Levenshtein distance $d$ is approximated (substitution-only) as:

$$P(\text{exactly } d) = \binom{L}{d} \cdot \frac{(k-1)^d}{k^L}$$

For ONT UMIs with a 3-letter alphabet ($k = 3$):

$$P(\text{exactly } d) = \binom{L}{d} \cdot \frac{2^d}{3^L}$$

The expected number of random pairs at distance $d$ among $N$ UMIs is:

$$E_{\text{false}}(d) = \binom{N}{2} \cdot P(\text{exactly } d) = \frac{N(N-1)}{2} \cdot \binom{L}{d} \cdot \frac{2^d}{3^L}$$

#### Per-distance false positive rate

After edge finding, the false positive rate at distance $d$ is:

$$\text{FPR}(d) = \frac{E_{\text{false}}(d)}{|\{e \in E : e.\text{dist} = d\}|}$$

If $\text{FPR}(d) > \alpha$ (default $\alpha = 0.05$), all edges at distance $d$ are discarded. The effective threshold is the highest distance that survives filtering.

#### Collision probabilities for ONT UMIs

For a typical ONT UMI with $L = 16$ bases over a 3-letter alphabet:

| Distance $d$ | $P(\text{exactly } d)$ | Interpretation |
|---|---|---|
| 1 | $7.4 \times 10^{-7}$ | $\sim 1$ in 1.3 million |
| 2 | $1.1 \times 10^{-5}$ | $\sim 1$ in 90,000 |
| 3 | $1.0 \times 10^{-4}$ | $\sim 1$ in 9,600 |

The expected number of false pairs grows rapidly with component size:

| $N$ (unique UMIs) | $d=1$ expected false | $d=2$ expected false | $d=3$ expected false |
|---|---|---|---|
| 100 | $< 0.1$ | 0.1 | 0.5 |
| 500 | 0.1 | 1.4 | 13 |
| 1,000 | 0.4 | 5.6 | 52 |
| 5,000 | 9.3 | 139 | 1,301 |
| 10,000 | 37 | 558 | 5,203 |

At $\alpha = 0.05$, the minimum number of real edges at distance $d$ needed to survive filtering is $E_{\text{false}}(d) / \alpha$:

| $N$ | $d=1$ min edges | $d=2$ min edges | $d=3$ min edges |
|---|---|---|---|
| 100 | 1 | 1 | 10 |
| 1,000 | 7 | 111 | 1,040 |
| 5,000 | 186 | 2,787 | 26,013 |
| 10,000 | 743 | 11,150 | 104,063 |

#### Effect of alphabet size

The 3-letter ONT alphabet has a much smaller sequence space than the standard 4-letter DNA alphabet ($3^{16} \approx 4.3 \times 10^7$ vs. $4^{16} \approx 4.3 \times 10^9$), making random collisions 20--67$\times$ more likely:

| Distance $d$ | 3-letter $P$ | 4-letter $P$ | Ratio |
|---|---|---|---|
| 1 | $7.4 \times 10^{-7}$ | $1.1 \times 10^{-8}$ | 66.5$\times$ |
| 2 | $1.1 \times 10^{-5}$ | $2.5 \times 10^{-7}$ | 44.3$\times$ |
| 3 | $1.0 \times 10^{-4}$ | $3.5 \times 10^{-6}$ | 29.6$\times$ |

With the default $\alpha = 0.05$ and approximately 1 real neighbor per UMI at $d = 3$, the adaptive threshold begins excluding $d = 3$ edges at around $N \approx 1{,}000$ unique UMIs.

### Clustering methods

All clustering methods operate on the same edge set (optionally filtered by the adaptive threshold). They differ in how edges are used to form clusters. In all methods, UMIs are pre-sorted by read count (descending), with ties broken lexicographically for deterministic output.

#### Connected (single-linkage)

Union-find is applied directly to all edges. Two UMIs are in the same cluster if they are connected by any path of edges, regardless of path length. This is equivalent to computing connected components of the UMI similarity graph.

**Limitation:** Single-linkage clustering suffers from *chaining*: if $d(A, B) \leq T$ and $d(B, C) \leq T$, then $A$, $B$, and $C$ form one cluster even if $d(A, C) \gg T$. In large components, this can produce mega-clusters containing thousands of UMIs.

#### Adjacency (default)

A greedy, center-based algorithm. UMIs are processed in count-descending order:

1. Each unassigned UMI becomes a cluster center.
2. All unassigned direct neighbors (within the edge set) join that cluster.
3. No further expansion occurs.

Each UMI is assigned exactly once, and merges are strictly one-hop. This eliminates chaining entirely: every member of a cluster is a direct neighbor of the center within the edit distance threshold.

**Limitation:** UMIs reachable only through intermediates are not clustered, even if the intermediates are close. If $A \to B$ (d=1) and $B \to C$ (d=1) but no edge $A \to C$ exists, then $C$ will not join $A$'s cluster.

#### Directional

Edges are pre-filtered by a PCR error count model before union-find. An edge between UMI $i$ (low count) and UMI $j$ (high count) survives only if:

$$\text{count}(i) \leq 2 \cdot \text{count}(j) \cdot \left(\frac{1}{4}\right)^{d(i,j)}$$

This models the expectation that PCR errors of a high-count UMI should appear at low frequency proportional to the per-base error rate raised to the power of the edit distance. After filtering, standard union-find is applied to the surviving edges.

This method is from UMI-tools (Smith et al., *Genome Research* 2017) and is calibrated for Illumina error rates ($\sim 0.1\%$). The assumed per-base error rate of $\frac{1}{4}$ is generally too conservative for ONT data where error rates are 5--15%.

**Properties:** Two equally-expressed UMIs do not merge (the count ratio test fails in both directions). Chaining can still occur on surviving edges.

#### Tiered (distance-attenuated BFS)

A novel method combining center selection from adjacency with multi-hop expansion at decreasing stringency. Starting from the highest-count unassigned UMI (the center), a BFS expands outward with the allowed edit distance decreasing by 1 at each hop:

$$\text{maxDist}(\text{hop}) = T - \text{hop}$$

where $T$ is the edit distance threshold. Expansion stops when $\text{maxDist} < 1$.

For $T = 3$:

| Hop | Max distance | Effect |
|-----|-------------|--------|
| 0 (center) | 3 | Neighbors up to $d = 3$ join |
| 1 | 2 | Their neighbors up to $d = 2$ join |
| 2 | 1 | Their neighbors up to $d = 1$ join |
| 3 | 0 | Stop |

**Properties:** Close UMIs ($d = 1$) can chain through multiple hops, capturing errors-of-errors. Distant UMIs ($d = 3$) are limited to direct neighbors of the center --- a chain $A \xrightarrow{d=3} B \xrightarrow{d=3} C$ is blocked because hop 2 requires $d \leq 2$.

Each complete BFS runs to completion before the next center is selected. Clusters are never merged after formation.

### Representative selection

The representative UMI for each cluster is the member with the highest read count. Ties are broken first by normalized string length (descending), then by lexicographic order (ascending). This ensures deterministic output across runs. The representative is always written in normalized (slash-separated) form.

### Molecule identifier assignment

Each UMI cluster receives a two-level molecule identifier of the form `mi_COMP.CLUST`, where `COMP` is a sequential component (read-overlap group) number and `CLUST` is the cluster index within that component. This encoding preserves information about which clusters were evaluated together during UMI clustering, enabling downstream analysis of clustering decisions.

## Discussion

### Method selection

The four clustering methods represent a spectrum of permissiveness:

| Method | Chaining | Multi-hop | Count-aware | Recommended for |
|--------|----------|-----------|-------------|-----------------|
| Connected | Unbounded | Yes | No | Baseline comparison |
| Adjacency | None | No | No | ONT data (default) |
| Directional | Limited | Yes | Yes | Illumina data |
| Tiered | Distance-bounded | Yes | No | ONT with errors-of-errors |

**Adjacency** is recommended as the default for ONT data because it avoids chaining while still merging UMIs that are direct neighbors within the threshold. Two high-count UMIs that are close (e.g., 2 edits apart) always merge, which handles the ONT error model where two well-amplified copies of the same molecule can have independent sequencing errors.

**Tiered** is appropriate when errors-of-errors are expected --- for example, when a UMI with 2 errors generates a PCR copy with 1 additional error, producing a sequence 3 edits from the original but reachable through an intermediate at 2 edits. The distance attenuation prevents false chaining at high distances while allowing this multi-hop recovery.

### HP-aware distance

The per-segment first-free model strikes a balance between three approaches:

1. **Standard Levenshtein**: penalizes HP errors at full cost, inflating distances between copies of the same molecule.
2. **Full HP compression**: distorts non-HP mutations (e.g., `AAGA` $\to$ `AGA` vs. `AAAA` $\to$ `A` gives distance 2 instead of 1).
3. **All HP indels free**: too lenient, allowing sequences with standard distance 12 to appear at HP-distance 3.

The first-free model tolerates the most common ONT error ($\pm 1$ HP length) while correctly penalizing substitutions and large HP collapses. The segment-level scope ensures that each 4-mer block of the UMI gets its own tolerance budget, preventing accumulation of free indels across the entire UMI.

### Adaptive threshold

The adaptive threshold addresses a fundamental tension in fixed-threshold clustering: a threshold that works for small components (where random collisions are negligible) produces excessive false edges in large components (where the quadratic growth of pairwise comparisons overwhelms the signal).

The per-distance FPR approach is conservative: at $\alpha = 0.05$, up to 1 in 20 edges at a given distance may be a random collision. This is permissive enough to avoid discarding real edges in small components, while aggressive enough to catch the random-collision problem in large ones. The 3-letter ONT alphabet makes this correction particularly important, as collision probabilities are 20--67$\times$ higher than for the standard 4-letter DNA alphabet.

### Splice junction matching

For RNA applications, positional overlap alone is insufficient for correct deduplication. Reads from different transcript isoforms can overlap extensively in genomic coordinates while representing distinct molecules. The junction matching requirement ensures that only reads with compatible splicing patterns are grouped, preventing cross-isoform UMI merging.

The pre-merging of adjacent junctions handles a practical issue with ONT alignments: small exons may be missed by the aligner in some reads but correctly identified in others. Without pre-merging, these reads would appear to have incompatible junction sets despite originating from the same molecule.

## Conclusions

Accurate UMI deduplication in ONT data requires error-tolerant methods that account for the technology's distinctive error profile. The combination of HP-aware edit distance, adaptive statistical thresholding, and multiple clustering strategies provides a flexible framework that can be tuned to the specific characteristics of each dataset. The addition of splice junction matching extends the pipeline to RNA applications where transcript isoform diversity would otherwise confound deduplication.

## References

- Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Research* 27(3):491--499, 2017.
