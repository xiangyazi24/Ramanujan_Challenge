# Q4279 exact candidate-first q6 leaf hunt

This directory is a computation artifact for the P3.2 q6 fibre question. It is **not** project source and lives only on `chatgpt-drop`.

## Build

The scanner uses only the C++ standard library:

```bash
g++ -O3 -march=native -std=c++20 -pthread -Wall -Wextra -Wpedantic \
  q4279_q6_leaf_hunt.cpp -o q4279_q6_leaf_hunt
```

The constructor syntax is brace-based (`Mont32 m{p}`), so the inherited `FastMod mod(uint32_t(p));` most-vexing-parse defect is absent. No Boost, FLINT, GMP, OpenMP, Python, or nonstandard header is required.

Example:

```bash
./q4279_q6_leaf_hunt \
  --qmin-exclusive 1000000 \
  --qmax 2000000 \
  --threads 4 \
  --block-width 100000 \
  --out results-1000000-2000000
```

## Exact enumeration

For every lower prime in the forced range, the program computes the zero set of the Apéry sequence modulo that prime using the division-free cleared recurrence

```text
B_(n+1) = (34n^3+51n^2+27n+5) B_n - n^6 B_(n-1) (mod p),
B_n = (n!)^3 b_n.
```

Because every index used is below `p`, the factorial is a unit. The program computes through `(p-1)/2` and materializes the other half by the proved reflection symmetry.

Every ordered lower-zero pair is projected by the exact q6 inverse:

```text
plus:  q=(7p-r+a)/6,      t=2r-a-p,
minus: q=p+(r+a+1)/7,    t=(r-6a-6)/7.
```

It then checks integrality, primality, the literal raw plus/minus windows, and the raw minus-first tie-break. The upper mark is tested last by exact membership `t in Z_q`.

The program enumerates the **raw two-row q6 overcarrier**. Archive-only residual/depth/primitive masks are not separately reconstructed in this release. This does not weaken a zero result: an actual filtered leaf is a raw leaf, so absence in the raw overcarrier implies absence after every further deletion. Candidate-size and margin statistics must, however, be read as statistics of the raw candidate graph.

## Built-in validation

Each run performs:

- 12,000 Montgomery multiplication cross-checks;
- full divided-recurrence zero-set comparisons for every prime `7 <= p <= 500`;
- exact equality between the candidate-first inverse projection and an independent direct quotient-row scan through `q=5000`;
- the `605` selected-state checkpoint through `q=5000`;
- reflection and no-consecutive-zero checks throughout the scanned lower range;
- inverse/window identities on every emitted candidate.

A nonzero failure aborts the run.

## Outputs

Each result directory contains:

```text
manifest.json                     run parameters, counts, checks, wall clock
summary.md                        finite result and closest candidate states
witnesses.txt                     witness ledger or explicit finite absence
block_stats.csv                   scale and margin statistics
q_stats.csv.gz                    one row per upper prime
candidate_occurrences.csv.gz      raw/assigned lower-complete occurrences
candidate_states_by_margin.csv.gz candidate states ordered by near-miss
prime_zero_counts.csv.gz          one zero-count row per computed prime
zero_records.csv.gz               full modular zero ledger
stdout.log / stderr-and-time.log  execution provenance
SHA256SUMS                        payload digests
```

`q4279_analyze.cpp` is a second standard-library program that derives aggregate histograms, quantiles, and scaling diagnostics from `q_stats.csv` and `block_stats.csv`.

## Interpretation boundary

Every statement produced by this directory is finite. A line saying no leaf was found in an exhausted interval is not a theorem of eventual or global emptiness. If no actual leaf occurs, the actual first mass in that interval is zero and ratios normalized by it are undefined.
