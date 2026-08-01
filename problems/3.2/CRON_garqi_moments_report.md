# Dyadic GARQI factorial-moment calibration

Date: 2026-08-01

## Result

For every dyadic (X=2^{13},ldots,2^{19}), let

\[
I_X=(X,2X],\qquad P_X=(X/2,2X],\qquad
H_X(n)=\#\{p\in P_X:n-p\in Z_p\}.
\]

The first intensity is stable:

\[
0.6717\leq \lambda_X\log X\leq 0.7238,
\qquad \lambda_X=S_X/X.
\]

The resolvable Poisson ratios do not drift upward.  For (k=2), all seven
ratios lie in ([0.927,1.160]).  For (k=3), they lie in
([0.412,2.648]) with no monotone trend; this larger variation is expected
because each interval contains only 1--11 values with (H=3).  Every scanned
interval has maximum (H=3), so (M_4,M_5,M_6) vanish identically and these
orders are not statistically resolved by this range.

### Intensity and maxima

| (X) | (S_X) | (lambda_X) | (lambda_X\log X) | (max H) | first argmax (n) | number of argmaxes |
|---:|---:|---:|---:|---:|---:|---:|
| 8,192 | 640 | 0.078125000 | 0.703977605 | 3 | 11,576 | 1 |
| 16,384 | 1,222 | 0.074584961 | 0.723776975 | 3 | 18,444 | 3 |
| 32,768 | 2,226 | 0.067932129 | 0.706304454 | 3 | 47,066 | 4 |
| 65,536 | 3,969 | 0.060562134 | 0.671655557 | 3 | 105,433 | 1 |
| 131,072 | 7,802 | 0.059524536 | 0.701407495 | 3 | 182,715 | 3 |
| 262,144 | 14,487 | 0.055263519 | 0.689503546 | 3 | 305,455 | 11 |
| 524,288 | 27,489 | 0.052431107 | 0.690507000 | 3 | 553,239 | 11 |

### Factorial moments

Each entry is

\[
M_k(X)\ \bigl(M_k(X)/(X\lambda_X^k)\bigr).
\]

| (X) | (k=1) | (k=2) | (k=3) | (k=4) | (k=5) | (k=6) |
|---:|---:|---:|---:|---:|---:|---:|
| 8,192 | 640 (1.000) | 58 (1.160) | 6 (1.536) | 0 (0) | 0 (0) | 0 (0) |
| 16,384 | 1,222 (1.000) | 100 (1.097) | 18 (2.648) | 0 (0) | 0 (0) | 0 (0) |
| 32,768 | 2,226 (1.000) | 154 (1.018) | 24 (2.336) | 0 (0) | 0 (0) | 0 (0) |
| 65,536 | 3,969 (1.000) | 238 (0.990) | 6 (0.412) | 0 (0) | 0 (0) | 0 (0) |
| 131,072 | 7,802 (1.000) | 462 (0.995) | 18 (0.651) | 0 (0) | 0 (0) | 0 (0) |
| 262,144 | 14,487 (1.000) | 742 (0.927) | 66 (1.492) | 0 (0) | 0 (0) | 0 (0) |
| 524,288 | 27,489 (1.000) | 1,458 (1.012) | 66 (0.873) | 0 (0) | 0 (0) | 0 (0) |

### Per-(X) verdicts

- (X=8{,}192): (k=2,3) are order one; no drift signal.
- (X=16{,}384): (k=3) rises to 2.648, but this is three (H=3)
  observations rather than a growing upper tail.
- (X=32{,}768): (k=2) moves closer to one and (k=3) falls; no monotone
  drift.
- (X=65{,}536): (k=2) remains at one while the single (H=3) observation
  makes (k=3) fluctuate downward.
- (X=131{,}072): both resolvable ratios remain order one and flat in scale.
- (X=262{,}144): the (k=3) ratio rebounds, but (k=2) remains near one;
  still no monotone drift.
- (X=524{,}288): the largest-(X) ratios are 1.012 and 0.873, returning close
  to the Poisson benchmark.

The data therefore support the AQI calibration at the fixed orders visible in
this sample.  They do not test the requested theoretical regime
(k\asymp\log X), because no value with (H\geq4) occurs.

## Implementation

The program is `CRON_garqi_moments.c`.  It scans every prime
(7\leq p\leq2^{20}).  Rather than divide at every recurrence step, the hot
kernel sets (c_n=(n!)^3b_n).  For (n<p), (n!) is a unit modulo (p), so
(b_n=0\pmod p) if and only if (c_n=0\pmod p), and

\[
c_{n+1}=(34n^3+51n^2+27n+5)c_n-n^6c_{n-1}\pmod p.
\]

The cubic coefficient and (n^6) are advanced by finite differences, and the
two products per step use 32-bit Montgomery reduction.  Eight pthread workers
take the largest remaining prime, retain zero pairs in private vectors, and
merge them after the scan.  The merged pairs are then scattered into the nine
hit arrays for the two reference sizes and seven dyadic sizes.  Progress and
ETA are printed to standard error about every five seconds.

The independent reference path uses the original (b_n) recurrence and an
explicit table of inverses modulo (p); it does not reuse the transformed hot
recurrence.

## Sanity gates and cross-checks

All gates passed.

1. The optimized zero-set kernel agreed position by position with the original
   inverse-table recurrence for (p=7,11,17,31,181,379,3137).
2. The mandatory Python-reference values were reproduced exactly:
   - (X=4000): (S_X=311), (max H=2), (M_2=18);
   - (X=8000): (S_X=622), (max H=3), (M_2=54), (M_3=6).
3. At (X=8192), 20 deterministic pseudorandom (n)'s (SplitMix64 seed
   `0x47524151492d3230`) agreed with direct tests of
   (5b_{n-p}\pmod p) for every (p\in(n/2,n]).  The sample happened to have
   (H=0) at all 20 points, so the scanner additionally checked the positive
   control (n=11576), obtaining (H=3) from both routes.
4. The threaded scan completed its exact expected count of
   41,162,092,072 recurrence steps.  The merged pair list was checked to be
   sorted, in range, and duplicate-free.
5. Below (10^6), this implementation found 78,462 zero pairs and exactly two
   midpoint-zero primes:

   ```text
   11, 3137
   ```

   Both the zero-pair total and the midpoint list agree exactly with the
   independently written sibling scan `CRON_zp_bigscan.c`.

## Performance and reproduction

The final run used an arm64 Mac mini, Apple clang 17.0.0, and eight threads.

| quantity | value |
|---|---:|
| primes scanned | 82,022 |
| recurrence steps | 41,162,092,072 |
| zero pairs through (p\leq2^{20}) | 81,950 |
| scan time | 159.948 s |
| total time including all gates | 160.175 s |
| scan throughput | (2.57346\times10^8) steps/s |

This is well below the 15-minute wall target.  The source SHA-256 for this run
was

```text
771cca004d358081ed447cd61db0e833c1817fdc1a1b86ae8b393cf17c3066d5
```

From the repository root:

```sh
cc -O3 -march=native -Wall -Wextra -Wpedantic \
  -o /tmp/CRON_garqi_moments problems/3.2/CRON_garqi_moments.c \
  -lpthread -lm
/tmp/CRON_garqi_moments 8
```
