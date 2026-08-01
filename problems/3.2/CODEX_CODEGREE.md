# Transfer-cocycle codegree protocol

This report records the exact output of
`python3 research/scripts/q32_codegree_protocol.py`. The run used Python
3.9.6, NumPy 1.26.4 for the exact integer bucket-Gram calculation, and the
version-independent SplitMix64 generator with seed `20260801` for the
annealed trajectories. The pure-stdlib bit-plane backend was separately
cross-checked against NumPy at p = 101 and 211.

## Mandatory gate

The gate passed at all 204 distinct primes used by the experiments: all five
Experiment-1 primes and the first 200 primes at least 1000, with 1009 in both
sets. At each prime:

- the transfer orbit's infinity visits for `1 ≤ n ≤ p−2` were identical to
  the zeros obtained directly from the Apéry recurrence modulo p;
- `x₁ = 5⁻¹ (mod p)`;
- every visit was followed by the forced restart `x_(n+1) = 0`;
- `b_(p−1−n) = b_n (mod p)`, and hence the visit set was exactly invariant
  under `n ↦ p−1−n`;
- the projective state identity
  `x_(p−1−n) x_(n+1) = −(n+1)⁶`
  was verified for every `0 ≤ n ≤ p−2`. Projective interpretation makes
  this identity valid when one factor is zero and the other is infinity.

The five Experiment-1 gate records were:

| p | x₁ | infinity visits |
|---:|---:|:---|
| 101 | 81 | `[]` |
| 211 | 169 | `[100,110]` |
| 401 | 321 | `[]` |
| 601 | 481 | `[]` |
| 1009 | 202 | `[422,586]` |

The script prints a `GATE VERIFIED` line, including the complete ordered
visit set, for every one of the 204 primes.

## Experiment 1: codegree exceptional locus

All `(p+1)²` pairs were scanned at every requested prime; p = 1009 was not
dropped. Infinity is represented projectively by `[1:0]`, while a finite x
is `[x:1]`. The exact flagged sets were:

| p | projective-boundary flags | finite diagonal flags | finite off-diagonal flags | total |
|---:|---:|:---|---:|---:|
| 101 | 203 (every pair with an infinity coordinate) | all `(x,x)` except x = 29 | 0 | 303 |
| 211 | 423 | all 211 | 0 | 634 |
| 401 | 803 | all 401 | 0 | 1204 |
| 601 | 1203 | all 601 | 0 | 1804 |
| 1009 | 2019 | all 1009 | 0 | 3028 |

Thus the compact descriptions above are the complete flagged sets, not
samples. At p = 101, `C₁₀₁(29,29) = 158`, so its deviation 57 is below
`6√101`; this is the lone threshold miss on the finite diagonal.

The codegree ranges were:

| p | finite diagonal Cₚ, min / mean / max | finite off-diagonal Cₚ, min / max |
|---:|:---|:---|
| 101 | 158 / 199.029703 / 248 | 60 / 136 |
| 211 | 374 / 417.033175 / 464 | 152 / 265 |
| 401 | 720 / 799.007481 / 910 | 315 / 483 |
| 601 | 1104 / 1195.011647 / 1326 | 497 / 711 |
| 1009 | 1884 / 2011.006938 / 2148 | 870 / 1149 |

The finite diagonal therefore has `Cₚ ≈ 2p`. The projective boundary is a
deterministic degeneracy: for finite x,
`Cₚ(x,∞) = Cₚ(∞,x) = 0`, while `Cₚ(∞,∞) = (p−1)²`.

Exact nullspace calculations over every Fₚ found the same two bilinear
relations at all five primes:

1. `XZ′ − ZX′ = 0`, the projective diagonal `x = x′`;
2. `ZZ′ = 0`, the reducible boundary `x = ∞` or `x′ = ∞`.

The second is not an affine correspondence. Among finite flagged pairs, no
relation `xx′ = c` or `x + x′ = c` had support at least four: the maximum
supports were two and one, respectively. There were no finite off-diagonal
flags from which another bilinear fit could arise.

The Apéry reflection does induce a reciprocal relation on states, but its
constant depends on time. Putting `m = p−1−n` and using
`b_(p−1−k) = b_k` gives

`x_m = −(n+1)³ b_(n+1)/b_n`,
and `x_(n+1) = (n+1)³ b_n/b_(n+1)`.

Their product is therefore `−(n+1)⁶`, not a fixed c. The scan finds **no
stable nontrivial affine exceptional relation beyond the diagonal**. The
only additional stable locus is the forced infinity boundary.

## Experiment 2: annealed versus ordered visits

The sample comprises the first 200 primes in the allowed range, from 1009
through 2503. One annealed trajectory of length p−1 was run at each prime.
The complete ordered/annealed record has SHA-256
`6a22c4905fba291e699d95e882c0ecc54a4df037e03b7c2091caddf856284004`.

### Count distributions

| visits m | ordered | annealed | Poisson(1), expected count out of 200 |
|---:|---:|---:|---:|
| 0 | 131 | 74 | 73.575888 |
| 1 | 0 | 70 | 73.575888 |
| 2 | 52 | 40 | 36.787944 |
| 3 | 0 | 14 | 12.262648 |
| 4 | 15 | 2 | 3.065662 |
| 5 | 0 | 0 | 0.613132 |
| 6 | 2 | 0 | 0.102189 |

No central visit at `n = (p−1)/2` occurred. Consequently `Hₚ = Vₚ/2` for
every sampled prime, with distribution:

| paired orbits h | observed | Poisson(1/2), expected count out of 200 |
|---:|---:|---:|
| 0 | 131 | 121.306132 |
| 1 | 52 | 60.653066 |
| 2 | 15 | 15.163266 |
| 3 | 2 | 2.527211 |

### Moments and parity

Here `(m)_k` denotes the falling factorial.

| ensemble | mean | variance | E[(m)₂] | E[(m)₃] |
|:---|---:|---:|---:|---:|
| ordered Vₚ | 0.880000 | 1.825600 | 1.720000 | 3.000000 |
| annealed | 1.000000 | 0.940000 | 0.940000 | 0.660000 |
| ordered pairs Hₚ | 0.440000 | 0.456400 | 0.210000 | 0.060000 |
| Poisson(1) | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| Poisson(1/2) | 0.500000 | 0.500000 | 0.250000 | 0.125000 |

All 200 ordered counts were even. In contrast, 116/200 = 0.580000 of the
annealed counts were even, close to the Poisson(1) value 0.567668. The
annealed total was 200, compared with
`Σₚ (p−1)/(p+1) = 199.753418`.

The annealed count frequencies, mean, variance, and second factorial moment
are consistent with Poisson(1) at this sample size. Its third factorial
moment is low because the realization contains no count above four. The
halved ordered counts are likewise broadly consistent with Poisson(1/2): the
mean is 0.44 and the first four frequency cells are close to their Poisson
expectations. Its third factorial moment is based on only two observations
with `Hₚ = 3`, so it is particularly noisy.

The raw ensembles do **not** agree at the tested resolution. Their means are
close, but the ordered orbit has an exact pairing signature (100% even versus
58% even) and correspondingly much larger variance and factorial moments.
This is a visible quenched-order gap.

### Reflection-pair gaps

Reflection symmetry was checked machine-exactly at all 200 primes. There
were 88 pairs. For the normalized gap `(p−1−2r)/(p−1)`, the distribution was:

| normalized gap | count | fraction |
|:---|---:|---:|
| [0.0, 0.1) | 10 | 0.113636 |
| [0.1, 0.2) | 9 | 0.102273 |
| [0.2, 0.3) | 7 | 0.079545 |
| [0.3, 0.4) | 10 | 0.113636 |
| [0.4, 0.5) | 11 | 0.125000 |
| [0.5, 0.6) | 10 | 0.113636 |
| [0.6, 0.7) | 8 | 0.090909 |
| [0.7, 0.8) | 2 | 0.022727 |
| [0.8, 0.9) | 13 | 0.147727 |
| [0.9, 1.0) | 8 | 0.090909 |

The smallest absolute reflection gap was 4, the largest was 2190, and the
mean normalized gap was 0.487570. No asymmetric visit occurred.

## Experiment 3: post-visit restart structure

Return times are consecutive infinity-visit spacings observed within each
length-p−1 trajectory. The deterministic next state after every visit was
checked to be zero.

| ensemble | completed return intervals | minimum | frequency of minimum | mean return | mean τ/(p+1) |
|:---|---:|---:|---:|---:|---:|
| ordered | 107 | 4 | 3 | 620.130841 | 0.350030 |
| annealed | 74 | 7 | 1 | 510.662162 | 0.308448 |

The short-return counts were:

| cutoff | ordered τ ≤ cutoff | annealed τ ≤ cutoff |
|---:|---:|---:|
| 2 | 0 | 0 |
| 5 | 3 | 0 |
| 10 | 7 | 1 |
| 25 | 9 | 2 |
| 50 | 15 | 3 |
| 100 | 19 | 9 |

The smallest ordered return time is therefore **4, occurring three times**.
The first ordered return-time frequencies were
`4:3, 8:2, 10:2, 11:2, 32:1, 36:1, 43:2, 48:2, 56:1, 64:1, 71:2`.
There was no return at time 2, but the three time-4 returns identify the
short-return danger zone requested by the protocol.

For comparison with the exponential scaling prediction:

| t | ordered empirical Pr(τ/(p+1) ≤ t) | annealed empirical | 1−exp(−t) |
|---:|---:|---:|---:|
| 0.10 | 0.205607 | 0.202703 | 0.095163 |
| 0.25 | 0.476636 | 0.432432 | 0.221199 |
| 0.50 | 0.719626 | 0.837838 | 0.393469 |
| 0.75 | 0.869159 | 0.986486 | 0.527633 |
| 1.00 | 1.000000 | 1.000000 | 0.632121 |

Both empirical CDFs lie well above the unconditioned exponential curve.
This comparison is descriptive: only returns completed before the length-p−1
horizon enter the sample, which preferentially retains short intervals. The
similar distortion in the annealed control confirms that finite-window
selection is substantial; these CDF values should not be read as an
uncensored geometric-law test.

## LIMITATIONS

- The exceptional-locus conclusion is an exact statement about the five
  tested finite fields, not a characteristic-zero proof. The algebraic
  search was restricted to constant product, constant sum, and general
  bilinear relations.
- The threshold `6√p` misses one finite diagonal point at p = 101, showing
  that threshold membership is not identical to algebraic-locus membership
  at the smallest scale.
- Experiment 2 uses a permitted subrange—200 primes from 1009 to 2503—not
  primes spread across the full interval up to 20,000.
- There is only one annealed trajectory per prime. Tail frequencies and
  third factorial moments therefore have appreciable Monte Carlo noise even
  though the seed and generator are fully reproducible.
- Return times are right-censored by the observation horizon, and only
  completed spacings are tabulated. A survival-analysis treatment with many
  independent, longer annealed paths would be needed for a quantitative
  geometric-law goodness-of-fit test.
- All finite-field arithmetic, gates, flag decisions, and relation fits are
  exact. Floating-point values enter only in displayed means and reference
  probabilities.
