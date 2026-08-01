# Dual low-rank tests for the Apéry projective clock

## Verdict

The algebraic test is negative in every prescribed window: all eight real
mirror-quotient matrices have the generic full column rank, including the four
near-square pressure matrices.  The observed nullity and the null-model
nullity are both zero, so the surplus is zero.  All 120 control windows are
also full column rank.  In the terminology of Q6523 §2.8, there is no candidate
anomaly, hence no structural or proof-grade syzygy.  This is cap evidence for
the tested mirror-saturated Stepanov windows, not a theorem outside those
windows.

The analytic test is not positive: neither the collision channel `Y` nor the
linear control `Z` has a bounded-rank plateau.  `Z` has the maximum possible
diagonal and Kronecker-Hankel rank in every one of its 3,552 paired
specializations.  `Y` reaches maximum rank at the largest windows for every
tested long twist when `p=3001,10007`, and almost does so for `p=1009`; however,
the translation twists `(a,b)=(1,±1)` have startup defects on the third-largest
window.  Consequently the literal Q6523 strong-negative condition—every long
twist on each of the last three scales—does not hold.  The §3 verdict is
therefore **intermediate growth with strong negative evidence**, not a
proof-grade cap.  Fixed-rank affine clock closure is strongly disfavored, but
the tested data do not justify excluding slower/transient, semilinear,
nonlinear, Möbius, or differently coordinated closures.

## Reproduction and checks

Commands used:

```text
sage -python CRON_lowrank_alg.py
sage -python CRON_lowrank_ana.py --prime 1009
sage -python CRON_lowrank_ana.py --prime 3001
sage -python CRON_lowrank_ana.py --prime 10007
```

Both programs emit JSON Lines, including every raw rank.  Before any production
case, the algebra script checks:

1. coefficient-by-coefficient agreement of `N_h` with `CRON_avggcd.py`;
2. `N_h(r)=0` iff the two computed projective orbit points agree, for every
   valid `(r,h)` in the `p=101, h<=8` test;
3. the required centered mirror parity and theoretical quotient degrees;
4. equality of the direct `G_h^M` matrix rank and the equivalent Hasse-jet
   matrix rank on a small squarefree example;
5. a hand-checkable one-row case.

The analytic script checks determinant zero against projective equality at
every valid `(r,h)`, checks and removes the forced even-`h` mirror event, checks
`Y(k,h)=Y(-k,h)` on a small instance, verifies the finite-field root has exact
order `p`, and directly recomputes sample Fourier sums with Python integers.
All ranks are Sage ranks over `GF(q)`; no complex FFT, SVD, or tolerance enters
the rank decision.  Fourier dot products are evaluated by binary64 matrix
multiplication only after asserting the worst possible integer sum is below
`2^53`; all products and partial sums are therefore exactly representable, and
direct integer spot checks are also enforced.

## 1. Algebraic mirror-quotient syzygy test

The pressure degree is the unique largest integer `A_pressure` satisfying
`Q=(A_pressure+1)(B_eff+1)<T`, as required by the Q6523 near-square pressure
rule.  The actual squarefree quotient degrees, not theoretical degrees, are
used in `T`.

### Raw real-matrix data

| p | H | M | window | A | B | T rows | Q cols | rank | nullity0 | nullity | surplus |
|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1009 | 12 | 1 | A | 21 | 3 | 96 | 88 | 88 | 0 | 0 | 0 |
| 1009 | 12 | 1 | pressure | 22 | 3 | 96 | 92 | 92 | 0 | 0 | 0 |
| 2003 | 20 | 1 | A | 41 | 5 | 280 | 252 | 252 | 0 | 0 | 0 |
| 2003 | 20 | 1 | pressure | 45 | 5 | 280 | 276 | 276 | 0 | 0 | 0 |
| 5003 | 16 | 3 | A | 94 | 4 | 528 | 475 | 475 | 0 | 0 | 0 |
| 5003 | 16 | 3 | pressure | 104 | 4 | 528 | 525 | 525 | 0 | 0 | 0 |
| 9001 | 32 | 2 | A | 146 | 8 | 1472 | 1323 | 1323 | 0 | 0 | 0 |
| 9001 | 32 | 2 | pressure | 162 | 8 | 1472 | 1467 | 1467 | 0 | 0 | 0 |

There were no squarefree degree drops and no accidental odd-`h` center roots
in any real case.

### Raw control data

Each entry below is `full-column-rank trials / trials`, separately for the A
and pressure windows.  Seeds are `1000..1004`, `1100..1104`, and `1200..1204`
for the three controls.

| p | same-degree random A / pressure | mirror-preserving random A / pressure | clock-shuffle A / pressure |
|---:|:---:|:---:|:---:|
| 1009 | 5/5 ; 5/5 | 5/5 ; 5/5 | 5/5 ; 5/5 |
| 2003 | 5/5 ; 5/5 | 5/5 ; 5/5 | 5/5 ; 5/5 |
| 5003 | 5/5 ; 5/5 | 5/5 ; 5/5 | 5/5 ; 5/5 |
| 9001 | 5/5 ; 5/5 | 5/5 ; 5/5 | 5/5 ; 5/5 |

The same-degree and mirror-preserving controls are distinct seeded trials.  In
the mirror quotient their ensembles have the same distribution; the latter is
also lifted to the prescribed even/odd polynomial in `z` and round-tripped by
an assertion before the quotient matrix is built.

### Q6523 §2.8 three-level decision

| level | required signal | observed decision |
|:---|:---|:---|
| Candidate anomaly | surplus at least 1 in `Q<T`, repeated across primes and rare in controls | **No:** surplus is 0 in all eight real windows |
| Structural syzygy | stable support/reconstruction over at least four primes and nested H, stable jets | **No candidate to reconstruct** |
| Proof-grade surplus | low-bidegree primitive generator or a shift module with growing surplus | **No:** every tested map is injective |

## 2. Clock-Fourier Hankel tests

### Finite-field and window data

| p | Hmax | auxiliary q values | dyadic L for h0=1 | h0=3 | h0=7 | primitive events through Hmax |
|---:|---:|:---|:---|:---|:---|---:|
| 1009 | 31 | 10091, 12109 | 2,4,8,16 | 2,4,8 | 2,4,8 | 26 |
| 3001 | 54 | 30011, 36013 | 2,4,8,16 | 2,4,8,16 | 2,4,8,16 | 60 |
| 10007 | 100 | 240169, 380267 | 2,4,8,16,32 | 2,4,8,16,32 | 2,4,8,16,32 | 103 |

The four seeds are `1, floor(p/7), floor(2p/7), floor(3p/7)`.  All 24 affine
twists specified by Q6523 were run.  The table below restricts only the compact
display to the 18 actual long-orbit twists at `s=8`, excluding the two stated
identity/mirror sanity twists and other short-orbit cases.  `full` counts
full-rank twists out of those 18.  Both auxiliary characteristics give every
displayed rank exactly; no specialization instability occurred.

### Raw `Y` growth curves at s=8

`diag` is the range across four seeds and the displayed long twists.
`bispectral` is `min..max / rmax` across the same twists.

| p | h0 | L | diag | bispectral | full |
|---:|---:|---:|:---:|:---:|:---:|
| 1009 | 1 | 2 | 1..1 | 4..8 / 16 | 0/18 |
| 1009 | 1 | 4 | 4..4 | 20..32 / 32 | 16/18 |
| 1009 | 1 | 8 | 8..8 | 46..64 / 64 | 16/18 |
| 1009 | 1 | 16 | 16..16 | 128..128 / 128 | 18/18 |
| 1009 | 3 | 2 | 2..2 | 8..16 / 16 | 16/18 |
| 1009 | 3 | 4 | 4..4 | 12..32 / 32 | 16/18 |
| 1009 | 3 | 8 | 8..8 | 46..64 / 64 | 16/18 |
| 1009 | 7 | 2 | 2..2 | 4..16 / 16 | 16/18 |
| 1009 | 7 | 4 | 2..2 | 4..16 / 32 | 0/18 |
| 1009 | 7 | 8 | 8..8 | 62..64 / 64 | 16/18 |
| 3001 | 1 | 2 | 2..2 | 8..16 / 16 | 16/18 |
| 3001 | 1 | 4 | 4..4 | 26..32 / 32 | 16/18 |
| 3001 | 1 | 8 | 8..8 | 64..64 / 64 | 18/18 |
| 3001 | 1 | 16 | 16..16 | 128..128 / 128 | 18/18 |
| 3001 | 3 | 2 | 2..2 | 6..16 / 16 | 16/18 |
| 3001 | 3 | 4 | 4..4 | 22..32 / 32 | 16/18 |
| 3001 | 3 | 8 | 8..8 | 64..64 / 64 | 18/18 |
| 3001 | 3 | 16 | 16..16 | 128..128 / 128 | 18/18 |
| 3001 | 7 | 2 | 2..2 | 4..16 / 16 | 16/18 |
| 3001 | 7 | 4 | 4..4 | 16..32 / 32 | 16/18 |
| 3001 | 7 | 8 | 8..8 | 60..64 / 64 | 16/18 |
| 3001 | 7 | 16 | 16..16 | 128..128 / 128 | 18/18 |
| 10007 | 1 | 2 | 2..2 | 4..16 / 16 | 16/18 |
| 10007 | 1 | 4 | 4..4 | 30..32 / 32 | 16/18 |
| 10007 | 1 | 8 | 8..8 | 64..64 / 64 | 18/18 |
| 10007 | 1 | 16 | 16..16 | 128..128 / 128 | 18/18 |
| 10007 | 1 | 32 | 32..32 | 256..256 / 256 | 18/18 |
| 10007 | 3 | 2 | 2..2 | 10..16 / 16 | 16/18 |
| 10007 | 3 | 4 | 4..4 | 32..32 / 32 | 18/18 |
| 10007 | 3 | 8 | 8..8 | 56..64 / 64 | 16/18 |
| 10007 | 3 | 16 | 16..16 | 128..128 / 128 | 18/18 |
| 10007 | 3 | 32 | 32..32 | 256..256 / 256 | 18/18 |
| 10007 | 7 | 2 | 2..2 | 8..16 / 16 | 16/18 |
| 10007 | 7 | 4 | 4..4 | 18..32 / 32 | 16/18 |
| 10007 | 7 | 8 | 8..8 | 38..64 / 64 | 16/18 |
| 10007 | 7 | 16 | 16..16 | 128..128 / 128 | 18/18 |
| 10007 | 7 | 32 | 32..32 | 256..256 / 256 | 18/18 |

The two non-full long twists in almost every affected row are precisely the
translations `(a,b)=(1,±1)`.  One raw representative is `(1,1)`:

| p | h0 | rank/rmax as L doubles |
|---:|---:|:---|
| 1009 | 1 | 4/16, 20/32, 46/64, 128/128 |
| 1009 | 3 | 8/16, 12/32, 46/64 |
| 1009 | 7 | 4/16, 4/32, 62/64 |
| 3001 | 1 | 8/16, 26/32, 64/64, 128/128 |
| 3001 | 3 | 6/16, 22/32, 64/64, 128/128 |
| 3001 | 7 | 4/16, 16/32, 60/64, 128/128 |
| 10007 | 1 | 4/16, 30/32, 64/64, 128/128, 256/256 |
| 10007 | 3 | 10/16, 32/32, 56/64, 128/128, 256/256 |
| 10007 | 7 | 8/16, 18/32, 38/64, 128/128, 256/256 |

At the largest admissible L, all long twists are full rank for every
`s=1,2,4,8` when `p=3001,10007`.  For `p=1009`, the only largest-window
defects are:

| h0 | L | s=1 | s=2 | s=4 | s=8 |
|---:|---:|:---:|:---:|:---:|:---:|
| 1 | 16 | 16/16 | 32/32 | 64/64 | 128/128 |
| 3 | 8 | 8/8 | 16/16 | 28..32/32 | 46..64/64 |
| 7 | 8 | 8/8 | 16/16 | 32/32 | 62..64/64 |

### Raw `Z` control result

For every `p`, both `q`, all `h0`, all admissible `L`, all 24 twists, every
effective `s in {1,2,4,8}`, and all four diagonal seeds:

```text
rank_diag = L
rank_bispectral = s_eff * L = rmax
rho = rho2 = 1                 (when rmax > 2)
```

Thus `Z` has 960/960, 1152/1152, and 1440/1440 full paired bispectral records
for `p=1009,3001,10007`, respectively.  `Y` has 718/960, 1094/1152, and
1371/1440 full paired records.  All 7,104 `q` pairs are specialization-stable.
The determinant weights in `Z` use their canonical residues `0,...,p-1`, as
produced by the `% p` arithmetic in the Q6523 Python skeleton, before the
integer Fourier sums are reduced modulo each auxiliary `q`.

### Q6523 §3 three-level decision

| level | required signal | observed decision |
|:---|:---|:---|
| Strong positive | rank at most 8 (or weakly 16) and stable under L-doubling, two q, at least four p, at least two long twists, low-rank Z, exact holdout prediction | **No:** ranks grow to 128 or 256 and `Z` is maximally ranked; the mandated run also has only three main primes |
| Strong negative / cap | linear rank growth for every long twist on all last three scales, stable across p,q,h0,seeds | **Not literally met by Y:** the small-p and translation startup windows fail the all-scales condition, although the final windows are full |
| Intermediate growth | neither a bounded plateau nor the uniform strong-negative threshold | **Yes:** this is the conservative verdict for Y; Z separately gives maximal-rank companion-level negative evidence |

## Scope of the conclusion

The data close the tested low-bidegree mirror-quotient syzygy windows and give
strong evidence against fixed-dimensional affine clock realizations of the
direct collision observable.  They do not prove an absolute no-go for all
hidden structures.  In particular, Q6523's stated caveats remain active:
nonlinear zero-set mechanisms, semilinear/Frobenius actions, different clock
coordinates, Möbius frequency maps, and growth strictly between bounded and
linear rank are not ruled out by these experiments.
