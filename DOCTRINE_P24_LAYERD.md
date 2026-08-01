# DOCTRINE — P2.4 Q⁻ Layer D/E closure

**Goal (one sentence):** drive `lean/RamanujanChallenge/Problem24QuadraticAlt.lean`
from 11 sorries down to 1 (only the frontier theorem
`alternatingQuadraticEulerTerm24_hasSum` may remain open).

## Established ground truth (verified, do not re-litigate)

- The whole certificate chain is numerically correct end to end, 25 digits:
  series = coefficient integral = −2·V·Dminus (IBP) = W0/H1/H2 form (Möbius)
  = six-integral combination = `alternatingQuadraticEulerValue24` = −0.062366096543714235.
- **All endpoint singularities on [0,1] are REMOVABLE** (numerically confirmed):
  - `F = −2·V(x)·J(−x)` → 0 at x→0⁺ (~x²log²x/2) and at x→1⁻ (V~(1−x)²/4 beats J~log²(1−x))
  - `A = −2(log x/(x(1+x)))·J(−x)` → 0 at x→0⁺ (~x·log x)
  - `B = (−2V)(−Dminus)` → 0 at x→1⁻ (V(1)=0 kills the log blow-up)
  So no "integrable singularity" machinery is needed — continuous extension suffices.
- `W0(1) = 0` exactly; the six I_ab all converge.
- Already banked this run: `logSq_mul_self_tendsto`, `log_mul_self_tendsto`,
  `slope_tendsto_of_hasDerivAt_zero`, `quadAltV_mul_self_tendsto`,
  `quadAltJneg_hasDerivAt_zero`, `quadAltJneg_zero`, `quadAltF_tendsto_zero_right`.

## Avenues, ranked

### (a) Right endpoint limit + FTC swap  — kills the `hFcont` stub
Prove `Tendsto F (𝓝[<] 1) (𝓝 0)`, then replace
`intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le` (needs ContinuousOn on the
closed interval) with `integral_eq_sub_of_hasDerivAt_of_tendsto` (needs the two
endpoint limits). Technique mirrors the left endpoint: reflect the atoms through
x ↦ 1−x, then split `F = (V/(1−x))·((1−x)·J(−x))`, both factors → 0.
**Terminal:** `hFcont` deleted and the file compiles, OR a written verdict that
the tendsto-FTC cannot be threaded (with the exact type mismatch).

### (b) Integrability of A, B, A+B — kills `hFint`, `hA_int`, `hB_int`
Each integrand extends continuously to [0,1] and Lean's junk conventions already
give the right endpoint VALUES (`0/0 = 0`, `V 1 = 0`). So: prove
`ContinuousOn` on `Icc 0 1` and apply `ContinuousOn.intervalIntegrable`.
Fallback if continuity at the endpoint resists: bound |A|,|B| by an explicit
constant near the bad end and use a bounded-measurable → integrable route.
**Terminal:** three stubs closed, OR a concrete failing goal state written down.

### (c) `quadAltSixIntegralLinear` (2 sorries)
Linearity/splitting of the six-integral combination. Expect these to be
integrability side conditions of the same family as (b); reuse (b)'s lemmas.
**Terminal:** closed, or the residual is shown to need a genuinely new bound.

### (d) `quadAltI11_eq_integral` (4 sorries) — Layer E
First endpoint evaluation. Only after (a)-(c).

## Fallbacks

- If the continuous-extension route for (b) fails, use
  `MeasureTheory.IntegrableOn` of a bounded measurable function on a finite
  measure set, bounding via the same asymptotic atoms.
- If a limit resists, the repo's `Problem26Cyclotomic` proves the same
  `log²x·log(1±x)` family (private — reuse the TECHNIQUE, not the lemma).

## Anti-zombie note

Before proving anything new: grep the repo. This run has already found TWO stubs
that were false gaps (`hcontM` derivable from a lemma proved 15 min earlier;
`quadAlt_dilog_landen` already in `Problem26WeightThree`). Assume more exist.

## Layer E targets — the six endpoint integrals, PSLQ-verified at 60 digits

Computed independently (mpmath, dps=60, PSLQ against the weight-4 basis
`{Li4(1/2), log⁴2, log²2·ζ2, log2·ζ3, ζ2²}`). These are the exact targets to
formalize; do not re-derive them, and do not accept a closed form that
disagrees with one of these.

```
I10 = -2 Li4(1/2) -  (1/12) log⁴2 -      log²2·ζ2 -  (7/4) log2·ζ3 + (1/10) ζ2²
I11 =                                                -(7/2) log2·ζ3 +  (3/4) ζ2²
I12 = -6 Li4(1/2) -   (1/4) log⁴2 +  3   log²2·ζ2 - (21/4) log2·ζ3 +  (9/5) ζ2²
I20 = -2 Li4(1/2) -  (1/12) log⁴2 +  (1/2) log²2·ζ2 - (7/4) log2·ζ3 +  (1/4) ζ2²
I21 = -6 Li4(1/2) -   (1/4) log⁴2 -  (3/2) log²2·ζ2 - (7/2) log2·ζ3 + (51/20) ζ2²
I22 = -6 Li4(1/2) -   (1/4) log⁴2 +  (3/2) log²2·ζ2 - (21/4) log2·ζ3 + (23/10) ζ2²
```

**`I11 = cubicLinearEulerValue24` exactly** — the repo's `Tplus` is
`-(7/2) log2 ζ3 + (3/4) ζ2²`, identical to the I11 row. Verified to 48 digits.
So the coefficient of `Tminus` in I11 is ZERO; an earlier reading of the
certificate that mixed `Tplus` and `Tminus` into this individual row was wrong —
that mixture belongs only to the final six-integral combination.

Consistency of the whole table, checked symbolically term by term on
`-2I10 - 2I11 + 2I12 + 4I20 + 6I21 - 5I22`:

| basis element | from the table | target |
| --- | --- | --- |
| Li4(1/2)   | 4 - 12 - 8 - 36 + 30 = -22 | -22 |
| log⁴2      | (2-6-4-18+15)/12 = -11/12  | -11/12 |
| log²2·ζ2   | 2+6+2-9-7.5 = -13/2        | -13/2 |
| log2·ζ3    | 3.5+7-10.5-7-21+26.25 = -7/4 | -7/4 |
| ζ2²        | -0.2-1.5+3.6+1+15.3-11.5 = 67/10 | 67/10 |

and numerically the combination agrees with `alternatingQuadraticEulerValue24`
to 57 digits. The certificate's endpoint table is therefore independently
confirmed, not merely transcribed.

Numerical note: the `/t` and `/(1-t)` rows need the substitution `t = e^{-u}`
(resp. `1-t = e^{-u}`) to reach PSLQ precision — direct quadrature on `[0,1]`
stalls around 14 digits near the log² endpoint and makes PSLQ report "no
relation" for I10, which is an artifact, not a fact about the integral.

## Layer E route: build on Problem24Euler's existing kernel library

`Problem24Euler.lean` already contains a dozen-plus weight-3/4 integral
evaluations built on a consistent template:

```
<name>Kernel24        the integrand
<name>Primitive24     its antiderivative
<name>Primitive24_tendsto_zero / _tendsto_one    the two ENDPOINT LIMITS
<name>Kernel24_intervalIntegrable
<name>Term24_hasSum_integral : HasSum term (∫ kernel)
<name>Integral24      : ∫ kernel = closed form
```

Kernels present include `minusLogSquareKernel24 = log(1-x)²/x`,
`plusLogSquareKernel24 = log(1+x)²/x`, `minusRadialKernel24 = log x log(1-x)/x`,
plus paired/cross/radial variants. **The repo already uses the endpoint-limit
form of the FTC** (`Primitive24_tendsto_zero/_one`) — the same move rediscovered
in Layer D. Follow this template rather than inventing a new one.

All of them are `private`. Reuse therefore needs de-privatization, for which
there is precedent in this file (five declarations in `Problem24Euler` were made
public so `Problem24QuadraticAlt` could use them).

Relation to Layer E: those kernels are weight-3; substituting `u = 1-t` turns
`I11` into `∫₀¹ [ζ2 - 2 Li₂((1-u)/2) - log²((1-u)/2)] · (-log u)/u du`, which is
weight-4 with a half-argument dilog. So they are building blocks, not matches —
but Layer E is extending an existing library, not starting from scratch.

`cubicLinearEulerTerm24_hasSum : HasSum cubicLinearEulerTerm24 cubicLinearEulerValue24`
is already proved (Problem24Euler.lean:3800) and `I11 = cubicLinearEulerValue24`
exactly, so `I11` reduces to finding an integral representation of the
cubic-linear sum that matches `I11`'s integrand — not to a fresh weight-4
evaluation. Start there.
