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

## I11: the derivation skeleton (independently confirmed, no MZV table needed)

An independent derivation reaches the same value PSLQ found,
`I11 = -(7/2) log2 ζ3 + π⁴/48 = -(7/2) log2 ζ3 + (3/4) ζ2²`, and — importantly
for formalization — **without** invoking `Li4(1/2)`, an alternating Euler-sum
table, or any multiple-zeta reduction. The route:

1. **One integration by parts**, using
   `d/dt Li₂(t/2) = -log(1-t/2)/t`, hence
   `W0'(t) = (2/t)·log((2-t)/t)`, together with `(H(t)²/2)' = H(t)/(1-t)`.
   This is exactly the `hprod'` / `g11` structure already written in
   `quadAltI11_eq_integral` — the Lean file's existing IBP is the right first
   step, not a detour.
   Reduces the problem to three logarithmic integrals.
2. Two of the three are immediate from nonnegative power-series expansions.
3. The third: introduce one parameter and apply **Tonelli** (nonnegativity makes
   the interchange free — this is the step that would otherwise need a
   dominated-convergence argument).
4. Two elementary integrations by parts reduce it to the alternating polylog
   series at `-1`.

Only standard input needed: `∑ H_n/n³ = (5/4) ζ4 = π⁴/72`. Note `ζ4 = π⁴/90 =
(2/5) ζ2²` is ordinary even-zeta normalisation, not an extra MZV input.

Implication for the Lean work: the existing `quadAltI11_eq_integral` skeleton
(IBP with `g11 = H1²/2`, endpoint limits already proved) is on the right track;
what remains after Codex closes its integrability side conditions is steps 2-4,
and the only genuinely new analytic ingredient is one Tonelli interchange on a
nonnegative double series.

## Layer E collapses to ONE weight-4 integral

The six endpoint evaluations are **not** six independent weight-4 problems.
Modulo the lower-product space `span_Q{L²ζ2, Lζ3, ζ4}` their span is
one-dimensional, generated by

```
K := ∫₀¹ log²x · log(1+x)/(1+x) dx
   = 4 Li₄(1/2) + (1/6) L⁴ - L²ζ2 + (7/2) Lζ3 - (15/4) ζ4
```

(closed form verified numerically to 51 digits). In the basis
`{K, L²ζ2, Lζ3, ζ4}` the table is:

```
I10 = -(1/2)K - (3/2)L²ζ2              - (13/8)ζ4
I11 =                      -(7/2)Lζ3   + (15/8)ζ4
I12 = -(3/2)K + (3/2)L²ζ2              -  (9/8)ζ4
I20 = -(1/2)K                          -  (5/4)ζ4
I21 = -(3/2)K -   3  L²ζ2  + (7/4)Lζ3  +  (3/4)ζ4
I22 = -(3/2)K                          +  (1/8)ζ4
```

`L⁴` cancels in every row. That is the structural reason for the grouping noted
earlier — `I10, I20` share `(-2, -1/12)` and `I12, I21, I22` share `(-6, -1/4)`
in the Li₄(1/2) basis precisely because those ratios are the ones that kill `L⁴`
when rewritten through `K`. The rows derived symbolically this way agree
coefficient-by-coefficient with the independent PSLQ output.

The combination is

```
-2I10 - 2I11 + 2I12 + 4I20 + 6I21 - 5I22
  = -(11/2)K - 12 L²ζ2 + (35/2) Lζ3 - (31/8) ζ4
```

which equals `alternatingQuadraticEulerValue24` exactly.

**Consequence for the Lean work.** Layer E needs exactly one hard analytic
result — the evaluation of `K` — plus rational arithmetic and the standard
constants. Two possible shapes:

- prove `K`'s closed form, then the six rows follow by elementary algebra; or
- skip `quadAltSixIntegralLinear`'s split entirely and prove the single combined
  identity, which needs only `K` with coefficient `-11/2`.

The second is likely cheaper and should be tried first. Either way the estimate
"six independent weight-4 evaluations" is wrong — it is one.

## …and K itself reduces to constants the repo has already proved

One step further. Expanding `log(1+x)/(1+x) = Σ_{n≥1} (-1)^{n+1} H_n xⁿ` and using
`∫₀¹ xⁿ log²x dx = 2/(n+1)³`, then shifting the index with
`H_n = H_{n+1} - 1/(n+1)`:

```
K = 2 ( η(4) - A ),    A := Σ_{n≥1} (-1)^{n-1} H_n / n³,   η(4) = (7/8) ζ4
```

and the alternating linear Euler sum `A` satisfies

```
A = Tminus + Tplus + (1/4) ζ2²
```

where `Tplus = cubicLinearEulerValue24` and
`Tminus = alternatingCubicLinearEulerValue24` — **both already proved in this
repo** (`cubicLinearEulerTerm24_hasSum`, and its alternating companion). Hence

```
K = 2 ( (7/8) ζ4 - Tminus - Tplus - (1/4) ζ2² )
```

All three relations verified numerically to 40 digits.

**So Layer E needs no new weight-4 evaluation at all.** The chain is

```
six-integral combination
  = -(11/2) K - 12 L²ζ2 + (35/2) Lζ3 - (31/8) ζ4        [rational algebra]
  = elementary function of Tplus, Tminus, ζ4, ζ2², Lζ3   [the K identity above]
  = alternatingQuadraticEulerValue24                     [already proved:
                                                          bridgeValue_eq_…]
```

The only genuinely new Lean content left is the **integral-to-series step for K**:
expand, integrate term by term, and reindex. That is one interchange on an
alternating series — dominated convergence with `Σ H_n/(n+1)³` as dominator, or
split into even/odd and use Tonelli twice.

Everything downstream of that is rational arithmetic over constants the repo
already owns. The earlier readings — "six independent weight-4 evaluations",
then "one weight-4 evaluation K" — were both overestimates; the true remaining
analytic content is a single termwise integration.

## Layer E build order (decided)

Keep the six named `I_ab` rows; do **not** collapse into one combined theorem as
the primary route. An independent architecture review makes the case: the
combined integrand is genuinely smoother at `t = 0`, but at `t = 1` there is no
extra cancellation among the six coefficients — the vanishing still comes from
`W0`'s double zero — and an explicit primitive of the combined bracket is just
the same linear combination of the six row primitives. So the combined route
compresses six statements into one large one without compressing the analytic
dependency graph, and is harder to debug and reuse.

Order:

1. a few shared master endpoint theorems, living in `Problem24Euler.lean`
   alongside the existing kernel library;
2. each `I_ab` derived in `Tminus`/`Tplus` normal form (this is where the
   `K = 2((7/8)ζ4 - Tminus - Tplus - (1/4)ζ2²)` reduction does its work — no new
   weight-4 analysis);
3. the fully expanded PSLQ forms by `unfold` + `ring`;
4. finish with `quadAltSixIntegralLinear` / `quadAltCoeffIntegral_eq_six`, both
   already proved and axiom-clean.

A direct combined identity is still worth adding afterwards as a regression
cross-check — but not as the main proof.

The one genuinely new analytic step in the whole of Layer E remains the termwise
integration behind `K` (expand `log(1+x)/(1+x)`, integrate `∫₀¹ xⁿ log²x = 2/(n+1)³`,
reindex with `H_n = H_{n+1} - 1/(n+1)`). Everything else is rational arithmetic
over constants this repo already owns.

## The one remaining analytic step, stated exactly

Numerically pinned (agreement to the truncation error of a 3000-term alternating
series with the adjacent-average correction, 1.5e-13):

```lean
∫ x in (0:ℝ)..1, (Real.log x)^2 * Real.log (1+x) / (1+x)
  = 2 * ∑' k : ℕ, (-1)^k * harmonicNumber (k+1) / ((k:ℝ)+2)^3
```

Nat-indexed from 0, no shift needed. Ingredients:
`log(1+x)/(1+x) = Σ_{n≥1} (-1)^{n+1} H_n xⁿ` and `∫₀¹ xⁿ (log x)² dx = 2/(n+1)³`.

The series is ALTERNATING, so Tonelli does not apply — the interchange needs
`MeasureTheory.hasSum_integral_of_summable_integral_norm` (the same lemma Layer C
already uses), whose hypotheses are per-term integrability plus summability of
`∫‖·‖`. Here `∫₀¹ |xⁿ (log x)² H_n| dx = 2H_n/(n+1)³`, and `Σ H_n/(n+1)³`
converges, so the norm-summability hypothesis is met.

Then `K = 2(η(4) - A)` with `A = Σ (-1)^{n-1} H_n/n³`, and
`A = Tminus + Tplus + (1/4)ζ2²` closes it against constants the repo owns.

## Definition audit (done — numerics verified against the Lean objects)

Every formula used in this session's numerical work was checked character by
character against the Lean definitions: `dilog` (which is literally
`∑' n, z^(n+1)/(n+1)^2`, i.e. Li₂), `W0`, `H1`, `H2`, `quadAltV`,
`quadAltDminus`, `quadAltMclosed` (all seven terms including
`- 2 * dilog ((1+x)/2)`), `quadAltJclosed`, `quadAltQclosed`, `quadAltR`,
`quadAltFclosed`, and `alternatingQuadraticEulerValue24`.

This matters because it is the error class a numerical check cannot catch: had
the Python used a different branch, normalisation or endpoint convention, every
digit could still have matched while the Lean statement said something else. It
did not — the 25-to-60-digit agreements above are agreements about the Lean
objects.

---

## Layer E status after the 2026-08-01 build-out

**Discharged.** `hK` is gone: `quadAltK_eq` proves
`K = (1/5) ζ₂² − 2·Tminus − 2·Tplus`, axiom-clean.  The route was not the
planned one.  `K`'s termwise step (`quadAltK_hasSum`) came from
`harmonicNumber_generating_hasSum`, which the repo already had — at `x = −t` it
IS the generating function of the `K` integrand.  Then `K`'s value came from two
more theorems the repo already had, buried as private steps inside the
cubic-linear constants: with `a n = (−1)^{n+1}H_{n+1}/(n+1)³` and
`b n = (−1)^{n+1}/(n+1)⁴`, the `K` summand is exactly `2(a(n+1) − b(n+1))`, and
`(a−b)(0) = 0` so the shift is free.  Two of the three ingredients were already
in the repo; both were found by grepping, not by deriving.

**Definition audit, complete and clean.** Checked character by character against
the numerics: `dilog z = Σ z^{n+1}/(n+1)²` (= Li₂, not Rogers),
`polylog4 z = Σ z^{n+1}/(n+1)⁴`, `W0 t = π²/6 − 2·dilog(t/2) − log(t/2)²`,
`H1 = −log(1−t)`, `H2 = −log(1−t/2)`, the six `I_ab`, and `quadAltK`.  All match.

**The frontier statement is true, and the sharp test is at index 2.** An
adversarial audit flagged that `parityRemainder24 n = H_n − 2A_n`, not the
alternating partial sum `A_n`.  It is not a defect: `parityRemainder24` is an
internal reduction object, not a problem-statement object (the competition asks
about `ΣΣ C(m,k)²H_k²/((m+1)²C(2m,m))`).  Decided numerically —
repo definition gives `tsum = −0.06236615…` against the target `−0.06236610…`,
the alternative gives `−0.23638…`.  The cheap discriminator the audit itself
proposed also lands on the repo: `alternatingQuadraticEulerTerm24 2 = 4/27`, not
`2/27`.  Indices 0 and 1 cannot tell the two apart — only index 2 and beyond can.

**Negative result worth knowing before formalizing I11.**  `H1(t)/(1−t) = Σ H_m tᵐ`
makes `I11 = Σ_m H_m ∫₀¹ W0(t) tᵐ dt` the natural series.  It is NOT termwise
equal to `Σ P_m/m³`:

```
m         H_m·∫₀¹W0·tᵐ        P_m/m³
1         -0.38629436         -1.0
2         -0.18172581          0.0625
3         -0.09709522          0.00617284
```

Both partial sums converge to `Tplus = −0.886852…`, so the identity is real, but
it needs an actual rearrangement — `congr_fun` against the repo's
`cubicLinearEulerTerm24_hasSum` will not close it.  Nor do the repo's `Tplus`
kernels help directly: they are trilogarithm-based
(`(ζ₃ − Li₃(x))/(1−x)` and `(ζ₃ − Li₃(−x))/(1+x)`), not `W0`-based.

**Remaining.** Six hypotheses, the rows `I10 … I22`.

## The rows are much cheaper than the Li₄-basis derivations suggest

The published analytic derivations all route through dilog/trilog functional
equations at `1/2` and produce `Li₄(1/2)` by way of an antiderivative — which in
Lean would need the derivative of `Li₄`, which neither this repo nor Mathlib
v4.29 has.  That is the expensive route and it is not necessary.

Substituting `u = t/2` FIRST collapses `W0`: `log²(t/2)` becomes `log²u`, so

```
W0(t) = ζ₂ − 2 Li₂(u) − log²u        (u = t/2)
```

and the remaining factor is an exact differential in every case:

```
H1(t)/t  dt = −log(1−t)/t dt  = d[Li₂(t)]
H2(t)/t  dt = −log(1−u)/u du  = d[Li₂(u)]
H2(t)/(2−t) dt = −log(1−u)/(1−u) du = d[log²(1−u)/2]
```

so the `Li₂` cross term integrates by the chain rule, `∫ Li₂ d Li₂ = Li₂²/2`,
with no series at all.  Two rows come out immediately, both verified as exact
symbolic identities against the recorded table at 41 digits:

```
I20 = ζ₂·Li₂(½) − Li₂(½)² − log²2·Li₂(½) − 2 log2·Li₃(½) − 2 Li₄(½)
I10 = ζ₂² − 2 V₂ − 2 ζ₄ − 2 log2·ζ₃ − log²2·ζ₂ ,   V₂ = Σ_{m≥1} H_m/(2^m m³)
```

`I20` needs no Euler sum whatsoever — only `Li₂(½)` (already proved here as
`quadAlt_dilog_half`), `Li₃(½)`, and `Li₄(½)`.  `I10` needs one, `V₂`, plus the
three elementary moments `Σ(1/n)(2/n³) = 2ζ₄`, `Σ(1/n)(1/n²) = ζ₃`,
`Σ 1/n² = ζ₂`.

The remaining three (`I12`, `I21`, `I22`) go the same way but land on the two
half-interval moments `J12 = ∫₀^{1/2} log u log²(1−u)/u du` and
`D3 = ∫₀^{1/2} log³(1−u)/u du`, whose series expansions produce `V₃` and `V₄`.
Both are confirmed to close in the basis — see
`scripts/p24_half_argument_basis.py`.  Consistency check passed: the direct
route gives `I22 = ζ₂L²/2 − Li₂(½)L² − L⁴/2 + J12 − D3`, and the earlier IBP
route gives `I22 = J12 − D3`; these agree exactly because
`Li₂(½) = ζ₂/2 − L²/2`.

**Order of attack, cheapest first: I20, I10, then I22, I12, I21.**

## The six rows collapse to ONE integral plus ONE row

The coefficient vector `(-2,-2,2,4,6,-5)` is not arbitrary — it comes out of the
Möbius substitution — and it carries structure that row-by-row evaluation throws
away.  Four of the six kernels are exact differentials:

```
H1/t     dt = d[Li₂(t)]        H1/(1-t) dt = d[H1²/2]
H2/t     dt = d[Li₂(t/2)]      H2/(2-t) dt = d[H2²/2]
```

The leftover `2 H1/(2-t) + 6 H2/(1-t)` is NOT exact, but it is exact *up to a
smaller multiple of one kernel*: using
`d[log(1-t) log(2-t)] = -log(2-t)/(1-t) - log(1-t)/(2-t)` and
`H2 = log 2 - log(2-t)`, the coefficient 6 drops to 4.  Writing

```
Φ(t) = -2 Li₂(t) - H1(t)² + 4 Li₂(t/2) - (5/2) H2(t)²
       + 2 log(1-t) log(2-t) + 2 log2 · H1(t)
```

the whole bracket satisfies, verified at 30 digits at six interior points,

```
G(t) = Φ'(t) + 4 H2(t)/(1-t)
```

Both boundary terms of `∫ W0 dΦ` vanish — `W0·Φ` is `3e-14` at `t = 1e-8` and
`6e-14` at `t = 1-1e-8`, decaying quadratically — because `W0` has a double zero
at `1` and `Φ = O(t)` at `0`.  Hence, with `W0' = -2R/t`,

```
S  :=  -2I10 - 2I11 + 2I12 + 4I20 + 6I21 - 5I22
    =  2 ∫₀¹ Φ(t) R(t)/t dt  +  4 I21
```

confirmed numerically to 30 digits against the target constant.

**Six obligations become two.**  And the surviving integral no longer contains
`W0`: `R(t) = log t - log(2-t)` is elementary, so every piece of `∫ Φ R/t` is
`(elementary or Li₂) × log / t` — exactly the shape the repo's termwise
template eats.  One of the six pieces, `2 log2 · ∫ H1 R/t`, is only weight three.

This is an alternative to the row-by-row route, not a replacement: keep whichever
lands first.  Row-by-row already has `I20` and `I10` in closed form and `I11`
under way.

### The Φ-route ledger, closed

Writing `P = ∫₀¹ Φ R/t = -2P1 - P2 + 4P3 - (5/2)P4 + 2P5 + 2L·P6` with

```
P1 = ∫₀¹ Li₂(t) R/t        P2 = ∫₀¹ H1² R/t       P3 = ∫₀¹ Li₂(t/2) R/t
P4 = ∫₀¹ H2² R/t           P5 = ∫₀¹ log(1-t)log(2-t) R/t
P6 = ∫₀¹ H1 R/t            (weight three)
```

every piece is in the basis `(Li₄(½), L⁴, L²ζ₂, Lζ₃, ζ₂²)`:

```
P1 = (-1,  -1/24, -1/2,  -7/8,   1/20)
P2 = ( 0,   0,     0,    -7/2,   3/4  )     = I11
P3 = (-1,  -1/24,  1/4,  -7/8,   1/8  )
P4 = (-6,  -1/4,   3/2, -21/4,  23/10 )     = I22
P5 = (-6,  -1/4,   9/4, -35/8,  87/40 )
P6 = -(3/2) L ζ₂                            (weight three)
```

`P2` and `P4` are `I11` and `I22` in their post-IBP form, which is a consistency
check, not new information.  `P5` was obtained by solving the linear relation and
then confirmed against its own direct quadrature to 16 digits — independent, not
circular.  `P1` likewise cross-checks to 17 digits.

`P1` has a short derivation worth keeping: the inner sum has a closed form,
`Σ_{m≥1} 1/(m²(m+n)) = ζ₂/n - H_n/n²`, so
`P1 = -ζ₄ - Lζ₃ + ζ₂ Li₂(½) - V₂` with `V₂ = Σ H_n/(2ⁿn³)`.

**Method warning, paid for the hard way.**  A first pass PSLQ'd `P1, P3, P5, P6`
straight from quadrature and reported NO RELATION for all four.  That was a
precision artifact: two different quadrature schemes for the same integral
disagreed at the 17th digit, so the input never had the ~40 digits a five-term
relation needs.  Recomputing by series — using the closed form above for the
inner sum — gave clean relations immediately.  This is the same trap already
recorded in UNDERSTANDING.md under "substitute before trusting PSLQ"; the lesson
generalises past substitution to *any* claim that a constant is not in a basis.
Cross-check two independent evaluations before believing a negative.

### …but the FULL collapse is circular — do not pursue it

The leftover `2H1/(2-t) + 6H2/(1-t)` IS exact, with primitive
`2 H1 H2 + 4L·H1 - 4 Li₂(t-1) - 2ζ₂`.  So the whole bracket has a primitive and
`S = ∫ W0 dΨ` for a single `Ψ`.  That looks like it collapses everything to one
integration by parts — and it does, straight back to where we started:

```
Ψ(t) = -quadAltJclosed(-t/(2-t))
```

i.e. the combined primitive is exactly the Möbius pullback of the Layer C
integrand.  Integrating by parts against `W0' = -2R/t` therefore reverses the
first IBP + Möbius step and returns the Layer C coefficient integral.  Nothing is
evaluated.

This is worth stating because the collapse is genuinely attractive and the
circularity is not visible until you identify `Ψ`.  The partial collapse recorded
above (`S = 2∫ΦR/t + 4 I21`) escapes it only because it stops short of full
exactness.

**Conclusion: row-by-row is the route.**  The useful collapsed statement would be
`sixIntegralCombination = -(11/2)K - 12 L²ζ₂ + (35/2) Lζ₃ - (31/20) ζ₂²` — but
that is not a consequence of exactness; it IS the whole weight-four reduction
written on one line.

### Cross-row relations: I12 and I21 need no independent evaluation

With `I10, I11, I20, I22` banked, two elementary relations (R12), (R21) obtained
from a Möbius difference give `I12` and `I21`.  They need one extra integral, a
GROUPED packet that is elementary in the unit-interval alphabet `{0, 1, -1}` — no
half-argument endpoint, no `Li₄` derivative.

Two warnings that came with it, both matching traps already hit here:
- **The packet must stay grouped.**  Splitting it under `1/(1-x)` produces two
  individually divergent integrals; only the sum is integrable.
- **Do not expand its three terms into series separately.**  That produces
  coloured depth-three Euler sums outside the four-item half-argument list, which
  would break the series-only plan.

### I10 reduces to a single evaluation, and V₂ is tied to the already-proved K

Splitting `W0` rather than integrating by parts,

```
I10 = ζ₂·∫₀¹ H1/t  -  2·∫₀¹ Li₂(t/2) H1/t  -  ∫₀¹ log²(t/2) H1/t
    = ζ₂·ζ₂  -  2 V₂  -  (2ζ₄ + 2Lζ₃ + L²ζ₂)
```

with `V₂ = ∫₀¹ Li₂(t/2) H1(t)/t dt = Σ_{m≥1} H_m/(2^m m³)`.  The first and third
integrals are termwise integrations against moments already proved here
(`∫₀¹ t^{n-1} = 1/n`, `∫₀¹ t^{n-1} log t = -1/n²`, `∫₀¹ t^{n-1} log²t = 2/n³`).
`V₂` needs `∫₀¹ t^{m-1}(-log(1-t)) dt = H_m/m` and then its own value.

That value is tied to `K`, which is already proved:

```
V₂ = (1/4) K + (1/4) L²ζ₂ - Lζ₃ + (17/40) ζ₂²
```

verified numerically.  Substituting it back reproduces the recorded K-normal-form
row `I10 = -(1/2)K - (3/2)L²ζ₂ - (13/8)ζ₄` exactly, which is a consistency check
on both.

**Do not mistake this for a free reduction.**  `K`'s series is alternating and
full-argument, `Σ(-1)ⁿ H_{n+1}/(n+2)³`; `V₂`'s is non-alternating and
half-argument, `Σ H_m/(2^m m³)`.  The relation between them is a genuine
weight-four identity, not a rearrangement.  So `I10` costs exactly one
evaluation — either `V₂` directly (half-interval log moments, per the analytic
derivation already in hand) or the `V₂ ↔ K` identity — whichever is cheaper in
Lean.  It does not cost zero, and it does not need any constant beyond `K`.

### Carried-hypothesis audit, done by instantiation rather than by reading

A hostile audit raised two dangers that `#print axioms` cannot see:

1. the coefficient-vector theorem might be about a *named* integral that is not
   literally the one the series bridge produces;
2. the seven-hypothesis conditional theorem might be vacuous, or hide an eighth
   obligation.

Both were closed by making Lean answer, not by reading the statements:

```lean
example : (∑' n, alternatingQuadraticEulerTerm24 n)
    = -2*I10 - 2*I11 + 2*I12 + 4*I20 + 6*I21 - 5*I22 :=
  quadAlt_tsum_eq_coeff_integral.trans quadAltCoeffIntegral_eq_six
```

`trans` typechecks, so the two theorems are about the *same* integral — Lean
unified them syntactically. Danger 1 is closed.

```lean
example (h10 : I10 = …) (h12 : I12 = …) (h21 : I21 = …) :
    sixIntegralCombination = bridgeValue :=
  quadAltSixIntegral_eq_bridgeValue quadAltK_eq h10 quadAltI11_eq h12
    quadAltI20_eq h21 quadAltI22_eq
```

Feeding in every hypothesis that is already proved leaves exactly three, and the
application typechecks. So the conditional theorem is not vacuous, there is no
hidden eighth obligation, and the endgame shape is settled: **I10, I12, I21 land
and the capstone closes.**

This is the right way to discharge the "carried hypotheses read and judged
satisfiable" clause in UNDERSTANDING.md — instantiate what you have and let the
kernel report the remainder, rather than reading the statement and forming an
opinion.

### The last two obligations contain no polylogarithm at all

`I12` and `I21` follow from `(R12)`, `(R21)`, which as single integrals are

```
A := ∫₀¹ W0 (H1-H2)(1/(2-t) + 1/t) dt
B := ∫₀¹ [ W0 H2 (1/(1-t) - 1/(2-t)) - 2 W0 (H1-H2)/t + (1/2) W0 H1/(1-t) ] dt
```

and their values are, exactly,

```
A = -η(2)² - ζ(4)          = -(13/20) ζ₂²
B = (1/2) η(2)² + 2 ζ(4)   =  (37/40) ζ₂²
```

with `η(2) = ζ₂/2`.  **No `Li₄(1/2)`, no `log⁴2`, no `log2·ζ₃`, no polylogarithm
of any kind** — only the alternating zeta value at 2 and `ζ(4) = π⁴/90`, the
latter already proved here as `shifted_zeta_four_hasSum24`.

The mechanism, which is why this is not luck: after the Möbius substitution `(A)`
becomes an alternating OFF-DIAGONAL SQUARE, so its value is `(Σ(-1)ⁿ/n²)²`
together with `ζ(4)`.  For `(B)` the same Möbius coordinate plus the
endpoint-reversing involution `x ↦ (1-x)/(1+x)` makes the one possible
half-endpoint weight-four period appear with the SAME coefficient in two
logarithmic terms, so it cancels before integration rather than after.

That distinction is the whole point for formalization: a proof that produces
`Li₄(1/2)` and cancels it needs machinery this repo does not have (no `Li₄`
derivative).  A proof that never produces it needs none.

Both values verified exactly (agreement at 1e-31 against the recorded rows) and
the integrals themselves to 15 digits.
