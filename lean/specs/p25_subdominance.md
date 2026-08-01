# Task: Prove catalanError is subdominant (close the last P2.5 sorry)

## The one sorry

In `RamanujanChallenge/Problem25Moment.lean` line 191:

```lean
theorem catalanError_over_denominator_tendsto_zero (j : Fin 3) :
    Filter.Tendsto (fun N => catalanError N j / (denominator N j : ℝ))
      Filter.atTop (nhds 0) := by
  sorry
```

This is the ONLY sorry blocking `problem25_solved : Problem25Claim`.

## What this means mathematically

The three challenge ratios `P_n/Q_n` converge to `commonLimit` (proved via Hilbert
metric contraction in Problem25Connection.lean). We need `commonLimit = catalanConstant`.

Equivalently: `catalanError N j / denominator N j = catalanConstant - P_n/Q_n → 0`.

This says the catalanError grows STRICTLY SLOWER than the denominator (subdominant).

## What's already proved

1. **catalanConstant_eq_integral** (Problem25Integral.lean, 0 sorry):
   `catalanConstant = ∫₀¹ (-log t)/(1+t²) dt`

2. **integral_neg_log_01** (Problem25Moment.lean):
   `∫₀¹ (-log t) dt = 1`

3. **Sign pattern at N=0** (Problem25.lean):
   - `positiveCatalanError_zero_zero_neg`: pCE 0 0 < 0
   - `positiveCatalanError_zero_one_pos`: 0 < pCE 0 1
   - `positiveCatalanError_zero_two_pos`: 0 < pCE 0 2

4. **Sign pattern at N=1** (Problem25SignN1.lean, 0 sorry):
   - `positiveCatalanError_one_zero_neg`: pCE 1 0 < 0
   - `positiveCatalanError_one_two_pos`: 0 < pCE 1 2

5. **Tight Catalan bounds** (Problem25TightBounds.lean, 0 sorry):
   - `catalan_tight_lower`: 1590511050/1736437500 < catalanConstant
   - `catalan_tight_upper`: catalanConstant < 21390206625/23352603750

6. **Projective contraction** (Problem25Connection.lean):
   - All three ratios converge to `commonLimit` at 2/3 geometric rate
   - `positiveRatio_envelope`: ratios are trapped in contracting envelopes
   - `problem25Claim_iff_commonLimit_eq_catalan`: Problem25Claim ↔ commonLimit = catalanConstant

7. **Recurrence** (Problem25.lean):
   - `positiveCatalanError_succ`: pCE evolves by the same positive matrix as pDen
   - `positiveCatalanError_eq`: pCE N j = catalanConstant * pDen N j - pNum N j
   - All positiveMatrix entries are positive (positiveMatrix_pos)

## Why this is hard

The error ratio `e_j = catalanError_j / denominator_j` satisfies:
`e_j(n+1) = Σ_i w_{i,j} · e_i(n)` (weighted average, weights > 0, sum to 1)

This proves `e_j → L` for some `L = catalanConstant - commonLimit`. But it CANNOT
prove `L = 0` — the weighted average converges but doesn't determine the limit.

Numerical evidence: sign pattern (-, +, +) persists through N=100+,
ε_{k+1}/ε_k → -0.12 (alternating decay). All consistent with L = 0.

## Possible approaches (ranked by promise)

### A. Prove the sign pattern persists at ALL N (squeeze argument)
If `∀ N, ∃ j, pCE N j < 0` AND `∀ N, ∃ j, pCE N j > 0`, then
G ∈ [lowerEnvelope N, upperEnvelope N] for all N, and since both → commonLimit,
G = commonLimit by squeeze.

To prove this: show the (-, +, +) sign pattern is an invariant of the positive
matrix multiplication. This requires showing that for each N, the j=0 column
weight of the negative entry dominates. Numerically verified through N=100.

### B. Operator factorization
The 3rd-order scalar ODE for the challenge has been extracted (degree-27 polynomial
coefficients in n). If it factors as L1 · L2 (with L1 first-order), then the error
satisfies L2 (lower order), giving subdominant growth automatically.

### C. Explicit second solution
If we can exhibit a SPECIFIC second solution S_n of the recurrence (with growth rate 1
or (17-12√2)^n) and show catalanError is proportional to it, then subdominance follows.
The initial conditions pin the proportionality constant.

### D. Wronskian / determinant argument
The Wronskian W(Q, P, E) = 0 trivially (E = G*Q - P). But a Casorati determinant
of Q, E, S (where S is known subdominant) might give useful information.

### E. Direct norm bound from integral representation
catalanError = ∫ kernel · R(t²) where kernel = (-log t)/(1+t²) and R is affine.
The integral = G*Q - P by construction (tautological). But if we can show |R(t²)|
is bounded by C · ρ^N on [0,1] for some ρ < 17+12√2, that gives subdominance.
This requires showing the remainder polynomial coefficients (Q-P and P) have
specific cancellation — essentially the Padé theory content.

## Build & verify
```bash
~/.elan/toolchains/leanprover--lean4---v4.29.0/bin/lake env lean RamanujanChallenge/Problem25Moment.lean
```
Expected: 0 sorry warnings (currently 1).

## Files you may modify
- `RamanujanChallenge/Problem25Moment.lean` (the sorry is here)
- You may create new helper files if needed

## Do NOT modify
- Problem25.lean, Problem25Connection.lean, Problem25Integral.lean,
  Problem25TightBounds.lean, Problem25SignN1.lean (all proven, 0 sorry)
