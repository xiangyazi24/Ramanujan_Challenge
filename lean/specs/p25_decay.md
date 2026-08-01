# Task: Prove catalanError subdominance for Problem 2.5

## Context
Working in `/Users/huangx/repos/Ramanujan_Challenge/lean/`.
Lean v4.29.0 + Mathlib v4.29.0.

In `RamanujanChallenge/Problem25Moment.lean` there is one remaining sorry:

```lean
theorem catalanError_over_denominator_tendsto_zero (j : Fin 3) :
    Filter.Tendsto (fun N => catalanError N j / (denominator N j : ℝ))
      Filter.atTop (nhds 0) := by
  sorry
```

## Key definitions (all in Problem25.lean)
- `catalanError N j = (denominator N j : ℝ) * catalanConstant - (numerator N j : ℝ)`
- `catalanConstant = ∑' n : ℕ, (-1 : ℝ) ^ n / (2 * (n : ℝ) + 1) ^ 2`
- `denominator` and `numerator` are integer sequences from a 3×3 matrix recurrence
- `challengeRatio N j = numerator N j / denominator N j`
- `positiveRatio`, `positiveDenominator`, `positiveNumerator` are sign-conjugated versions
- Already proved: `positiveMatrix_pos`, `positiveDenominator_pos`, `positiveNumerator_pos`

## Key existing theorems
In Problem25Connection.lean:
- `challengeRatio_tendsto_common j`: all three ratios → `commonLimit`
- `lowerEnvelope_tendsto_common`: lowerEnvelope → commonLimit
- `upperEnvelope_tendsto_common`: upperEnvelope → commonLimit
- `positiveRatio_envelope n j`: lowerEnvelope n ≤ positiveRatio n j ≤ upperEnvelope n
- `problem25Claim_iff_commonLimit_eq_catalan`: Problem25Claim ↔ commonLimit = catalanConstant

In Problem25.lean:
- `positiveCatalanError_zero_zero_neg`: positiveCatalanError 0 0 < 0
- `positiveCatalanError_zero_one_pos`: 0 < positiveCatalanError 0 1
- `positiveCatalanError_zero_two_pos`: 0 < positiveCatalanError 0 2
- `catalan_lower_bound`: 32972/36000 < catalanConstant
- `catalan_upper_bound`: catalanConstant < 30921/33750
- `catalanMagnitude_antitone.alternating_series_le_tendsto`
- `catalanMagnitude_antitone.tendsto_le_alternating_series`

## Approach: Numerical squeeze

The challenge ratios converge to `commonLimit`, and we need to show `commonLimit = catalanConstant`.

Strategy: prove tighter Catalan bounds (using more terms of the alternating series), then show G is in the contracting envelope at step N=1.

### Step 1: Tighter Catalan bounds
Prove `catalan_tight_lower` and `catalan_tight_upper` using k=130 (260 terms):
```lean
theorem catalan_tight_lower : (X : ℝ) / Y < catalanConstant := by
  calc (X : ℝ) / Y < ∑ i ∈ Finset.range (2 * 130), (-1 : ℝ) ^ i * catalanMagnitude i := by
    set_option maxHeartbeats 0 in norm_num [catalanMagnitude, Finset.sum_range_succ]
  _ ≤ catalanConstant :=
    catalanMagnitude_antitone.alternating_series_le_tendsto
      catalanMagnitude_partialSum_tendsto 130
```
The rational X/Y should be chosen so that X/Y > lowerEnvelope 1.

### Step 2: Show G ∈ envelope at N=1
Compute the three positiveRatio values at N=1 exactly (from `approximants_one` or by extending it). Show:
- `catalan_tight_lower` > min ratio at N=1 (= lowerEnvelope 1)
- `catalan_tight_upper` < max ratio at N=1 (= upperEnvelope 1)

### Step 3: Squeeze
Since G ∈ [lowerEnvelope n, upperEnvelope n] for n=0 (proved) and n=1 (step 2),
and [lowerEnvelope 1, upperEnvelope 1] ⊂ [lowerEnvelope 0, upperEnvelope 0],
and both envelopes contract to commonLimit:
lowerEnvelope 1 ≤ G ≤ upperEnvelope 1
lowerEnvelope 1 ≤ commonLimit ≤ upperEnvelope 1
Width ≤ 4.5e-6 at N=1

This gives |G - commonLimit| ≤ 4.5e-6 but NOT exact equality.

### Alternative: Prove at more steps
If G ∈ envelope at N=0,1,2,...,K for some K, and the envelope width → 0,
then G = commonLimit. But showing G ∈ envelope at each step requires
the sign pattern (-, +, +) of positiveCatalanError to persist, which
requires computing the signs at each step.

### Your job
1. First, try to compute positiveCatalanError signs at N=1 (by computing approximants_one and using the existing Catalan bounds). If you can show the sign pattern persists at N=1, try N=2.

2. If the sign pattern persists for enough steps that the envelope width drops below the precision of the Catalan bounds, you can squeeze G = commonLimit.

3. The key identity: catalanError N j / denominator N j = catalanConstant - challengeRatio N j. So if challengeRatio → catalanConstant, then the ratio → 0.

4. To PROVE challengeRatio → catalanConstant: you need commonLimit = catalanConstant. So this IS the gap.

5. Most promising: prove that positiveCatalanError has both signs at EVERY step (by induction, using the matrix structure). If min(positiveCatalanError n j) < 0 < max(positiveCatalanError n j) for all n, then catalanConstant ∈ [lowerEnvelope n, upperEnvelope n] → catalanConstant = commonLimit.

## Build
```bash
~/.elan/toolchains/leanprover--lean4---v4.29.0/bin/lake env lean RamanujanChallenge/Problem25Moment.lean
```

## Rules
- You may modify ONLY Problem25Moment.lean
- Target: replace the sorry in `catalanError_over_denominator_tendsto_zero`
- If you can't prove it fully, try to reduce the sorry to something smaller
- Use `set_option maxHeartbeats 0` freely
