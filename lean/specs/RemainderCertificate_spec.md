# Codex Spec: RemainderCertificate.lean — Close All Sorries

## Target
Close all sorry's in `RamanujanChallenge/RemainderCertificate.lean`. There are 2 sorry'd theorems:
1. `tendsto_div_of_remainder_tendsto_zero` — the master convergence lemma
2. `tendsto_zero_of_geometric_bound` — geometric decay implies convergence to 0

The recurrence-linearity lemmas (`remainder_satisfies_recurrence_*`) already close by `ring`.

## Working directory
`~/repos/Ramanujan_Challenge/lean/`

## Strategy

### `tendsto_div_of_remainder_tendsto_zero`

Given: q n ≠ 0 eventually, and (q n * C - p n) → 0.
Prove: p n / q n → C.

**Route:** Rewrite p n / q n = C - (q n * C - p n) / q n. Then:
- The second term → 0/q n. But we need q n to be eventually bounded away from 0... 
  Actually no: we just need q n ≠ 0 eventually and the numerator → 0.
- Use `Filter.Tendsto.div` or `Filter.Tendsto.div_atTop` from Mathlib.
- Actually cleaner: `p n / q n = C - (q n * C - p n) / q n`. 
  We need `(q n * C - p n) / q n → 0`. This follows from:
  - If r n → 0, and we can show |r n / q n| ≤ |r n| (when |q n| ≥ 1), then done.
  - But actually for the general case: use `Tendsto.div` from Mathlib's filter library.
    `Tendsto r atTop (nhds 0)` and `q n ≠ 0` eventually gives us what we need IF
    we know something about q n (it doesn't tend to 0).
  
**Simpler rewrite approach:**
For integer-valued q (which is the case in all our problems), q n ≠ 0 implies |q n| ≥ 1.
But the lemma is stated for ℝ. So use the general approach:

Actually, the cleanest proof: the hypothesis says q n * C - p n → 0.
So p n = q n * C - (q n * C - p n). 
Thus p n / q n = C - (q n * C - p n) / q n.

We need: (q n * C - p n) / q n → 0.

Use: if f → 0 and g is eventually nonzero AND g doesn't tend to 0... 
Actually this is NOT true in general (g could oscillate between large and small values).

**Better approach for our specific use case:** Add a hypothesis that |q n| → ∞ 
(which is true for all our problems since q n are integer sequences going to infinity).
OR: use the fact that p n / q n = C - r n / q n, and if r n → 0 and |q n| ≥ 1 eventually,
then |r n / q n| ≤ |r n| → 0.

**Recommended:** Strengthen the lemma slightly: add `(hq_bound : ∀ᶠ n in atTop, 1 ≤ |q n|)`.
This is trivially satisfied for all our integer-valued sequences.
Then |r n / q n| ≤ |r n| / 1 = |r n| → 0.

Or keep the current signature and use `Tendsto.div`:
Mathlib has `Filter.Tendsto.div` which requires the denominator filter to tend to some nonzero value.
If we don't know where q tends, we need a different approach.

**Best approach:** Just require `Filter.Tendsto (fun n => |q n|) atTop atTop` (q diverges).
Then use `tendsto_div_nhds_zero_of_tendsto_nhds_zero_of_tendsto_atTop`:
  if numerator → 0 and |denominator| → ∞, then ratio → 0.

Actually, look for `isLittleO_one_iff` or `Asymptotics.IsLittleO` in Mathlib.
Or just: r n = o(q n) since r n → 0 and q n → ∞ gives r n / q n → 0.

Use `Filter.Tendsto.div_atTop` if available, or construct manually.

### `tendsto_zero_of_geometric_bound`

Given: |r (n+1)| ≤ ρ |r n| for n ≥ N, with 0 ≤ ρ < 1.
Prove: r n → 0.

**Route:**
1. By induction: |r (N+k)| ≤ |r N| · ρ^k for all k.
2. |r N| · ρ^k → 0 since ρ < 1 (use `tendsto_pow_atTop_nhds_zero_of_lt_one`).
3. Squeeze: |r n| ≤ |r N| · ρ^(n-N) → 0.
4. Use `squeeze_zero_norm` or `tendsto_of_tendsto_of_tendsto_of_le_of_le`.

## Hard rules
- No sorry, no axiom, no native_decide
- Keep the existing `ring`-proved recurrence lemmas unchanged
- If the current signature of `tendsto_div_of_remainder_tendsto_zero` is too weak,
  you may strengthen it (add hypotheses like `Tendsto (|q ·|) atTop atTop`)
  BUT update the docstring to reflect the change.

## Verification
```bash
export PATH="$HOME/.elan/bin:$PATH"
lake env lean RamanujanChallenge/RemainderCertificate.lean 2>&1
```
Must compile with no errors and no sorry warnings.

## Stall protocol
If stuck, deliver what compiles + precise stall report.
