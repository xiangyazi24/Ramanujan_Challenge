# Codex Spec: Problem 2.3 (π + e) — Close All Sorries

## Target
Close all 3 sorry'd theorems in `RamanujanChallenge/Problem23.lean`:
1. `lambertRatio_tendsto` — B_m/A_m → π/4
2. `derangementRatio_tendsto` — m!/D_m → e
3. `problem23_limit` — P_m/Q_m → π + e

## Assume-success frame
Assume this target is achievable. Your job is to find the path, not to judge whether one exists.

## Working directory
`~/repos/Ramanujan_Challenge/lean/`

## Strategy

### Part 1: Lambert CF convergence (B_m/A_m → π/4)

The Lambert continued fraction `B_m/A_m` with `X_m = (2m+1)X_{m-1} + m²X_{m-2}` converges to π/4.

**Route: Euler's CF for arctan(1).**
The classical identity: `arctan(x) = x/(1 + x²/(3 + (2x)²/(5 + (3x)²/(7 + ...))))`.
At x=1: `π/4 = 1/(1 + 1/(3 + 4/(5 + 9/(7 + ...))))`.

This is equivalent to our Lambert recurrence: the convergents of this CF satisfy
exactly `X_m = (2m+1)X_{m-1} + m²X_{m-2}`.

**Proof approach:**
1. Define the partial convergent ratio `B_m/A_m` as a function ℕ → ℝ.
2. Prove that the Lambert recurrence equivalences (the odd-CF wallis-type identity):
   `A_{m+1}·B_m - A_m·B_{m+1} = (-1)^{m+1} · (m!)²` by induction.
3. This gives `B_{m+1}/A_{m+1} - B_m/A_m = (-1)^{m+1}·(m!)²/(A_m·A_{m+1})`.
4. Show A_m grows like (2m)!!·C via crude bound: A_m ≥ (2m+1)!! (from recurrence positivity).
5. The alternating series with terms decaying geometrically converges.
6. The limit equals arctan(1) = π/4. Use `Real.arctan_one` from Mathlib.

**Alternative simpler route:** Prove by induction that
`B_m/A_m = ∑_{k=0}^{m-1} (-1)^k / (2k+1) + error_m`
where the partial sums of the Leibniz series `∑ (-1)^k/(2k+1)` are exactly the CF convergents.
Then the limit is `∑' k, (-1)^k/(2k+1) = arctan(1) = π/4`.

For the infinite series = π/4 step, you may need to prove this from the Taylor series of arctan.
Mathlib has `Real.arctan_one : arctan 1 = π / 4` and `Real.hasSum_pow_mul_geometric_of_abs_lt_one`
and related series tools.

### Part 2: Derangement ratio (m!/D_m → e)

The derangement numbers D_m satisfy D_0=1, D_1=0, D_m=(m-1)(D_{m-1}+D_{m-2}).

**Key identity:** D_m = m! · ∑_{k=0}^{m} (-1)^k/k!

Prove this by induction from the recurrence. Then:
- m!/D_m = 1/(∑_{k=0}^m (-1)^k/k!)
- The partial sums ∑_{k=0}^m (-1)^k/k! → e^{-1} (this is the Taylor series of exp(-1))
- Therefore m!/D_m → 1/e^{-1} = e = Real.exp 1

For the exp series convergence, use Mathlib's `Real.hasSum_exp` or
`NormedSpace.exp_hasFPowerSeriesOnBall` to extract that
∑' k, (-1)^k/k! = exp(-1).

Use `Filter.Tendsto.inv₀` with `Real.exp_pos` to invert.

### Part 3: Combining (P_m/Q_m → π + e)

From the file:
```
challengeP m = 4 * lambertB m * derangement m + lambertA m * m!
challengeQ m = lambertA m * derangement m
```

So: `challengeP m / challengeQ m = 4 * (lambertB m / lambertA m) + m! / derangement m`

Use `Filter.Tendsto.add` after proving the decomposition holds eventually
(need lambertA m ≠ 0 and derangement m ≠ 0 eventually).

Prove `lambertA m > 0` for all m by induction (base cases + positive recurrence coefficients).
Prove `derangement m > 0` for m ≥ 2 by induction.

Then apply Parts 1 and 2 with `Tendsto.add` and `Tendsto.const_mul` (for the factor of 4).

## Hard rules
- No sorry, no axiom, no native_decide
- Only kernel `decide` for small computations
- Use `norm_num` for numeric verification, `ring` for polynomial identities
- Use `field_simp` before `ring` when clearing denominators
- Line length ≤ 100 characters

## What does NOT count as done
- A sorry anywhere in the file
- Using `native_decide` (injects ofReduceBool axiom)
- A theorem whose hypothesis IS the conclusion in disguise
- Leaving the limit as a hypothesis rather than proving it

## Verification
```bash
export PATH="$HOME/.elan/bin:$PATH"
lake env lean RamanujanChallenge/Problem23.lean 2>&1
```
Must compile with no errors. Then:
```bash
lake env lean -c - <<'EOF'
import RamanujanChallenge.Problem23
#print axioms lambertRatio_tendsto
#print axioms derangementRatio_tendsto
#print axioms problem23_limit
EOF
```
Must show only {propext, Classical.choice, Quot.sound}.

## Stall protocol
If stuck on a specific lemma, deliver what compiles + precise stall report naming the exact goal state.
