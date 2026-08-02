# Ramanujan Challenge — Problem 2.5

**Catalan's constant G as a 3x3 CMF limit**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

## What is in this package

- solution.pdf / solution.tex — mathematical proof (11 pages)
- lean/ — partial Lean 4 formalization (~3200 lines, 1 sorry)
- README.md — this file

## Summary

We prove that the 3x3 CMF in Problem 2.5 converges to Catalan's constant G.
The proof proceeds by:
1. Exact moment lift: the error G·Q-P has an integral representation
   via the kernel (-log t)/(1+t²).
2. Catalan integral: G = (1/2)∫₀¹ K(k) dk via the Brafman identity.
3. Delannoy basis decomposition: Q̂_N = Σ f(k)·B(N,k) with f(k) > 0.
4. Epsilon decay: h(k) = g(k)-G·f(k) has an exact integral representation
   via the inverse Delannoy transform, giving h(k) = O((-1/8)^k).
5. Splitting argument: head (polynomial) + tail (bounded by ε_k) → P̂/Q̂ → G.

## Lean formalization status

~3200 lines across 8 files. Main theorem: `problem25_solved : Problem25Claim`.

### Fully proved (0 sorry, axiom-clean)

- Catalan integral identity: G = ∫₀¹ (-log t)/(1+t²) dt (214 lines)
- Projective contraction: all 3 ratios → commonLimit at 2/3 rate (757 lines)
- Delannoy basis decomposition, uniqueness, F/G coefficients (143 lines)
- Delannoy positivity f(k) > 0 and |ε_k| < 10⁻⁵ for k ≤ 10 (159 lines)
- Sign pattern at N=0 and N=1 via tight alternating-series bounds (53+39 lines)
- Full chain: positiveCatalanError_brackets → G ∈ envelope → squeeze
  → commonLimit = G → error → 0 → Problem25Claim (553 lines)

### 1 sorry

`positiveCatalanError_brackets`: for all N, the Catalan error has both
a non-positive and a non-negative entry among the three columns.

This corresponds to Theorem 18 in solution.pdf (the epsilon-decay theorem).
The mathematical proof uses the exact integral representation of h(k)
via the inverse Delannoy transform. The Lean formalization of this step
requires:
(a) the k-recurrence for f and g (order 8, established by ore_algebra),
(b) the inverse Delannoy closed form for B⁻¹(k,N), and
(c) the integral bound on the finite sum Ψ_k.
All other steps from this sorry to Problem25Claim are fully proved.

### Axioms

All non-sorry theorems: {propext, Classical.choice, Quot.sound}.
