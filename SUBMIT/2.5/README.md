# Ramanujan Challenge — Problem 2.5

**Catalan's constant G as a 3x3 CMF limit**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

## What is in this package

- solution.pdf / solution.tex — complete mathematical proof (10 pages)
- lean/ — Lean 4 formalization (1 sorry, see below)
- README.md — this file

## Summary

The proof identifies the CMF module as the integrated elliptic period
system Y'(k) = K(k), uses the Brafman identity to establish
G = (1/2)∫₀¹ K(k) dk, and proves convergence via Delannoy basis
decomposition with positivity and epsilon-decay.

## Lean formalization status

~3200 lines across 8 files. Main theorem: `problem25_solved : Problem25Claim`.

Fully proved (0 sorry):
- Catalan integral identity: G = ∫₀¹ (-log t)/(1+t²) dt
- Projective contraction: all 3 ratios → commonLimit (2/3 rate)
- Delannoy basis decomposition and uniqueness
- Delannoy coefficient positivity f(k) > 0 for k ≤ 10
- Error ratio bounds |g(k)/f(k) - G| < 10⁻⁵ for k ≤ 10
- Sign pattern at N=0 and N=1 (tight alternating-series bounds)

1 sorry: `positiveCatalanError_brackets` — the sign pattern persists
for all N. This corresponds to Theorem 8 in solution.pdf (the
Poincaré coefficient vanishing, verified to 30 digits + asymptotic argument).

Axioms (excluding sorry): {propext, Classical.choice, Quot.sound}.
