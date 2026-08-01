# DOCTRINE — Ramanujan Challenge Lean Formalization (2026-08-01)

## Goal
Close the Lean formalization gaps for P2.5, P2.6, P2.7 — in that order.
Each problem has deep infrastructure (0 sorry); each gap is one analytic identity.

## Toolchain
Lean v4.29.0 + Mathlib v4.29.0. Local build (uisai2 down).
Lake: `/Users/huangx/.elan/toolchains/leanprover--lean4---v4.29.0/bin/lake`

## P2.5 — Catalan's constant G
Gap: `commonLimit = catalanConstant` (Problem25Connection.lean:748).

**(a) Moment formula route** (paper §2):
  - Moment identities DONE (Problem25Moment.lean, 0 sorry).
  - Next: Catalan integral G = -∫₀¹ log(t)/(1+t²) dt.
  - Then: remainder integral G·Q-P = ∫[-log(t)/(1+t²)]·R(t²) dt.
  - Then: decay → commonLimit = G.

**(b) Torsion + numerical** (P3.1-style rational reconstruction).

**(c) Direct Brafman** (needs K(k) — heavy).

## P2.6 — ζ(2) + ζ(3)
Gap: nested weight-3 inverse-binomial evaluation.
**(a)** Existing Cyclotomic/Generating/WeightThree/Nested files (~4700 lines).

## P2.7 — ζ(2) + ζ(3) (4-term Zudilin)
Gap: source normalization = Barnes contour integral.
**(a)** Barnes integral representation via residue calculus.

## Status
- P3.1: DONE (SUBMIT, 0 sorry)
- P2.8: DONE (SUBMIT, 0 sorry)
- P2.5: moment formulas proved, avenue (a) active
- P2.6: ~6200 lines, gap isolated
- P2.7: ~2800 lines, Barnes integral → 0 proved
