# MUM Barrier: Fixed-Point Obstruction for R_RTCRN

**Core research problem — relates to Ripple OPEN_PROBLEMS.md §3.3**

## The Problem
For Apéry-type constants (ζ(3), 1/π, etc.), the generating function ODE
has x=0 as a singular point. When parameterizing x=t, the PIVP system
has all RHS = 0 at t=0 (fixed point). The trajectory doesn't start.

## Existing Work (Ripple/notes/apery-pivp.tex)
Xiang already explored: start from x₀ ∈ (0,1)∩Q with N-term partial
sum F_N(x₀) as rational IC. Error ≤ C(x₀/4)^N. But this has issues
(details TBD — Xiang says "走不通").

## Fable's Analysis (2026-07-13)

### Route 1: Rescaled Variables (one-stage PIVP)
- Set U = x² · F'' (absorbs the 1/x² singularity)
- dx/dt = 1 (NOT zero at t=0, so x=0 is NOT a fixed point of the full system)
- dU/dt = polynomial in (x, F, U, ...) — no denominators
- U(0) = 0 but dU/dt|₀ = q_d'(0)·a₁ ≠ 0 → U starts moving immediately
- Key: the FULL system (x, F, U, ...) has dx/dt = 1 ≠ 0, so t=0 is not a fixed point

### Route 2: Connection-Point with Rational Approximation (two-stage)
- Stage 1: PIVP computes F_N(x₀) (partial sum, rational)
- Stage 2: PIVP integrates from x₀ with approximate IC, flows to x=1
- Error from IC approximation propagates with bounded amplification
- No need for EXACT F(x₀), only rational approximation

### Route 3: Time Reparameterization
- Instead of dx/dt = 1, use dx/dt = x(1-x) or dx/dt = (1-x)
- Then x: 0 → 1 takes infinite time (t → ∞)
- All variables bounded, exponential convergence
- Still polynomial RHS

## ChatGPT's Analysis (Q4721)
- The Apéry ODE at x=0 is NOT full MUM: local exponents are {0, 0, 1/2}
  (not {0, 0, 0}). There's a double logarithmic sector + one x^{1/2} sector.
- The general approach: extract a rational recurrence with exponential
  error bound, avoiding the ODE singular-endpoint issue entirely.

## Xiang's Concern
F(x₀) at any rational x₀ ∈ (0,1) is an infinite series sum — NOT rational.
So it can't be an exact PIVP initial condition. The partial-sum approximation
F_N(x₀) IS rational but introduces error.

## Open Question
Does the error from F_N(x₀) ≈ F(x₀) propagate boundedly through the
ODE integration from x₀ to x=1? Xiang's apery-pivp.tex claims yes
(amplification factor ~1) but says the overall approach "走不通" — 
need to understand exactly why.

## Connection to Ramanujan Challenge
This session's 2.6 result provides a worked example: the variation-of-constants
sum is computing a period of a Fuchsian ODE, and the central-binomial
identities are the discrete encoding of that period integral.
