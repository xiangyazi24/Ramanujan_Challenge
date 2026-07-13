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

## Concrete Computation (2026-07-13): dx/dt=1 DOES NOT WORK

With x=t (linear clock), the Apéry system has:
- F(0)=0, G₁(0)=F'(0)=1/2, G₂(0)=F''(0)=-1/12
- dG₂/dt = [1-(2+t)G₁-t(10+3t)G₂] / [t²(4+t)]
- Numerator at t=0: 1-2·(1/2)-0 = 0 (one-fold zero)
- Denominator at t=0: 4t² (two-fold zero)
- So dG₂/dt ~ 1/(8t) → DIVERGES

Frobenius deflation U=t·G₂: U starts moving (dU/dt|₀=1/24≠0),
but recovering G₂=U/t reintroduces 1/t in dG₁/dt = G₂ = U/t.

ROOT CAUSE: F'(0) = 1/2 ≠ 0 (forced by the inhomogeneous term "1").
If all ICs were zero, deflation could work recursively. But the
non-zero Frobenius coefficient at the non-zero indicial root makes
the singularity irreducible by polynomial state-variable substitution.

## Open Question (refined)
Is the regular singular point with F'(0)≠0 an ABSOLUTE barrier for
first-order polynomial PIVP? Or is there a non-obvious variable
substitution (beyond simple Frobenius deflation) that resolves it?

## Connection to Ramanujan Challenge
This session's 2.6 result provides a worked example: the variation-of-constants
sum is computing a period of a Fuchsian ODE, and the central-binomial
identities are the discrete encoding of that period integral.

## BREAKTHROUGH (ChatGPT Q4724): Apéry ODE is COMPLETELY REDUCIBLE

The double root at ρ=0 is NOT an irreducible MUM block — it's just θ² = (constant)·(log).
The homogeneous operator factors into THREE first-order factors:
  L_hom = (D + p₃/q₃)·(D + 1/x)·D

This means:
- Differentiating once (Y=F') removes the double root: exponents become {-1, -1/2}
- The analytic inhomogeneous branch has a first-order hypergeometric recurrence
- aₙ = (-1)^{n-1}/(n³·C(2n,n)) with ratio < 1/4 (geometric convergence)

**Constructive PIVP route:** compile the rational partial-sum recurrence directly.
No need to fight the singular point at all.

**The general lesson:** "double indicial root" ≠ "irreducible". Check Ore factorizability
BEFORE declaring MUM. For ζ(3), the operator is completely reducible.
The TRUE MUM barrier arises only for irreducible blocks (e.g., some Calabi-Yau operators).

## FABLE DEEP ANALYSIS: Recurrence = Fixed Point, GF = Unstable Manifold

### Minimal 2D system (Ore factorization L = M₁·D₁·D)
```
dx/dt = x(4+x)
dJ/dt = 1 - (2+x)J     where J = xF'' + F'
```
Fixed point: (0, 1/2). Jacobian eigenvalues: 4 (unstable), -2 (stable).
Clean saddle — no degeneracy.

### The beautiful reframing
- The RECURRENCE for aₙ IS the fixed-point condition (RHS vanishes order by order at x=0)
- The GENERATING FUNCTION F(x) IS the unstable manifold of the saddle
- ε-perturbation IS the natural parameterization along the unstable manifold
- ICs at x=ε: J(ε) = 1/2 - ε/12 + ε²/24 - ... (all RATIONAL coefficients)

### Eigenvalue analysis for 4D system (x, J, I, F)
Eigenvalues: 4, -2, 0, 0 (saddle + neutral center)
Unstable manifold is 1-dimensional, parameterized by x.

### Why the fixed point is inescapable (4 independent proofs)
(a) Picard-Lindelöf uniqueness at Lipschitz fixed point
(b) Tangent proportionality → proportionality factor absorbs pole → vanishes
(c) Taylor coefficients of RHS = recurrence relations (vanish by design)
(d) Every blow-up produces another 0/0 (infinite regress for resonant indices)

### Paper-ready theorem
Every holonomic GF at a regular singular point:
(i) No PIVP from x=0 (fixed point = recurrence)
(ii) Ore factorization → minimal saddle system, GF = unstable manifold
(iii) ε-perturbation gives rational ICs to arbitrary precision
(iv) Apéry ζ(3): 2-species CRN suffices

## XIANG'S PRECISION OBJECTION (2026-07-13)

The ε-perturbation with N-term truncation gives J(ε) rational with O(ε^N) error.
But RTCRN requires the trajectory LIMIT to be EXACTLY the target constant.
Truncation error propagates through the ODE → the limit is NOT exactly ζ(3).

This is NOT an engineering issue — it's a theoretical precision requirement.
A single fixed PIVP must have its trajectory converge to ζ(3) EXACTLY.

Current status: generating function direct encoding at a regular singular
point REMAINS OPEN in the exact RTCRN sense. The ratio and integral routes
work because they avoid the singular point entirely.

## POISSON TRANSFORM ROUTE (ChatGPT Q4725) — potential resolution

Z*(t) = e^{-t} · Σₙ Sₙ · tⁿ/n!  where Sₙ = partial sums of aₖ

Properties:
- No truncation, no ε, contains the FULL infinite sequence
- Z*(t) → ζ(3) as t→∞, error ≤ C·e^{-(1-q)t} (EXACT limit)
- Z* is D-finite → differentially algebraic → has GPAC/PIVP representation
- All Taylor coefficients at t=0 are RATIONAL (= Sₙ/n!)
- With t=r², the ODE for Z* becomes polynomial

Open question: can we write the EXPLICIT small-dimensional polynomial
system, and is it non-singular at t=0?

Shannon-Pour-El gives existence but not a printed vector field.
The companion ODE from the D-finite structure may still have
regular singular point at t=0 (same issue recycled).

Xiang's verdict: pending.
