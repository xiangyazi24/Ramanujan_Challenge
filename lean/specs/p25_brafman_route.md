# Task: Close P2.5 via Brafman identity + elliptic K integral

## Goal
Prove `catalanError_over_denominator_tendsto_zero` in Problem25Moment.lean
by establishing the Brafman integral representation, which gives subdominance
structurally (not numerically).

## The proof route (from problems/2.5/proof.tex §4-6)

### Step 1: Brafman identity for Delannoy squares
The central Delannoy numbers D_n = P_n(3) = 1, 3, 13, 63, 321, ...
satisfy (n+1)D_{n+1} = 3(2n+1)D_n - nD_{n-1}.

Brafman's theorem (1951): The generating function of D_n² is
  Σ_{n≥0} D_n² z^n = (2/(π(1-z))) · K(4√(2z)/(1-z))
for |z| < ρ = 17-12√2, where K is the complete elliptic integral of the first kind:
  K(k) = ∫_0^{π/2} dθ/√(1 - k²sin²θ) = (π/2) · ₂F₁(1/2, 1/2; 1; k²)

This follows from Clausen's identity ₂F₁(a,b;a+b+1/2;z)² = ₃F₂(2a,2b,a+b;2a+2b,a+b+1/2;z)
applied to the Legendre function P_n(x) = ₂F₁(-n,n+1;1;(1-x)/2) at x = 3.

### Step 2: Catalan integral via Brafman
Define the substitution k(z) = 4√(2z)/(1-z), mapping [0,ρ] → [0,1].
Then:
  G = (π√2/2) ∫_0^ρ ((1+z)/(√z·(1-z))) · F(z) dz
    = (1/2) ∫_0^1 K(k) dk

where the second equality uses k'(z) = 2√2(1+z)/(√z·(1-z)²) and the Brafman identity.

The final step: (1/2)∫_0^1 K(k) dk = ∫_0^1 arctan(t)/t dt = G
uses Fubini and the integral representation of arctan.

### Step 3: Beukers-type error integral
From the Brafman representation, the Delannoy basis expansion gives:
  Q̂_N = Σ_k f(k)·B(N,k)  where B(N,k) = 2^k·C(2k,k)·C(N,k)·C(N+k,k)

The error has the integral representation:
  G·Q̂_N - P̂_N = ∫_0^1 (x²(-log x)/(1+x²)) · P_N^Leg(√(1-8x²))² dx

where P_N^Leg is the Legendre polynomial. The integrand is bounded by
C · n · σ^N with σ = (15+4√14)/(17+12√2) < 1, giving exponential decay.

However, the CMF in Problem 2.5 uses a DIFFERENT, faster construction with
rate ρ ≈ 0.0294 (vs σ ≈ 0.882). The connection between the CMF and the
Delannoy model goes through the elliptic K module structure.

### Step 4: Convergence from positivity + rate bound
Once we have:
  |G·Q̂_N - P̂_N| ≤ C · rate^N · |Q̂_N|
with rate < 1, dividing by Q̂_N gives the desired convergence to 0.

## Ripple infrastructure available

The repo ~/repos/Ripple has (Lean v4.29):
- `ellipticK` defined (Number/Hypergeometric/PeriodBridge.lean:112)
- Clausen's identity (Number/Hypergeometric/Clausen.lean)
- ₂F₁ hypergeometric functions and their properties
- Period bridge infrastructure connecting ₂F₁ to ellipticK

## What needs to be built

1. **Delannoy numbers in Lean**: D_n = Σ_k C(n,k)²·2^k. Verify they satisfy
   the order-2 recurrence. Show D_n² satisfies the order-3 recurrence
   (shared Poincaré roots with Problem 2.5 CMF).

2. **Brafman identity**: Σ D_n² z^n = (2/(π(1-z)))·K(4√(2z)/(1-z)).
   This needs Clausen's identity from Ripple + the specific parameter values
   for Legendre polynomials at x=3.

3. **Catalan integral from Brafman**: G = (1/2)∫₀¹ K(k) dk.
   Uses the substitution k(z) and integration.

4. **Error bound**: The integral representation gives
   |G·Q_N - P_N| = O(ρ^N · Q_N) with ρ < 1.

## Build & verify
```bash
~/.elan/toolchains/leanprover--lean4---v4.29.0/bin/lake env lean RamanujanChallenge/Problem25Moment.lean
```
Target: 0 sorry.

## Key references
- Brafman (1951): "Generating functions of Jacobi and related polynomials"
- proof.tex §4-6 in problems/2.5/
- ~/repos/Ripple/Ripple/Number/Hypergeometric/PeriodBridge.lean
- ~/repos/Ripple/Ripple/Number/Hypergeometric/Clausen.lean

## Strategy notes
- Start by reading Ripple's Clausen.lean and PeriodBridge.lean to understand
  what's available
- The Brafman identity is the KEY bridge — it connects the combinatorial
  Delannoy world to the analytic elliptic integral world
- The error bound follows once the integral representation is established
- Don't try to factor the recurrence operator — it's irreducible (proved in
  proof.tex Theorem 4). The proof goes AROUND the operator, through integrals.

## Files to create/modify
- May create new files in RamanujanChallenge/ for the Brafman/Delannoy layer
- Target: replace the sorry in Problem25Moment.lean line 191
- May import from ~/repos/Ripple if needed (add as lake dependency or copy)
