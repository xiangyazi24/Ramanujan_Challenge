# Task: Close P2.5 via Delannoy basis decomposition (proof.tex §5-6)

## The one sorry
`catalanError_over_denominator_tendsto_zero` in Problem25Moment.lean line 191.

## Proof route (from problems/2.5/proof.tex Theorems 6-8)

The proof does NOT use operator factorization or Brafman gauge (those fail —
nullity 0 confirmed). It works entirely within the Delannoy basis.

### Architecture

**Layer 1: Delannoy basis B(N,k)**
```
B(N,k) = 2^k · C(2k,k) · C(N,k) · C(N+k,k)
```
This is lower-triangular: B(N,k) = 0 for k > N, B(k,k) = 2^k · C(2k,k)^2 > 0.

**Layer 2: Decomposition coefficients f(k), g(k)**
By triangular inversion:
```
Q̂_N = Σ_{k=0}^N f(k) · B(N,k)     (denominator, first column)
P̂_N = Σ_{k=0}^N g(k) · B(N,k)     (numerator, first column)
```
f(k), g(k) are RATIONAL sequences computable by forward inversion.

**Layer 3: Positivity f(k) > 0**
Verified exactly for k = 0, ..., 400.
For k > 400: Poincaré asymptotics of the k-recurrence gives f(k) ~ Ak + B
with A > 0, correction O((-1/8)^k) < 10^{-300}.

**Layer 4: Error splitting**
```
P̂_N - G·Q̂_N = Σ_k f(k)·ε_k·B(N,k)
```
where ε_k = g(k)/f(k) - G.

Split at parameter K:
- Head (k ≤ K): |Σ_{k≤K}| ≤ C_K · N^{2K}, negligible vs Q̂_N ~ λ^N
- Tail (k > K): |Σ_{k>K}| ≤ max_{k>K}|ε_k| · Q̂_N (by positivity of f,B)

**Layer 5: ε_k → 0 (THE KEY STEP)**
h(k) = g(k) - G·f(k) satisfies the LCLM(L_f, L_g) operator of order 9.
Poincaré roots: (ξ-1)²(ξ+1/8)^7.

The ξ=1 modes contribute h_dominant(k) = Ak + B (polynomial growth).
The ξ=-1/8 modes contribute h_sub(k) = O((-1/8)^k) (exponential decay).

If the ξ=1 coefficients vanish: ε_k = h(k)/f(k) = O((-1/8)^k / k) → 0.

**How to prove the ξ=1 coefficients vanish:**
The coefficients A, B in h_dominant = Ak + B are determined by h(0),...,h(8).
Since h(k) = g(k) - G·f(k), these are linear in G:
  A = A_g - G·A_f
  B = B_g - G·B_f

where A_f, B_f, A_g, B_g are the ξ=1 asymptotic coefficients of f and g
(computable from f(k)/k → A_f and g(k)/k → A_g as k → ∞).

The claim A = B = 0 reduces to G = A_g/A_f = B_g/B_f.

Strategy to verify: compute A_f, A_g numerically to 100+ digits from the
k-recurrence, verify A_g/A_f = G to 100+ digits. Then the Lean proof:
1. Compute f(0),...,f(K) and g(0),...,g(K) exactly (rational)
2. Use the Catalan integral identity (already proved) to bound G
3. Bound |ε_k| for k ≤ K using rational G bounds
4. For k > K, use the k-recurrence Poincaré asymptotics

## What to build in Lean

### File 1: Problem25Delannoy.lean
- Define B(N,k) = 2^k · C(2k,k) · C(N,k) · C(N+k,k)
- Prove B(N,k) > 0, B(k,k) > 0, B(N,k) = 0 for k > N
- Define f(k) and g(k) by triangular inversion from Q̂ and P̂
- Prove f(k) > 0 for k = 0,...,K (norm_num, with set_option maxRecDepth/maxHeartbeats)
- Prove the decomposition Q̂_N = Σ f(k)·B(N,k) and P̂_N = Σ g(k)·B(N,k)

### File 2: Problem25EpsilonDecay.lean
- State ε_k = g(k)/f(k) - catalanConstant
- Prove |ε_k| ≤ bound_k for k ≤ K (using catalan_tight_lower/upper)
- State the k-recurrence and its Poincaré polynomial (ξ-1)²(ξ+1/8)^7
- Prove the asymptotic bound for k > K

### File 3: Problem25Close.lean (or modify Problem25Moment.lean)
- The splitting argument: head + tail bound
- Close catalanError_over_denominator_tendsto_zero
- Close problem25_solved : Problem25Claim

## Existing infrastructure
- Problem25.lean: Q̂, P̂ sequences, catalanConstant, all recurrence machinery
- Problem25Integral.lean: catalanConstant = ∫(-log t)/(1+t²)dt (0 sorry)
- Problem25TightBounds.lean: tight rational bounds on G (0 sorry)
- Problem25Connection.lean: projective contraction, commonLimit (0 sorry)

## Key Python scripts (for reference values)
- scripts/p25_epsilon_rate.py: computes ε_k and verifies rate → -1/8
- scripts/p25_solution_decomposition.py: computes f(k), g(k), decomposes in Delannoy basis

## Build
```bash
~/.elan/toolchains/leanprover--lean4---v4.29.0/bin/lake env lean RamanujanChallenge/Problem25Moment.lean
```
Target: 0 sorry.

## CRITICAL: what makes this a REAL proof (not numerical verification)

The paper proof says "verified to 40+ digits." That's not a proof.
The REAL proof must show the ξ=1 coefficient is EXACTLY zero.

Approach: the ξ=1 coefficient of h(k) = g(k) - G·f(k) is
  c₁ = Σ_{j=0}^{8} α_j · (g(j) - G·f(j))
where α_j are Casorati inversion coefficients (rational, from the k-recurrence).

This equals c₁ = (Σ α_j·g(j)) - G·(Σ α_j·f(j)) = B₁ - G·A₁.

c₁ = 0 iff G = B₁/A₁. The Lean proof must:
1. Compute A₁ and B₁ exactly (rational, from the k-recurrence)
2. Show catalanConstant = B₁/A₁

Step 2 can use catalanConstant_eq_integral + the integral evaluation
∫(-log t)/(1+t²)dt = B₁/A₁.

This is the genuine mathematical content: an identity between the
Catalan integral and a specific ratio arising from the k-recurrence.
