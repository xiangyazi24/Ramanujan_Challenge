# Ratio Framework: From Apéry to General Constants

**Note for future work — not part of the challenge submission.**

## The Pattern (observed across 6+ problems)

Given a constant L (like ζ(3), γ, π+e, G, ζ(2)+ζ(3), √10005/π):

1. **Two sequences, same recurrence, different initial conditions**
   - a_n, b_n satisfy the same polynomial-coefficient recurrence
   - a₀, a₁, ... (integer/rational, "denominator" sequence)
   - b₀, b₁, ... (rational, "numerator" sequence)
   - b_n/a_n → L as n → ∞

2. **Generating functions → ODE**
   - A(z) = Σ aₙ zⁿ satisfies a homogeneous ODE
   - B(z) = Σ bₙ zⁿ satisfies an INHOMOGENEOUS ODE with a constant RHS
   - The constant (e.g., 6 for Apéry ζ(3)) encodes the target value

3. **Ratio ODE (the key step)**
   - The ratio R(z) = B(z)/A(z) satisfies a polynomial ODE
   - R(z) → L as z → singularity (Poincaré-Perron)
   - The ODE can be realized as a bounded polynomial dynamical system

4. **Convergence rate from Poincaré roots**
   - The characteristic polynomial of the recurrence determines the rate
   - |b_n/a_n - L| = O(ρⁿ) where ρ = |recessive root / dominant root|

## Instances from the Ramanujan Challenge

| Problem | Target L | Recurrence order | Poincaré roots | ρ (rate) |
|---------|----------|-----------------|----------------|----------|
| 2.2 | γ | 4-term (order 3) | triple root (r-1)³ | subexponential |
| 2.3 | π+e | 5-term (order 4) | {-1±√2, 1±√2} | geometric |
| 2.6 | ζ(2)+ζ(3) | 3-term (order 2) | 1, 1/4 | polynomial (dominant!) |
| 2.7 | ζ(2)+ζ(3) | 4-term (order 3) | irred cubic F(μ) | geometric |
| 2.8 | √10005/π | 5-term (order 4) | 64R dominant | ~14 digits/term |
| Apéry | ζ(3) | 3-term (order 2) | (1±√2)⁴ | ~1.5 digits/term |

## Xiang's Apéry Adaptation (infsup.com/math/apery-adaptation-zeta3/)

The blog post establishes:
- ζ(3) ∈ R_RTCRN (computable by bounded polynomial ODE with rational ICs)
- The ratio B''/A'' satisfies a polynomial ODE
- The inhomogeneous constant 6 is load-bearing

## Generalization Plan

### Step 1: Abstract the Lean framework
Currently in Ripple: Apéry-specific (AperySequences.lean, AperyCertificate.lean,
AperyGeneratingFunction.lean, PoincaréPerron.lean).

Generalize to: `RatioConvergence.lean` — given ANY polynomial-coefficient
recurrence with two solutions a_n, b_n and Poincaré analysis showing
b_n/a_n → L, produce L as a computable real number in R_RTCRN.

### Step 2: Instantiate for each constant
- γ: Aptekarev recurrence (gauge (n!)², triple root)
- π+e: Lambert×derangement (Ore factorization L=M·P)
- G (Catalan): silver-ratio recurrence (Sym²(t²-6t+1))
- ζ(2)+ζ(3): Ore-factorizable order-2 (both solutions hypergeometric)
- √10005/π: Chudnovsky CMF (Poincaré root 64R)

### Step 3: The universal recipe
```
constant → integral representation → parametric kernel K_n
→ creative telescoping → polynomial recurrence
→ generating function ODE → ratio ODE → polynomial dynamical system
→ R_RTCRN membership
```

Each arrow is a functor in the holonomic/D-finite category.
The CMF (Conservative Matrix Field) is the matrix-valued version.

### Step 4: Connect to the hierarchy
This fits into Xiang's research program:
- GPAC ⊃ CRN ⊃ PP ⊃ NAP (the Chomsky hierarchy of analog computation)
- R_RTCRN is the set of reals computable by rate-independent CRNs
- Each Ramanujan Challenge constant proved to be in R_RTCRN via this framework
  extends the known membership of the hierarchy

## References
- Xiang's blog: infsup.com/math/apery-adaptation-zeta3/
- Ripple project: ~/repos/Ripple (8400+ lines Lean 4)
- Ramanujan Challenge: ~/repos/Ramanujan_Challenge (72 commits)
