# Impossibility Theorem: Direct GF Encoding at Regular Singular Points

**Theorem (Fable, 2026-07-13).** Let F(x) be holonomic with minimal ODE
L[F] = g(x) having a regular singular point at x=0. There is NO polynomial
PIVP dy/dt = P(y), y(0) = y₀ ∈ Qᵐ, with P(y₀) ≠ 0, such that some
polynomial projection recovers F(t).

**Proof sketch:**
1. Tangent constraint: P(y) ∝ T(t) with proportionality factor f
2. Pole absorption: F⁽ᵈ⁾ = N(x)/q_d(x) has pole → f must vanish at x=0
3. Fixed point: P(y₀) = 0, Picard-Lindelöf uniqueness → y(t) ≡ y₀

**Verified for:**
- Apéry ODE: x²(4+x)F''' + x(10+3x)F'' + (2+x)F' = 1
- Bessel J₀: t²F'' + tF' + t²F = 0

**Why ratio works:** Riccati R = B/A has R' = W/A², absorbing the 1/q_d pole.
This is algebraic desingularization at the projective bundle level.

**Constructive resolution:** Two-stage PIVP (series accumulation + regular-point
integration) is optimal. Proves: periods of Fuchsian ODEs ∈ R_RTCRN.
