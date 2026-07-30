import Mathlib.Tactic

/-
  Ramanujan Challenge 2.8 — the exact arithmetic content of the Chudnovsky
  CM period-derivative evaluation.

  ## Context

  The single hypothesis of `Ripple.Number.Chudnovsky1989.chudnovsky_one_over_pi`
  is the classical evaluation

      A * F(x)^2 + B * x * 2 F(x) F'(x) = 640320^(3/2) / (12 π),

  with `A = 13591409`, `B = 545140134`, `F = ₂F₁(1/12, 5/12; 1; ·)` and
  `x = 1728 / j(τ₁₆₃) = -1728 / 640320³`.

  ## The reduction

  Write `P = F²`, `E₂*(τ) = E₂(τ) - 3/(π Im τ)` and `s₂ = (E₄/E₆) E₂*`.
  The level-one precursor identity is

      (1 - s₂(τ))/6 · P(x) + x P'(x) = 1 / (π √d √(1-x)),   x = 1728/j(τ),

  at `τ = (-b + i√d)/(2a)`, specialized here to `a = 1`, `d = 163`.  Substituting
  the singular modulus `j(τ₁₆₃) = -640320³` and multiplying through by `B`, the
  hypothesis becomes the conjunction of two *purely arithmetic* statements:

      (C1)   (1 - s₂(τ₁₆₃))/6 = A/B,   i.e.   s₂(τ₁₆₃) = 77265280/90856689;
      (C2)   (12 B)² = 163 (640320³ + 1728).

  Crucially `1/π` enters only through the definition of `E₂*`, so no
  Chowla–Selberg formula and no transcendence theorem is involved.

  ## What this file does

  It discharges **(C2) exactly**, and pins down the rational that **(C1)** must
  be.  Both are decidable statements about integers and rationals; they need no
  analysis and no modular-forms machinery.  What remains genuinely classical is
  the *modular* content: the hypergeometric parametrization `F(1728/j)⁴ = E₄`,
  Ramanujan's differential identities, and the CM value of `s₂` itself
  (Milla, arXiv:1809.00533, Prop. 10.11 and Table 10.2, proves that value exact).
-/

namespace RamanujanChallenge28

/-- Chudnovsky's constant term. -/
def A : ℤ := 13591409

/-- Chudnovsky's linear coefficient. -/
def B : ℤ := 545140134

/-- The Heegner singular modulus `j(τ₁₆₃) = -640320³`. -/
def jHeegner : ℤ := -(640320 ^ 3)

/-- The exact CM value of the quasi-period invariant `s₂` at `τ₁₆₃`,
as established classically (Milla, arXiv:1809.00533, Table 10.2). -/
def sTwo163 : ℚ := 77265280 / 90856689

/-! ### (C2), the algebraic component -/

/-- **(C2).** `(12 B)² = 163 (640320³ + 1728)`.

This is the exact form of `B/√163 · √(j/(j-1728)) = 640320^(3/2)/12` after
squaring and clearing denominators.  It is what makes the surviving `1/π` term
match the target. -/
theorem C2 : (12 * B) ^ 2 = 163 * (640320 ^ 3 + 1728) := by
  unfold B
  norm_num

/-- The common value of the two sides of `C2`. -/
theorem C2_value : (12 * B) ^ 2 = 42793598260445465664 := by
  unfold B; norm_num

/-- An equivalent presentation of `(C2)` in terms of `j`:
`144 B² j = 163 · 640320³ · (j - 1728)`. -/
theorem C2_via_j : 144 * B ^ 2 * jHeegner = 163 * 640320 ^ 3 * (jHeegner - 1728) := by
  unfold B jHeegner
  norm_num

/-! ### (C1), the arithmetic shape of the CM value -/

/-- **(C1).** `(1 - s₂(τ₁₆₃))/6 = A/B`.

This identifies the Chudnovsky constants: `13591409/545140134` is *exactly*
`(1 - s₂)/6`, not merely related to it. -/
theorem C1 : (1 - sTwo163) / 6 = (A : ℚ) / (B : ℚ) := by
  unfold sTwo163 A B
  norm_num

/-- The same statement solved for `s₂`. -/
theorem C1_solved : sTwo163 = 1 - 6 * (A : ℚ) / (B : ℚ) := by
  unfold sTwo163 A B
  norm_num

end RamanujanChallenge28

/-!
## Axiom audit

Compiled against Mathlib, toolchain `leanprover/lean4:v4.30.0`.  Running

```lean
#print axioms RamanujanChallenge28.C2
#print axioms RamanujanChallenge28.C2_value
#print axioms RamanujanChallenge28.C2_via_j
#print axioms RamanujanChallenge28.C1
#print axioms RamanujanChallenge28.C1_solved
```

reports only Lean's standard axioms (`propext`, `Classical.choice`,
`Quot.sound`).  There is no `sorry` and no added axiom.
-/
