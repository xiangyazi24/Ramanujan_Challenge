/-
  Ramanujan Challenge Problem 2.1: a polynomial continued fraction for 6/(3−π)

  The challenge asks to prove

      a₀ + b₁/(a₁ + b₂/(a₂ + b₃/(a₃ + ⋯)))  =  6/(3−π)

  for  aₙ = −220n³ − 484n² − 301n − 42,  bₙ = 4n²(2n+1)²(5n−4)(5n+6).

  THE STRUCTURE.  Let α(n) = 220n³ − 176n² − 7n + 5 and β(n) = bₙ.  Cohen's
  encyclopedic dictionary of polynomial continued fractions (arXiv:2607.06581),
  Entry 5.3.22, records

      π = 3 + 6/(α(1) + β(1)/(α(2) + β(2)/(α(3) + ⋯)))
        = 3 + 6/(42 + 396/(1047 + 38400/(4340 + ⋯))).

  Two elementary facts connect the two:

    (1) an index shift:   aₙ = −α(n+1)   and   bₙ = β(n)   (polynomial identities);
    (2) a sign flip:      negating every partial denominator of a continued
        fraction, keeping the partial numerators, negates every convergent.

  So the challenge continued fraction is the sign-flip of Cohen's tail
  T = 6/(π−3), and its value is −T = 6/(3−π).

  Fact (2) is proved here at the level of CONVERGENTS, which avoids all
  questions about convergence of tails: we show outright that

      P̃ k = (−1)^k · P k      and      Q̃ k = −(−1)^k · Q k,

  hence P̃ k / Q̃ k = −(P k / Q k) for every k, so the two continued fractions
  converge or diverge together and their values are negatives.

  WHAT IS PROVED HERE vs. WHAT IS CITED.  (1) and (2) are proved below, and the
  main theorem is derived from them.  Cohen's Entry 5.3.22 — the value of the
  π-continued fraction itself — is the single classical input, and appears as an
  explicit hypothesis of `problem21_pcf_value`, so the dependency is visible in
  the statement rather than hidden in an axiom.

  Reference: Xiang Huang, "Solution to Ramanujan Challenge Problem 2.1", 2026.
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

open Filter Topology

noncomputable section

namespace RamanujanChallenge.P21

/-! ## Coefficients -/

/-- The challenge's partial denominators `aₙ`. -/
def a21 (n : ℕ) : ℝ := -220 * (n : ℝ) ^ 3 - 484 * (n : ℝ) ^ 2 - 301 * (n : ℝ) - 42

/-- The challenge's partial numerators `bₙ`. -/
def b21 (n : ℕ) : ℝ :=
  4 * (n : ℝ) ^ 2 * (2 * (n : ℝ) + 1) ^ 2 * (5 * (n : ℝ) - 4) * (5 * (n : ℝ) + 6)

/-- Cohen Entry 5.3.22: partial denominators `α(n)`. -/
def alphaC (n : ℕ) : ℝ := 220 * (n : ℝ) ^ 3 - 176 * (n : ℝ) ^ 2 - 7 * (n : ℝ) + 5

/-- Cohen Entry 5.3.22: partial numerators `β(n)`; identical to `b21`. -/
def betaC (n : ℕ) : ℝ :=
  4 * (n : ℝ) ^ 2 * (2 * (n : ℝ) + 1) ^ 2 * (5 * (n : ℝ) - 4) * (5 * (n : ℝ) + 6)

/-! ## The index-shift identities -/

theorem shift_a (n : ℕ) : a21 n = -alphaC (n + 1) := by
  simp only [a21, alphaC]; push_cast; ring

theorem shift_b (n : ℕ) : b21 n = betaC n := rfl

/-- Sanity values against Cohen's displayed convergents 42, 396, 1047, 38400. -/
theorem alphaC_one : alphaC 1 = 42 := by norm_num [alphaC]
theorem alphaC_two : alphaC 2 = 1047 := by norm_num [alphaC]
theorem betaC_one : betaC 1 = 396 := by norm_num [betaC]
theorem betaC_two : betaC 2 = 38400 := by norm_num [betaC]
theorem a21_zero : a21 0 = -42 := by norm_num [a21]
theorem a21_one : a21 1 = -1047 := by norm_num [a21]
theorem b21_one : b21 1 = 396 := by norm_num [b21]

/-! ## Continued-fraction convergents

For a continued fraction `c₀ + d₁/(c₁ + d₂/(c₂ + ⋯))` the classical convergent
recursions are

  P₋₁ = 1, P₀ = c₀, Pₙ = cₙ Pₙ₋₁ + dₙ Pₙ₋₂,
  Q₋₁ = 0, Q₀ = 1,  Qₙ = cₙ Qₙ₋₁ + dₙ Qₙ₋₂.

We index by `k = n + 1` so that `k : ℕ`; thus `cfP c d 0 = P₋₁` and
`cfP c d (n+1) = Pₙ`. -/

def cfP (c d : ℕ → ℝ) : ℕ → ℝ
  | 0 => 1
  | 1 => c 0
  | (n + 2) => c (n + 1) * cfP c d (n + 1) + d (n + 1) * cfP c d n

def cfQ (c d : ℕ → ℝ) : ℕ → ℝ
  | 0 => 0
  | 1 => 1
  | (n + 2) => c (n + 1) * cfQ c d (n + 1) + d (n + 1) * cfQ c d n

/-! ## The sign-flip lemma

Negating all partial denominators multiplies `Pₖ` by `(−1)^k` and `Qₖ` by
`−(−1)^k`, hence negates every convergent. -/

theorem cfP_neg (c d : ℕ → ℝ) :
    ∀ k, cfP (fun n => -c n) d k = (-1 : ℝ) ^ k * cfP c d k
  | 0 => by simp [cfP]
  | 1 => by simp [cfP]
  | (n + 2) => by
      have h1 := cfP_neg c d (n + 1)
      have h2 := cfP_neg c d n
      simp only [cfP, h1, h2]
      ring

theorem cfQ_neg (c d : ℕ → ℝ) :
    ∀ k, cfQ (fun n => -c n) d k = -((-1 : ℝ) ^ k * cfQ c d k)
  | 0 => by simp [cfQ]
  | 1 => by simp [cfQ]
  | (n + 2) => by
      have h1 := cfQ_neg c d (n + 1)
      have h2 := cfQ_neg c d n
      simp only [cfQ, h1, h2]
      ring

/-- Every convergent of the sign-flipped continued fraction is the negative of
the corresponding convergent of the original.  No convergence hypothesis is
needed: this is an identity for each `k`. -/
theorem cf_neg_convergent (c d : ℕ → ℝ) (k : ℕ) :
    cfP (fun n => -c n) d k / cfQ (fun n => -c n) d k
      = -(cfP c d k / cfQ c d k) := by
  rw [cfP_neg, cfQ_neg]
  have hs : ((-1 : ℝ) ^ k) ≠ 0 := pow_ne_zero _ (by norm_num)
  rcases eq_or_ne (cfQ c d k) 0 with hq | hq
  · simp [hq]
  · field_simp

/-! ## The two continued fractions of this problem

`cohenC`/`cohenD` present Cohen's tail
`T = α(1) + β(1)/(α(2) + β(2)/(α(3) + ⋯))`, and `a21`/`b21` present the
challenge's `a₀ + b₁/(a₁ + b₂/(a₂ + ⋯))`.  By `shift_a`, the challenge's
denominator sequence is exactly the negation of Cohen's. -/

def cohenC (n : ℕ) : ℝ := alphaC (n + 1)
def cohenD (n : ℕ) : ℝ := betaC n

theorem challenge_is_neg_of_cohen : a21 = fun n => -cohenC n := by
  funext n; exact shift_a n

/-- The challenge's convergents are the negatives of Cohen's tail convergents. -/
theorem challenge_convergent_eq (k : ℕ) :
    cfP a21 b21 k / cfQ a21 b21 k = -(cfP cohenC cohenD k / cfQ cohenC cohenD k) := by
  have hb : b21 = cohenD := by funext n; exact shift_b n
  rw [hb, challenge_is_neg_of_cohen]
  exact cf_neg_convergent cohenC cohenD k

/-! ## Main theorem

The hypothesis is exactly Cohen's Entry 5.3.22, in the form `π = 3 + 6/T`,
i.e. `T = 6/(π−3)`. -/

theorem problem21_pcf_value
    (hCohen : Tendsto (fun k => cfP cohenC cohenD k / cfQ cohenC cohenD k) atTop
      (𝓝 (6 / (Real.pi - 3)))) :
    Tendsto (fun k => cfP a21 b21 k / cfQ a21 b21 k) atTop (𝓝 (6 / (3 - Real.pi))) := by
  have hval : (6 : ℝ) / (3 - Real.pi) = -(6 / (Real.pi - 3)) := by
    rw [← neg_sub Real.pi 3, div_neg]
  rw [hval]
  exact hCohen.neg.congr fun k => (challenge_convergent_eq k).symm

end RamanujanChallenge.P21

end
