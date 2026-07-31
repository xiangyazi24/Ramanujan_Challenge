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

  The target identity in Cohen's Entry 5.3.22 is not assumed below.  Its
  continued fraction is proved to be the exact even contraction of the earlier
  period-two Apéry fraction

      3 + 1/(9 - 36/(18 + 18/(38 - 400/(55 + ⋯)))).

  This includes the equivalence multipliers, identities for both continuants,
  positivity/nonvanishing, and the limit transfer through `x ↦ 6/(x-3)`.
  The cofinal odd diagonal is then proved to converge to `π` by an elementary
  moment integral.  A rational derivative certificate gives its three-term
  recurrence, two explicit antiderivatives give the initial values, and a
  geometric integrand bound makes the moment error tend to zero.  Thus
  `problem21_pcf_value` has no external hypothesis.

  Reference: Xiang Huang, "Solution to Ramanujan Challenge Problem 2.1", 2026.
-/
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
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

/-! ## The period-two Apéry fraction behind Cohen's entry

Cohen's entry is not used as a black box below.  It is the even contraction
of the following period-two continued fraction for `π`:

    3 + 1/(9 - 36/(18 + 18/(38 - 400/(55 + 150/(87 - ⋯))))).

The next definitions use `n / 2` only to package the two polynomial branches
as ordinary coefficient sequences. -/

def aperyC21 (n : ℕ) : ℝ :=
  if n = 0 then 3
  else if n % 2 = 0 then
    10 * ((n / 2 : ℕ) : ℝ) ^ 2 + 7 * ((n / 2 : ℕ) : ℝ) + 1
  else
    10 * ((n / 2 : ℕ) : ℝ) ^ 2 + 19 * ((n / 2 : ℕ) : ℝ) + 9

def aperyD21 (n : ℕ) : ℝ :=
  if n = 1 then 1
  else if n % 2 = 0 then
    -4 * ((n / 2 : ℕ) : ℝ) ^ 2 * (2 * ((n / 2 : ℕ) : ℝ) + 1) ^ 2
  else
    ((n / 2 : ℕ) : ℝ) * (((n / 2 : ℕ) : ℝ) + 1) *
      (2 * ((n / 2 : ℕ) : ℝ) + 1) ^ 2

@[simp] theorem aperyC21_zero : aperyC21 0 = 3 := by simp [aperyC21]

@[simp] theorem aperyC21_odd (n : ℕ) :
    aperyC21 (2 * n + 1) = 10 * (n : ℝ) ^ 2 + 19 * (n : ℝ) + 9 := by
  have hdiv : (2 * n + 1) / 2 = n := by omega
  have hmod : (2 * n + 1) % 2 = 1 := by omega
  simp [aperyC21, hdiv, hmod]

@[simp] theorem aperyC21_even_succ (n : ℕ) :
    aperyC21 (2 * n + 2) =
      10 * ((n : ℝ) + 1) ^ 2 + 7 * ((n : ℝ) + 1) + 1 := by
  have hdiv : (2 * n + 2) / 2 = n + 1 := by omega
  have hmod : (2 * n + 2) % 2 = 0 := by omega
  simp [aperyC21, hdiv, hmod]

@[simp] theorem aperyC21_odd_succ (n : ℕ) :
    aperyC21 (2 * n + 3) =
      10 * ((n : ℝ) + 1) ^ 2 + 19 * ((n : ℝ) + 1) + 9 := by
  rw [show 2 * n + 3 = 2 * (n + 1) + 1 by omega, aperyC21_odd]
  push_cast
  ring

@[simp] theorem aperyC21_even_succ_succ (n : ℕ) :
    aperyC21 (2 * n + 4) =
      10 * ((n : ℝ) + 2) ^ 2 + 7 * ((n : ℝ) + 2) + 1 := by
  rw [show 2 * n + 4 = 2 * (n + 1) + 2 by omega, aperyC21_even_succ]
  push_cast
  ring

@[simp] theorem aperyD21_one : aperyD21 1 = 1 := by simp [aperyD21]

@[simp] theorem aperyD21_odd_succ (n : ℕ) :
    aperyD21 (2 * n + 3) =
      ((n : ℝ) + 1) * ((n : ℝ) + 2) * (2 * (n : ℝ) + 3) ^ 2 := by
  have hdiv : (2 * n + 3) / 2 = n + 1 := by omega
  have hmod : (2 * n + 3) % 2 = 1 := by omega
  simp [aperyD21, hdiv, hmod]
  ring

@[simp] theorem aperyD21_even_succ (n : ℕ) :
    aperyD21 (2 * n + 2) =
      -4 * ((n : ℝ) + 1) ^ 2 * (2 * (n : ℝ) + 3) ^ 2 := by
  have hdiv : (2 * n + 2) / 2 = n + 1 := by omega
  have hmod : (2 * n + 2) % 2 = 0 := by omega
  have hne : 2 * n + 2 ≠ 1 := by omega
  simp [aperyD21, hdiv, hmod, hne]
  left
  ring

@[simp] theorem aperyD21_even_succ_succ (n : ℕ) :
    aperyD21 (2 * n + 4) =
      -4 * ((n : ℝ) + 2) ^ 2 * (2 * (n : ℝ) + 5) ^ 2 := by
  rw [show 2 * n + 4 = 2 * (n + 1) + 2 by omega, aperyD21_even_succ]
  push_cast
  ring

/-! The even contraction is equivalent to the continued fraction obtained by
adjoining `3 + 6/(...)` to Cohen's tail. -/

def wholeC21 : ℕ → ℝ
  | 0 => 3
  | n + 1 => cohenC n

def wholeD21 : ℕ → ℝ
  | 0 => 0
  | 1 => 6
  | n + 2 => cohenD (n + 1)

/-- Equivalence multipliers for the even contraction. -/
def contractScale21 : ℕ → ℝ
  | 0 => 1
  | n + 1 =>
      contractScale21 n * (5 * (n : ℝ) + 1) /
        (((n : ℝ) + 1) * (2 * (n : ℝ) + 3))

theorem contractScale21_pos : ∀ n, 0 < contractScale21 n
  | 0 => by simp [contractScale21]
  | n + 1 => by
      have hn := contractScale21_pos n
      rw [contractScale21]
      positivity

private theorem contraction_c_coefficient (n : ℕ) :
    (10 * ((n : ℝ) + 1) ^ 2 + 7 * ((n : ℝ) + 1) + 1) *
        alphaC (n + 2) * contractScale21 (n + 1) =
      contractScale21 (n + 2) *
        ((10 * ((n : ℝ) + 1) ^ 2 + 7 * ((n : ℝ) + 1) + 1) *
            ((10 * ((n : ℝ) + 2) ^ 2 + 7 * ((n : ℝ) + 2) + 1) *
                (10 * ((n : ℝ) + 1) ^ 2 + 19 * ((n : ℝ) + 1) + 9) -
              4 * ((n : ℝ) + 2) ^ 2 * (2 * (n : ℝ) + 5) ^ 2) +
          (10 * ((n : ℝ) + 2) ^ 2 + 7 * ((n : ℝ) + 2) + 1) *
            (((n : ℝ) + 1) * ((n : ℝ) + 2) * (2 * (n : ℝ) + 3) ^ 2)) := by
  simp only [contractScale21, alphaC]
  push_cast
  field_simp
  ring

private theorem contraction_d_coefficient (n : ℕ) :
    (10 * ((n : ℝ) + 1) ^ 2 + 7 * ((n : ℝ) + 1) + 1) *
        betaC (n + 1) * contractScale21 n =
      -contractScale21 (n + 2) *
        (10 * ((n : ℝ) + 2) ^ 2 + 7 * ((n : ℝ) + 2) + 1) *
        (((n : ℝ) + 1) * ((n : ℝ) + 2) * (2 * (n : ℝ) + 3) ^ 2) *
        (-4 * ((n : ℝ) + 1) ^ 2 * (2 * (n : ℝ) + 3) ^ 2) := by
  simp only [contractScale21, betaC]
  push_cast
  field_simp
  ring

private theorem apery_scaled_even_recurrence
    (u : ℕ → ℝ)
    (hu : ∀ m, u (m + 2) = aperyC21 (m + 1) * u (m + 1) + aperyD21 (m + 1) * u m)
    (n : ℕ) :
    contractScale21 (n + 2) * u (2 * n + 5) =
      alphaC (n + 2) * contractScale21 (n + 1) * u (2 * n + 3) +
        betaC (n + 1) * contractScale21 n * u (2 * n + 1) := by
  have hrec := hu (2 * n + 1)
  have hodd := hu (2 * n + 2)
  have heven := hu (2 * n + 3)
  have hrec' :
      u (2 * n + 3) =
        (10 * ((n : ℝ) + 1) ^ 2 + 7 * ((n : ℝ) + 1) + 1) * u (2 * n + 2) +
          (-4 * ((n : ℝ) + 1) ^ 2 * (2 * (n : ℝ) + 3) ^ 2) * u (2 * n + 1) := by
    rw [show (2 * n + 1) + 2 = 2 * n + 3 by omega,
      show (2 * n + 1) + 1 = 2 * n + 2 by omega] at hrec
    simpa only [aperyC21_even_succ, aperyD21_even_succ] using hrec
  have hodd' :
      u (2 * n + 4) =
        (10 * ((n : ℝ) + 1) ^ 2 + 19 * ((n : ℝ) + 1) + 9) * u (2 * n + 3) +
          (((n : ℝ) + 1) * ((n : ℝ) + 2) * (2 * (n : ℝ) + 3) ^ 2) *
            u (2 * n + 2) := by
    rw [show (2 * n + 2) + 2 = 2 * n + 4 by omega,
      show (2 * n + 2) + 1 = 2 * n + 3 by omega] at hodd
    simpa only [aperyC21_odd_succ, aperyD21_odd_succ] using hodd
  have heven' :
      u (2 * n + 5) =
        (10 * ((n : ℝ) + 2) ^ 2 + 7 * ((n : ℝ) + 2) + 1) * u (2 * n + 4) +
          (-4 * ((n : ℝ) + 2) ^ 2 * (2 * (n : ℝ) + 5) ^ 2) * u (2 * n + 3) := by
    rw [show (2 * n + 3) + 2 = 2 * n + 5 by omega,
      show (2 * n + 3) + 1 = 2 * n + 4 by omega] at heven
    simpa only [aperyC21_even_succ_succ, aperyD21_even_succ_succ] using heven
  have hA :
      (10 * ((n : ℝ) + 1) ^ 2 + 7 * ((n : ℝ) + 1) + 1) ≠ 0 := by
    positivity
  have hy :
      u (2 * n + 2) =
        (u (2 * n + 3) -
            (-4 * ((n : ℝ) + 1) ^ 2 * (2 * (n : ℝ) + 3) ^ 2) * u (2 * n + 1)) /
          (10 * ((n : ℝ) + 1) ^ 2 + 7 * ((n : ℝ) + 1) + 1) := by
    apply (eq_div_iff hA).2
    linarith [hrec']
  have hc := contraction_c_coefficient n
  have hd := contraction_d_coefficient n
  rw [heven', hodd', hy]
  field_simp [hA] at *
  linear_combination -hc * u (2 * n + 3) - hd * u (2 * n + 1)

theorem wholeP21_eq_scaled_even : ∀ n,
    cfP wholeC21 wholeD21 (n + 1) =
      contractScale21 n * cfP aperyC21 aperyD21 (2 * n + 1)
  | 0 => by norm_num [cfP, wholeC21, contractScale21, aperyC21]
  | 1 => by
      norm_num [cfP, wholeC21, wholeD21, cohenC, alphaC, cohenD, betaC,
        contractScale21, aperyC21, aperyD21]
  | n + 2 => by
      have ih1 := wholeP21_eq_scaled_even (n + 1)
      have ih0 := wholeP21_eq_scaled_even n
      rw [show (n + 2) + 1 = (n + 1) + 2 by omega, cfP, ih1, ih0]
      have hrec := apery_scaled_even_recurrence
        (u := cfP aperyC21 aperyD21) (fun _ => rfl) n
      convert hrec.symm using 1
      all_goals simp only [wholeC21, wholeD21, cohenC, cohenD]
      all_goals ring_nf

theorem wholeQ21_eq_scaled_even : ∀ n,
    cfQ wholeC21 wholeD21 (n + 1) =
      contractScale21 n * cfQ aperyC21 aperyD21 (2 * n + 1)
  | 0 => by norm_num [cfQ, wholeC21, contractScale21]
  | 1 => by
      norm_num [cfQ, wholeC21, wholeD21, cohenC, alphaC, cohenD, betaC,
        contractScale21, aperyC21, aperyD21]
  | n + 2 => by
      have ih1 := wholeQ21_eq_scaled_even (n + 1)
      have ih0 := wholeQ21_eq_scaled_even n
      rw [show (n + 2) + 1 = (n + 1) + 2 by omega, cfQ, ih1, ih0]
      have hrec := apery_scaled_even_recurrence
        (u := cfQ aperyC21 aperyD21) (fun _ => rfl) n
      convert hrec.symm using 1
      all_goals simp only [wholeC21, wholeD21, cohenC, cohenD]
      all_goals ring_nf

/-- The convergents of `wholeC21`/`wholeD21` are exactly the even convergents
of the period-two Apéry fraction. -/
theorem whole_convergent_eq_apery_even (n : ℕ) :
    cfP wholeC21 wholeD21 (n + 1) / cfQ wholeC21 wholeD21 (n + 1) =
      cfP aperyC21 aperyD21 (2 * n + 1) / cfQ aperyC21 aperyD21 (2 * n + 1) := by
  rw [wholeP21_eq_scaled_even, wholeQ21_eq_scaled_even]
  have hs := ne_of_gt (contractScale21_pos n)
  rcases eq_or_ne (cfQ aperyC21 aperyD21 (2 * n + 1)) 0 with hq | hq
  · simp [hq]
  · field_simp

/-! ## The integral seed for the odd diagonal

Put

`Jₙ = 2 ∫₀¹ x^(2n+2) (1-x²)^(2n+1) / (1+x²)^(2n+2) dx`.

Thus `Jₙ` is the elementary hyperbolic moment usually denoted `I_(2n+2)`
after the substitution `x = tanh (t/2)`.  The first two values are
`J₀ = π - 3` and `J₁ = 22 - 7π`.  A rational certificate proves the
three-term recurrence whose two independent rational coefficient sequences
are exactly the odd continuants above.  Finally the integrand is bounded by a
geometric sequence, so `Jₙ → 0`. -/

private def momentBase21 (x : ℝ) : ℝ :=
  2 * x ^ 2 * (1 - x ^ 2) / (1 + x ^ 2) ^ 2

private def momentRatio21 (x : ℝ) : ℝ :=
  x ^ 2 * (1 - x ^ 2) ^ 2 / (1 + x ^ 2) ^ 2

private def momentIntegrand21 (n : ℕ) (x : ℝ) : ℝ :=
  momentBase21 x * momentRatio21 x ^ n

private def moment21 (n : ℕ) : ℝ :=
  ∫ x in (0 : ℝ)..1, momentIntegrand21 n x

private def certificatePoly21 (m u : ℝ) : ℝ :=
  (-10 * m ^ 2 - 12 * m) +
    (30 * m ^ 2 + 51 * m + 14) * u +
    (60 * m ^ 2 + 82 * m + 18) * u ^ 2 +
    (-30 * m ^ 2 - 41 * m - 6) * u ^ 3 +
    (-10 * m ^ 2 - 12 * m - 2) * u ^ 4

private def certificate21 (m x : ℝ) : ℝ :=
  x * (1 - x ^ 2) * certificatePoly21 m (x ^ 2) / (1 + x ^ 2) ^ 3

private def certificateDeriv21 (m x : ℝ) : ℝ :=
  (50 * m ^ 2 * x ^ 12 + 170 * m ^ 2 * x ^ 10 + 90 * m ^ 2 * x ^ 8 -
      660 * m ^ 2 * x ^ 6 + 30 * m ^ 2 * x ^ 4 + 170 * m ^ 2 * x ^ 2 -
      10 * m ^ 2 +
    60 * m * x ^ 12 + 219 * m * x ^ 10 + 138 * m * x ^ 8 -
      892 * m * x ^ 6 - 34 * m * x ^ 4 + 249 * m * x ^ 2 - 12 * m +
    10 * x ^ 12 + 34 * x ^ 10 + 12 * x ^ 8 - 172 * x ^ 6 -
      22 * x ^ 4 + 42 * x ^ 2) / (1 + x ^ 2) ^ 4

private theorem hasDerivAt_certificatePoly21 (m u : ℝ) :
    HasDerivAt (certificatePoly21 m)
      ((30 * m ^ 2 + 51 * m + 14) +
        2 * (60 * m ^ 2 + 82 * m + 18) * u +
        3 * (-30 * m ^ 2 - 41 * m - 6) * u ^ 2 +
        4 * (-10 * m ^ 2 - 12 * m - 2) * u ^ 3) u := by
  change HasDerivAt
    (fun v =>
      (-10 * m ^ 2 - 12 * m) +
        (30 * m ^ 2 + 51 * m + 14) * v +
        (60 * m ^ 2 + 82 * m + 18) * v ^ 2 +
        (-30 * m ^ 2 - 41 * m - 6) * v ^ 3 +
        (-10 * m ^ 2 - 12 * m - 2) * v ^ 4)
    _ u
  have h0 := hasDerivAt_const u (-10 * m ^ 2 - 12 * m)
  have h1 :=
    (hasDerivAt_const u (30 * m ^ 2 + 51 * m + 14)).mul
      (hasDerivAt_id u)
  have h2 :=
    (hasDerivAt_const u (60 * m ^ 2 + 82 * m + 18)).mul
      ((hasDerivAt_id u).pow 2)
  have h3 :=
    (hasDerivAt_const u (-30 * m ^ 2 - 41 * m - 6)).mul
      ((hasDerivAt_id u).pow 3)
  have h4 :=
    (hasDerivAt_const u (-10 * m ^ 2 - 12 * m - 2)).mul
      ((hasDerivAt_id u).pow 4)
  convert ((((h0.add h1).add h2).add h3).add h4) using 1 <;>
    simp only [id_eq] <;> ring

private theorem moment_certificate_derivative (n : ℕ) (x : ℝ) :
    HasDerivAt
      (fun y => momentIntegrand21 n y * certificate21 ((n : ℝ) + 1) y)
      (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * (5 * (n : ℝ) + 6) *
          momentIntegrand21 (n + 2) x +
        alphaC (n + 2) * momentIntegrand21 (n + 1) x -
        2 * ((n : ℝ) + 1) * (2 * (n : ℝ) + 3) * (5 * (n : ℝ) + 11) *
          momentIntegrand21 n x) x := by
  have hden : 1 + x ^ 2 ≠ 0 := by positivity
  have hx : HasDerivAt (fun y : ℝ => y) 1 x := hasDerivAt_id x
  have hsq : HasDerivAt (fun y : ℝ => y ^ 2) (2 * x) x := by
    convert hx.pow 2 using 1 <;> ring
  have hsub : HasDerivAt (fun y : ℝ => 1 - y ^ 2) (-2 * x) x := by
    convert (hasDerivAt_const x 1).sub hsq using 1 <;> ring
  have hadd : HasDerivAt (fun y : ℝ => 1 + y ^ 2) (2 * x) x := by
    convert (hasDerivAt_const x 1).add hsq using 1 <;> ring
  have hbase :=
    (((hasDerivAt_const x (2 : ℝ)).mul hsq).mul hsub).div
      (hadd.pow 2) (pow_ne_zero 2 hden)
  have hratio :=
    (hsq.mul (hsub.pow 2)).div (hadd.pow 2) (pow_ne_zero 2 hden)
  have hpoly :=
    (hasDerivAt_certificatePoly21 ((n : ℝ) + 1) (x ^ 2)).comp x hsq
  have hcert :=
    ((hx.mul hsub).mul hpoly).div (hadd.pow 3) (pow_ne_zero 3 hden)
  have hbase' : HasDerivAt momentBase21
      (-4 * x * (3 * x ^ 2 - 1) / (1 + x ^ 2) ^ 3) x := by
    convert hbase using 1
    simp
    field_simp [hden]
    ring
  have hratio' : HasDerivAt momentRatio21
      (2 * x * (x - 1) * (x + 1) * (x ^ 4 + 4 * x ^ 2 - 1) /
        (1 + x ^ 2) ^ 3) x := by
    convert hratio using 1
    simp
    field_simp [hden]
    ring
  have hcert' : HasDerivAt (certificate21 ((n : ℝ) + 1))
      (certificateDeriv21 ((n : ℝ) + 1) x) x := by
    convert hcert using 1
    simp [certificatePoly21, certificateDeriv21]
    field_simp [hden]
    ring
  have htotal := (hbase'.mul (hratio'.pow n)).mul hcert'
  convert htotal using 1
  cases n with
  | zero =>
      norm_num [momentIntegrand21, momentBase21, momentRatio21, certificate21,
        certificatePoly21, certificateDeriv21, alphaC]
      field_simp [hden]
      ring
  | succ k =>
      simp [momentIntegrand21, pow_succ, alphaC, momentBase21, momentRatio21,
        certificate21, certificatePoly21, certificateDeriv21]
      field_simp [hden]
      ring

private theorem continuous_momentBase21 : Continuous momentBase21 := by
  unfold momentBase21
  exact
    (((continuous_const.mul (continuous_id.pow 2)).mul
      (continuous_const.sub (continuous_id.pow 2))).div
      ((continuous_const.add (continuous_id.pow 2)).pow 2)
      (fun x => by positivity))

private theorem continuous_momentRatio21 : Continuous momentRatio21 := by
  unfold momentRatio21
  exact
    (((continuous_id.pow 2).mul
      ((continuous_const.sub (continuous_id.pow 2)).pow 2)).div
      ((continuous_const.add (continuous_id.pow 2)).pow 2)
      (fun x => by positivity))

private theorem continuous_momentIntegrand21 (n : ℕ) :
    Continuous (momentIntegrand21 n) := by
  exact continuous_momentBase21.mul (continuous_momentRatio21.pow n)

private theorem moment21_recurrence (n : ℕ) :
    2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * (5 * (n : ℝ) + 6) *
          moment21 (n + 2) +
        alphaC (n + 2) * moment21 (n + 1) -
        2 * ((n : ℝ) + 1) * (2 * (n : ℝ) + 3) * (5 * (n : ℝ) + 11) *
          moment21 n = 0 := by
  let A : ℝ :=
    2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * (5 * (n : ℝ) + 6)
  let B : ℝ := alphaC (n + 2)
  let C : ℝ :=
    2 * ((n : ℝ) + 1) * (2 * (n : ℝ) + 3) * (5 * (n : ℝ) + 11)
  have hA : IntervalIntegrable
      (fun x => A * momentIntegrand21 (n + 2) x) MeasureTheory.volume 0 1 :=
    ((continuous_momentIntegrand21 (n + 2)).intervalIntegrable 0 1).const_mul A
  have hB : IntervalIntegrable
      (fun x => B * momentIntegrand21 (n + 1) x) MeasureTheory.volume 0 1 :=
    ((continuous_momentIntegrand21 (n + 1)).intervalIntegrable 0 1).const_mul B
  have hC : IntervalIntegrable
      (fun x => C * momentIntegrand21 n x) MeasureTheory.volume 0 1 :=
    ((continuous_momentIntegrand21 n).intervalIntegrable 0 1).const_mul C
  have hAB := hA.add hB
  have hABC := hAB.sub hC
  have hftc :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt
      (f := fun x =>
        momentIntegrand21 n x * certificate21 ((n : ℝ) + 1) x)
      (f' := fun x =>
        A * momentIntegrand21 (n + 2) x +
          B * momentIntegrand21 (n + 1) x -
          C * momentIntegrand21 n x)
      (fun x _ => by
        simpa only [A, B, C] using moment_certificate_derivative n x)
      hABC
  rw [intervalIntegral.integral_sub hAB hC,
    intervalIntegral.integral_add hA hB,
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul] at hftc
  simpa [A, B, C, moment21, momentIntegrand21, momentBase21, certificate21] using hftc

private def momentPrimitiveZero21 (x : ℝ) : ℝ :=
  -2 * x - 2 * x / (1 + x ^ 2) + 4 * Real.arctan x

private theorem hasDerivAt_momentPrimitiveZero21 (x : ℝ) :
    HasDerivAt momentPrimitiveZero21 (momentIntegrand21 0 x) x := by
  have hden : 1 + x ^ 2 ≠ 0 := by positivity
  have hx : HasDerivAt (fun y : ℝ => y) 1 x := hasDerivAt_id x
  have hsq : HasDerivAt (fun y : ℝ => y ^ 2) (2 * x) x := by
    convert hx.pow 2 using 1 <;> ring
  have hadd : HasDerivAt (fun y : ℝ => 1 + y ^ 2) (2 * x) x := by
    convert (hasDerivAt_const x 1).add hsq using 1 <;> ring
  have hlin := (hasDerivAt_const x (-2 : ℝ)).mul hx
  have hfrac :=
    ((hasDerivAt_const x (-2 : ℝ)).mul hx).div hadd hden
  have hatan :=
    (hasDerivAt_const x (4 : ℝ)).mul (Real.hasDerivAt_arctan x)
  convert (hlin.add hfrac).add hatan using 1
  · funext y
    simp [momentPrimitiveZero21]
    ring
  · simp [momentIntegrand21, momentBase21]
    field_simp [hden]
    ring

private theorem moment21_zero : moment21 0 = Real.pi - 3 := by
  rw [moment21]
  calc
    (∫ x in (0 : ℝ)..1, momentIntegrand21 0 x) =
        momentPrimitiveZero21 1 - momentPrimitiveZero21 0 := by
      exact intervalIntegral.integral_eq_sub_of_hasDerivAt
        (fun x _ => hasDerivAt_momentPrimitiveZero21 x)
        ((continuous_momentIntegrand21 0).intervalIntegrable 0 1)
    _ = Real.pi - 3 := by
      simp [momentPrimitiveZero21, Real.arctan_one]
      ring

private def momentPrimitiveOne21 (x : ℝ) : ℝ :=
  -2 * x ^ 3 / 3 + 14 * x +
      (66 * x ^ 5 + 100 * x ^ 3 + 42 * x) / (3 * (1 + x ^ 2) ^ 3) -
    28 * Real.arctan x

private theorem hasDerivAt_momentPrimitiveOne21 (x : ℝ) :
    HasDerivAt momentPrimitiveOne21 (momentIntegrand21 1 x) x := by
  have hden : 1 + x ^ 2 ≠ 0 := by positivity
  have hx : HasDerivAt (fun y : ℝ => y) 1 x := hasDerivAt_id x
  have hsq : HasDerivAt (fun y : ℝ => y ^ 2) (2 * x) x := by
    convert hx.pow 2 using 1 <;> ring
  have hadd : HasDerivAt (fun y : ℝ => 1 + y ^ 2) (2 * x) x := by
    convert (hasDerivAt_const x 1).add hsq using 1 <;> ring
  have hcubic := (hasDerivAt_const x (-2 / 3 : ℝ)).mul (hx.pow 3)
  have hlinear := (hasDerivAt_const x (14 : ℝ)).mul hx
  have hnum :=
    (((hasDerivAt_const x (66 : ℝ)).mul (hx.pow 5)).add
      ((hasDerivAt_const x (100 : ℝ)).mul (hx.pow 3))).add
      ((hasDerivAt_const x (42 : ℝ)).mul hx)
  have hdenom :=
    (hasDerivAt_const x (3 : ℝ)).mul (hadd.pow 3)
  have hdenom_ne : (3 : ℝ) * (1 + x ^ 2) ^ 3 ≠ 0 :=
    mul_ne_zero (by norm_num) (pow_ne_zero 3 hden)
  have hquot := hnum.div hdenom hdenom_ne
  have hatan :=
    (hasDerivAt_const x (28 : ℝ)).mul (Real.hasDerivAt_arctan x)
  convert ((hcubic.add hlinear).add hquot).sub hatan using 1
  · funext y
    simp [momentPrimitiveOne21]
    ring
  · simp [momentIntegrand21, momentBase21, momentRatio21]
    field_simp [hden]
    ring

private theorem moment21_one : moment21 1 = 22 - 7 * Real.pi := by
  rw [moment21]
  calc
    (∫ x in (0 : ℝ)..1, momentIntegrand21 1 x) =
        momentPrimitiveOne21 1 - momentPrimitiveOne21 0 := by
      exact intervalIntegral.integral_eq_sub_of_hasDerivAt
        (fun x _ => hasDerivAt_momentPrimitiveOne21 x)
        ((continuous_momentIntegrand21 1).intervalIntegrable 0 1)
    _ = 22 - 7 * Real.pi := by
      simp [momentPrimitiveOne21, Real.arctan_one]
      ring

private theorem moment_pointwise_bounds {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    0 ≤ momentBase21 x ∧ momentBase21 x ≤ 2 ∧
      0 ≤ momentRatio21 x ∧ momentRatio21 x ≤ (1 : ℝ) / 4 := by
  have hx2 : 0 ≤ x ^ 2 := sq_nonneg x
  have hx2le : x ^ 2 ≤ 1 := by
    nlinarith [mul_nonneg hx0 (sub_nonneg.mpr hx1)]
  have hsub : 0 ≤ 1 - x ^ 2 := sub_nonneg.mpr hx2le
  have hsuble : 1 - x ^ 2 ≤ 1 := by linarith
  have hden : 0 < 1 + x ^ 2 := by positivity
  have hprod : x ^ 2 * (1 - x ^ 2) ≤ 1 := by
    calc
      x ^ 2 * (1 - x ^ 2) ≤ 1 * (1 - x ^ 2) :=
        mul_le_mul_of_nonneg_right hx2le hsub
      _ ≤ 1 := by linarith
  have hbase0 : 0 ≤ momentBase21 x := by
    simp only [momentBase21]
    positivity
  have hbase2 : momentBase21 x ≤ 2 := by
    rw [momentBase21, div_le_iff₀ (pow_pos hden 2)]
    have hden1 : 1 ≤ (1 + x ^ 2) ^ 2 := by nlinarith
    nlinarith
  have hratio0 : 0 ≤ momentRatio21 x := by
    simp only [momentRatio21]
    positivity
  have hsubsq : (1 - x ^ 2) ^ 2 ≤ 1 := by
    nlinarith [mul_nonneg hsub (sub_nonneg.mpr hsuble)]
  have hprod2 : x ^ 2 * (1 - x ^ 2) ^ 2 ≤ x ^ 2 :=
    mul_le_of_le_one_right hx2 hsubsq
  have hamgm : 4 * x ^ 2 ≤ (1 + x ^ 2) ^ 2 := by
    nlinarith [sq_nonneg (1 - x ^ 2)]
  have hratio4 : momentRatio21 x ≤ (1 : ℝ) / 4 := by
    rw [momentRatio21, div_le_iff₀ (pow_pos hden 2)]
    nlinarith
  exact ⟨hbase0, hbase2, hratio0, hratio4⟩

private theorem moment21_norm_bound (n : ℕ) :
    ‖moment21 n‖ ≤ 2 * ((1 : ℝ) / 4) ^ n := by
  have hpoint :
      ∀ x ∈ Set.uIoc (0 : ℝ) 1, ‖momentIntegrand21 n x‖ ≤
        2 * ((1 : ℝ) / 4) ^ n := by
    intro x hx
    rw [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1), Set.mem_Ioc] at hx
    obtain ⟨hb0, hb2, hr0, hr4⟩ :=
      moment_pointwise_bounds hx.1.le hx.2
    rw [Real.norm_eq_abs, abs_of_nonneg]
    · exact mul_le_mul hb2 (pow_le_pow_left₀ hr0 hr4 n)
        (pow_nonneg hr0 n) (by norm_num)
    · exact mul_nonneg hb0 (pow_nonneg hr0 n)
  simpa [moment21] using
    (intervalIntegral.norm_integral_le_of_norm_le_const hpoint)

private theorem moment21_tendsto_zero :
    Tendsto moment21 atTop (𝓝 0) := by
  have hpow :
      Tendsto (fun n : ℕ => ((1 : ℝ) / 4) ^ n) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)
  have hgeom :
      Tendsto (fun n : ℕ => 2 * ((1 : ℝ) / 4) ^ n) atTop (𝓝 0) := by
    simpa using tendsto_const_nhds.mul hpow
  rw [Metric.tendsto_atTop] at hgeom ⊢
  intro ε hε
  obtain ⟨N, hN⟩ := hgeom ε hε
  refine ⟨N, fun n hn => ?_⟩
  rw [dist_zero_right]
  exact lt_of_le_of_lt (moment21_norm_bound n) (by
    simpa [dist_zero_right] using hN n hn)

theorem wholeP21_eq_tail : ∀ k,
    cfP wholeC21 wholeD21 (k + 1) =
      3 * cfP cohenC cohenD k + 6 * cfQ cohenC cohenD k
  | 0 => by norm_num [cfP, cfQ, wholeC21]
  | 1 => by
      norm_num [cfP, cfQ, wholeC21, wholeD21, cohenC, alphaC, cohenD, betaC]
  | n + 2 => by
      have ih1 := wholeP21_eq_tail (n + 1)
      have ih0 := wholeP21_eq_tail n
      rw [show (n + 2) + 1 = (n + 1) + 2 by omega, cfP, ih1, ih0]
      simp only [wholeC21, wholeD21, cohenC, cohenD, cfQ, cfP]
      ring

theorem wholeQ21_eq_tail : ∀ k,
    cfQ wholeC21 wholeD21 (k + 1) = cfP cohenC cohenD k
  | 0 => by norm_num [cfP, cfQ, wholeC21]
  | 1 => by
      norm_num [cfP, cfQ, wholeC21, wholeD21, cohenC, alphaC]
  | n + 2 => by
      have ih1 := wholeQ21_eq_tail (n + 1)
      have ih0 := wholeQ21_eq_tail n
      rw [show (n + 2) + 1 = (n + 1) + 2 by omega, cfQ, ih1, ih0]
      simp only [wholeC21, wholeD21, cohenC, cohenD, cfP]

private theorem cohenC_pos (n : ℕ) : 0 < cohenC n := by
  have hpoly :
      cohenC n =
        220 * (n : ℝ) ^ 3 + 484 * (n : ℝ) ^ 2 + 301 * (n : ℝ) + 42 := by
    simp only [cohenC, alphaC]
    push_cast
    ring
  rw [hpoly]
  positivity

private theorem cohenD_nonneg (n : ℕ) : 0 ≤ cohenD n := by
  rcases n with _ | n
  · norm_num [cohenD, betaC]
  · have hn : (1 : ℝ) ≤ (n + 1 : ℕ) := by norm_num
    have hfac : 0 ≤ 5 * ((n + 1 : ℕ) : ℝ) - 4 := by linarith
    simp only [cohenD, betaC]
    positivity

private theorem cohenD_pos {n : ℕ} (hn : 0 < n) : 0 < cohenD n := by
  have hnr : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hfac : 0 < 5 * (n : ℝ) - 4 := by linarith
  simp only [cohenD, betaC]
  positivity

private theorem cohenP_pos : ∀ k, 0 < cfP cohenC cohenD k
  | 0 => by simp [cfP]
  | 1 => by simpa [cfP] using cohenC_pos 0
  | n + 2 => by
      rw [cfP]
      have hc := cohenC_pos (n + 1)
      have hd := cohenD_pos (Nat.succ_pos n)
      have hp1 := cohenP_pos (n + 1)
      have hp0 := cohenP_pos n
      positivity

private theorem cohenQ_nonneg : ∀ k, 0 ≤ cfQ cohenC cohenD k
  | 0 => by simp [cfQ]
  | 1 => by simp [cfQ]
  | n + 2 => by
      rw [cfQ]
      exact add_nonneg
        (mul_nonneg (cohenC_pos (n + 1)).le (cohenQ_nonneg (n + 1)))
        (mul_nonneg (cohenD_nonneg (n + 1)) (cohenQ_nonneg n))

private theorem cohenQ_pos_succ : ∀ k, 0 < cfQ cohenC cohenD (k + 1)
  | 0 => by simp [cfQ]
  | n + 1 => by
      rw [show (n + 1) + 1 = n + 2 by omega, cfQ]
      exact add_pos_of_pos_of_nonneg
        (mul_pos (cohenC_pos (n + 1)) (cohenQ_pos_succ n))
        (mul_nonneg (cohenD_nonneg (n + 1)) (cohenQ_nonneg n))

/-! ## Identification of the moments with the continuant error -/

private def scaleDen21 (n : ℕ) : ℝ :=
  2 * ((n : ℝ) + 1) * (2 * (n : ℝ) + 3) * (5 * (n : ℝ) + 1)

private def momentScale21 : ℕ → ℝ
  | 0 => 1
  | n + 1 => -momentScale21 n / scaleDen21 n

private theorem scaleDen21_pos (n : ℕ) : 0 < scaleDen21 n := by
  simp only [scaleDen21]
  positivity

private theorem momentScale21_A (n : ℕ) :
    2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * (5 * (n : ℝ) + 6) *
        momentScale21 (n + 2) =
      -momentScale21 (n + 1) := by
  simp [momentScale21, scaleDen21]
  field_simp
  ring

private theorem momentScale21_beta (n : ℕ) :
    betaC (n + 1) * momentScale21 (n + 1) =
      -(2 * ((n : ℝ) + 1) * (2 * (n : ℝ) + 3) *
          (5 * (n : ℝ) + 11) * momentScale21 n) := by
  simp [momentScale21, scaleDen21, betaC]
  field_simp
  ring

private def wholeLinearError21 (n : ℕ) : ℝ :=
  Real.pi * cfQ wholeC21 wholeD21 (n + 1) -
    cfP wholeC21 wholeD21 (n + 1)

private theorem wholeLinearError21_recurrence (n : ℕ) :
    wholeLinearError21 (n + 2) =
      alphaC (n + 2) * wholeLinearError21 (n + 1) +
        betaC (n + 1) * wholeLinearError21 n := by
  unfold wholeLinearError21
  rw [show n + 2 + 1 = (n + 1) + 2 by omega, cfQ, cfP]
  simp only [wholeC21, wholeD21, cohenC, cohenD]
  ring

private def wholeError21 (n : ℕ) : ℝ :=
  momentScale21 n * wholeLinearError21 n

private theorem wholeError21_recurrence (n : ℕ) :
    2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * (5 * (n : ℝ) + 6) *
          wholeError21 (n + 2) +
        alphaC (n + 2) * wholeError21 (n + 1) -
        2 * ((n : ℝ) + 1) * (2 * (n : ℝ) + 3) * (5 * (n : ℝ) + 11) *
          wholeError21 n = 0 := by
  let A : ℝ :=
    2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * (5 * (n : ℝ) + 6)
  let C : ℝ :=
    2 * ((n : ℝ) + 1) * (2 * (n : ℝ) + 3) * (5 * (n : ℝ) + 11)
  have hsA : A * momentScale21 (n + 2) = -momentScale21 (n + 1) := by
    simpa only [A] using momentScale21_A n
  have hsB :
      betaC (n + 1) * momentScale21 (n + 1) =
        -(C * momentScale21 n) := by
    simpa only [C, mul_assoc] using momentScale21_beta n
  have hc1 : A * momentScale21 (n + 2) + momentScale21 (n + 1) = 0 := by
    rw [hsA]
    ring
  have hc0 :
      A * momentScale21 (n + 2) * betaC (n + 1) -
          C * momentScale21 n = 0 := by
    rw [hsA]
    calc
      -momentScale21 (n + 1) * betaC (n + 1) -
          C * momentScale21 n =
          -(betaC (n + 1) * momentScale21 (n + 1)) -
            C * momentScale21 n := by ring
      _ = 0 := by rw [hsB]; ring
  rw [wholeError21, wholeLinearError21_recurrence]
  change A * (momentScale21 (n + 2) *
      (alphaC (n + 2) * wholeLinearError21 (n + 1) +
        betaC (n + 1) * wholeLinearError21 n)) +
    alphaC (n + 2) *
      (momentScale21 (n + 1) * wholeLinearError21 (n + 1)) -
    C * (momentScale21 n * wholeLinearError21 n) = 0
  calc
    _ =
        (A * momentScale21 (n + 2) + momentScale21 (n + 1)) *
            alphaC (n + 2) * wholeLinearError21 (n + 1) +
          (A * momentScale21 (n + 2) * betaC (n + 1) -
            C * momentScale21 n) * wholeLinearError21 n := by ring
    _ = 0 := by rw [hc1, hc0]; ring

private theorem moment21_eq_wholeError21 : ∀ n, moment21 n = wholeError21 n
  | 0 => by
      rw [moment21_zero]
      norm_num [wholeError21, wholeLinearError21, momentScale21, cfP, cfQ,
        wholeC21]
  | 1 => by
      rw [moment21_one]
      norm_num [wholeError21, wholeLinearError21, momentScale21, scaleDen21,
        cfP, cfQ, wholeC21, wholeD21, cohenC, alphaC]
      ring
  | n + 2 => by
      have ih0 := moment21_eq_wholeError21 n
      have ih1 := moment21_eq_wholeError21 (n + 1)
      have hm := moment21_recurrence n
      have he := wholeError21_recurrence n
      rw [ih0, ih1] at hm
      let A : ℝ :=
        2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * (5 * (n : ℝ) + 6)
      have hA : A ≠ 0 := by
        dsimp [A]
        positivity
      apply (mul_left_cancel₀ hA)
      change A * moment21 (n + 2) = A * wholeError21 (n + 2)
      dsimp [A] at hm he ⊢
      linarith

private theorem wholeQ21_pos (n : ℕ) :
    0 < cfQ wholeC21 wholeD21 (n + 1) := by
  rw [wholeQ21_eq_tail]
  exact cohenP_pos n

private theorem wholeQ21_nonneg : ∀ n, 0 ≤ cfQ wholeC21 wholeD21 n
  | 0 => by simp [cfQ]
  | n + 1 => (wholeQ21_pos n).le

private theorem wholeD21_nonneg : ∀ n, 0 ≤ wholeD21 n
  | 0 => by simp [wholeD21]
  | 1 => by simp [wholeD21]
  | n + 2 => by
      simpa only [wholeD21] using cohenD_nonneg (n + 1)

private theorem scaleDen21_le_alpha (n : ℕ) :
    scaleDen21 n ≤ alphaC (n + 1) := by
  have hn : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
  have hfac :
      0 ≤ (5 * (n : ℝ) + 4) *
        (40 * (n : ℝ) ^ 2 + 54 * (n : ℝ) + 9) := by
    positivity
  simp only [scaleDen21, alphaC]
  push_cast
  nlinarith

private theorem wholeQ21_step_lower (n : ℕ) :
    scaleDen21 n * cfQ wholeC21 wholeD21 (n + 1) ≤
      cfQ wholeC21 wholeD21 (n + 2) := by
  have hq := (wholeQ21_pos n).le
  have hα := scaleDen21_le_alpha n
  have hdq :
      0 ≤ wholeD21 (n + 1) * cfQ wholeC21 wholeD21 n :=
    mul_nonneg (wholeD21_nonneg (n + 1)) (wholeQ21_nonneg n)
  calc
    scaleDen21 n * cfQ wholeC21 wholeD21 (n + 1) ≤
        alphaC (n + 1) * cfQ wholeC21 wholeD21 (n + 1) :=
      mul_le_mul_of_nonneg_right hα hq
    _ ≤ alphaC (n + 1) * cfQ wholeC21 wholeD21 (n + 1) +
        wholeD21 (n + 1) * cfQ wholeC21 wholeD21 n :=
      le_add_of_nonneg_right hdq
    _ = cfQ wholeC21 wholeD21 (n + 2) := by
      rw [cfQ]
      simp only [wholeC21, cohenC]

private theorem scaledDen21_mono (n : ℕ) :
    |momentScale21 n| * cfQ wholeC21 wholeD21 (n + 1) ≤
      |momentScale21 (n + 1)| * cfQ wholeC21 wholeD21 (n + 2) := by
  have hD := scaleDen21_pos n
  have hstep := wholeQ21_step_lower n
  rw [momentScale21, abs_div, abs_neg, abs_of_pos hD]
  calc
    |momentScale21 n| * cfQ wholeC21 wholeD21 (n + 1) =
        (|momentScale21 n| / scaleDen21 n) *
          (scaleDen21 n * cfQ wholeC21 wholeD21 (n + 1)) := by
      field_simp [ne_of_gt hD]
    _ ≤ (|momentScale21 n| / scaleDen21 n) *
          cfQ wholeC21 wholeD21 (n + 2) :=
      mul_le_mul_of_nonneg_left hstep (div_nonneg (abs_nonneg _) hD.le)

private theorem scaledDen21_ge_one : ∀ n,
    1 ≤ |momentScale21 n| * cfQ wholeC21 wholeD21 (n + 1)
  | 0 => by norm_num [momentScale21, cfQ]
  | n + 1 => (scaledDen21_ge_one n).trans (scaledDen21_mono n)

private theorem whole21_convergent_error_bound (n : ℕ) :
    |Real.pi -
        cfP wholeC21 wholeD21 (n + 1) /
          cfQ wholeC21 wholeD21 (n + 1)| ≤
      |moment21 n| := by
  let p := cfP wholeC21 wholeD21 (n + 1)
  let q := cfQ wholeC21 wholeD21 (n + 1)
  have hq : 0 < q := by
    simpa only [q] using wholeQ21_pos n
  have hlinear : Real.pi * q - p = q * (Real.pi - p / q) := by
    field_simp [ne_of_gt hq]
  have herr : moment21 n =
      momentScale21 n * (Real.pi * q - p) := by
    rw [moment21_eq_wholeError21]
    rfl
  have habs :
      |moment21 n| =
        (|momentScale21 n| * q) * |Real.pi - p / q| := by
    rw [herr, hlinear, abs_mul, abs_mul, abs_of_pos hq]
    ring
  rw [habs]
  calc
    |Real.pi - p / q| = 1 * |Real.pi - p / q| := by ring
    _ ≤ (|momentScale21 n| * q) * |Real.pi - p / q| :=
      mul_le_mul_of_nonneg_right
        (by simpa only [q] using scaledDen21_ge_one n)
        (abs_nonneg _)

private theorem whole21_shift_converges :
    Tendsto
      (fun n => cfP wholeC21 wholeD21 (n + 1) /
        cfQ wholeC21 wholeD21 (n + 1))
      atTop (𝓝 Real.pi) := by
  have hm := moment21_tendsto_zero
  rw [Metric.tendsto_atTop] at hm ⊢
  intro ε hε
  obtain ⟨N, hN⟩ := hm ε hε
  refine ⟨N, fun n hn => ?_⟩
  have hj : |moment21 n| < ε := by
    simpa [Real.dist_eq] using hN n hn
  have he := lt_of_le_of_lt (whole21_convergent_error_bound n) hj
  simpa [Real.dist_eq, abs_sub_comm] using he

/-- The whole contracted continued fraction converges to `π`, with no
external convergence hypothesis. -/
theorem whole21_converges :
    Tendsto (fun k => cfP wholeC21 wholeD21 k / cfQ wholeC21 wholeD21 k)
      atTop (𝓝 Real.pi) :=
  (tendsto_add_atTop_iff_nat 1).mp whole21_shift_converges

/-- The odd diagonal of the period-two Apéry fraction converges to `π`. -/
theorem apery21_odd_diagonal_converges :
    Tendsto
      (fun n => cfP aperyC21 aperyD21 (2 * n + 1) /
        cfQ aperyC21 aperyD21 (2 * n + 1))
      atTop (𝓝 Real.pi) := by
  exact whole21_shift_converges.congr fun n =>
    whole_convergent_eq_apery_even n

theorem cohen_convergent_eq_whole_mobius (k : ℕ) :
    cfP cohenC cohenD (k + 1) / cfQ cohenC cohenD (k + 1) =
      6 / (cfP wholeC21 wholeD21 (k + 2) / cfQ wholeC21 wholeD21 (k + 2) - 3) := by
  rw [wholeP21_eq_tail (k + 1), wholeQ21_eq_tail (k + 1)]
  have hp := ne_of_gt (cohenP_pos (k + 1))
  have hq := ne_of_gt (cohenQ_pos_succ k)
  field_simp
  ring

/-! The following conditional transfer lemmas are retained as reusable
interfaces.  The integral argument above has already proved their diagonal
hypothesis unconditionally. -/

theorem whole21_converges_of_apery_diagonal
    (hAperyDiagonal : Tendsto
      (fun n => cfP aperyC21 aperyD21 (2 * n + 1) /
        cfQ aperyC21 aperyD21 (2 * n + 1))
      atTop (𝓝 Real.pi)) :
    Tendsto (fun k => cfP wholeC21 wholeD21 k / cfQ wholeC21 wholeD21 k)
      atTop (𝓝 Real.pi) := by
  have hshift :
      Tendsto
        (fun n => cfP wholeC21 wholeD21 (n + 1) / cfQ wholeC21 wholeD21 (n + 1))
        atTop (𝓝 Real.pi) :=
    hAperyDiagonal.congr fun n => (whole_convergent_eq_apery_even n).symm
  exact (tendsto_add_atTop_iff_nat 1).mp hshift

theorem whole21_converges_of_apery
    (hApery : Tendsto
      (fun k => cfP aperyC21 aperyD21 k / cfQ aperyC21 aperyD21 k)
      atTop (𝓝 Real.pi)) :
    Tendsto (fun k => cfP wholeC21 wholeD21 k / cfQ wholeC21 wholeD21 k)
      atTop (𝓝 Real.pi) := by
  have hindex : Tendsto (fun n : ℕ => 2 * n + 1) atTop atTop := by
    rw [tendsto_atTop]
    intro b
    exact eventually_atTop.2 ⟨b, fun n hn => by omega⟩
  exact whole21_converges_of_apery_diagonal (hApery.comp hindex)

/-- Cohen's Entry 5.3.22 follows algebraically from the period-two Apéry
continued fraction.  In fact only its cofinal odd-indexed diagonal is needed,
rather than convergence of the full period-two sequence. -/
theorem cohen_converges_of_apery_diagonal
    (hAperyDiagonal : Tendsto
      (fun n => cfP aperyC21 aperyD21 (2 * n + 1) /
        cfQ aperyC21 aperyD21 (2 * n + 1))
      atTop (𝓝 Real.pi)) :
    Tendsto (fun k => cfP cohenC cohenD k / cfQ cohenC cohenD k)
      atTop (𝓝 (6 / (Real.pi - 3))) := by
  have hwhole := whole21_converges_of_apery_diagonal hAperyDiagonal
  have hwholeShift :
      Tendsto
        (fun k => cfP wholeC21 wholeD21 (k + 2) / cfQ wholeC21 wholeD21 (k + 2))
        atTop (𝓝 Real.pi) :=
    hwhole.comp (tendsto_add_atTop_nat 2)
  have hpi : Real.pi - 3 ≠ 0 := ne_of_gt (sub_pos.mpr Real.pi_gt_three)
  have hmobius :
      Tendsto
        (fun k => 6 /
          (cfP wholeC21 wholeD21 (k + 2) / cfQ wholeC21 wholeD21 (k + 2) - 3))
        atTop (𝓝 (6 / (Real.pi - 3))) :=
    tendsto_const_nhds.div (hwholeShift.sub_const 3) hpi
  have hshift :
      Tendsto
        (fun k => cfP cohenC cohenD (k + 1) / cfQ cohenC cohenD (k + 1))
        atTop (𝓝 (6 / (Real.pi - 3))) :=
    hmobius.congr fun k => (cohen_convergent_eq_whole_mobius k).symm
  exact (tendsto_add_atTop_iff_nat 1).mp hshift

theorem cohen_converges_of_apery
    (hApery : Tendsto
      (fun k => cfP aperyC21 aperyD21 k / cfQ aperyC21 aperyD21 k)
      atTop (𝓝 Real.pi)) :
    Tendsto (fun k => cfP cohenC cohenD k / cfQ cohenC cohenD k)
      atTop (𝓝 (6 / (Real.pi - 3))) := by
  have hindex : Tendsto (fun n : ℕ => 2 * n + 1) atTop atTop := by
    rw [tendsto_atTop]
    intro b
    exact eventually_atTop.2 ⟨b, fun n hn => by omega⟩
  exact cohen_converges_of_apery_diagonal (hApery.comp hindex)

theorem cohen21_converges :
    Tendsto (fun k => cfP cohenC cohenD k / cfQ cohenC cohenD k)
      atTop (𝓝 (6 / (Real.pi - 3))) :=
  cohen_converges_of_apery_diagonal apery21_odd_diagonal_converges

/-! ## Main theorem

The first theorem records the old one-line reduction from Cohen's entry.  The
main theorem then discharges that entry by the verified integral certificate
and even-contraction argument above. -/

/-- The literal assertion of Problem 2.1. -/
def Problem21Claim : Prop :=
  Tendsto (fun k => cfP a21 b21 k / cfQ a21 b21 k)
    atTop (𝓝 (6 / (3 - Real.pi)))

theorem problem21_pcf_value_of_cohen
    (hCohen : Tendsto (fun k => cfP cohenC cohenD k / cfQ cohenC cohenD k) atTop
      (𝓝 (6 / (Real.pi - 3)))) :
    Tendsto (fun k => cfP a21 b21 k / cfQ a21 b21 k) atTop (𝓝 (6 / (3 - Real.pi))) := by
  have hval : (6 : ℝ) / (3 - Real.pi) = -(6 / (Real.pi - 3)) := by
    rw [← neg_sub Real.pi 3, div_neg]
  rw [hval]
  exact hCohen.neg.congr fun k => (challenge_convergent_eq k).symm

theorem problem21_pcf_value_of_apery_diagonal
    (hAperyDiagonal : Tendsto
      (fun n => cfP aperyC21 aperyD21 (2 * n + 1) /
        cfQ aperyC21 aperyD21 (2 * n + 1))
      atTop (𝓝 Real.pi)) :
    Tendsto (fun k => cfP a21 b21 k / cfQ a21 b21 k)
      atTop (𝓝 (6 / (3 - Real.pi))) :=
  problem21_pcf_value_of_cohen
    (cohen_converges_of_apery_diagonal hAperyDiagonal)

/-- Problem 2.1, proved unconditionally. -/
theorem problem21_pcf_value : Problem21Claim := by
  unfold Problem21Claim
  exact problem21_pcf_value_of_cohen cohen21_converges

/-- Backwards-compatible wrapper accepting convergence of the full period-two
sequence. -/
theorem problem21_pcf_value_of_full_apery
    (_hApery : Tendsto
      (fun k => cfP aperyC21 aperyD21 k / cfQ aperyC21 aperyD21 k)
      atTop (𝓝 Real.pi)) :
    Tendsto (fun k => cfP a21 b21 k / cfQ a21 b21 k)
      atTop (𝓝 (6 / (3 - Real.pi))) := by
  exact problem21_pcf_value

end RamanujanChallenge.P21

end
