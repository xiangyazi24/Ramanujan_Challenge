import Mathlib.Data.Rat.Cast.Order
import Mathlib.Tactic

/-!
# Rational reconstruction from a denominator bound

The last step of the Ramanujan Challenge 3.1 proof is:

* the quantity `Re[Delta R]/pi^2` is known to be RATIONAL with denominator
  dividing `Q = 2040` (from the Merkurjev--Suslin torsion order of the two
  endpoint fields);
* it is computed numerically to `301` digits;
* therefore it equals `-4/85`.

This file proves the general principle that makes that step rigorous: two
distinct rationals with denominator at most `Q` are separated by at least
`1/Q^2`, so an approximation of error less than `1/Q^2` determines the rational
uniquely.

Nothing here is specific to the knot problem; the statement is reusable.

Main results:

* `Rat.inv_den_mul_den_le_abs_sub` — distinct rationals are separated by the
  reciprocal product of their denominators;
* `Rat.eq_of_den_le_of_abs_sub_lt_inv_sq` — the `1/Q^2` separation;
* `rat_reconstruct` — the reconstruction principle;
* `rat_reconstruct_half` — the `1/(2Q^2)` form actually used, which also gives
  the two-sided margin used in the write-up.
-/

namespace Rat

/-- Distinct rationals are separated by the reciprocal product of their
denominators. -/
theorem inv_den_mul_den_le_abs_sub {r s : ℚ} (hrs : r ≠ s) :
    1 / ((r.den : ℝ) * (s.den : ℝ)) ≤ |(r : ℝ) - (s : ℝ)| := by
  have hrden : (0 : ℝ) < r.den := by exact_mod_cast r.den_pos
  have hsden : (0 : ℝ) < s.den := by exact_mod_cast s.den_pos
  have hcross : r.num * (s.den : ℤ) - s.num * (r.den : ℤ) ≠ 0 := by
    intro h
    exact hrs (Rat.eq_iff_mul_eq_mul.mpr (sub_eq_zero.mp h))
  have hcross_abs :
      (1 : ℝ) ≤ |((r.num * (s.den : ℤ) - s.num * (r.den : ℤ) : ℤ) : ℝ)| := by
    exact_mod_cast Int.one_le_abs hcross
  rw [Rat.cast_def, Rat.cast_def]
  have hfrac :
      (r.num : ℝ) / (r.den : ℝ) - (s.num : ℝ) / (s.den : ℝ) =
        ((r.num * (s.den : ℤ) - s.num * (r.den : ℤ) : ℤ) : ℝ) /
          ((r.den : ℝ) * (s.den : ℝ)) := by
    push_cast
    field_simp
  rw [hfrac, abs_div, abs_of_pos (mul_pos hrden hsden)]
  exact (div_le_div_iff_of_pos_right (mul_pos hrden hsden)).2 hcross_abs

/-- Two rationals of denominator at most `Q` cannot be closer than `1/Q^2`. -/
theorem eq_of_den_le_of_abs_sub_lt_inv_sq {r s : ℚ} {Q : ℕ}
    (hQ : 0 < Q) (hr : r.den ≤ Q) (hs : s.den ≤ Q)
    (hclose : |(r : ℝ) - (s : ℝ)| < 1 / (Q : ℝ) ^ 2) :
    r = s := by
  by_contra hrs
  have hrden : (0 : ℝ) < r.den := by exact_mod_cast r.den_pos
  have hsden : (0 : ℝ) < s.den := by exact_mod_cast s.den_pos
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  have hrQ : (r.den : ℝ) ≤ Q := by exact_mod_cast hr
  have hsQ : (s.den : ℝ) ≤ Q := by exact_mod_cast hs
  have hprod : (r.den : ℝ) * (s.den : ℝ) ≤ (Q : ℝ) ^ 2 := by
    calc (r.den : ℝ) * (s.den : ℝ) ≤ (Q : ℝ) * Q :=
          mul_le_mul hrQ hsQ hsden.le hQr.le
      _ = (Q : ℝ) ^ 2 := by ring
  have hinv : 1 / (Q : ℝ) ^ 2 ≤ 1 / ((r.den : ℝ) * (s.den : ℝ)) :=
    (one_div_le_one_div (pow_pos hQr 2) (mul_pos hrden hsden)).2 hprod
  exact (not_lt_of_ge (hinv.trans (Rat.inv_den_mul_den_le_abs_sub hrs))) hclose

end Rat

/-- **Rational reconstruction.** A real number known to be rational with
denominator at most `Q`, and known to within `1/Q^2` of a rational `q` of
denominator at most `Q`, equals `q`. -/
theorem rat_reconstruct {x : ℝ} {q : ℚ} {Q : ℕ}
    (hQ : 0 < Q) (hq : q.den ≤ Q)
    (hx : ∃ r : ℚ, r.den ≤ Q ∧ x = (r : ℝ))
    (happrox : |x - (q : ℝ)| < 1 / (Q : ℝ) ^ 2) :
    x = (q : ℝ) := by
  rcases hx with ⟨r, hr, rfl⟩
  have hrq : r = q := Rat.eq_of_den_le_of_abs_sub_lt_inv_sq hQ hr hq happrox
  simp [hrq]

/-- The `1/(2Q^2)` form, which is what the numerical certificate supplies. -/
theorem rat_reconstruct_half {x : ℝ} {q : ℚ} {Q : ℕ}
    (hQ : 0 < Q) (hq : q.den ≤ Q)
    (hx : ∃ r : ℚ, r.den ≤ Q ∧ x = (r : ℝ))
    (happrox : |x - (q : ℝ)| < 1 / (2 * (Q : ℝ) ^ 2)) :
    x = (q : ℝ) := by
  refine rat_reconstruct hQ hq hx (lt_of_lt_of_le happrox ?_)
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  have h2 : (0 : ℝ) < (Q : ℝ) ^ 2 := pow_pos hQr 2
  have hle : (Q : ℝ) ^ 2 ≤ 2 * (Q : ℝ) ^ 2 := by linarith
  exact one_div_le_one_div_of_le h2 hle

/-!
## The instance used in Problem 3.1

`Q = 2040 = lcm (w₂ (F_alpha), w₂ (F_beta)) = lcm 120 408`, and the target
rational is `-4/85`.  Note `(-4/85 : ℚ).den = 85 ≤ 2040`, and indeed `85 ∣ 2040`
with quotient `24`.
-/

example : ((-4 : ℚ) / 85).den = 85 := by norm_num

example : (2040 : ℕ) % 85 = 0 := by norm_num

/-- The concrete reconstruction step of Problem 3.1: if the regulator quotient is
rational with denominator at most `2040`, and lies within `1/(2·2040²)` of
`-4/85`, then it *is* `-4/85`. -/
theorem regulator_quotient_eq
    {x : ℝ}
    (hx : ∃ r : ℚ, r.den ≤ 2040 ∧ x = (r : ℝ))
    (happrox : |x - ((-4 : ℚ) / 85 : ℚ)| < 1 / (2 * (2040 : ℝ) ^ 2)) :
    x = ((-4 : ℚ) / 85 : ℚ) := by
  refine rat_reconstruct_half (by norm_num) ?_ hx happrox
  norm_num

/-!
## Robustness of the reconstruction

A referee objection: the passage "the extended Bloch class is torsion of order
dividing `m`" ⟹ "the Rogers value lies in `(1/m)π²ℤ`" carries convention-dependent
factors — 6, 12, 24 — so the bound `Q = 2040` might really be `6·2040` or
`24·2040`.

**In this instance the objection is harmless, and the theorem below says exactly
why.**  The numerical certificate has error `1.1·10⁻³⁰¹`, while the separation of
rationals with denominator `≤ Q` is `1/Q²`.  So the reconstruction survives *any*
denominator bound up to about `1.7·10¹⁵⁰`; a stray factor of 6 or 24, or even of
`10¹⁴⁶`, changes nothing.

The theorem is stated with the generous ceiling `Q ≤ 10¹⁵⁰` so that no
normalization question has to be settled to believe the conclusion.
-/

/-- **The reconstruction is insensitive to the normalization of the denominator
bound.**  For *any* `Q` with `85 ≤ Q ≤ 10¹⁵⁰`, a quantity known to be rational
with denominator at most `Q`, and known to agree with `-4/85` to `301` digits,
equals `-4/85`. -/
theorem regulator_quotient_eq_robust {x : ℝ} {Q : ℕ}
    (hQ : 85 ≤ Q) (hQle : Q ≤ 10 ^ 150)
    (hx : ∃ r : ℚ, r.den ≤ Q ∧ x = (r : ℝ))
    (happrox : |x - ((-4 : ℚ) / 85 : ℚ)| < 1 / (10 : ℝ) ^ 301) :
    x = ((-4 : ℚ) / 85 : ℚ) := by
  have hQpos : 0 < Q := lt_of_lt_of_le (by norm_num) hQ
  have hQr : (0 : ℝ) < (Q : ℝ) := by exact_mod_cast hQpos
  refine rat_reconstruct hQpos ?_ hx (lt_of_lt_of_le happrox ?_)
  · rw [show ((-4 : ℚ) / 85).den = 85 from by norm_num]
    exact hQ
  · have h1 : (Q : ℝ) ≤ (10 : ℝ) ^ 150 := by exact_mod_cast hQle
    have hmul : (Q : ℝ) * (Q : ℝ) ≤ (10 : ℝ) ^ 150 * (10 : ℝ) ^ 150 :=
      mul_self_le_mul_self hQr.le h1
    have h2 : (Q : ℝ) ^ 2 ≤ (10 : ℝ) ^ 301 := by
      calc (Q : ℝ) ^ 2 = (Q : ℝ) * (Q : ℝ) := by ring
        _ ≤ (10 : ℝ) ^ 150 * (10 : ℝ) ^ 150 := hmul
        _ = (10 : ℝ) ^ 300 := by rw [← pow_add]
        _ ≤ (10 : ℝ) ^ 301 := pow_le_pow_right₀ (by norm_num) (by norm_num)
    exact one_div_le_one_div_of_le (pow_pos hQr 2) h2

/-- The concrete instance with the *stated* bound `Q = 2040`, and with the bound
inflated by the largest conventional factor anyone might insist on (`24`), both
give the same answer. -/
example {x : ℝ}
    (hx : ∃ r : ℚ, r.den ≤ 24 * 2040 ∧ x = (r : ℝ))
    (happrox : |x - ((-4 : ℚ) / 85 : ℚ)| < 1 / (10 : ℝ) ^ 301) :
    x = ((-4 : ℚ) / 85 : ℚ) :=
  regulator_quotient_eq_robust (by norm_num) (by norm_num) hx happrox
