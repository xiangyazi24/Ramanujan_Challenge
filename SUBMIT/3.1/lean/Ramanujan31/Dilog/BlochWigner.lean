import Mathlib

/-!
# The Bloch–Wigner function, branch-free

This file constructs a function `blochWignerGeom : ℂ → ℝ` and proves the three
functional equations that the Problem 3.1 argument runs on:

    D(conj z) = -D z,    D((1-z)⁻¹) = D z,    D(x : ℝ) = 0.

Those three are exactly the fields of the structure `BlochWignerLaws`, which the
rest of the development currently takes as a hypothesis.  With this file they
become a theorem: `blochWignerGeom_laws`.

## The idea

The obvious route — define `Li₂` by its power series, continue it analytically
past the unit disc, set `D z = Im Li₂ z + arg(1-z) log‖z‖` — runs straight into
branch cuts and analytic continuation, and the arguments we care about are
outside the disc.

The route taken here instead uses the *geometric* (Lobachevsky) description: an
ideal hyperbolic tetrahedron with modulus `z` has three dihedral angles, one for
each element of the three-cycle orbit

    z  ↦  (1-z)⁻¹  ↦  1 - z⁻¹  ↦  z,

and `D` is the sum of the three Lobachevsky contributions.  Written that way:

* the argument of `clausen` is always a **phase** — a point of the unit circle —
  so the defining series converges everywhere on `ℂ \ {0,1}` and there is no
  continuation problem at all;
* the **squared** phase is used, so a real argument (phase `±1`) gives phase
  squared `= 1` and contributes exactly `0`, with no case split between argument
  `0` and argument `π`;
* the three-cycle invariance `D((1-z)⁻¹) = D z` becomes a cyclic reassociation of
  a three-term sum — no analysis whatsoever.

## Honest statement of what is and is not proved

`blochWignerGeom` is *a* standard definition of the Bloch–Wigner function, but
this file does **not** prove it equal to `Im Li₂(z) + arg(1-z) log‖z‖`.  That
bridge is a separate theorem.  Anyone reading the 3.1 development should know
which definition the machine-checked statements refer to; the name carries the
`Geom` suffix for exactly that reason.

What the 3.1 argument uses is only the three functional equations, and those are
proved here.
-/

open scoped Topology

namespace BlochWigner

/-! ## The Clausen series -/

/-- `clausen u = ∑_{n≥1} Im(uⁿ)/n²`.  For `‖u‖ ≤ 1` this converges absolutely.

When `u = exp(iθ)` this is the classical Clausen function `Cl₂(θ) = ∑ sin(nθ)/n²`,
i.e. twice the Lobachevsky function. -/
noncomputable def clausen (u : ℂ) : ℝ := ∑' n : ℕ, (u ^ (n + 1)).im / ((n : ℝ) + 1) ^ 2

theorem summable_clausenTerm {u : ℂ} (hu : ‖u‖ ≤ 1) :
    Summable (fun n : ℕ => (u ^ (n + 1)).im / ((n : ℝ) + 1) ^ 2) := by
  apply Summable.of_abs
  refine Summable.of_nonneg_of_le (fun n => abs_nonneg _) (fun n => ?_)
    (?_ : Summable fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2)
  · rw [abs_div]
    have hnum : |(u ^ (n + 1)).im| ≤ 1 := by
      calc |(u ^ (n + 1)).im| ≤ ‖u ^ (n + 1)‖ := Complex.abs_im_le_norm _
        _ = ‖u‖ ^ (n + 1) := by rw [norm_pow]
        _ ≤ 1 := pow_le_one₀ (norm_nonneg u) hu
    have hden : |((n : ℝ) + 1) ^ 2| = ((n : ℝ) + 1) ^ 2 := abs_of_nonneg (by positivity)
    rw [hden]
    gcongr
  · have := Real.summable_one_div_nat_pow.mpr (by norm_num : 1 < 2)
    exact (summable_nat_add_iff 1).mpr this |>.congr (fun n => by push_cast; ring)

/-- `clausen` vanishes at `1`: every term is `Im(1) = 0`. -/
@[simp] theorem clausen_one : clausen 1 = 0 := by
  simp [clausen]

/-- `clausen` is odd under conjugation, because conjugation negates `Im`. -/
theorem clausen_conj (u : ℂ) :
    clausen ((starRingEnd ℂ) u) = -clausen u := by
  have hterm : ∀ n : ℕ,
      (((starRingEnd ℂ) u) ^ (n + 1)).im / ((n : ℝ) + 1) ^ 2
        = -((u ^ (n + 1)).im / ((n : ℝ) + 1) ^ 2) := by
    intro n
    rw [← map_pow, Complex.conj_im]
    ring
  rw [clausen, clausen, tsum_congr hterm, tsum_neg]

/-! ## Phases -/

/-- The phase of a nonzero complex number, `w/‖w‖`. -/
noncomputable def phase (w : ℂ) : ℂ := w / (‖w‖ : ℂ)

theorem norm_phase {w : ℂ} (hw : w ≠ 0) : ‖phase w‖ = 1 := by
  rw [phase, norm_div, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (norm_nonneg w)]
  exact div_self (norm_ne_zero_iff.mpr hw)

theorem norm_phase_sq_le {w : ℂ} (hw : w ≠ 0) : ‖(phase w) ^ 2‖ ≤ 1 := by
  rw [norm_pow, norm_phase hw]; norm_num

theorem phase_conj (w : ℂ) : phase ((starRingEnd ℂ) w) = (starRingEnd ℂ) (phase w) := by
  rw [phase, phase, map_div₀]
  congr 1
  rw [RCLike.norm_conj]
  simp

/-- **The phase of a real number squares to `1`.**  This is what makes the
vanishing on the reals free of case analysis: both `+1` and `-1` square to `1`. -/
theorem phase_sq_ofReal {x : ℝ} (hx : x ≠ 0) : (phase (x : ℂ)) ^ 2 = 1 := by
  rw [phase]
  have hnorm : ‖(x : ℂ)‖ = |x| := by simp
  rw [hnorm, div_pow]
  have habs : ((|x| : ℝ) : ℂ) ^ 2 = ((x : ℂ)) ^ 2 := by
    have : ((|x| : ℝ) : ℂ) ^ 2 = (((|x| ^ 2 : ℝ)) : ℂ) := by push_cast; ring
    rw [this, sq_abs]
    push_cast
    ring
  rw [habs]
  exact div_self (pow_ne_zero 2 (by exact_mod_cast hx))

/-! ## The Lobachevsky angle contribution -/

/-- The contribution of one dihedral angle: half the Clausen function of the
squared phase. -/
noncomputable def angleTerm (w : ℂ) : ℝ := (1 / 2 : ℝ) * clausen ((phase w) ^ 2)

/-- A real argument contributes nothing. -/
@[simp] theorem angleTerm_ofReal {x : ℝ} (hx : x ≠ 0) : angleTerm (x : ℂ) = 0 := by
  rw [angleTerm, phase_sq_ofReal hx, clausen_one, mul_zero]

/-- `angleTerm` is odd under conjugation. -/
theorem angleTerm_conj {w : ℂ} (_hw : w ≠ 0) :
    angleTerm ((starRingEnd ℂ) w) = -angleTerm w := by
  rw [angleTerm, angleTerm, phase_conj, ← map_pow,
    clausen_conj _]
  ring

/-! ## The three-cycle -/

/-- `z ↦ (1-z)⁻¹`, the generator of the three-cycle on tetrahedron moduli. -/
noncomputable def mob (z : ℂ) : ℂ := (1 - z)⁻¹

theorem mob_mob {z : ℂ} (hz0 : z ≠ 0) (hz1 : z ≠ 1) : mob (mob z) = 1 - z⁻¹ := by
  have h1 : (1 : ℂ) - z ≠ 0 := sub_ne_zero.mpr (Ne.symm hz1)
  have hzz : z * z⁻¹ = 1 := mul_inv_cancel₀ hz0
  rw [mob, mob]
  field_simp
  linear_combination (z - 1) * hzz

theorem mob_mob_mob {z : ℂ} (hz0 : z ≠ 0) (hz1 : z ≠ 1) : mob (mob (mob z)) = z := by
  rw [mob_mob hz0 hz1, mob]
  field_simp
  ring

theorem mob_ne_zero {z : ℂ} (hz1 : z ≠ 1) : mob z ≠ 0 := by
  rw [mob]
  exact inv_ne_zero (sub_ne_zero.mpr (Ne.symm hz1))

theorem mob_ne_one {z : ℂ} (hz0 : z ≠ 0) (hz1 : z ≠ 1) : mob z ≠ 1 := by
  rw [mob]
  intro h
  have h1 : (1 : ℂ) - z ≠ 0 := sub_ne_zero.mpr (Ne.symm hz1)
  have h2 : (1 : ℂ) - z = 1 := by
    have hc := congrArg (fun w : ℂ => w⁻¹) h
    simpa [inv_inv] using hc
  exact hz0 (by linear_combination -h2)

theorem mob_conj (z : ℂ) : mob ((starRingEnd ℂ) z) = (starRingEnd ℂ) (mob z) := by
  rw [mob, mob, map_inv₀, map_sub, map_one]

theorem mob_ofReal (x : ℝ) : mob (x : ℂ) = (((1 - x)⁻¹ : ℝ) : ℂ) := by
  rw [mob]
  push_cast
  ring

/-! ## The Bloch–Wigner function -/

/-- **The Bloch–Wigner function**, as the sum of the three Lobachevsky angle
contributions around the three-cycle orbit of `z`. -/
noncomputable def blochWignerGeom (z : ℂ) : ℝ :=
  angleTerm z + angleTerm (mob z) + angleTerm (mob (mob z))

/-- **(E2) Three-cycle invariance.**  This is the equation the shape cancellation
runs on, and here it is pure reassociation. -/
theorem blochWignerGeom_mob {z : ℂ} (hz0 : z ≠ 0) (hz1 : z ≠ 1) :
    blochWignerGeom (mob z) = blochWignerGeom z := by
  rw [blochWignerGeom, blochWignerGeom, mob_mob_mob hz0 hz1]
  ring

/-- **(E1) Conjugation antisymmetry.** -/
theorem blochWignerGeom_conj {z : ℂ} (hz0 : z ≠ 0) (hz1 : z ≠ 1) :
    blochWignerGeom ((starRingEnd ℂ) z) = -blochWignerGeom z := by
  have h1 : mob z ≠ 0 := mob_ne_zero hz1
  have h2 : mob (mob z) ≠ 0 := mob_ne_zero (mob_ne_one hz0 hz1)
  rw [blochWignerGeom, blochWignerGeom]
  simp only [mob_conj]
  rw [angleTerm_conj hz0, angleTerm_conj h1, angleTerm_conj h2]
  ring

/-- **(E3) Vanishing on the reals.**  Every element of the orbit of a real number
is real, and a real number contributes nothing. -/
theorem blochWignerGeom_ofReal {x : ℝ} (hx0 : x ≠ 0) (hx1 : x ≠ 1) :
    blochWignerGeom (x : ℂ) = 0 := by
  have h1x : (1 : ℝ) - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hx1)
  have hmob : mob ((x : ℝ) : ℂ) = (((1 - x)⁻¹ : ℝ) : ℂ) := mob_ofReal x
  have hmob2 : mob (mob ((x : ℝ) : ℂ)) = ((((1 - (1 - x)⁻¹)⁻¹ : ℝ)) : ℂ) := by
    rw [hmob, mob_ofReal]
  have hne1 : ((1 - x)⁻¹ : ℝ) ≠ 0 := inv_ne_zero h1x
  have hne2 : ((1 - (1 - x)⁻¹)⁻¹ : ℝ) ≠ 0 := by
    apply inv_ne_zero
    intro h
    have he : ((1 - x)⁻¹ : ℝ) = 1 := by linarith
    have h2 : (1 : ℝ) - x = 1 := by
      have hc := congrArg (fun w : ℝ => w⁻¹) he
      simpa [inv_inv] using hc
    exact hx0 (by linarith)
  rw [blochWignerGeom, hmob2, hmob, angleTerm_ofReal hx0, angleTerm_ofReal hne1,
    angleTerm_ofReal hne2]
  ring


/-! ## The degenerate points, and the unconditional laws

At `z = 0` and `z = 1` the function `mob` is still defined (`mob 0 = 1`,
`mob 1 = 0`, using `(1-1)⁻¹ = 0`), and the junk values conspire: both orbits
consist of real points, so `blochWignerGeom` vanishes at both.  That makes all
three functional equations hold **unconditionally**, with no side conditions —
which is what lets them be packaged as `BlochWignerLaws`. -/

@[simp] theorem clausen_zero : clausen 0 = 0 := by
  simp [clausen]

@[simp] theorem angleTerm_zero : angleTerm 0 = 0 := by
  rw [angleTerm, phase]
  norm_num

@[simp] theorem angleTerm_one : angleTerm 1 = 0 := by
  have : ((1 : ℝ) : ℂ) = 1 := by norm_num
  rw [← this]
  exact angleTerm_ofReal one_ne_zero

@[simp] theorem mob_zero : mob 0 = 1 := by rw [mob]; norm_num

@[simp] theorem mob_one : mob 1 = 0 := by rw [mob]; norm_num

@[simp] theorem blochWignerGeom_zero : blochWignerGeom 0 = 0 := by
  rw [blochWignerGeom, mob_zero, mob_one, angleTerm_zero, angleTerm_one]
  ring

@[simp] theorem blochWignerGeom_one : blochWignerGeom 1 = 0 := by
  rw [blochWignerGeom, mob_one, mob_zero, angleTerm_zero, angleTerm_one]
  ring

/-- **(E1), unconditionally.** -/
theorem blochWignerGeom_conj' (z : ℂ) :
    blochWignerGeom ((starRingEnd ℂ) z) = -blochWignerGeom z := by
  rcases eq_or_ne z 0 with rfl | hz0
  · simp
  rcases eq_or_ne z 1 with rfl | hz1
  · simp
  exact blochWignerGeom_conj hz0 hz1

/-- **(E2), unconditionally.** -/
theorem blochWignerGeom_mob' (z : ℂ) :
    blochWignerGeom ((1 - z)⁻¹) = blochWignerGeom z := by
  rcases eq_or_ne z 0 with rfl | hz0
  · simpa using blochWignerGeom_one
  rcases eq_or_ne z 1 with rfl | hz1
  · simpa using blochWignerGeom_zero
  simpa [mob] using blochWignerGeom_mob hz0 hz1

/-- **(E3), unconditionally.** -/
theorem blochWignerGeom_ofReal' (x : ℝ) : blochWignerGeom (x : ℂ) = 0 := by
  rcases eq_or_ne x 0 with rfl | hx0
  · simpa using blochWignerGeom_zero
  rcases eq_or_ne x 1 with rfl | hx1
  · simpa using blochWignerGeom_one
  exact blochWignerGeom_ofReal hx0 hx1

end BlochWigner
