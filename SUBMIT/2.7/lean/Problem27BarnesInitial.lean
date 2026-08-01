import RamanujanChallenge.Problem27BarnesPole

open Filter Set MeasureTheory Asymptotics Topology

noncomputable section

namespace RamanujanChallenge.P27

private theorem barnesPolePoint_ne_zero27 (m : ℕ) (y : ℝ) :
    barnesPolePoint27 m y ≠ 0 := by
  intro h
  have hre := congrArg Complex.re h
  simp [barnesPolePoint27] at hre
  linarith

private theorem integral_y_sechSq27 :
    (∫ y : ℝ, (y : ℂ) * sechSq27 y) = 0 := by
  let f : ℝ → ℂ := fun y => (y : ℂ) * sechSq27 y
  have hneg := (Measure.measurePreserving_neg (volume : Measure ℝ)).integral_comp
    (Homeomorph.neg ℝ).measurableEmbedding f
  have hfun : (fun y : ℝ => f (-y)) = fun y => -f y := by
    funext y
    unfold f sechSq27
    rw [show Real.pi * -y = -(Real.pi * y) by ring]
    rw [Real.cosh_neg]
    push_cast
    ring
  rw [hfun, MeasureTheory.integral_neg] at hneg
  change (∫ y : ℝ, (y : ℂ) * sechSq27 y) = 0
  change -(∫ y : ℝ, (y : ℂ) * sechSq27 y) =
    ∫ y : ℝ, (y : ℂ) * sechSq27 y at hneg
  linear_combination -hneg / 2

private theorem barnesPhi_zero27 (s : ℂ) (hs : -(1 / 2 : ℝ) ≤ s.re) :
    zudilinBarnesPhi27 0 s = barnesPoleBlock27 (s + 1) := by
  have hd : s + 1 ≠ 0 := by
    simpa [zudilinBarnesDenominatorFactor27] using
      (zudilinBarnesDenominatorFactor_ne_zero27 (n := 0) (k := 0) hs)
  rw [zudilinBarnesPhi27,
    zudilinBarnesFPrime_eq_mul_logDerivative27 hs]
  simp [zudilinBarnesF27, zudilinBarnesLogDerivative27,
    zudilinBarnesNumeratorFactor27, zudilinBarnesDenominatorFactor27,
    barnesPoleBlock27, Finset.prod_range_succ, Finset.sum_range_succ]
  field_simp [hd]
  ring

private theorem barnesPhi_one27 (s : ℂ) (hs : -(1 / 2 : ℝ) ≤ s.re) :
    zudilinBarnesPhi27 1 s =
      (s + 2) - 9 / 2 - barnesPoleBlock27 (s + 2) +
        8 * barnesPoleBlock27 (s + 3) := by
  have hn : s + 1 ≠ 0 := by
    simpa [zudilinBarnesNumeratorFactor27] using
      (zudilinBarnesNumeratorFactor_ne_zero27 (k := 0) hs)
  have hd0 : s + 2 ≠ 0 := by
    simpa [zudilinBarnesDenominatorFactor27] using
      (zudilinBarnesDenominatorFactor_ne_zero27 (n := 1) (k := 0) hs)
  have hd1 : s + 3 ≠ 0 := by
    simpa [zudilinBarnesDenominatorFactor27] using
      (zudilinBarnesDenominatorFactor_ne_zero27 (n := 1) (k := 1) hs)
  have hd1' : 3 + s ≠ 0 := by simpa [add_comm] using hd1
  rw [zudilinBarnesPhi27,
    zudilinBarnesFPrime_eq_mul_logDerivative27 hs]
  simp only [zudilinBarnesF27, zudilinBarnesLogDerivative27,
    zudilinBarnesNumeratorFactor27, zudilinBarnesDenominatorFactor27,
    Finset.prod_range_succ, Finset.prod_range_zero, Finset.sum_range_succ,
    Finset.sum_range_zero, Nat.factorial_one, Nat.cast_one, one_pow, one_mul]
  norm_num
  change
    ((s + 1) ^ 3 / ((s + 2) * (s + 3)) -
      ((s + 1) ^ 3 / ((s + 2) * (s + 3))) *
        (3 * (s + 1)⁻¹ - ((s + 2)⁻¹ + (s + 3)⁻¹)) / 2) = _
  unfold barnesPoleBlock27
  field_simp [hn, hd0, hd1, hd1']
  ring

private theorem barnesPhi_two27 (s : ℂ) (hs : -(1 / 2 : ℝ) ≤ s.re) :
    zudilinBarnesPhi27 2 s =
      ((1 / 4 : ℂ) * (s + 3) ^ 3 - (27 / 8 : ℂ) * (s + 3) ^ 2 +
        (79 / 4 : ℂ) * (s + 3) - 547 / 8) +
        barnesPoleBlock27 (s + 3) - 54 * barnesPoleBlock27 (s + 4) +
        216 * barnesPoleBlock27 (s + 5) := by
  have hn0 : s + 1 ≠ 0 := by
    simpa [zudilinBarnesNumeratorFactor27] using
      (zudilinBarnesNumeratorFactor_ne_zero27 (k := 0) hs)
  have hn1 : s + 2 ≠ 0 := by
    simpa [zudilinBarnesNumeratorFactor27] using
      (zudilinBarnesNumeratorFactor_ne_zero27 (k := 1) hs)
  have hd0 : s + 3 ≠ 0 := by
    simpa [zudilinBarnesDenominatorFactor27] using
      (zudilinBarnesDenominatorFactor_ne_zero27 (n := 2) (k := 0) hs)
  have hd1 : s + 4 ≠ 0 := by
    simpa [zudilinBarnesDenominatorFactor27] using
      (zudilinBarnesDenominatorFactor_ne_zero27 (n := 2) (k := 1) hs)
  have hd2 : s + 5 ≠ 0 := by
    simpa [zudilinBarnesDenominatorFactor27] using
      (zudilinBarnesDenominatorFactor_ne_zero27 (n := 2) (k := 2) hs)
  have hn1' : 2 + s ≠ 0 := by simpa [add_comm] using hn1
  have hd0' : 3 + s ≠ 0 := by simpa [add_comm] using hd0
  have hd1' : 4 + s ≠ 0 := by simpa [add_comm] using hd1
  have hd2' : 5 + s ≠ 0 := by simpa [add_comm] using hd2
  have hpoly : 20 + s * 9 + s ^ 2 ≠ 0 := by
    have hmul := mul_ne_zero hd1' hd2'
    convert hmul using 1 <;> ring
  rw [zudilinBarnesPhi27,
    zudilinBarnesFPrime_eq_mul_logDerivative27 hs]
  simp only [zudilinBarnesF27, zudilinBarnesLogDerivative27,
    zudilinBarnesNumeratorFactor27, zudilinBarnesDenominatorFactor27,
    Finset.prod_range_succ, Finset.prod_range_zero, Finset.sum_range_succ,
    Finset.sum_range_zero, Nat.factorial_succ, Nat.factorial_one,
    Nat.cast_ofNat, Nat.cast_one, one_mul]
  norm_num
  change
    ((((s + 1) * (s + 2)) ^ 3 /
        (4 * ((s + 3) * (s + 4) * (s + 5)))) -
      (((s + 1) * (s + 2)) ^ 3 /
        (4 * ((s + 3) * (s + 4) * (s + 5)))) *
        (3 * ((s + 1)⁻¹ + (s + 2)⁻¹) -
          ((s + 3)⁻¹ + (s + 4)⁻¹ + (s + 5)⁻¹)) / 2) = _
  unfold barnesPoleBlock27
  field_simp [hn0, hn1, hd0, hd1, hd2, hn1', hd0', hd1', hd2', hpoly]
  ring

private theorem barnesPhi_zero_line27 (y : ℝ) :
    zudilinBarnesPhi27 0 (zudilinBarnesLine27 y) =
      barnesPoleBlock27 (barnesPolePoint27 0 y) := by
  rw [barnesPhi_zero27 _ (by simp)]
  apply congrArg barnesPoleBlock27
  unfold zudilinBarnesLine27 barnesPolePoint27
  push_cast
  ring

private theorem barnesPhi_one_line27 (y : ℝ) :
    zudilinBarnesPhi27 1 (zudilinBarnesLine27 y) =
      (-3 : ℂ) + (y : ℂ) * Complex.I -
        barnesPoleBlock27 (barnesPolePoint27 1 y) +
          8 * barnesPoleBlock27 (barnesPolePoint27 2 y) := by
  rw [barnesPhi_one27 _ (by simp)]
  simp [zudilinBarnesLine27, barnesPolePoint27]
  ring

private theorem barnesPhi_two_line27 (y : ℝ) :
    zudilinBarnesPhi27 2 (zudilinBarnesLine27 y) =
      ((1 / 4 : ℂ) * (barnesPolePoint27 2 y) ^ 3 -
        (27 / 8 : ℂ) * (barnesPolePoint27 2 y) ^ 2 +
        (79 / 4 : ℂ) * barnesPolePoint27 2 y - 547 / 8) +
        barnesPoleBlock27 (barnesPolePoint27 2 y) -
        54 * barnesPoleBlock27 (barnesPolePoint27 3 y) +
        216 * barnesPoleBlock27 (barnesPolePoint27 4 y) := by
  rw [barnesPhi_two27 _ (by simp)]
  simp [zudilinBarnesLine27, barnesPolePoint27]
  ring

private theorem barnesRealIntegrand_zero27 (y : ℝ) :
    zudilinBarnesRealIntegrand27 0 y =
      sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 0 y) := by
  rw [zudilinBarnesRealIntegrand27, barnesPhi_zero_line27]
  unfold sechSq27
  ring

theorem zudilinBarnesErrorIntegral_zero27 :
    zudilinBarnesErrorIntegral27 0 =
      (((Real.pi ^ 2 / 6 + zeta3 : ℝ)) : ℂ) := by
  unfold zudilinBarnesErrorIntegral27
  rw [show (fun y : ℝ => zudilinBarnesRealIntegrand27 0 y) =
      fun y => sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 0 y) by
    funext y
    exact barnesRealIntegrand_zero27 y]
  rw [integral_barnesPoleBlock27, barnesPoleTail_zero27]

private theorem normalized_integral_linear27 :
    ((Real.pi / 2 : ℝ) : ℂ) *
        ∫ y : ℝ, ((-3 : ℂ) + (y : ℂ) * Complex.I) * sechSq27 y =
      -3 := by
  have h0 : Integrable (fun y : ℝ => (-3 : ℂ) * sechSq27 y) := by
    simpa using (integrable_pow_sechSq27 0).const_mul (-3 : ℂ)
  have h1 : Integrable
      (fun y : ℝ => ((y : ℂ) * sechSq27 y) * Complex.I) := by
    simpa using (integrable_pow_sechSq27 1).mul_const Complex.I
  rw [show (fun y : ℝ => ((-3 : ℂ) + (y : ℂ) * Complex.I) * sechSq27 y) =
      fun y => (-3 : ℂ) * sechSq27 y +
        ((y : ℂ) * sechSq27 y) * Complex.I by
    funext y
    ring]
  have hm : (∫ y : ℝ, (-3 : ℂ) * sechSq27 y) =
      (-3 : ℂ) * ∫ y : ℝ, sechSq27 y :=
    MeasureTheory.integral_const_mul (-3 : ℂ) sechSq27
  have hi : (∫ y : ℝ, ((y : ℂ) * sechSq27 y) * Complex.I) =
      (∫ y : ℝ, (y : ℂ) * sechSq27 y) * Complex.I :=
    MeasureTheory.integral_mul_const Complex.I
      (fun y : ℝ => (y : ℂ) * sechSq27 y)
  rw [MeasureTheory.integral_add h0 h1, hm, hi,
    integral_sechSq27, integral_y_sechSq27]
  push_cast
  field_simp [Real.pi_ne_zero]
  ring

private theorem barnesRealIntegrand_one27 (y : ℝ) :
    zudilinBarnesRealIntegrand27 1 y =
      (((-3 : ℂ) + (y : ℂ) * Complex.I) * sechSq27 y) -
        sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 1 y) +
        8 * (sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 2 y)) := by
  rw [zudilinBarnesRealIntegrand27, barnesPhi_one_line27]
  unfold sechSq27
  ring

theorem zudilinBarnesErrorIntegral_one27 :
    zudilinBarnesErrorIntegral27 1 =
      (((7 * (Real.pi ^ 2 / 6 + zeta3) - 20 : ℝ)) : ℂ) := by
  let p : ℝ → ℂ := fun y =>
    ((-3 : ℂ) + (y : ℂ) * Complex.I) * sechSq27 y
  let b1 : ℝ → ℂ := fun y =>
    sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 1 y)
  let b2 : ℝ → ℂ := fun y =>
    sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 2 y)
  have hp : Integrable p := by
    unfold p
    have h0 : Integrable (fun y : ℝ => (-3 : ℂ) * sechSq27 y) := by
      simpa using (integrable_pow_sechSq27 0).const_mul (-3 : ℂ)
    have h1 : Integrable
        (fun y : ℝ => ((y : ℂ) * sechSq27 y) * Complex.I) := by
      simpa using (integrable_pow_sechSq27 1).mul_const Complex.I
    apply (h0.add h1).congr
    filter_upwards with y
    change (-3 : ℂ) * sechSq27 y +
      ((y : ℂ) * sechSq27 y) * Complex.I = p y
    unfold p
    ring_nf
  have hb1 : Integrable b1 := integrable_sechSq_mul_barnesPoleBlock27 1
  have hb2 : Integrable b2 := integrable_sechSq_mul_barnesPoleBlock27 2
  unfold zudilinBarnesErrorIntegral27
  rw [show (fun y : ℝ => zudilinBarnesRealIntegrand27 1 y) =
      fun y => p y - b1 y + 8 * b2 y by
    funext y
    exact barnesRealIntegrand_one27 y]
  have hsplit :
      (∫ y : ℝ, p y - b1 y + 8 * b2 y) =
        (∫ y : ℝ, p y) - (∫ y : ℝ, b1 y) +
          8 * (∫ y : ℝ, b2 y) := by
    have hsub : (∫ y : ℝ, (p - b1) y) =
        (∫ y : ℝ, p y) - ∫ y : ℝ, b1 y := by
      simpa only [Pi.sub_apply] using MeasureTheory.integral_sub hp hb1
    have hmul : (∫ y : ℝ, (8 : ℂ) * b2 y) =
        8 * ∫ y : ℝ, b2 y :=
      MeasureTheory.integral_const_mul (8 : ℂ) b2
    calc
      (∫ y : ℝ, p y - b1 y + 8 * b2 y) =
          ∫ y : ℝ, (p - b1) y + (fun y => (8 : ℂ) * b2 y) y := rfl
      _ = (∫ y : ℝ, (p - b1) y) +
          ∫ y : ℝ, (8 : ℂ) * b2 y :=
        MeasureTheory.integral_add (hp.sub hb1) (hb2.const_mul 8)
      _ = (∫ y : ℝ, p y) - (∫ y : ℝ, b1 y) +
          8 * (∫ y : ℝ, b2 y) := by
        rw [hsub, hmul]
  rw [hsplit]
  have hpval := normalized_integral_linear27
  have hb1val := integral_barnesPoleBlock27 1
  have hb2val := integral_barnesPoleBlock27 2
  rw [barnesPoleTail_one27] at hb1val
  rw [barnesPoleTail_two27] at hb2val
  push_cast at hpval hb1val hb2val ⊢
  linear_combination hpval - hb1val + 8 * hb2val

private theorem integral_cube_sechSq27 :
    (∫ y : ℝ, (y : ℂ) ^ 3 * sechSq27 y) = 0 := by
  let f : ℝ → ℂ := fun y => (y : ℂ) ^ 3 * sechSq27 y
  have hneg := (Measure.measurePreserving_neg (volume : Measure ℝ)).integral_comp
    (Homeomorph.neg ℝ).measurableEmbedding f
  have hfun : (fun y : ℝ => f (-y)) = fun y => -f y := by
    funext y
    unfold f sechSq27
    rw [show Real.pi * -y = -(Real.pi * y) by ring]
    rw [Real.cosh_neg]
    push_cast
    ring
  rw [hfun, MeasureTheory.integral_neg] at hneg
  change (∫ y : ℝ, (y : ℂ) ^ 3 * sechSq27 y) = 0
  change -(∫ y : ℝ, (y : ℂ) ^ 3 * sechSq27 y) =
    ∫ y : ℝ, (y : ℂ) ^ 3 * sechSq27 y at hneg
  linear_combination -hneg / 2

private theorem barnesPolynomial_two_expand27 (y : ℝ) :
    ((1 / 4 : ℂ) * (barnesPolePoint27 2 y) ^ 3 -
        (27 / 8 : ℂ) * (barnesPolePoint27 2 y) ^ 2 +
        (79 / 4 : ℂ) * barnesPolePoint27 2 y - 547 / 8) * sechSq27 y =
      (sechSq27 y * (-579 / 16 : ℂ) +
        ((y : ℂ) * sechSq27 y) * ((121 / 16 : ℂ) * Complex.I)) +
      (((y : ℂ) ^ 2 * sechSq27 y) * (3 / 2 : ℂ) +
        ((y : ℂ) ^ 3 * sechSq27 y) * (-(Complex.I) / 4)) := by
  unfold barnesPolePoint27
  push_cast
  ring_nf
  rw [Complex.I_sq, Complex.I_pow_three]
  ring

private theorem integrable_barnesPolynomial_two27 :
    Integrable (fun y : ℝ =>
      ((1 / 4 : ℂ) * (barnesPolePoint27 2 y) ^ 3 -
          (27 / 8 : ℂ) * (barnesPolePoint27 2 y) ^ 2 +
          (79 / 4 : ℂ) * barnesPolePoint27 2 y - 547 / 8) * sechSq27 y) := by
  have h0 : Integrable (fun y : ℝ =>
      sechSq27 y * (-579 / 16 : ℂ)) := by
    simpa [Complex.real_smul] using
      (integrable_pow_sechSq27 0).mul_const (-579 / 16 : ℂ)
  have h1 : Integrable (fun y : ℝ =>
      ((y : ℂ) * sechSq27 y) * ((121 / 16 : ℂ) * Complex.I)) := by
    simpa [Complex.real_smul] using
      (integrable_pow_sechSq27 1).mul_const ((121 / 16 : ℂ) * Complex.I)
  have h2 : Integrable (fun y : ℝ =>
      ((y : ℂ) ^ 2 * sechSq27 y) * (3 / 2 : ℂ)) := by
    simpa [Complex.real_smul] using
      (integrable_pow_sechSq27 2).mul_const (3 / 2 : ℂ)
  have h3 : Integrable (fun y : ℝ =>
      ((y : ℂ) ^ 3 * sechSq27 y) * (-(Complex.I) / 4)) := by
    simpa [Complex.real_smul] using
      (integrable_pow_sechSq27 3).mul_const (-(Complex.I) / 4)
  apply ((h0.add h1).add (h2.add h3)).congr
  filter_upwards with y
  rw [barnesPolynomial_two_expand27]
  simp only [Pi.add_apply]

private theorem normalized_integral_barnesPolynomial_two27 :
    ((Real.pi / 2 : ℝ) : ℂ) *
        ∫ y : ℝ,
          ((1 / 4 : ℂ) * (barnesPolePoint27 2 y) ^ 3 -
            (27 / 8 : ℂ) * (barnesPolePoint27 2 y) ^ 2 +
            (79 / 4 : ℂ) * barnesPolePoint27 2 y - 547 / 8) * sechSq27 y =
      (-577 / 16 : ℂ) := by
  let f0 : ℝ → ℂ := fun y => sechSq27 y * (-579 / 16 : ℂ)
  let f1 : ℝ → ℂ := fun y =>
    ((y : ℂ) * sechSq27 y) * ((121 / 16 : ℂ) * Complex.I)
  let f2 : ℝ → ℂ := fun y =>
    ((y : ℂ) ^ 2 * sechSq27 y) * (3 / 2 : ℂ)
  let f3 : ℝ → ℂ := fun y =>
    ((y : ℂ) ^ 3 * sechSq27 y) * (-(Complex.I) / 4)
  have h0 : Integrable f0 := by
    simpa [f0] using (integrable_pow_sechSq27 0).mul_const (-579 / 16 : ℂ)
  have h1 : Integrable f1 :=
    by simpa [f1, Complex.real_smul] using
      (integrable_pow_sechSq27 1).mul_const ((121 / 16 : ℂ) * Complex.I)
  have h2 : Integrable f2 :=
    by simpa [f2, Complex.real_smul] using
      (integrable_pow_sechSq27 2).mul_const (3 / 2 : ℂ)
  have h3 : Integrable f3 :=
    by simpa [f3, Complex.real_smul] using
      (integrable_pow_sechSq27 3).mul_const (-(Complex.I) / 4)
  rw [show (fun y : ℝ =>
      ((1 / 4 : ℂ) * (barnesPolePoint27 2 y) ^ 3 -
        (27 / 8 : ℂ) * (barnesPolePoint27 2 y) ^ 2 +
        (79 / 4 : ℂ) * barnesPolePoint27 2 y - 547 / 8) * sechSq27 y) =
      fun y => (f0 + f1) y + (f2 + f3) y by
    funext y
    exact barnesPolynomial_two_expand27 y]
  rw [MeasureTheory.integral_add (h0.add h1) (h2.add h3)]
  have h01 : (∫ y : ℝ, (f0 + f1) y) =
      (∫ y : ℝ, f0 y) + ∫ y : ℝ, f1 y := by
    simpa only [Pi.add_apply] using MeasureTheory.integral_add h0 h1
  have h23 : (∫ y : ℝ, (f2 + f3) y) =
      (∫ y : ℝ, f2 y) + ∫ y : ℝ, f3 y := by
    simpa only [Pi.add_apply] using MeasureTheory.integral_add h2 h3
  rw [h01, h23]
  unfold f0 f1 f2 f3
  have hi0 : (∫ y : ℝ, sechSq27 y * (-579 / 16 : ℂ)) =
      (∫ y : ℝ, sechSq27 y) * (-579 / 16 : ℂ) :=
    MeasureTheory.integral_mul_const (-579 / 16 : ℂ) sechSq27
  have hi1 : (∫ y : ℝ, ((y : ℂ) * sechSq27 y) *
      ((121 / 16 : ℂ) * Complex.I)) =
      (∫ y : ℝ, (y : ℂ) * sechSq27 y) *
        ((121 / 16 : ℂ) * Complex.I) :=
    MeasureTheory.integral_mul_const ((121 / 16 : ℂ) * Complex.I)
      (fun y : ℝ => (y : ℂ) * sechSq27 y)
  have hi2 : (∫ y : ℝ, ((y : ℂ) ^ 2 * sechSq27 y) * (3 / 2 : ℂ)) =
      (∫ y : ℝ, (y : ℂ) ^ 2 * sechSq27 y) * (3 / 2 : ℂ) :=
    MeasureTheory.integral_mul_const (3 / 2 : ℂ)
      (fun y : ℝ => (y : ℂ) ^ 2 * sechSq27 y)
  have hi3 : (∫ y : ℝ, ((y : ℂ) ^ 3 * sechSq27 y) *
      (-(Complex.I) / 4)) =
      (∫ y : ℝ, (y : ℂ) ^ 3 * sechSq27 y) * (-(Complex.I) / 4) :=
    MeasureTheory.integral_mul_const (-(Complex.I) / 4)
      (fun y : ℝ => (y : ℂ) ^ 3 * sechSq27 y)
  rw [hi0, hi1, hi2, hi3,
    integral_sechSq27, integral_y_sechSq27, integral_sq_sechSq27,
    integral_cube_sechSq27]
  push_cast
  field_simp [Real.pi_ne_zero]
  ring

private theorem barnesRealIntegrand_two27 (y : ℝ) :
    zudilinBarnesRealIntegrand27 2 y =
      (((1 / 4 : ℂ) * (barnesPolePoint27 2 y) ^ 3 -
          (27 / 8 : ℂ) * (barnesPolePoint27 2 y) ^ 2 +
          (79 / 4 : ℂ) * barnesPolePoint27 2 y - 547 / 8) * sechSq27 y) +
        sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 2 y) -
        54 * (sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 3 y)) +
        216 * (sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 4 y)) := by
  rw [zudilinBarnesRealIntegrand27, barnesPhi_two_line27]
  unfold sechSq27
  ring

theorem zudilinBarnesErrorIntegral_two27 :
    zudilinBarnesErrorIntegral27 2 =
      (((163 * (Real.pi ^ 2 / 6 + zeta3) - 7425 / 16 : ℝ)) : ℂ) := by
  let p : ℝ → ℂ := fun y =>
    ((1 / 4 : ℂ) * (barnesPolePoint27 2 y) ^ 3 -
      (27 / 8 : ℂ) * (barnesPolePoint27 2 y) ^ 2 +
      (79 / 4 : ℂ) * barnesPolePoint27 2 y - 547 / 8) * sechSq27 y
  let b2 : ℝ → ℂ := fun y =>
    sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 2 y)
  let b3 : ℝ → ℂ := fun y =>
    sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 3 y)
  let b4 : ℝ → ℂ := fun y =>
    sechSq27 y * barnesPoleBlock27 (barnesPolePoint27 4 y)
  have hp : Integrable p := integrable_barnesPolynomial_two27
  have hb2 : Integrable b2 := integrable_sechSq_mul_barnesPoleBlock27 2
  have hb3 : Integrable b3 := integrable_sechSq_mul_barnesPoleBlock27 3
  have hb4 : Integrable b4 := integrable_sechSq_mul_barnesPoleBlock27 4
  unfold zudilinBarnesErrorIntegral27
  rw [show (fun y : ℝ => zudilinBarnesRealIntegrand27 2 y) =
      fun y => p y + b2 y - 54 * b3 y + 216 * b4 y by
    funext y
    exact barnesRealIntegrand_two27 y]
  have hsplit :
      (∫ y : ℝ, p y + b2 y - 54 * b3 y + 216 * b4 y) =
        (∫ y : ℝ, p y) + (∫ y : ℝ, b2 y) -
          54 * (∫ y : ℝ, b3 y) + 216 * (∫ y : ℝ, b4 y) := by
    have hpb : (∫ y : ℝ, p y + b2 y) =
        (∫ y : ℝ, p y) + ∫ y : ℝ, b2 y :=
      MeasureTheory.integral_add hp hb2
    have h3 : (∫ y : ℝ, (54 : ℂ) * b3 y) =
        54 * ∫ y : ℝ, b3 y :=
      MeasureTheory.integral_const_mul (54 : ℂ) b3
    have h4 : (∫ y : ℝ, (216 : ℂ) * b4 y) =
        216 * ∫ y : ℝ, b4 y :=
      MeasureTheory.integral_const_mul (216 : ℂ) b4
    have hsub : (∫ y : ℝ, p y + b2 y - 54 * b3 y) =
        (∫ y : ℝ, p y + b2 y) - ∫ y : ℝ, 54 * b3 y :=
      MeasureTheory.integral_sub (hp.add hb2) (hb3.const_mul 54)
    calc
      (∫ y : ℝ, p y + b2 y - 54 * b3 y + 216 * b4 y) =
          (∫ y : ℝ, p y + b2 y - 54 * b3 y) +
            ∫ y : ℝ, 216 * b4 y :=
        MeasureTheory.integral_add ((hp.add hb2).sub (hb3.const_mul 54))
          (hb4.const_mul 216)
      _ = ((∫ y : ℝ, p y + b2 y) - ∫ y : ℝ, 54 * b3 y) +
            ∫ y : ℝ, 216 * b4 y := by
        rw [hsub]
      _ = (∫ y : ℝ, p y) + (∫ y : ℝ, b2 y) -
          54 * (∫ y : ℝ, b3 y) + 216 * (∫ y : ℝ, b4 y) := by
        rw [hpb, h3, h4]
  rw [hsplit]
  have hpval := normalized_integral_barnesPolynomial_two27
  have hb2val := integral_barnesPoleBlock27 2
  have hb3val := integral_barnesPoleBlock27 3
  have hb4val := integral_barnesPoleBlock27 4
  rw [barnesPoleTail_two27] at hb2val
  rw [barnesPoleTail_three27] at hb3val
  rw [barnesPoleTail_four27] at hb4val
  push_cast at hpval hb2val hb3val hb4val ⊢
  linear_combination hpval + hb2val - 54 * hb3val + 216 * hb4val

end RamanujanChallenge.P27
