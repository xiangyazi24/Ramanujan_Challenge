import RamanujanChallenge.Q6427HorizontalFinal
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic

open Filter Set Topology MeasureTheory
open scoped Interval Real BigOperators

noncomputable section

namespace RamanujanChallenge.P27.Q6427

/-! ## Integrability of the polynomial-exponential majorant -/

theorem integrableOn_pow_mul_exp_neg_mul27_final
    (p : ℕ) {a : ℝ} (ha : 0 < a) :
    IntegrableOn (fun t : ℝ => t ^ p * Real.exp (-(a * t))) (Ioi 0) := by
  have hbase :
      IntegrableOn (fun x : ℝ => Real.exp (-x) * x ^ p) (Ioi 0) := by
    convert Real.GammaIntegral_convergent
      (s := (p : ℝ) + 1) (by positivity) using 1
    ext x
    rw [show (p : ℝ) + 1 - 1 = p by ring, Real.rpow_natCast]
  have hbase0 :
      IntegrableOn (fun x : ℝ => Real.exp (-x) * x ^ p) (Ioi (a * 0)) := by
    simpa only [mul_zero] using hbase
  have hscaled :
      IntegrableOn
        (fun t : ℝ => Real.exp (-(a * t)) * (a * t) ^ p) (Ioi 0) := by
    exact (integrableOn_Ioi_comp_mul_left_iff
      (fun x : ℝ => Real.exp (-x) * x ^ p) 0 ha).2 hbase0
  have hc := hscaled.const_mul ((a ^ p)⁻¹)
  refine hc.congr ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  rw [mul_pow]
  field_simp [ha.ne']

theorem integrableOn_ctHorizontalMajorant27_final (n : ℕ) :
    IntegrableOn (ctHorizontalMajorant27Final n) (Ioi 0) := by
  let C : ℝ :=
    ((2 * ((n : ℝ) + 2)) ^ (3 * n) /
      (1 / 2 : ℝ) ^ (n + 1)) * (16 * Real.pi ^ 2)
  have hbase := integrableOn_pow_mul_exp_neg_mul27_final
    (3 * n) Real.two_pi_pos
  have hc := hbase.const_mul C
  refine hc.congr ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  unfold ctHorizontalMajorant27Final
  have hE : Real.exp (Real.pi * t) ^ 2 = Real.exp (2 * Real.pi * t) := by
    rw [pow_two, ← Real.exp_add]
    congr 1
    ring
  rw [hE]
  simp only [div_eq_mul_inv]
  rw [← Real.exp_neg]
  dsimp [C]
  ring

def ctAbsMajorant27Final (n : ℕ) (y : ℝ) : ℝ :=
  ctHorizontalMajorant27Final n |y|

theorem integrable_ctAbsMajorant27_final (n : ℕ) :
    Integrable (ctAbsMajorant27Final n) := by
  have hpos0 := integrableOn_ctHorizontalMajorant27_final n
  have hpos : IntegrableOn (ctAbsMajorant27Final n) (Ioi 0) := by
    refine hpos0.congr ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with y hy
    have hy0 : 0 < y := hy
    simp [ctAbsMajorant27Final, abs_of_pos hy0]
  have hneg : IntegrableOn (ctAbsMajorant27Final n) (Iic 0) := by
    rw [← Measure.map_neg_eq_self (volume : Measure ℝ)]
    let e : MeasurableEmbedding (fun x : ℝ => -x) :=
      (Homeomorph.neg ℝ).measurableEmbedding
    rw [e.integrableOn_map_iff]
    simp_rw [Function.comp_def, ctAbsMajorant27Final, abs_neg,
      neg_preimage, neg_Iic, neg_zero]
    exact Iff.mpr integrableOn_Ici_iff_integrableOn_Ioi hpos
  have hu : IntegrableOn (ctAbsMajorant27Final n) (Iic 0 ∪ Ioi 0) :=
    hneg.union hpos
  simpa [Iic_union_Ioi, IntegrableOn] using hu

/-! ## Continuity and vertical integrability of the removable extension -/

theorem ctExtension_differentiableAt27_final
    {n m : ℕ} (hm1 : 1 ≤ m) {t : ℂ}
    (ht : t ∈ halfIntegerStrip (m : ℤ)) :
    DifferentiableAt ℂ (ctExtension27 n m) t := by
  have hremoved : DifferentiableAt ℂ (ctRemoved27 n m) t := by
    unfold ctRemoved27
    apply DifferentiableAt.div
    · exact ctNumeratorWithout_differentiableAt27 n m t
    · exact (differentiableAt_const _).mul
        (ctPoleProduct_differentiableAt27 (n + 1) t)
    · apply mul_ne_zero
      · exact pow_ne_zero _ (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n))
      · exact ctPoleProduct_ne_zero_on_strip27 hm1 ht
  have hnum : DifferentiableAt ℂ
      (fun z => (Real.pi : ℂ) ^ 2 * (z - (m : ℂ)) * ctRemoved27 n m z) t :=
    ((differentiableAt_const _).mul
      (differentiableAt_id.sub_const (m : ℂ))).mul hremoved
  have hden : DifferentiableAt ℂ (fun z => sineSlope (m : ℤ) z ^ 2) t :=
    (sineSlope_differentiable (m : ℤ) t).pow 2
  simpa [ctExtension27] using hnum.div hden
    (pow_ne_zero 2 (sineSlope_ne_zero_on_strip (m : ℤ) ht))

theorem continuous_ctExtension_vertical27_final
    {n m : ℕ} (hm1 : 1 ≤ m) {x : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2)) :
    Continuous (fun y : ℝ => ctExtension27 n m (verticalPoint x y)) := by
  rw [continuous_iff_continuousAt]
  intro y
  have ht : verticalPoint x y ∈ halfIntegerStrip (m : ℤ) := by
    simpa [verticalPoint, halfIntegerStrip] using hx
  exact (ctExtension_differentiableAt27_final hm1 ht).continuousAt.comp
    (by fun_prop : ContinuousAt (verticalPoint x) y)

theorem integrable_ctExtension_vertical27_final
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) {x : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hxm : x ≠ (m : ℝ)) :
    Integrable (fun y : ℝ => ctExtension27 n m (verticalPoint x y)) := by
  let f : ℝ → ℂ := fun y => ctExtension27 n m (verticalPoint x y)
  let g : ℝ → ℝ := ctAbsMajorant27Final n
  have hg : Integrable g := integrable_ctAbsMajorant27_final n
  have hfcont : Continuous f := by
    simpa [f] using continuous_ctExtension_vertical27_final hm1 hx
  have hcentral : IntegrableOn f (Set.Icc (-1 : ℝ) 1) :=
    hfcont.continuousOn.integrableOn_compact isCompact_Icc
  have hpos : IntegrableOn f (Ioi 1) := by
    have hgpos : IntegrableOn g (Ioi 1) := hg.integrableOn
    apply hgpos.mono'
    · exact hfcont.aestronglyMeasurable.restrict
    · filter_upwards [ae_restrict_mem measurableSet_Ioi] with y hy
      have hygt : 1 < y := hy
      have hy1 : 1 ≤ y := hygt.le
      have hy0 : 0 < y := zero_lt_one.trans hygt
      have hmem : verticalPoint x y ∈ halfIntegerStrip (m : ℤ) := by
        simpa [verticalPoint, halfIntegerStrip] using hx
      have hne : verticalPoint x y ≠ (m : ℂ) := by
        intro h
        have hre := congrArg Complex.re h
        simp [verticalPoint] at hre
        exact hxm hre
      rw [show f y = ctIntegrand27 n (verticalPoint x y) by
        dsimp [f]
        symm
        exact ctIntegrand_eq_extension27 hm1 hmn hmem hne]
      have hbound := norm_ctIntegrand_top_le_final hm1 hmn hx hy1
      simpa [g, ctAbsMajorant27Final, verticalPoint, abs_of_pos hy0] using hbound
  have hneg : IntegrableOn f (Iio (-1)) := by
    have hgneg : IntegrableOn g (Iio (-1)) := hg.integrableOn
    apply hgneg.mono'
    · exact hfcont.aestronglyMeasurable.restrict
    · filter_upwards [ae_restrict_mem measurableSet_Iio] with y hy
      have hylt : y < -1 := hy
      have hy0 : y < 0 := hylt.trans (by norm_num)
      have hT : 1 ≤ -y := by linarith
      have hmem : verticalPoint x y ∈ halfIntegerStrip (m : ℤ) := by
        simpa [verticalPoint, halfIntegerStrip] using hx
      have hne : verticalPoint x y ≠ (m : ℂ) := by
        intro h
        have hre := congrArg Complex.re h
        simp [verticalPoint] at hre
        exact hxm hre
      rw [show f y = ctIntegrand27 n (verticalPoint x y) by
        dsimp [f]
        symm
        exact ctIntegrand_eq_extension27 hm1 hmn hmem hne]
      have hpoint : verticalPoint x y =
          (x : ℂ) - ((-y : ℝ) : ℂ) * Complex.I := by
        apply Complex.ext <;> simp [verticalPoint]
      rw [hpoint]
      have hbound := norm_ctIntegrand_bottom_le_final hm1 hmn hx hT
      simpa [g, ctAbsMajorant27Final, abs_of_neg hy0] using hbound
  have hall : IntegrableOn f
      (Iio (-1 : ℝ) ∪ Set.Icc (-1 : ℝ) 1 ∪ Ioi 1) :=
    (hneg.union hcentral).union hpos
  have hcover : Iio (-1 : ℝ) ∪ Set.Icc (-1 : ℝ) 1 ∪ Ioi 1 = Set.univ := by
    ext y
    simp only [Set.mem_union, Set.mem_Iio, Set.mem_Icc, Set.mem_Ioi,
      Set.mem_univ, iff_true]
    by_cases hlow : y < -1
    · exact Or.inl hlow
    by_cases hhigh : 1 < y
    · exact Or.inr (Or.inr hhigh)
    · exact Or.inr (Or.inl ⟨le_of_not_gt hlow, le_of_not_gt hhigh⟩)
  rw [hcover] at hall
  simpa [IntegrableOn] using hall

theorem integrable_ctIntegrand_vertical27_final
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) {x : ℝ}
    (hx : x ∈ Set.Icc ((m : ℝ) - 1 / 2) ((m : ℝ) + 1 / 2))
    (hxm : x ≠ (m : ℝ)) :
    Integrable (fun y : ℝ => ctIntegrand27 n (verticalPoint x y)) := by
  have hext := integrable_ctExtension_vertical27_final hm1 hmn hx hxm
  apply hext.congr
  filter_upwards with y
  symm
  apply ctIntegrand_eq_extension27 hm1 hmn
  · simpa [verticalPoint, halfIntegerStrip] using hx
  · intro h
    have hre := congrArg Complex.re h
    simp [verticalPoint] at hre
    exact hxm hre

/-! ## Complete one-strip and repeated-shift theorems -/

def ctVerticalIntegral27Final (n : ℕ) (x : ℝ) : ℂ :=
  ∫ y : ℝ, ctIntegrand27 n (verticalPoint x y)

theorem ctR_one_strip_shift_complete27_final
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    ctVerticalIntegral27Final n ((m : ℝ) - 1 / 2) =
      ctVerticalIntegral27Final n ((m : ℝ) + 1 / 2) := by
  unfold ctVerticalIntegral27Final
  apply ctR_one_strip_shift27 hm1 hmn
  · apply integrable_ctIntegrand_vertical27_final hm1 hmn
    · constructor <;> linarith
    · linarith
  · apply integrable_ctIntegrand_vertical27_final hm1 hmn
    · constructor <;> linarith
    · linarith
  · exact ctHorizontal_top_tendsto_zero27_final hm1 hmn
  · exact ctHorizontal_bottom_tendsto_zero27_final hm1 hmn

theorem ctR_repeated_shift27_final
    (n a d : ℕ) (had : a + d ≤ n) :
    ctVerticalIntegral27Final n ((a : ℝ) + 1 / 2) =
      ctVerticalIntegral27Final n (((a + d : ℕ) : ℝ) + 1 / 2) := by
  revert had
  induction d with
  | zero => simp
  | succ d ih =>
      intro had
      have had' : a + d ≤ n := by omega
      calc
        ctVerticalIntegral27Final n ((a : ℝ) + 1 / 2) =
            ctVerticalIntegral27Final n (((a + d : ℕ) : ℝ) + 1 / 2) :=
          ih had'
        _ = ctVerticalIntegral27Final n (((a + d + 1 : ℕ) : ℝ) + 1 / 2) := by
          have hone := ctR_one_strip_shift_complete27_final
            (n := n) (m := a + d + 1) (by omega) (by omega)
          convert hone using 1 <;> push_cast <;> ring_nf
        _ = ctVerticalIntegral27Final n (((a + (d + 1) : ℕ) : ℝ) + 1 / 2) := by
          simp [Nat.add_assoc]

theorem ctR_shift_to_native27_final
    {n a : ℕ} (ha : a ≤ n) :
    ctVerticalIntegral27Final n ((a : ℝ) + 1 / 2) =
      ctVerticalIntegral27Final n ((n : ℝ) + 1 / 2) := by
  have h := ctR_repeated_shift27_final n a (n - a) (by omega)
  simpa [Nat.add_sub_of_le ha] using h

#print axioms ctR_one_strip_shift_complete27_final
#print axioms ctR_repeated_shift27_final
#print axioms ctR_shift_to_native27_final

end RamanujanChallenge.P27.Q6427
