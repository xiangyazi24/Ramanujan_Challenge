import RamanujanChallenge.Problem25DualCertificateDivergence
import RamanujanChallenge.Problem25DualCertificateRow1
import RamanujanChallenge.Problem25DualCertificateRow2

noncomputable section

namespace RamanujanChallenge.P25

open MeasureTheory Set Filter Topology

set_option maxRecDepth 10000

private theorem cube_ae_bounds' :
    ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      0 < x.1 ∧ x.1 ≤ 1 ∧
      0 < x.2.1 ∧ x.2.1 ≤ 1 ∧
      0 < x.2.2 ∧ x.2.2 ≤ 1 := by
  have hmem : ∀ᵐ x : ℝ × (ℝ × ℝ) ∂cubeMeasure,
      x ∈ Ioc (0 : ℝ) 1 ×ˢ (Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1) := by
    rw [Measure.ae_prod_mem_iff_ae_ae_mem
      (measurableSet_Ioc.prod (measurableSet_Ioc.prod measurableSet_Ioc))]
    filter_upwards [unit_ae_bounds] with p hp
    have hqv : ∀ᵐ y : ℝ × ℝ ∂unitMeasure.prod unitMeasure,
        y ∈ Ioc (0 : ℝ) 1 ×ˢ Ioc (0 : ℝ) 1 := by
      rw [Measure.ae_prod_mem_iff_ae_ae_mem
        (measurableSet_Ioc.prod measurableSet_Ioc)]
      filter_upwards [unit_ae_bounds] with q hq
      filter_upwards [unit_ae_bounds] with v hv
      exact ⟨hq, hv⟩
    filter_upwards [hqv] with y hy
    exact ⟨hp, hy⟩
  filter_upwards [hmem] with x hx
  exact ⟨hx.1.1, hx.1.2, hx.2.1.1, hx.2.1.2, hx.2.2.1, hx.2.2.2⟩

@[simp] theorem dualCertEval_NextNumPoly (n p q v : ℝ) (j : Fin 3) :
    dualCertEval (dualCertNextNumPoly j) n p q v =
      dualCertNextNum n j p q v := by
  fin_cases j <;>
    simp [dualCertNextNumPoly, dualCertNextNum, dualCertSnumPoly,
      dualCertDPoly, dualCertSnum, dualCertD, Matrix.cons_val_two]

@[simp] theorem dualCertEval_CurNumPoly (n p q v : ℝ) (i : Fin 3) :
    dualCertEval (dualCertCurNumPoly i) n p q v =
      dualCertCurNum n i p q v := by
  fin_cases i <;>
    simp [dualCertCurNumPoly, dualCertCurNum, dualCertDPoly, dualCertD,
      Matrix.cons_val_two]

private theorem dualCertD_eq_dualD (p q v : ℝ) :
    dualCertD p q v = dualD p q v := rfl

private theorem dualCertNext_weight_ae (n : ℕ) (j : Fin 3) :
    dualCertWeightedIntegrand n (dualCertNextNumPoly j) =ᵐ[cubeMeasure]
      ![rawMomentIntegrand (n + 1) 0 0 0 0,
        fun x => (2 * ((n : ℝ) + 3)) *
          rawMomentIntegrand (n + 1) 0 0 1 1 x,
        fun x => -((n : ℝ) + 3) *
            rawMomentIntegrand (n + 1) 0 0 1 1 x +
          (2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7)) *
            rawMomentIntegrand (n + 1) 0 0 2 2 x] j := by
  fin_cases j
  all_goals
    filter_upwards [cube_ae_bounds'] with x hx
    rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
    have hD : dualD x.1 x.2.1 x.2.2 ≠ 0 := by
      apply ne_of_gt
      unfold dualD
      positivity
    unfold dualCertWeightedIntegrand rawMomentIntegrand
    simp only [dualCertEval_NextNumPoly]
    simp [dualCertNextNum, Matrix.cons_val_two]
    simp only [dualCertSnum]
    rw [dualCertD_eq_dualD]
  all_goals
    rw [show 2 * (n + 1) + 6 = (2 * n + 6) + 2 by omega,
      show 2 * (n + 1) + 5 = (2 * n + 5) + 2 by omega,
      show 2 * (n + 1) + 3 = (2 * n + 3) + 2 by omega,
      show 2 * n + 4 + 4 = (2 * n + 6) + 2 by omega,
      pow_add, pow_succ]
    field_simp [hD]
    ring

private theorem dualCertCur_weight_ae (n : ℕ) (i : Fin 3) :
    dualCertWeightedIntegrand n (dualCertCurNumPoly i) =ᵐ[cubeMeasure]
      ![rawMomentIntegrand n 0 0 0 0,
        fun x => (2 * ((n : ℝ) + 2)) * rawMomentIntegrand n 0 0 1 1 x,
        fun x => -((n : ℝ) + 2) * rawMomentIntegrand n 0 0 1 1 x +
          (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5)) *
            rawMomentIntegrand n 0 0 2 2 x] i := by
  fin_cases i
  all_goals
    filter_upwards [cube_ae_bounds'] with x hx
    rcases hx with ⟨hp, hp1, hq, hq1, hv, hv1⟩
    have hD : dualD x.1 x.2.1 x.2.2 ≠ 0 := by
      apply ne_of_gt
      unfold dualD
      positivity
    unfold dualCertWeightedIntegrand rawMomentIntegrand
    simp only [dualCertEval_CurNumPoly]
    simp [dualCertCurNum, Matrix.cons_val_two]
    rw [dualCertD_eq_dualD]
  all_goals
    rw [show 2 * n + 4 + 4 = (2 * n + 4) + 4 by omega, pow_add]
    field_simp [hD]
    ring

theorem dualCertNext_integral (n : ℕ) (j : Fin 3) :
    (∫ x, dualCertWeightedIntegrand n (dualCertNextNumPoly j) x ∂cubeMeasure) =
      dualVector (n + 1) j := by
  rw [integral_congr_ae (dualCertNext_weight_ae n j)]
  fin_cases j
  · rfl
  · change (∫ x, (2 * ((n : ℝ) + 3)) *
        rawMomentIntegrand (n + 1) 0 0 1 1 x ∂cubeMeasure) =
      2 * (((n + 1 : ℕ) : ℝ) + 2) * dualMoment (n + 1) 0 0 1 1
    rw [MeasureTheory.integral_const_mul]
    change 2 * ((n : ℝ) + 3) * dualMoment (n + 1) 0 0 1 1 = _
    congr 2
    push_cast
    ring
  · change (∫ x, -((n : ℝ) + 3) *
          rawMomentIntegrand (n + 1) 0 0 1 1 x +
        (2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7)) *
          rawMomentIntegrand (n + 1) 0 0 2 2 x ∂cubeMeasure) =
      -(((n + 1 : ℕ) : ℝ) + 2) * dualMoment (n + 1) 0 0 1 1 +
        2 * (((n + 1 : ℕ) : ℝ) + 2) *
          (2 * ((n + 1 : ℕ) : ℝ) + 5) * dualMoment (n + 1) 0 0 2 2
    have h1 := rawMomentIntegrand_integrable (n + 1) 0 0 1 1 (by omega)
    have h2 := rawMomentIntegrand_integrable (n + 1) 0 0 2 2 (by omega)
    rw [MeasureTheory.integral_add (h1.const_mul _) (h2.const_mul _),
      MeasureTheory.integral_const_mul, MeasureTheory.integral_const_mul]
    change -((n : ℝ) + 3) * dualMoment (n + 1) 0 0 1 1 +
        2 * ((n : ℝ) + 3) * (2 * (n : ℝ) + 7) *
          dualMoment (n + 1) 0 0 2 2 = _
    push_cast
    ring

theorem dualCertCur_integral (n : ℕ) (i : Fin 3) :
    (∫ x, dualCertWeightedIntegrand n (dualCertCurNumPoly i) x ∂cubeMeasure) =
      dualVector n i := by
  rw [integral_congr_ae (dualCertCur_weight_ae n i)]
  fin_cases i
  · rfl
  · change (∫ x, (2 * ((n : ℝ) + 2)) *
        rawMomentIntegrand n 0 0 1 1 x ∂cubeMeasure) =
      2 * ((n : ℝ) + 2) * dualMoment n 0 0 1 1
    rw [MeasureTheory.integral_const_mul]
    rfl
  · change (∫ x, -((n : ℝ) + 2) * rawMomentIntegrand n 0 0 1 1 x +
        (2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5)) *
          rawMomentIntegrand n 0 0 2 2 x ∂cubeMeasure) =
      -((n : ℝ) + 2) * dualMoment n 0 0 1 1 +
        2 * ((n : ℝ) + 2) * (2 * (n : ℝ) + 5) * dualMoment n 0 0 2 2
    have h1 := rawMomentIntegrand_integrable n 0 0 1 1 (by omega)
    have h2 := rawMomentIntegrand_integrable n 0 0 2 2 (by omega)
    rw [MeasureTheory.integral_add (h1.const_mul _) (h2.const_mul _),
      MeasureTheory.integral_const_mul, MeasureTheory.integral_const_mul]
    rfl

private theorem dualCert_row_adjoint (n : ℕ) (i : Fin 3)
    (Pp Pq Pv : DualCertPoly)
    (hidentity : ∀ p q v : ℝ,
      (4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)) *
          ((∑ j : Fin 3, (positiveMatrix (n : ℤ) i j : ℝ) *
              dualCertEval (dualCertNextNumPoly j) (n : ℝ) p q v) -
            dualCertLambda (n : ℝ) *
              dualCertEval (dualCertCurNumPoly i) (n : ℝ) p q v) =
        dualCertEval (dualCertOpPPoly Pp) (n : ℝ) p q v +
          dualCertEval (dualCertOpQPoly Pq) (n : ℝ) p q v +
          dualCertEval (dualCertOpVPoly Pv) (n : ℝ) p q v) :
    ∑ j : Fin 3, (positiveMatrix (n : ℤ) i j : ℝ) * dualVector (n + 1) j =
      dualCertLambda (n : ℝ) * dualVector n i := by
  let δ : ℝ := 4 * (2 * (n : ℝ) + 3) * ((n : ℝ) + 2)
  let a0 : ℝ := positiveMatrix (n : ℤ) i 0
  let a1 : ℝ := positiveMatrix (n : ℤ) i 1
  let a2 : ℝ := positiveMatrix (n : ℤ) i 2
  let L : ℝ := dualCertLambda (n : ℝ)
  let f0 := dualCertWeightedIntegrand n (dualCertNextNumPoly 0)
  let f1 := dualCertWeightedIntegrand n (dualCertNextNumPoly 1)
  let f2 := dualCertWeightedIntegrand n (dualCertNextNumPoly 2)
  let fc := dualCertWeightedIntegrand n (dualCertCurNumPoly i)
  let fp := dualCertWeightedIntegrand n (dualCertOpPPoly Pp)
  let fq := dualCertWeightedIntegrand n (dualCertOpQPoly Pq)
  let fv := dualCertWeightedIntegrand n (dualCertOpVPoly Pv)
  have hf0 : Integrable f0 cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  have hf1 : Integrable f1 cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  have hf2 : Integrable f2 cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  have hfc : Integrable fc cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  have hfp : Integrable fp cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  have hfq : Integrable fq cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  have hfv : Integrable fv cubeMeasure := dualCertWeightedIntegrand_integrable _ _
  have hpoint : ∀ x : ℝ × (ℝ × ℝ),
      δ * (a0 * f0 x + a1 * f1 x + a2 * f2 x - L * fc x) =
        fp x + fq x + fv x := by
    intro x
    have h := hidentity x.1 x.2.1 x.2.2
    simp only [Fin.sum_univ_three] at h
    unfold δ a0 a1 a2 L f0 f1 f2 fc fp fq fv
    unfold dualCertWeightedIntegrand
    linear_combination (rawMomentIntegrand n 0 0 0 4 x) * h
  have hsum : Integrable (fun x => a0 * f0 x + a1 * f1 x + a2 * f2 x)
      cubeMeasure :=
    ((hf0.const_mul a0).add (hf1.const_mul a1)).add (hf2.const_mul a2)
  have hsum_integral :
      (∫ x, a0 * f0 x + a1 * f1 x + a2 * f2 x ∂cubeMeasure) =
        a0 * (∫ x, f0 x ∂cubeMeasure) +
          a1 * (∫ x, f1 x ∂cubeMeasure) +
          a2 * (∫ x, f2 x ∂cubeMeasure) := by
    calc
      _ = (∫ x, a0 * f0 x + a1 * f1 x ∂cubeMeasure) +
          ∫ x, a2 * f2 x ∂cubeMeasure := by
        exact MeasureTheory.integral_add
          ((hf0.const_mul a0).add (hf1.const_mul a1)) (hf2.const_mul a2)
      _ = _ := by
        rw [MeasureTheory.integral_add (hf0.const_mul a0) (hf1.const_mul a1),
          MeasureTheory.integral_const_mul, MeasureTheory.integral_const_mul,
          MeasureTheory.integral_const_mul]
  have hleft :
      (∫ x, δ * (a0 * f0 x + a1 * f1 x + a2 * f2 x - L * fc x)
          ∂cubeMeasure) =
        δ * ((a0 * dualVector (n + 1) 0 + a1 * dualVector (n + 1) 1 +
          a2 * dualVector (n + 1) 2) - L * dualVector n i) := by
    rw [MeasureTheory.integral_const_mul,
      MeasureTheory.integral_sub hsum (hfc.const_mul L), hsum_integral,
      MeasureTheory.integral_const_mul,
      dualCertNext_integral, dualCertNext_integral, dualCertNext_integral,
      dualCertCur_integral]
  have hzero :
      (∫ x, δ * (a0 * f0 x + a1 * f1 x + a2 * f2 x - L * fc x)
          ∂cubeMeasure) = 0 := by
    calc
      _ = ∫ x, fp x + fq x + fv x ∂cubeMeasure :=
        integral_congr_ae (Filter.Eventually.of_forall hpoint)
      _ = (∫ x, fp x ∂cubeMeasure) + (∫ x, fq x ∂cubeMeasure) +
          (∫ x, fv x ∂cubeMeasure) := by
        calc
          _ = (∫ x, fp x + fq x ∂cubeMeasure) +
              ∫ x, fv x ∂cubeMeasure :=
            MeasureTheory.integral_add (hfp.add hfq) hfv
          _ = _ := by rw [MeasureTheory.integral_add hfp hfq]
      _ = 0 := by
        rw [dualCertOpP_cube_integral_zero, dualCertOpQ_cube_integral_zero,
          dualCertOpV_cube_integral_zero]
        ring
  rw [hleft] at hzero
  have hδ : δ ≠ 0 := by
    unfold δ
    positivity
  have hinner :
      a0 * dualVector (n + 1) 0 + a1 * dualVector (n + 1) 1 +
        a2 * dualVector (n + 1) 2 - L * dualVector n i = 0 :=
    (mul_eq_zero.mp hzero).resolve_left hδ
  simp only [Fin.sum_univ_three]
  change a0 * dualVector (n + 1) 0 + a1 * dualVector (n + 1) 1 +
    a2 * dualVector (n + 1) 2 = L * dualVector n i
  linarith

theorem dualVector_adjoint (n : ℕ) (i : Fin 3) :
    ∑ j : Fin 3, (positiveMatrix (n : ℤ) i j : ℝ) * dualVector (n + 1) j =
      dualCertLambda (n : ℝ) * dualVector n i := by
  fin_cases i
  · apply dualCert_row_adjoint n 0 dualCertPp0Poly dualCertPq0Poly dualCertPv0Poly
    intro p q v
    simpa [Fin.sum_univ_three, dualCertPp0, dualCertPq0, dualCertPv0,
      dualCertDPp0, dualCertDPq0, dualCertDPv0] using
      dualCert_identity_zero n p q v
  · apply dualCert_row_adjoint n 1 dualCertPp1Poly dualCertPq1Poly dualCertPv1Poly
    intro p q v
    simpa [Fin.sum_univ_three, dualCertPp1, dualCertPq1, dualCertPv1,
      dualCertDPp1, dualCertDPq1, dualCertDPv1] using
      dualCert_identity_one n p q v
  · apply dualCert_row_adjoint n 2 dualCertPp2Poly dualCertPq2Poly dualCertPv2Poly
    intro p q v
    simpa [Fin.sum_univ_three, dualCertPp2, dualCertPq2, dualCertPv2,
      dualCertDPp2, dualCertDPq2, dualCertDPv2] using
      dualCert_identity_two n p q v

end RamanujanChallenge.P25

end
