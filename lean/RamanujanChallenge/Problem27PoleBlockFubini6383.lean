import RamanujanChallenge.Problem27PoleBlockFourier6383

open Filter Set MeasureTheory Topology
open scoped Topology

noncomputable section

namespace RamanujanChallenge.P27.Q6383

def doubleKernel (m : ℕ) (t y : ℝ) : ℂ :=
  (laplaceWeight m t : ℂ) * fourierKernel t y

theorem integrable_doubleKernel (m : ℕ) :
    Integrable (Function.uncurry (doubleKernel m))
      ((volume.restrict (Ioi 0)).prod volume) := by
  have hmajor := (integrableOn_laplaceWeight m).norm.mul_prod integrable_sechSq
  apply hmajor.mono'
  · have hlap : Continuous fun p : ℝ × ℝ =>
        (laplaceWeight m p.1 : ℂ) := by
      unfold laplaceWeight
      fun_prop
    have hphase : Continuous fun p : ℝ × ℝ =>
        Complex.exp (-(p.2 * p.1 : ℝ) * Complex.I) := by
      fun_prop
    have hden : Continuous fun p : ℝ × ℝ =>
        (Real.cosh (Real.pi * p.2) : ℂ) ^ 2 := by
      fun_prop
    have hden0 : ∀ p : ℝ × ℝ,
        (Real.cosh (Real.pi * p.2) : ℂ) ^ 2 ≠ 0 := fun p =>
      pow_ne_zero 2 (Complex.ofReal_ne_zero.mpr (Real.cosh_pos _).ne')
    have hcont : Continuous fun p : ℝ × ℝ =>
        (laplaceWeight m p.1 : ℂ) *
          (Complex.exp (-(p.2 * p.1 : ℝ) * Complex.I) /
            (Real.cosh (Real.pi * p.2) : ℂ) ^ 2) :=
      hlap.mul (hphase.div hden hden0)
    change AEStronglyMeasurable
      (fun p : ℝ × ℝ =>
        (laplaceWeight m p.1 : ℂ) *
          (Complex.exp (-(p.2 * p.1 : ℝ) * Complex.I) /
            (Real.cosh (Real.pi * p.2) : ℂ) ^ 2))
      ((volume.restrict (Ioi 0)).prod volume)
    exact hcont.aestronglyMeasurable
  · filter_upwards with p
    rcases p with ⟨t, y⟩
    change ‖(laplaceWeight m t : ℂ) * fourierKernel t y‖ ≤
      ‖laplaceWeight m t‖ * sechSq y
    rw [norm_mul, Complex.norm_real]
    apply mul_le_mul_of_nonneg_left _ (norm_nonneg _)
    unfold fourierKernel sechSq
    rw [norm_div, Complex.norm_exp, norm_pow, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos (Real.cosh_pos _)]
    simp

theorem integral_doubleKernel_t (m : ℕ) (y : ℝ) :
    (∫ t : ℝ in Ioi 0, doubleKernel m t y) =
      poleBlock m y / (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
  calc
    (∫ t : ℝ in Ioi 0, doubleKernel m t y) =
        ∫ t : ℝ in Ioi 0,
          (((1 + t / 2 : ℝ) : ℂ) *
            Complex.exp (-polePoint m y * (t : ℂ))) /
              (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
          apply integral_congr_ae
          filter_upwards with t
          unfold doubleKernel laplaceWeight fourierKernel polePoint
          push_cast
          have hexp :
              Complex.exp (-(↑(poleAbscissa m) * (t : ℂ))) *
                  Complex.exp (-((y : ℂ) * (t : ℂ)) * Complex.I) =
                Complex.exp (-(↑(poleAbscissa m) + (y : ℂ) * Complex.I) *
                  (t : ℂ)) := by
            rw [← Complex.exp_add]
            congr 2
            ring
          calc
            (1 + (t : ℂ) / 2) *
                  Complex.exp (-(↑(poleAbscissa m) * (t : ℂ))) *
                    (Complex.exp (-((y : ℂ) * (t : ℂ)) * Complex.I) /
                      Complex.cosh ((Real.pi : ℂ) * (y : ℂ)) ^ 2) =
                (1 + (t : ℂ) / 2) *
                  (Complex.exp (-(↑(poleAbscissa m) * (t : ℂ))) *
                    Complex.exp (-((y : ℂ) * (t : ℂ)) * Complex.I)) /
                      Complex.cosh ((Real.pi : ℂ) * (y : ℂ)) ^ 2 := by ring
            _ = (1 + (t : ℂ) / 2) *
                  Complex.exp (-(↑(poleAbscissa m) + (y : ℂ) * Complex.I) *
                    (t : ℂ)) /
                    Complex.cosh ((Real.pi : ℂ) * (y : ℂ)) ^ 2 := by rw [hexp]
    _ = (∫ t : ℝ in Ioi 0,
          ((1 + t / 2 : ℝ) : ℂ) *
            Complex.exp (-polePoint m y * (t : ℂ))) /
              (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
          exact MeasureTheory.integral_div _ _
    _ = poleBlock m y /
          (Real.cosh (Real.pi * y) : ℂ) ^ 2 := by
          rw [integral_laplace_poleBlock]
          · rfl
          · simpa using poleAbscissa_pos m

theorem integral_doubleKernel_y (m : ℕ) {t : ℝ} (ht : 0 < t) :
    (∫ y : ℝ, doubleKernel m t y) =
      (laplaceWeight m t : ℂ) *
        ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
  calc
    (∫ y : ℝ, doubleKernel m t y) =
        (laplaceWeight m t : ℂ) * ∫ y : ℝ, fourierKernel t y := by
          exact MeasureTheory.integral_const_mul _ _
    _ = _ := by rw [integral_fourierKernel ht]

theorem integral_poleBlock_eq_fourierLaplace (m : ℕ) :
    (∫ y : ℝ,
      poleBlock m y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
      ∫ t : ℝ in Ioi 0,
        (laplaceWeight m t : ℂ) *
          ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
  have hswap := MeasureTheory.integral_integral_swap
    (integrable_doubleKernel m)
  calc
    (∫ y : ℝ,
        poleBlock m y / (Real.cosh (Real.pi * y) : ℂ) ^ 2) =
        ∫ y : ℝ, ∫ t : ℝ in Ioi 0, doubleKernel m t y := by
          apply integral_congr_ae
          filter_upwards with y
          exact (integral_doubleKernel_t m y).symm
    _ = ∫ t : ℝ in Ioi 0, ∫ y : ℝ, doubleKernel m t y := hswap.symm
    _ = ∫ t : ℝ in Ioi 0,
        (laplaceWeight m t : ℂ) *
          ((t / (Real.pi * Real.sinh (t / 2)) : ℝ) : ℂ) := by
          apply integral_congr_ae
          filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
          exact integral_doubleKernel_y m ht

#print axioms integral_poleBlock_eq_fourierLaplace

end RamanujanChallenge.P27.Q6383
