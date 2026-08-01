import RamanujanChallenge.Problem27BarnesContour

/-!
# Problem 2.7: recurrence for the Barnes errors

This file iterates the proved one-strip contour shifts, moves every Barnes
integral to the common line `re t = 1/2`, and integrates the exact creative
telescoping identity.
-/

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27

def ctRPhiVerticalIntegral27 (n : ℕ) (x : ℝ) : ℂ :=
  ∫ y : ℝ, ctRKernelRaw27 n (verticalPoint27 x y)

def ctSPhiVerticalIntegral27 (n : ℕ) (x : ℝ) : ℂ :=
  ∫ y : ℝ, ctSKernelRaw27 n (verticalPoint27 x y)

theorem ctRPhiVerticalIntegral_one_strip27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    ctRPhiVerticalIntegral27 n ((m : ℝ) - 1 / 2) =
      ctRPhiVerticalIntegral27 n ((m : ℝ) + 1 / 2) := by
  simpa only [ctRPhiVerticalIntegral27] using
    (ctRKernelRaw_one_strip27 (n := n) (m := m) hm1 hmn)

theorem ctSPhiVerticalIntegral_one_strip27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) :
    ctSPhiVerticalIntegral27 n ((m : ℝ) - 1 / 2) =
      ctSPhiVerticalIntegral27 n ((m : ℝ) + 1 / 2) := by
  simpa only [ctSPhiVerticalIntegral27] using
    (ctSKernelRaw_one_strip27 (n := n) (m := m) hm1 hmn)

theorem ctRPhiVerticalIntegral_shift_down27
    (n m : ℕ) (hmn : m ≤ n) :
    ctRPhiVerticalIntegral27 n ((m : ℝ) + 1 / 2) =
      ctRPhiVerticalIntegral27 n (1 / 2) := by
  induction m with
  | zero => norm_num
  | succ m ih =>
      have hm' : m ≤ n := by omega
      have hstrip := ctRPhiVerticalIntegral_one_strip27
        (n := n) (m := m + 1) (by omega) hmn
      calc
        ctRPhiVerticalIntegral27 n
            (((m + 1 : ℕ) : ℝ) + 1 / 2) =
          ctRPhiVerticalIntegral27 n
            ((m : ℝ) + 1 / 2) := by
              convert hstrip.symm using 1 <;> push_cast <;> ring
        _ = ctRPhiVerticalIntegral27 n (1 / 2) := ih hm'

theorem ctRPhiVerticalIntegral_native_to_half27 (n : ℕ) :
    ctRPhiVerticalIntegral27 n ((n : ℝ) + 1 / 2) =
      ctRPhiVerticalIntegral27 n (1 / 2) :=
  ctRPhiVerticalIntegral_shift_down27 n n le_rfl

private theorem verticalPoint_native_eq_translate27
    (n : ℕ) (y : ℝ) :
    verticalPoint27 ((n : ℝ) + 1 / 2) y =
      zudilinBarnesLine27 y + (((n + 1 : ℕ) : ℂ)) := by
  unfold verticalPoint27 zudilinBarnesLine27
  push_cast
  ring

theorem zudilinBarnesFixedLineRaw_eq_ctRPhiNative27 (n : ℕ) :
    (∫ y : ℝ,
      zudilinBarnesPhi27 n (zudilinBarnesLine27 y) *
        zudilinBarnesSquaredSineKernel27
          (zudilinBarnesLine27 y)) =
      ctRPhiVerticalIntegral27 n ((n : ℝ) + 1 / 2) := by
  unfold ctRPhiVerticalIntegral27 ctRKernelRaw27
  apply integral_congr_ae
  filter_upwards with y
  rw [zudilinBarnesPhi_eq_ctRPhi_translate27]
  have hp := verticalPoint_native_eq_translate27 n y
  rw [hp, squaredSineKernel_add_nat27]

theorem zudilinBarnesErrorIntegral_eq_ctRPhiNative27 (n : ℕ) :
    zudilinBarnesErrorIntegral27 n =
      (1 / (2 * (Real.pi : ℂ))) *
        ctRPhiVerticalIntegral27 n ((n : ℝ) + 1 / 2) := by
  rw [zudilinBarnesErrorIntegral_eq_fixedLine27]
  congr 1
  exact zudilinBarnesFixedLineRaw_eq_ctRPhiNative27 n

theorem zudilinBarnesErrorIntegral_eq_ctRPhiHalf27 (n : ℕ) :
    zudilinBarnesErrorIntegral27 n =
      (1 / (2 * (Real.pi : ℂ))) *
        ctRPhiVerticalIntegral27 n (1 / 2) := by
  rw [zudilinBarnesErrorIntegral_eq_ctRPhiNative27,
    ctRPhiVerticalIntegral_native_to_half27]

theorem ctSPhiVerticalIntegral_threeHalf_eq_half27 (k : ℕ) :
    ctSPhiVerticalIntegral27 (k + 2) (3 / 2) =
      ctSPhiVerticalIntegral27 (k + 2) (1 / 2) := by
  have h := (ctSPhiVerticalIntegral_one_strip27
    (n := k + 2) (m := 1) (by norm_num) (by omega)).symm
  convert h using 1 <;> norm_num

theorem ctSPhi_certificate_integral_shift27 (k : ℕ) :
    (∫ y : ℝ,
      ctSPhi27 (k + 2) (verticalPoint27 (1 / 2) y + 1) *
        zudilinBarnesSquaredSineKernel27
          (verticalPoint27 (1 / 2) y)) =
      ∫ y : ℝ,
        ctSPhi27 (k + 2) (verticalPoint27 (1 / 2) y) *
          zudilinBarnesSquaredSineKernel27
            (verticalPoint27 (1 / 2) y) := by
  calc
    _ = ctSPhiVerticalIntegral27 (k + 2) (3 / 2) := by
      unfold ctSPhiVerticalIntegral27 ctSKernelRaw27
      apply integral_congr_ae
      filter_upwards with y
      have hp : verticalPoint27 (1 / 2) y + 1 =
          verticalPoint27 (3 / 2) y := by
        unfold verticalPoint27
        push_cast
        ring
      rw [← squaredSineKernel_add_one27
        (verticalPoint27 (1 / 2) y), hp]
    _ = ctSPhiVerticalIntegral27 (k + 2) (1 / 2) :=
      ctSPhiVerticalIntegral_threeHalf_eq_half27 k
    _ = _ := by rfl

theorem integrable_ctRKernelRaw_half27 (n : ℕ) :
    Integrable (fun y : ℝ =>
      ctRKernelRaw27 n (verticalPoint27 (1 / 2) y)) := by
  convert integrable_ctRKernelRaw_left27 n 1 (by norm_num) using 1 <;>
    norm_num

theorem integrable_ctSKernelRaw_half27 (n : ℕ) :
    Integrable (fun y : ℝ =>
      ctSKernelRaw27 n (verticalPoint27 (1 / 2) y)) := by
  convert integrable_ctSKernelRaw_left27 n 1 (by norm_num) using 1 <;>
    norm_num

theorem integrable_ctSPhi_certificate_shift27 (k : ℕ) :
    Integrable (fun y : ℝ =>
      ctSPhi27 (k + 2) (verticalPoint27 (1 / 2) y + 1) *
        zudilinBarnesSquaredSineKernel27
          (verticalPoint27 (1 / 2) y)) := by
  apply (integrable_ctSKernelRaw_right27 (k + 2) 1).congr
  filter_upwards with y
  have hp : verticalPoint27 (1 / 2) y + 1 =
      verticalPoint27 (3 / 2) y := by
    unfold verticalPoint27
    push_cast
    ring
  unfold ctSKernelRaw27
  rw [← squaredSineKernel_add_one27 (verticalPoint27 (1 / 2) y), hp]
  norm_num

theorem ctRPhiVerticalIntegral_telescoper27 (k : ℕ) :
    (ctAlpha27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 3) (1 / 2)
      - (ctBeta27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 2) (1 / 2)
      + (ctGamma27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 1) (1 / 2)
      - (ctDelta27 k : ℂ) * ctRPhiVerticalIntegral27 k (1 / 2) = 0 := by
  let t : ℝ → ℂ := fun y => verticalPoint27 (1 / 2) y
  let ker : ℝ → ℂ := fun y =>
    zudilinBarnesSquaredSineKernel27 (t y)
  let r : ℕ → ℝ → ℂ := fun n y => ctRKernelRaw27 n (t y)
  let s1 : ℝ → ℂ := fun y => ctSPhi27 (k + 2) (t y + 1) * ker y
  let s0 : ℝ → ℂ := fun y => ctSKernelRaw27 (k + 2) (t y)
  have hpoint : ∀ y : ℝ,
      (ctAlpha27 k : ℂ) * r (k + 3) y
        - (ctBeta27 k : ℂ) * r (k + 2) y
        + (ctGamma27 k : ℂ) * r (k + 1) y
        - (ctDelta27 k : ℂ) * r k y = s1 y - s0 y := by
    intro y
    have htmem : t y ∈ ctClosedStrip27 1 := by
      dsimp only [t]
      norm_num [ctClosedStrip27, verticalPoint27]
    have hpole : ctPoleProduct27 (k + 4) (t y) ≠ 0 :=
      ctPoleProduct_ne_zero_on_strip27 (m := 1) (M := k + 4)
        (by norm_num) htmem
    have htel := ctPhi_telescoper_step27 k (t y) hpole
    dsimp only [r, s1, s0, ker]
    unfold ctRKernelRaw27 ctSKernelRaw27
    linear_combination htel *
      zudilinBarnesSquaredSineKernel27 (t y)
  have h3 : Integrable (r (k + 3)) := by
    simpa only [r, t] using integrable_ctRKernelRaw_half27 (k + 3)
  have h2 : Integrable (r (k + 2)) := by
    simpa only [r, t] using integrable_ctRKernelRaw_half27 (k + 2)
  have h1 : Integrable (r (k + 1)) := by
    simpa only [r, t] using integrable_ctRKernelRaw_half27 (k + 1)
  have h0 : Integrable (r k) := by
    simpa only [r, t] using integrable_ctRKernelRaw_half27 k
  have hs1 : Integrable s1 := by
    simpa only [s1, ker, t] using integrable_ctSPhi_certificate_shift27 k
  have hs0 : Integrable s0 := by
    simpa only [s0, t] using integrable_ctSKernelRaw_half27 (k + 2)
  have ha3 := h3.const_mul (ctAlpha27 k : ℂ)
  have ha2 := h2.const_mul (ctBeta27 k : ℂ)
  have ha1 := h1.const_mul (ctGamma27 k : ℂ)
  have ha0 := h0.const_mul (ctDelta27 k : ℂ)
  have h32 := ha3.sub ha2
  have h321 := h32.add ha1
  have hint :
      (∫ y : ℝ,
        (ctAlpha27 k : ℂ) * r (k + 3) y
          - (ctBeta27 k : ℂ) * r (k + 2) y
          + (ctGamma27 k : ℂ) * r (k + 1) y
          - (ctDelta27 k : ℂ) * r k y) =
        ∫ y : ℝ, s1 y - s0 y := by
    apply integral_congr_ae
    filter_upwards with y
    exact hpoint y
  have hlinR :
      (∫ y : ℝ,
        (ctAlpha27 k : ℂ) * r (k + 3) y
          - (ctBeta27 k : ℂ) * r (k + 2) y
          + (ctGamma27 k : ℂ) * r (k + 1) y
          - (ctDelta27 k : ℂ) * r k y) =
        (ctAlpha27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 3) (1 / 2)
          - (ctBeta27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 2) (1 / 2)
          + (ctGamma27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 1) (1 / 2)
          - (ctDelta27 k : ℂ) * ctRPhiVerticalIntegral27 k (1 / 2) := by
    have eouter := MeasureTheory.integral_sub h321 ha0
    have eadd := MeasureTheory.integral_add h32 ha1
    have esub := MeasureTheory.integral_sub ha3 ha2
    simp only [Pi.sub_apply, Pi.add_apply] at eouter eadd esub
    have e3 := MeasureTheory.integral_const_mul (μ := volume)
      (ctAlpha27 k : ℂ) (r (k + 3))
    have e2 := MeasureTheory.integral_const_mul (μ := volume)
      (ctBeta27 k : ℂ) (r (k + 2))
    have e1 := MeasureTheory.integral_const_mul (μ := volume)
      (ctGamma27 k : ℂ) (r (k + 1))
    have e0 := MeasureTheory.integral_const_mul (μ := volume)
      (ctDelta27 k : ℂ) (r k)
    rw [eouter, eadd, esub]
    calc
      _ =
          (((ctAlpha27 k : ℂ) * (∫ y : ℝ, r (k + 3) y) -
              (ctBeta27 k : ℂ) * (∫ y : ℝ, r (k + 2) y)) +
            (ctGamma27 k : ℂ) * (∫ y : ℝ, r (k + 1) y)) -
            (ctDelta27 k : ℂ) * (∫ y : ℝ, r k y) := by
              exact congrArg₂ (· - ·)
                (congrArg₂ (· + ·) (congrArg₂ (· - ·) e3 e2) e1) e0
      _ = _ := by rfl
  have hlinS :
      (∫ y : ℝ, s1 y - s0 y) =
        (∫ y : ℝ, s1 y) - ∫ y : ℝ, s0 y :=
    MeasureTheory.integral_sub hs1 hs0
  have hcert : (∫ y : ℝ, s1 y) = ∫ y : ℝ, s0 y := by
    simpa only [s1, s0, ker, t, ctSKernelRaw27] using
      ctSPhi_certificate_integral_shift27 k
  calc
    (ctAlpha27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 3) (1 / 2)
        - (ctBeta27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 2) (1 / 2)
        + (ctGamma27 k : ℂ) * ctRPhiVerticalIntegral27 (k + 1) (1 / 2)
        - (ctDelta27 k : ℂ) * ctRPhiVerticalIntegral27 k (1 / 2) =
      ∫ y : ℝ,
        (ctAlpha27 k : ℂ) * r (k + 3) y
          - (ctBeta27 k : ℂ) * r (k + 2) y
          + (ctGamma27 k : ℂ) * r (k + 1) y
          - (ctDelta27 k : ℂ) * r k y := hlinR.symm
    _ = ∫ y : ℝ, s1 y - s0 y := hint
    _ = (∫ y : ℝ, s1 y) - ∫ y : ℝ, s0 y := hlinS
    _ = 0 := by rw [hcert, sub_self]

def ZudilinSatisfiesRecC27 (u : ℕ → ℂ) : Prop :=
  ∀ k : ℕ,
    (ctAlpha27 k : ℂ) * u (k + 3)
      - (ctBeta27 k : ℂ) * u (k + 2)
      + (ctGamma27 k : ℂ) * u (k + 1)
      - (ctDelta27 k : ℂ) * u k = 0

theorem ctRPhiVerticalIntegral_recC27 :
    ZudilinSatisfiesRecC27
      (fun n => ctRPhiVerticalIntegral27 n (1 / 2)) :=
  ctRPhiVerticalIntegral_telescoper27

theorem zudilinBarnesErrorIntegral_recC27 :
    ZudilinSatisfiesRecC27 zudilinBarnesErrorIntegral27 := by
  intro k
  rw [zudilinBarnesErrorIntegral_eq_ctRPhiHalf27,
    zudilinBarnesErrorIntegral_eq_ctRPhiHalf27,
    zudilinBarnesErrorIntegral_eq_ctRPhiHalf27,
    zudilinBarnesErrorIntegral_eq_ctRPhiHalf27]
  have h := ctRPhiVerticalIntegral_telescoper27 k
  calc
    (ctAlpha27 k : ℂ) *
          ((1 / (2 * (Real.pi : ℂ))) *
            ctRPhiVerticalIntegral27 (k + 3) (1 / 2))
        - (ctBeta27 k : ℂ) *
          ((1 / (2 * (Real.pi : ℂ))) *
            ctRPhiVerticalIntegral27 (k + 2) (1 / 2))
        + (ctGamma27 k : ℂ) *
          ((1 / (2 * (Real.pi : ℂ))) *
            ctRPhiVerticalIntegral27 (k + 1) (1 / 2))
        - (ctDelta27 k : ℂ) *
          ((1 / (2 * (Real.pi : ℂ))) *
            ctRPhiVerticalIntegral27 k (1 / 2)) =
      (1 / (2 * (Real.pi : ℂ))) *
        ((ctAlpha27 k : ℂ) *
            ctRPhiVerticalIntegral27 (k + 3) (1 / 2)
          - (ctBeta27 k : ℂ) *
            ctRPhiVerticalIntegral27 (k + 2) (1 / 2)
          + (ctGamma27 k : ℂ) *
            ctRPhiVerticalIntegral27 (k + 1) (1 / 2)
          - (ctDelta27 k : ℂ) *
            ctRPhiVerticalIntegral27 k (1 / 2)) := by ring
    _ = 0 := by rw [h, mul_zero]

end RamanujanChallenge.P27
