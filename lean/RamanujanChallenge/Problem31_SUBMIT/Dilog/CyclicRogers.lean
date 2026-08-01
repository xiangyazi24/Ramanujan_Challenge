import Ramanujan31.Dilog.FiveTerm

/-!
# Cyclic sine cross-ratios and Rogers five-term relations

For four consecutive positive angular gaps summing to `π`, `sineCross` is
the cross-ratio of the corresponding four cyclically ordered points.  This
file proves its positivity, its complement symmetry, and the five-point
relation used by the finite `17`-gon endpoint certificate.
-/

namespace Real

/-- The positive cross-ratio associated with four cyclic angular gaps. -/
noncomputable def sineCross (A B C D : ℝ) : ℝ :=
  sin A * sin C / (sin (A + B) * sin (B + C))

/-- Ptolemy's identity in sine coordinates. -/
theorem sine_ptolemy
    {A B C D : ℝ} (hsum : A + B + C + D = Real.pi) :
    sin (A + B) * sin (B + C) =
      sin A * sin C + sin B * sin D := by
  have hD : sin D = sin (A + B + C) := by
    calc
      sin D = sin (Real.pi - (A + B + C)) := by
        congr 1
        linarith
      _ = sin (A + B + C) := Real.sin_pi_sub _
  rw [hD]
  rw [show A + B + C = (A + B) + C by ring]
  rw [Real.sin_add A B, Real.sin_add B C,
    Real.sin_add (A + B) C, Real.sin_add A B, Real.cos_add A B]
  linear_combination
    (Real.sin A * Real.sin C) * (Real.sin_sq_add_cos_sq B)

private theorem sin_eq_of_add_eq_pi
    {x y : ℝ} (h : x + y = Real.pi) :
    sin x = sin y := by
  calc
    sin x = sin (Real.pi - y) := by rw [show x = Real.pi - y by linarith]
    _ = sin y := Real.sin_pi_sub y

/-- A cyclic sine cross-ratio lies strictly between zero and one. -/
theorem sineCross_mem_Ioo
    {A B C D : ℝ}
    (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hsum : A + B + C + D = Real.pi) :
    0 < sineCross A B C D ∧ sineCross A B C D < 1 := by
  have hAB0 : 0 < A + B := by linarith
  have hBC0 : 0 < B + C := by linarith
  have hApi : A < Real.pi := by linarith
  have hBpi : B < Real.pi := by linarith
  have hCpi : C < Real.pi := by linarith
  have hDpi : D < Real.pi := by linarith
  have hABpi : A + B < Real.pi := by linarith
  have hBCpi : B + C < Real.pi := by linarith
  have hsA := Real.sin_pos_of_pos_of_lt_pi hA hApi
  have hsB := Real.sin_pos_of_pos_of_lt_pi hB hBpi
  have hsC := Real.sin_pos_of_pos_of_lt_pi hC hCpi
  have hsD := Real.sin_pos_of_pos_of_lt_pi hD hDpi
  have hsAB := Real.sin_pos_of_pos_of_lt_pi hAB0 hABpi
  have hsBC := Real.sin_pos_of_pos_of_lt_pi hBC0 hBCpi
  have hden : 0 < sin (A + B) * sin (B + C) := mul_pos hsAB hsBC
  unfold sineCross
  constructor
  · exact div_pos (mul_pos hsA hsC) hden
  · rw [div_lt_one hden, sine_ptolemy hsum]
    nlinarith [mul_pos hsB hsD]

/-- Rotating the four cyclic gaps by one replaces the cross-ratio by its
Euler complement. -/
theorem sineCross_rot_one
    {A B C D : ℝ}
    (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hsum : A + B + C + D = Real.pi) :
    sineCross B C D A = 1 - sineCross A B C D := by
  have hsAB : 0 < sin (A + B) := by
    exact Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsBC : 0 < sin (B + C) := by
    exact Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hCD : sin (C + D) = sin (A + B) := by
    exact (sin_eq_of_add_eq_pi (by linarith)).symm
  unfold sineCross
  rw [hCD]
  field_simp [ne_of_gt hsAB, ne_of_gt hsBC]
  nlinarith [sine_ptolemy hsum]

/-- Rotating the gaps by two leaves the positive cross-ratio unchanged. -/
theorem sineCross_rot_two
    {A B C D : ℝ} (hsum : A + B + C + D = Real.pi) :
    sineCross C D A B = sineCross A B C D := by
  have hCD : sin (C + D) = sin (A + B) :=
    (sin_eq_of_add_eq_pi (by linarith)).symm
  have hDA : sin (D + A) = sin (B + C) :=
    (sin_eq_of_add_eq_pi (by linarith)).symm
  unfold sineCross
  rw [hCD, hDA]
  ring

/-- The reflection preserving the two numerator gaps leaves the cross-ratio
unchanged. -/
theorem sineCross_reflect
    {A B C D : ℝ} (hsum : A + B + C + D = Real.pi) :
    sineCross A D C B = sineCross A B C D := by
  have hAD : sin (A + D) = sin (B + C) :=
    sin_eq_of_add_eq_pi (by linarith)
  have hDC : sin (D + C) = sin (A + B) :=
    sin_eq_of_add_eq_pi (by linarith)
  unfold sineCross
  rw [hAD, hDC]
  ring

/-- The alternating five-face relation for five positive cyclic gaps. -/
theorem rogers_sineCross_five
    {A B C D E : ℝ}
    (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D) (hE : 0 < E)
    (hsum : A + B + C + D + E = Real.pi) :
    rogers (sineCross B C D (E + A))
      - rogers (sineCross (A + B) C D E)
      + rogers (sineCross A (B + C) D E)
      - rogers (sineCross A B (C + D) E)
      + rogers (sineCross A B C (D + E)) = 0 := by
  let z0 := sineCross B C D (E + A)
  let z1 := sineCross (A + B) C D E
  let z2 := sineCross A (B + C) D E
  let z3 := sineCross A B (C + D) E
  let z4 := sineCross A B C (D + E)
  have hz0 := sineCross_mem_Ioo
    (A := B) (B := C) (C := D) (D := E + A)
    hB hC hD (by linarith) (by linarith)
  have hz1 := sineCross_mem_Ioo
    (A := A + B) (B := C) (C := D) (D := E)
    (by linarith) hC hD hE (by linarith)
  have hz2 := sineCross_mem_Ioo
    (A := A) (B := B + C) (C := D) (D := E)
    hA (by linarith) hD hE (by linarith)
  have hz3 := sineCross_mem_Ioo
    (A := A) (B := B) (C := C + D) (D := E)
    hA hB (by linarith) hE (by linarith)
  have hz4 := sineCross_mem_Ioo
    (A := A) (B := B) (C := C) (D := D + E)
    hA hB hC (by linarith) (by linarith)

  have hsA : sin A ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi hA (by linarith)
  have hsB : sin B ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi hB (by linarith)
  have hsC : sin C ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi hC (by linarith)
  have hsD : sin D ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi hD (by linarith)
  have hsE : sin E ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi hE (by linarith)
  have hsAB : sin (A + B) ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsBC : sin (B + C) ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsCD : sin (C + D) ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsDE : sin (D + E) ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsABC : sin (A + B + C) ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsBCD : sin (B + C + D) ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hsCDE : sin (C + D + E) ≠ 0 := ne_of_gt <|
    Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)

  have hp : z1 * z3 = z2 := by
    dsimp [z1, z2, z3, sineCross]
    ring_nf
    field_simp [hsA, hsB, hsC, hsD, hsE, hsAB, hsBC, hsCD, hsDE,
      hsABC, hsBCD, hsCDE]

  have hcomp3 :
      sineCross B (C + D) E A = 1 - z3 := by
    dsimp [z3]
    exact sineCross_rot_one hA hB (by linarith) hE (by linarith)
  have hcomp2 :
      sineCross (B + C) D E A = 1 - z2 := by
    dsimp [z2]
    exact sineCross_rot_one hA (by linarith) hD hE (by linarith)
  have hcomp1 :
      sineCross C D E (A + B) = 1 - z1 := by
    dsimp [z1]
    exact sineCross_rot_one (by linarith) hC hD hE (by linarith)

  have hAB_CDE : sin (A + B) = sin (C + D + E) :=
    sin_eq_of_add_eq_pi (by linarith)
  have hABC_DE : sin (A + B + C) = sin (D + E) :=
    sin_eq_of_add_eq_pi (by linarith)
  have hBC_DEA : sin (B + C) = sin (D + E + A) :=
    sin_eq_of_add_eq_pi (by linarith)

  have hq : z1 * (1 - z3) / (1 - z1 * z3) = z0 := by
    rw [hp, ← hcomp3, ← hcomp2]
    dsimp [z0, z1, sineCross]
    rw [hAB_CDE, hABC_DE]
    rw [show B + (C + D) = B + C + D by ring]
    field_simp [hsA, hsB, hsC, hsD, hsE, hsAB, hsBC, hsCD, hsDE,
      hsABC, hsBCD, hsCDE]

  have hr : z3 * (1 - z1) / (1 - z1 * z3) = z4 := by
    rw [hp, ← hcomp1, ← hcomp2]
    dsimp [z3, z4, sineCross]
    rw [hBC_DEA]
    rw [show B + (C + D) = B + C + D by ring]
    field_simp [hsA, hsB, hsC, hsD, hsE, hsAB, hsBC, hsCD, hsDE,
      hsABC, hsBCD, hsCDE]

  have hf :
      rogers z1 + rogers z3 =
        rogers (z1 * z3)
          + rogers (z1 * (1 - z3) / (1 - z1 * z3))
          + rogers (z3 * (1 - z1) / (1 - z1 * z3)) :=
    rogers_five_term hz1.1 hz1.2 hz3.1 hz3.2
  rw [hq, hr, hp] at hf
  dsimp [z0, z1, z2, z3, z4] at hf ⊢
  linarith

end Real
