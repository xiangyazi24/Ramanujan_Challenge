import Ramanujan31.APolynomial

/-!
# The real Rogers shape chamber

On the official real component the auxiliary chart root is negative.  This
file proves that, once this sign is known, all four real Rogers arguments lie
strictly between zero and one.  The inequalities use only the cleared
holonomy equation and `0 < X < 1`, `X ≤ L`.
-/

/-- First real Rogers argument obtained from the negative shape. -/
noncomputable def rogersArgumentA (r : ℝ) : ℝ := 1 - 1 / r ^ 2

/-- Second real Rogers argument. -/
noncomputable def rogersArgumentB (r : ℝ) : ℝ := chartT r / r

/-- Third real Rogers argument. -/
noncomputable def rogersArgumentC (X r : ℝ) : ℝ := chartT r * X / r

/-- Fourth real Rogers argument. -/
noncomputable def rogersArgumentD (X r : ℝ) : ℝ :=
  (chartT r - r * X) / chartT r

theorem rogersArgumentB_eq_of_holonomy
    {X L r : ℝ} (hX : 0 < X) (hL : 0 < L) (hr : r < -1)
    (hH : holonomyEquation X L r = 0) :
    rogersArgumentB r = X * (L + X) / (L + X ^ 3) := by
  have hr0 : r ≠ 0 := by linarith
  have hsum : L + X ^ 3 ≠ 0 := by positivity
  unfold rogersArgumentB holonomyEquation at *
  field_simp [hr0, hsum]
  linear_combination hH

theorem rogersArgumentA_mem_Ioo
    {r : ℝ} (hr : r < -1) :
    0 < rogersArgumentA r ∧ rogersArgumentA r < 1 := by
  have hr0 : r ≠ 0 := by linarith
  have hr2 : 1 < r ^ 2 := by nlinarith
  unfold rogersArgumentA
  constructor
  · have hinvlt : 1 / r ^ 2 < 1 := by
      exact (div_lt_one (sq_pos_of_ne_zero hr0)).mpr hr2
    linarith
  · have hinv : 0 < 1 / r ^ 2 := by positivity
    linarith

theorem rogersArgumentB_bounds
    {X L r : ℝ} (hX : 0 < X) (hX1 : X < 1) (hXL : X ≤ L)
    (hr : r < -1) (hH : holonomyEquation X L r = 0) :
    X < rogersArgumentB r ∧ rogersArgumentB r < 1 := by
  have hL : 0 < L := lt_of_lt_of_le hX hXL
  have hsum : 0 < L + X ^ 3 := by positivity
  have hB := rogersArgumentB_eq_of_holonomy hX hL hr hH
  rw [hB]
  constructor
  · rw [lt_div_iff₀ hsum]
    have hcube : X ^ 3 < X := by
      have hX2 : X ^ 2 < 1 := by nlinarith
      nlinarith
    nlinarith
  · rw [div_lt_one hsum]
    have hX2 : X ^ 2 < X := by nlinarith
    nlinarith

theorem rogersArgumentC_eq
    (X r : ℝ) :
    rogersArgumentC X r = rogersArgumentB r * X := by
  unfold rogersArgumentC rogersArgumentB
  ring

theorem rogersArgumentD_eq
    {X r : ℝ} (hr : r < -1) :
    rogersArgumentD X r = 1 - X / rogersArgumentB r := by
  have hr0 : r ≠ 0 := by linarith
  have ht : chartT r ≠ 0 := by
    unfold chartT
    nlinarith
  unfold rogersArgumentD rogersArgumentB
  field_simp [hr0, ht]

/-- All four transformed real shapes lie in the open unit interval. -/
theorem rogersArguments_mem_Ioo
    {X L r : ℝ} (hX : 0 < X) (hX1 : X < 1) (hXL : X ≤ L)
    (hr : r < -1) (hH : holonomyEquation X L r = 0) :
    (0 < rogersArgumentA r ∧ rogersArgumentA r < 1)
      ∧ (0 < rogersArgumentB r ∧ rogersArgumentB r < 1)
      ∧ (0 < rogersArgumentC X r ∧ rogersArgumentC X r < 1)
      ∧ (0 < rogersArgumentD X r ∧ rogersArgumentD X r < 1) := by
  have hA := rogersArgumentA_mem_Ioo hr
  have hBbounds := rogersArgumentB_bounds hX hX1 hXL hr hH
  have hB0 : 0 < rogersArgumentB r := lt_trans hX hBbounds.1
  have hB : 0 < rogersArgumentB r ∧ rogersArgumentB r < 1 :=
    ⟨hB0, hBbounds.2⟩
  have hCeq := rogersArgumentC_eq X r
  have hC : 0 < rogersArgumentC X r ∧ rogersArgumentC X r < 1 := by
    rw [hCeq]
    constructor
    · positivity
    · nlinarith [mul_lt_mul_of_pos_right hB.2 hX]
  have hDeq := rogersArgumentD_eq (X := X) hr
  have hD : 0 < rogersArgumentD X r ∧ rogersArgumentD X r < 1 := by
    rw [hDeq]
    have hratio0 : 0 < X / rogersArgumentB r := div_pos hX hB0
    have hratio1 : X / rogersArgumentB r < 1 := by
      exact (div_lt_one hB0).mpr hBbounds.1
    constructor <;> linarith
  exact ⟨hA, hB, hC, hD⟩
