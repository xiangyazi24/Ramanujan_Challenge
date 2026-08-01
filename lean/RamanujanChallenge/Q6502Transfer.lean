import RamanujanChallenge.Problem22Concentration

open Filter Topology Real
open scoped BigOperators

noncomputable section

namespace RamanujanChallenge.P22.Q6502

/-- Uniform harmonic-minus-log remainder using `log (m+1)`, valid also at `m=0`. -/
theorem abs_harmonic_sub_log_succ_sub_gamma_le (m : ℕ) :
    |((harmonic m : ℚ) : ℝ) - Real.log ((m : ℝ) + 1) -
        Real.eulerMascheroniConstant| ≤
      1 / ((m : ℝ) + 1) := by
  have hlo := Real.eulerMascheroniSeq_lt_eulerMascheroniConstant m
  have hhi := Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' (m + 1)
  simp only [Real.eulerMascheroniSeq, Real.eulerMascheroniSeq',
    Nat.succ_ne_zero, if_false, harmonic_succ, Rat.cast_add,
    Rat.cast_inv, Rat.cast_natCast, Rat.cast_one, Nat.cast_add,
    Nat.cast_one] at hlo hhi
  have hneg :
      ((harmonic m : ℚ) : ℝ) - Real.log ((m : ℝ) + 1) -
          Real.eulerMascheroniConstant ≤ 0 := by
    linarith
  rw [abs_of_nonpos hneg]
  have hbound :
      -(((harmonic m : ℚ) : ℝ) - Real.log ((m : ℝ) + 1) -
          Real.eulerMascheroniConstant) <
        (((m : ℝ) + 1)⁻¹) := by
    linarith
  simpa only [one_div] using hbound.le

/-- Lipschitz bound for the logarithm on a positive interval. -/
theorem abs_log_sub_log_le_abs_sub_div_min
    {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    |Real.log a - Real.log b| ≤ |a - b| / min a b := by
  rcases le_total a b with hab | hba
  · have hlog : Real.log a ≤ Real.log b := Real.log_le_log ha hab
    rw [abs_of_nonpos (sub_nonpos.mpr hlog), abs_of_nonpos (sub_nonpos.mpr hab),
      min_eq_left hab]
    have hratio : 0 < b / a := div_pos hb ha
    have hle := Real.log_le_sub_one_of_pos hratio
    rw [Real.log_div hb.ne' ha.ne'] at hle
    have ha0 : a ≠ 0 := ha.ne'
    field_simp [ha0] at hle ⊢
    nlinarith
  · have hlog : Real.log b ≤ Real.log a := Real.log_le_log hb hba
    rw [abs_of_nonneg (sub_nonneg.mpr hlog), abs_of_nonneg (sub_nonneg.mpr hba),
      min_eq_right hba]
    have hratio : 0 < a / b := div_pos ha hb
    have hle := Real.log_le_sub_one_of_pos hratio
    rw [Real.log_div ha.ne' hb.ne'] at hle
    have hb0 : b ≠ 0 := hb.ne'
    field_simp [hb0] at hle ⊢
    nlinarith

/-- The logarithmic saddle observable with shifted positive arguments. -/
def rivoalLogSaddle22 (n k : ℕ) : ℝ :=
  3 * Real.log ((k : ℝ) + 1) -
    2 * Real.log (((n - k : ℕ) : ℝ) + 1)

/-- Exact harmonic-minus-log decomposition, with no `log 0` branch. -/
theorem rivoalHarmonic_error_le_remainders_log
    (n k : ℕ) (hk : k ≤ n) :
    |rivoalRealHarmonicValue22 n k - Real.eulerMascheroniConstant| ≤
      3 * (1 / ((k : ℝ) + 1)) +
        2 * (1 / (((n - k : ℕ) : ℝ) + 1)) +
          |rivoalLogSaddle22 n k| := by
  let rk : ℝ :=
    ((harmonic k : ℚ) : ℝ) - Real.log ((k : ℝ) + 1) -
      Real.eulerMascheroniConstant
  let rj : ℝ :=
    ((harmonic (n - k) : ℚ) : ℝ) -
      Real.log (((n - k : ℕ) : ℝ) + 1) -
        Real.eulerMascheroniConstant
  have hkrem : |rk| ≤ 1 / ((k : ℝ) + 1) := by
    simpa only [rk] using abs_harmonic_sub_log_succ_sub_gamma_le k
  have hjrem : |rj| ≤ 1 / (((n - k : ℕ) : ℝ) + 1) := by
    simpa only [rj] using abs_harmonic_sub_log_succ_sub_gamma_le (n - k)
  have hid :
      rivoalRealHarmonicValue22 n k - Real.eulerMascheroniConstant =
        3 * rk - 2 * rj + rivoalLogSaddle22 n k := by
    simp only [rivoalRealHarmonicValue22, rivoalHarmonicKernel22,
      rivoalLogSaddle22, rk, rj]
    push_cast
    ring
  rw [hid]
  calc
    |3 * rk - 2 * rj + rivoalLogSaddle22 n k| ≤
        |3 * rk - 2 * rj| + |rivoalLogSaddle22 n k| :=
      abs_add_le _ _
    _ ≤ (|3 * rk| + |2 * rj|) + |rivoalLogSaddle22 n k| := by
      exact add_le_add_right (abs_sub_le _ _) _
    _ = (3 * |rk| + 2 * |rj|) + |rivoalLogSaddle22 n k| := by
      rw [abs_mul, abs_mul]
      norm_num
    _ ≤ 3 * (1 / ((k : ℝ) + 1)) +
          2 * (1 / (((n - k : ℕ) : ℝ) + 1)) +
            |rivoalLogSaddle22 n k| := by
      gcongr

/-- Crude global envelope used only on the small bad set. -/
theorem rivoalHarmonic_error_le_envelope
    (n k : ℕ) (hk : k ≤ n) :
    |rivoalRealHarmonicValue22 n k - Real.eulerMascheroniConstant| ≤
      5 + 5 * Real.log ((n : ℝ) + 1) := by
  have hbase := rivoalHarmonic_error_le_remainders_log n k hk
  have hkrec : 1 / ((k : ℝ) + 1) ≤ 1 := by
    apply (div_le_iff₀ (by positivity : 0 < (k : ℝ) + 1)).2
    linarith
  have hjrec : 1 / (((n - k : ℕ) : ℝ) + 1) ≤ 1 := by
    apply (div_le_iff₀ (by positivity :
      0 < ((n - k : ℕ) : ℝ) + 1)).2
    linarith
  have hkarg : (k : ℝ) + 1 ≤ (n : ℝ) + 1 := by
    exact_mod_cast Nat.add_le_add_right hk 1
  have hjarg : ((n - k : ℕ) : ℝ) + 1 ≤ (n : ℝ) + 1 := by
    exact_mod_cast Nat.add_le_add_right (Nat.sub_le n k) 1
  have hlogk0 : 0 ≤ Real.log ((k : ℝ) + 1) :=
    Real.log_nonneg (by positivity)
  have hlogj0 : 0 ≤ Real.log (((n - k : ℕ) : ℝ) + 1) :=
    Real.log_nonneg (by positivity)
  have hlogk := Real.log_le_log (by positivity : 0 < (k : ℝ) + 1) hkarg
  have hlogj := Real.log_le_log
    (by positivity : 0 < ((n - k : ℕ) : ℝ) + 1) hjarg
  have hlogs :
      |rivoalLogSaddle22 n k| ≤ 5 * Real.log ((n : ℝ) + 1) := by
    unfold rivoalLogSaddle22
    calc
      |3 * Real.log ((k : ℝ) + 1) -
          2 * Real.log (((n - k : ℕ) : ℝ) + 1)| ≤
        |3 * Real.log ((k : ℝ) + 1)| +
          |2 * Real.log (((n - k : ℕ) : ℝ) + 1)| := abs_sub_le _ _
      _ = 3 * Real.log ((k : ℝ) + 1) +
          2 * Real.log (((n - k : ℕ) : ℝ) + 1) := by
        rw [abs_mul, abs_mul, abs_of_nonneg hlogk0, abs_of_nonneg hlogj0]
        norm_num
      _ ≤ 5 * Real.log ((n : ℝ) + 1) := by
        nlinarith
  linarith

#print axioms abs_harmonic_sub_log_succ_sub_gamma_le
#print axioms abs_log_sub_log_le_abs_sub_div_min
#print axioms rivoalHarmonic_error_le_remainders_log
#print axioms rivoalHarmonic_error_le_envelope

end RamanujanChallenge.P22.Q6502
