import RamanujanChallenge.Problem22Moment
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.NumberTheory.Harmonic.Bounds

/-!
# Problem 2.2: from cubic-saddle concentration to harmonic concentration
-/

noncomputable section

open Filter Topology Real
open scoped BigOperators

namespace RamanujanChallenge.P22

/-- The error in the classical approximation `H_m - log m -> gamma`. -/
def harmonicRemainder22 (m : ℕ) : ℝ :=
  (harmonic m : ℝ) - Real.log (m : ℝ) - eulerMascheroniConstant

theorem harmonicRemainder22_nonneg
    {m : ℕ} (hm : 1 ≤ m) :
    0 ≤ harmonicRemainder22 m := by
  have h := Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' m
  simp only [Real.eulerMascheroniSeq', if_neg (Nat.one_le_iff_ne_zero.mp hm)] at h
  unfold harmonicRemainder22
  linarith

/-- Explicit one-sided error bound for harmonic numbers. -/
theorem harmonicRemainder22_le_inv
    {m : ℕ} (hm : 1 ≤ m) :
    harmonicRemainder22 m ≤ 1 / (m : ℝ) := by
  have hmreal : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hlower := Real.eulerMascheroniSeq_lt_eulerMascheroniConstant m
  have hupper := Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' m
  simp only [Real.eulerMascheroniSeq, Real.eulerMascheroniSeq',
    if_neg (Nat.one_le_iff_ne_zero.mp hm)] at hlower hupper
  have hgap :
      harmonicRemainder22 m ≤
        Real.log ((m : ℝ) + 1) - Real.log (m : ℝ) := by
    unfold harmonicRemainder22
    linarith
  have hlog :
      Real.log ((m : ℝ) + 1) - Real.log (m : ℝ) ≤ 1 / (m : ℝ) := by
    rw [← Real.log_div (by positivity) hmreal.ne']
    have hpos : 0 < ((m : ℝ) + 1) / (m : ℝ) := div_pos (by positivity) hmreal
    calc
      Real.log (((m : ℝ) + 1) / (m : ℝ)) ≤
          ((m : ℝ) + 1) / (m : ℝ) - 1 :=
        Real.log_le_sub_one_of_pos hpos
      _ = 1 / (m : ℝ) := by field_simp; ring
  exact hgap.trans hlog

theorem abs_harmonicRemainder22_le_inv
    {m : ℕ} (hm : 1 ≤ m) :
    |harmonicRemainder22 m| ≤ 1 / (m : ℝ) := by
  rw [abs_of_nonneg (harmonicRemainder22_nonneg hm)]
  exact harmonicRemainder22_le_inv hm

/-- Algebraic decomposition of the sampled harmonic value into two classical
harmonic remainders and the logarithmic saddle coordinate. -/
theorem rivoalHarmonicValue22_decompose
    {n k : ℕ} (_hk : k ≤ n) :
    rivoalRealHarmonicValue22 n k - eulerMascheroniConstant =
      3 * harmonicRemainder22 k - 2 * harmonicRemainder22 (n - k) +
        (3 * Real.log (k : ℝ) - 2 * Real.log ((n - k : ℕ) : ℝ)) := by
  simp only [rivoalRealHarmonicValue22, rivoalHarmonicKernel22,
    Rat.cast_sub, Rat.cast_mul, Rat.cast_ofNat, harmonicRemainder22]
  ring

/-- Logarithmic saddle coordinate as one logarithm. -/
theorem rivoalLogSaddle22_eq
    {n k : ℕ} (hk0 : 1 ≤ k) (hkn : k < n) :
    3 * Real.log (k : ℝ) - 2 * Real.log ((n - k : ℕ) : ℝ) =
      Real.log (((k : ℝ) ^ 3) / (((n - k : ℕ) : ℝ) ^ 2)) := by
  have hkpos : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk0
  have hdiff : 0 < n - k := Nat.sub_pos_of_lt hkn
  have hdiffpos : (0 : ℝ) < ((n - k : ℕ) : ℝ) := by exact_mod_cast hdiff
  rw [Real.log_div (pow_ne_zero _ hkpos.ne') (pow_ne_zero _ hdiffpos.ne'),
    Real.log_pow, Real.log_pow]
  norm_num

/-- A relative error at most `1/2` gives a robust logarithm bound. -/
theorem abs_log_le_two_mul_abs_sub_one22
    {x : ℝ} (hx : 0 < x) (hnear : |x - 1| ≤ 1 / 2) :
    |Real.log x| ≤ 2 * |x - 1| := by
  have hxlower : 1 / 2 ≤ x := by
    rw [abs_le] at hnear
    linarith
  rcases le_total x 1 with hx1 | h1x
  · have hlog_nonpos : Real.log x ≤ 0 := Real.log_nonpos hx.le hx1
    rw [abs_of_nonpos hlog_nonpos, abs_of_nonpos (sub_nonpos.mpr hx1)]
    have h := Real.log_le_sub_one_of_pos (inv_pos.mpr hx)
    rw [Real.log_inv] at h
    calc
      -Real.log x ≤ x⁻¹ - 1 := h
      _ = (1 - x) / x := by field_simp [hx.ne']
      _ ≤ 2 * (1 - x) := by
        apply (div_le_iff₀ hx).2
        have hprod : 0 ≤ (1 - x) * (2 * x - 1) :=
          mul_nonneg (sub_nonneg.mpr hx1) (by linarith)
        nlinarith
      _ = 2 * -(x - 1) := by ring
  · have hlog_nonneg : 0 ≤ Real.log x := Real.log_nonneg h1x
    rw [abs_of_nonneg hlog_nonneg, abs_of_nonneg (sub_nonneg.mpr h1x)]
    calc
      Real.log x ≤ x - 1 := Real.log_le_sub_one_of_pos hx
      _ ≤ 2 * (x - 1) := by nlinarith

/-- Geometry forced by a small normalized cubic-saddle error. -/
theorem rivoalGoodSaddleGeometry22
    {δ : ℝ} (_hδ0 : 0 ≤ δ) (hδ : δ ≤ 1 / 16)
    {n k : ℕ} (hn : 16 ≤ n) (hkn : k ≤ n)
    (hgood : |rivoalSaddleError22 n k| ≤ δ * (n : ℝ) ^ 2) :
    1 ≤ k ∧ k < n ∧
      (3 / 16 : ℝ) * (n : ℝ) ^ 2 ≤ (k : ℝ) ^ 3 ∧
      (n : ℝ) / 2 ≤ (n : ℝ) - (k : ℝ) := by
  have hn16 : (16 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hn0 : (0 : ℝ) < (n : ℝ) := by linarith
  have hk0 : (0 : ℝ) ≤ (k : ℝ) := by positivity
  have hcast : (k : ℝ) ≤ (n : ℝ) := by exact_mod_cast hkn
  have hdiff0 : 0 ≤ (n : ℝ) - (k : ℝ) := sub_nonneg.mpr hcast
  have hdiff_le : (n : ℝ) - (k : ℝ) ≤ (n : ℝ) := by linarith
  rw [abs_le] at hgood
  have hδn : δ * (n : ℝ) ^ 2 ≤ (1 / 16 : ℝ) * (n : ℝ) ^ 2 := by
    exact mul_le_mul_of_nonneg_right hδ (sq_nonneg _)
  have hn3 : 16 * (n : ℝ) ^ 2 ≤ (n : ℝ) ^ 3 := by
    nlinarith [mul_nonneg (sq_nonneg (n : ℝ)) (sub_nonneg.mpr hn16)]
  have hhalf : (k : ℝ) ≤ (n : ℝ) / 2 := by
    by_contra hnot
    have hkhalf : (n : ℝ) / 2 < (k : ℝ) := lt_of_not_ge hnot
    have hcub : ((n : ℝ) / 2) ^ 3 ≤ (k : ℝ) ^ 3 := by gcongr
    have hdiffsq : ((n : ℝ) - (k : ℝ)) ^ 2 ≤ (n : ℝ) ^ 2 := by gcongr
    simp only [rivoalSaddleError22] at hgood
    nlinarith [sq_pos_of_pos hn0]
  have hdiffhalf : (n : ℝ) / 2 ≤ (n : ℝ) - (k : ℝ) := by linarith
  have hdiffsq_lower : (n : ℝ) ^ 2 / 4 ≤ ((n : ℝ) - (k : ℝ)) ^ 2 := by
    nlinarith [sq_nonneg ((n : ℝ) - (k : ℝ) - (n : ℝ) / 2)]
  have hkcube :
      (3 / 16 : ℝ) * (n : ℝ) ^ 2 ≤ (k : ℝ) ^ 3 := by
    simp only [rivoalSaddleError22] at hgood
    nlinarith
  have hk_ne : k ≠ 0 := by
    intro hkz
    subst k
    norm_num at hkcube
    nlinarith [sq_pos_of_pos hn0]
  have hk1 : 1 ≤ k := Nat.one_le_iff_ne_zero.mpr hk_ne
  have h2cast : (2 * k : ℕ) ≤ n := by
    exact_mod_cast (show 2 * (k : ℝ) ≤ (n : ℝ) by linarith)
  have hklt : k < n := by omega
  exact ⟨hk1, hklt, hkcube, hdiffhalf⟩

/-- On the cubic-saddle good set the sampled harmonic value is uniformly close
to Euler's constant. -/
theorem rivoalGoodSaddleHarmonicError22
    {δ : ℝ} (hδ0 : 0 ≤ δ) (hδ : δ ≤ 1 / 16)
    {M n k : ℕ} (hM : 1 ≤ M) (hn : 16 ≤ n) (hkn : k ≤ n)
    (hMn : 16 * (M : ℝ) ^ 3 ≤ 3 * (n : ℝ) ^ 2)
    (hgood : |rivoalSaddleError22 n k| ≤ δ * (n : ℝ) ^ 2) :
    |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant| ≤
      5 / (M : ℝ) + 8 * δ := by
  obtain ⟨hk1, hklt, hkcube, hdiffhalf⟩ :=
    rivoalGoodSaddleGeometry22 hδ0 hδ hn hkn hgood
  have hMreal : (0 : ℝ) < (M : ℝ) := by exact_mod_cast hM
  have hkreal : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk1
  have hdiffnat : 1 ≤ n - k := Nat.sub_pos_iff_lt.mpr hklt
  have hdiffreal : (0 : ℝ) < ((n - k : ℕ) : ℝ) := by exact_mod_cast hdiffnat
  have hMk : M ≤ k := by
    by_contra hnot
    have hkM : k < M := Nat.lt_of_not_ge hnot
    have hcubelt : (k : ℝ) ^ 3 < (M : ℝ) ^ 3 := by
      gcongr
    nlinarith
  have hMdiff : M ≤ n - k := by
    have hcastsub : (((n - k : ℕ) : ℝ)) = (n : ℝ) - (k : ℝ) := by
      exact Nat.cast_sub hkn
    exact_mod_cast (show (M : ℝ) ≤ ((n - k : ℕ) : ℝ) by
      rw [hcastsub]
      have hMhalf : (M : ℝ) ≤ (n : ℝ) / 2 := by
        by_contra hnot
        have hlt : (n : ℝ) / 2 < (M : ℝ) := lt_of_not_ge hnot
        have hcub : ((n : ℝ) / 2) ^ 3 < (M : ℝ) ^ 3 := by gcongr
        have hn16 : (16 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
        have hn3 : 16 * (n : ℝ) ^ 2 ≤ (n : ℝ) ^ 3 := by
          nlinarith [mul_nonneg (sq_nonneg (n : ℝ)) (sub_nonneg.mpr hn16)]
        nlinarith
      linarith)
  have hrk : |harmonicRemainder22 k| ≤ 1 / (M : ℝ) := by
    exact (abs_harmonicRemainder22_le_inv hk1).trans
      (one_div_le_one_div_of_le hMreal (by exact_mod_cast hMk))
  have hrd : |harmonicRemainder22 (n - k)| ≤ 1 / (M : ℝ) := by
    exact (abs_harmonicRemainder22_le_inv hdiffnat).trans
      (one_div_le_one_div_of_le hMreal (by exact_mod_cast hMdiff))
  let x : ℝ := (k : ℝ) ^ 3 / (((n - k : ℕ) : ℝ) ^ 2)
  have hx : 0 < x := by
    dsimp [x]
    positivity
  have hcastsub : (((n - k : ℕ) : ℝ)) = (n : ℝ) - (k : ℝ) := by
    exact Nat.cast_sub hkn
  have hdenlower : (n : ℝ) ^ 2 / 4 ≤ (((n - k : ℕ) : ℝ) ^ 2) := by
    rw [hcastsub]
    nlinarith [sq_nonneg ((n : ℝ) - (k : ℝ) - (n : ℝ) / 2)]
  have hxsub : |x - 1| =
      |rivoalSaddleError22 n k| / (((n - k : ℕ) : ℝ) ^ 2) := by
    have hdenne : (((n - k : ℕ) : ℝ) ^ 2) ≠ 0 :=
      pow_ne_zero _ hdiffreal.ne'
    have hxalg : x - 1 =
        rivoalSaddleError22 n k / (((n - k : ℕ) : ℝ) ^ 2) := by
      dsimp [x]
      rw [div_sub_one hdenne]
      simp only [rivoalSaddleError22, hcastsub]
    rw [hxalg, abs_div, abs_of_pos (sq_pos_of_pos hdiffreal)]
  have hxnear : |x - 1| ≤ 4 * δ := by
    rw [hxsub]
    apply (div_le_iff₀ (sq_pos_of_pos hdiffreal)).2
    have hscaled := mul_le_mul_of_nonneg_left hdenlower (show 0 ≤ 4 * δ by positivity)
    nlinarith
  have hxhalf : |x - 1| ≤ 1 / 2 :=
    hxnear.trans (by nlinarith)
  have hlog : |Real.log x| ≤ 8 * δ := by
    calc
      |Real.log x| ≤ 2 * |x - 1| :=
        abs_log_le_two_mul_abs_sub_one22 hx hxhalf
      _ ≤ 8 * δ := by nlinarith
  rw [rivoalHarmonicValue22_decompose hkn,
    rivoalLogSaddle22_eq hk1 hklt]
  change |3 * harmonicRemainder22 k - 2 * harmonicRemainder22 (n - k) +
      Real.log x| ≤ _
  calc
    _ ≤ 3 * |harmonicRemainder22 k| +
        2 * |harmonicRemainder22 (n - k)| + |Real.log x| := by
          calc
            _ ≤ |3 * harmonicRemainder22 k -
                  2 * harmonicRemainder22 (n - k)| + |Real.log x| :=
              abs_add_le _ _
            _ ≤ (|3 * harmonicRemainder22 k| +
                  |-2 * harmonicRemainder22 (n - k)|) + |Real.log x| := by
              gcongr
              simpa [sub_eq_add_neg] using
                (abs_add_le (3 * harmonicRemainder22 k)
                  (-2 * harmonicRemainder22 (n - k)))
            _ = _ := by simp [abs_mul]
    _ ≤ 3 * (1 / (M : ℝ)) + 2 * (1 / (M : ℝ)) + 8 * δ := by
          exact add_le_add (add_le_add (mul_le_mul_of_nonneg_left hrk (by norm_num))
            (mul_le_mul_of_nonneg_left hrd (by norm_num))) hlog
    _ = 5 / (M : ℝ) + 8 * δ := by ring

/-- A coarse uniform bound used only on the complement of the saddle set. -/
def rivoalHarmonicEnvelope22 (n : ℕ) : ℝ :=
  6 + 5 * Real.log (n : ℝ)

theorem harmonic_nonneg_real22 (m : ℕ) :
    0 ≤ (harmonic m : ℝ) := by
  have hone : (1 : ℝ) ≤ ((m + 1 : ℕ) : ℝ) := by
    exact_mod_cast (show 1 ≤ m + 1 by omega)
  exact (Real.log_nonneg hone).trans (log_add_one_le_harmonic m)

theorem harmonic_le_one_add_log_of_le22
    {m n : ℕ} (hn : 1 ≤ n) (hmn : m ≤ n) :
    (harmonic m : ℝ) ≤ 1 + Real.log (n : ℝ) := by
  by_cases hm0 : m = 0
  · subst m
    simp only [harmonic_zero, Rat.cast_zero]
    have hnreal : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
    positivity
  · have hmreal : (0 : ℝ) < (m : ℝ) := by
      exact_mod_cast (Nat.pos_of_ne_zero hm0)
    have hmcast : (m : ℝ) ≤ (n : ℝ) := by exact_mod_cast hmn
    exact (harmonic_le_one_add_log m).trans
      (add_le_add (le_refl 1) (Real.log_le_log hmreal hmcast))

theorem rivoalHarmonicError22_le_envelope
    {n k : ℕ} (hn : 1 ≤ n) (hkn : k ≤ n) :
    |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant| ≤
      rivoalHarmonicEnvelope22 n := by
  have hk0 := harmonic_nonneg_real22 k
  have hd0 := harmonic_nonneg_real22 (n - k)
  have hk := harmonic_le_one_add_log_of_le22 hn hkn
  have hd := harmonic_le_one_add_log_of_le22 hn (Nat.sub_le n k)
  have hlog : 0 ≤ Real.log (n : ℝ) := by
    exact Real.log_nonneg (by exact_mod_cast hn)
  have hγ0 : 0 < eulerMascheroniConstant :=
    (by linarith [Real.one_half_lt_eulerMascheroniConstant])
  have hγ1 : eulerMascheroniConstant < 1 :=
    Real.eulerMascheroniConstant_lt_two_thirds.trans (by norm_num)
  simp only [rivoalRealHarmonicValue22, rivoalHarmonicKernel22,
    Rat.cast_sub, Rat.cast_mul, Rat.cast_ofNat,
    rivoalHarmonicEnvelope22]
  rw [abs_le]
  constructor <;> nlinarith

/-- Off the saddle set, the normalized square saddle error is at least one. -/
theorem one_le_rivoalSaddleError22_sq_div
    {δ : ℝ} (hδ : 0 < δ) {n k : ℕ} (hn : 1 ≤ n)
    (hbad : δ * (n : ℝ) ^ 2 < |rivoalSaddleError22 n k|) :
    1 ≤ (rivoalSaddleError22 n k) ^ 2 /
      (δ ^ 2 * (n : ℝ) ^ 4) := by
  have hnreal : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn
  have hscale : 0 < δ * (n : ℝ) ^ 2 := by positivity
  have hsquare :
      (δ * (n : ℝ) ^ 2) ^ 2 ≤
        |rivoalSaddleError22 n k| ^ 2 := by
    nlinarith [sq_nonneg
      (|rivoalSaddleError22 n k| - δ * (n : ℝ) ^ 2)]
  rw [sq_abs] at hsquare
  apply (le_div_iff₀ (show 0 < δ ^ 2 * (n : ℝ) ^ 4 by positivity)).2
  convert hsquare using 1
  ring

/-- Pointwise good-set/bad-set decomposition. -/
theorem rivoalHarmonicError22_pointwise_le
    {δ : ℝ} (hδ : 0 < δ) (hδle : δ ≤ 1 / 16)
    {M n k : ℕ} (hM : 1 ≤ M) (hn : 16 ≤ n) (hkn : k ≤ n)
    (hMn : 16 * (M : ℝ) ^ 3 ≤ 3 * (n : ℝ) ^ 2) :
    |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant| ≤
      (5 / (M : ℝ) + 8 * δ) +
        rivoalHarmonicEnvelope22 n *
          ((rivoalSaddleError22 n k) ^ 2 /
            (δ ^ 2 * (n : ℝ) ^ 4)) := by
  have hn1 : 1 ≤ n := hn.trans' (by norm_num)
  have henv : 0 ≤ rivoalHarmonicEnvelope22 n := by
    unfold rivoalHarmonicEnvelope22
    have : 0 ≤ Real.log (n : ℝ) :=
      Real.log_nonneg (by exact_mod_cast hn1)
    positivity
  have hmain : 0 ≤ 5 / (M : ℝ) + 8 * δ := by positivity
  by_cases hgood :
      |rivoalSaddleError22 n k| ≤ δ * (n : ℝ) ^ 2
  · exact (rivoalGoodSaddleHarmonicError22 hδ.le hδle hM hn hkn hMn hgood).trans
      (le_add_of_nonneg_right
        (mul_nonneg henv (div_nonneg (sq_nonneg _) (by positivity))))
  · have hbad : δ * (n : ℝ) ^ 2 < |rivoalSaddleError22 n k| :=
      lt_of_not_ge hgood
    have hratio := one_le_rivoalSaddleError22_sq_div hδ hn1 hbad
    calc
      |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant| ≤
          rivoalHarmonicEnvelope22 n :=
        rivoalHarmonicError22_le_envelope hn1 hkn
      _ ≤ rivoalHarmonicEnvelope22 n *
          ((rivoalSaddleError22 n k) ^ 2 /
            (δ ^ 2 * (n : ℝ) ^ 4)) :=
        (by simpa using mul_le_mul_of_nonneg_left hratio henv)
      _ ≤ (5 / (M : ℝ) + 8 * δ) +
          rivoalHarmonicEnvelope22 n *
            ((rivoalSaddleError22 n k) ^ 2 /
              (δ ^ 2 * (n : ℝ) ^ 4)) :=
        le_add_of_nonneg_left hmain

/-- Quantitative harmonic concentration bound obtained from the second moment. -/
theorem rivoalWeightedHarmonicError22_le
    {δ : ℝ} (hδ : 0 < δ) (hδle : δ ≤ 1 / 16)
    {M n : ℕ} (hM : 1 ≤ M) (hn : 16 ≤ n)
    (hMn : 16 * (M : ℝ) ^ 3 ≤ 3 * (n : ℝ) ^ 2) :
    rivoalWeightedHarmonicError22 n ≤
      5 / (M : ℝ) + 8 * δ +
        81 * rivoalHarmonicEnvelope22 n /
          (δ ^ 2 * Real.sqrt (n : ℝ)) := by
  have hn1 : 1 ≤ n := hn.trans' (by norm_num)
  have hnreal : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn1
  have hsqrt : 0 < Real.sqrt (n : ℝ) := Real.sqrt_pos.2 hnreal
  have hsqrt_sq : Real.sqrt (n : ℝ) ^ 2 = (n : ℝ) :=
    Real.sq_sqrt hnreal.le
  have hmass := rivoalWeightMass22_pos n
  have henv : 0 ≤ rivoalHarmonicEnvelope22 n := by
    unfold rivoalHarmonicEnvelope22
    have : 0 ≤ Real.log (n : ℝ) :=
      Real.log_nonneg (by exact_mod_cast hn1)
    positivity
  have hcoef :
      0 ≤ rivoalHarmonicEnvelope22 n /
        (δ ^ 2 * (n : ℝ) ^ 4) := by positivity
  have hpoint : ∀ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k *
          |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant| ≤
        rivoalRealWeight22 n k *
          ((5 / (M : ℝ) + 8 * δ) +
            rivoalHarmonicEnvelope22 n *
              ((rivoalSaddleError22 n k) ^ 2 /
                (δ ^ 2 * (n : ℝ) ^ 4))) := by
    intro k hk
    exact mul_le_mul_of_nonneg_left
      (rivoalHarmonicError22_pointwise_le hδ hδle hM hn
        (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)) hMn)
      (rivoalRealWeight22_nonneg n k)
  have hsum :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant|) ≤
        (5 / (M : ℝ) + 8 * δ) * rivoalWeightMass22 n +
          (rivoalHarmonicEnvelope22 n /
            (δ ^ 2 * (n : ℝ) ^ 4)) *
            (∑ k ∈ Finset.range (n + 1),
              rivoalRealWeight22 n k * (rivoalSaddleError22 n k) ^ 2) := by
    calc
      _ ≤ ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            ((5 / (M : ℝ) + 8 * δ) +
              rivoalHarmonicEnvelope22 n *
                ((rivoalSaddleError22 n k) ^ 2 /
                  (δ ^ 2 * (n : ℝ) ^ 4))) :=
        Finset.sum_le_sum hpoint
      _ = _ := by
        calc
          _ = ∑ k ∈ Finset.range (n + 1),
              (rivoalRealWeight22 n k * (5 / (M : ℝ) + 8 * δ) +
                (rivoalHarmonicEnvelope22 n /
                  (δ ^ 2 * (n : ℝ) ^ 4)) *
                  (rivoalRealWeight22 n k *
                    (rivoalSaddleError22 n k) ^ 2)) := by
                apply Finset.sum_congr rfl
                intro k hk
                field_simp [hδ.ne', hnreal.ne']
          _ = (∑ k ∈ Finset.range (n + 1),
                rivoalRealWeight22 n k * (5 / (M : ℝ) + 8 * δ)) +
              ∑ k ∈ Finset.range (n + 1),
                (rivoalHarmonicEnvelope22 n /
                  (δ ^ 2 * (n : ℝ) ^ 4)) *
                  (rivoalRealWeight22 n k *
                    (rivoalSaddleError22 n k) ^ 2) := by
                rw [Finset.sum_add_distrib]
          _ = _ := by
                rw [← Finset.sum_mul, rivoalRealWeight22_sum_eq_mass]
                rw [← Finset.mul_sum]
                ring
  have hmoment := rivoalSaddleSecondMoment22_le hn1
  have hsum' :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            |rivoalRealHarmonicValue22 n k - eulerMascheroniConstant|) ≤
        (5 / (M : ℝ) + 8 * δ) * rivoalWeightMass22 n +
          (rivoalHarmonicEnvelope22 n /
            (δ ^ 2 * (n : ℝ) ^ 4)) *
            (81 * (n : ℝ) ^ 3 * Real.sqrt (n : ℝ) *
              rivoalWeightMass22 n) := by
    exact hsum.trans (add_le_add (le_refl _)
      (mul_le_mul_of_nonneg_left hmoment hcoef))
  rw [rivoalWeightedHarmonicError22]
  change _ / rivoalWeightMass22 n ≤ _
  apply (div_le_iff₀ hmass).2
  apply hsum'.trans_eq
  rw [← hsqrt_sq]
  rw [Real.sqrt_sq hsqrt.le]
  field_simp [hδ.ne', hsqrt.ne', hmass.ne']

theorem rivoalWeightedHarmonicError22_nonneg (n : ℕ) :
    0 ≤ rivoalWeightedHarmonicError22 n := by
  rw [rivoalWeightedHarmonicError22]
  apply div_nonneg
  · apply Finset.sum_nonneg
    intro k hk
    exact mul_nonneg (rivoalRealWeight22_nonneg n k) (abs_nonneg _)
  · exact (rivoalWeightMass22_pos n).le

theorem log_div_sqrt_nat_tendsto_zero22 :
    Tendsto
      (fun n : ℕ => Real.log (n : ℝ) / Real.sqrt (n : ℝ))
      atTop (nhds 0) := by
  have hreal :=
    (isLittleO_log_rpow_atTop
      (show (0 : ℝ) < 1 / 2 by norm_num)).tendsto_div_nhds_zero
  have hnat := hreal.comp tendsto_natCast_atTop_atTop
  simpa [Real.sqrt_eq_rpow] using hnat

theorem rivoalHarmonicEnvelope22_div_sqrt_tendsto_zero :
    Tendsto
      (fun n : ℕ =>
        rivoalHarmonicEnvelope22 n / Real.sqrt (n : ℝ))
      atTop (nhds 0) := by
  have hsqrtTop :
      Tendsto (fun n : ℕ => Real.sqrt (n : ℝ)) atTop atTop :=
    Real.tendsto_sqrt_atTop.comp tendsto_natCast_atTop_atTop
  have hconst :
      Tendsto (fun n : ℕ => (6 : ℝ) / Real.sqrt (n : ℝ))
        atTop (nhds 0) :=
    tendsto_const_nhds.div_atTop hsqrtTop
  have hfive : Tendsto (fun _ : ℕ => (5 : ℝ)) atTop (nhds 5) :=
    tendsto_const_nhds
  have hlog := hfive.mul log_div_sqrt_nat_tendsto_zero22
  simpa [rivoalHarmonicEnvelope22, add_div, mul_div_assoc] using
    hconst.add hlog

theorem rivoalHarmonicBadBound22_tendsto_zero
    {δ : ℝ} (_hδ : 0 < δ) :
    Tendsto
      (fun n : ℕ =>
        81 * rivoalHarmonicEnvelope22 n /
          (δ ^ 2 * Real.sqrt (n : ℝ)))
      atTop (nhds 0) := by
  have hc :
      Tendsto (fun _ : ℕ => (81 : ℝ) / δ ^ 2)
        atTop (nhds ((81 : ℝ) / δ ^ 2)) :=
    tendsto_const_nhds
  have h := hc.mul rivoalHarmonicEnvelope22_div_sqrt_tendsto_zero
  simpa [div_eq_mul_inv, mul_inv, mul_assoc, mul_left_comm, mul_comm] using h

/-- The analytic condition left open in `Problem22` follows from the finite
Stein identity and the quantitative saddle estimate above. -/
theorem rivoalHarmonicConcentrationClaim22 :
    RivoalHarmonicConcentrationClaim22 := by
  rw [RivoalHarmonicConcentrationClaim22, Metric.tendsto_atTop]
  intro ε hε
  let δ : ℝ := min (1 / 16) (ε / 32)
  have hδ : 0 < δ := by
    dsimp [δ]
    exact lt_min (by norm_num) (by positivity)
  have hδle : δ ≤ 1 / 16 := by
    dsimp [δ]
    exact min_le_left _ _
  have hδε : δ ≤ ε / 32 := by
    dsimp [δ]
    exact min_le_right _ _
  have h8δ : 8 * δ ≤ ε / 4 := by linarith
  obtain ⟨M, hMgt⟩ := exists_nat_gt (20 / ε)
  have hMreal : (0 : ℝ) < (M : ℝ) := by
    exact lt_trans (by positivity : (0 : ℝ) < 20 / ε) hMgt
  have hM : 1 ≤ M := by exact_mod_cast hMreal
  have h20 : 20 < (M : ℝ) * ε := (div_lt_iff₀ hε).mp hMgt
  have h5M : 5 / (M : ℝ) < ε / 4 := by
    apply (div_lt_iff₀ hMreal).2
    nlinarith
  have hbad := (tendsto_order.1
    (rivoalHarmonicBadBound22_tendsto_zero hδ)).2 (ε / 2) (by linarith)
  rw [eventually_atTop] at hbad
  obtain ⟨N, hN⟩ := hbad
  refine ⟨max N (max 16 (16 * M ^ 3)), ?_⟩
  intro n hnlarge
  have hnN : N ≤ n := by omega
  have hn16 : 16 ≤ n := by omega
  have hMnNat : 16 * M ^ 3 ≤ n := by omega
  have hMnCast : 16 * (M : ℝ) ^ 3 ≤ (n : ℝ) := by
    exact_mod_cast hMnNat
  have hnreal : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn16.trans' (by norm_num)
  have hMn : 16 * (M : ℝ) ^ 3 ≤ 3 * (n : ℝ) ^ 2 := by
    calc
      _ ≤ (n : ℝ) := hMnCast
      _ ≤ 3 * (n : ℝ) ^ 2 := by nlinarith [sq_nonneg ((n : ℝ) - 1)]
  have hquant := rivoalWeightedHarmonicError22_le hδ hδle hM hn16 hMn
  have hbadn := hN n hnN
  rw [Real.dist_eq, sub_zero,
    abs_of_nonneg (rivoalWeightedHarmonicError22_nonneg n)]
  exact hquant.trans_lt (by linarith)

/-- Unconditional closure of Problem 2.2. -/
theorem problem22_solved : Problem22Claim :=
  problem22_of_harmonic_concentration rivoalHarmonicConcentrationClaim22

end RamanujanChallenge.P22

end
