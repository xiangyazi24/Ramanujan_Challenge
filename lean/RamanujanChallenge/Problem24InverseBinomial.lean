import RamanujanChallenge.Problem24
import Mathlib.Analysis.Normed.Group.Tannery

/-! WZ proofs of the inverse-central-binomial sums of weight four. -/

set_option maxHeartbeats 100000

open scoped BigOperators Topology
open Filter Finset

namespace RamanujanChallenge

private theorem hasSum_of_diagonal_wz
    (F G : ℕ → ℕ → ℝ) (s : ℝ)
    (hwz : ∀ n k, F (n + 1) k - F n k = G n (k + 1) - G n k)
    (hsource : HasSum (F 0) s)
    (htarget : Summable (fun n => G n 0))
    (hFboundary : Tendsto (fun N => ∑ k ∈ range N, F N k) atTop (nhds 0))
    (hGboundary : Tendsto (fun N => ∑ n ∈ range N, G n N) atTop (nhds 0)) :
    HasSum (fun n => G n 0) s := by
  have hrectangle (N : ℕ) :
      (∑ n ∈ range N, G n 0) =
        (∑ k ∈ range N, F 0 k) - (∑ k ∈ range N, F N k) +
          ∑ n ∈ range N, G n N := by
    have hdouble :
        (∑ k ∈ range N, (F N k - F 0 k)) =
          ∑ n ∈ range N, (G n N - G n 0) := by
      calc
        (∑ k ∈ range N, (F N k - F 0 k)) =
            ∑ k ∈ range N, ∑ n ∈ range N,
              (F (n + 1) k - F n k) := by
                apply sum_congr rfl
                intro k hk
                exact (Finset.sum_range_sub (fun n => F n k) N).symm
        _ = ∑ k ∈ range N, ∑ n ∈ range N,
              (G n (k + 1) - G n k) := by
                apply sum_congr rfl
                intro k hk
                apply sum_congr rfl
                intro n hn
                exact hwz n k
        _ = ∑ n ∈ range N, ∑ k ∈ range N,
              (G n (k + 1) - G n k) := by
                rw [sum_comm]
        _ = ∑ n ∈ range N, (G n N - G n 0) := by
                apply sum_congr rfl
                intro n hn
                exact Finset.sum_range_sub (fun k => G n k) N
    rw [sum_sub_distrib, sum_sub_distrib] at hdouble
    linarith
  have hlimit :
      Tendsto (fun N => ∑ n ∈ range N, G n 0) atTop (nhds s) := by
    have hright :=
      (hsource.tendsto_sum_nat.sub hFboundary).add hGboundary
    simpa using
      hright.congr' (Eventually.of_forall fun N => (hrectangle N).symm)
  have htargetLimit := htarget.hasSum.tendsto_sum_nat
  have hvalue : ∑' n, G n 0 = s := tendsto_nhds_unique htargetLimit hlimit
  simpa [hvalue] using htarget.hasSum

private noncomputable def leshBase (n k : ℕ) : ℝ :=
  (-1 : ℝ) ^ k * (k.factorial : ℝ) * (n.factorial : ℝ) ^ 2 /
    ((2 * n + k + 1).factorial : ℝ)

private noncomputable def leshCoeff (n : ℕ) : ℝ :=
  (n.factorial : ℝ) ^ 2 / ((2 * n + 1).factorial : ℝ)

private theorem leshCoeff_nonneg (n : ℕ) : 0 ≤ leshCoeff n := by
  unfold leshCoeff
  positivity

private theorem leshCoeff_le (n : ℕ) :
    leshCoeff n ≤ 1 / (n + 1 : ℝ) ^ 2 := by
  have hchooseNat :
      Nat.choose (2 * n + 1) n * n.factorial * (n + 1).factorial =
        (2 * n + 1).factorial := by
    simpa [show 2 * n + 1 - n = n + 1 by omega] using
      Nat.choose_mul_factorial_mul_factorial
        (n := 2 * n + 1) (k := n) (by omega : n ≤ 2 * n + 1)
  have hchooseLower : n + 1 ≤ Nat.choose (2 * n + 1) n := by
    have hmono := Nat.choose_le_choose n (show n + 1 ≤ 2 * n + 1 by omega)
    simpa using hmono
  have hchooseR := congrArg (fun m : ℕ => (m : ℝ)) hchooseNat
  have hlowerR : (n + 1 : ℝ) ≤ Nat.choose (2 * n + 1) n := by
    exact_mod_cast hchooseLower
  rw [Nat.factorial_succ] at hchooseR
  push_cast at hchooseR
  have heq : leshCoeff n =
      1 / ((n + 1 : ℝ) * Nat.choose (2 * n + 1) n) := by
    unfold leshCoeff
    field_simp [Nat.factorial_ne_zero,
      Nat.choose_ne_zero (by omega : n ≤ 2 * n + 1)]
    nlinarith [hchooseR]
  rw [heq]
  apply one_div_le_one_div_of_le (by positivity)
  have hn1 : (0 : ℝ) ≤ n + 1 := by positivity
  simpa [pow_two] using mul_le_mul_of_nonneg_left hlowerR hn1

private theorem leshBase_norm_le (n k : ℕ) :
    ‖leshBase n k‖ ≤ leshCoeff n / (k + 1 : ℝ) := by
  have hchooseNat :
      Nat.choose (2 * n + k + 1) k * k.factorial *
          (2 * n + 1).factorial = (2 * n + k + 1).factorial := by
    simpa [show 2 * n + k + 1 - k = 2 * n + 1 by omega] using
      Nat.choose_mul_factorial_mul_factorial
        (n := 2 * n + k + 1) (k := k)
        (by omega : k ≤ 2 * n + k + 1)
  have hchooseLower : k + 1 ≤ Nat.choose (2 * n + k + 1) k := by
    have hmono := Nat.choose_le_choose k
      (show k + 1 ≤ 2 * n + k + 1 by omega)
    simpa using hmono
  have hchooseR := congrArg (fun m : ℕ => (m : ℝ)) hchooseNat
  have hlowerR : (k + 1 : ℝ) ≤ Nat.choose (2 * n + k + 1) k := by
    exact_mod_cast hchooseLower
  push_cast at hchooseR
  have heq : ‖leshBase n k‖ =
      leshCoeff n / Nat.choose (2 * n + k + 1) k := by
    unfold leshBase leshCoeff
    rw [Real.norm_eq_abs, abs_div, abs_mul, abs_mul, abs_pow,
      abs_pow, abs_neg, abs_one, one_pow,
      abs_of_nonneg (show 0 ≤ (k.factorial : ℝ) by positivity),
      abs_of_nonneg (show 0 ≤ (n.factorial : ℝ) by positivity),
      abs_of_nonneg
        (show 0 ≤ ((2 * n + k + 1).factorial : ℝ) by positivity)]
    field_simp [Nat.factorial_ne_zero,
      Nat.choose_ne_zero (by omega : k ≤ 2 * n + k + 1)]
    nlinarith [hchooseR]
  rw [heq]
  exact div_le_div_of_nonneg_left (leshCoeff_nonneg n) (by positivity) hlowerR

private noncomputable def leshFCoef (n k : ℕ) : ℝ :=
  1 / (n + k + 1 : ℝ) ^ 3 -
    harmonicSquare24 n / (n + k + 1 : ℝ)

private noncomputable def leshGCoef (n k : ℕ) : ℝ :=
  let N : ℝ := n + k + 1
  let A : ℝ := 3 * (n + 1 : ℝ) ^ 2 + (k : ℝ) ^ 2 +
    4 * k * (n + 1)
  1 / (2 * (2 * n + k + 2 : ℝ)) *
    (1 / N ^ 2 + A / N ^ 4 - A * harmonicSquare24 n / N ^ 2)

private noncomputable def leshF (n k : ℕ) : ℝ :=
  leshBase n k * leshFCoef n k

private noncomputable def leshG (n k : ℕ) : ℝ :=
  leshBase n k * leshGCoef n k

private theorem leshBase_succ_n (n k : ℕ) :
    leshBase (n + 1) k = leshBase n k *
      ((n + 1 : ℝ) ^ 2 /
        ((2 * n + k + 2 : ℝ) * (2 * n + k + 3 : ℝ))) := by
  unfold leshBase
  have hfn : (2 * (n + 1) + k + 1).factorial =
      (2 * n + k + 3) * (2 * n + k + 2) *
        (2 * n + k + 1).factorial := by
    rw [show 2 * (n + 1) + k + 1 = (2 * n + k + 2) + 1 by omega,
      Nat.factorial_succ,
      show 2 * n + k + 2 = (2 * n + k + 1) + 1 by omega,
      Nat.factorial_succ]
    ring
  rw [hfn, Nat.factorial_succ]
  push_cast
  field_simp [Nat.factorial_ne_zero]

private theorem leshBase_succ_k (n k : ℕ) :
    leshBase n (k + 1) = leshBase n k *
      (-(k + 1 : ℝ) / (2 * n + k + 2 : ℝ)) := by
  unfold leshBase
  have hfg : (2 * n + (k + 1) + 1).factorial =
      (2 * n + k + 2) * (2 * n + k + 1).factorial := by
    rw [show 2 * n + (k + 1) + 1 = (2 * n + k + 1) + 1 by omega,
      Nat.factorial_succ]
  rw [hfg, Nat.factorial_succ, pow_succ]
  push_cast
  field_simp [Nat.factorial_ne_zero]

private theorem lesh_rational_wz (n k : ℕ) :
    ((n + 1 : ℝ) ^ 2 /
        ((2 * n + k + 2 : ℝ) * (2 * n + k + 3 : ℝ))) *
          leshFCoef (n + 1) k - leshFCoef n k =
      (-(k + 1 : ℝ) / (2 * n + k + 2 : ℝ)) *
          leshGCoef n (k + 1) - leshGCoef n k := by
  unfold leshFCoef leshGCoef
  rw [harmonicSquare24_succ]
  push_cast
  field_simp
  ring

private theorem lesh_wz (n k : ℕ) :
    leshF (n + 1) k - leshF n k = leshG n (k + 1) - leshG n k := by
  unfold leshF leshG
  rw [leshBase_succ_n, leshBase_succ_k]
  calc
    leshBase n k *
          ((n + 1 : ℝ) ^ 2 /
            ((2 * n + k + 2 : ℝ) * (2 * n + k + 3 : ℝ))) *
          leshFCoef (n + 1) k - leshBase n k * leshFCoef n k =
        leshBase n k *
          (((n + 1 : ℝ) ^ 2 /
              ((2 * n + k + 2 : ℝ) * (2 * n + k + 3 : ℝ))) *
            leshFCoef (n + 1) k - leshFCoef n k) := by ring
    _ = leshBase n k *
          ((-(k + 1 : ℝ) / (2 * n + k + 2 : ℝ)) *
            leshGCoef n (k + 1) -
              leshGCoef n k) := by rw [lesh_rational_wz]
    _ = leshBase n k * (-(k + 1 : ℝ) / (2 * n + k + 2 : ℝ)) *
          leshGCoef n (k + 1) - leshBase n k * leshGCoef n k := by ring

private theorem leshF_zero (k : ℕ) :
    leshF 0 k = (-1 : ℝ) ^ k / (k + 1 : ℝ) ^ 4 := by
  unfold leshF leshBase leshFCoef
  rw [show (2 * 0 + k + 1).factorial = (k + 1) * k.factorial by
    simpa [Nat.add_comm] using Nat.factorial_succ k]
  push_cast
  simp [harmonicSquare24]
  field_simp [Nat.factorial_ne_zero]

private theorem leshBase_zero_eq (n : ℕ) :
    leshBase n 0 =
      2 * (n + 1 : ℝ) * inverseCentralCoefficient24 n := by
  have hchooseNat :
      Nat.choose (2 * (n + 1)) (n + 1) * (n + 1).factorial *
          (n + 1).factorial = (2 * (n + 1)).factorial := by
    convert Nat.add_choose_mul_factorial_mul_factorial (n + 1) (n + 1) using 1 <;>
      ring
  have hchooseR := congrArg (fun m : ℕ => (m : ℝ)) hchooseNat
  have htwice : 2 * (n + 1) = 2 * n + 2 := by omega
  rw [htwice] at hchooseR
  have hfact : (2 * n + 2).factorial =
      (2 * n + 2) * (2 * n + 1).factorial := by
    rw [show 2 * n + 2 = (2 * n + 1) + 1 by omega, Nat.factorial_succ]
  rw [hfact, Nat.factorial_succ] at hchooseR
  push_cast at hchooseR
  have hrel :
      (Nat.choose (2 * n + 2) (n + 1) : ℝ) * (n + 1) *
          (n.factorial : ℝ) ^ 2 =
        2 * ((2 * n + 1).factorial : ℝ) := by
    have hn : (0 : ℝ) < n + 1 := by positivity
    ring_nf at hchooseR ⊢
    nlinarith [hchooseR]
  have hchooseNe :
      (Nat.choose (2 * n + 2) (n + 1) : ℝ) ≠ 0 := by
    exact_mod_cast Nat.choose_ne_zero (by omega : n + 1 ≤ 2 * n + 2)
  unfold leshBase inverseCentralCoefficient24
  simp only [pow_zero, Nat.factorial_zero, Nat.cast_one, one_mul, add_zero]
  rw [show 2 * (n + 1) = 2 * n + 2 by omega]
  field_simp [Nat.factorial_ne_zero, hchooseNe]
  ring_nf at hrel ⊢
  nlinarith [hrel]

private theorem leshG_zero (n : ℕ) :
    leshG n 0 = leshchinerWeightFourTerm24 n := by
  unfold leshG
  rw [leshBase_zero_eq]
  unfold leshGCoef leshchinerWeightFourTerm24
    inverseCentralFourthCoefficient24
  simp only [Nat.cast_zero, add_zero, mul_zero]
  field_simp
  ring

private theorem leshCoeff_summable : Summable leshCoeff := by
  apply shifted_zeta_two_hasSum.summable.of_norm_bounded
  intro n
  rw [Real.norm_eq_abs, abs_of_nonneg (leshCoeff_nonneg n)]
  exact leshCoeff_le n

private theorem leshCoeff_tendsto_zero :
    Tendsto leshCoeff atTop (nhds 0) :=
  leshCoeff_summable.tendsto_atTop_zero

private theorem leshF_norm_le (n k : ℕ) :
    ‖leshF n k‖ ≤
      (1 + Real.pi ^ 2 / 6) * leshCoeff n / (k + 1 : ℝ) ^ 2 := by
  have hN : (0 : ℝ) < n + k + 1 := by positivity
  have hk1 : (0 : ℝ) < k + 1 := by positivity
  have hkN : (k + 1 : ℝ) ≤ n + k + 1 := by
    linarith
  have hH : 0 ≤ harmonicSquare24 n := by
    unfold harmonicSquare24
    positivity
  have hpow : (k + 1 : ℝ) ≤ (n + k + 1 : ℝ) ^ 3 := by
    calc
      (k + 1 : ℝ) ≤ n + k + 1 := hkN
      _ = (n + k + 1 : ℝ) * 1 := by ring
      _ ≤ (n + k + 1 : ℝ) * (n + k + 1) ^ 2 := by
        apply mul_le_mul_of_nonneg_left _ hN.le
        nlinarith [sq_nonneg ((n + k + 1 : ℝ) - 1)]
      _ = (n + k + 1 : ℝ) ^ 3 := by ring
  have hcoef : ‖leshFCoef n k‖ ≤
      (1 + Real.pi ^ 2 / 6) / (k + 1 : ℝ) := by
    unfold leshFCoef
    rw [Real.norm_eq_abs]
    calc
      |1 / (n + k + 1 : ℝ) ^ 3 -
          harmonicSquare24 n / (n + k + 1 : ℝ)| ≤
          1 / (n + k + 1 : ℝ) ^ 3 +
            harmonicSquare24 n / (n + k + 1 : ℝ) := by
              rw [abs_sub_le_iff]
              constructor <;> nlinarith [show
                0 ≤ 1 / (n + k + 1 : ℝ) ^ 3 by positivity,
                show 0 ≤ harmonicSquare24 n / (n + k + 1 : ℝ) by positivity]
      _ ≤ 1 / (k + 1 : ℝ) + (Real.pi ^ 2 / 6) / (k + 1 : ℝ) := by
        apply add_le_add
        · exact one_div_le_one_div_of_le hk1 hpow
        · calc
            harmonicSquare24 n / (n + k + 1 : ℝ) ≤
                (Real.pi ^ 2 / 6) / (n + k + 1 : ℝ) := by
                  exact div_le_div_of_nonneg_right
                    (harmonicSquare24_le_zeta_two n) hN.le
            _ ≤ (Real.pi ^ 2 / 6) / (k + 1 : ℝ) :=
              div_le_div_of_nonneg_left (by positivity) hk1 hkN
      _ = (1 + Real.pi ^ 2 / 6) / (k + 1 : ℝ) := by ring
  unfold leshF
  rw [norm_mul]
  calc
    ‖leshBase n k‖ * ‖leshFCoef n k‖ ≤
        (leshCoeff n / (k + 1 : ℝ)) *
          ((1 + Real.pi ^ 2 / 6) / (k + 1 : ℝ)) :=
      mul_le_mul (leshBase_norm_le n k) hcoef (norm_nonneg _)
        (div_nonneg (leshCoeff_nonneg n) hk1.le)
    _ = (1 + Real.pi ^ 2 / 6) * leshCoeff n /
          (k + 1 : ℝ) ^ 2 := by
      field_simp [ne_of_gt hk1]

private theorem leshG_norm_le (n k : ℕ) :
    ‖leshG n k‖ ≤
      (4 + 3 * (Real.pi ^ 2 / 6)) * leshCoeff n /
        (2 * (k + 1 : ℝ) ^ 2) := by
  let N : ℝ := n + k + 1
  let A : ℝ := 3 * (n + 1 : ℝ) ^ 2 + (k : ℝ) ^ 2 +
    4 * k * (n + 1)
  have hN : 1 ≤ N := by
    unfold N
    linarith
  have hNpos : 0 < N := lt_of_lt_of_le (by norm_num) hN
  have hk1 : (0 : ℝ) < k + 1 := by positivity
  have hkN : (k + 1 : ℝ) ≤ N := by
    unfold N
    linarith
  have hA : 0 ≤ A := by
    unfold A
    positivity
  have hAle : A ≤ 3 * N ^ 2 := by
    unfold A N
    nlinarith [mul_nonneg (show (0 : ℝ) ≤ n + 1 by positivity)
      (show (0 : ℝ) ≤ k by positivity)]
  have hN2 : (1 : ℝ) ≤ N ^ 2 := by nlinarith [sq_nonneg (N - 1)]
  have hAdiv : A / N ^ 2 ≤ 3 := by
    rw [div_le_iff₀ (sq_pos_of_pos hNpos)]
    exact hAle
  have hH : 0 ≤ harmonicSquare24 n := by
    unfold harmonicSquare24
    positivity
  have hbracket :
      ‖1 / N ^ 2 + A / N ^ 4 - A * harmonicSquare24 n / N ^ 2‖ ≤
        4 + 3 * (Real.pi ^ 2 / 6) := by
    rw [Real.norm_eq_abs]
    calc
      |1 / N ^ 2 + A / N ^ 4 - A * harmonicSquare24 n / N ^ 2| ≤
          (1 / N ^ 2 + A / N ^ 4) +
            A * harmonicSquare24 n / N ^ 2 := by
              rw [abs_sub_le_iff]
              constructor <;> nlinarith [show 0 ≤ 1 / N ^ 2 by positivity,
                show 0 ≤ A / N ^ 4 by positivity,
                show 0 ≤ A * harmonicSquare24 n / N ^ 2 by positivity]
      _ ≤ 1 + 3 + 3 * (Real.pi ^ 2 / 6) := by
        gcongr
        · simpa using one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < 1) hN2
        · calc
            A / N ^ 4 = (A / N ^ 2) / N ^ 2 := by field_simp
            _ ≤ 3 / N ^ 2 := by gcongr
            _ ≤ 3 := by
              rw [div_le_iff₀ (sq_pos_of_pos hNpos)]
              nlinarith
        · calc
            A * harmonicSquare24 n / N ^ 2 =
                (A / N ^ 2) * harmonicSquare24 n := by field_simp
            _ ≤ 3 * harmonicSquare24 n :=
              mul_le_mul_of_nonneg_right hAdiv hH
            _ ≤ 3 * (Real.pi ^ 2 / 6) := by
              gcongr
              exact harmonicSquare24_le_zeta_two n
      _ = 4 + 3 * (Real.pi ^ 2 / 6) := by ring
  have hden : (k + 1 : ℝ) ≤ 2 * n + k + 2 := by
    linarith
  have hcoef : ‖leshGCoef n k‖ ≤
      (4 + 3 * (Real.pi ^ 2 / 6)) / (2 * (k + 1 : ℝ)) := by
    unfold leshGCoef
    change ‖1 / (2 * (2 * n + k + 2 : ℝ)) *
      (1 / N ^ 2 + A / N ^ 4 - A * harmonicSquare24 n / N ^ 2)‖ ≤ _
    rw [norm_mul, Real.norm_eq_abs,
      abs_of_nonneg (by positivity : (0 : ℝ) ≤ 1 / (2 * (2 * n + k + 2 : ℝ)))]
    calc
      1 / (2 * (2 * n + k + 2 : ℝ)) *
          ‖1 / N ^ 2 + A / N ^ 4 - A * harmonicSquare24 n / N ^ 2‖ ≤
          1 / (2 * (2 * n + k + 2 : ℝ)) *
            (4 + 3 * (Real.pi ^ 2 / 6)) :=
        mul_le_mul_of_nonneg_left hbracket (by positivity)
      _ ≤ (4 + 3 * (Real.pi ^ 2 / 6)) / (2 * (k + 1 : ℝ)) := by
        have hinv := one_div_le_one_div_of_le
          (show (0 : ℝ) < 2 * (k + 1 : ℝ) by positivity)
          (show 2 * (k + 1 : ℝ) ≤ 2 * (2 * n + k + 2 : ℝ) by nlinarith)
        calc
          1 / (2 * (2 * n + k + 2 : ℝ)) *
              (4 + 3 * (Real.pi ^ 2 / 6)) ≤
              1 / (2 * (k + 1 : ℝ)) *
                (4 + 3 * (Real.pi ^ 2 / 6)) :=
            mul_le_mul_of_nonneg_right hinv (by positivity)
          _ = (4 + 3 * (Real.pi ^ 2 / 6)) / (2 * (k + 1 : ℝ)) := by ring
  unfold leshG
  rw [norm_mul]
  calc
    ‖leshBase n k‖ * ‖leshGCoef n k‖ ≤
        (leshCoeff n / (k + 1 : ℝ)) *
          ((4 + 3 * (Real.pi ^ 2 / 6)) / (2 * (k + 1 : ℝ))) :=
      mul_le_mul (leshBase_norm_le n k) hcoef (norm_nonneg _)
        (div_nonneg (leshCoeff_nonneg n) hk1.le)
    _ = (4 + 3 * (Real.pi ^ 2 / 6)) * leshCoeff n /
          (2 * (k + 1 : ℝ) ^ 2) := by
      field_simp [ne_of_gt hk1]

private noncomputable def bbbBase (n k : ℕ) : ℝ :=
  (n.factorial : ℝ) ^ 4 * (k.factorial : ℝ) ^ 2 /
    (((2 * n).factorial : ℝ) * ((n + k + 1).factorial : ℝ) ^ 2)

private noncomputable def bbbCoeff (n : ℕ) : ℝ :=
  (n.factorial : ℝ) ^ 2 / ((2 * n).factorial : ℝ)

private theorem bbbCoeff_nonneg (n : ℕ) : 0 ≤ bbbCoeff n := by
  unfold bbbCoeff
  positivity

private theorem bbbCoeff_succ (n : ℕ) :
    bbbCoeff (n + 1) = bbbCoeff n *
      ((n + 1 : ℝ) ^ 2 / ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ))) := by
  unfold bbbCoeff
  have htwo : (2 * (n + 1)).factorial =
      (2 * n + 2) * (2 * n + 1) * (2 * n).factorial := by
    rw [show 2 * (n + 1) = (2 * n + 1) + 1 by omega,
      Nat.factorial_succ, Nat.factorial_succ]
    ring
  rw [htwo, Nat.factorial_succ]
  push_cast
  field_simp [Nat.factorial_ne_zero]

private theorem bbbCoeff_le_one (n : ℕ) : bbbCoeff n ≤ 1 := by
  have hchooseNat : Nat.choose (2 * n) n * n.factorial * n.factorial =
      (2 * n).factorial := by
    simpa [show 2 * n - n = n by omega] using
      Nat.choose_mul_factorial_mul_factorial
        (n := 2 * n) (k := n) (by omega : n ≤ 2 * n)
  have hchooseR := congrArg (fun m : ℕ => (m : ℝ)) hchooseNat
  push_cast at hchooseR
  unfold bbbCoeff
  rw [div_le_one (by positivity : (0 : ℝ) < (2 * n).factorial)]
  have hchoose : (1 : ℝ) ≤ Nat.choose (2 * n) n := by
    exact_mod_cast Nat.choose_pos (by omega : n ≤ 2 * n)
  have hchoose0 : (0 : ℝ) ≤ Nat.choose (2 * n) n :=
    le_trans (by norm_num) hchoose
  nlinarith [mul_nonneg hchoose0 (sq_nonneg (n.factorial : ℝ))]

private theorem bbbCoeff_summable : Summable bbbCoeff := by
  apply summable_of_ratio_norm_eventually_le (r := (1 / 2 : ℝ)) (by norm_num)
  exact Eventually.of_forall fun n => by
    rw [bbbCoeff_succ]
    simp only [norm_mul, Real.norm_eq_abs,
      abs_of_nonneg (bbbCoeff_nonneg n)]
    have hratio :
        (n + 1 : ℝ) ^ 2 /
            ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ)) ≤ 1 / 2 := by
      rw [div_le_iff₀ (by positivity :
        (0 : ℝ) < (2 * n + 1 : ℝ) * (2 * n + 2 : ℝ))]
      nlinarith [show (0 : ℝ) ≤ n by positivity]
    rw [abs_of_nonneg (by positivity :
      (0 : ℝ) ≤ (n + 1 : ℝ) ^ 2 /
        ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ)))]
    simpa [mul_comm] using
      mul_le_mul_of_nonneg_left hratio (bbbCoeff_nonneg n)

private theorem bbbCoeff_tendsto_zero :
    Tendsto bbbCoeff atTop (nhds 0) :=
  bbbCoeff_summable.tendsto_atTop_zero

private theorem bbbBase_norm_le (n k : ℕ) :
    ‖bbbBase n k‖ ≤ bbbCoeff n / (k + 1 : ℝ) ^ 2 := by
  have hchooseNat :
      Nat.choose (n + k + 1) n * n.factorial * (k + 1).factorial =
        (n + k + 1).factorial := by
    simpa [show n + k + 1 - n = k + 1 by omega] using
      Nat.choose_mul_factorial_mul_factorial
        (n := n + k + 1) (k := n) (by omega : n ≤ n + k + 1)
  have hchooseR := congrArg (fun m : ℕ => (m : ℝ)) hchooseNat
  rw [Nat.factorial_succ] at hchooseR
  push_cast at hchooseR
  have hsq := congrArg (fun x : ℝ => x ^ 2) hchooseR
  have hchoose : (1 : ℝ) ≤ Nat.choose (n + k + 1) n := by
    exact_mod_cast Nat.choose_pos (by omega : n ≤ n + k + 1)
  have heq : ‖bbbBase n k‖ =
      bbbCoeff n /
        ((Nat.choose (n + k + 1) n : ℝ) ^ 2 * (k + 1 : ℝ) ^ 2) := by
    unfold bbbBase bbbCoeff
    rw [Real.norm_eq_abs, abs_div, abs_mul, abs_pow, abs_pow,
      abs_mul, abs_pow,
      abs_of_nonneg (show 0 ≤ (n.factorial : ℝ) by positivity),
      abs_of_nonneg (show 0 ≤ (k.factorial : ℝ) by positivity),
      abs_of_nonneg (show 0 ≤ ((2 * n).factorial : ℝ) by positivity),
      abs_of_nonneg
        (show 0 ≤ ((n + k + 1).factorial : ℝ) by positivity)]
    field_simp [Nat.factorial_ne_zero,
      Nat.choose_ne_zero (by omega : n ≤ n + k + 1)]
    ring_nf at hsq ⊢
    nlinarith [hsq]
  rw [heq]
  apply div_le_div_of_nonneg_left (bbbCoeff_nonneg n) (by positivity)
  have hc2 : (1 : ℝ) ≤ (Nat.choose (n + k + 1) n : ℝ) ^ 2 := by
    nlinarith [hchoose,
      sq_nonneg ((Nat.choose (n + k + 1) n : ℝ) - 1)]
  calc
    (k + 1 : ℝ) ^ 2 = 1 * (k + 1 : ℝ) ^ 2 := by ring
    _ ≤ (Nat.choose (n + k + 1) n : ℝ) ^ 2 * (k + 1 : ℝ) ^ 2 :=
      mul_le_mul_of_nonneg_right hc2 (sq_nonneg _)

private noncomputable def bbbFCoef (n k : ℕ) : ℝ :=
  harmonicSquare24 (n + k + 1) - harmonicSquare24 k -
    4 * harmonicSquare24 n

private noncomputable def bbbGCoef (n k : ℕ) : ℝ :=
  (n + 1 : ℝ) * (3 * n + 3 + 2 * k : ℝ) /
      ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ)) *
    bbbFCoef n k

private noncomputable def bbbF (n k : ℕ) : ℝ :=
  bbbBase n k * bbbFCoef n k

private noncomputable def bbbG (n k : ℕ) : ℝ :=
  bbbBase n k * bbbGCoef n k

private theorem bbbBase_succ_n (n k : ℕ) :
    bbbBase (n + 1) k = bbbBase n k *
      ((n + 1 : ℝ) ^ 4 /
        ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ) *
          (n + k + 2 : ℝ) ^ 2)) := by
  unfold bbbBase
  have htwo : (2 * (n + 1)).factorial =
      (2 * n + 2) * (2 * n + 1) * (2 * n).factorial := by
    rw [show 2 * (n + 1) = (2 * n + 1) + 1 by omega,
      Nat.factorial_succ,
      show 2 * n + 1 = 2 * n + 1 by rfl,
      Nat.factorial_succ]
    ring
  have hsum : (n + 1 + k + 1).factorial =
      (n + k + 2) * (n + k + 1).factorial := by
    rw [show n + 1 + k + 1 = (n + k + 1) + 1 by omega,
      Nat.factorial_succ]
  rw [htwo, hsum, Nat.factorial_succ]
  push_cast
  field_simp [Nat.factorial_ne_zero]

private theorem bbbBase_succ_k (n k : ℕ) :
    bbbBase n (k + 1) = bbbBase n k *
      ((k + 1 : ℝ) ^ 2 / (n + k + 2 : ℝ) ^ 2) := by
  unfold bbbBase
  have hsum : (n + (k + 1) + 1).factorial =
      (n + k + 2) * (n + k + 1).factorial := by
    rw [show n + (k + 1) + 1 = (n + k + 1) + 1 by omega,
      Nat.factorial_succ]
  rw [hsum, Nat.factorial_succ]
  push_cast
  field_simp [Nat.factorial_ne_zero]

private theorem bbb_rational_wz (n k : ℕ) :
    ((n + 1 : ℝ) ^ 4 /
        ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ) *
          (n + k + 2 : ℝ) ^ 2)) * bbbFCoef (n + 1) k -
      bbbFCoef n k =
      ((k + 1 : ℝ) ^ 2 / (n + k + 2 : ℝ) ^ 2) *
        bbbGCoef n (k + 1) - bbbGCoef n k := by
  unfold bbbGCoef bbbFCoef
  simp only [show n + 1 + k + 1 = n + k + 1 + 1 by omega,
    show n + (k + 1) + 1 = n + k + 1 + 1 by omega]
  rw [harmonicSquare24_succ (n + k + 1),
    harmonicSquare24_succ k, harmonicSquare24_succ n]
  push_cast
  field_simp
  ring

private theorem bbb_wz (n k : ℕ) :
    bbbF (n + 1) k - bbbF n k = bbbG n (k + 1) - bbbG n k := by
  unfold bbbF bbbG
  rw [bbbBase_succ_n, bbbBase_succ_k]
  calc
    bbbBase n k *
          ((n + 1 : ℝ) ^ 4 /
            ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ) *
              (n + k + 2 : ℝ) ^ 2)) * bbbFCoef (n + 1) k -
        bbbBase n k * bbbFCoef n k =
      bbbBase n k *
        (((n + 1 : ℝ) ^ 4 /
          ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ) *
            (n + k + 2 : ℝ) ^ 2)) * bbbFCoef (n + 1) k -
          bbbFCoef n k) := by ring
    _ = bbbBase n k *
        (((k + 1 : ℝ) ^ 2 / (n + k + 2 : ℝ) ^ 2) *
          bbbGCoef n (k + 1) - bbbGCoef n k) := by rw [bbb_rational_wz]
    _ = bbbBase n k *
          ((k + 1 : ℝ) ^ 2 / (n + k + 2 : ℝ) ^ 2) *
            bbbGCoef n (k + 1) - bbbBase n k * bbbGCoef n k := by ring

private theorem bbbF_zero (k : ℕ) :
    bbbF 0 k = 1 / (k + 1 : ℝ) ^ 4 := by
  unfold bbbF bbbBase bbbFCoef
  rw [harmonicSquare24_succ]
  rw [show (0 + k + 1).factorial = (k + 1) * k.factorial by
    simpa [Nat.add_comm] using Nat.factorial_succ k]
  push_cast
  simp [harmonicSquare24]
  field_simp [Nat.factorial_ne_zero]

private theorem bbbG_zero (n : ℕ) :
    bbbG n 0 = bbbWeightFourTerm24 n := by
  have hbase : bbbBase n 0 =
      leshBase n 0 * (2 * n + 1 : ℝ) / (n + 1 : ℝ) ^ 2 := by
    unfold bbbBase leshBase
    simp only [Nat.factorial_zero, Nat.cast_one, one_pow, mul_one,
      pow_zero, one_mul, add_zero]
    rw [show (2 * n + 1).factorial =
      (2 * n + 1) * (2 * n).factorial by
        rw [show 2 * n + 1 = 2 * n + 1 by rfl, Nat.factorial_succ]]
    rw [Nat.factorial_succ]
    push_cast
    field_simp [Nat.factorial_ne_zero]
  unfold bbbG
  rw [hbase, leshBase_zero_eq]
  unfold bbbGCoef bbbFCoef bbbWeightFourTerm24
    inverseCentralFourthCoefficient24
  rw [harmonicSquare24_succ]
  simp only [harmonicSquare24, sum_range_zero, Nat.cast_zero, sub_zero]
  field_simp
  ring

private theorem bbbF_norm_le (n k : ℕ) :
    ‖bbbF n k‖ ≤
      (Real.pi ^ 2) * bbbCoeff n / (k + 1 : ℝ) ^ 2 := by
  have hnonneg (m : ℕ) : 0 ≤ harmonicSquare24 m := by
    unfold harmonicSquare24
    positivity
  have hcoef : ‖bbbFCoef n k‖ ≤ Real.pi ^ 2 := by
    unfold bbbFCoef
    rw [Real.norm_eq_abs]
    calc
      |harmonicSquare24 (n + k + 1) - harmonicSquare24 k -
          4 * harmonicSquare24 n| ≤
          |harmonicSquare24 (n + k + 1) - harmonicSquare24 k| +
            |4 * harmonicSquare24 n| := abs_sub _ _
      _ ≤ (|harmonicSquare24 (n + k + 1)| + |harmonicSquare24 k|) +
            |4 * harmonicSquare24 n| := by gcongr; exact abs_sub _ _
      _ = harmonicSquare24 (n + k + 1) + harmonicSquare24 k +
            4 * harmonicSquare24 n := by
              rw [abs_of_nonneg (hnonneg _), abs_of_nonneg (hnonneg _),
                abs_mul, abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 4),
                abs_of_nonneg (hnonneg _)]
      _ ≤ Real.pi ^ 2 / 6 + Real.pi ^ 2 / 6 +
            4 * (Real.pi ^ 2 / 6) := by
              gcongr <;> apply harmonicSquare24_le_zeta_two
      _ = Real.pi ^ 2 := by ring
  unfold bbbF
  rw [norm_mul]
  calc
    ‖bbbBase n k‖ * ‖bbbFCoef n k‖ ≤
        (bbbCoeff n / (k + 1 : ℝ) ^ 2) * Real.pi ^ 2 :=
      mul_le_mul (bbbBase_norm_le n k) hcoef (norm_nonneg _)
        (div_nonneg (bbbCoeff_nonneg n) (sq_nonneg _))
    _ = Real.pi ^ 2 * bbbCoeff n / (k + 1 : ℝ) ^ 2 := by ring

private theorem bbbG_norm_le (n k : ℕ) :
    ‖bbbG n k‖ ≤
      (3 * Real.pi ^ 2) * bbbCoeff n / (k + 1 : ℝ) := by
  have hk1 : (0 : ℝ) < k + 1 := by positivity
  have hnonneg (m : ℕ) : 0 ≤ harmonicSquare24 m := by
    unfold harmonicSquare24
    positivity
  have hFcoef : ‖bbbFCoef n k‖ ≤ Real.pi ^ 2 := by
    unfold bbbFCoef
    rw [Real.norm_eq_abs]
    calc
      |harmonicSquare24 (n + k + 1) - harmonicSquare24 k -
          4 * harmonicSquare24 n| ≤
          |harmonicSquare24 (n + k + 1) - harmonicSquare24 k| +
            |4 * harmonicSquare24 n| := abs_sub _ _
      _ ≤ (|harmonicSquare24 (n + k + 1)| + |harmonicSquare24 k|) +
            |4 * harmonicSquare24 n| := by gcongr; exact abs_sub _ _
      _ = harmonicSquare24 (n + k + 1) + harmonicSquare24 k +
            4 * harmonicSquare24 n := by
              rw [abs_of_nonneg (hnonneg _), abs_of_nonneg (hnonneg _),
                abs_mul, abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 4),
                abs_of_nonneg (hnonneg _)]
      _ ≤ Real.pi ^ 2 / 6 + Real.pi ^ 2 / 6 +
            4 * (Real.pi ^ 2 / 6) := by
              gcongr <;> apply harmonicSquare24_le_zeta_two
      _ = Real.pi ^ 2 := by ring
  have hprefNonneg :
      0 ≤ (n + 1 : ℝ) * (3 * n + 3 + 2 * k : ℝ) /
        ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ)) := by positivity
  have hpref :
      (n + 1 : ℝ) * (3 * n + 3 + 2 * k : ℝ) /
          ((2 * n + 1 : ℝ) * (2 * n + 2 : ℝ)) ≤
        3 * (k + 1 : ℝ) := by
    rw [div_le_iff₀ (by positivity :
      (0 : ℝ) < (2 * n + 1 : ℝ) * (2 * n + 2 : ℝ))]
    nlinarith [mul_nonneg (show (0 : ℝ) ≤ n by positivity)
      (show (0 : ℝ) ≤ k by positivity)]
  have hGcoef : ‖bbbGCoef n k‖ ≤ 3 * (k + 1 : ℝ) * Real.pi ^ 2 := by
    unfold bbbGCoef
    rw [norm_mul, Real.norm_eq_abs, abs_of_nonneg hprefNonneg]
    exact mul_le_mul hpref hFcoef (norm_nonneg _) (by positivity)
  unfold bbbG
  rw [norm_mul]
  calc
    ‖bbbBase n k‖ * ‖bbbGCoef n k‖ ≤
        (bbbCoeff n / (k + 1 : ℝ) ^ 2) *
          (3 * (k + 1 : ℝ) * Real.pi ^ 2) :=
      mul_le_mul (bbbBase_norm_le n k) hGcoef (norm_nonneg _)
        (div_nonneg (bbbCoeff_nonneg n) (sq_nonneg _))
    _ = (3 * Real.pi ^ 2) * bbbCoeff n / (k + 1 : ℝ) := by
      field_simp [ne_of_gt hk1]

private theorem alternating_zeta_four_hasSum24 :
    HasSum (fun k : ℕ => (-1 : ℝ) ^ k / (k + 1 : ℝ) ^ 4)
      ((7 / 8 : ℝ) * (Real.pi ^ 4 / 90)) := by
  let f : ℕ → ℝ := fun k => 1 / (k + 1 : ℝ) ^ 4
  let a : ℕ → ℝ := fun k => (-1 : ℝ) ^ k * f k
  have hf : HasSum f (Real.pi ^ 4 / 90) := by
    simpa [f] using shifted_zeta_four_hasSum24
  have haSummable : Summable a := by
    apply hf.summable.of_norm_bounded
    intro k
    simp [a, f, Real.norm_eq_abs, abs_of_pos (by positivity : (0 : ℝ) < k + 1)]
  have ha := haSummable.hasSum
  have hpaired := ha.pair_consecutive24
  have hodd :
      HasSum (fun m => f (2 * m + 1))
        ((1 / 16 : ℝ) * (Real.pi ^ 4 / 90)) := by
    convert hf.mul_left (1 / 16 : ℝ) using 1
    funext m
    simp [f]
    field_simp
    ring
  have hdiff :
      HasSum (fun m => f (2 * m) - f (2 * m + 1))
        ((7 / 8 : ℝ) * (Real.pi ^ 4 / 90)) := by
    convert hf.pair_consecutive24.sub (hodd.mul_left 2) using 1
    · funext m
      ring
    · ring
  have hpairEq :
      (fun m => a (2 * m) + a (2 * m + 1)) =
        fun m => f (2 * m) - f (2 * m + 1) := by
    funext m
    simp [a, pow_add, f]
    ring
  rw [hpairEq] at hpaired
  have hvalue : (∑' k, a k) = (7 / 8 : ℝ) * (Real.pi ^ 4 / 90) :=
    hpaired.unique hdiff
  rw [hvalue] at ha
  convert ha using 1
  funext k
  simp [a, f, div_eq_mul_inv]

private theorem inverseBinomial_tendsto_zero_of_norm24
    {u : ℕ → ℝ} (h : Tendsto (fun n => ‖u n‖) atTop (nhds 0)) :
    Tendsto u atTop (nhds 0) := by
  rw [NormedAddGroup.tendsto_nhds_zero]
  intro ε hε
  simpa using (NormedAddGroup.tendsto_nhds_zero.mp h ε hε)

private theorem leshF_boundary :
    Tendsto (fun N => ∑ k ∈ range N, leshF N k) atTop (nhds 0) := by
  let f : ℕ → ℕ → ℝ := fun N k => if k < N then leshF N k else 0
  let bound : ℕ → ℝ := fun k =>
    (1 + Real.pi ^ 2 / 6) * (1 / (k + 1 : ℝ) ^ 2)
  have hboundSummable : Summable bound := by
    unfold bound
    simpa [bound] using
      shifted_zeta_two_hasSum.summable.mul_left (1 + Real.pi ^ 2 / 6)
  have hpoint (k : ℕ) : Tendsto (fun N => f N k) atTop (nhds 0) := by
    have hnorm : Tendsto (fun N => ‖leshF N k‖) atTop (nhds 0) := by
      refine squeeze_zero (fun _ => norm_nonneg _)
        (fun N => leshF_norm_le N k) ?_
      convert (leshCoeff_tendsto_zero.const_mul
          ((1 + Real.pi ^ 2 / 6) / (k + 1 : ℝ) ^ 2)) using 1
      · funext N
        ring
      · ring
    have hterm : Tendsto (fun N => leshF N k) atTop (nhds 0) :=
      inverseBinomial_tendsto_zero_of_norm24 hnorm
    apply hterm.congr'
    filter_upwards [eventually_gt_atTop k] with N hN
    simp [f, hN]
  have hbounded : Filter.Eventually (fun n => ∀ k, ‖f n k‖ ≤ bound k) atTop := by
    exact Eventually.of_forall fun n k => by
      by_cases hk : k < n
      · simp only [f, hk, if_true]
        calc
          ‖leshF n k‖ ≤
              (1 + Real.pi ^ 2 / 6) * leshCoeff n /
                (k + 1 : ℝ) ^ 2 := leshF_norm_le n k
          _ ≤ bound k := by
            unfold bound
            have hc : leshCoeff n ≤ 1 := by
              exact (leshCoeff_le n).trans (by
                rw [div_le_one (by positivity)]
                nlinarith [show (1 : ℝ) ≤ n + 1 by norm_num])
            calc
              (1 + Real.pi ^ 2 / 6) * leshCoeff n /
                    (k + 1 : ℝ) ^ 2 ≤
                  (1 + Real.pi ^ 2 / 6) * 1 /
                    (k + 1 : ℝ) ^ 2 :=
                div_le_div_of_nonneg_right
                  (mul_le_mul_of_nonneg_left hc (by positivity)) (sq_nonneg _)
              _ = (1 + Real.pi ^ 2 / 6) *
                    (1 / (k + 1 : ℝ) ^ 2) := by ring
      · simp [f, hk, bound]
        positivity
  have ht := tendsto_tsum_of_dominated_convergence
    hboundSummable hpoint hbounded
  have hfinite (N : ℕ) :
      (∑' k, f N k) = ∑ k ∈ range N, leshF N k := by
    rw [tsum_eq_sum (s := range N) (by
      intro k hk
      simp only [mem_range, not_lt] at hk
      simp [f, hk])]
    apply sum_congr rfl
    intro k hk
    simp [f, mem_range.mp hk]
  convert ht using 1
  · funext N
    exact (hfinite N).symm
  · simp

private theorem leshG_boundary :
    Tendsto (fun N => ∑ n ∈ range N, leshG n N) atTop (nhds 0) := by
  let f : ℕ → ℕ → ℝ := fun N n => if n < N then leshG n N else 0
  let bound : ℕ → ℝ := fun n =>
    ((4 + 3 * (Real.pi ^ 2 / 6)) / 2) * leshCoeff n
  have hboundSummable : Summable bound := by
    unfold bound
    exact leshCoeff_summable.mul_left ((4 + 3 * (Real.pi ^ 2 / 6)) / 2)
  have hpoint (n : ℕ) : Tendsto (fun N => f N n) atTop (nhds 0) := by
    have hnorm : Tendsto (fun N => ‖leshG n N‖) atTop (nhds 0) := by
      refine squeeze_zero (fun _ => norm_nonneg _) (fun N => leshG_norm_le n N) ?_
      have hzero := (tendsto_one_div_add_atTop_nhds_zero_nat
        (𝕜 := ℝ)).pow 2
      convert hzero.const_mul
          (((4 + 3 * (Real.pi ^ 2 / 6)) * leshCoeff n) / 2) using 1
      · funext N
        field_simp
      · ring
    have hterm : Tendsto (fun N => leshG n N) atTop (nhds 0) :=
      inverseBinomial_tendsto_zero_of_norm24 hnorm
    apply hterm.congr'
    filter_upwards [eventually_gt_atTop n] with N hN
    simp [f, hN]
  have hbounded : Filter.Eventually (fun N => ∀ n, ‖f N n‖ ≤ bound n) atTop := by
    exact Eventually.of_forall fun N n => by
      by_cases hn : n < N
      · simp only [f, hn, if_true]
        calc
          ‖leshG n N‖ ≤
              (4 + 3 * (Real.pi ^ 2 / 6)) * leshCoeff n /
                (2 * (N + 1 : ℝ) ^ 2) := leshG_norm_le n N
          _ ≤ bound n := by
            unfold bound
            have hsq : (1 : ℝ) ≤ (N + 1 : ℝ) ^ 2 := by
              nlinarith [show (1 : ℝ) ≤ N + 1 by norm_num]
            have hc : 0 ≤ (4 + 3 * (Real.pi ^ 2 / 6)) * leshCoeff n := by
              exact mul_nonneg (by positivity) (leshCoeff_nonneg n)
            calc
              (4 + 3 * (Real.pi ^ 2 / 6)) * leshCoeff n /
                    (2 * (N + 1 : ℝ) ^ 2) ≤
                  (4 + 3 * (Real.pi ^ 2 / 6)) * leshCoeff n / 2 :=
                div_le_div_of_nonneg_left hc (by norm_num)
                  (by nlinarith : (2 : ℝ) ≤ 2 * (N + 1 : ℝ) ^ 2)
              _ = (4 + 3 * (Real.pi ^ 2 / 6)) / 2 * leshCoeff n := by ring
      · simp only [f, hn, if_false, norm_zero]
        unfold bound
        exact mul_nonneg (by positivity) (leshCoeff_nonneg n)
  have ht := tendsto_tsum_of_dominated_convergence
    hboundSummable hpoint hbounded
  have hfinite (N : ℕ) :
      (∑' n, f N n) = ∑ n ∈ range N, leshG n N := by
    rw [tsum_eq_sum (s := range N) (by
      intro n hn
      simp only [mem_range, not_lt] at hn
      simp [f, hn])]
    apply sum_congr rfl
    intro n hn
    simp [f, mem_range.mp hn]
  convert ht using 1
  · funext N
    exact (hfinite N).symm
  · simp

private theorem bbbF_boundary :
    Tendsto (fun N => ∑ k ∈ range N, bbbF N k) atTop (nhds 0) := by
  let f : ℕ → ℕ → ℝ := fun N k => if k < N then bbbF N k else 0
  let bound : ℕ → ℝ := fun k => Real.pi ^ 2 / (k + 1 : ℝ) ^ 2
  have hboundSummable : Summable bound := by
    unfold bound
    simpa [bound] using
      shifted_zeta_two_hasSum.summable.mul_left (Real.pi ^ 2)
  have hpoint (k : ℕ) : Tendsto (fun N => f N k) atTop (nhds 0) := by
    have hnorm : Tendsto (fun N => ‖bbbF N k‖) atTop (nhds 0) := by
      refine squeeze_zero (fun _ => norm_nonneg _) (fun N => bbbF_norm_le N k) ?_
      convert bbbCoeff_tendsto_zero.const_mul
          (Real.pi ^ 2 / (k + 1 : ℝ) ^ 2) using 1
      · funext N
        ring
      · ring
    have hterm : Tendsto (fun N => bbbF N k) atTop (nhds 0) :=
      inverseBinomial_tendsto_zero_of_norm24 hnorm
    apply hterm.congr'
    filter_upwards [eventually_gt_atTop k] with N hN
    simp [f, hN]
  have hbounded : Filter.Eventually (fun N => ∀ k, ‖f N k‖ ≤ bound k) atTop := by
    exact Eventually.of_forall fun N k => by
      by_cases hk : k < N
      · simp only [f, hk, if_true]
        calc
          ‖bbbF N k‖ ≤ Real.pi ^ 2 * bbbCoeff N /
              (k + 1 : ℝ) ^ 2 := bbbF_norm_le N k
          _ ≤ bound k := by
            unfold bound
            exact div_le_div_of_nonneg_right
              (by simpa using
                mul_le_mul_of_nonneg_left (bbbCoeff_le_one N) (sq_nonneg Real.pi))
              (sq_nonneg _)
      · simp [f, hk, bound]
        positivity
  have ht := tendsto_tsum_of_dominated_convergence
    hboundSummable hpoint hbounded
  have hfinite (N : ℕ) :
      (∑' k, f N k) = ∑ k ∈ range N, bbbF N k := by
    rw [tsum_eq_sum (s := range N) (by
      intro k hk
      simp only [mem_range, not_lt] at hk
      simp [f, hk])]
    apply sum_congr rfl
    intro k hk
    simp [f, mem_range.mp hk]
  convert ht using 1
  · funext N
    exact (hfinite N).symm
  · simp

private theorem bbbG_boundary :
    Tendsto (fun N => ∑ n ∈ range N, bbbG n N) atTop (nhds 0) := by
  let f : ℕ → ℕ → ℝ := fun N n => if n < N then bbbG n N else 0
  let bound : ℕ → ℝ := fun n => 3 * Real.pi ^ 2 * bbbCoeff n
  have hboundSummable : Summable bound := by
    unfold bound
    exact bbbCoeff_summable.mul_left (3 * Real.pi ^ 2)
  have hpoint (n : ℕ) : Tendsto (fun N => f N n) atTop (nhds 0) := by
    have hnorm : Tendsto (fun N => ‖bbbG n N‖) atTop (nhds 0) := by
      refine squeeze_zero (fun _ => norm_nonneg _) (fun N => bbbG_norm_le n N) ?_
      have hzero := tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)
      convert hzero.const_mul (3 * Real.pi ^ 2 * bbbCoeff n) using 1
      · funext N
        ring
      · ring
    have hterm : Tendsto (fun N => bbbG n N) atTop (nhds 0) :=
      inverseBinomial_tendsto_zero_of_norm24 hnorm
    apply hterm.congr'
    filter_upwards [eventually_gt_atTop n] with N hN
    simp [f, hN]
  have hbounded : Filter.Eventually (fun N => ∀ n, ‖f N n‖ ≤ bound n) atTop := by
    exact Eventually.of_forall fun N n => by
      by_cases hn : n < N
      · simp only [f, hn, if_true]
        calc
          ‖bbbG n N‖ ≤ 3 * Real.pi ^ 2 * bbbCoeff n /
              (N + 1 : ℝ) := bbbG_norm_le n N
          _ ≤ bound n := by
            unfold bound
            have hc : 0 ≤ 3 * Real.pi ^ 2 * bbbCoeff n :=
              mul_nonneg (by positivity) (bbbCoeff_nonneg n)
            exact div_le_self hc (by norm_num)
      · simp only [f, hn, if_false, norm_zero]
        unfold bound
        exact mul_nonneg (by positivity) (bbbCoeff_nonneg n)
  have ht := tendsto_tsum_of_dominated_convergence
    hboundSummable hpoint hbounded
  have hfinite (N : ℕ) :
      (∑' n, f N n) = ∑ n ∈ range N, bbbG n N := by
    rw [tsum_eq_sum (s := range N) (by
      intro n hn
      simp only [mem_range, not_lt] at hn
      simp [f, hn])]
    apply sum_congr rfl
    intro n hn
    simp [f, mem_range.mp hn]
  convert ht using 1
  · funext N
    exact (hfinite N).symm
  · simp

theorem leshchinerWeightFourTerm24_hasSum :
    HasSum leshchinerWeightFourTerm24
      ((7 / 8) * (Real.pi ^ 4 / 90)) := by
  have hsource :
      HasSum (leshF 0) ((7 / 8) * (Real.pi ^ 4 / 90)) := by
    convert alternating_zeta_four_hasSum24 using 1
    funext k
    exact leshF_zero k
  have h : HasSum (fun n => leshG n 0)
      ((7 / 8) * (Real.pi ^ 4 / 90)) := hasSum_of_diagonal_wz leshF leshG
    ((7 / 8) * (Real.pi ^ 4 / 90)) lesh_wz hsource
    (summable_leshchinerWeightFourTerm24.congr fun n => (leshG_zero n).symm)
    leshF_boundary leshG_boundary
  exact HasSum.congr_fun h fun n => (leshG_zero n).symm

theorem bbbWeightFourTerm24_hasSum :
    HasSum bbbWeightFourTerm24 (Real.pi ^ 4 / 90) := by
  have hsource : HasSum (bbbF 0) (Real.pi ^ 4 / 90) := by
    convert shifted_zeta_four_hasSum24 using 1
    funext k
    exact bbbF_zero k
  have h : HasSum (fun n => bbbG n 0) (Real.pi ^ 4 / 90) :=
    hasSum_of_diagonal_wz bbbF bbbG (Real.pi ^ 4 / 90)
    bbb_wz hsource
    (summable_bbbWeightFourTerm24.congr fun n => (bbbG_zero n).symm)
    bbbF_boundary bbbG_boundary
  exact HasSum.congr_fun h fun n => (bbbG_zero n).symm

theorem inverseCentralQuartic_hasSum :
    HasSum (fun n : ℕ => 1 / (((n : ℝ) + 1) ^ 4 *
        (Nat.choose (2 * (n + 1)) (n + 1) : ℝ)))
      (17 * Real.pi ^ 4 / 3240) := by
  have h := inverseCentralFourthCoefficient24_of_wz
    leshchinerWeightFourTerm24_hasSum bbbWeightFourTerm24_hasSum
  convert h using 1
  · funext n
    unfold inverseCentralFourthCoefficient24 inverseCentralCoefficient24
    have hc : (Nat.choose (2 * (n + 1)) (n + 1) : ℝ) ≠ 0 := by
      exact_mod_cast Nat.choose_ne_zero (by omega : n + 1 ≤ 2 * (n + 1))
    field_simp [hc]
  · ring

theorem inverseCentralHarmonicSquare_hasSum :
    HasSum (fun n : ℕ => harmonicSquare24 n /
        (((n : ℝ) + 1) ^ 2 *
          (Nat.choose (2 * (n + 1)) (n + 1) : ℝ)))
      (5 * Real.pi ^ 4 / 9720) := by
  have h : HasSum
      (fun n => harmonicSquare24 n * inverseCentralCoefficient24 n)
      ((1 / 9 : ℝ) *
        (3 * ((17 / 36 : ℝ) * (Real.pi ^ 4 / 90)) -
          Real.pi ^ 4 / 90)) := by
    have hU := inverseCentralFourthCoefficient24_of_wz
      leshchinerWeightFourTerm24_hasSum bbbWeightFourTerm24_hasSum
    convert ((hU.mul_left 3).sub
      bbbWeightFourTerm24_hasSum).mul_left (1 / 9 : ℝ) using 1
    · funext n
      unfold bbbWeightFourTerm24
      ring
  convert h using 1
  · funext n
    unfold inverseCentralCoefficient24
    ring
  · ring

end RamanujanChallenge
