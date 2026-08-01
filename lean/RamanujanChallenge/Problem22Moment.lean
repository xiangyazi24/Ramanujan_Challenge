import RamanujanChallenge.Problem22Concentration

/-!
# Problem 2.2: finite moment bounds for the Rivoal weights

This module turns the finite Stein identity from `Problem22Concentration` into
quantitative concentration at the cubic saddle

`k^3 = (n - k)^2`.

All estimates are finite polynomial inequalities; no asymptotic theorem about
the recurrence is used.
-/

noncomputable section

open Filter Topology Real
open scoped BigOperators

namespace RamanujanChallenge.P22

/-- The real total mass of the Rivoal weights. -/
def rivoalWeightMass22 (n : ℕ) : ℝ :=
  ((rivoalExplicitQ22 n : ℚ) : ℝ)

theorem rivoalWeightMass22_pos (n : ℕ) :
    0 < rivoalWeightMass22 n := by
  unfold rivoalWeightMass22
  exact_mod_cast rivoalExplicitQ22_pos n

theorem rivoalRealWeight22_sum_eq_mass (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1), rivoalRealWeight22 n k) =
      rivoalWeightMass22 n := by
  exact rivoalRealWeight22_sum n

theorem rivoalSaddleError22_succ (n k : ℕ) :
    rivoalSaddleError22 n (k + 1) - rivoalSaddleError22 n k =
      3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ) := by
  simp only [rivoalSaddleError22, Nat.cast_add, Nat.cast_one]
  ring

theorem rivoalBirth_sub_death22 (n k : ℕ) :
    rivoalBirth22 n k - rivoalDeath22 n k =
      2 * ((n : ℝ) - (k : ℝ)) ^ 2 -
        (2 * (n : ℝ) + (k : ℝ)) * rivoalSaddleError22 n k := by
  simp only [rivoalBirth22, rivoalDeath22, rivoalSaddleError22]
  ring

/-- The Stein identity with the constant test function: total birth and death
rates agree. -/
theorem rivoalBirthDeathBalance22 (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * rivoalDeath22 n k) =
      ∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * rivoalBirth22 n k := by
  have h := rivoalWeightStein22 n (fun _ => (1 : ℝ))
  simp only [sub_self, mul_zero, zero_add, mul_one] at h
  have h' :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalBirth22 n k) -
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalDeath22 n k) = 0 := by
    rw [← Finset.sum_sub_distrib]
    calc
      _ = ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (rivoalBirth22 n k - rivoalDeath22 n k) := by
              apply Finset.sum_congr rfl
              intro k hk
              ring
      _ = 0 := h
  exact (sub_eq_zero.mp h').symm

/-- On the finite support, the death rate controls the third moment. -/
theorem rivoalDeath_lower22
    {n k : ℕ} (_hk : k < n + 1) :
    2 * (n : ℝ) * (k : ℝ) ^ 3 ≤ rivoalDeath22 n k := by
  have hk0 : (0 : ℝ) ≤ k := by positivity
  simp only [rivoalDeath22]
  nlinarith [mul_nonneg hk0 (pow_nonneg hk0 3)]

/-- On the finite support and for `n ≥ 1`, the birth rate is at most
`5 n^3`. -/
theorem rivoalBirth_upper22
    {n k : ℕ} (hn : 1 ≤ n) (hk : k < n + 1) :
    rivoalBirth22 n k ≤ 5 * (n : ℝ) ^ 3 := by
  have hkn : k ≤ n := by omega
  have hn0 : (0 : ℝ) ≤ n := by positivity
  have hk0 : (0 : ℝ) ≤ k := by positivity
  have hcast : (k : ℝ) ≤ (n : ℝ) := by exact_mod_cast hkn
  have hdiff : 0 ≤ (n : ℝ) - (k : ℝ) := sub_nonneg.mpr hcast
  have hdiff_le : (n : ℝ) - (k : ℝ) ≤ (n : ℝ) := by linarith
  have hlin : 2 * (n : ℝ) + (k : ℝ) + 2 ≤ 5 * (n : ℝ) := by
    have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
    linarith
  simp only [rivoalBirth22]
  calc
    (2 * (n : ℝ) + (k : ℝ) + 2) * ((n : ℝ) - (k : ℝ)) ^ 2 ≤
        (5 * (n : ℝ)) * ((n : ℝ) - (k : ℝ)) ^ 2 := by
          gcongr
    _ ≤ (5 * (n : ℝ)) * (n : ℝ) ^ 2 := by
          gcongr
    _ = 5 * (n : ℝ) ^ 3 := by ring

/-- Cubic moment bound under the normalized Rivoal weights. -/
theorem rivoalWeightedCube22_le
    {n : ℕ} (hn : 1 ≤ n) :
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * (k : ℝ) ^ 3) ≤
      (5 / 2 : ℝ) * (n : ℝ) ^ 2 * rivoalWeightMass22 n := by
  have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
  have hdeath :
      2 * (n : ℝ) *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * (k : ℝ) ^ 3) ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalDeath22 n k := by
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum
    intro k hk
    calc
      2 * (n : ℝ) * (rivoalRealWeight22 n k * (k : ℝ) ^ 3) =
          rivoalRealWeight22 n k * (2 * (n : ℝ) * (k : ℝ) ^ 3) := by ring
      _ ≤ rivoalRealWeight22 n k * rivoalDeath22 n k :=
        mul_le_mul_of_nonneg_left
          (rivoalDeath_lower22 (Finset.mem_range.mp hk))
          (rivoalRealWeight22_nonneg n k)
  have hbirth :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalBirth22 n k) ≤
        5 * (n : ℝ) ^ 3 * rivoalWeightMass22 n := by
    calc
      _ ≤ ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * (5 * (n : ℝ) ^ 3) := by
            apply Finset.sum_le_sum
            intro k hk
            exact mul_le_mul_of_nonneg_left
              (rivoalBirth_upper22 hn (Finset.mem_range.mp hk))
              (rivoalRealWeight22_nonneg n k)
      _ = 5 * (n : ℝ) ^ 3 * rivoalWeightMass22 n := by
            rw [← Finset.sum_mul, rivoalRealWeight22_sum_eq_mass]
            ring
  have hcombined := hdeath.trans
    ((rivoalBirthDeathBalance22 n).le.trans hbirth)
  calc
    _ ≤ (5 * (n : ℝ) ^ 3 * rivoalWeightMass22 n) /
        (2 * (n : ℝ)) := by
          apply (le_div_iff₀ (show 0 < 2 * (n : ℝ) by positivity)).2
          convert hcombined using 1
          ring
    _ = (5 / 2 : ℝ) * (n : ℝ) ^ 2 * rivoalWeightMass22 n := by
          field_simp

/-- Elementary interpolation inequality used to turn a third moment into a
second moment. -/
theorem sq_le_sq_add_cube_div22
    {x a : ℝ} (hx : 0 ≤ x) (ha : 0 < a) :
    x ^ 2 ≤ a ^ 2 + x ^ 3 / a := by
  rcases le_total x a with hxa | hax
  · have hsquare : x ^ 2 ≤ a ^ 2 := by nlinarith
    exact hsquare.trans (le_add_of_nonneg_right (div_nonneg (pow_nonneg hx 3) ha.le))
  · have hcube : 0 ≤ x ^ 3 / a := div_nonneg (pow_nonneg hx 3) ha.le
    have hmain : x ^ 2 ≤ x ^ 3 / a := by
      apply (le_div_iff₀ ha).2
      nlinarith [mul_nonneg (sq_nonneg x) (sub_nonneg.mpr hax)]
    exact hmain.trans (le_add_of_nonneg_left (sq_nonneg a))

/-- Quadratic moment bound under the Rivoal weights. -/
theorem rivoalWeightedSq22_le
    {n : ℕ} (hn : 1 ≤ n) :
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * (k : ℝ) ^ 2) ≤
      (7 / 2 : ℝ) * (n : ℝ) * Real.sqrt (n : ℝ) *
        rivoalWeightMass22 n := by
  have hnreal : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn
  have hsqrt : 0 < Real.sqrt (n : ℝ) := Real.sqrt_pos.2 hnreal
  have hsqrt_sq : Real.sqrt (n : ℝ) ^ 2 = (n : ℝ) :=
    Real.sq_sqrt hnreal.le
  have hsqrt_one : (1 : ℝ) ≤ Real.sqrt (n : ℝ) := by
    rw [← Real.sqrt_one]
    exact Real.sqrt_le_sqrt (by exact_mod_cast hn)
  have hpoint : ∀ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k * (k : ℝ) ^ 2 ≤
        rivoalRealWeight22 n k *
          ((n : ℝ) + (k : ℝ) ^ 3 / Real.sqrt (n : ℝ)) := by
    intro k hk
    have hinterp := sq_le_sq_add_cube_div22
      (show (0 : ℝ) ≤ (k : ℝ) by positivity) hsqrt
    rw [hsqrt_sq] at hinterp
    exact mul_le_mul_of_nonneg_left hinterp
      (rivoalRealWeight22_nonneg n k)
  calc
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * (k : ℝ) ^ 2) ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            ((n : ℝ) + (k : ℝ) ^ 3 / Real.sqrt (n : ℝ)) := by
          exact Finset.sum_le_sum hpoint
    _ = (n : ℝ) * rivoalWeightMass22 n +
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * (k : ℝ) ^ 3) /
            Real.sqrt (n : ℝ) := by
          simp_rw [mul_add]
          rw [Finset.sum_add_distrib]
          rw [← Finset.sum_mul, rivoalRealWeight22_sum_eq_mass]
          rw [Finset.sum_div]
          apply congrArg₂ (· + ·)
          · ring
          · apply Finset.sum_congr rfl
            intro k hk
            ring
    _ ≤ (n : ℝ) * rivoalWeightMass22 n +
        ((5 / 2 : ℝ) * (n : ℝ) ^ 2 * rivoalWeightMass22 n) /
          Real.sqrt (n : ℝ) := by
          gcongr
          exact rivoalWeightedCube22_le hn
    _ = (n : ℝ) * rivoalWeightMass22 n +
        (5 / 2 : ℝ) * (n : ℝ) * Real.sqrt (n : ℝ) *
          rivoalWeightMass22 n := by
          congr 1
          field_simp [hsqrt.ne']
          rw [hsqrt_sq]
          ring
    _ ≤ (7 / 2 : ℝ) * (n : ℝ) * Real.sqrt (n : ℝ) *
        rivoalWeightMass22 n := by
          have hmass := (rivoalWeightMass22_pos n).le
          have hbase :
              (n : ℝ) * rivoalWeightMass22 n ≤
                (n : ℝ) * Real.sqrt (n : ℝ) * rivoalWeightMass22 n := by
            calc
              (n : ℝ) * rivoalWeightMass22 n =
                  ((n : ℝ) * rivoalWeightMass22 n) * 1 := by ring
              _ ≤ ((n : ℝ) * rivoalWeightMass22 n) * Real.sqrt (n : ℝ) := by
                    gcongr
              _ = _ := by ring
          linarith

/-- Exact second-moment identity at the cubic saddle. -/
theorem rivoalSaddleSecondMomentExact22 (n : ℕ) :
    (∑ k ∈ Finset.range (n + 1),
      rivoalRealWeight22 n k *
        ((2 * (n : ℝ) + (k : ℝ)) * rivoalSaddleError22 n k ^ 2)) =
      ∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          (2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
            rivoalBirth22 n k *
              (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ))) := by
  have h := rivoalWeightStein22 n (rivoalSaddleError22 n)
  rw [← sub_eq_zero, ← Finset.sum_sub_distrib]
  calc
    _ = - (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k *
          (rivoalBirth22 n k *
              (rivoalSaddleError22 n (k + 1) - rivoalSaddleError22 n k) +
            (rivoalBirth22 n k - rivoalDeath22 n k) *
              rivoalSaddleError22 n k)) := by
          rw [← Finset.sum_neg_distrib]
          apply Finset.sum_congr rfl
          intro k hk
          rw [rivoalSaddleError22_succ, rivoalBirth_sub_death22]
          ring
    _ = 0 := by rw [h, neg_zero]

/-- Young's inequality in the exact form used in the saddle estimate. -/
theorem two_mul_le_mul_sq_add_sq_div22
    {a b c : ℝ} (hc : 0 < c) :
    2 * b * a ≤ c * a ^ 2 + b ^ 2 / c := by
  field_simp [hc.ne']
  nlinarith [sq_nonneg (c * a - b)]

theorem rivoalSaddleIncrement_upper22
    {n k : ℕ} (hk : k < n + 1) :
    3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ) ≤
      4 * (k : ℝ) ^ 2 + 2 * (n : ℝ) := by
  have hkn : k ≤ n := by omega
  have hk_sq : (k : ℝ) ≤ (k : ℝ) ^ 2 := by
    rcases k with _ | k
    · norm_num
    · have : (1 : ℝ) ≤ (k + 1 : ℕ) := by exact_mod_cast Nat.succ_le_succ (Nat.zero_le k)
      nlinarith
  linarith

theorem rivoalBirthIncrement_upper22
    {n k : ℕ} (hn : 1 ≤ n) (hk : k < n + 1) :
    rivoalBirth22 n k *
        (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)) ≤
      20 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2 + 10 * (n : ℝ) ^ 4 := by
  have hbirth := rivoalBirth_upper22 hn hk
  have hinc := rivoalSaddleIncrement_upper22 hk
  have hbirth0 : 0 ≤ rivoalBirth22 n k := by
    simp only [rivoalBirth22]
    positivity
  have hinc0 :
      0 ≤ 3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ) := by positivity
  calc
    _ ≤ (5 * (n : ℝ) ^ 3) *
        (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)) := by
          exact mul_le_mul_of_nonneg_right hbirth hinc0
    _ ≤ (5 * (n : ℝ) ^ 3) *
        (4 * (k : ℝ) ^ 2 + 2 * (n : ℝ)) := by
          exact mul_le_mul_of_nonneg_left hinc (by positivity)
    _ = 20 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2 + 10 * (n : ℝ) ^ 4 := by
          ring

/-- Explicit quadratic concentration at the cubic saddle. -/
theorem rivoalSaddleSecondMoment22_le
    {n : ℕ} (hn : 1 ≤ n) :
    (∑ k ∈ Finset.range (n + 1),
        rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) ≤
      81 * (n : ℝ) ^ 3 * Real.sqrt (n : ℝ) *
        rivoalWeightMass22 n := by
  have hnreal : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn
  have hsqrt_one : (1 : ℝ) ≤ Real.sqrt (n : ℝ) := by
    rw [← Real.sqrt_one]
    exact Real.sqrt_le_sqrt (by exact_mod_cast hn)
  have hleft :
      2 * (n : ℝ) *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) ≤
        ∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            ((2 * (n : ℝ) + (k : ℝ)) *
              rivoalSaddleError22 n k ^ 2) := by
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum
    intro k hk
    have hw := rivoalRealWeight22_nonneg n k
    have hg : 0 ≤ rivoalSaddleError22 n k ^ 2 := sq_nonneg _
    calc
      2 * (n : ℝ) *
          (rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) =
          rivoalRealWeight22 n k *
            (2 * (n : ℝ) * rivoalSaddleError22 n k ^ 2) := by ring
      _ ≤ rivoalRealWeight22 n k *
          ((2 * (n : ℝ) + (k : ℝ)) *
            rivoalSaddleError22 n k ^ 2) := by
            apply mul_le_mul_of_nonneg_left _ hw
            apply mul_le_mul_of_nonneg_right _ hg
            exact le_add_of_nonneg_right (Nat.cast_nonneg k)
  have hright :
      (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k *
            (2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
              rivoalBirth22 n k *
                (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)))) ≤
        (n : ℝ) *
            (∑ k ∈ Finset.range (n + 1),
              rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) +
          20 * (n : ℝ) ^ 3 *
            (∑ k ∈ Finset.range (n + 1),
              rivoalRealWeight22 n k * (k : ℝ) ^ 2) +
          11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n := by
    calc
      _ ≤ ∑ k ∈ Finset.range (n + 1),
          (rivoalRealWeight22 n k *
              ((n : ℝ) * rivoalSaddleError22 n k ^ 2) +
            rivoalRealWeight22 n k *
              (20 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2) +
            rivoalRealWeight22 n k * (11 * (n : ℝ) ^ 4)) := by
              apply Finset.sum_le_sum
              intro k hk
              have hklt : k < n + 1 := Finset.mem_range.mp hk
              have hkn : k ≤ n := by omega
              have hcast : (k : ℝ) ≤ (n : ℝ) := by exact_mod_cast hkn
              have hdiff0 : 0 ≤ (n : ℝ) - (k : ℝ) := sub_nonneg.mpr hcast
              have hdiff_le : (n : ℝ) - (k : ℝ) ≤ (n : ℝ) := by linarith
              have hyoung := two_mul_le_mul_sq_add_sq_div22
                (a := rivoalSaddleError22 n k)
                (b := ((n : ℝ) - (k : ℝ)) ^ 2) hnreal
              have hfour :
                  (((n : ℝ) - (k : ℝ)) ^ 2) ^ 2 / (n : ℝ) ≤
                    (n : ℝ) ^ 3 := by
                apply (div_le_iff₀ hnreal).2
                have hp : ((n : ℝ) - (k : ℝ)) ^ 4 ≤ (n : ℝ) ^ 4 := by
                  gcongr
                nlinarith
              have hfirst :
                  2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k ≤
                    (n : ℝ) * rivoalSaddleError22 n k ^ 2 +
                      (n : ℝ) ^ 3 := by
                calc
                  _ ≤ (n : ℝ) * rivoalSaddleError22 n k ^ 2 +
                      (((n : ℝ) - (k : ℝ)) ^ 2) ^ 2 / (n : ℝ) := hyoung
                  _ ≤ _ := add_le_add_right hfour _
              have hsecond := rivoalBirthIncrement_upper22 hn
                (Finset.mem_range.mp hk)
              have hcore :
                  2 * ((n : ℝ) - (k : ℝ)) ^ 2 * rivoalSaddleError22 n k +
                      rivoalBirth22 n k *
                        (3 * (k : ℝ) ^ 2 + (k : ℝ) + 2 * (n : ℝ)) ≤
                    (n : ℝ) * rivoalSaddleError22 n k ^ 2 +
                      20 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2 +
                      11 * (n : ℝ) ^ 4 := by
                have hn4 : (n : ℝ) ^ 3 ≤ (n : ℝ) ^ 4 := by
                  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
                  nlinarith [mul_nonneg (pow_nonneg hnreal.le 3) (sub_nonneg.mpr hn1)]
                nlinarith
              calc
                _ ≤ rivoalRealWeight22 n k *
                    ((n : ℝ) * rivoalSaddleError22 n k ^ 2 +
                      20 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2 +
                      11 * (n : ℝ) ^ 4) :=
                    mul_le_mul_of_nonneg_left hcore
                      (rivoalRealWeight22_nonneg n k)
                _ = _ := by ring
      _ = _ := by
            rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
            have h1 :
                (∑ k ∈ Finset.range (n + 1),
                    rivoalRealWeight22 n k *
                      ((n : ℝ) * rivoalSaddleError22 n k ^ 2)) =
                  (n : ℝ) *
                    (∑ k ∈ Finset.range (n + 1),
                      rivoalRealWeight22 n k *
                        rivoalSaddleError22 n k ^ 2) := by
              rw [Finset.mul_sum]
              apply Finset.sum_congr rfl
              intro k hk
              ring
            have h2 :
                (∑ k ∈ Finset.range (n + 1),
                    rivoalRealWeight22 n k *
                      (20 * (n : ℝ) ^ 3 * (k : ℝ) ^ 2)) =
                  20 * (n : ℝ) ^ 3 *
                    (∑ k ∈ Finset.range (n + 1),
                      rivoalRealWeight22 n k * (k : ℝ) ^ 2) := by
              rw [Finset.mul_sum]
              apply Finset.sum_congr rfl
              intro k hk
              ring
            have h3 :
                (∑ k ∈ Finset.range (n + 1),
                    rivoalRealWeight22 n k * (11 * (n : ℝ) ^ 4)) =
                  11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n := by
              rw [← Finset.sum_mul, rivoalRealWeight22_sum_eq_mass]
              ring
            rw [h1, h2, h3]
  have hchain := hleft.trans
    ((rivoalSaddleSecondMomentExact22 n).le.trans hright)
  have hmove :
      (n : ℝ) *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) ≤
        20 * (n : ℝ) ^ 3 *
            (∑ k ∈ Finset.range (n + 1),
              rivoalRealWeight22 n k * (k : ℝ) ^ 2) +
          11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n := by
    nlinarith
  have hsq := rivoalWeightedSq22_le hn
  have hsqterm :
      20 * (n : ℝ) ^ 3 *
          (∑ k ∈ Finset.range (n + 1),
            rivoalRealWeight22 n k * (k : ℝ) ^ 2) ≤
        70 * (n : ℝ) ^ 4 * Real.sqrt (n : ℝ) *
          rivoalWeightMass22 n := by
    calc
      _ ≤ 20 * (n : ℝ) ^ 3 *
          ((7 / 2 : ℝ) * (n : ℝ) * Real.sqrt (n : ℝ) *
            rivoalWeightMass22 n) := by
              exact mul_le_mul_of_nonneg_left hsq (by positivity)
      _ = _ := by ring
  have hconst :
      11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n ≤
        11 * (n : ℝ) ^ 4 * Real.sqrt (n : ℝ) *
          rivoalWeightMass22 n := by
    have hmass := (rivoalWeightMass22_pos n).le
    calc
      11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n =
          (11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n) * 1 := by ring
      _ ≤ (11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n) *
          Real.sqrt (n : ℝ) := by gcongr
      _ = _ := by ring
  have htotal :
      20 * (n : ℝ) ^ 3 *
            (∑ k ∈ Finset.range (n + 1),
              rivoalRealWeight22 n k * (k : ℝ) ^ 2) +
          11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n ≤
        81 * (n : ℝ) ^ 4 * Real.sqrt (n : ℝ) *
          rivoalWeightMass22 n := by
    linarith
  apply le_of_mul_le_mul_left _ hnreal
  calc
    (n : ℝ) *
        (∑ k ∈ Finset.range (n + 1),
          rivoalRealWeight22 n k * rivoalSaddleError22 n k ^ 2) ≤
      20 * (n : ℝ) ^ 3 *
            (∑ k ∈ Finset.range (n + 1),
              rivoalRealWeight22 n k * (k : ℝ) ^ 2) +
          11 * (n : ℝ) ^ 4 * rivoalWeightMass22 n := hmove
    _ ≤ 81 * (n : ℝ) ^ 4 * Real.sqrt (n : ℝ) *
          rivoalWeightMass22 n := htotal
    _ = (n : ℝ) *
        (81 * (n : ℝ) ^ 3 * Real.sqrt (n : ℝ) *
          rivoalWeightMass22 n) := by ring

end RamanujanChallenge.P22

end
