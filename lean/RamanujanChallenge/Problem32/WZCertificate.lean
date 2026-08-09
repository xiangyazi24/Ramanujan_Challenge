/-
  Problem 3.2 — WZ Certificate for the Apéry recurrence.

  The Apéry numbers b_n = Σ_{k=0}^n C(n,k)²C(n+k,k)² satisfy
    (n+2)³ b_{n+2} = P(n+1) b_{n+1} - (n+1)³ b_n
  where P(m) = 34m³ + 51m² + 27m + 5.

  Proof: Zeilberger's WZ certificate in division-free form.
  The anti-difference is G(n,k+1) = 4(2n+3)(4n²+12n-2k²-k+9) · f(n+1,k),
  with G(n,0) = 0 by convention and G(n,n+3) = 0 since C(n+1,n+2) = 0.
  The WZ equation -(n+1)³f(n,k) + P(n+1)f(n+1,k) - (n+2)³f(n+2,k) = G(n,k+1) - G(n,k)
  telescopes, proving the recurrence.
-/
import RamanujanChallenge.Problem32.AperyDef
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Tactic

noncomputable section

open Finset

/-! ## The summand (Apéry term) -/

def aperyTerm (n k : ℕ) : ℤ :=
  (Nat.choose n k : ℤ) ^ 2 * (Nat.choose (n + k) k : ℤ) ^ 2

theorem aperyB_eq_sum (n : ℕ) :
    aperyB n = ∑ k ∈ Finset.range (n + 1), aperyTerm n k := by
  simp only [aperyB, aperyTerm]

/-! ## Key binomial ratio identities (over ℤ, no division)

  These express shifted binomial coefficients in multiplicative form.
-/

-- (n+1-k) · C(n+1, k) = (n+1) · C(n, k)
theorem choose_succ_mul (n k : ℕ) (_hk : k ≤ n) :
    (n + 1 - k) * Nat.choose (n + 1) k = (n + 1) * Nat.choose n k := by
  simpa [Nat.mul_comm] using (Nat.choose_mul_succ_eq n k).symm

/-- Moving the upper Apéry index by one produces the square of the two
linear factors in the hypergeometric quotient, with no division. -/
private theorem aperyTerm_succ (n k : ℕ) (hk : k ≤ n) :
    ((n + 1 - k : ℕ) : ℤ) ^ 2 * aperyTerm (n + 1) k =
      ((n + 1 + k : ℕ) : ℤ) ^ 2 * aperyTerm n k := by
  have hC : ((n + 1 - k : ℕ) : ℤ) * (Nat.choose (n + 1) k : ℤ) =
      ((n + 1 : ℕ) : ℤ) * (Nat.choose n k : ℤ) := by
    exact_mod_cast choose_succ_mul n k hk
  have hDnat :
      (n + 1) * Nat.choose (n + 1 + k) k =
        (n + 1 + k) * Nat.choose (n + k) k := by
    have h := (Nat.choose_mul_succ_eq (n + k) k).symm
    have hsub : n + k + 1 - k = n + 1 := by omega
    rw [hsub] at h
    simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm, Nat.mul_comm] using h
  have hD : ((n + 1 : ℕ) : ℤ) * (Nat.choose (n + 1 + k) k : ℤ) =
      ((n + 1 + k : ℕ) : ℤ) * (Nat.choose (n + k) k : ℤ) := by
    exact_mod_cast hDnat
  unfold aperyTerm
  calc
    ((n + 1 - k : ℕ) : ℤ) ^ 2 *
          ((Nat.choose (n + 1) k : ℤ) ^ 2 *
            (Nat.choose (n + 1 + k) k : ℤ) ^ 2) =
        (((n + 1 - k : ℕ) : ℤ) * (Nat.choose (n + 1) k : ℤ)) ^ 2 *
          (Nat.choose (n + 1 + k) k : ℤ) ^ 2 := by ring
    _ = (((n + 1 : ℕ) : ℤ) * (Nat.choose n k : ℤ)) ^ 2 *
          (Nat.choose (n + 1 + k) k : ℤ) ^ 2 := by rw [hC]
    _ = (Nat.choose n k : ℤ) ^ 2 *
          (((n + 1 : ℕ) : ℤ) * (Nat.choose (n + 1 + k) k : ℤ)) ^ 2 := by ring
    _ = (Nat.choose n k : ℤ) ^ 2 *
          (((n + 1 + k : ℕ) : ℤ) * (Nat.choose (n + k) k : ℤ)) ^ 2 := by rw [hD]
    _ = ((n + 1 + k : ℕ) : ℤ) ^ 2 *
          ((Nat.choose n k : ℤ) ^ 2 * (Nat.choose (n + k) k : ℤ) ^ 2) := by ring

/-- Moving the summation index down by one gives the other hypergeometric
quotient in a division-free form. -/
private theorem aperyTerm_k_step (n k : ℕ) (hk0 : 1 ≤ k) (hkn : k ≤ n) :
    (k : ℤ) ^ 4 * aperyTerm n k =
      ((n + 1 - k : ℕ) : ℤ) ^ 2 * ((n + k : ℕ) : ℤ) ^ 2 *
        aperyTerm n (k - 1) := by
  have hCnat :
      k * Nat.choose n k = (n + 1 - k) * Nat.choose n (k - 1) := by
    have h := Nat.choose_succ_right_eq n (k - 1)
    have hkpred : k - 1 + 1 = k := by omega
    have hnsub : n - (k - 1) = n + 1 - k := by omega
    rw [hkpred, hnsub] at h
    simpa [Nat.mul_comm] using h
  have hDnat :
      k * Nat.choose (n + k) k =
        (n + k) * Nat.choose (n + (k - 1)) (k - 1) := by
    have h := (Nat.add_one_mul_choose_eq (n + k - 1) (k - 1)).symm
    have htop : n + k - 1 + 1 = n + k := by omega
    have hidx : k - 1 + 1 = k := by omega
    have hprev : n + k - 1 = n + (k - 1) := by omega
    rw [htop, hidx, hprev] at h
    simpa [Nat.mul_comm] using h
  have hC : (k : ℤ) * (Nat.choose n k : ℤ) =
      ((n + 1 - k : ℕ) : ℤ) * (Nat.choose n (k - 1) : ℤ) := by
    exact_mod_cast hCnat
  have hD : (k : ℤ) * (Nat.choose (n + k) k : ℤ) =
      ((n + k : ℕ) : ℤ) * (Nat.choose (n + (k - 1)) (k - 1) : ℤ) := by
    exact_mod_cast hDnat
  unfold aperyTerm
  calc
    (k : ℤ) ^ 4 *
          ((Nat.choose n k : ℤ) ^ 2 * (Nat.choose (n + k) k : ℤ) ^ 2) =
        ((k : ℤ) * (Nat.choose n k : ℤ)) ^ 2 *
          ((k : ℤ) * (Nat.choose (n + k) k : ℤ)) ^ 2 := by ring
    _ = (((n + 1 - k : ℕ) : ℤ) * (Nat.choose n (k - 1) : ℤ)) ^ 2 *
          (((n + k : ℕ) : ℤ) *
            (Nat.choose (n + (k - 1)) (k - 1) : ℤ)) ^ 2 := by rw [hC, hD]
    _ = ((n + 1 - k : ℕ) : ℤ) ^ 2 * ((n + k : ℕ) : ℤ) ^ 2 *
          ((Nat.choose n (k - 1) : ℤ) ^ 2 *
            (Nat.choose (n + (k - 1)) (k - 1) : ℤ) ^ 2) := by ring

/-- The quotient needed for the backward WZ flux combines one move in each
index; the intermediate factor `n + k` cancels without division. -/
private theorem aperyTerm_cross (n k : ℕ) (hk0 : 1 ≤ k) (hkn : k ≤ n) :
    (k : ℤ) ^ 4 * aperyTerm n k =
      ((n + 1 - k : ℕ) : ℤ) ^ 2 * ((n + 2 - k : ℕ) : ℤ) ^ 2 *
        aperyTerm (n + 1) (k - 1) := by
  have hkprev : k - 1 ≤ n := by omega
  have hs := aperyTerm_succ n (k - 1) hkprev
  have hleft : n + 1 - (k - 1) = n + 2 - k := by omega
  have hright : n + 1 + (k - 1) = n + k := by omega
  rw [hleft, hright] at hs
  calc
    (k : ℤ) ^ 4 * aperyTerm n k =
        ((n + 1 - k : ℕ) : ℤ) ^ 2 * ((n + k : ℕ) : ℤ) ^ 2 *
          aperyTerm n (k - 1) := aperyTerm_k_step n k hk0 hkn
    _ = ((n + 1 - k : ℕ) : ℤ) ^ 2 *
          (((n + k : ℕ) : ℤ) ^ 2 * aperyTerm n (k - 1)) := by ring
    _ = ((n + 1 - k : ℕ) : ℤ) ^ 2 *
          (((n + 2 - k : ℕ) : ℤ) ^ 2 * aperyTerm (n + 1) (k - 1)) := by rw [hs]
    _ = ((n + 1 - k : ℕ) : ℤ) ^ 2 * ((n + 2 - k : ℕ) : ℤ) ^ 2 *
          aperyTerm (n + 1) (k - 1) := by ring

/-- Two consecutive central binomial coefficients, in the division-free form
needed at the upper WZ boundary. -/
private theorem choose_central_succ (m : ℕ) :
    (m + 1) * Nat.choose (2 * m + 2) (m + 1) =
      2 * (2 * m + 1) * Nat.choose (2 * m) m := by
  have h1 := Nat.add_one_mul_choose_eq (2 * m) m
  have hsymm := Nat.choose_symm_half m
  have h1' :
      (2 * m + 1) * Nat.choose (2 * m) m =
        (m + 1) * Nat.choose (2 * m + 1) m := by
    rw [← hsymm]
    simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm, Nat.mul_comm] using h1
  have h2 := Nat.add_one_mul_choose_eq (2 * m + 1) m
  have hcancel :
      (m + 1) * (2 * Nat.choose (2 * m + 1) m) =
        (m + 1) * Nat.choose (2 * m + 2) (m + 1) := by
    calc
      _ = (2 * m + 2) * Nat.choose (2 * m + 1) m := by ring
      _ = Nat.choose (2 * m + 2) (m + 1) * (m + 1) := by
        simpa [Nat.add_assoc] using h2
      _ = _ := by ring
  have hdiag :
      2 * Nat.choose (2 * m + 1) m = Nat.choose (2 * m + 2) (m + 1) :=
    Nat.mul_left_cancel (by omega) hcancel
  calc
    (m + 1) * Nat.choose (2 * m + 2) (m + 1) =
        (m + 1) * (2 * Nat.choose (2 * m + 1) m) := by rw [hdiag]
    _ = 2 * ((m + 1) * Nat.choose (2 * m + 1) m) := by ring
    _ = 2 * ((2 * m + 1) * Nat.choose (2 * m) m) := by rw [h1']
    _ = _ := by ring

/-- The near-diagonal summand at the first upper boundary. -/
private theorem aperyTerm_upper (n : ℕ) :
    aperyTerm (n + 2) (n + 1) =
      (2 * (n : ℤ) + 3) ^ 2 * aperyTerm (n + 1) (n + 1) := by
  have h := Nat.add_one_mul_choose_eq (2 * n + 2) (n + 1)
  have hsymm := Nat.choose_symm_half (n + 1)
  have hidx : 2 * (n + 1) + 1 = 2 * n + 3 := by omega
  rw [hidx] at hsymm
  have hsymm' :
      Nat.choose (2 * n + 3) (n + 2) = Nat.choose (2 * n + 3) (n + 1) := by
    exact hsymm
  have hidx' : 2 * n + 2 + 1 = 2 * n + 3 := by omega
  rw [hidx'] at h
  have hnat :
      (n + 2) * Nat.choose (2 * n + 3) (n + 1) =
        (2 * n + 3) * Nat.choose (2 * n + 2) (n + 1) := by
    calc
      _ = Nat.choose (2 * n + 3) (n + 2) * (n + 2) := by rw [hsymm']; ring
      _ = _ := by
        simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using h.symm
  have hint :
      ((n + 2 : ℕ) : ℤ) * (Nat.choose (2 * n + 3) (n + 1) : ℤ) =
        ((2 * n + 3 : ℕ) : ℤ) * (Nat.choose (2 * n + 2) (n + 1) : ℤ) := by
    exact_mod_cast hnat
  simp only [aperyTerm, Nat.choose_succ_self_right, Nat.choose_self, Nat.cast_one, one_pow,
    one_mul]
  have htop1 : n + 2 + (n + 1) = 2 * n + 3 := by omega
  have htop0 : n + 1 + (n + 1) = 2 * n + 2 := by omega
  rw [htop1, htop0]
  push_cast at hint ⊢
  calc
    _ = (((n : ℤ) + 2) * (Nat.choose (2 * n + 3) (n + 1) : ℤ)) ^ 2 := by ring
    _ = ((2 * (n : ℤ) + 3) * (Nat.choose (2 * n + 2) (n + 1) : ℤ)) ^ 2 := by rw [hint]
    _ = _ := by ring

/-- The last diagonal summand and the WZ flux have the same central-binomial
scale at the final upper boundary. -/
private theorem aperyTerm_top (n : ℕ) :
    ((n : ℤ) + 2) ^ 2 * aperyTerm (n + 2) (n + 2) =
      4 * (2 * (n : ℤ) + 3) ^ 2 * aperyTerm (n + 1) (n + 1) := by
  have hnat := choose_central_succ (n + 1)
  have hint :
      ((n + 2 : ℕ) : ℤ) * (Nat.choose (2 * n + 4) (n + 2) : ℤ) =
        (2 * (2 * n + 3) : ℕ) * (Nat.choose (2 * n + 2) (n + 1) : ℤ) := by
    exact_mod_cast hnat
  simp only [aperyTerm, Nat.choose_self, Nat.cast_one, one_pow, one_mul]
  have htop2 : n + 2 + (n + 2) = 2 * n + 4 := by omega
  have htop1 : n + 1 + (n + 1) = 2 * n + 2 := by omega
  rw [htop2, htop1]
  push_cast at hint ⊢
  calc
    _ = (((n : ℤ) + 2) * (Nat.choose (2 * n + 4) (n + 2) : ℤ)) ^ 2 := by ring
    _ = (2 * (2 * (n : ℤ) + 3) *
          (Nat.choose (2 * n + 2) (n + 1) : ℤ)) ^ 2 := by rw [hint]
    _ = _ := by ring

/-- The penultimate summand is a quarter of the diagonal summand after the
obvious `(n+1)^2` factor is removed. -/
private theorem aperyTerm_penultimate (n : ℕ) :
    4 * aperyTerm (n + 1) n =
      ((n : ℤ) + 1) ^ 2 * aperyTerm (n + 1) (n + 1) := by
  have h := Nat.add_one_mul_choose_eq (2 * n + 1) n
  have hcancel :
      (n + 1) * (2 * Nat.choose (2 * n + 1) n) =
        (n + 1) * Nat.choose (2 * n + 2) (n + 1) := by
    calc
      _ = (2 * n + 2) * Nat.choose (2 * n + 1) n := by ring
      _ = Nat.choose (2 * n + 2) (n + 1) * (n + 1) := by
        simpa [Nat.add_assoc] using h
      _ = _ := by ring
  have hnat :
      2 * Nat.choose (2 * n + 1) n = Nat.choose (2 * n + 2) (n + 1) :=
    Nat.mul_left_cancel (by omega) hcancel
  have hint :
      2 * (Nat.choose (2 * n + 1) n : ℤ) =
        (Nat.choose (2 * n + 2) (n + 1) : ℤ) := by exact_mod_cast hnat
  simp only [aperyTerm, Nat.choose_succ_self_right, Nat.choose_self, Nat.cast_one, one_pow,
    one_mul]
  have htop0 : n + 1 + n = 2 * n + 1 := by omega
  have htop1 : n + 1 + (n + 1) = 2 * n + 2 := by omega
  rw [htop0, htop1]
  push_cast at hint ⊢
  calc
    _ = ((2 * (Nat.choose (2 * n + 1) n : ℤ)) ^ 2) * ((n : ℤ) + 1) ^ 2 := by ring
    _ = (Nat.choose (2 * n + 2) (n + 1) : ℤ) ^ 2 * ((n : ℤ) + 1) ^ 2 := by rw [hint]
    _ = _ := by ring

-- k · C(n, k) = (n - k + 1) · C(n, k-1)  [absorption]
-- (n+k+1) · C(n+k, k) = (n+1) · C(n+k+1, k+1) / (k+1) ... complicated
-- We may need the identity in "multiplied" form directly.

/-! ## The WZ polynomial identity

  After clearing all denominators and dividing by f(n+1,k),
  the WZ equation reduces to a polynomial identity in ℤ[n,k].
  This is verified by `ring`.
-/

-- The polynomial identity (degree 7, 28 terms):
-- LHS_poly(n,k) = RHS_poly(n,k) where both are polynomials in n, k
-- This is the core of the WZ proof.
theorem wz_polynomial_identity (n k : ℤ) :
    -(n+1)^3 * (n+1-k)^2 * (n+2-k)^2
    + (2*n+3) * (17*n^2+51*n+39) * (n+1+k)^2 * (n+2-k)^2
    - (n+2)^3 * (n+2+k)^2 * (n+1+k)^2
    = 4*(2*n+3) * (4*n^2+12*n-2*k^2-k+9) * (n+1+k)^2 * (n+2-k)^2
    - 4*(2*n+3) * (4*n^2+12*n-2*k^2+3*k+8) * k^4 := by
  ring

/-- The Zeilberger anti-difference.  At `k+1` it is the certificate
polynomial times the middle summand `f(n+1,k)`; the value at zero is the
vanishing lower boundary. -/
private def wzFlux (n : ℕ) : ℕ → ℤ
  | 0 => 0
  | k + 1 =>
      4 * (2 * (n : ℤ) + 3) *
        (4 * (n : ℤ) ^ 2 + 12 * (n : ℤ) - 2 * (k : ℤ) ^ 2 - (k : ℤ) + 9) *
          aperyTerm (n + 1) k

/-- The pointwise WZ identity away from the two upper boundary terms. -/
private theorem wz_step_interior (n k : ℕ) (hk : k ≤ n) :
    -((n : ℤ) + 1) ^ 3 * aperyTerm n k +
        (2 * (n : ℤ) + 3) * (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) *
          aperyTerm (n + 1) k -
        ((n : ℤ) + 2) ^ 3 * aperyTerm (n + 2) k =
      wzFlux n (k + 1) - wzFlux n k := by
  have hs1 := aperyTerm_succ n k hk
  have hs2 := aperyTerm_succ (n + 1) k (by omega)
  push_cast [Nat.cast_sub (by omega : k ≤ n + 1)] at hs1
  push_cast [Nat.cast_sub (by omega : k ≤ n + 2)] at hs2
  have hmid :
      (((n : ℤ) + 1 - (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2) *
          aperyTerm (n + 1) k =
        ((n : ℤ) + 1 + (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2 *
          aperyTerm n k := by
    calc
      _ = ((n : ℤ) + 2 - (k : ℤ)) ^ 2 *
            (((n : ℤ) + 1 - (k : ℤ)) ^ 2 * aperyTerm (n + 1) k) := by ring
      _ = ((n : ℤ) + 2 - (k : ℤ)) ^ 2 *
            (((n : ℤ) + 1 + (k : ℤ)) ^ 2 * aperyTerm n k) := by rw [hs1]
      _ = _ := by ring
  have htop :
      (((n : ℤ) + 1 - (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2) *
          aperyTerm (n + 2) k =
        ((n : ℤ) + 2 + (k : ℤ)) ^ 2 * ((n : ℤ) + 1 + (k : ℤ)) ^ 2 *
          aperyTerm n k := by
    calc
      _ = ((n : ℤ) + 1 - (k : ℤ)) ^ 2 *
            (((n : ℤ) + 2 - (k : ℤ)) ^ 2 * aperyTerm (n + 2) k) := by ring
      _ = ((n : ℤ) + 1 - (k : ℤ)) ^ 2 *
            (((n : ℤ) + 2 + (k : ℤ)) ^ 2 * aperyTerm (n + 1) k) := by
              rw [hs2]
              ring
      _ = ((n : ℤ) + 2 + (k : ℤ)) ^ 2 *
            (((n : ℤ) + 1 - (k : ℤ)) ^ 2 * aperyTerm (n + 1) k) := by ring
      _ = ((n : ℤ) + 2 + (k : ℤ)) ^ 2 *
            (((n : ℤ) + 1 + (k : ℤ)) ^ 2 * aperyTerm n k) := by rw [hs1]
      _ = _ := by ring
  have hApos : 0 < (n : ℤ) + 1 - (k : ℤ) := by omega
  have hBpos : 0 < (n : ℤ) + 2 - (k : ℤ) := by omega
  have hDne :
      ((n : ℤ) + 1 - (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2 ≠ 0 := by
    positivity
  apply mul_left_cancel₀ hDne
  calc
    (((n : ℤ) + 1 - (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2) *
          (-((n : ℤ) + 1) ^ 3 * aperyTerm n k +
            (2 * (n : ℤ) + 3) * (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) *
              aperyTerm (n + 1) k -
            ((n : ℤ) + 2) ^ 3 * aperyTerm (n + 2) k) =
        (-((n : ℤ) + 1) ^ 3 *
              ((n : ℤ) + 1 - (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2 +
          (2 * (n : ℤ) + 3) * (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) *
              ((n : ℤ) + 1 + (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2 -
          ((n : ℤ) + 2) ^ 3 * ((n : ℤ) + 2 + (k : ℤ)) ^ 2 *
              ((n : ℤ) + 1 + (k : ℤ)) ^ 2) * aperyTerm n k := by
            linear_combination
              ((2 * (n : ℤ) + 3) *
                (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39)) * hmid -
              ((n : ℤ) + 2) ^ 3 * htop
    _ = (4 * (2 * (n : ℤ) + 3) *
              (4 * (n : ℤ) ^ 2 + 12 * (n : ℤ) - 2 * (k : ℤ) ^ 2 - (k : ℤ) + 9) *
                ((n : ℤ) + 1 + (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2 -
          4 * (2 * (n : ℤ) + 3) *
              (4 * (n : ℤ) ^ 2 + 12 * (n : ℤ) - 2 * (k : ℤ) ^ 2 + 3 * (k : ℤ) + 8) *
                (k : ℤ) ^ 4) * aperyTerm n k := by
          rw [wz_polynomial_identity]
    _ = (((n : ℤ) + 1 - (k : ℤ)) ^ 2 * ((n : ℤ) + 2 - (k : ℤ)) ^ 2) *
          (wzFlux n (k + 1) - wzFlux n k) := by
      cases k with
      | zero =>
          simp [wzFlux, aperyTerm]
          ring
      | succ j =>
          have hc := aperyTerm_cross n (j + 1) (by omega) (by omega)
          push_cast [Nat.cast_sub (by omega : j + 1 ≤ n + 1),
            Nat.cast_sub (by omega : j + 1 ≤ n + 2)] at hc
          have hjn : j ≤ n := by omega
          have hjn1 : j ≤ n + 1 := by omega
          push_cast [Nat.cast_sub hjn, Nat.cast_sub hjn1] at hc
          simp only [wzFlux]
          push_cast at hmid hc ⊢
          linear_combination
            -(4 * (2 * (n : ℤ) + 3) *
                (4 * (n : ℤ) ^ 2 + 12 * (n : ℤ) -
                  2 * ((j : ℤ) + 1) ^ 2 - ((j : ℤ) + 1) + 9)) * hmid -
              (4 * (2 * (n : ℤ) + 3) *
                (4 * (n : ℤ) ^ 2 + 12 * (n : ℤ) -
                  2 * ((j : ℤ) + 1) ^ 2 + 3 * ((j : ℤ) + 1) + 8)) * hc

/-- The first of the two upper boundary terms in the WZ telescope. -/
private theorem wz_step_upper (n : ℕ) :
    -((n : ℤ) + 1) ^ 3 * aperyTerm n (n + 1) +
        (2 * (n : ℤ) + 3) * (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) *
          aperyTerm (n + 1) (n + 1) -
        ((n : ℤ) + 2) ^ 3 * aperyTerm (n + 2) (n + 1) =
      wzFlux n (n + 2) - wzFlux n (n + 1) := by
  have hz : aperyTerm n (n + 1) = 0 := by
    unfold aperyTerm
    rw [Nat.choose_eq_zero_of_lt (by omega : n < n + 1)]
    simp
  have hu := aperyTerm_upper n
  have hp := aperyTerm_penultimate n
  rw [hz]
  have hn2 : n + 2 = (n + 1) + 1 := by omega
  rw [hn2]
  simp only [wzFlux]
  push_cast at hu hp ⊢
  linear_combination
    -((n : ℤ) + 2) ^ 3 * hu +
      ((2 * (n : ℤ) + 3) *
        (4 * (n : ℤ) ^ 2 + 12 * (n : ℤ) - 2 * (n : ℤ) ^ 2 - (n : ℤ) + 9)) * hp

/-- The final upper boundary term; the outgoing flux is zero because the
middle binomial coefficient is already beyond its diagonal. -/
private theorem wz_step_top (n : ℕ) :
    -((n : ℤ) + 1) ^ 3 * aperyTerm n (n + 2) +
        (2 * (n : ℤ) + 3) * (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) *
          aperyTerm (n + 1) (n + 2) -
        ((n : ℤ) + 2) ^ 3 * aperyTerm (n + 2) (n + 2) =
      wzFlux n (n + 3) - wzFlux n (n + 2) := by
  have hz0 : aperyTerm n (n + 2) = 0 := by
    unfold aperyTerm
    rw [Nat.choose_eq_zero_of_lt (by omega : n < n + 2)]
    simp
  have hz1 : aperyTerm (n + 1) (n + 2) = 0 := by
    unfold aperyTerm
    rw [Nat.choose_eq_zero_of_lt (by omega : n + 1 < n + 2)]
    simp
  have ht := aperyTerm_top n
  rw [hz0, hz1]
  have hflux3 : wzFlux n (n + 3) = 0 := by
    rw [show n + 3 = (n + 2) + 1 by omega]
    simp only [wzFlux]
    rw [hz1]
    ring
  have hflux2 : wzFlux n (n + 2) =
      4 * (2 * (n : ℤ) + 3) *
        (4 * (n : ℤ) ^ 2 + 12 * (n : ℤ) -
          2 * ((n + 1 : ℕ) : ℤ) ^ 2 - ((n + 1 : ℕ) : ℤ) + 9) *
          aperyTerm (n + 1) (n + 1) := by
    rw [show n + 2 = (n + 1) + 1 by omega]
    simp only [wzFlux]
  rw [hflux3, hflux2]
  push_cast at ht ⊢
  linear_combination -((n : ℤ) + 2) * ht

/-- The three pointwise cases together cover exactly the range used in the
finite WZ telescope. -/
private theorem wz_step (n k : ℕ) (hk : k < n + 3) :
    -((n : ℤ) + 1) ^ 3 * aperyTerm n k +
        (2 * (n : ℤ) + 3) * (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) *
          aperyTerm (n + 1) k -
        ((n : ℤ) + 2) ^ 3 * aperyTerm (n + 2) k =
      wzFlux n (k + 1) - wzFlux n k := by
  by_cases hkn : k ≤ n
  · exact wz_step_interior n k hkn
  · have hkcases : k = n + 1 ∨ k = n + 2 := by omega
    rcases hkcases with rfl | rfl
    · exact wz_step_upper n
    · exact wz_step_top n

/-- Extending the defining Apéry sum by one place adds only a zero binomial
term. -/
private theorem sum_aperyTerm_pad_one (n : ℕ) :
    ∑ k ∈ Finset.range (n + 2), aperyTerm n k = aperyB n := by
  rw [show n + 2 = (n + 1) + 1 by omega, Finset.sum_range_succ]
  have hz : aperyTerm n (n + 1) = 0 := by
    unfold aperyTerm
    rw [Nat.choose_eq_zero_of_lt (by omega : n < n + 1)]
    simp
  rw [hz, add_zero]
  exact (aperyB_eq_sum n).symm

/-- Extending by two places likewise leaves the Apéry sum unchanged. -/
private theorem sum_aperyTerm_pad_two (n : ℕ) :
    ∑ k ∈ Finset.range (n + 3), aperyTerm n k = aperyB n := by
  rw [show n + 3 = (n + 2) + 1 by omega, Finset.sum_range_succ]
  have hz : aperyTerm n (n + 2) = 0 := by
    unfold aperyTerm
    rw [Nat.choose_eq_zero_of_lt (by omega : n < n + 2)]
    simp
  rw [hz, add_zero]
  exact sum_aperyTerm_pad_one n

/-! ## The recurrence theorem

  Using the WZ certificate, we prove:
  (n+2)³ b_{n+2} = P(n+1) b_{n+1} - (n+1)³ b_n  for all n ≥ 0.

  This is equivalent to `aperyB_recurrence_int` with index shifted by 1.
-/

theorem aperyB_recurrence_shifted (n : ℕ) :
    ((n : ℤ) + 2) ^ 3 * aperyB (n + 2) =
      (34 * ((n : ℤ) + 1) ^ 3 + 51 * ((n : ℤ) + 1) ^ 2 + 27 * ((n : ℤ) + 1) + 5) * aperyB (n + 1) -
      ((n : ℤ) + 1) ^ 3 * aperyB n := by
  have hsum :
      (∑ k ∈ Finset.range (n + 3),
        (-((n : ℤ) + 1) ^ 3 * aperyTerm n k +
            (2 * (n : ℤ) + 3) *
                (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) * aperyTerm (n + 1) k -
            ((n : ℤ) + 2) ^ 3 * aperyTerm (n + 2) k)) =
        ∑ k ∈ Finset.range (n + 3), (wzFlux n (k + 1) - wzFlux n k) := by
    apply Finset.sum_congr rfl
    intro k hk
    exact wz_step n k (Finset.mem_range.mp hk)
  have htel :
      (∑ k ∈ Finset.range (n + 3), (wzFlux n (k + 1) - wzFlux n k)) = 0 := by
    rw [Finset.sum_range_sub]
    have hz : aperyTerm (n + 1) (n + 2) = 0 := by
      unfold aperyTerm
      rw [Nat.choose_eq_zero_of_lt (by omega : n + 1 < n + 2)]
      simp
    rw [show n + 3 = (n + 2) + 1 by omega]
    simp only [wzFlux]
    rw [hz]
    ring
  have hleft :
      (∑ k ∈ Finset.range (n + 3),
        (-((n : ℤ) + 1) ^ 3 * aperyTerm n k +
            (2 * (n : ℤ) + 3) *
                (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) * aperyTerm (n + 1) k -
            ((n : ℤ) + 2) ^ 3 * aperyTerm (n + 2) k)) =
        -((n : ℤ) + 1) ^ 3 * aperyB n +
          (2 * (n : ℤ) + 3) *
              (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) * aperyB (n + 1) -
          ((n : ℤ) + 2) ^ 3 * aperyB (n + 2) := by
    rw [Finset.sum_sub_distrib, Finset.sum_add_distrib]
    simp only [← Finset.mul_sum]
    rw [sum_aperyTerm_pad_two n, sum_aperyTerm_pad_one (n + 1), aperyB_eq_sum (n + 2)]
  have hclosed :
      -((n : ℤ) + 1) ^ 3 * aperyB n +
          (2 * (n : ℤ) + 3) *
              (17 * (n : ℤ) ^ 2 + 51 * (n : ℤ) + 39) * aperyB (n + 1) -
          ((n : ℤ) + 2) ^ 3 * aperyB (n + 2) = 0 := by
    rw [← hleft, hsum, htel]
  linear_combination -hclosed

end
