import RamanujanChallenge.Problem27BarnesShift

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real

noncomputable section

namespace RamanujanChallenge.P27

/-! ## Removable triple zeros on the contour-shift strips -/

/-- The closed vertical strip with real part in
`[m - 1/2, m + 1/2]`. -/
def ctClosedStrip27 (m : ℕ) : Set ℂ :=
  {t | (m : ℝ) - 1 / 2 ≤ t.re ∧ t.re ≤ (m : ℝ) + 1 / 2}

private theorem ctComplex_ne_zero_of_re_pos27
    {z : ℂ} (hz : 0 < z.re) : z ≠ 0 := by
  intro h
  rw [h] at hz
  norm_num at hz

/-- Every denominator factor `t+j`, `j ≥ 0`, is nonzero on a strip
centered at a positive integer. -/
private theorem ctPoleFactor_ne_zero_on_strip27
    {m j : ℕ} (hm : 1 ≤ m) {t : ℂ}
    (ht : t ∈ ctClosedStrip27 m) :
    t + (j : ℂ) ≠ 0 := by
  change (m : ℝ) - 1 / 2 ≤ t.re ∧
    t.re ≤ (m : ℝ) + 1 / 2 at ht
  apply ctComplex_ne_zero_of_re_pos27
  have hm' : (1 : ℝ) ≤ (m : ℝ) := by
    exact_mod_cast hm
  have hj' : (0 : ℝ) ≤ (j : ℝ) := Nat.cast_nonneg j
  have hpos : 0 < t.re + (j : ℝ) := by
    linarith [ht.1]
  simpa using hpos

/-- The complete pole product is nonzero throughout a positive strip. -/
theorem ctPoleProduct_ne_zero_on_strip27
    {m M : ℕ} (hm : 1 ≤ m) {t : ℂ}
    (ht : t ∈ ctClosedStrip27 m) :
    ctPoleProduct27 M t ≠ 0 := by
  unfold ctPoleProduct27
  rw [Finset.prod_ne_zero_iff]
  intro j hj
  exact ctPoleFactor_ne_zero_on_strip27 hm ht

/-! ### The deleted numerator factor -/

/-- The numerator product with the factor indexed by `m-1` deleted.
There is no division by `t-m`. -/
def ctNumeratorErase27 (n m : ℕ) (t : ℂ) : ℂ :=
  ∏ r ∈ (Finset.range n).erase (m - 1),
    (t - (((r + 1 : ℕ) : ℂ))) ^ 3

/-- Exact finite-product factorization.  It is a polynomial identity and
therefore holds also at `t=m`. -/
theorem ctNumerator_factor_erase27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) (t : ℂ) :
    ctNumerator27 n t =
      (t - (m : ℂ)) ^ 3 * ctNumeratorErase27 n m t := by
  have hidx : m - 1 ∈ Finset.range n :=
    Finset.mem_range.mpr (by omega)
  have hprod :=
    Finset.mul_prod_erase (Finset.range n)
      (fun r : ℕ => (t - (((r + 1 : ℕ) : ℂ))) ^ 3) hidx
  simpa only [ctNumerator27, ctNumeratorErase27,
    Nat.sub_add_cancel hm1] using hprod.symm

/-- `R_n` after deleting its triple zero at `m`. -/
def ctRQuotient27 (n m : ℕ) (t : ℂ) : ℂ :=
  ctNumeratorErase27 n m t /
    (((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t)

/-- Global factorization of `R_n`; no condition on `t`, and in particular
no condition `t ≠ m`. -/
theorem ctR_eq_cube_mul_quotient27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) (t : ℂ) :
    ctR27 n t =
      (t - (m : ℂ)) ^ 3 * ctRQuotient27 n m t := by
  unfold ctR27 ctRQuotient27
  rw [ctNumerator_factor_erase27 hm1 hmn t]
  simp only [div_eq_mul_inv]
  ring

private theorem ctNumeratorErase_differentiableAt27
    (n m : ℕ) (t : ℂ) :
    DifferentiableAt ℂ (ctNumeratorErase27 n m) t := by
  unfold ctNumeratorErase27
  exact DifferentiableAt.fun_finset_prod
    (u := (Finset.range n).erase (m - 1)) fun r _ =>
      (differentiableAt_id.sub_const _).fun_pow 3

/-- Stronger pointwise form used to compute the actual derivative at the
boundary of the closed strip. -/
theorem ctRQuotient_differentiableAt_on_strip27
    {n m : ℕ} (hm : 1 ≤ m) {t : ℂ}
    (ht : t ∈ ctClosedStrip27 m) :
    DifferentiableAt ℂ (ctRQuotient27 n m) t := by
  unfold ctRQuotient27
  apply DifferentiableAt.div
  · exact ctNumeratorErase_differentiableAt27 n m t
  · exact (differentiableAt_const _).mul
      (ctPoleProduct_differentiableAt27 (n + 1) t)
  · exact mul_ne_zero
      (pow_ne_zero _
        (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n)))
      (ctPoleProduct_ne_zero_on_strip27 hm ht)

/-- Holomorphy of the deleted-factor quotient on the whole closed strip. -/
theorem ctRQuotient_differentiableOn_strip27
    {n m : ℕ} (hm : 1 ≤ m) :
    DifferentiableOn ℂ (ctRQuotient27 n m) (ctClosedStrip27 m) := by
  intro t ht
  exact (ctRQuotient_differentiableAt_on_strip27
    (n := n) (m := m) hm ht).differentiableWithinAt

/-- The actual complex derivative of the deleted-factor quotient. -/
def ctRQuotientPrime27 (n m : ℕ) (t : ℂ) : ℂ :=
  deriv (ctRQuotient27 n m) t

/-! ### The analogous deleted factor for the certificate -/

/-- `S_n` after deleting the triple zero at `m`.  The hypothesis needed later
is `1 ≤ m < n`, since `S_n` contains `ctNumerator27 (n-1)`. -/
def ctSQuotient27 (n m : ℕ) (t : ℂ) : ℂ :=
  ctQhatC27 (n : ℂ) (t - (n : ℂ)) *
      ctNumeratorErase27 (n - 1) m t /
    (((n.factorial : ℂ) ^ 2) * ctPoleProduct27 (n + 1) t)

/-- Global triple-zero factorization for the cancelled certificate. -/
theorem ctS_eq_cube_mul_quotient27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) (t : ℂ) :
    ctS27 n t =
      (t - (m : ℂ)) ^ 3 * ctSQuotient27 n m t := by
  unfold ctS27 ctSQuotient27
  rw [ctNumerator_factor_erase27
    (n := n - 1) (m := m) hm1 (by omega) t]
  simp only [div_eq_mul_inv]
  ring

/-- Pointwise holomorphy of the certificate quotient on the strip. -/
theorem ctSQuotient_differentiableAt_on_strip27
    {n m : ℕ} (hm : 1 ≤ m) {t : ℂ}
    (ht : t ∈ ctClosedStrip27 m) :
    DifferentiableAt ℂ (ctSQuotient27 n m) t := by
  unfold ctSQuotient27
  apply DifferentiableAt.div
  · apply DifferentiableAt.mul
    · unfold ctQhatC27
      fun_prop
    · exact ctNumeratorErase_differentiableAt27 (n - 1) m t
  · exact (differentiableAt_const _).mul
      (ctPoleProduct_differentiableAt27 (n + 1) t)
  · exact mul_ne_zero
      (pow_ne_zero _
        (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n)))
      (ctPoleProduct_ne_zero_on_strip27 hm ht)

/-- Holomorphy of the certificate quotient on the closed strip. -/
theorem ctSQuotient_differentiableOn_strip27
    {n m : ℕ} (hm : 1 ≤ m) :
    DifferentiableOn ℂ (ctSQuotient27 n m) (ctClosedStrip27 m) := by
  intro t ht
  exact (ctSQuotient_differentiableAt_on_strip27
    (n := n) (m := m) hm ht).differentiableWithinAt

/-- The actual complex derivative of the certificate quotient. -/
def ctSQuotientPrime27 (n m : ℕ) (t : ℂ) : ℂ :=
  deriv (ctSQuotient27 n m) t

/-! ### One generic product-rule calculation -/

/-- Actual derivative of `(t-m)^3 Q(t)`. -/
private theorem ctDeriv_cube_mul27
    (Q : ℂ → ℂ) (m : ℕ) (t : ℂ)
    (hQ : DifferentiableAt ℂ Q t) :
    deriv (fun z : ℂ => (z - (m : ℂ)) ^ 3 * Q z) t =
      3 * (t - (m : ℂ)) ^ 2 * Q t +
        (t - (m : ℂ)) ^ 3 * deriv Q t := by
  have hcube : HasDerivAt
      (fun z : ℂ => (z - (m : ℂ)) ^ 3)
      (3 * (t - (m : ℂ)) ^ 2) t := by
    convert ((hasDerivAt_id t).sub_const (m : ℂ)).fun_pow 3 using 1 <;>
      norm_num <;> ring
  exact (hcube.mul hQ.hasDerivAt).deriv

/-- The common algebraic identity behind both removable `Phi` factors. -/
private theorem ctPhi_cube_mul27
    (F Q : ℂ → ℂ) (m : ℕ) (t : ℂ)
    (hF : ∀ z : ℂ, F z = (z - (m : ℂ)) ^ 3 * Q z)
    (hQ : DifferentiableAt ℂ Q t) :
    F t - deriv F t / 2 =
      (t - (m : ℂ)) ^ 2 *
        ((t - (m : ℂ)) * Q t -
          (3 * Q t + (t - (m : ℂ)) * deriv Q t) / 2) := by
  have hfun : F = fun z : ℂ => (z - (m : ℂ)) ^ 3 * Q z :=
    funext hF
  rw [hfun, ctDeriv_cube_mul27 Q m t hQ]
  ring

/-- Exact order-two factorization of `ctRPhi27` on the strip.  It uses the
actual derivative of the quotient and no logarithmic derivative. -/
theorem ctRPhi_factor_on_strip27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) {t : ℂ}
    (ht : t ∈ ctClosedStrip27 m) :
    ctRPhi27 n t =
      (t - (m : ℂ)) ^ 2 *
        ((t - (m : ℂ)) * ctRQuotient27 n m t -
          (3 * ctRQuotient27 n m t +
            (t - (m : ℂ)) * ctRQuotientPrime27 n m t) / 2) := by
  simpa only [ctRPhi27, ctRQuotientPrime27] using
    ctPhi_cube_mul27
      (ctR27 n) (ctRQuotient27 n m) m t
      (fun z => ctR_eq_cube_mul_quotient27 hm1 hmn z)
      (ctRQuotient_differentiableAt_on_strip27 hm1 ht)

/-- Exact order-two factorization of the certificate `Phi` factor. -/
theorem ctSPhi_factor_on_strip27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) {t : ℂ}
    (ht : t ∈ ctClosedStrip27 m) :
    ctSPhi27 n t =
      (t - (m : ℂ)) ^ 2 *
        ((t - (m : ℂ)) * ctSQuotient27 n m t -
          (3 * ctSQuotient27 n m t +
            (t - (m : ℂ)) * ctSQuotientPrime27 n m t) / 2) := by
  simpa only [ctSPhi27, ctSQuotientPrime27] using
    ctPhi_cube_mul27
      (ctS27 n) (ctSQuotient27 n m) m t
      (fun z => ctS_eq_cube_mul_quotient27 hm1 hmn z)
      (ctSQuotient_differentiableAt_on_strip27 hm1 ht)

#print axioms ctNumerator_factor_erase27
#print axioms ctR_eq_cube_mul_quotient27
#print axioms ctRQuotient_differentiableOn_strip27
#print axioms ctRPhi_factor_on_strip27
#print axioms ctS_eq_cube_mul_quotient27
#print axioms ctSQuotient_differentiableOn_strip27
#print axioms ctSPhi_factor_on_strip27

end RamanujanChallenge.P27
