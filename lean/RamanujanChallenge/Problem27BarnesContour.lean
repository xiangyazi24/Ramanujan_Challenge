import RamanujanChallenge.Problem27BarnesTripleZero
import RamanujanChallenge.Problem27BarnesMoments

/-!
# Problem 2.7: quantitative contour shifts

This file supplies the estimates omitted in the source paper's Barnes-integral
argument.  The functions used on a strip are genuine holomorphic extensions at
the crossed integer; no totalized quotient is asserted to be continuous there.
-/

open Filter Set MeasureTheory Topology
open scoped BigOperators Interval Real ComplexConjugate

noncomputable section

namespace RamanujanChallenge.P27

/-! ## Polynomial growth on the pole-free half-plane -/

def ctRightHalfPlane27 : Set ℂ := {z | 0 < z.re}

theorem isOpen_ctRightHalfPlane27 : IsOpen ctRightHalfPlane27 := by
  exact isOpen_lt continuous_const Complex.continuous_re

private theorem ctPoleFactor_ne_zero_on_right27
    {j : ℕ} {z : ℂ} (hz : z ∈ ctRightHalfPlane27) :
    z + (j : ℂ) ≠ 0 := by
  intro h
  have hre := congrArg Complex.re h
  change 0 < z.re at hz
  simp only [Complex.add_re, Complex.natCast_re, Complex.zero_re] at hre
  have hj : (0 : ℝ) ≤ (j : ℝ) := Nat.cast_nonneg j
  linarith

private theorem ctPoleProduct_ne_zero_on_right27
    {M : ℕ} {z : ℂ} (hz : z ∈ ctRightHalfPlane27) :
    ctPoleProduct27 M z ≠ 0 := by
  unfold ctPoleProduct27
  rw [Finset.prod_ne_zero_iff]
  intro j hj
  exact ctPoleFactor_ne_zero_on_right27 hz

private theorem ctClosedStrip_subset_right27
    {m : ℕ} (hm : 1 ≤ m) :
    ctClosedStrip27 m ⊆ ctRightHalfPlane27 := by
  intro z hz
  change (m : ℝ) - 1 / 2 ≤ z.re ∧
    z.re ≤ (m : ℝ) + 1 / 2 at hz
  change 0 < z.re
  have hm' : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
  linarith

/-- `f` and its first derivative have polynomial growth uniformly on
`re z ≥ 1/2`, and `f` is holomorphic on the open right half-plane. -/
def HasCTPolyGrowth27 (f : ℂ → ℂ) : Prop :=
  DifferentiableOn ℂ f ctRightHalfPlane27 ∧
    ∃ C : ℝ, ∃ d : ℕ, 0 ≤ C ∧
      ∀ z : ℂ, 1 / 2 ≤ z.re →
        ‖f z‖ + ‖deriv f z‖ ≤ C * (1 + ‖z‖) ^ d

private theorem ctCore_mem_right27 {z : ℂ} (hz : 1 / 2 ≤ z.re) :
    z ∈ ctRightHalfPlane27 := by
  change 0 < z.re
  linarith

private theorem hasCTPolyGrowth_const27 (c : ℂ) :
    HasCTPolyGrowth27 (fun _ : ℂ => c) := by
  refine ⟨differentiableOn_const c, ‖c‖, 0, norm_nonneg _, ?_⟩
  intro z hz
  simp

private theorem hasCTPolyGrowth_id27 :
    HasCTPolyGrowth27 (fun z : ℂ => z) := by
  refine ⟨differentiableOn_id, 2, 1, by norm_num, ?_⟩
  intro z hz
  rw [show deriv (fun w : ℂ => w) z = 1 by simp]
  norm_num
  nlinarith [norm_nonneg z]

private theorem HasCTPolyGrowth27.add
    {f g : ℂ → ℂ}
    (hf : HasCTPolyGrowth27 f) (hg : HasCTPolyGrowth27 g) :
    HasCTPolyGrowth27 (fun z => f z + g z) := by
  rcases hf with ⟨hfd, Cf, df, hCf, hf⟩
  rcases hg with ⟨hgd, Cg, dg, hCg, hg⟩
  refine ⟨hfd.add hgd, Cf + Cg, df + dg, add_nonneg hCf hCg, ?_⟩
  intro z hz
  have hzR := ctCore_mem_right27 hz
  have hfa : DifferentiableAt ℂ f z :=
    hfd.differentiableAt (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  have hga : DifferentiableAt ℂ g z :=
    hgd.differentiableAt (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  have hH : (1 : ℝ) ≤ 1 + ‖z‖ := by linarith [norm_nonneg z]
  have hdf : (1 + ‖z‖) ^ df ≤ (1 + ‖z‖) ^ (df + dg) :=
    pow_le_pow_right₀ hH (Nat.le_add_right df dg)
  have hdg : (1 + ‖z‖) ^ dg ≤ (1 + ‖z‖) ^ (df + dg) :=
    pow_le_pow_right₀ hH (Nat.le_add_left dg df)
  rw [show deriv (fun w => f w + g w) z = deriv f z + deriv g z by
    simpa using deriv_add hfa hga]
  calc
    ‖f z + g z‖ + ‖deriv f z + deriv g z‖ ≤
        (‖f z‖ + ‖g z‖) + (‖deriv f z‖ + ‖deriv g z‖) :=
      add_le_add (norm_add_le _ _) (norm_add_le _ _)
    _ = (‖f z‖ + ‖deriv f z‖) +
        (‖g z‖ + ‖deriv g z‖) := by ring
    _ ≤ Cf * (1 + ‖z‖) ^ df + Cg * (1 + ‖z‖) ^ dg :=
      add_le_add (hf z hz) (hg z hz)
    _ ≤ Cf * (1 + ‖z‖) ^ (df + dg) +
        Cg * (1 + ‖z‖) ^ (df + dg) :=
      add_le_add (mul_le_mul_of_nonneg_left hdf hCf)
        (mul_le_mul_of_nonneg_left hdg hCg)
    _ = (Cf + Cg) * (1 + ‖z‖) ^ (df + dg) := by ring

private theorem HasCTPolyGrowth27.mul
    {f g : ℂ → ℂ}
    (hf : HasCTPolyGrowth27 f) (hg : HasCTPolyGrowth27 g) :
    HasCTPolyGrowth27 (fun z => f z * g z) := by
  rcases hf with ⟨hfd, Cf, df, hCf, hf⟩
  rcases hg with ⟨hgd, Cg, dg, hCg, hg⟩
  refine ⟨hfd.mul hgd, Cf * Cg, df + dg, mul_nonneg hCf hCg, ?_⟩
  intro z hz
  have hzR := ctCore_mem_right27 hz
  have hfa : DifferentiableAt ℂ f z :=
    hfd.differentiableAt (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  have hga : DifferentiableAt ℂ g z :=
    hgd.differentiableAt (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  have hfz := hf z hz
  have hgz := hg z hz
  rw [show deriv (fun w => f w * g w) z =
      deriv f z * g z + f z * deriv g z by
    simpa using deriv_mul hfa hga]
  calc
    ‖f z * g z‖ + ‖deriv f z * g z + f z * deriv g z‖ ≤
        ‖f z‖ * ‖g z‖ +
          (‖deriv f z‖ * ‖g z‖ + ‖f z‖ * ‖deriv g z‖) := by
      rw [norm_mul]
      apply add_le_add (le_refl _)
      calc
        ‖deriv f z * g z + f z * deriv g z‖ ≤
            ‖deriv f z * g z‖ + ‖f z * deriv g z‖ := norm_add_le _ _
        _ = ‖deriv f z‖ * ‖g z‖ + ‖f z‖ * ‖deriv g z‖ := by
          rw [norm_mul, norm_mul]
    _ ≤ (‖f z‖ + ‖deriv f z‖) *
        (‖g z‖ + ‖deriv g z‖) := by
      nlinarith [norm_nonneg (deriv f z), norm_nonneg (deriv g z),
        norm_nonneg (f z), norm_nonneg (g z)]
    _ ≤ (Cf * (1 + ‖z‖) ^ df) *
        (Cg * (1 + ‖z‖) ^ dg) := by
      exact mul_le_mul hfz hgz (by positivity) (by positivity)
    _ = (Cf * Cg) * (1 + ‖z‖) ^ (df + dg) := by
      rw [pow_add]
      ring

private theorem HasCTPolyGrowth27.pow
    {f : ℂ → ℂ} (hf : HasCTPolyGrowth27 f) (n : ℕ) :
    HasCTPolyGrowth27 (fun z => (f z) ^ n) := by
  induction n with
  | zero => simpa using hasCTPolyGrowth_const27 1
  | succ n ih =>
      simpa [pow_succ] using ih.mul hf

private theorem hasCTPolyGrowth_affine_sub27 (c : ℂ) :
    HasCTPolyGrowth27 (fun z : ℂ => z - c) := by
  simpa [sub_eq_add_neg] using
    hasCTPolyGrowth_id27.add (hasCTPolyGrowth_const27 (-c))

private theorem hasCTPolyGrowth_finset_prod27
    {ι : Type*} [DecidableEq ι] (s : Finset ι) (f : ι → ℂ → ℂ)
    (hf : ∀ i ∈ s, HasCTPolyGrowth27 (f i)) :
    HasCTPolyGrowth27 (fun z => ∏ i ∈ s, f i z) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using hasCTPolyGrowth_const27 1
  | @insert a s ha ih =>
      have hfa := hf a (Finset.mem_insert_self a s)
      have hfs : ∀ i ∈ s, HasCTPolyGrowth27 (f i) := by
        intro i hi
        exact hf i (Finset.mem_insert_of_mem hi)
      simpa [Finset.prod_insert, ha] using hfa.mul (ih hfs)

private theorem norm_inv_poleFactor_le_two27
    (j : ℕ) {z : ℂ} (hz : 1 / 2 ≤ z.re) :
    ‖(z + (j : ℂ))⁻¹‖ ≤ 2 := by
  rw [norm_inv]
  have hre : 1 / 2 ≤ (z + (j : ℂ)).re := by
    have hj : (0 : ℝ) ≤ (j : ℝ) := Nat.cast_nonneg j
    simp only [Complex.add_re, Complex.natCast_re]
    linarith
  have hn : 1 / 2 ≤ ‖z + (j : ℂ)‖ :=
    hre.trans (Complex.re_le_norm _)
  have hnpos : 0 < ‖z + (j : ℂ)‖ := lt_of_lt_of_le (by norm_num) hn
  apply (inv_le_iff_one_le_mul₀' hnpos).2
  nlinarith

private theorem hasCTPolyGrowth_inv_poleFactor27 (j : ℕ) :
    HasCTPolyGrowth27 (fun z : ℂ => (z + (j : ℂ))⁻¹) := by
  have hdiff : DifferentiableOn ℂ
      (fun z : ℂ => (z + (j : ℂ))⁻¹) ctRightHalfPlane27 := by
    intro z hz
    exact (((hasDerivAt_id z).add_const (j : ℂ)).inv
      (ctPoleFactor_ne_zero_on_right27 hz)).differentiableAt.differentiableWithinAt
  refine ⟨hdiff, 6, 0, by norm_num, ?_⟩
  intro z hz
  have hzR := ctCore_mem_right27 hz
  have hne := ctPoleFactor_ne_zero_on_right27 (j := j) hzR
  have hderiv : HasDerivAt (fun w : ℂ => (w + (j : ℂ))⁻¹)
      (-1 / (z + (j : ℂ)) ^ 2) z := by
    simpa using ((hasDerivAt_id z).add_const (j : ℂ)).inv hne
  have hinv := norm_inv_poleFactor_le_two27 j hz
  have hbase : 0 < ‖z + (j : ℂ)‖ := norm_pos_iff.mpr hne
  have hprime : ‖deriv (fun w : ℂ => (w + (j : ℂ))⁻¹) z‖ ≤ 4 := by
    rw [hderiv.deriv, neg_div, norm_neg, norm_div, norm_one, norm_pow]
    rw [div_le_iff₀ (sq_pos_of_pos hbase)]
    have hn : (1 / 2 : ℝ) ≤ ‖z + (j : ℂ)‖ := by
      have hre : 1 / 2 ≤ (z + (j : ℂ)).re := by
        have hj : (0 : ℝ) ≤ (j : ℝ) := Nat.cast_nonneg j
        simp only [Complex.add_re, Complex.natCast_re]
        linarith
      exact hre.trans (Complex.re_le_norm _)
    nlinarith
  norm_num
  rw [norm_inv] at hinv
  linarith

private theorem hasCTPolyGrowth_poleInvProduct27 (N : ℕ) :
    HasCTPolyGrowth27
      (fun z : ℂ => ∏ j ∈ Finset.range N, (z + (j : ℂ))⁻¹) := by
  apply hasCTPolyGrowth_finset_prod27
  intro j hj
  exact hasCTPolyGrowth_inv_poleFactor27 j

private theorem hasCTPolyGrowth_numeratorErase27 (n m : ℕ) :
    HasCTPolyGrowth27 (ctNumeratorErase27 n m) := by
  classical
  unfold ctNumeratorErase27
  apply hasCTPolyGrowth_finset_prod27
  intro r hr
  exact (hasCTPolyGrowth_affine_sub27 (((r + 1 : ℕ) : ℂ))).pow 3

private theorem hasCTPolyGrowth_numerator27 (n : ℕ) :
    HasCTPolyGrowth27 (ctNumerator27 n) := by
  classical
  unfold ctNumerator27
  apply hasCTPolyGrowth_finset_prod27
  intro r hr
  exact (hasCTPolyGrowth_affine_sub27 (((r + 1 : ℕ) : ℂ))).pow 3

private theorem ctRQuotient_alt27 (n m : ℕ) :
    ctRQuotient27 n m = fun z : ℂ =>
      (((n.factorial : ℂ) ^ 2)⁻¹ * ctNumeratorErase27 n m z) *
        ∏ j ∈ Finset.range (n + 1), (z + (j : ℂ))⁻¹ := by
  funext z
  unfold ctRQuotient27 ctPoleProduct27
  rw [div_eq_mul_inv, mul_inv, Finset.prod_inv_distrib]
  ring

theorem hasCTPolyGrowth_RQuotient27 (n m : ℕ) :
    HasCTPolyGrowth27 (ctRQuotient27 n m) := by
  rw [ctRQuotient_alt27]
  exact (((hasCTPolyGrowth_const27 (((n.factorial : ℂ) ^ 2)⁻¹)).mul
    (hasCTPolyGrowth_numeratorErase27 n m)).mul
      (hasCTPolyGrowth_poleInvProduct27 (n + 1)))

private theorem ctR_alt27 (n : ℕ) :
    ctR27 n = fun z : ℂ =>
      (((n.factorial : ℂ) ^ 2)⁻¹ * ctNumerator27 n z) *
        ∏ j ∈ Finset.range (n + 1), (z + (j : ℂ))⁻¹ := by
  funext z
  unfold ctR27 ctPoleProduct27
  rw [div_eq_mul_inv, mul_inv, Finset.prod_inv_distrib]
  ring

theorem hasCTPolyGrowth_R27 (n : ℕ) :
    HasCTPolyGrowth27 (ctR27 n) := by
  rw [ctR_alt27]
  exact (((hasCTPolyGrowth_const27 (((n.factorial : ℂ) ^ 2)⁻¹)).mul
    (hasCTPolyGrowth_numerator27 n)).mul
      (hasCTPolyGrowth_poleInvProduct27 (n + 1)))

private theorem hasCTPolyGrowth_QhatShift27 (n : ℕ) :
    HasCTPolyGrowth27 (fun z : ℂ =>
      ctQhatC27 (n : ℂ) (z - (n : ℂ))) := by
  let x : ℂ → ℂ := fun z => z - (n : ℂ)
  have hx : HasCTPolyGrowth27 x := hasCTPolyGrowth_affine_sub27 (n : ℂ)
  have h0 := hasCTPolyGrowth_const27 (ctQ0C27 (n : ℂ))
  have h1 := (hasCTPolyGrowth_const27 (ctQ1C27 (n : ℂ))).mul hx
  have h2 := (hasCTPolyGrowth_const27 (ctQ2C27 (n : ℂ))).mul (hx.pow 2)
  have h3 := (hasCTPolyGrowth_const27 (ctQ3C27 (n : ℂ))).mul (hx.pow 3)
  have h4 := (hasCTPolyGrowth_const27 (ctQ4C27 (n : ℂ))).mul (hx.pow 4)
  have h5 := (hasCTPolyGrowth_const27 (ctQ5C27 (n : ℂ))).mul (hx.pow 5)
  have h6 := (hasCTPolyGrowth_const27 (ctQ6C27 (n : ℂ))).mul (hx.pow 6)
  simpa only [ctQhatC27, x] using
    ((((((h0.add h1).add h2).add h3).add h4).add h5).add h6)

private theorem ctSQuotient_alt27 (n m : ℕ) :
    ctSQuotient27 n m = fun z : ℂ =>
      (((n.factorial : ℂ) ^ 2)⁻¹ *
        (ctQhatC27 (n : ℂ) (z - (n : ℂ)) *
          ctNumeratorErase27 (n - 1) m z)) *
        ∏ j ∈ Finset.range (n + 1), (z + (j : ℂ))⁻¹ := by
  funext z
  unfold ctSQuotient27 ctPoleProduct27
  rw [div_eq_mul_inv, mul_inv, Finset.prod_inv_distrib]
  ring

theorem hasCTPolyGrowth_SQuotient27 (n m : ℕ) :
    HasCTPolyGrowth27 (ctSQuotient27 n m) := by
  rw [ctSQuotient_alt27]
  have hnum :=
    (hasCTPolyGrowth_const27 (((n.factorial : ℂ) ^ 2)⁻¹)).mul
      ((hasCTPolyGrowth_QhatShift27 n).mul
        (hasCTPolyGrowth_numeratorErase27 (n - 1) m))
  exact hnum.mul (hasCTPolyGrowth_poleInvProduct27 (n + 1))

private theorem ctS_alt27 (n : ℕ) :
    ctS27 n = fun z : ℂ =>
      (((n.factorial : ℂ) ^ 2)⁻¹ *
        (ctQhatC27 (n : ℂ) (z - (n : ℂ)) *
          ctNumerator27 (n - 1) z)) *
        ∏ j ∈ Finset.range (n + 1), (z + (j : ℂ))⁻¹ := by
  funext z
  unfold ctS27 ctPoleProduct27
  rw [div_eq_mul_inv, mul_inv, Finset.prod_inv_distrib]
  ring

theorem hasCTPolyGrowth_S27 (n : ℕ) :
    HasCTPolyGrowth27 (ctS27 n) := by
  rw [ctS_alt27]
  have hnum :=
    (hasCTPolyGrowth_const27 (((n.factorial : ℂ) ^ 2)⁻¹)).mul
      ((hasCTPolyGrowth_QhatShift27 n).mul
        (hasCTPolyGrowth_numerator27 (n - 1)))
  exact hnum.mul (hasCTPolyGrowth_poleInvProduct27 (n + 1))

/-! ## Genuine removable extensions on a unit strip -/

theorem ctRQuotient_differentiableOn_right27 (n m : ℕ) :
    DifferentiableOn ℂ (ctRQuotient27 n m) ctRightHalfPlane27 :=
  (hasCTPolyGrowth_RQuotient27 n m).1

theorem ctRQuotientPrime_differentiableOn_right27 (n m : ℕ) :
    DifferentiableOn ℂ (ctRQuotientPrime27 n m) ctRightHalfPlane27 := by
  unfold ctRQuotientPrime27
  exact (ctRQuotient_differentiableOn_right27 n m).deriv
    isOpen_ctRightHalfPlane27

theorem ctSQuotient_differentiableOn_right27 (n m : ℕ) :
    DifferentiableOn ℂ (ctSQuotient27 n m) ctRightHalfPlane27 :=
  (hasCTPolyGrowth_SQuotient27 n m).1

theorem ctSQuotientPrime_differentiableOn_right27 (n m : ℕ) :
    DifferentiableOn ℂ (ctSQuotientPrime27 n m) ctRightHalfPlane27 := by
  unfold ctSQuotientPrime27
  exact (ctSQuotient_differentiableOn_right27 n m).deriv
    isOpen_ctRightHalfPlane27

def ctRKernelBracket27 (n m : ℕ) (z : ℂ) : ℂ :=
  (z - (m : ℂ)) * ctRQuotient27 n m z -
    (3 * ctRQuotient27 n m z +
      (z - (m : ℂ)) * ctRQuotientPrime27 n m z) / 2

def ctSKernelBracket27 (n m : ℕ) (z : ℂ) : ℂ :=
  (z - (m : ℂ)) * ctSQuotient27 n m z -
    (3 * ctSQuotient27 n m z +
      (z - (m : ℂ)) * ctSQuotientPrime27 n m z) / 2

def ctRKernelExt27 (n m : ℕ) (z : ℂ) : ℂ :=
  ctRKernelBracket27 n m z *
    ((Real.pi : ℂ) / sineSlope27 (m : ℤ) z) ^ 2

def ctSKernelExt27 (n m : ℕ) (z : ℂ) : ℂ :=
  ctSKernelBracket27 n m z *
    ((Real.pi : ℂ) / sineSlope27 (m : ℤ) z) ^ 2

def ctRKernelRaw27 (n : ℕ) (z : ℂ) : ℂ :=
  ctRPhi27 n z * zudilinBarnesSquaredSineKernel27 z

def ctSKernelRaw27 (n : ℕ) (z : ℂ) : ℂ :=
  ctSPhi27 n z * zudilinBarnesSquaredSineKernel27 z

private theorem ctRKernelBracket_differentiableAt27
    {n m : ℕ} (hm : 1 ≤ m) {z : ℂ}
    (hz : z ∈ ctClosedStrip27 m) :
    DifferentiableAt ℂ (ctRKernelBracket27 n m) z := by
  have hzR := ctClosedStrip_subset_right27 hm hz
  have hQ : DifferentiableAt ℂ (ctRQuotient27 n m) z :=
    (ctRQuotient_differentiableOn_right27 n m).differentiableAt
      (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  have hQ' : DifferentiableAt ℂ (ctRQuotientPrime27 n m) z :=
    (ctRQuotientPrime_differentiableOn_right27 n m).differentiableAt
      (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  unfold ctRKernelBracket27
  fun_prop

private theorem ctSKernelBracket_differentiableAt27
    {n m : ℕ} (hm : 1 ≤ m) {z : ℂ}
    (hz : z ∈ ctClosedStrip27 m) :
    DifferentiableAt ℂ (ctSKernelBracket27 n m) z := by
  have hzR := ctClosedStrip_subset_right27 hm hz
  have hQ : DifferentiableAt ℂ (ctSQuotient27 n m) z :=
    (ctSQuotient_differentiableOn_right27 n m).differentiableAt
      (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  have hQ' : DifferentiableAt ℂ (ctSQuotientPrime27 n m) z := by
    exact (ctSQuotientPrime_differentiableOn_right27 n m).differentiableAt
      (isOpen_ctRightHalfPlane27.mem_nhds hzR)
  unfold ctSKernelBracket27
  fun_prop

private theorem sineSlope_ne_zero_on_ctStrip27
    {m : ℕ} {z : ℂ} (hz : z ∈ ctClosedStrip27 m) :
    sineSlope27 (m : ℤ) z ≠ 0 := by
  apply sineSlope_ne_zero_on_strip27 (m : ℤ)
  simpa [ctClosedStrip27, stripLeft27, stripRight27] using hz

theorem ctRKernelExt_differentiableOn_strip27
    {n m : ℕ} (hm : 1 ≤ m) :
    DifferentiableOn ℂ (ctRKernelExt27 n m) (ctClosedStrip27 m) := by
  intro z hz
  apply DifferentiableAt.differentiableWithinAt
  unfold ctRKernelExt27
  apply DifferentiableAt.mul
  · exact ctRKernelBracket_differentiableAt27 hm hz
  · apply DifferentiableAt.fun_pow
    apply DifferentiableAt.div
    · exact differentiableAt_const _
    · exact sineSlope_differentiable27 (m : ℤ) z
    · exact sineSlope_ne_zero_on_ctStrip27 hz

theorem ctSKernelExt_differentiableOn_strip27
    {n m : ℕ} (hm : 1 ≤ m) :
    DifferentiableOn ℂ (ctSKernelExt27 n m) (ctClosedStrip27 m) := by
  intro z hz
  apply DifferentiableAt.differentiableWithinAt
  unfold ctSKernelExt27
  apply DifferentiableAt.mul
  · exact ctSKernelBracket_differentiableAt27 hm hz
  · apply DifferentiableAt.fun_pow
    apply DifferentiableAt.div
    · exact differentiableAt_const _
    · exact sineSlope_differentiable27 (m : ℤ) z
    · exact sineSlope_ne_zero_on_ctStrip27 hz

theorem ctRKernelRaw_eq_ext27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n)
    {z : ℂ} (hz : z ∈ ctClosedStrip27 m)
    (hzm : z ≠ (m : ℂ)) :
    ctRKernelRaw27 n z = ctRKernelExt27 n m z := by
  have hs := sineSlope_ne_zero_on_ctStrip27 hz
  have hsub : z - (m : ℂ) ≠ 0 := sub_ne_zero.mpr hzm
  rw [ctRKernelRaw27, ctRKernelExt27,
    ctRPhi_factor_on_strip27 hm1 hmn hz]
  unfold zudilinBarnesSquaredSineKernel27
  rw [show Complex.sin ((Real.pi : ℂ) * z) = sinePi27 z by rfl]
  rw [sinePi_eq_sub_mul_sineSlope27 (m : ℤ) z]
  rw [show (((m : ℤ) : ℂ)) = (m : ℂ) by norm_cast]
  unfold ctRKernelBracket27
  field_simp [hsub, hs]

theorem ctSKernelRaw_eq_ext27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n)
    {z : ℂ} (hz : z ∈ ctClosedStrip27 m)
    (hzm : z ≠ (m : ℂ)) :
    ctSKernelRaw27 n z = ctSKernelExt27 n m z := by
  have hs := sineSlope_ne_zero_on_ctStrip27 hz
  have hsub : z - (m : ℂ) ≠ 0 := sub_ne_zero.mpr hzm
  rw [ctSKernelRaw27, ctSKernelExt27,
    ctSPhi_factor_on_strip27 hm1 hmn hz]
  unfold zudilinBarnesSquaredSineKernel27
  rw [show Complex.sin ((Real.pi : ℂ) * z) = sinePi27 z by rfl]
  rw [sinePi_eq_sub_mul_sineSlope27 (m : ℤ) z]
  rw [show (((m : ℤ) : ℂ)) = (m : ℂ) by norm_cast]
  unfold ctSKernelBracket27
  field_simp [hsub, hs]

/-! ## Bounds for the raw kernels away from the removable point -/

private theorem ctPhi_poly_bound_of_growth27
    {F : ℂ → ℂ} (hF : HasCTPolyGrowth27 F) :
    ∃ C : ℝ, ∃ d : ℕ, 0 ≤ C ∧
      ∀ z : ℂ, 1 / 2 ≤ z.re →
        ‖F z - deriv F z / 2‖ ≤ C * (1 + ‖z‖) ^ d := by
  rcases hF.2 with ⟨C, d, hC, hbound⟩
  refine ⟨C, d, hC, ?_⟩
  intro z hz
  calc
    ‖F z - deriv F z / 2‖ ≤ ‖F z‖ + ‖deriv F z / 2‖ :=
      norm_sub_le _ _
    _ = ‖F z‖ + ‖deriv F z‖ / 2 := by
      rw [norm_div]
      norm_num
    _ ≤ ‖F z‖ + ‖deriv F z‖ := by
      linarith [norm_nonneg (deriv F z)]
    _ ≤ C * (1 + ‖z‖) ^ d := hbound z hz

theorem ctRPhi_poly_bound27 (n : ℕ) :
    ∃ C : ℝ, ∃ d : ℕ, 0 ≤ C ∧
      ∀ z : ℂ, 1 / 2 ≤ z.re →
        ‖ctRPhi27 n z‖ ≤ C * (1 + ‖z‖) ^ d := by
  simpa only [ctRPhi27] using
    ctPhi_poly_bound_of_growth27 (hasCTPolyGrowth_R27 n)

theorem ctSPhi_poly_bound27 (n : ℕ) :
    ∃ C : ℝ, ∃ d : ℕ, 0 ≤ C ∧
      ∀ z : ℂ, 1 / 2 ≤ z.re →
        ‖ctSPhi27 n z‖ ≤ C * (1 + ‖z‖) ^ d := by
  simpa only [ctSPhi27] using
    ctPhi_poly_bound_of_growth27 (hasCTPolyGrowth_S27 n)

theorem sinePi_norm_sq27 (z : ℂ) :
    ‖sinePi27 z‖ ^ 2 =
      Real.sin (Real.pi * z.re) ^ 2 +
        Real.sinh (Real.pi * z.im) ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  have harg : (Real.pi : ℂ) * z =
      ((Real.pi * z.re : ℝ) : ℂ) +
        ((Real.pi * z.im : ℝ) : ℂ) * Complex.I := by
    apply Complex.ext <;> simp
  unfold sinePi27
  rw [harg, Complex.sin_add_mul_I]
  simp only [← Complex.ofReal_sin, ← Complex.ofReal_cos,
    ← Complex.ofReal_cosh, ← Complex.ofReal_sinh,
    Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im,
    Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
  norm_num
  ring_nf
  rw [Real.cosh_sq]
  have htrig := Real.sin_sq_add_cos_sq (Real.pi * z.re)
  ring_nf at htrig ⊢
  nlinarith

theorem abs_sinh_le_norm_sinePi27 (z : ℂ) :
    |Real.sinh (Real.pi * z.im)| ≤ ‖sinePi27 z‖ := by
  apply (sq_le_sq₀ (abs_nonneg _) (norm_nonneg _)).mp
  rw [sq_abs, sinePi_norm_sq27]
  nlinarith [sq_nonneg (Real.sin (Real.pi * z.re))]

private theorem sinh_lower_exp_quarter27 {T : ℝ} (hT : 1 ≤ T) :
    Real.exp (Real.pi * T) / 4 ≤ Real.sinh (Real.pi * T) := by
  have ha : 1 ≤ Real.pi * T := by
    nlinarith [Real.two_le_pi]
  have hnonneg : 0 ≤ Real.pi * T := le_trans (by norm_num) ha
  have hneg : Real.exp (-(Real.pi * T)) ≤ 1 :=
    Real.exp_le_one_iff.mpr (by linarith)
  have hpos : 2 ≤ Real.exp (Real.pi * T) := by
    linarith [Real.add_one_le_exp (Real.pi * T)]
  rw [Real.sinh_eq]
  nlinarith

theorem norm_squaredSineKernel_le_exp27
    {z : ℂ} (hz : 1 ≤ |z.im|) :
    ‖zudilinBarnesSquaredSineKernel27 z‖ ≤
      16 * Real.pi ^ 2 *
        Real.exp (-2 * Real.pi * |z.im|) := by
  let T : ℝ := |z.im|
  have hT : 1 ≤ T := hz
  have hsinh := sinh_lower_exp_quarter27 hT
  have hsinhpos : 0 < Real.sinh (Real.pi * T) := by
    exact Real.sinh_pos_iff.mpr (mul_pos Real.pi_pos (lt_of_lt_of_le zero_lt_one hT))
  have hsin : Real.exp (Real.pi * T) / 4 ≤ ‖sinePi27 z‖ := by
    calc
      Real.exp (Real.pi * T) / 4 ≤ Real.sinh (Real.pi * T) := hsinh
      _ = |Real.sinh (Real.pi * z.im)| := by
        rw [Real.abs_sinh]
        simp only [T, abs_mul, abs_of_pos Real.pi_pos]
      _ ≤ ‖sinePi27 z‖ := abs_sinh_le_norm_sinePi27 z
  have hdenpos : 0 < ‖sinePi27 z‖ :=
    lt_of_lt_of_le (div_pos (Real.exp_pos _) (by norm_num)) hsin
  have hquot : Real.pi / ‖sinePi27 z‖ ≤
      Real.pi / (Real.exp (Real.pi * T) / 4) := by
    exact div_le_div_of_nonneg_left Real.pi_pos.le
      (div_pos (Real.exp_pos _) (by norm_num)) hsin
  unfold zudilinBarnesSquaredSineKernel27
  rw [norm_pow, norm_div, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos Real.pi_pos]
  calc
    (Real.pi / ‖Complex.sin ((Real.pi : ℂ) * z)‖) ^ 2 ≤
        (Real.pi / (Real.exp (Real.pi * T) / 4)) ^ 2 :=
      pow_le_pow_left₀ (div_nonneg Real.pi_pos.le hdenpos.le) hquot 2
    _ = 16 * Real.pi ^ 2 * Real.exp (-2 * Real.pi * T) := by
      rw [show -2 * Real.pi * T = -(Real.pi * T) + -(Real.pi * T) by ring,
        Real.exp_add, Real.exp_neg]
      field_simp [ne_of_gt (Real.exp_pos (Real.pi * T))]
      ring
    _ = 16 * Real.pi ^ 2 * Real.exp (-2 * Real.pi * |z.im|) := rfl

/-! ## The squared-sine kernel on all half-integer lines -/

theorem squaredSineKernel_add_one27 (z : ℂ) :
    zudilinBarnesSquaredSineKernel27 (z + 1) =
      zudilinBarnesSquaredSineKernel27 z := by
  unfold zudilinBarnesSquaredSineKernel27
  rw [show (Real.pi : ℂ) * (z + 1) =
      (Real.pi : ℂ) * z + (Real.pi : ℂ) by ring,
    Complex.sin_add_pi]
  ring

theorem squaredSineKernel_add_nat27 (z : ℂ) (n : ℕ) :
    zudilinBarnesSquaredSineKernel27 (z + (n : ℂ)) =
      zudilinBarnesSquaredSineKernel27 z := by
  induction n with
  | zero => simp
  | succ n ih =>
      calc
        zudilinBarnesSquaredSineKernel27 (z + ((n + 1 : ℕ) : ℂ)) =
            zudilinBarnesSquaredSineKernel27 ((z + (n : ℂ)) + 1) := by
              congr 1
              push_cast
              ring
        _ = zudilinBarnesSquaredSineKernel27 (z + (n : ℂ)) :=
          squaredSineKernel_add_one27 _
        _ = zudilinBarnesSquaredSineKernel27 z := ih

theorem squaredSineKernel_left_halfLine27 (m : ℕ) (y : ℝ) :
    zudilinBarnesSquaredSineKernel27
        (verticalPoint27 ((m : ℝ) - 1 / 2) y) =
      (Real.pi : ℂ) ^ 2 * sechSq27 y := by
  have hpoint : verticalPoint27 ((m : ℝ) - 1 / 2) y =
      zudilinBarnesLine27 y + (m : ℂ) := by
    unfold verticalPoint27 zudilinBarnesLine27
    push_cast
    ring
  rw [hpoint, squaredSineKernel_add_nat27,
    zudilinBarnesSquaredSineKernel_line27]
  unfold sechSq27
  ring

theorem squaredSineKernel_right_halfLine27 (m : ℕ) (y : ℝ) :
    zudilinBarnesSquaredSineKernel27
        (verticalPoint27 ((m : ℝ) + 1 / 2) y) =
      (Real.pi : ℂ) ^ 2 * sechSq27 y := by
  have hpoint : verticalPoint27 ((m : ℝ) + 1 / 2) y =
      zudilinBarnesLine27 y + ((m + 1 : ℕ) : ℂ) := by
    unfold verticalPoint27 zudilinBarnesLine27
    push_cast
    ring
  rw [hpoint, squaredSineKernel_add_nat27,
    zudilinBarnesSquaredSineKernel_line27]
  unfold sechSq27
  ring

/-! ## Integrability on the vertical boundaries -/

private theorem ctRPhi_differentiableOn_right27 (n : ℕ) :
    DifferentiableOn ℂ (ctRPhi27 n) ctRightHalfPlane27 := by
  unfold ctRPhi27
  have hR := (hasCTPolyGrowth_R27 n).1
  exact hR.sub ((hR.deriv isOpen_ctRightHalfPlane27).div_const 2)

private theorem ctSPhi_differentiableOn_right27 (n : ℕ) :
    DifferentiableOn ℂ (ctSPhi27 n) ctRightHalfPlane27 := by
  unfold ctSPhi27
  have hS := (hasCTPolyGrowth_S27 n).1
  exact hS.sub ((hS.deriv isOpen_ctRightHalfPlane27).div_const 2)

private theorem continuous_sechSq27 : Continuous sechSq27 := by
  unfold sechSq27
  apply Continuous.div continuous_const
    ((Complex.ofRealCLM.continuous.comp
      (Real.continuous_cosh.comp (continuous_const.mul continuous_id))).pow 2)
  intro y
  exact pow_ne_zero _
    (Complex.ofReal_ne_zero.mpr (Real.cosh_pos _).ne')

private theorem integrable_one_add_abs_pow_mul_norm_sechSq27 (d : ℕ) :
    Integrable (fun y : ℝ => (1 + |y|) ^ d * ‖sechSq27 y‖) := by
  have h0 : Integrable (fun y : ℝ => ‖sechSq27 y‖) := by
    simpa using (integrable_pow_sechSq27 0).norm
  have hd : Integrable (fun y : ℝ => |y| ^ d * ‖sechSq27 y‖) := by
    simpa [norm_smul, Real.norm_eq_abs, abs_pow] using
      (integrable_pow_sechSq27 d).norm
  let g : ℝ → ℝ := fun y =>
    (2 : ℝ) ^ (d - 1) * (‖sechSq27 y‖ + |y| ^ d * ‖sechSq27 y‖)
  have hg : Integrable g := by
    exact (h0.add hd).const_mul ((2 : ℝ) ^ (d - 1))
  apply hg.mono'
  · apply Continuous.aestronglyMeasurable
    exact ((continuous_const.add continuous_abs).pow d).mul
      continuous_sechSq27.norm
  · filter_upwards with y
    have hp := add_pow_le (by norm_num : (0 : ℝ) ≤ 1) (abs_nonneg y) d
    calc
      ‖(1 + |y|) ^ d * ‖sechSq27 y‖‖ =
          (1 + |y|) ^ d * ‖sechSq27 y‖ := by
        rw [Real.norm_eq_abs, abs_of_nonneg]
        positivity
      _ ≤ (2 : ℝ) ^ (d - 1) * (1 ^ d + |y| ^ d) * ‖sechSq27 y‖ :=
        mul_le_mul_of_nonneg_right hp (norm_nonneg _)
      _ = g y := by
        dsimp only [g]
        simp
        ring

private theorem norm_verticalPoint_le_abs_add27 (x y : ℝ) :
    ‖verticalPoint27 x y‖ ≤ |x| + |y| := by
  unfold verticalPoint27
  calc
    ‖(x : ℂ) + (y : ℂ) * Complex.I‖ ≤
        ‖(x : ℂ)‖ + ‖(y : ℂ) * Complex.I‖ := norm_add_le _ _
    _ = |x| + |y| := by simp

private theorem integrable_phi_mul_sechSq27
    {Phi : ℂ → ℂ}
    (hPhiDiff : DifferentiableOn ℂ Phi ctRightHalfPlane27)
    (hPhiBound : ∃ C : ℝ, ∃ d : ℕ, 0 ≤ C ∧
      ∀ z : ℂ, 1 / 2 ≤ z.re →
        ‖Phi z‖ ≤ C * (1 + ‖z‖) ^ d)
    {x : ℝ} (hx : 1 / 2 ≤ x) :
    Integrable (fun y : ℝ =>
      Phi (verticalPoint27 x y) * ((Real.pi : ℂ) ^ 2 * sechSq27 y)) := by
  rcases hPhiBound with ⟨C, d, hC, hbound⟩
  let K : ℝ := C * (1 + |x|) ^ d * Real.pi ^ 2
  have hK : 0 ≤ K := by
    dsimp only [K]
    positivity
  have hmajor : Integrable
      (fun y : ℝ => K * ((1 + |y|) ^ d * ‖sechSq27 y‖)) :=
    (integrable_one_add_abs_pow_mul_norm_sechSq27 d).const_mul K
  apply hmajor.mono'
  · apply Continuous.aestronglyMeasurable
    apply Continuous.mul
    · apply continuous_iff_continuousAt.mpr
      intro y
      have hp : verticalPoint27 x y ∈ ctRightHalfPlane27 := by
        change 0 < (verticalPoint27 x y).re
        simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
          Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
          Complex.I_im]
        norm_num
        linarith
      exact (hPhiDiff.differentiableAt
        (isOpen_ctRightHalfPlane27.mem_nhds hp)).continuousAt.comp
          (by
            unfold verticalPoint27
            fun_prop)
    · exact (continuous_const.pow 2).mul continuous_sechSq27
  · filter_upwards with y
    have hpointRe : 1 / 2 ≤ (verticalPoint27 x y).re := by
      simpa [verticalPoint27] using hx
    have hnorm := norm_verticalPoint_le_abs_add27 x y
    have hbase : 1 + ‖verticalPoint27 x y‖ ≤
        (1 + |x|) * (1 + |y|) := by
      nlinarith [abs_nonneg x, abs_nonneg y]
    have hpow : (1 + ‖verticalPoint27 x y‖) ^ d ≤
        ((1 + |x|) * (1 + |y|)) ^ d := by
      exact pow_le_pow_left₀ (by positivity) hbase d
    have hPhi : ‖Phi (verticalPoint27 x y)‖ ≤
        C * ((1 + |x|) * (1 + |y|)) ^ d :=
      (hbound _ hpointRe).trans
        (mul_le_mul_of_nonneg_left hpow hC)
    rw [norm_mul, norm_mul, norm_pow, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos Real.pi_pos]
    calc
      ‖Phi (verticalPoint27 x y)‖ *
          (Real.pi ^ 2 * ‖sechSq27 y‖) ≤
          (C * ((1 + |x|) * (1 + |y|)) ^ d) *
            (Real.pi ^ 2 * ‖sechSq27 y‖) := by
        exact mul_le_mul_of_nonneg_right hPhi (by positivity)
      _ = K * ((1 + |y|) ^ d * ‖sechSq27 y‖) := by
        dsimp only [K]
        rw [mul_pow]
        ring

theorem integrable_ctRKernelRaw_left27
    (n m : ℕ) (hm1 : 1 ≤ m) :
    Integrable (fun y : ℝ =>
      ctRKernelRaw27 n (verticalPoint27 ((m : ℝ) - 1 / 2) y)) := by
  have hx : 1 / 2 ≤ (m : ℝ) - 1 / 2 := by
    have hm1' : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm1
    linarith
  have h := integrable_phi_mul_sechSq27
    (ctRPhi_differentiableOn_right27 n) (ctRPhi_poly_bound27 n) hx
  simpa only [ctRKernelRaw27, squaredSineKernel_left_halfLine27] using h

theorem integrable_ctRKernelRaw_right27 (n m : ℕ) :
    Integrable (fun y : ℝ =>
      ctRKernelRaw27 n (verticalPoint27 ((m : ℝ) + 1 / 2) y)) := by
  have hx : 1 / 2 ≤ (m : ℝ) + 1 / 2 := by
    have hm0 : (0 : ℝ) ≤ (m : ℝ) := Nat.cast_nonneg m
    linarith
  have h := integrable_phi_mul_sechSq27
    (ctRPhi_differentiableOn_right27 n) (ctRPhi_poly_bound27 n) hx
  simpa only [ctRKernelRaw27, squaredSineKernel_right_halfLine27] using h

theorem integrable_ctSKernelRaw_left27
    (n m : ℕ) (hm1 : 1 ≤ m) :
    Integrable (fun y : ℝ =>
      ctSKernelRaw27 n (verticalPoint27 ((m : ℝ) - 1 / 2) y)) := by
  have hx : 1 / 2 ≤ (m : ℝ) - 1 / 2 := by
    have hm1' : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm1
    linarith
  have h := integrable_phi_mul_sechSq27
    (ctSPhi_differentiableOn_right27 n) (ctSPhi_poly_bound27 n) hx
  simpa only [ctSKernelRaw27, squaredSineKernel_left_halfLine27] using h

theorem integrable_ctSKernelRaw_right27 (n m : ℕ) :
    Integrable (fun y : ℝ =>
      ctSKernelRaw27 n (verticalPoint27 ((m : ℝ) + 1 / 2) y)) := by
  have hx : 1 / 2 ≤ (m : ℝ) + 1 / 2 := by
    have hm0 : (0 : ℝ) ≤ (m : ℝ) := Nat.cast_nonneg m
    linarith
  have h := integrable_phi_mul_sechSq27
    (ctSPhi_differentiableOn_right27 n) (ctSPhi_poly_bound27 n) hx
  simpa only [ctSKernelRaw27, squaredSineKernel_right_halfLine27] using h

theorem integrable_ctRKernelExt_left27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Integrable (fun y : ℝ =>
      ctRKernelExt27 n m (verticalPoint27 ((m : ℝ) - 1 / 2) y)) := by
  apply (integrable_ctRKernelRaw_left27 n m hm1).congr
  filter_upwards with y
  have hz : verticalPoint27 ((m : ℝ) - 1 / 2) y ∈ ctClosedStrip27 m := by
    norm_num [ctClosedStrip27, verticalPoint27]
    linarith
  have hne : verticalPoint27 ((m : ℝ) - 1 / 2) y ≠ (m : ℂ) := by
    intro h
    have hre := congrArg Complex.re h
    simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
      Complex.I_im, Complex.natCast_re] at hre
    norm_num at hre
  exact ctRKernelRaw_eq_ext27 hm1 hmn hz hne

theorem integrable_ctRKernelExt_right27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Integrable (fun y : ℝ =>
      ctRKernelExt27 n m (verticalPoint27 ((m : ℝ) + 1 / 2) y)) := by
  apply (integrable_ctRKernelRaw_right27 n m).congr
  filter_upwards with y
  have hz : verticalPoint27 ((m : ℝ) + 1 / 2) y ∈ ctClosedStrip27 m := by
    norm_num [ctClosedStrip27, verticalPoint27]
    linarith
  have hne : verticalPoint27 ((m : ℝ) + 1 / 2) y ≠ (m : ℂ) := by
    intro h
    have hre := congrArg Complex.re h
    simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
      Complex.I_im, Complex.natCast_re] at hre
    norm_num at hre
  exact ctRKernelRaw_eq_ext27 hm1 hmn hz hne

theorem integrable_ctSKernelExt_left27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) :
    Integrable (fun y : ℝ =>
      ctSKernelExt27 n m (verticalPoint27 ((m : ℝ) - 1 / 2) y)) := by
  apply (integrable_ctSKernelRaw_left27 n m hm1).congr
  filter_upwards with y
  have hz : verticalPoint27 ((m : ℝ) - 1 / 2) y ∈ ctClosedStrip27 m := by
    norm_num [ctClosedStrip27, verticalPoint27]
    linarith
  have hne : verticalPoint27 ((m : ℝ) - 1 / 2) y ≠ (m : ℂ) := by
    intro h
    have hre := congrArg Complex.re h
    simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
      Complex.I_im, Complex.natCast_re] at hre
    norm_num at hre
  exact ctSKernelRaw_eq_ext27 hm1 hmn hz hne

theorem integrable_ctSKernelExt_right27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) :
    Integrable (fun y : ℝ =>
      ctSKernelExt27 n m (verticalPoint27 ((m : ℝ) + 1 / 2) y)) := by
  apply (integrable_ctSKernelRaw_right27 n m).congr
  filter_upwards with y
  have hz : verticalPoint27 ((m : ℝ) + 1 / 2) y ∈ ctClosedStrip27 m := by
    norm_num [ctClosedStrip27, verticalPoint27]
    linarith
  have hne : verticalPoint27 ((m : ℝ) + 1 / 2) y ≠ (m : ℂ) := by
    intro h
    have hre := congrArg Complex.re h
    simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
      Complex.I_im, Complex.natCast_re] at hre
    norm_num at hre
  exact ctSKernelRaw_eq_ext27 hm1 hmn hz hne

/-! ## Vanishing horizontal boundaries -/

private theorem rawKernel_horizontal_bound27
    {Phi : ℂ → ℂ} {C : ℝ} {d m : ℕ}
    (hC : 0 ≤ C)
    (hPhi : ∀ z : ℂ, 1 / 2 ≤ z.re →
      ‖Phi z‖ ≤ C * (1 + ‖z‖) ^ d)
    (hm1 : 1 ≤ m) {T x u : ℝ}
    (hT : 1 ≤ T)
    (hx : x ∈ [[(m : ℝ) - 1 / 2, (m : ℝ) + 1 / 2]])
    (hu : |u| = T) :
    ‖Phi (verticalPoint27 x u) *
        zudilinBarnesSquaredSineKernel27 (verticalPoint27 x u)‖ ≤
      (C * (16 * Real.pi ^ 2) * ((m : ℝ) + 3) ^ d) *
        (1 * T) ^ d * Real.exp (-(1 * T)) := by
  have hab : (m : ℝ) - 1 / 2 ≤ (m : ℝ) + 1 / 2 := by linarith
  rw [uIcc_of_le hab] at hx
  have hm1' : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm1
  have hx0 : 0 ≤ x := by linarith [hx.1]
  have hpointRe : 1 / 2 ≤ (verticalPoint27 x u).re := by
    simpa [verticalPoint27] using (show (1 / 2 : ℝ) ≤ x by linarith [hx.1])
  have hnorm0 := norm_verticalPoint_le_abs_add27 x u
  have hnorm : ‖verticalPoint27 x u‖ ≤ x + T := by
    rw [abs_of_nonneg hx0, hu] at hnorm0
    exact hnorm0
  have hprod : 0 ≤ (T - 1) * ((m : ℝ) + 2) :=
    mul_nonneg (by linarith) (by positivity)
  have hbase : 1 + ‖verticalPoint27 x u‖ ≤ ((m : ℝ) + 3) * T := by
    nlinarith [hx.2]
  have hpow : (1 + ‖verticalPoint27 x u‖) ^ d ≤
      (((m : ℝ) + 3) * T) ^ d :=
    pow_le_pow_left₀ (by positivity) hbase d
  have hPhi' : ‖Phi (verticalPoint27 x u)‖ ≤
      C * (((m : ℝ) + 3) * T) ^ d :=
    (hPhi _ hpointRe).trans (mul_le_mul_of_nonneg_left hpow hC)
  have hkernel : ‖zudilinBarnesSquaredSineKernel27
      (verticalPoint27 x u)‖ ≤
      16 * Real.pi ^ 2 * Real.exp (-2 * Real.pi * T) := by
    have him : |(verticalPoint27 x u).im| = T := by
      simpa [verticalPoint27] using hu
    simpa only [him] using
      (norm_squaredSineKernel_le_exp27 (z := verticalPoint27 x u)
        (by simpa only [him] using hT))
  have hexp : Real.exp (-2 * Real.pi * T) ≤ Real.exp (-T) := by
    apply Real.exp_le_exp.mpr
    nlinarith [Real.two_le_pi]
  rw [norm_mul]
  calc
    ‖Phi (verticalPoint27 x u)‖ *
        ‖zudilinBarnesSquaredSineKernel27 (verticalPoint27 x u)‖ ≤
        (C * (((m : ℝ) + 3) * T) ^ d) *
          (16 * Real.pi ^ 2 * Real.exp (-2 * Real.pi * T)) := by
      exact mul_le_mul hPhi' hkernel (norm_nonneg _) (by positivity)
    _ ≤ (C * (((m : ℝ) + 3) * T) ^ d) *
          (16 * Real.pi ^ 2 * Real.exp (-T)) := by
      gcongr
    _ = (C * (16 * Real.pi ^ 2) * ((m : ℝ) + 3) ^ d) *
          (1 * T) ^ d * Real.exp (-(1 * T)) := by
      rw [mul_pow]
      ring

private theorem ctRKernelRaw_horizontal_top_tendsto27
    (n : ℕ) {m : ℕ} (hm1 : 1 ≤ m) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctRKernelRaw27 n ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rcases ctRPhi_poly_bound27 n with ⟨C, d, hC, hPhi⟩
  apply horizontalIntegral_tendsto_zero_of_pow_exp_bound27
    (C := C * (16 * Real.pi ^ 2) * ((m : ℝ) + 3) ^ d)
    (c := 1) (d := d) (by positivity) (by norm_num)
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  intro x hx
  simpa only [ctRKernelRaw27, verticalPoint27] using
    rawKernel_horizontal_bound27 hC hPhi hm1 hT hx (abs_of_nonneg (by linarith))

private theorem ctRKernelRaw_horizontal_bottom_tendsto27
    (n : ℕ) {m : ℕ} (hm1 : 1 ≤ m) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctRKernelRaw27 n ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rcases ctRPhi_poly_bound27 n with ⟨C, d, hC, hPhi⟩
  let F : ℂ → ℂ := fun z => ctRKernelRaw27 n (conj z)
  have htop : Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        F ((x : ℂ) + (T : ℂ) * Complex.I)) atTop (𝓝 0) := by
    apply horizontalIntegral_tendsto_zero_of_pow_exp_bound27
      (C := C * (16 * Real.pi ^ 2) * ((m : ℝ) + 3) ^ d)
      (c := 1) (d := d) (by positivity) (by norm_num)
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
    intro x hx
    have hu : |-T| = T := by rw [abs_neg, abs_of_nonneg (by linarith)]
    have hconj : conj ((x : ℂ) + (T : ℂ) * Complex.I) =
        verticalPoint27 x (-T) := by
      unfold verticalPoint27
      simp only [map_add, map_mul, Complex.conj_ofReal,
        Complex.conj_I, Complex.ofReal_neg]
      ring
    change ‖ctRKernelRaw27 n (conj ((x : ℂ) + (T : ℂ) * Complex.I))‖ ≤ _
    rw [hconj]
    simpa only [ctRKernelRaw27] using
      rawKernel_horizontal_bound27 hC hPhi hm1 hT hx hu
  have hconj (x T : ℝ) : conj ((x : ℂ) + (T : ℂ) * Complex.I) =
      (x : ℂ) - (T : ℂ) * Complex.I := by
    simp only [map_add, map_mul, Complex.conj_ofReal, Complex.conj_I]
    ring
  simpa only [F, hconj] using htop

private theorem ctSKernelRaw_horizontal_top_tendsto27
    (n : ℕ) {m : ℕ} (hm1 : 1 ≤ m) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctSKernelRaw27 n ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rcases ctSPhi_poly_bound27 n with ⟨C, d, hC, hPhi⟩
  apply horizontalIntegral_tendsto_zero_of_pow_exp_bound27
    (C := C * (16 * Real.pi ^ 2) * ((m : ℝ) + 3) ^ d)
    (c := 1) (d := d) (by positivity) (by norm_num)
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  intro x hx
  simpa only [ctSKernelRaw27, verticalPoint27] using
    rawKernel_horizontal_bound27 hC hPhi hm1 hT hx (abs_of_nonneg (by linarith))

private theorem ctSKernelRaw_horizontal_bottom_tendsto27
    (n : ℕ) {m : ℕ} (hm1 : 1 ≤ m) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctSKernelRaw27 n ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  rcases ctSPhi_poly_bound27 n with ⟨C, d, hC, hPhi⟩
  let F : ℂ → ℂ := fun z => ctSKernelRaw27 n (conj z)
  have htop : Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        F ((x : ℂ) + (T : ℂ) * Complex.I)) atTop (𝓝 0) := by
    apply horizontalIntegral_tendsto_zero_of_pow_exp_bound27
      (C := C * (16 * Real.pi ^ 2) * ((m : ℝ) + 3) ^ d)
      (c := 1) (d := d) (by positivity) (by norm_num)
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
    intro x hx
    have hu : |-T| = T := by rw [abs_neg, abs_of_nonneg (by linarith)]
    have hconj : conj ((x : ℂ) + (T : ℂ) * Complex.I) =
        verticalPoint27 x (-T) := by
      unfold verticalPoint27
      simp only [map_add, map_mul, Complex.conj_ofReal,
        Complex.conj_I, Complex.ofReal_neg]
      ring
    change ‖ctSKernelRaw27 n (conj ((x : ℂ) + (T : ℂ) * Complex.I))‖ ≤ _
    rw [hconj]
    simpa only [ctSKernelRaw27] using
      rawKernel_horizontal_bound27 hC hPhi hm1 hT hx hu
  have hconj (x T : ℝ) : conj ((x : ℂ) + (T : ℂ) * Complex.I) =
      (x : ℂ) - (T : ℂ) * Complex.I := by
    simp only [map_add, map_mul, Complex.conj_ofReal, Complex.conj_I]
    ring
  simpa only [F, hconj] using htop

private theorem ctRKernelExt_horizontal_top_tendsto27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctRKernelExt27 n m ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  apply Filter.Tendsto.congr' ?_ (ctRKernelRaw_horizontal_top_tendsto27 n hm1)
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  apply intervalIntegral.integral_congr
  intro x hx
  have hz : (x : ℂ) + (T : ℂ) * Complex.I ∈ ctClosedStrip27 m := by
    have hab : (m : ℝ) - 1 / 2 ≤ (m : ℝ) + 1 / 2 := by linarith
    rw [uIcc_of_le hab] at hx
    simp only [ctClosedStrip27, Complex.add_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
    norm_num
    exact ⟨by linarith [hx.1], hx.2⟩
  have hne : (x : ℂ) + (T : ℂ) * Complex.I ≠ (m : ℂ) := by
    intro h
    have him := congrArg Complex.im h
    simp only [Complex.add_im, Complex.mul_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im,
      Complex.natCast_im] at him
    norm_num at him
    linarith
  exact ctRKernelRaw_eq_ext27 hm1 hmn hz hne

private theorem ctRKernelExt_horizontal_bottom_tendsto27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctRKernelExt27 n m ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  apply Filter.Tendsto.congr' ?_ (ctRKernelRaw_horizontal_bottom_tendsto27 n hm1)
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  apply intervalIntegral.integral_congr
  intro x hx
  have hz : (x : ℂ) - (T : ℂ) * Complex.I ∈ ctClosedStrip27 m := by
    have hab : (m : ℝ) - 1 / 2 ≤ (m : ℝ) + 1 / 2 := by linarith
    rw [uIcc_of_le hab] at hx
    simp only [ctClosedStrip27, Complex.sub_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
    norm_num
    exact ⟨by linarith [hx.1], hx.2⟩
  have hne : (x : ℂ) - (T : ℂ) * Complex.I ≠ (m : ℂ) := by
    intro h
    have him := congrArg Complex.im h
    simp only [Complex.sub_im, Complex.mul_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im,
      Complex.natCast_im] at him
    norm_num at him
    linarith
  exact ctRKernelRaw_eq_ext27 hm1 hmn hz hne

private theorem ctSKernelExt_horizontal_top_tendsto27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctSKernelExt27 n m ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  apply Filter.Tendsto.congr' ?_ (ctSKernelRaw_horizontal_top_tendsto27 n hm1)
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  apply intervalIntegral.integral_congr
  intro x hx
  have hz : (x : ℂ) + (T : ℂ) * Complex.I ∈ ctClosedStrip27 m := by
    have hab : (m : ℝ) - 1 / 2 ≤ (m : ℝ) + 1 / 2 := by linarith
    rw [uIcc_of_le hab] at hx
    simp only [ctClosedStrip27, Complex.add_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
    norm_num
    exact ⟨by linarith [hx.1], hx.2⟩
  have hne : (x : ℂ) + (T : ℂ) * Complex.I ≠ (m : ℂ) := by
    intro h
    have him := congrArg Complex.im h
    simp only [Complex.add_im, Complex.mul_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im,
      Complex.natCast_im] at him
    norm_num at him
    linarith
  exact ctSKernelRaw_eq_ext27 hm1 hmn hz hne

private theorem ctSKernelExt_horizontal_bottom_tendsto27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) :
    Tendsto
      (fun T : ℝ => ∫ x in (m : ℝ) - 1 / 2..(m : ℝ) + 1 / 2,
        ctSKernelExt27 n m ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  apply Filter.Tendsto.congr' ?_ (ctSKernelRaw_horizontal_bottom_tendsto27 n hm1)
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with T hT
  apply intervalIntegral.integral_congr
  intro x hx
  have hz : (x : ℂ) - (T : ℂ) * Complex.I ∈ ctClosedStrip27 m := by
    have hab : (m : ℝ) - 1 / 2 ≤ (m : ℝ) + 1 / 2 := by linarith
    rw [uIcc_of_le hab] at hx
    simp only [ctClosedStrip27, Complex.sub_re, Complex.mul_re,
      Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
    norm_num
    exact ⟨by linarith [hx.1], hx.2⟩
  have hne : (x : ℂ) - (T : ℂ) * Complex.I ≠ (m : ℂ) := by
    intro h
    have him := congrArg Complex.im h
    simp only [Complex.sub_im, Complex.mul_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im,
      Complex.natCast_im] at him
    norm_num at him
    linarith
  exact ctSKernelRaw_eq_ext27 hm1 hmn hz hne

theorem ctRKernelRaw_one_strip27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m ≤ n) :
    (∫ y : ℝ,
      ctRKernelRaw27 n (verticalPoint27 ((m : ℝ) - 1 / 2) y)) =
      ∫ y : ℝ,
        ctRKernelRaw27 n (verticalPoint27 ((m : ℝ) + 1 / 2) y) := by
  have hext := verticalIntegral_eq_of_horizontal_tendsto27
    (F := ctRKernelExt27 n m)
    (a := (m : ℝ) - 1 / 2) (b := (m : ℝ) + 1 / 2)
    (by linarith)
    (ctRKernelExt_differentiableOn_strip27 hm1)
    (integrable_ctRKernelExt_left27 hm1 hmn)
    (integrable_ctRKernelExt_right27 hm1 hmn)
    (ctRKernelExt_horizontal_top_tendsto27 hm1 hmn)
    (ctRKernelExt_horizontal_bottom_tendsto27 hm1 hmn)
  calc
    (∫ y : ℝ,
      ctRKernelRaw27 n (verticalPoint27 ((m : ℝ) - 1 / 2) y)) =
        ∫ y : ℝ,
          ctRKernelExt27 n m (verticalPoint27 ((m : ℝ) - 1 / 2) y) := by
      apply integral_congr_ae
      filter_upwards with y
      apply ctRKernelRaw_eq_ext27 hm1 hmn
      · norm_num [ctClosedStrip27, verticalPoint27]
        linarith
      · intro h
        have hre := congrArg Complex.re h
        simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
          Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
          Complex.I_im, Complex.natCast_re] at hre
        norm_num at hre
    _ = ∫ y : ℝ,
          ctRKernelExt27 n m (verticalPoint27 ((m : ℝ) + 1 / 2) y) := hext
    _ = ∫ y : ℝ,
          ctRKernelRaw27 n (verticalPoint27 ((m : ℝ) + 1 / 2) y) := by
      apply integral_congr_ae
      filter_upwards with y
      symm
      apply ctRKernelRaw_eq_ext27 hm1 hmn
      · norm_num [ctClosedStrip27, verticalPoint27]
        linarith
      · intro h
        have hre := congrArg Complex.re h
        simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
          Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
          Complex.I_im, Complex.natCast_re] at hre
        norm_num at hre

theorem ctSKernelRaw_one_strip27
    {n m : ℕ} (hm1 : 1 ≤ m) (hmn : m < n) :
    (∫ y : ℝ,
      ctSKernelRaw27 n (verticalPoint27 ((m : ℝ) - 1 / 2) y)) =
      ∫ y : ℝ,
        ctSKernelRaw27 n (verticalPoint27 ((m : ℝ) + 1 / 2) y) := by
  have hext := verticalIntegral_eq_of_horizontal_tendsto27
    (F := ctSKernelExt27 n m)
    (a := (m : ℝ) - 1 / 2) (b := (m : ℝ) + 1 / 2)
    (by linarith)
    (ctSKernelExt_differentiableOn_strip27 hm1)
    (integrable_ctSKernelExt_left27 hm1 hmn)
    (integrable_ctSKernelExt_right27 hm1 hmn)
    (ctSKernelExt_horizontal_top_tendsto27 hm1 hmn)
    (ctSKernelExt_horizontal_bottom_tendsto27 hm1 hmn)
  calc
    (∫ y : ℝ,
      ctSKernelRaw27 n (verticalPoint27 ((m : ℝ) - 1 / 2) y)) =
        ∫ y : ℝ,
          ctSKernelExt27 n m (verticalPoint27 ((m : ℝ) - 1 / 2) y) := by
      apply integral_congr_ae
      filter_upwards with y
      apply ctSKernelRaw_eq_ext27 hm1 hmn
      · norm_num [ctClosedStrip27, verticalPoint27]
        linarith
      · intro h
        have hre := congrArg Complex.re h
        simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
          Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
          Complex.I_im, Complex.natCast_re] at hre
        norm_num at hre
    _ = ∫ y : ℝ,
          ctSKernelExt27 n m (verticalPoint27 ((m : ℝ) + 1 / 2) y) := hext
    _ = ∫ y : ℝ,
          ctSKernelRaw27 n (verticalPoint27 ((m : ℝ) + 1 / 2) y) := by
      apply integral_congr_ae
      filter_upwards with y
      symm
      apply ctSKernelRaw_eq_ext27 hm1 hmn
      · norm_num [ctClosedStrip27, verticalPoint27]
        linarith
      · intro h
        have hre := congrArg Complex.re h
        simp only [verticalPoint27, Complex.add_re, Complex.mul_re,
          Complex.ofReal_re, Complex.ofReal_im, Complex.I_re,
          Complex.I_im, Complex.natCast_re] at hre
        norm_num at hre

end RamanujanChallenge.P27
