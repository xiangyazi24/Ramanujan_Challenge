import RamanujanChallenge.Problem25DualCertificateSparseRow0

noncomputable section

namespace RamanujanChallenge.P25

def spCoeffEval (n : ℝ) : List ℤ → ℝ
  | [] => 0
  | c :: cs => c + n * spCoeffEval n cs

theorem spCoeffEval_normalize (n : ℝ) (cs : List ℤ) :
    spCoeffEval n (spCoeffNormalize cs) = spCoeffEval n cs := by
  induction cs with
  | nil => rfl
  | cons c cs ih =>
      cases h : spCoeffNormalize cs with
      | nil =>
          rw [h] at ih
          simp only [spCoeffEval] at ih
          split_ifs with hc
          · simp [spCoeffNormalize, h, spCoeffEval, hc, ← ih]
          · simp [spCoeffNormalize, h, spCoeffEval, hc, ← ih]
      | cons d ds =>
          rw [h] at ih
          simp only [spCoeffEval] at ih
          simp only [spCoeffNormalize, h, spCoeffEval]
          rw [ih]

theorem spCoeffEval_addRaw (n : ℝ) (as bs : List ℤ) :
    spCoeffEval n (spCoeffAddRaw as bs) =
      spCoeffEval n as + spCoeffEval n bs := by
  induction as generalizing bs with
  | nil => simp [spCoeffAddRaw, spCoeffEval]
  | cons a as ih =>
      cases bs with
      | nil => simp [spCoeffAddRaw, spCoeffEval]
      | cons b bs =>
          simp [spCoeffAddRaw, spCoeffEval, ih]
          ring

@[simp] theorem spCoeffEval_add (n : ℝ) (as bs : List ℤ) :
    spCoeffEval n (spCoeffAdd as bs) =
      spCoeffEval n as + spCoeffEval n bs := by
  rw [spCoeffAdd, spCoeffEval_normalize, spCoeffEval_addRaw]

theorem spCoeffEval_map_neg (n : ℝ) (as : List ℤ) :
    spCoeffEval n (as.map (-·)) = -spCoeffEval n as := by
  induction as with
  | nil => simp [spCoeffEval]
  | cons a as ih =>
      simp [spCoeffEval, ih]
      ring

@[simp] theorem spCoeffEval_neg (n : ℝ) (as : List ℤ) :
    spCoeffEval n (spCoeffNeg as) = -spCoeffEval n as := by
  rw [spCoeffNeg, spCoeffEval_normalize, spCoeffEval_map_neg]

@[simp] theorem spCoeffEval_scale (n : ℝ) (a : ℤ) (bs : List ℤ) :
    spCoeffEval n (spCoeffScale a bs) = a * spCoeffEval n bs := by
  induction bs with
  | nil => simp [spCoeffScale, spCoeffEval]
  | cons b bs ih =>
      simp only [spCoeffScale, List.map_cons, spCoeffEval]
      rw [ih]
      push_cast
      ring

theorem spCoeffEval_mulRaw (n : ℝ) (as bs : List ℤ) :
    spCoeffEval n (spCoeffMulRaw as bs) =
      spCoeffEval n as * spCoeffEval n bs := by
  induction as with
  | nil => simp [spCoeffMulRaw, spCoeffEval]
  | cons a as ih =>
      simp [spCoeffMulRaw, spCoeffEval_addRaw, spCoeffEval_scale,
        spCoeffEval, ih]
      ring

@[simp] theorem spCoeffEval_mul (n : ℝ) (as bs : List ℤ) :
    spCoeffEval n (spCoeffMul as bs) =
      spCoeffEval n as * spCoeffEval n bs := by
  rw [spCoeffMul, spCoeffEval_normalize, spCoeffEval_mulRaw]

def sparseTermEval (n p q v : ℝ) (t : SparseTerm) : ℝ :=
  spCoeffEval n t.nCoeffs * p ^ t.pExp * q ^ t.qExp * v ^ t.vExp

def sparseEval (P : SparsePoly) (n p q v : ℝ) : ℝ :=
  (P.terms.map (sparseTermEval n p q v)).sum

private theorem sparseExpCompare_eq_iff (a b : SparseTerm) :
    dualCertExpCompare (sparseExp a) (sparseExp b) = .eq ↔
      sparseExp a = sparseExp b := by
  constructor
  · intro h
    by_contra hab
    simp only [dualCertExpCompare, if_neg hab] at h
    split at h <;> simp_all
  · intro h
    simp [dualCertExpCompare, h]

private theorem sparseTermEval_add_coeff
    (n p q v : ℝ) (a b : SparseTerm) (h : sparseExp a = sparseExp b) :
    sparseTermEval n p q v { a with nCoeffs := spCoeffAdd a.nCoeffs b.nCoeffs } =
      sparseTermEval n p q v a + sparseTermEval n p q v b := by
  rcases a with ⟨ac, ap, aq, av⟩
  rcases b with ⟨bc, bp, bq, bv⟩
  simp only [sparseExp] at h
  injection h with _ hp hq hv
  subst bp
  subst bq
  subst bv
  simp [sparseTermEval]
  ring

private theorem sparseEval_insert (n p q v : ℝ) (t : SparseTerm) :
    ∀ ts : List SparseTerm,
      ((sparseInsert t ts).map (sparseTermEval n p q v)).sum =
        sparseTermEval n p q v t +
          (ts.map (sparseTermEval n p q v)).sum := by
  intro ts
  induction ts with
  | nil =>
      simp only [sparseInsert]
      split_ifs with h
      · simp [sparseTermEval, h, spCoeffEval]
      · simp
  | cons u us ih =>
      simp only [sparseInsert]
      cases hcmp : dualCertExpCompare (sparseExp t) (sparseExp u) with
      | lt =>
          simp only [hcmp]
          split_ifs with h
          · simp [sparseTermEval, h, spCoeffEval]
          · simp
      | eq =>
          simp only [hcmp]
          have hexp : sparseExp t = sparseExp u :=
            (sparseExpCompare_eq_iff t u).mp hcmp
          split_ifs with hsum
          · have heval :
                sparseTermEval n p q v t + sparseTermEval n p q v u = 0 := by
              rw [← sparseTermEval_add_coeff n p q v t u hexp]
              simp [sparseTermEval, hsum, spCoeffEval]
            simp only [List.map_cons, List.sum_cons]
            linarith
          · simp only [List.map_cons, List.sum_cons]
            rw [sparseTermEval_add_coeff n p q v t u hexp]
            ring
      | gt =>
          simp only [hcmp, List.map_cons, List.sum_cons, ih]
          ring

theorem sparseEval_normalize (P : SparsePoly) (n p q v : ℝ) :
    sparseEval (sparseNormalize P) n p q v = sparseEval P n p q v := by
  rcases P with ⟨ts⟩
  unfold sparseNormalize sparseEval
  induction ts with
  | nil => simp
  | cons t ts ih =>
      simp only [List.foldr_cons, List.map_cons, List.sum_cons]
      rw [sparseEval_insert, ih]

@[simp] theorem sparseEval_const (z : ℤ) (n p q v : ℝ) :
    sparseEval (SparsePoly.const z) n p q v = z := by
  rw [SparsePoly.const, sparseEval_normalize]
  by_cases hz : z = 0
  · simp [sparseEval, sparseTermEval, spCoeffEval, hz]
  · simp [sparseEval, sparseTermEval, spCoeffEval, hz]

@[simp] theorem sparseEval_n (n p q v : ℝ) :
    sparseEval SparsePoly.n n p q v = n := by
  simp [SparsePoly.n, sparseEval, sparseTermEval, spCoeffEval]

@[simp] theorem sparseEval_p (n p q v : ℝ) :
    sparseEval SparsePoly.p n p q v = p := by
  simp [SparsePoly.p, sparseEval, sparseTermEval, spCoeffEval]

@[simp] theorem sparseEval_q (n p q v : ℝ) :
    sparseEval SparsePoly.q n p q v = q := by
  simp [SparsePoly.q, sparseEval, sparseTermEval, spCoeffEval]

@[simp] theorem sparseEval_v (n p q v : ℝ) :
    sparseEval SparsePoly.v n p q v = v := by
  simp [SparsePoly.v, sparseEval, sparseTermEval, spCoeffEval]

@[simp] theorem sparseEval_ofNat (k : ℕ) (n p q v : ℝ) :
    sparseEval (OfNat.ofNat k : SparsePoly) n p q v = k := by
  rw [sparse_ofNat_eq, sparseEval_const]
  norm_num

@[simp] theorem sparseEval_add (P Q : SparsePoly) (n p q v : ℝ) :
    sparseEval (P + Q) n p q v =
      sparseEval P n p q v + sparseEval Q n p q v := by
  rw [sparse_add_eq, sparseEval_normalize]
  simp [sparseEval]

@[simp] theorem sparseEval_neg (P : SparsePoly) (n p q v : ℝ) :
    sparseEval (-P) n p q v = -sparseEval P n p q v := by
  rw [sparse_neg_eq, sparseEval_normalize]
  rcases P with ⟨ts⟩
  induction ts with
  | nil => simp [sparseEval]
  | cons t ts ih =>
      simp only [sparseEval, List.map_cons, List.sum_cons] at ih ⊢
      simp [sparseTermEval, ih]
      ring

@[simp] theorem sparseEval_sub (P Q : SparsePoly) (n p q v : ℝ) :
    sparseEval (P - Q) n p q v =
      sparseEval P n p q v - sparseEval Q n p q v := by
  rw [sparse_sub_eq, sparseEval_add, sparseEval_neg]
  ring

private theorem sparseTermEval_mul (n p q v : ℝ) (a b : SparseTerm) :
    sparseTermEval n p q v (sparseMulTerm a b) =
      sparseTermEval n p q v a * sparseTermEval n p q v b := by
  rcases a with ⟨ac, ap, aq, av⟩
  rcases b with ⟨bc, bp, bq, bv⟩
  simp [sparseMulTerm, sparseTermEval, pow_add]
  ring

private theorem sparseEval_map_mul_left
    (a : SparseTerm) (qs : List SparseTerm) (n p q v : ℝ) :
    ((qs.map (sparseMulTerm a)).map (sparseTermEval n p q v)).sum =
      sparseTermEval n p q v a *
        (qs.map (sparseTermEval n p q v)).sum := by
  induction qs with
  | nil => simp
  | cons b qs ih =>
      simp only [List.map_cons, List.sum_cons, ih, sparseTermEval_mul]
      ring

private theorem sparseEval_mul_raw (P Q : SparsePoly) (n p q v : ℝ) :
    ((P.terms.flatMap fun a => Q.terms.map (sparseMulTerm a)).map
        (sparseTermEval n p q v)).sum =
      sparseEval P n p q v * sparseEval Q n p q v := by
  rcases P with ⟨ps⟩
  rcases Q with ⟨qs⟩
  simp only [sparseEval]
  induction ps with
  | nil => simp
  | cons a ps ih =>
      simp only [List.flatMap_cons, List.map_append, List.sum_append,
        List.map_cons, List.sum_cons, ih]
      rw [sparseEval_map_mul_left]
      ring

@[simp] theorem sparseEval_mul (P Q : SparsePoly) (n p q v : ℝ) :
    sparseEval (P * Q) n p q v =
      sparseEval P n p q v * sparseEval Q n p q v := by
  rw [sparse_mul_eq, sparseEval_normalize]
  exact sparseEval_mul_raw P Q n p q v

@[simp] theorem sparseEval_pow (P : SparsePoly) (k : ℕ) (n p q v : ℝ) :
    sparseEval (P ^ k) n p q v = sparseEval P n p q v ^ k := by
  induction k with
  | zero => simp
  | succ k ih =>
      change sparseEval ((P ^ k) * P) n p q v =
        sparseEval P n p q v ^ k * sparseEval P n p q v
      rw [sparseEval_mul, ih]

end RamanujanChallenge.P25

end
