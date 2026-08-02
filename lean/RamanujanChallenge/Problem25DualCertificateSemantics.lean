import RamanujanChallenge.Problem25DualCertificateCommon
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Pow

noncomputable section

namespace RamanujanChallenge.P25

private theorem dualCertPDerivMonomial_hasDerivAt
    (t : DualCertMonomial) (n q v p : ℝ) :
    HasDerivAt (fun x => dualCertMonomialEval n x q v t)
      (((dualCertPDerivMonomial t).map
        (dualCertMonomialEval n p q v)).sum) p := by
  rcases t with ⟨c, ⟨ne, pe, qe, ve⟩⟩
  cases pe with
  | zero =>
      simpa [dualCertPDerivMonomial, dualCertMonomialEval] using
        (hasDerivAt_const p ((c : ℝ) * n ^ ne * 1 * q ^ qe * v ^ ve))
  | succ pe =>
      have h := (hasDerivAt_pow (pe + 1) p).const_mul
        ((c : ℝ) * n ^ ne * q ^ qe * v ^ ve)
      convert h using 1 <;> (try funext x) <;>
        simp [dualCertPDerivMonomial, dualCertMonomialEval, Nat.cast_add,
          Nat.cast_one, Nat.succ_eq_add_one] <;> ring

private theorem dualCertQDerivMonomial_hasDerivAt
    (t : DualCertMonomial) (n p v q : ℝ) :
    HasDerivAt (fun x => dualCertMonomialEval n p x v t)
      (((dualCertQDerivMonomial t).map
        (dualCertMonomialEval n p q v)).sum) q := by
  rcases t with ⟨c, ⟨ne, pe, qe, ve⟩⟩
  cases qe with
  | zero =>
      simpa [dualCertQDerivMonomial, dualCertMonomialEval] using
        (hasDerivAt_const q ((c : ℝ) * n ^ ne * p ^ pe * 1 * v ^ ve))
  | succ qe =>
      have h := (hasDerivAt_pow (qe + 1) q).const_mul
        ((c : ℝ) * n ^ ne * p ^ pe * v ^ ve)
      convert h using 1 <;> (try funext x) <;>
        simp [dualCertQDerivMonomial, dualCertMonomialEval, Nat.cast_add,
          Nat.cast_one, Nat.succ_eq_add_one] <;> ring

private theorem dualCertVDerivMonomial_hasDerivAt
    (t : DualCertMonomial) (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertMonomialEval n p q x t)
      (((dualCertVDerivMonomial t).map
        (dualCertMonomialEval n p q v)).sum) v := by
  rcases t with ⟨c, ⟨ne, pe, qe, ve⟩⟩
  cases ve with
  | zero =>
      simpa [dualCertVDerivMonomial, dualCertMonomialEval] using
        (hasDerivAt_const v ((c : ℝ) * n ^ ne * p ^ pe * q ^ qe * 1))
  | succ ve =>
      have h := (hasDerivAt_pow (ve + 1) v).const_mul
        ((c : ℝ) * n ^ ne * p ^ pe * q ^ qe)
      convert h using 1 <;> (try funext x) <;>
        simp [dualCertVDerivMonomial, dualCertMonomialEval, Nat.cast_add,
          Nat.cast_one, Nat.succ_eq_add_one] <;> ring

theorem dualCertEval_p_hasDerivAt (P : DualCertPoly) (n q v p : ℝ) :
    HasDerivAt (fun x => dualCertEval P n x q v)
      (dualCertEval (dualCertPDeriv P) n p q v) p := by
  rcases P with ⟨ts⟩
  simp only [dualCertEval, dualCertPDeriv]
  induction ts with
  | nil => simpa using hasDerivAt_const p (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, List.flatMap_cons, List.map_append,
        List.sum_append]
      exact (dualCertPDerivMonomial_hasDerivAt t n q v p).add ih

theorem dualCertEval_q_hasDerivAt (P : DualCertPoly) (n p v q : ℝ) :
    HasDerivAt (fun x => dualCertEval P n p x v)
      (dualCertEval (dualCertQDeriv P) n p q v) q := by
  rcases P with ⟨ts⟩
  simp only [dualCertEval, dualCertQDeriv]
  induction ts with
  | nil => simpa using hasDerivAt_const q (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, List.flatMap_cons, List.map_append,
        List.sum_append]
      exact (dualCertQDerivMonomial_hasDerivAt t n p v q).add ih

theorem dualCertEval_v_hasDerivAt (P : DualCertPoly) (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertEval P n p q x)
      (dualCertEval (dualCertVDeriv P) n p q v) v := by
  rcases P with ⟨ts⟩
  simp only [dualCertEval, dualCertVDeriv]
  induction ts with
  | nil => simpa using hasDerivAt_const v (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, List.flatMap_cons, List.map_append,
        List.sum_append]
      exact (dualCertVDerivMonomial_hasDerivAt t n p q v).add ih

macro "dual_cert_eval_num" name:ident value:num : command =>
  `(@[simp] private theorem $name (n p q v : ℝ) :
      dualCertEval ($value : DualCertPoly) n p q v = ($value : ℝ) := by
    change dualCertEval (dualCertConst ($value : ℤ)) n p q v = _
    simp)

dual_cert_eval_num dualCertEval_zero' 0
dual_cert_eval_num dualCertEval_one' 1
dual_cert_eval_num dualCertEval_two' 2
dual_cert_eval_num dualCertEval_three' 3
dual_cert_eval_num dualCertEval_four' 4
dual_cert_eval_num dualCertEval_five' 5
dual_cert_eval_num dualCertEval_six' 6
dual_cert_eval_num dualCertEval_seven' 7
dual_cert_eval_num dualCertEval_eight' 8
dual_cert_eval_num dualCertEval_nine' 9
dual_cert_eval_num dualCertEval_ten' 10

@[simp] theorem dualCertEval_DPoly (n p q v : ℝ) :
    dualCertEval dualCertDPoly n p q v = dualCertD p q v := by
  simp [dualCertDPoly, dualCertD]

@[simp] theorem dualCertEval_SnumPoly (n p q v : ℝ) :
    dualCertEval dualCertSnumPoly n p q v = dualCertSnum p q v := by
  simp [dualCertSnumPoly, dualCertSnum]

@[simp] theorem dualCertEval_OpPPoly (P : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (dualCertOpPPoly P) n p q v =
      dualCertOpP n (dualCertEval P n p q v)
        (dualCertEval (dualCertPDeriv P) n p q v) p q v := by
  simp [dualCertOpPPoly, dualCertOpP] <;> ring

@[simp] theorem dualCertEval_OpQPoly (P : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (dualCertOpQPoly P) n p q v =
      dualCertOpQ n (dualCertEval P n p q v)
        (dualCertEval (dualCertQDeriv P) n p q v) p q v := by
  simp [dualCertOpQPoly, dualCertOpQ] <;> ring

@[simp] theorem dualCertEval_OpVPoly (P : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (dualCertOpVPoly P) n p q v =
      dualCertOpV n (dualCertEval P n p q v)
        (dualCertEval (dualCertVDeriv P) n p q v) p q v := by
  simp [dualCertOpVPoly, dualCertOpV] <;> ring

end RamanujanChallenge.P25

end
