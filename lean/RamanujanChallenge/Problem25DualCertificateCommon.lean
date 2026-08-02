import RamanujanChallenge.Problem25
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Mul

noncomputable section

namespace RamanujanChallenge.P25

structure DualCertTerm where
  nCoeffs : List ℤ
  pExp : ℕ
  qExp : ℕ
  vExp : ℕ
  deriving DecidableEq

structure DualCertExp where
  nExp : ℕ
  pExp : ℕ
  qExp : ℕ
  vExp : ℕ
  deriving DecidableEq, Repr

structure DualCertMonomial where
  coeff : ℤ
  exp : DualCertExp
  deriving DecidableEq, Repr

structure DualCertPoly where
  terms : List DualCertMonomial
  deriving DecidableEq, Repr

def dualCertConst (z : ℤ) : DualCertPoly :=
  ⟨[{ coeff := z, exp := ⟨0, 0, 0, 0⟩ }]⟩

def dualCertN : DualCertPoly :=
  ⟨[{ coeff := 1, exp := ⟨1, 0, 0, 0⟩ }]⟩

def dualCertP : DualCertPoly :=
  ⟨[{ coeff := 1, exp := ⟨0, 1, 0, 0⟩ }]⟩

def dualCertQ : DualCertPoly :=
  ⟨[{ coeff := 1, exp := ⟨0, 0, 1, 0⟩ }]⟩

def dualCertV : DualCertPoly :=
  ⟨[{ coeff := 1, exp := ⟨0, 0, 0, 1⟩ }]⟩

def dualCertAdd (P Q : DualCertPoly) : DualCertPoly :=
  ⟨P.terms ++ Q.terms⟩

def dualCertNeg (P : DualCertPoly) : DualCertPoly :=
  ⟨P.terms.map fun t => { t with coeff := -t.coeff }⟩

def dualCertMulMonomial (a b : DualCertMonomial) : DualCertMonomial :=
  { coeff := a.coeff * b.coeff
    exp :=
      { nExp := a.exp.nExp + b.exp.nExp
        pExp := a.exp.pExp + b.exp.pExp
        qExp := a.exp.qExp + b.exp.qExp
        vExp := a.exp.vExp + b.exp.vExp } }

def dualCertMul (P Q : DualCertPoly) : DualCertPoly :=
  ⟨P.terms.flatMap fun a => Q.terms.map (dualCertMulMonomial a)⟩

def dualCertPow (P : DualCertPoly) : ℕ → DualCertPoly
  | 0 => dualCertConst 1
  | k + 1 => dualCertMul (dualCertPow P k) P

instance : OfNat DualCertPoly n where
  ofNat := dualCertConst n

instance : Add DualCertPoly := ⟨dualCertAdd⟩
instance : Neg DualCertPoly := ⟨dualCertNeg⟩
instance : Sub DualCertPoly := ⟨fun P Q => dualCertAdd P (dualCertNeg Q)⟩
instance : Mul DualCertPoly := ⟨dualCertMul⟩
instance : Pow DualCertPoly ℕ := ⟨dualCertPow⟩

def dualCertOfNCoeffs (pExp qExp vExp : ℕ) : ℕ → List ℤ → List DualCertMonomial
  | _, [] => []
  | nExp, c :: cs =>
      { coeff := c, exp := ⟨nExp, pExp, qExp, vExp⟩ } ::
        dualCertOfNCoeffs pExp qExp vExp (nExp + 1) cs

def dualCertOfTerm (t : DualCertTerm) : List DualCertMonomial :=
  dualCertOfNCoeffs t.pExp t.qExp t.vExp 0 t.nCoeffs

def dualCertOfTerms (ts : List DualCertTerm) : DualCertPoly :=
  ⟨ts.flatMap dualCertOfTerm⟩

def dualCertPDerivMonomial (t : DualCertMonomial) : List DualCertMonomial :=
  if h : t.exp.pExp = 0 then [] else
    [{ coeff := t.coeff * t.exp.pExp,
       exp := { t.exp with pExp := t.exp.pExp - 1 } }]

def dualCertQDerivMonomial (t : DualCertMonomial) : List DualCertMonomial :=
  if h : t.exp.qExp = 0 then [] else
    [{ coeff := t.coeff * t.exp.qExp,
       exp := { t.exp with qExp := t.exp.qExp - 1 } }]

def dualCertVDerivMonomial (t : DualCertMonomial) : List DualCertMonomial :=
  if h : t.exp.vExp = 0 then [] else
    [{ coeff := t.coeff * t.exp.vExp,
       exp := { t.exp with vExp := t.exp.vExp - 1 } }]

def dualCertPDeriv (P : DualCertPoly) : DualCertPoly :=
  ⟨P.terms.flatMap dualCertPDerivMonomial⟩

def dualCertQDeriv (P : DualCertPoly) : DualCertPoly :=
  ⟨P.terms.flatMap dualCertQDerivMonomial⟩

def dualCertVDeriv (P : DualCertPoly) : DualCertPoly :=
  ⟨P.terms.flatMap dualCertVDerivMonomial⟩

def dualCertExpLT (a b : DualCertExp) : Bool :=
  if a.nExp < b.nExp then true else if b.nExp < a.nExp then false else
  if a.pExp < b.pExp then true else if b.pExp < a.pExp then false else
  if a.qExp < b.qExp then true else if b.qExp < a.qExp then false else
  a.vExp < b.vExp

def dualCertExpCompare (a b : DualCertExp) : Ordering :=
  if a = b then .eq else if dualCertExpLT a b then .lt else .gt

def dualCertInsert (t : DualCertMonomial) :
    List DualCertMonomial → List DualCertMonomial
  | [] => if t.coeff = 0 then [] else [t]
  | u :: us =>
      match dualCertExpCompare t.exp u.exp with
      | .lt => if t.coeff = 0 then u :: us else t :: u :: us
      | .eq =>
          let c := t.coeff + u.coeff
          if c = 0 then us else { coeff := c, exp := t.exp } :: us
      | .gt => u :: dualCertInsert t us

def dualCertNormalize (P : DualCertPoly) : DualCertPoly :=
  ⟨P.terms.foldr dualCertInsert []⟩

def dualCertDPoly : DualCertPoly :=
  dualCertP * dualCertQ * (1 + dualCertV ^ 2) + 2 * dualCertV

def dualCertSnumPoly : DualCertPoly :=
  dualCertP ^ 2 * dualCertQ ^ 2 * (1 - dualCertP ^ 2) *
    (1 - dualCertQ ^ 2) * dualCertV ^ 2

def dualCertNextNumPoly (j : Fin 3) : DualCertPoly :=
  ![dualCertSnumPoly * dualCertDPoly ^ 2,
    dualCertSnumPoly * 2 * (dualCertN + 3) * dualCertV * dualCertDPoly,
    dualCertSnumPoly *
      (-(dualCertN + 3) * dualCertV * dualCertDPoly +
        2 * (dualCertN + 3) * (2 * dualCertN + 7) * dualCertV ^ 2)] j

def dualCertCurNumPoly (i : Fin 3) : DualCertPoly :=
  ![dualCertDPoly ^ 4,
    2 * (dualCertN + 2) * dualCertV * dualCertDPoly ^ 3,
    (-(dualCertN + 2) * dualCertV * dualCertDPoly +
        2 * (dualCertN + 2) * (2 * dualCertN + 5) * dualCertV ^ 2) *
      dualCertDPoly ^ 2] i

def dualCertLambdaPoly : DualCertPoly :=
  (dualCertN + 1) * (dualCertN + 2) ^ 2 * (dualCertN + 3) ^ 2 *
    (2 * dualCertN + 7) ^ 2

def dualCertDeltaPoly : DualCertPoly :=
  4 * (2 * dualCertN + 3) * (dualCertN + 2)

def dualCertOpPPoly (P : DualCertPoly) : DualCertPoly :=
  dualCertDPoly *
      (dualCertP * (1 - dualCertP ^ 2) * dualCertPDeriv P +
        ((2 * dualCertN + 7) - (4 * dualCertN + 9) * dualCertP ^ 2) * P) -
    (2 * dualCertN + 7) * dualCertP * (1 - dualCertP ^ 2) *
      (dualCertQ * (1 + dualCertV ^ 2)) * P

def dualCertOpQPoly (P : DualCertPoly) : DualCertPoly :=
  dualCertDPoly *
      (dualCertQ * (1 - dualCertQ ^ 2) * dualCertQDeriv P +
        ((2 * dualCertN + 6) - (4 * dualCertN + 10) * dualCertQ ^ 2) * P) -
    (2 * dualCertN + 7) * dualCertQ * (1 - dualCertQ ^ 2) *
      (dualCertP * (1 + dualCertV ^ 2)) * P

def dualCertOpVPoly (P : DualCertPoly) : DualCertPoly :=
  dualCertDPoly *
      (dualCertV * (1 - dualCertV ^ 2) * dualCertVDeriv P +
        ((2 * dualCertN + 4) - (2 * dualCertN + 6) * dualCertV ^ 2) * P) -
    (2 * dualCertN + 7) * dualCertV * (1 - dualCertV ^ 2) *
      (2 * dualCertP * dualCertQ * dualCertV + 2) * P

def dualCertMonomialEval (n p q v : ℝ) (t : DualCertMonomial) : ℝ :=
  (t.coeff : ℝ) * n ^ t.exp.nExp * p ^ t.exp.pExp *
    q ^ t.exp.qExp * v ^ t.exp.vExp

def dualCertEval (P : DualCertPoly) (n p q v : ℝ) : ℝ :=
  (P.terms.map (dualCertMonomialEval n p q v)).sum

private theorem dualCertExpCompare_eq_iff (a b : DualCertExp) :
    dualCertExpCompare a b = .eq ↔ a = b := by
  constructor
  · intro h
    by_contra hab
    simp only [dualCertExpCompare, if_neg hab] at h
    split at h <;> simp_all
  · intro h
    subst b
    simp [dualCertExpCompare]

private theorem dualCertMonomialEval_add_coeff
    (n p q v : ℝ) (a b : DualCertMonomial) (h : a.exp = b.exp) :
    dualCertMonomialEval n p q v { coeff := a.coeff + b.coeff, exp := a.exp } =
      dualCertMonomialEval n p q v a + dualCertMonomialEval n p q v b := by
  rcases a with ⟨ac, ae⟩
  rcases b with ⟨bc, be⟩
  change ae = be at h
  subst be
  unfold dualCertMonomialEval
  push_cast
  ring

private theorem dualCertEval_insert (n p q v : ℝ) (t : DualCertMonomial) :
    ∀ ts : List DualCertMonomial,
      ((dualCertInsert t ts).map (dualCertMonomialEval n p q v)).sum =
        dualCertMonomialEval n p q v t +
          (ts.map (dualCertMonomialEval n p q v)).sum := by
  intro ts
  induction ts with
  | nil =>
      simp only [dualCertInsert]
      split_ifs with h
      · simp [dualCertMonomialEval, h]
      · simp
  | cons u us ih =>
      simp only [dualCertInsert]
      cases hcmp : dualCertExpCompare t.exp u.exp with
      | lt =>
          simp only [hcmp]
          split_ifs with h
          · simp [dualCertMonomialEval, h]
          · simp
      | eq =>
          simp only [hcmp]
          have hexp : t.exp = u.exp :=
            (dualCertExpCompare_eq_iff t.exp u.exp).mp hcmp
          split_ifs with hsum
          · have heval :
                dualCertMonomialEval n p q v t +
                    dualCertMonomialEval n p q v u = 0 := by
              rw [← dualCertMonomialEval_add_coeff n p q v t u hexp]
              simp [dualCertMonomialEval, hsum]
            simp only [List.map_cons, List.sum_cons]
            linarith
          · simp only [List.map_cons, List.sum_cons]
            rw [dualCertMonomialEval_add_coeff n p q v t u hexp]
            ring
      | gt =>
          simp only [hcmp, List.map_cons, List.sum_cons, ih]
          ring

theorem dualCertEval_normalize (P : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (dualCertNormalize P) n p q v = dualCertEval P n p q v := by
  rcases P with ⟨ts⟩
  unfold dualCertNormalize dualCertEval
  induction ts with
  | nil => simp
  | cons t ts ih =>
      simp only [List.foldr_cons, List.map_cons, List.sum_cons]
      rw [dualCertEval_insert, ih]

@[simp] theorem dualCertEval_const (z : ℤ) (n p q v : ℝ) :
    dualCertEval (dualCertConst z) n p q v = z := by
  simp [dualCertEval, dualCertConst, dualCertMonomialEval]

@[simp] theorem dualCertEval_N (n p q v : ℝ) :
    dualCertEval dualCertN n p q v = n := by
  simp [dualCertEval, dualCertN, dualCertMonomialEval]

@[simp] theorem dualCertEval_P (n p q v : ℝ) :
    dualCertEval dualCertP n p q v = p := by
  simp [dualCertEval, dualCertP, dualCertMonomialEval]

@[simp] theorem dualCertEval_Q (n p q v : ℝ) :
    dualCertEval dualCertQ n p q v = q := by
  simp [dualCertEval, dualCertQ, dualCertMonomialEval]

@[simp] theorem dualCertEval_V (n p q v : ℝ) :
    dualCertEval dualCertV n p q v = v := by
  simp [dualCertEval, dualCertV, dualCertMonomialEval]

@[simp] theorem dualCertEval_ofNat (k : ℕ) (n p q v : ℝ) :
    dualCertEval (OfNat.ofNat k : DualCertPoly) n p q v = k := by
  change dualCertEval (dualCertConst k) n p q v = k
  simp

@[simp] theorem dualCertEval_add (P Q : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (P + Q) n p q v =
      dualCertEval P n p q v + dualCertEval Q n p q v := by
  simp [HAdd.hAdd, Add.add, dualCertAdd, dualCertEval]

@[simp] theorem dualCertEval_neg (P : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (-P) n p q v = -dualCertEval P n p q v := by
  rcases P with ⟨ts⟩
  change dualCertEval (dualCertNeg ⟨ts⟩) n p q v = -dualCertEval ⟨ts⟩ n p q v
  unfold dualCertNeg dualCertEval
  induction ts with
  | nil => simp
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, ih]
      unfold dualCertMonomialEval
      push_cast
      ring

@[simp] theorem dualCertEval_sub (P Q : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (P - Q) n p q v =
      dualCertEval P n p q v - dualCertEval Q n p q v := by
  change dualCertEval (dualCertAdd P (dualCertNeg Q)) n p q v = _
  rw [show dualCertAdd P (dualCertNeg Q) = P + (-Q) by rfl]
  rw [dualCertEval_add, dualCertEval_neg]
  ring

private theorem dualCertMonomialEval_mul (n p q v : ℝ)
    (a b : DualCertMonomial) :
    dualCertMonomialEval n p q v (dualCertMulMonomial a b) =
      dualCertMonomialEval n p q v a * dualCertMonomialEval n p q v b := by
  unfold dualCertMonomialEval dualCertMulMonomial
  push_cast
  simp only [pow_add]
  ring

@[simp] theorem dualCertEval_mul (P Q : DualCertPoly) (n p q v : ℝ) :
    dualCertEval (P * Q) n p q v =
      dualCertEval P n p q v * dualCertEval Q n p q v := by
  rcases P with ⟨ps⟩
  rcases Q with ⟨qs⟩
  change dualCertEval (dualCertMul ⟨ps⟩ ⟨qs⟩) n p q v = _
  unfold dualCertMul dualCertEval
  induction ps with
  | nil => simp
  | cons a ps ih =>
      simp only [List.flatMap_cons, List.map_append, List.sum_append,
        List.map_map, ih]
      have hmap :
          List.map (dualCertMonomialEval n p q v ∘ dualCertMulMonomial a) qs =
            List.map (fun b => dualCertMonomialEval n p q v a *
              dualCertMonomialEval n p q v b) qs := by
        apply List.map_congr_left
        intro b hb
        exact dualCertMonomialEval_mul n p q v a b
      rw [hmap]
      rw [List.sum_map_mul_left]
      simp only [List.map_cons, List.sum_cons]
      ring

@[simp] theorem dualCertEval_pow (P : DualCertPoly) (k : ℕ) (n p q v : ℝ) :
    dualCertEval (P ^ k) n p q v = dualCertEval P n p q v ^ k := by
  change dualCertEval (dualCertPow P k) n p q v = _
  induction k with
  | zero => simp [dualCertPow]
  | succ k ih =>
      rw [dualCertPow]
      rw [show dualCertMul (dualCertPow P k) P = dualCertPow P k * P by rfl]
      rw [dualCertEval_mul, ih, pow_succ]

private theorem dualCertMonomial_hasDerivAt_p (t : DualCertMonomial)
    (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertMonomialEval n x q v t)
      ((dualCertPDerivMonomial t).map
        (dualCertMonomialEval n p q v)).sum p := by
  by_cases h : t.exp.pExp = 0
  · simp only [dualCertPDerivMonomial, dif_pos h, List.map_nil, List.sum_nil]
    convert hasDerivAt_const p (dualCertMonomialEval n p q v t) using 1
    funext x
    unfold dualCertMonomialEval
    simp [h]
  · simp only [dualCertPDerivMonomial, dif_neg h, List.map_cons, List.map_nil,
      List.sum_cons, List.sum_nil, add_zero]
    let C : ℝ := (t.coeff : ℝ) * n ^ t.exp.nExp *
      q ^ t.exp.qExp * v ^ t.exp.vExp
    have hd := (hasDerivAt_id p).pow t.exp.pExp
    convert hd.const_mul C using 1
    · funext x
      unfold dualCertMonomialEval C
      simp only [Pi.pow_apply, id_eq]
      ring
    · unfold dualCertMonomialEval C
      simp only [Pi.pow_apply, id_eq]
      push_cast
      ring

theorem dualCertEval_hasDerivAt_p (P : DualCertPoly) (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertEval P n x q v)
      (dualCertEval (dualCertPDeriv P) n p q v) p := by
  rcases P with ⟨ts⟩
  unfold dualCertEval dualCertPDeriv
  induction ts with
  | nil => simpa using hasDerivAt_const p (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, List.flatMap_cons,
        List.map_append, List.sum_append]
      exact (dualCertMonomial_hasDerivAt_p t n p q v).add ih

private theorem dualCertMonomial_hasDerivAt_q (t : DualCertMonomial)
    (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertMonomialEval n p x v t)
      ((dualCertQDerivMonomial t).map
        (dualCertMonomialEval n p q v)).sum q := by
  by_cases h : t.exp.qExp = 0
  · simp only [dualCertQDerivMonomial, dif_pos h, List.map_nil, List.sum_nil]
    convert hasDerivAt_const q (dualCertMonomialEval n p q v t) using 1
    funext x
    unfold dualCertMonomialEval
    simp [h]
  · simp only [dualCertQDerivMonomial, dif_neg h, List.map_cons, List.map_nil,
      List.sum_cons, List.sum_nil, add_zero]
    let C : ℝ := (t.coeff : ℝ) * n ^ t.exp.nExp *
      p ^ t.exp.pExp * v ^ t.exp.vExp
    have hd := (hasDerivAt_id q).pow t.exp.qExp
    convert hd.const_mul C using 1
    · funext x
      unfold dualCertMonomialEval C
      simp only [Pi.pow_apply, id_eq]
      ring
    · unfold dualCertMonomialEval C
      simp only [Pi.pow_apply, id_eq]
      push_cast
      ring

theorem dualCertEval_hasDerivAt_q (P : DualCertPoly) (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertEval P n p x v)
      (dualCertEval (dualCertQDeriv P) n p q v) q := by
  rcases P with ⟨ts⟩
  unfold dualCertEval dualCertQDeriv
  induction ts with
  | nil => simpa using hasDerivAt_const q (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, List.flatMap_cons,
        List.map_append, List.sum_append]
      exact (dualCertMonomial_hasDerivAt_q t n p q v).add ih

private theorem dualCertMonomial_hasDerivAt_v (t : DualCertMonomial)
    (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertMonomialEval n p q x t)
      ((dualCertVDerivMonomial t).map
        (dualCertMonomialEval n p q v)).sum v := by
  by_cases h : t.exp.vExp = 0
  · simp only [dualCertVDerivMonomial, dif_pos h, List.map_nil, List.sum_nil]
    convert hasDerivAt_const v (dualCertMonomialEval n p q v t) using 1
    funext x
    unfold dualCertMonomialEval
    simp [h]
  · simp only [dualCertVDerivMonomial, dif_neg h, List.map_cons, List.map_nil,
      List.sum_cons, List.sum_nil, add_zero]
    let C : ℝ := (t.coeff : ℝ) * n ^ t.exp.nExp *
      p ^ t.exp.pExp * q ^ t.exp.qExp
    have hd := (hasDerivAt_id v).pow t.exp.vExp
    convert hd.const_mul C using 1
    unfold dualCertMonomialEval C
    simp only [Pi.pow_apply, id_eq]
    push_cast
    ring

theorem dualCertEval_hasDerivAt_v (P : DualCertPoly) (n p q v : ℝ) :
    HasDerivAt (fun x => dualCertEval P n p q x)
      (dualCertEval (dualCertVDeriv P) n p q v) v := by
  rcases P with ⟨ts⟩
  unfold dualCertEval dualCertVDeriv
  induction ts with
  | nil => simpa using hasDerivAt_const v (0 : ℝ)
  | cons t ts ih =>
      simp only [List.map_cons, List.sum_cons, List.flatMap_cons,
        List.map_append, List.sum_append]
      exact (dualCertMonomial_hasDerivAt_v t n p q v).add ih

def dualCertD (p q v : ℝ) : ℝ :=
  p * q * (1 + v ^ 2) + 2 * v

def dualCertSnum (p q v : ℝ) : ℝ :=
  p ^ 2 * q ^ 2 * (1 - p ^ 2) * (1 - q ^ 2) * v ^ 2

def dualCertNextNum (n : ℝ) (j : Fin 3) (p q v : ℝ) : ℝ :=
  ![dualCertSnum p q v * dualCertD p q v ^ 2,
    dualCertSnum p q v * 2 * (n + 3) * v * dualCertD p q v,
    dualCertSnum p q v *
      (-(n + 3) * v * dualCertD p q v +
        2 * (n + 3) * (2 * n + 7) * v ^ 2)] j

def dualCertCurNum (n : ℝ) (i : Fin 3) (p q v : ℝ) : ℝ :=
  ![dualCertD p q v ^ 4,
    2 * (n + 2) * v * dualCertD p q v ^ 3,
    (-(n + 2) * v * dualCertD p q v +
        2 * (n + 2) * (2 * n + 5) * v ^ 2) * dualCertD p q v ^ 2] i

def dualCertLambda (n : ℝ) : ℝ :=
  (n + 1) * (n + 2) ^ 2 * (n + 3) ^ 2 * (2 * n + 7) ^ 2

def dualCertOpP (n P dP p q v : ℝ) : ℝ :=
  dualCertD p q v *
      (p * (1 - p ^ 2) * dP + ((2 * n + 7) - (4 * n + 9) * p ^ 2) * P) -
    (2 * n + 7) * p * (1 - p ^ 2) * (q * (1 + v ^ 2)) * P

def dualCertOpQ (n P dQ p q v : ℝ) : ℝ :=
  dualCertD p q v *
      (q * (1 - q ^ 2) * dQ + ((2 * n + 6) - (4 * n + 10) * q ^ 2) * P) -
    (2 * n + 7) * q * (1 - q ^ 2) * (p * (1 + v ^ 2)) * P

def dualCertOpV (n P dV p q v : ℝ) : ℝ :=
  dualCertD p q v *
      (v * (1 - v ^ 2) * dV + ((2 * n + 4) - (2 * n + 6) * v ^ 2) * P) -
    (2 * n + 7) * v * (1 - v ^ 2) * (2 * p * q * v + 2) * P

end RamanujanChallenge.P25

end
