import RamanujanChallenge.Problem25DualCertificateRow0

namespace RamanujanChallenge.P25

def spCoeffNormalize : List ℤ → List ℤ
  | [] => []
  | c :: cs =>
      match spCoeffNormalize cs with
      | [] => if c = 0 then [] else [c]
      | ds => c :: ds

def spCoeffAddRaw : List ℤ → List ℤ → List ℤ
  | [], bs => bs
  | as, [] => as
  | a :: as, b :: bs => (a + b) :: spCoeffAddRaw as bs

def spCoeffAdd (as bs : List ℤ) : List ℤ :=
  spCoeffNormalize (spCoeffAddRaw as bs)

def spCoeffNeg (as : List ℤ) : List ℤ :=
  spCoeffNormalize (as.map (-·))

def spCoeffScale (a : ℤ) (bs : List ℤ) : List ℤ :=
  bs.map (a * ·)

def spCoeffMulRaw : List ℤ → List ℤ → List ℤ
  | [], _ => []
  | a :: as, bs =>
      spCoeffAddRaw (spCoeffScale a bs) (0 :: spCoeffMulRaw as bs)

def spCoeffMul (as bs : List ℤ) : List ℤ :=
  spCoeffNormalize (spCoeffMulRaw as bs)

structure SparseTerm where
  nCoeffs : List ℤ
  pExp : ℕ
  qExp : ℕ
  vExp : ℕ
  deriving DecidableEq, Repr

structure SparsePoly where
  terms : List SparseTerm
  deriving DecidableEq, Repr

def sparseTermOfDual (t : DualCertTerm) : SparseTerm :=
  ⟨spCoeffNormalize t.nCoeffs, t.pExp, t.qExp, t.vExp⟩

def sparseExp (t : SparseTerm) : DualCertExp :=
  ⟨0, t.pExp, t.qExp, t.vExp⟩

def sparseInsert (t : SparseTerm) : List SparseTerm → List SparseTerm
  | [] => if t.nCoeffs = [] then [] else [t]
  | u :: us =>
      match dualCertExpCompare (sparseExp t) (sparseExp u) with
      | .lt => if t.nCoeffs = [] then u :: us else t :: u :: us
      | .eq =>
          let cs := spCoeffAdd t.nCoeffs u.nCoeffs
          if cs = [] then us else { t with nCoeffs := cs } :: us
      | .gt => u :: sparseInsert t us

def sparseNormalize (P : SparsePoly) : SparsePoly :=
  ⟨P.terms.foldr sparseInsert []⟩

def SparsePoly.const (z : ℤ) : SparsePoly :=
  sparseNormalize ⟨[⟨[z], 0, 0, 0⟩]⟩

def SparsePoly.n : SparsePoly := ⟨[⟨[0, 1], 0, 0, 0⟩]⟩
def SparsePoly.p : SparsePoly := ⟨[⟨[1], 1, 0, 0⟩]⟩
def SparsePoly.q : SparsePoly := ⟨[⟨[1], 0, 1, 0⟩]⟩
def SparsePoly.v : SparsePoly := ⟨[⟨[1], 0, 0, 1⟩]⟩

def SparsePoly.add (P Q : SparsePoly) : SparsePoly :=
  sparseNormalize ⟨P.terms ++ Q.terms⟩

def SparsePoly.neg (P : SparsePoly) : SparsePoly :=
  sparseNormalize ⟨P.terms.map fun t => { t with nCoeffs := spCoeffNeg t.nCoeffs }⟩

def sparseMulTerm (a b : SparseTerm) : SparseTerm :=
  ⟨spCoeffMul a.nCoeffs b.nCoeffs,
    a.pExp + b.pExp, a.qExp + b.qExp, a.vExp + b.vExp⟩

def SparsePoly.mul (P Q : SparsePoly) : SparsePoly :=
  sparseNormalize ⟨P.terms.flatMap fun a => Q.terms.map (sparseMulTerm a)⟩

def SparsePoly.pow (P : SparsePoly) : ℕ → SparsePoly
  | 0 => SparsePoly.const 1
  | k + 1 => SparsePoly.mul (SparsePoly.pow P k) P

instance : OfNat SparsePoly n where ofNat := SparsePoly.const n
instance : Add SparsePoly := ⟨SparsePoly.add⟩
instance : Neg SparsePoly := ⟨SparsePoly.neg⟩
instance : Sub SparsePoly := ⟨fun P Q => P + -Q⟩
instance : Mul SparsePoly := ⟨SparsePoly.mul⟩
instance : Pow SparsePoly ℕ := ⟨SparsePoly.pow⟩

@[simp] theorem sparse_ofNat_eq (n : ℕ) :
    (OfNat.ofNat n : SparsePoly) = SparsePoly.const n := rfl

@[simp] theorem sparse_add_eq (P Q : SparsePoly) :
    P + Q = sparseNormalize ⟨P.terms ++ Q.terms⟩ := rfl

@[simp] theorem sparse_neg_eq (P : SparsePoly) :
    -P = sparseNormalize
      ⟨P.terms.map fun t => { t with nCoeffs := spCoeffNeg t.nCoeffs }⟩ := rfl

@[simp] theorem sparse_sub_eq (P Q : SparsePoly) : P - Q = P + -Q := rfl

@[simp] theorem sparse_mul_eq (P Q : SparsePoly) :
    P * Q = sparseNormalize
      ⟨P.terms.flatMap fun a => Q.terms.map (sparseMulTerm a)⟩ := rfl

@[simp] theorem sparse_pow_zero (P : SparsePoly) :
    P ^ (0 : ℕ) = SparsePoly.const 1 := rfl

@[simp] theorem sparse_pow_succ (P : SparsePoly) (k : ℕ) :
    P ^ (k + 1) = (P ^ k) * P := rfl

def SparsePoly.ofDualTerms (ts : List DualCertTerm) : SparsePoly :=
  sparseNormalize ⟨ts.map sparseTermOfDual⟩

def spD : SparsePoly :=
  SparsePoly.p * SparsePoly.q * (1 + SparsePoly.v ^ 2) + 2 * SparsePoly.v

def spSnum : SparsePoly :=
  SparsePoly.p ^ 2 * SparsePoly.q ^ 2 * (1 - SparsePoly.p ^ 2) *
    (1 - SparsePoly.q ^ 2) * SparsePoly.v ^ 2

def spNext (j : Fin 3) : SparsePoly :=
  ![spSnum * spD ^ 2,
    spSnum * 2 * (SparsePoly.n + 3) * SparsePoly.v * spD,
    spSnum *
      (-(SparsePoly.n + 3) * SparsePoly.v * spD +
        2 * (SparsePoly.n + 3) * (2 * SparsePoly.n + 7) * SparsePoly.v ^ 2)] j

def spCur (i : Fin 3) : SparsePoly :=
  ![spD ^ 4,
    2 * (SparsePoly.n + 2) * SparsePoly.v * spD ^ 3,
    (-(SparsePoly.n + 2) * SparsePoly.v * spD +
        2 * (SparsePoly.n + 2) * (2 * SparsePoly.n + 5) * SparsePoly.v ^ 2) *
      spD ^ 2] i

def spLambda : SparsePoly :=
  (SparsePoly.n + 1) * (SparsePoly.n + 2) ^ 2 * (SparsePoly.n + 3) ^ 2 *
    (2 * SparsePoly.n + 7) ^ 2

def spDelta : SparsePoly := 4 * (2 * SparsePoly.n + 3) * (SparsePoly.n + 2)

def spDerivP (P : SparsePoly) : SparsePoly :=
  sparseNormalize ⟨P.terms.map fun t =>
    if t.pExp = 0 then { t with nCoeffs := [] }
    else { t with nCoeffs := spCoeffScale t.pExp t.nCoeffs, pExp := t.pExp - 1 }⟩

def spDerivQ (P : SparsePoly) : SparsePoly :=
  sparseNormalize ⟨P.terms.map fun t =>
    if t.qExp = 0 then { t with nCoeffs := [] }
    else { t with nCoeffs := spCoeffScale t.qExp t.nCoeffs, qExp := t.qExp - 1 }⟩

def spDerivV (P : SparsePoly) : SparsePoly :=
  sparseNormalize ⟨P.terms.map fun t =>
    if t.vExp = 0 then { t with nCoeffs := [] }
    else { t with nCoeffs := spCoeffScale t.vExp t.nCoeffs, vExp := t.vExp - 1 }⟩

def spOpP (P : SparsePoly) : SparsePoly :=
  spD *
      (SparsePoly.p * (1 - SparsePoly.p ^ 2) * spDerivP P +
        ((2 * SparsePoly.n + 7) - (4 * SparsePoly.n + 9) * SparsePoly.p ^ 2) * P) -
    (2 * SparsePoly.n + 7) * SparsePoly.p * (1 - SparsePoly.p ^ 2) *
      (SparsePoly.q * (1 + SparsePoly.v ^ 2)) * P

def spOpQ (P : SparsePoly) : SparsePoly :=
  spD *
      (SparsePoly.q * (1 - SparsePoly.q ^ 2) * spDerivQ P +
        ((2 * SparsePoly.n + 6) - (4 * SparsePoly.n + 10) * SparsePoly.q ^ 2) * P) -
    (2 * SparsePoly.n + 7) * SparsePoly.q * (1 - SparsePoly.q ^ 2) *
      (SparsePoly.p * (1 + SparsePoly.v ^ 2)) * P

def spOpV (P : SparsePoly) : SparsePoly :=
  spD *
      (SparsePoly.v * (1 - SparsePoly.v ^ 2) * spDerivV P +
        ((2 * SparsePoly.n + 4) - (2 * SparsePoly.n + 6) * SparsePoly.v ^ 2) * P) -
    (2 * SparsePoly.n + 7) * SparsePoly.v * (1 - SparsePoly.v ^ 2) *
      (2 * SparsePoly.p * SparsePoly.q * SparsePoly.v + 2) * P

def spM00 : SparsePoly :=
  (2 * SparsePoly.n + 5) * (SparsePoly.n + 3) ^ 2 *
    (136 * SparsePoly.n ^ 4 + 1424 * SparsePoly.n ^ 3 + 5548 * SparsePoly.n ^ 2 +
      9551 * SparsePoly.n + 6141)

def spM01 : SparsePoly :=
  384 * SparsePoly.n ^ 6 + 6384 * SparsePoly.n ^ 5 + 44168 * SparsePoly.n ^ 4 +
    162698 * SparsePoly.n ^ 3 + 336377 * SparsePoly.n ^ 2 + 369933 * SparsePoly.n +
      169011

def spM02 : SparsePoly :=
  480 * SparsePoly.n ^ 4 + 4980 * SparsePoly.n ^ 3 + 19210 * SparsePoly.n ^ 2 +
    32690 * SparsePoly.n + 20730

def spPp0 : SparsePoly := SparsePoly.ofDualTerms dualCertPp0Terms
def spPq0 : SparsePoly := SparsePoly.ofDualTerms dualCertPq0Terms
def spPv0 : SparsePoly := SparsePoly.ofDualTerms dualCertPv0Terms

set_option maxHeartbeats 0 in
set_option maxRecDepth 1000000 in
theorem dualCert_sparse_identity_zero :
    spDelta * (spM00 * spNext 0 + spM01 * spNext 1 + spM02 * spNext 2 -
      spLambda * spCur 0) =
    spOpP spPp0 + spOpQ spPq0 + spOpV spPv0 := by
  simp (config := { maxSteps := 10000000 })
    [spDelta, spM00, spM01, spM02, spNext, spCur, spLambda,
    spOpP, spOpQ, spOpV, spD, spSnum, spDerivP, spDerivQ, spDerivV,
    spPp0, spPq0, spPv0, SparsePoly.ofDualTerms,
    dualCertPp0Terms, dualCertPq0Terms, dualCertPv0Terms,
    SparsePoly.const, SparsePoly.n, SparsePoly.p, SparsePoly.q, SparsePoly.v,
    sparseMulTerm, sparseNormalize, sparseInsert, sparseExp, sparseTermOfDual,
    spCoeffMul, spCoeffMulRaw, spCoeffScale, spCoeffNeg,
    spCoeffAdd, spCoeffAddRaw, spCoeffNormalize, dualCertExpCompare,
    dualCertExpLT]

end RamanujanChallenge.P25
