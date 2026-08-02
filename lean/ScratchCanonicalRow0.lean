import RamanujanChallenge.Problem25DualCertificateRow0

namespace RamanujanChallenge.P25

structure CPoly where
  val : DualCertPoly

def CPoly.const (z : ℤ) : CPoly := ⟨dualCertConst z⟩
def CPoly.n : CPoly := ⟨dualCertN⟩
def CPoly.p : CPoly := ⟨dualCertP⟩
def CPoly.q : CPoly := ⟨dualCertQ⟩
def CPoly.v : CPoly := ⟨dualCertV⟩
def CPoly.add (P Q : CPoly) : CPoly := ⟨dualCertNormalize (P.val + Q.val)⟩
def CPoly.neg (P : CPoly) : CPoly := ⟨dualCertNormalize (-P.val)⟩
def CPoly.mul (P Q : CPoly) : CPoly := ⟨dualCertNormalize (P.val * Q.val)⟩
def CPoly.pow (P : CPoly) : ℕ → CPoly
  | 0 => CPoly.const 1
  | k + 1 => CPoly.mul (CPoly.pow P k) P

instance : OfNat CPoly n where ofNat := CPoly.const n
instance : Add CPoly := ⟨CPoly.add⟩
instance : Neg CPoly := ⟨CPoly.neg⟩
instance : Sub CPoly := ⟨fun P Q => P + -Q⟩
instance : Mul CPoly := ⟨CPoly.mul⟩
instance : Pow CPoly ℕ := ⟨CPoly.pow⟩

def cD : CPoly :=
  CPoly.p * CPoly.q * (1 + CPoly.v ^ 2) + 2 * CPoly.v

def cSnum : CPoly :=
  CPoly.p ^ 2 * CPoly.q ^ 2 * (1 - CPoly.p ^ 2) *
    (1 - CPoly.q ^ 2) * CPoly.v ^ 2

def cNext (j : Fin 3) : CPoly :=
  ![cSnum * cD ^ 2,
    cSnum * 2 * (CPoly.n + 3) * CPoly.v * cD,
    cSnum *
      (-(CPoly.n + 3) * CPoly.v * cD +
        2 * (CPoly.n + 3) * (2 * CPoly.n + 7) * CPoly.v ^ 2)] j

def cCur (i : Fin 3) : CPoly :=
  ![cD ^ 4,
    2 * (CPoly.n + 2) * CPoly.v * cD ^ 3,
    (-(CPoly.n + 2) * CPoly.v * cD +
        2 * (CPoly.n + 2) * (2 * CPoly.n + 5) * CPoly.v ^ 2) *
      cD ^ 2] i

def cLambda : CPoly :=
  (CPoly.n + 1) * (CPoly.n + 2) ^ 2 * (CPoly.n + 3) ^ 2 *
    (2 * CPoly.n + 7) ^ 2

def cDelta : CPoly := 4 * (2 * CPoly.n + 3) * (CPoly.n + 2)

def cOpP (P dP : CPoly) : CPoly :=
  cD *
      (CPoly.p * (1 - CPoly.p ^ 2) * dP +
        ((2 * CPoly.n + 7) - (4 * CPoly.n + 9) * CPoly.p ^ 2) * P) -
    (2 * CPoly.n + 7) * CPoly.p * (1 - CPoly.p ^ 2) *
      (CPoly.q * (1 + CPoly.v ^ 2)) * P

def cOpQ (P dQ : CPoly) : CPoly :=
  cD *
      (CPoly.q * (1 - CPoly.q ^ 2) * dQ +
        ((2 * CPoly.n + 6) - (4 * CPoly.n + 10) * CPoly.q ^ 2) * P) -
    (2 * CPoly.n + 7) * CPoly.q * (1 - CPoly.q ^ 2) *
      (CPoly.p * (1 + CPoly.v ^ 2)) * P

def cOpV (P dV : CPoly) : CPoly :=
  cD *
      (CPoly.v * (1 - CPoly.v ^ 2) * dV +
        ((2 * CPoly.n + 4) - (2 * CPoly.n + 6) * CPoly.v ^ 2) * P) -
    (2 * CPoly.n + 7) * CPoly.v * (1 - CPoly.v ^ 2) *
      (2 * CPoly.p * CPoly.q * CPoly.v + 2) * P

def cM00 : CPoly :=
  (2 * CPoly.n + 5) * (CPoly.n + 3) ^ 2 *
    (136 * CPoly.n ^ 4 + 1424 * CPoly.n ^ 3 + 5548 * CPoly.n ^ 2 +
      9551 * CPoly.n + 6141)

def cM01 : CPoly :=
  384 * CPoly.n ^ 6 + 6384 * CPoly.n ^ 5 + 44168 * CPoly.n ^ 4 +
    162698 * CPoly.n ^ 3 + 336377 * CPoly.n ^ 2 + 369933 * CPoly.n + 169011

def cM02 : CPoly :=
  480 * CPoly.n ^ 4 + 4980 * CPoly.n ^ 3 + 19210 * CPoly.n ^ 2 +
    32690 * CPoly.n + 20730

def cPp0 : CPoly := ⟨dualCertPp0Poly⟩
def cPq0 : CPoly := ⟨dualCertPq0Poly⟩
def cPv0 : CPoly := ⟨dualCertPv0Poly⟩
def cDPp0 : CPoly := ⟨dualCertPDeriv dualCertPp0Poly⟩
def cDPq0 : CPoly := ⟨dualCertQDeriv dualCertPq0Poly⟩
def cDPv0 : CPoly := ⟨dualCertVDeriv dualCertPv0Poly⟩

set_option maxHeartbeats 0 in
set_option maxRecDepth 100000 in
example :
    (cDelta * (cM00 * cNext 0 + cM01 * cNext 1 + cM02 * cNext 2 -
      cLambda * cCur 0)).val =
    (cOpP cPp0 cDPp0 + cOpQ cPq0 cDPq0 + cOpV cPv0 cDPv0).val := by
  rfl

end RamanujanChallenge.P25
