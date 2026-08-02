import RamanujanChallenge.Problem25DualCertificateCommon

#check RamanujanChallenge.P25.instOfNatDualCertPoly
#check RamanujanChallenge.P25.dualCertEval_ofNat
#print RamanujanChallenge.P25.dualCertEval_ofNat

open RamanujanChallenge.P25
example (n p q v : ℝ) : dualCertEval (4 : DualCertPoly) n p q v = 4 := by simp
