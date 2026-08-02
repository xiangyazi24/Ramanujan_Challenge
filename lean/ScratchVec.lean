import RamanujanChallenge.Problem25DualCertificateCommon

namespace RamanujanChallenge.P25

example (j : Fin 3) (a b c : ℝ) : ![a,b,c] j = ![a,b,c] j := by
  fin_cases j
  all_goals simp

example (j : Fin 3) (a b c : ℝ) : ![a,b,c] j = a := by
  fin_cases j
  · simp
  · sorry
  · sorry

end RamanujanChallenge.P25
