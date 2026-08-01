/-
  Axiom audit for the Problem 2.3 development.  Build the library first,
  then:   lake env lean AxiomCheck.lean
-/
import Ramanujan23.Problem23

open RamanujanChallenge.P23

#print axioms tensor_rec
#print axioms factorial_isDerRec
#print axioms challengeQ_rec
#print axioms challengeP_rec
#print axioms C0_ne_zero
#print axioms eq_of_satisfiesRec
#print axioms ratio_split
#print axioms lambertB_div_lambertA_tendsto_pi_div_four
#print axioms factorial_div_derang_tendsto_exp_one
#print axioms problem23_pi_add_e

-- The challenge's initial values, computed from the definitions:
#eval (challengeQ 0, challengeQ 1, challengeQ 2, challengeQ 3)   -- (1, 0, 4, 48)
#eval (challengeP 0, challengeP 1, challengeP 2, challengeP 3)   -- (1, 1, 20, 296)
-- and the first values the recurrence produces:
#eval (challengeQ 4, challengeQ 5)                               -- (1836, 97680)
#eval (challengeP 4, challengeP 5)                               -- (10656, 573344)
