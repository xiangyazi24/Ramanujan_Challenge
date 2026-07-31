/-
  Axiom audit for the Problem 2.1 development.  Build the library first,
  then:   lake env lean AxiomCheck.lean
-/
import Ramanujan21.Problem21

open RamanujanChallenge.P21

#print axioms shift_a
#print axioms shift_b
#print axioms cfP_neg
#print axioms cfQ_neg
#print axioms cf_neg_convergent
#print axioms challenge_convergent_eq
#print axioms problem21_pcf_value

#check @problem21_pcf_value
