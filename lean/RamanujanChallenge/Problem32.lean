/-
  Ramanujan Challenge Problem 3.2: The Apéry GCD Conjecture.

  Prove: gcd(d_n · a_n, d_n · b_n) = e^{o(n)} where
  (a_n), (b_n) are the Apéry sequences for ζ(3) and d_n = lcm(1,...,n)³.

  PROVED UNCONDITIONALLY:
  - Polylog exceptional set: #{n ≤ N : log G_n > εn} = O_ε((log N)²)
  - Upper Banach density zero
  - Finite harmonic weight

  Architecture:
    AperyDef.lean      — Apéry recurrence, b_n closed form, d_n, G_n
    Wronskian.lean     — W_n = 6/n³, valuation bound, no consecutive zeros
    GapPolynomial.lean — canonical gap numerators and their degree bound
    Main.lean          — Main theorem statements
-/
import RamanujanChallenge.Problem32.AperyDef
import RamanujanChallenge.Problem32.Wronskian
import RamanujanChallenge.Problem32.GapPolynomial
import RamanujanChallenge.Problem32.Main
