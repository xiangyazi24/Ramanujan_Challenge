# Problem 3.1: Knot Integral for π² (CONJECTURE)

## Statement

Let A_{7_2}(M,L) be the A-polynomial for the prime knot 7_2.
(Full polynomial given in paper — degree 22 in M, degree 5 in L.)

Let α ≈ 0.349269 be the real root of A_{7_2}(α, α^{1/2}) = 0 closest to 0.349269.
Let β ≈ 0.406813 be the real root of A_{7_2}(β, β) = 0 closest to 0.406813.
Let y = y(x) be the curve satisfying A(x, y(x)) = 0 for α ≤ x ≤ β,
positive and decreasing with y'(x) ≤ -2.

Prove:
```
4π²/85 = ∫_α^β (log x · dy/y - log y · dx/x)
```

## Context
Conjectural closed form of a Godbillon-Vey type Knot invariant (Khoi 2008).

## Approach Ideas
- This is from knot theory — quite different from the other problems
- The A-polynomial of 7_2 is a classical object in knot theory
- The integral is a "Mahler measure" type object
- Connection to Bloch-Wigner dilogarithm and hyperbolic volume
- Check: Khoi (2008) "On the Integral of log x dy/y - log y dx/x"
- Numerical verification first: compute the integral to high precision
- This may require algebraic geometry tools more than number theory
- Probably the hardest problem — leave for later
