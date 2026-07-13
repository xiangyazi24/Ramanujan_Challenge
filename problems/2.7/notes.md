# Problem 2.7: Efficient Four-Term Recurrence for ζ(2) + ζ(3)

## Statement

Define polynomials A_n, B_n, C_n, D_n (degree ~9 each) and recurrence:
```
u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2},  n ≥ 2
```

Two solutions p_n, q_n with large integer initial conditions.

Prove: lim p_n/q_n = ζ(2) + ζ(3)

## Numerical verification: TODO (large coefficients, needs mpmath with high precision — uisai2)

## Approach Ideas
- This is a 4-term version of 2.6 — likely more efficient convergence
- The very large initial values suggest this comes from a CMF decomposition
- Connection to Apéry's original approach for ζ(3)
- The polynomials A, B, C, D may factor in revealing ways
- Check if this recurrence is related to 2.6 by recurrence algebra
