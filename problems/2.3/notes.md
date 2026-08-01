# Problem 2.3: π + e as Apéry Limit

## Statement

5-term recurrence (order 4):
```
0 = (-n³ + 2n² + 7n + 3) u_n
  + (n+2)(2n⁴ + n³ - 26n² - 48n - 19) u_{n-1}
  + (n+2)(n⁶ + 9n⁵ + 8n⁴ - 87n³ - 249n² - 234n - 68) u_{n-2}
  + (n+1)²(n+2)(2n⁵ + 3n⁴ - 13n³ - 21n² + 4) u_{n-3}
  - n³(n+1)²(n+2)(n³ + n² - 8n - 11) u_{n-4}
```

Initial values:
- p_{-3}=1, p_{-2}=1, p_{-1}=20, p_0=296
- q_{-3}=1, q_{-2}=0, q_{-1}=4,  q_0=48

Prove: lim p_n/q_n = π + e

## Numerical verification: ✅ (matches to 50 digits at N=80)

## Key Findings

### Poincaré roots: {-1±√2, 1±√2} at two factorial scales
- (n/e)^{2n} modes: roots -1±√2 (double-factorial)
- (n/e)^n modes: roots 1±√2 (single-factorial)

### Ore algebra analysis (uisai2 Sage, 44 min computation)
- **REDUCIBLE but NOT COMPLETELY REDUCIBLE** (not semisimple)
- One 2D invariant subspace: Sol(P) annihilates (n/e)^{2n} pair
- L = Q·P as product factorization (left factor Q, right factor P)
- But NO complementary 2D invariant subspace — the (n/e)^n solutions
  are ENTANGLED with the double-factorial part
- Complement has non-rational (factorially-growing) coefficients

### LCLM test results
- LCLM of two order-2 operators in Q[n]<S>: FAILS (ChatGPT Q4639)
- Product factorization L = Q·P: EXISTS (uisai2 ore_algebra)
- But P and Q don't cleanly separate π from e

### Superseded false Lambert lead
- The debunked lead used the different recurrence implemented in
  `verify_decomposition.py`; its ratio tends to a Bessel-value quotient and its
  initial values do not match.
- The final proof uses `X_m = (2m+1)X_{m-1}+m²X_{m-2}`. Its tensor identity and
  all initial values are kernel-checked in `lean/RamanujanChallenge/Problem23.lean`.
  A positive moment representation now proves its ratio tends to `π/4`, so the
  final theorem is fully unconditional.

## Approach Ideas (updated)
- The non-semisimple structure suggests a Jordan-block phenomenon
- π+e is genuinely entangled at the operator level
- Need: CMF embedding or non-obvious integral representation
- Possible: the proof uses the PRODUCT structure L=Q·P differently
  from a simple decomposition
