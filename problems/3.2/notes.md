# Problem 3.2: Optimality of Apéry's Irrationality Measure for ζ(3) (CONJECTURE)

## Statement

The sequences (a_n), (b_n) satisfy Apéry's recurrence:
```
(n+1)³ u_{n+1} - (34n³ + 51n² + 27n + 5) u_n + n³ u_{n-1} = 0,  n ≥ 1
```

Initial values: a_0=0, a_1=6, b_0=1, b_1=5.

Define d_n := lcm(1,...,n)³. Note d_n·a_n, d_n·b_n ∈ Z (Apéry 1979).

Prove: gcd(d_n·a_n, d_n·b_n) = e^{o(n)}.

## Context
The a_n, b_n are Apéry's sequences for ζ(3) irrationality.
The linear form a_n - ζ(3)·b_n → 0 gives the irrationality measure.
If gcd(d_n·a_n, d_n·b_n) = e^{o(n)}, then Apéry's bound μ(ζ(3)) ≤ 1 + ... is optimal.

## Approach Ideas
- This is a deep number-theoretic conjecture
- Connected to p-adic properties of Apéry numbers
- Check: supercongruences for Apéry numbers (Beukers, Stienstra-Beukers)
- The gcd condition is related to absence of "unexpected" common factors
- Computationally: verify for many n that the gcd is small
- Theoretical: may need p-adic analysis for each prime p
- Connection to HolonomicCRN: Apéry recurrence is the canonical holonomic example
- This is likely open for a reason — very hard
