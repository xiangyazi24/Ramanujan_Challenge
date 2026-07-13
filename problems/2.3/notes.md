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

## Approach Ideas
- π + e appearing together is very unusual — no known closed-form integral for this
- The recurrence is order 4 — significantly harder
- This is likely the hardest problem in Section 2
- Check: does the recurrence factor into two order-2 recurrences (one for π, one for e)?
- Differential Galois theory perspective
- Possible connection to periods of mixed motives
