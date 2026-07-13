# Problem 2.2: Euler's γ as Apéry Limit

## Statement

4-term recurrence (order 3):
```
0 = (-8n³ - 51n² - 105n - 68) u_n
  + (24n⁵ + 337n⁴ + 1833n³ + 4818n² + 6092n + 2928) u_{n-1}
  - (n+2)(n+3)(24n⁵ + 273n⁴ + 1150n³ + 2154n² + 1635n + 268) u_{n-2}
  + (n+1)(n+2)⁴(n+3)(8n³ + 75n² + 231n + 232) u_{n-3}
```

Initial values:
- p_{-3}=0, p_{-2}=7, p_{-1}=179
- q_{-3}=1, q_{-2}=12, q_{-1}=306

Prove: lim p_n/q_n = γ (Euler-Mascheroni constant)

## Numerical verification: ✅ (matches to 40 digits at N=100)

## Approach Ideas
- Apéry-limit style: find integral representation for p_n - γ q_n
- The recurrence is order 3 — check if it's a known Apéry-like recurrence
- Connection to Chamberland et al. (2021) "Apéry Limits: Experiments and Proofs"
- Check Pilehrood-Pilehrood work on γ
- Creative telescoping / Zeilberger
