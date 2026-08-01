# Problem 2.2: Euler's γ as an Apéry Limit

## Status (2026-08-01)

Unconditional. The previous Aptekarev-identification route was false and has
been replaced. The Lean theorem is
`RamanujanChallenge.P22.problem22_solved`; it has 0 `sorry` and its axiom audit
is `{propext, Classical.choice, Quot.sound}`.

## Exact route

1. Shift the challenge sequences by three indices.
2. Prove an exact first-order Ore transform from Rivoal's order-3 recurrence:

   ```text
   C_n = (T_{n+1} + (n+1)(3n+4)T_n)/(8n+11).
   ```

   The transformed initial triples are exactly `(0,7,179)` and `(1,12,306)`.
3. Prove finite hypergeometric formulas for the unscaled Rivoal solutions:

   ```text
   w(n,k) = choose(n,k)^2/k!
   B_n = Σ_k (2n+k+1) w(n,k)
   A_n = Σ_k w(n,k)((2n+k+1)(3H_k-2H_{n-k})-1).
   ```

   The recurrence identities and all boundary terms are discharged by exact WZ
   certificates in `Problem22.lean`.
4. Normalize the positive weights `W(n,k)=(2n+k+1)w(n,k)`. Their adjacent
   ratio gives the finite birth-death identity

   ```text
   W(n,k) B(n,k) = W(n,k+1) D(n,k+1),
   B(n,k) = (2n+k+2)(n-k)^2,
   D(n,k) = (2n+k)k^3.
   ```

5. Set `G(n,k)=k^3-(n-k)^2`. The Stein identity gives

   ```text
   E[G(n,K)^2] ≤ 81 n^3 sqrt(n).
   ```

6. On `|G|≤δn²`, the cubic saddle forces both `k` and `n-k` large and
   `k³/(n-k)²` close to one. Using
   `0≤H_m-log(m)-γ≤1/m` gives error at most `5/M+8δ`.
7. Off this set, use the envelope `6+5 log n` and the second moment. The final
   quantitative estimate is

   ```text
   weighted harmonic error
     ≤ 5/M + 8δ + 81(6+5 log n)/(δ² sqrt(n)).
   ```

   First choose `δ`, then `M`, then let `n→∞`.
8. The rational correction is at most `1/(2n+1)`. The Ore transform preserves
   the limit because its quotient is a positive weighted average of two
   adjacent Rivoal quotients.

## Files

- `lean/RamanujanChallenge/Problem22.lean`: recurrence, WZ sums, Ore transform,
  limit transfer.
- `lean/RamanujanChallenge/Problem22Concentration.lean`: adjacent ratio and
  finite Stein identity.
- `lean/RamanujanChallenge/Problem22Moment.lean`: polynomial moment bounds.
- `lean/RamanujanChallenge/Problem22Harmonic.lean`: harmonic transfer and
  unconditional theorem.
- `lean/P22Audit.lean`: axiom audit.
- `problems/2.2/proof.tex`: matching standalone proof.
