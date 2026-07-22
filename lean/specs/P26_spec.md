# Codex Spec: Problem 2.6 (ζ(2)+ζ(3)) — Close All Sorries

## Target
Close all sorry's in `RamanujanChallenge/Problem26.lean`:
1. `recessiveRatio_limit` — v_{n+1}/v_n → 1/4
2. `zeta2_eq` — ζ(2) = π²/6
3. `problem26_identity` — the main theorem

## Working directory
`~/repos/Ramanujan_Challenge/lean/`

## Mathematical proof (remainder-certificate route)

The three-term recurrence A(n)u_n = B(n)u_{n-1} - C(n)u_{n-2} has:
- Poincaré roots 1 and 1/4 (already proved in the file: `poincare_poly_roots`)
- Recessive solution v_n = 1120·((n+2)!)²/(2n+6)! with ratio → 1/4
- Dominant solution u_n from variation of constants

### Strategy for `zeta2_eq`
Mathlib has `hasSum_zeta_two : HasSum (fun n => 1/(n+1)^2) (π²/6)`.
Bridge from `HasSum` to `tsum` using `HasSum.tsum_eq`.
The exact statement to close: `∑' n, 1/(↑n+1)^2 = π²/6`.

### Strategy for `recessiveRatio_limit`
The ratio is (n+3)²/(2(n+4)(2n+7)). As n → ∞ this is ~ n²/(4n²) = 1/4.
Use `Filter.Tendsto.div` with polynomial asymptotics.
Specifically: write as ((n+3)/(2n+7)) · ((n+3)/(2(n+4))), both → 1/2, product → 1/4.
Or directly: clear to rational function, use `Filter.Tendsto.ratCast` or manual manipulation.

### Strategy for `problem26_identity`
This is the hardest theorem. The remainder-certificate approach:

1. **Define the recessive solution explicitly:**
   v_n := 1120 * ((n+2)!)² / (2n+6)!
   Prove v satisfies the recurrence (by clearing factorials, reducing to polynomial identity).

2. **Casorati determinant:**
   W_n := u_n·v_{n+1} - u_{n+1}·v_n (Wronskian-type).
   For a 2nd-order recurrence, W_{n+1}/W_n = C(n+1)/A(n+1) (known, computable).
   This gives W_n in closed form.

3. **Variation of constants:**
   The dominant solution u_n can be written as u_n = v_n · ∑_{j≥n} f(j),
   where f(j) = 1/(v_j·v_{j+1}·(C(j+1)/A(j+1))...) — exact closed form.

4. **Remainder:**
   r_n = q_n·(ζ(2)+ζ(3)) - p_n = (exact tail sum involving v_n).
   Bound |r_n| by geometric decay (v ratio → 1/4 < 1).

5. **Apply `tendsto_div_of_remainder_tendsto_zero`** from RemainderCertificate.lean.

## Alternative simpler approach
Instead of the full Casorati machinery, you can:
1. Define p_n, q_n by the recurrence with initial values from the file.
2. Define r_n = q_n·C - p_n where C = ζ(2)+ζ(3).
3. Prove r_n satisfies the SAME recurrence (by `ring` — this is in RemainderCertificate.lean).
4. Compute r_0, r_1 from initial values (these will be specific rational numbers involving π²/6 and ζ(3)).
5. Show |r_{n+1}/r_n| → 1/4 (the recessive root), since r_n IS the recessive solution.
6. Geometric decay + remainder-certificate convergence.

The key challenge is: computing r_0, r_1 requires knowing the exact value of ζ(2)+ζ(3),
which is circular UNLESS we can verify the initial values numerically and then use the
recurrence structure to bootstrap.

## What the file already has
- Recurrence coefficients `coeff_A`, `coeff_B`, `coeff_C` (ℤ → ℤ)
- `poincare_poly_roots` — roots are 1 and 1/4 (PROVED)
- `recessiveRatio` — the ratio function (ℕ → ℚ)
- Initial values u₁, u₂ as rationals

## Hard rules
- No sorry, no axiom, no native_decide
- Use `hasSum_zeta_two` from Mathlib for ζ(2) = π²/6
- For ζ(3), define as `∑' n : ℕ, 1 / (↑n + 1) ^ 3` and relate to `riemannZeta 3`

## Verification
```bash
export PATH="$HOME/.elan/bin:$PATH"
lake env lean RamanujanChallenge/Problem26.lean 2>&1
```

## Stall protocol
If the full main theorem is too hard, close `zeta2_eq` and `recessiveRatio_limit` first,
then report the exact stall point on `problem26_identity`.
