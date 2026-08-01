# Ramanujan Challenge — Status (audited 2026-07-30)

**Deadline:** August 1, 2026, 23:59 UTC

**Submission set as of this audit:** 2.1, 2.3, 2.8, 3.1 — packaged under
`SUBMIT/`, combined archive `SUBMIT/dist/ramanujan-huang.zip` (2.0 MB).

## ⚠️ The old ✅ column was not accurate

Until 2026-07-30 this file marked 2.1–2.8 and 3.1 all "✅ unconditional".
Auditing the actual `.tex` files showed that several "complete proofs" contain
asserted steps. Two examples that were caught by reading the source:

- **2.2** — the decisive step, that the shifted operator equals Aptekarev's
  published operator, was written as *"can be checked by comparing sufficiently
  many evaluations"*. It was not done. Convergence was likewise asserted "by
  Poincaré–Perron theory adapted to the resonant case".
- **2.6** — the write-up itself says the identity is *"verified numerically to
  39 digits"*, and the "complete algebraic proof" is a sketch ending in
  *"reduces to a sum that decomposes as V + W"*.

The Lean side had the same problem in worse form; see `lean/UNDERSTANDING.md`.

The column below now records what has been **verified this session** versus what
is inherited and **not re-audited**. Nothing is marked done on the strength of an
earlier session's claim.

| Problem | Topic | Method | Audited? | State |
|---------|-------|--------|----------|-------|
| 2.1 | PCF → π | sign-flip of Cohen's Entry 5.3.22 | ✅ audited | **DONE — packaged in `SUBMIT/2.1/`** |
| 2.2 | γ Apéry | claimed = Aptekarev after shift m=n+3 | ✅ audited | **the identification is FALSE — see below** |
| 2.3 | π+e | tensor product: Lambert ⊗ derangement | ✅ audited | **DONE — packaged in `SUBMIT/2.3/`** |
| 2.4 | harmonic+polylog | weight-4 HPL symbolic summation | not re-audited | inherited claim |
| 2.5 | Catalan CMF | Delannoy decomposition + k-recurrence | not re-audited | inherited claim |
| 2.6 | ζ(2)+ζ(3) | GF ODE connection formula | ✅ audited | **gap: numerical + proof sketch** |
| 2.7 | 4-term ζ(2)+ζ(3) | rational gauge transfer from Zudilin | not re-audited | inherited claim |
| 2.8 | √10005/π | Chudnovsky in CMF disguise | ✅ (earlier) | **DONE — packaged in `SUBMIT/2.8/`** |
| 3.1 | knot π² | A-polynomial / Mahler measure | ✅ (earlier) | **DONE — packaged in `SUBMIT/3.1/`** |

### 2.2 — the identification is wrong

The write-up asserted that the challenge's initial values `(0,7,179)`/`(1,12,306)`
"are precisely the initial values of the Aptekarev rational approximants". They
are not. Aptekarev's published values (Hessami Pilehrood & Hessami Pilehrood,
arXiv:1010.1420, eq. (3)) are `p̃ = (0,2,31)`, `q̃ = (1,3,50)`, and his recurrence
is `(16n-15)y_{n+1} = (128n³+40n²-82n-45)y_n - n²(256n³-240n²+64n-7)y_{n-1} +
n²(n-1)²(16n+1)y_{n-2}` — coefficient degrees 1,3,5,5, not the challenge's
3,5,7,9. Rivoal's order-3 γ-recurrence (2009) does not match either: its
sequences are `Q = 1, 7, 65/2, 727/6, …` against the challenge's
`1, 12, 306, 13056, …`.

What *is* verified: the challenge's own sequences are integral and `p_n/q_n → γ`
(27.7 digits at n=60, subexponential — the Aptekarev-type signature). Also
checked and ruled out: neither `p` nor `q` satisfies any order-2 recurrence with
polynomial coefficients of degree ≤ 10, so the order-3 operator does not factor
the way 2.3's did. The right reference is probably Aptekarev–Tulyakov,
*Four-term recurrence relations for γ-forms*, which we have not obtained.
**2.2 is not submittable as it stands.**

### 2.1 — what closed it

The challenge PCF is the sign-flip of the tail of Cohen's Entry 5.3.22
(arXiv:2607.06581), which we retrieved and confirmed **verbatim**, including its
displayed quotients 42, 396, 1047, 38400, 4340. Two elementary proved steps: the
index shift `a_n = -α(n+1)`, `b_n = β(n)`, and a sign-flip lemma proved *at the
level of convergents* (`P̃_n = (-1)^{n+1}P_n`, `Q̃_n = (-1)^n Q_n`) so that no
tail-convergence question arises. Lean: 0 sorry, standard axioms only.

### 2.3 — what closed it

Order 4 = 2 × 2. The challenge operator annihilates every product
`X_{n+2} Y_{n+3}` of a Lambert-recurrence solution with a derangement-recurrence
solution; `m!` satisfies the *same* recurrence as `D_m`, which is why π and e
appear together and additively. Hence the **exact** splitting
`p_n/q_n = 4·B_{n+2}/A_{n+2} + (n+3)!/D_{n+3}`. The Lambert value is now proved
from a positive moment representation and geometric remainder bound. Lean: 0
sorry, standard axioms only, no `native_decide`; the main theorem is fully
unconditional.

## P3.2 — gcd(d_n a_n, d_n b_n) = e^{o(n)} for Apéry sequences

**Status:** ADDRESSED (10 pages). Two-tier proof + structural theory.

### Unconditional results
1. **Z(p) = O(p^{2/3}) for ALL primes** — effective bound, no excluded primes
   - No consecutive zeros of b_j mod p (backward induction from b_0=1)
   - Gap-h polynomial C_h(m) has degree 3(h-1), leading coeff U_{h-1}(17) (Chebyshev)
   - **Nonvanishing lemma** (two-point proof): N_h(-1) = b_{h-1}·((h-1)!)³,
     N_h(-2) = -5·b_{h-2}·((h-2)!)³ → if N_h ≡ 0 mod p, consecutive zeros → ⊥
   - Effective: Z(p) ≤ p/H + 3(H-1)(H-2)/2 + 1, optimize H = (p/3)^{1/3}
   - Reflection law: N_h(-m-h-1) = (-1)^{h-1} N_h(m)
2. **Content restriction:** p | cont(N_h) ⟹ h ≥ 2p (for p ≥ 7)
3. **Resultant coprimality:** Res(N_2, N_3) = -5^6 (coprime mod p for ALL p ≥ 7)
4. **Restart lemma:** N_r(x₀)=0, N_{r+1}(x₀)≠0 ⟹ N_{r+d}(x₀) = N_{r+1}(x₀)·N_d(x₀+r)
5. **Column O(H^{2/3}) bound:** #{h: N_h(x₀)≡0, h∈I} ≤ (3^{4/3}/2)H^{2/3}+O(H^{1/3})
6. **SL₂-tiling:** Dodgson relation τ_{h+1}(x)τ_{h+1}(x+1) - τ_{h+2}(x)τ_h(x+1) = 1
7. **Split-prime obstruction:** R_p(H)=O(H) uniformly is FALSE (Chebotarev)
8. **Density-1 unconditional:** G_n = e^{o(n)} for density-1 of n
   - First-moment argument over dyadic intervals using Z(p) = o(p)

### Conditional result (under Hypothesis Z)
5. Hypothesis Z-bar (average Z(p) = O(1)): log G_n = O(√n) for a set of n of natural
   density 1 (Theorem thm:main).  NOTE (2026-07-31): this is a density-1 statement, NOT a
   bound for all n; the earlier "for ALL n" phrasing here was wrong.  No upgrade to all n can
   follow from Z-bar plus the proved structure of Z_p: sets with |S_p| <= 2, reflection
   symmetry, no consecutive elements and bounded average can align at a single N and give
   T(N) ≍ N (explicit example S_p = {N-p, p-1-(N-p)}, verified at N = 20000).  The missing
   input is the uniform diagonal-discrepancy bound (DA); see
   research/working_notes/Q32_SESSION_2026-07-31_RESULTS.tex.

### Computational evidence (p ≤ 10^5, 9590 primes)
- Z(p) ∈ {0, 1, 2, 4, 6, 8, 10}, max = 10, mean ≈ 0.990
- P(Z=0) ≈ 61% ≈ e^{-1/2}; pair count Z/2 ∼ Poisson(1/2)
- Z(p) = 1 at exactly 2 non-ordinary primes (p | a_p(f), f = 8.4.a.a)
- Non-CM modular form → no systematic forced zeros (contrast: ζ(2) Apéry is CM)

### Key scripts
- `scripts/p32_zp_extended.py` — Z(p) for p ≤ 10^5
- `problems/3.2/proof.tex` — 10-page proof (compiles clean)
