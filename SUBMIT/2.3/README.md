# Ramanujan Challenge — Problem 2.3

**The sum π + e as an Apéry limit**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

---

## What is in this package

```
solution.pdf / solution.tex   human-readable proof (REQUIRED ARTIFACT), 7 pages
lean/                         self-contained Lean 4 formalization
verify.py                     exact-arithmetic certificates, reproducible
README.md                     this file
```

## Summary of the result

The challenge gives an order-four recurrence and two of its solutions `p_n`,
`q_n`, and asks for `lim p_n/q_n = π + e`.

The proof is structural rather than asymptotic. **Order 4 = 2 × 2**: the
challenge operator is the tensor product of two order-two systems,

| system | recurrence | solutions | ratio |
|---|---|---|---|
| Lambert `L` | `X_m = (2m+1)X_{m-1} + m² X_{m-2}` | `A_m`, `B_m` | `B_m/A_m → π/4` |
| derangement `D` | `Y_m = (m-1)(Y_{m-1} + Y_{m-2})` | `D_m`, `m!` | `m!/D_m → e` |

The observation that makes the problem work is that **`m!` satisfies the same
recurrence as the derangement numbers** — `(m-1)((m-1)! + (m-2)!) = m!` — so the
`D`-system carries both the derangement numbers and the factorials. That is why
π and e occur together here, and why they occur *additively*.

The main theorem is that *every* product `X_{n+2} Y_{n+3}` of a solution of `L`
with a solution of `D` satisfies the challenge recurrence. This is an identity
with four free initial values, not a statement about particular sequences. The
challenge's own `p_n, q_n` are then identified as

```
q_n = A_{n+2} D_{n+3}
p_n = 4 B_{n+2} D_{n+3} + A_{n+2} (n+3)!
```

(the closed forms satisfy the recurrence, match the eight prescribed initial
values, and the leading coefficient `c_0(n) = -n³+2n²+7n+3` never vanishes at a
positive integer, so the solution is unique). Dividing gives the **exact**
splitting

```
p_n / q_n  =  4 · B_{n+2}/A_{n+2}  +  (n+3)!/D_{n+3}   →   4·(π/4) + e  =  π + e.
```

No Poincaré–Perron or Birkhoff asymptotics enter the proof of the limit; the
splitting is an identity valid term by term.

## The Lean formalization, and exactly what it proves

`lean/` is a **self-contained** Lake project: `Ramanujan23/Problem23.lean`
imports only `Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic` and
`Mathlib.Combinatorics.Derangements.Exponential`, and nothing outside `lean/`.

**Build status.**

```
Build completed successfully (2071 jobs).   0 errors, 0 sorry
```

against the pinned toolchain `leanprover/lean4:v4.30.0` and Mathlib
`c5ea00351c28…` (tag `v4.30.0`), fixed in `lean-toolchain` and
`lake-manifest.json`.

*Scope of what we verified:* the project was built as a standalone Lake package
with these pins, from these sources alone. We did not additionally rebuild it
from scratch on a fresh machine with a cold Mathlib cache; the Mathlib build
artifacts were reused from the same pinned revision. To reproduce from nothing:

```bash
cd lean
lake exe cache get      # fetches the matching Mathlib build
lake build
lake env lean AxiomCheck.lean
```

**What is proved.**

| Lean name | statement |
|---|---|
| `tensor_rec` | for **any** `L`-solution `X` and **any** `D`-solution `Y`, `X_m·Y_m` satisfies the challenge recurrence |
| `factorial_isDerRec` | `m!` solves the derangement recurrence |
| `challengeQ_rec`, `challengeP_rec` | the closed forms solve the challenge recurrence |
| `challengeQ_zero … challengeP_three` | all eight prescribed initial values |
| `C0_ne_zero` | `c_0(n) ≠ 0` for every integer `n ≥ 1` |
| `eq_of_satisfiesRec` | a solution over a char-0 field is determined by its first four values |
| `ratio_split` | the exact splitting identity above |
| `factorial_div_derang_tendsto_exp_one` | `m!/D_m → e` |
| `problem23_pi_add_e` | the main theorem, given the Lambert hypothesis |

**Axiom audit.** `#print axioms` on the verified build (`AxiomCheck.lean`)
reports **only the three standard axioms**, and less for the algebraic core:

| theorem | axioms |
|---|---|
| `factorial_isDerRec` | `propext` |
| `tensor_rec`, `challengeQ_rec`, `challengeP_rec` | `propext, Quot.sound` |
| `C0_ne_zero`, `eq_of_satisfiesRec`, `ratio_split` | `propext, Classical.choice, Quot.sound` |
| `factorial_div_derang_tendsto_exp_one` | `propext, Classical.choice, Quot.sound` |
| `problem23_pi_add_e` | `propext, Classical.choice, Quot.sound` |

There is **no `native_decide`** anywhere in this development, and no `sorry`.

**The one hypothesis that remains.** The top-level theorem is stated
**conditionally** on the classical value of Lambert's continued fraction:

```lean
theorem problem23_pi_add_e
    (hLambert : Tendsto (fun m => (lambertB m : ℝ) / (lambertA m : ℝ))
                  atTop (𝓝 (Real.pi / 4))) :
    Tendsto (fun m => (challengeP m : ℝ) / (challengeQ m : ℝ))
      atTop (𝓝 (Real.pi + Real.exp 1))
```

`lambertA`, `lambertB` are the continuant sequences of

```
arctan z = z / (1 + z²/(3 + 4z²/(5 + 9z²/(7 + …)))),      π/4 = arctan 1,
```

whose partial numerators `m²` and partial denominators `2m+1` are exactly the
coefficients of the `L`-recurrence. The convergence and the value are due to
Lambert (1761) and Legendre (1794); see Wall, *Analytic Theory of Continued
Fractions* §92, or Jones–Thron Ch. 6, for modern proofs via Gauss's continued
fraction for `₂F₁` ratios. We deliberately kept this as an explicit hypothesis
rather than an axiom or a `sorry`, so that the dependency is visible in the
statement itself.

What the submission *does* establish unconditionally about that limit
(§5 of `solution.pdf`): from the Casorati determinant
`A_m B_{m-1} - A_{m-1} B_m = (-1)^{m+1} (m!)²` we get

```
B_m/A_m  =  1 + Σ_{j=1}^{m} (-1)^j (j!)² / (A_j A_{j-1}),
```

an alternating series whose terms are proved to decrease strictly to 0 with
ratio ≤ 1/5. So the *existence* of the limit is unconditional and elementary;
what is imported from the classical literature is its *value*.

The `e` half is fully formal: our `derang` is proved equal to Mathlib's
`numDerangements`, and the limit comes from Mathlib's
`numDerangements_tendsto_inv_e`.

## Reproducing the numerics

```bash
python3 verify.py
```

`sympy` and `mpmath` are optional; the script says which checks it skipped.
It performs three checks:

1. **Faithfulness.** Forward-solves the challenge recurrence over ℚ *from the
   problem statement's initial values alone*, confirms every value is an
   integer, and confirms agreement with the closed forms for `n = -3 … 40`.
   First new values: `q_1 = 1836`, `q_2 = 97680`, `p_1 = 10656`, `p_2 = 573344`.
2. **The theorem, not an instance.** Re-verifies the tensor identity
   symbolically with four *free* initial values, so that agreement of one
   particular solution cannot be mistaken for the general identity.
3. **The limit.** At `n = 50`, `p_n/q_n` agrees with
   `π + e = 5.859874482048838473822930854632…` to 40 decimal places
   (error `1.4e-41`), consistent with the geometric rate `(3-2√2)^n`.

## References

- J. H. Lambert, *Mémoire sur quelques propriétés remarquables des quantités
  transcendantes circulaires et logarithmiques*, Hist. Acad. Roy. Sci. Berlin
  **17** (1761/1768), 265–322.
- A.-M. Legendre, *Éléments de géométrie*, Paris, 1794; Note IV.
- H. S. Wall, *Analytic Theory of Continued Fractions*, Van Nostrand, 1948.
- W. B. Jones and W. J. Thron, *Continued Fractions: Analytic Theory and
  Applications*, Addison–Wesley, 1980.
- S. Khrushchev, *Orthogonal Polynomials and Continued Fractions: From Euler's
  Point of View*, Cambridge University Press, 2008.
