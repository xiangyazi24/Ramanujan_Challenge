# Ramanujan Challenge — Problem 2.8

**Very fast rational approximation of √10005/π**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

---

## What is in this package

```
solution.pdf / solution.tex   human-readable proof (REQUIRED ARTIFACT)
lean/                         self-contained Lean 4 formalization
verify.py                     numerical certificates, reproducible
README.md                     this file
```

## Summary of the result

The challenge presents a conservative-matrix-field (CMF) recurrence and asks to
prove that its convergents `P_{N,j}/Q_{N,j}` tend to `√10005/π` for `j = 1,2,3,4`.

The proof has three layers:

1. **Parameter identification.** The evaluation point is identified as the
   Chudnovsky CM singular modulus of discriminant `-163`: the parameter
   `R = 1 - j(τ₁₆₃)/1728`, the cubic `w = (2k+1)(6k+1)(6k+5)`, and
   `640320 = 8² · 10005`.
2. **Poincaré root bridge.** The scalar order-4 recurrence extracted from the
   matrix field has Poincaré root `ρ = 64(R-1) = 64 · 640320³/1728`, matching the
   Chudnovsky hypergeometric ratio `h_{k+1}/h_k → -1/(R-1)` at the same point.
   Hence the CMF convergents share the limit and the geometric rate of the
   Chudnovsky partial sums.
3. **Classical evaluation.** The Chudnovsky series equals `√10005/π`
   (Chudnovsky–Chudnovsky 1988; Borwein–Borwein 1987).

Layers 1 and 2 are the contribution of this submission. Layer 3 is a published
theorem and is cited as such.

## The Lean formalization, and exactly what it proves

`lean/` contains a **self-contained** development (27 modules, ~30 000 lines)
that formalizes Layer 3 down to a single isolated analytic evaluation. The
dependency closure was computed from `Ripple.Number.Chudnovsky1989` and isolated
from our larger `Ripple` project; **self-containedness was then verified by
unpacking this directory alone on a clean machine and building it**, which
reported `Build completed successfully (3678 jobs)` with zero errors. Nothing
outside `lean/` is referenced.

**Build status.** The development builds cleanly:

```
Build completed successfully (3675 jobs).   0 errors, 0 sorry
```

against the pinned toolchain `leanprover/lean4:v4.29.0` and Mathlib
`8a178386ffc0f5fef0b77738bb5449d50efeea95` (tag `v4.29.0`), both fixed in
`lean-toolchain` and `lakefile.toml`. *(It does not build against Mathlib
v4.30.0: `ModularForm.delta` was deprecated in favour of
`ModularForm.discriminant` on 2026-03-23 and the level-one group presentation
changed. The pin above is the version the development was verified against.)*

**Axiom audit — please read this, it is not uniform across the theorems.**
`#print axioms` on the verified build reports:

| theorem | axioms |
|---|---|
| `class_number_neg_163_eq_one` | `propext, Classical.choice, Quot.sound` |
| `chudnovsky_one_over_pi` | `propext, Classical.choice, Quot.sound` |
| `a_eq_3F2_coeff` | `propext, Classical.choice, Quot.sound` |
| `kleinJ_heegnerTau163_eq_heegnerJ163Target_unconditional` | the three standard axioms **plus five `native_decide` axioms** |

So:

- The class-number-one theorem, the `₃F₂` bridge, and the top-level conditional
  Chudnovsky theorem are **fully kernel-checked**.
- **`j(τ₁₆₃) = -640320³` is proved, but its proof invokes `native_decide`** at
  five points (degree and diagonal-coefficient facts about the sparse `Φ₄₁`
  terms, the level-41 Sturm bound, and a cofactor tail estimate). `native_decide`
  compiles the decision procedure to native code and trusts the compiler and
  runtime, so it enlarges the trusted base relative to kernel reduction. It is
  not a `sorry`, and Mathlib-style `decide` is infeasible on certificate arrays
  of this size, but we state the dependency explicitly rather than claim a
  kernel-checked result. *Mathlib does not contain this evaluation at all.*

**Also proved (kernel-checked):** the Chudnovsky coefficient recurrence, the
identification `a_k = ₃F₂` coefficients, absolute summability, the reduction of
the Chudnovsky linear series to a derivative series of `₃F₂`, and **Clausen's
identity** converting that `₃F₂` into the square of a Gauss `₂F₁`.

**The one hypothesis that remains.** The top-level theorem

```lean
theorem chudnovsky_one_over_pi :
    Hypergeometric.chudnovskyCM163GaussDerivativeCombination
        = (640320 : ℝ)^(3/2) / (12 * Real.pi) →
    (640320 : ℝ)^(3/2) / (12 * Real.pi)
      = ∑' k : ℕ, (-1)^k * a k * (13591409 + 545140134 * k) / (640320 : ℝ)^(3*k)
```

is stated **conditionally** on the classical CM period-derivative evaluation

```
13591409 · F(x)²  +  545140134 · x · 2 F(x) F'(x)  =  640320^(3/2) / (12π),
      F = ₂F₁(1/12, 5/12; 1; ·),      x = -1728/640320³.
```

We do not formalize this last evaluation; it is the analytic content of the
cited Chudnovsky theorem. Everything else in the chain is formal. We deliberately
kept it as an explicit hypothesis rather than an axiom or a `sorry`, so that the
dependency is visible in the statement itself.

`verify.py` checks this hypothesis numerically to **130 decimal digits**
(relative error `1.4e-131`).

### Reduction of the remaining hypothesis (`cm_reduction.py`)

We further reduce that single hypothesis to three named classical inputs and one
exact rational identity. Write `tau = tau₁₆₃ = (1 + i√163)/2`, `z = 1728/j(tau)`.
The classical inversion `F = E₄^{1/4}` together with

```
dz/dtau = 2πi · z · E₆/E₄        (from j' = -2πi j E₆/E₄)
E₄'     = (2πi/3)(E₂E₄ - E₆)     (Ramanujan)
```

gives `z · 2FF' = (1/6) E₄^{1/2} (E₂E₄/E₆ - 1)`, hence

```
A F² + B z · 2FF' = E₄^{1/2} [ A - B/6 + (B/6) E₂E₄/E₆ ].
```

Substituting the non-holomorphic completion `E₂ = E₂* + 3/(π Im tau)` — **the
only point at which 1/π enters** — and writing `s₂ = E₂* E₄/E₆`:

```
= E₄^{1/2} [ (A - B/6 + (B/6) s₂) + (B/(2π Im tau)) E₄/E₆ ].
```

The bracket's algebraic part must vanish, and the surviving term must match the
target. So the hypothesis is equivalent to the conjunction of

- **(C1)** `s₂(tau₁₆₃) = 1 - 6A/B = 77265280/90856689`;
- **(C2)** `B/√163 · E₄^{3/2}/E₆ = 640320^{3/2}/12`, which since
  `E₄³/E₆² = j/(j-1728)` is purely algebraic in `j`.

**(C2) is discharged exactly.** In rational arithmetic

```
B²/163 · j/(j-1728)  =  640320³/144  =  1823176476672000,
```

an identity between rational numbers, hence decidable — no numerics involved.

**(C1)** is the classical rationality of `s₂` at a class-number-one CM point,
together with its value; `cm_reduction.py` confirms it against `q`-expansions to
80 digits. This is a genuine test rather than a tautology: the naive
approximation `1 - 6/(π√163) = 0.8504082731872**41179...**` already departs from
`s₂ = 0.8504082731872**38861...**` in the 13th decimal, so the `q`-corrections
are doing real work.

**Consequence.** No Chowla–Selberg input is required. What remains classical is
only (i) the hypergeometric–modular inversion `F = E₄^{1/4}`, (ii) Ramanujan's
derivative formulas, and (iii) the value of `s₂` at `tau₁₆₃`.

## Building the Lean development

```bash
cd lean
lake exe cache get      # fetches the matching Mathlib build
lake build Ripple.Number.Chudnovsky1989
```

Toolchain is pinned in `lean-toolchain` (`leanprover/lean4:v4.29.0`) and Mathlib
in `lakefile.toml` (`v4.29.0`, commit `8a178386...`). No files outside this
directory are required.

Verified twice on a clean Linux host:

* the full upstream project at this pin: `Build completed successfully (3675 jobs)`;
* **this directory alone, unpacked from the submitted archive**:
  `Build completed successfully (3678 jobs)`.

Both with zero errors and zero `sorry`.

To inspect the axiom dependencies of the singular-modulus theorem:

```lean
import Ripple.Number.Modular.CMEvaluation163
open Ripple.Number.Modular
#print axioms kleinJ_heegnerTau163_eq_heegnerJ163Target_unconditional
```

## Reproducing the numerics

```bash
python3 verify.py
```

Requires `mpmath`. Prints the 130-digit check of the CM evaluation, the Poincaré
root comparison, and the CMF convergent rate.

## References

- D. V. Chudnovsky and G. V. Chudnovsky, *Approximations and complex
  multiplication according to Ramanujan*, in Ramanujan Revisited, Academic Press
  (1988), 375–472.
- J. M. Borwein and P. B. Borwein, *Pi and the AGM*, Wiley (1987).
- S. Ramanujan, *Modular equations and approximations to π*, Quart. J. Math. 45
  (1914), 350–372.
