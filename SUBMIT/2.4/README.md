# Problem 2.4 — harmonic numbers, polylogarithms and zeta values

**Status: proved, unconditionally.**

```lean
theorem problem24_unconditional : problem24Statement
-- depends on axioms: [propext, Classical.choice, Quot.sound]
```

## Contents

- `solution.pdf` / `solution.tex` — the mathematical write-up (5 pages).
- `lean/RamanujanChallenge/` — four Lean 4 source files.

## What is proved

The double sum is reduced to seven scalar series of weight at most four, and all
seven are evaluated. Writing
`P_m = H_m + 2·Σ_{k<m} (-1)^{k+1}/(k+1)`, so that `P_{2m} = 2H_m − H_{2m}` is the
harmonic remainder of the creative-telescoping certificate:

| series | value |
|---|---|
| `Σ (P_m² − H_m^{(2)})/m²` | `20 Li₄(½) + (5/6)L⁴ + 7L²Z₂ − (59/10)Z₂²` |
| `Σ (-1)^m (P_m² − H_m^{(2)})/m²` | `−22 Li₄(½) − (11/12)L⁴ − (13/2)L²Z₂ − (7/4)LZ₃ + (67/10)Z₂²` |
| `Σ P_m/m³` | classical cubic-linear |
| `Σ (-1)^m P_m/m³` | classical alternating cubic-linear |
| shifted linear Euler sum | classical |
| Leshchiner boundary sum | `(7/8)ζ(4)` |
| Borwein–Bradley–Broadhurst boundary sum | `ζ(4)` |

`L = log 2`, `Z_r = ζ(r)`. The certificate
`problem24_of_euler_and_classical` takes these seven as hypotheses; supplying all
seven gives `problem24Statement` with no hypothesis left.

Three of the evaluations appear to be new in this form: the two weight-four
quadratic Euler sums in the first two rows, and the reduction of the two
inverse-central-binomial boundary sums to the pair
`Σ n⁻⁴C(2n,n)⁻¹ = 17π⁴/3240` (classical, Comtet) and
`Σ H_{n-1}^{(2)} n⁻²C(2n,n)⁻¹ = 5π⁴/9720`.

The two quadratic sums are not related by a sign change. The coefficient
integral evaluates the generating function `Q(x) = 2J(x)/(1−x)` at `−x` in the
alternating case, keeping the argument in `(−1,0)`; the non-alternating case runs
it into the pole at `x = 1`. Half the pole cancels against `−log x ∼ 1−x`, but
`J` grows like `log²(1−x)`, so that direction carries an integrable `log²`
singularity the other does not. Its evaluation goes by parts with
`U(x) = log²x/2 + Li₂(1−x)` against `D⁺(x) = J'(x)`, and the `log2·ζ(3)` term
cancels identically against the generator's contribution — which is why the
non-alternating value carries no `ζ(3)`.

## Verification

Toolchain `leanprover/lean4:v4.29.0`, Mathlib pinned in `lake-manifest.json`.
The four files contain **no `sorry`** and introduce **no axiom**. Every theorem
named above depends only on `propext`, `Classical.choice`, `Quot.sound`.

The Lean sources are provided as the development they were written in. Building
them requires the accompanying `lakefile.toml`, `lean-toolchain` and
`lake-manifest.json` in `lean/`, and `Problem24QuadraticAlt.lean` additionally
imports the project's `Problem26*` and `Dilogarithm` modules, which are not
reproduced here. The claims above were checked in the full development, not in
this excerpt.

Independently of the formalization, the six endpoint values were obtained by
three mutually independent analytic derivations and cross-checked against
high-precision quadrature; the two Möbius relations were verified exactly over
`ℚ` against those values before being formalized.

## Scope note

The write-up's Stage 1 — the creative-telescoping certificate and the closed form
for the inner sum — is presented as a derivation, not as a formalized result.
The formalization begins at the seven scalar evaluations and the certificate that
assembles them.
