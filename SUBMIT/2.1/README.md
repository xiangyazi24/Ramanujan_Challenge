# Ramanujan Challenge — Problem 2.1

**Polynomial continued fraction for 6/(3−π)**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

---

## What is in this package

```
solution.pdf / solution.tex   human-readable proof (REQUIRED ARTIFACT), 4 pages
lean/                         self-contained Lean 4 formalization
verify.py                     exact-arithmetic certificates, reproducible
README.md                     this file
```

## Summary of the result

With `a_n = -220n³ - 484n² - 301n - 42` and `b_n = 4n²(2n+1)²(5n-4)(5n+6)`, the
challenge asks for

```
a_0 + b_1/(a_1 + b_2/(a_2 + b_3/(a_3 + …)))  =  6/(3-π).
```

The solution is short. Set

```
α(n) = 220n³ - 176n² - 7n + 5,     β(n) = 4n²(2n+1)²(5n-4)(5n+6).
```

**Step 1 (index shift, proved).** `a_n = -α(n+1)` and `b_n = β(n)`, both
polynomial identities.

**Step 2 (sign flip, proved).** If a continued fraction has convergents
`P_n/Q_n`, then negating every partial denominator while keeping every partial
numerator produces convergents

```
P̃_n = (-1)^{n+1} P_n,      Q̃_n = (-1)^n Q_n,      so   P̃_n/Q̃_n = -(P_n/Q_n).
```

This is proved **at the level of convergents**, by induction, with no
convergence hypothesis at all — which is what makes the argument airtight.
Sign manipulations of continued fractions are usually done on tails, where they
are delicate; here nothing is assumed and the two fractions converge or diverge
together by construction.

**Step 3 (the classical input, cited).** Cohen's *Continued Fractions of
Polynomial Type: Theory and Encyclopedic Dictionary* (arXiv:2607.06581),
Entry 5.3.22:

```
π = 3 + 6/(α(1) + β(1)/(α(2) + β(2)/(α(3) + …)))
  = 3 + 6/(42 + 396/(1047 + 38400/(4340 + …))).
```

By Step 1 the challenge's fraction is exactly the sign-flip of the tail
`T = 6/(π-3)`, so by Step 2 its value is `-T = 6/(3-π)`. ∎

## What we verified about the cited entry

We did not take Entry 5.3.22 on trust. We retrieved the source
(arXiv:2607.06581) and confirmed the entry verbatim — it is recorded there as

```
[()->Pi,[3,220*n^3-176*n^2-7*n+5],
         [6,4*n^2*(2*n+1)^2*(5*n-4)*(5*n+6)]]
```

with the displayed expansion beginning `42, 396, 1047, 38400, 4340`. `verify.py`
re-derives each of those five numbers from the formulas, and confirms the stated
value numerically: the truncation at n = 61 gives π to within `5.4e-129`.

## The Lean formalization, and exactly what it proves

`lean/` is a **self-contained** Lake project; `Ramanujan21/Problem21.lean`
imports only `Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic` and
nothing outside `lean/`.

**Build status.**

```
Build completed successfully (1915 jobs).   0 errors, 0 sorry
```

against `leanprover/lean4:v4.30.0` with Mathlib pinned in `lake-manifest.json`.
As with our other packages, this was verified as a standalone Lake project from
these sources; we did not additionally rebuild from a cold Mathlib cache on a
fresh machine. To reproduce:

```bash
cd lean
lake exe cache get
lake build
lake env lean AxiomCheck.lean
```

| Lean name | statement |
|---|---|
| `shift_a`, `shift_b` | `a_n = -α(n+1)`, `b_n = β(n)` |
| `alphaC_one … b21_one` | the displayed values 42, 396, 1047, 38400 |
| `cfP`, `cfQ` | the classical convergent recursions |
| `cfP_neg`, `cfQ_neg` | `P̃_k = (-1)^k P_k`, `Q̃_k = -(-1)^k Q_k` (shifted indexing) |
| `cf_neg_convergent` | `P̃_k/Q̃_k = -(P_k/Q_k)` **for every k**, no hypotheses |
| `challenge_convergent_eq` | the challenge's convergents are minus Cohen's |
| `problem21_pcf_value` | the main theorem, given Cohen's entry |

**Axiom audit.** `#print axioms` (see `AxiomCheck.lean`) reports only
`propext, Classical.choice, Quot.sound` for every theorem above. No
`native_decide`, no `sorry`.

**The one hypothesis that remains.** Entry 5.3.22 is an explicit hypothesis:

```lean
theorem problem21_pcf_value
    (hCohen : Tendsto (fun k => cfP cohenC cohenD k / cfQ cohenC cohenD k)
                atTop (𝓝 (6 / (Real.pi - 3)))) :
    Tendsto (fun k => cfP a21 b21 k / cfQ a21 b21 k)
      atTop (𝓝 (6 / (3 - Real.pi)))
```

Note that `a21`, `b21` are the challenge's own coefficients and `cfP`/`cfQ` the
classical convergent recursion, so the conclusion is a statement about the
challenge's actual continued fraction, not about existentially quantified
sequences.

## Reproducing the numerics

```bash
python3 verify.py
```

`mpmath` is optional (only the two limit checks need it). The script uses the
*same* recursion as the Lean file, so it doubles as a transcription check on the
formalization. Output at n = 61: Cohen's fraction agrees with π to `5.4e-129`,
and the challenge's agrees with `6/(3-π) = -42.3750798355862746187580309154…`
to `1.6e-126`, i.e. the golden-ratio rate of `10 log₁₀ φ = 2.0899…` digits per
term that Cohen records for this entry.

## References

- H. Cohen, *Continued Fractions of Polynomial Type: Theory and Encyclopedic
  Dictionary*, arXiv:2607.06581, 2026. Entry 5.3.22.
- H. Cohen, *A Database of Continued Fractions of Polynomial Type*,
  arXiv:2409.06086, 2024.
- W. B. Jones and W. J. Thron, *Continued Fractions: Analytic Theory and
  Applications*, Addison–Wesley, 1980.
