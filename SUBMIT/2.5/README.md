# Ramanujan Challenge — Problem 2.5

**Catalan's constant G as a 3x3 CMF limit**

Submitter: Xiang Huang (University of Illinois Springfield)
Contact: xhuan5@uis.edu

## What is in this package

- solution.pdf / solution.tex — mathematical proof (16 pages)
- independent_miller.tex — independent elementary Miller proof included by
  solution.tex
- lean/ — Lean 4 source snapshot for the independent proof
- README.md — this file

## Summary

We prove that all three column ratios of the 3x3 CMF in Problem 2.5
converge to Catalan's constant G. The paper now records two exact routes:

1. A Meijer-G inverse-transpose trajectory, reconstructed after release of
   the official solutions.
2. A mathematically independent route using only the printed matrix,
   elementary positive triple integrals, rational polynomial certificates,
   and a contracting Miller cone.

The second route constructs a fast positive adjoint solution, proves its
matrix recurrence by exact creative telescoping, selects it by a uniform
backward contraction, and connects it to the literal third column by an
exact rational terminal certificate. A separate positive-cocycle squeeze
then gives the same limit for all three columns.

## Lean formalization status

The canonical development is under `../../lean/RamanujanChallenge/`.
The `lean/` directory here mirrors the 17-file dependency closure of the
final proof together with four earlier Delannoy auxiliary modules
(21 files, 10,934 lines). Main theorem:
`problem25_solved : Problem25Claim`.

### Fully proved

- Catalan integral identity
- positivity and 2/3 projective contraction for all three ratios
- triple-moment integrability, positivity, and the 64^(-n) bound
- exact initial denominator and Catalan pairings
- kernel-checked sparse divergence certificates for the adjoint recurrence
- rational companion identities and coefficient bounds
- invariant Miller box and 1/8 contraction
- fast-solution selection by the defect-growth argument
- exact third-column terminal certificate
- `commonLimit = catalanConstant`, normalized error convergence, and the
  final challenge claim

### Axioms

There are 0 `sorry` and 0 added axioms.
The final theorem depends only on
`{propext, Classical.choice, Quot.sound}`.

## Verification

From the repository's `lean/` directory:

```bash
lake env lean RamanujanChallenge/Problem25Moment.lean
```
