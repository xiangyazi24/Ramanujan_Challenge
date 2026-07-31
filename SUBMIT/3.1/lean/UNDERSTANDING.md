# Problem 3.1 — unconditional formalization map

This document tracks the proof that Lean must ultimately expose.  It is not a
list of claims delegated to the paper: every box below must end as a theorem
whose axiom audit contains only Lean/Mathlib foundations.

## Headline statement

The final theorem must use the exact polynomial printed in the official
problem.  It must characterize:

- `α` as the unique root near `0.349269` of `A(α, sqrt α) = 0`;
- `β` as the unique root near `0.406813` of `A(β, β) = 0`;
- `y` as the positive differentiable A-polynomial branch on `[α, β]` with
  `y' ≤ -2`;
- the official interval integral, with value `4 * π^2 / 85`.

No theorem whose hypotheses already contain a regulator value, a torsion
denominator, a numerical approximation of the desired answer, an auxiliary
shape lift, or the target integral qualifies as the headline theorem.

## Proof route

The submitted paper uses geometric volume, extended Bloch groups,
Merkurjev--Suslin, and numerical rational reconstruction.  The unconditional
Lean proof instead uses the same four-tetrahedron chart but proves the needed
variation and endpoint values directly:

1. Set `X=M²`, `t=1-r²`.  The chart equations are

   ```
   H = (L+X³)t - Xr(L+X) = 0
   E = (t-r)(Xt-r)(t-rX) - Xr⁴t = 0.
   ```

2. Elimination of `r` gives exactly the official A-polynomial
   `aPolynomialX X L`.

3. For

   ```
   T=t, U=r/t, V=r/(tX), W=t/(t-rX),
   ```

   direct logarithmic differentiation gives

   ```
   d(R(T)+R(U)+R(V)+R(W))
     = -(log M dlog L - log L dlog M).
   ```

4. At the two real endpoints, transform the four principal real regulator
   terms to Rogers arguments in `(0,1)`:

   ```
   a = 1 - 1/r²
   b = t/r
   c = tX/r
   d = (t-rX)/t.
   ```

   Their exact relations are

   ```
   a² = b²(1-a)
   c = b²(1-d)
   c = bd(1-a)(1-b)(1-c)
   X = c/b
   L = c²(1-c)/(b³d).
   ```

   The endpoint constraints are `X=L⁴` at `α` and `X=L²` at `β`.

5. The remaining endpoint identities are

   ```
   R(aα)+R(bα)+R(cα)+R(dα) = π²/10
   R(aβ)+R(bβ)+R(cβ)+R(dβ) = 5π²/34.
   ```

   Their difference is `4π²/85`.

## Checked atoms

- [x] The exact official polynomial is defined in `APolynomial.lean`.
- [x] The cleared holonomy and edge equations are defined there.
- [x] The four-shape Rogers differential is checked by `ring` in
      `RegulatorDifferential.lean`, including the sign and the factor from
      `X=M²`.
- [x] On `E=H=0`, `0<X<1`, `r<-1`, the linear-subresultant coefficient is
      strictly negative.  `SubresultantSign.lean` proves this from two exact
      polynomial certificates.
- [x] The transformed endpoint relations above have been independently checked
      in exact endpoint number fields by Sage.

## Open atoms

- [x] `Elimination.lean` checks a denominator-free forward certificate
      `H=E=0 → A=0`.
- [x] `ReverseElimination.lean` checks exact linear-subresultant,
      pseudo-division, and Bézout certificates.  Thus every positive point with
      `0 < X < 1` on `A=0` has an explicitly reconstructed chart root satisfying
      both `H=0` and `E=0`, with no nonvanishing hypothesis.
- [ ] Prove rational isolating intervals and uniqueness for `α` and `β`.
- [x] `Dilog/Rogers.lean` constructs the real Rogers function and proves its
      derivative and Euler complement identity.  `Dilog/FiveTerm.lean` proves
      Abel's five-term identity on the open unit square.
- [x] `Dilog/SpecialValues.lean` derives `R(ρ)=π²/10` for
      `ρ²=1-ρ`, and verifies the explicit positive golden-ratio conjugate.
- [x] Exact number-field computation shows that the alpha endpoint field
      already contains `√5`; at the distinguished root,
      `ρ = -q¹¹+2q¹⁰-2q⁹+3q⁸-3q⁷+5q⁶-4q⁵+4q⁴-3q³+3q²-2q+1`.
- [x] `Dilog/EndpointTransform.lean` replays the 22-relation certificate
      reducing the beta chart sum to the five standard `π/17` Rogers
      arguments.  Its remaining hypotheses are explicit algebraic
      cross-identifications and open-unit inequalities, not regulator values.
- [ ] Produce a finite, machine-readable five-term/Y-system certificate for
      both endpoint identities.
- [x] `ShapeChamber.lean` proves all four real Rogers arguments lie in
      `(0,1)` from the negative chart root, `0<X<1`, `X≤L`, and the exact
      holonomy equation.  Along the official decreasing branch, `X≤L`
      follows from `y(x)≥x`; selecting the negative reconstructed root remains.
- [ ] Apply the fundamental theorem of calculus and identify the interval
      integral with the endpoint difference.
- [ ] State the exact official theorem and audit its assumptions and axioms.

## Non-solutions

The legacy `MainTheorem.regulator_value` consumes `torsion` and `numeric`
hypotheses.  It remains useful historical work but cannot be used in the
unconditional headline theorem.  Likewise, the old
`lean/RamanujanChallenge/Problem31.lean` proves a value chosen to be the target
and must not be imported.
