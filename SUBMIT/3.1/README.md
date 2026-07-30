# Ramanujan Challenge — Problem 3.1

**Claim.** Along the arc of the `7₂` A-polynomial curve between the two endpoints
identified in `solution.pdf`,

```
∫ (log x · dy/y − log y · dx/x) = 4π²/85.
```

Contents:

| path | what it is |
|---|---|
| `solution.pdf` / `solution.tex` | the write-up (5 pages) |
| `lean/` | a self-contained Lean 4 project checking the algebraic spine |
| `scripts/` | 15 Sage/Python scripts reproducing every number in the paper |

---

## The proof in one column

1. **Reduction.** Khoi's variation formula: `I = (GV(ρ_α) − GV(ρ_β))/4`.
2. **Beta endpoint, exact.** `S³₋₁(7₂) = Σ(2,3,17)`, maximal Fuchsian;
   Brooks–Goldman gives `GV(ρ_β) = 242π²/51`.  The claim becomes
   `GV(ρ_α) = 74π²/15`.
3. **Neumann–Zagier.** `I = −ΔR`, with `R` the Rogers dilogarithm sum over the
   four tetrahedron shapes.
4. **Rationality (the new ingredient).** Both endpoint minimal polynomials are
   palindromic; the trace polynomials are totally real with all but one/two roots
   in `[−2,2]`; hence every non-real embedding puts the eigenvalue on the unit
   circle, the chart parameter `u` is real, two shapes are real and the other two
   cancel — so the Borel regulator vanishes at **every** embedding, both extended
   Bloch classes are **torsion**, and `Re[ΔR] ∈ π²ℚ`.
5. **Denominator bound.** Merkurjev–Suslin: `|K₃^ind(F)_tors| = w₂(F)`.
   `w₂(F_α) = 120`, `w₂(F_β) = 408`, `Q = lcm = 2040`.
6. **Pin it.** `Re[ΔR]/π² = −4/85` to 301 digits; two distinct rationals of
   denominator `≤ 2040` differ by at least `1/2040² ≈ 2.4·10⁻⁷`.  Note `85 ∣ 2040`.

---

## Lean

`lean/` is a standalone Lake project.

```
cd lean
lake exe cache get
lake build
```

* **Toolchain** `leanprover/lean4:v4.29.0`; **Mathlib** pinned to `v4.29.0`.
* Verified build: `Build completed successfully (3295 jobs)`, zero errors.
* **No `sorry`. No `native_decide`.**  All 37 audited declarations depend only on
  `[propext, Classical.choice, Quot.sound]` — the four `ShapeField` certificates
  need only `[propext, Quot.sound]`.

### Module map

| module | content |
|---|---|
| `Ramanujan31/RatReconstruct.lean` | rational reconstruction from a denominator bound — step 6 |
| `Ramanujan31/ChartSymmetry.lean` | `u(1/a) = u(a)`; the palindromic decompositions `f = a^d·g(a+1/a)` for both endpoints |
| `Ramanujan31/ShapeField.lean` | algebraic certificate that `1 + 4u²` is a **square in the endpoint field** — so the shapes are defined over `F`, not a quadratic extension |
| `Ramanujan31/TraceRoots.lean` | the trace polynomials are totally real, with **exact** root counts in `(−2,2)` |
| `Ramanujan31/UnitCircle.lean` | real trace in `[−2,2]` **⟺** `‖a‖ = 1` |
| `Ramanujan31/ShapeCancellation.lean` | the Bloch–Wigner four-shape cancellation |
| `Ramanujan31/MainTheorem.lean` | the composition, plus the final step |

### Headline theorems

`MainTheorem.four_shape_sum_vanishes_of_trace_real`
: at any embedding whose trace is real with `|w| ≤ 2`, the Bloch–Wigner sum of
  the four shapes vanishes.  This is step 4 at a single embedding, assembled from
  the four modules above.

`MainTheorem.regulator_value`
: given the torsion denominator bound and the numerical certificate, the value
  **is** `−4/85`.  This is step 6.

`TraceRoots.gAlpha_totally_real` / `gBeta_totally_real`
: `gα` splits with five roots in `(−2,2)` and one in `(2, 5/2)`; `gβ` splits with
  six in `(−2,2)` and two in `(2,3)`.  Each statement carries an **exhaustion
  clause** (`∀ x, eval x = 0 → x = r₁ ∨ …`), so the counts are exact, not lower
  bounds — that is what rules out a non-real embedding escaping the unit circle.

### Two referee objections, answered inside Lean

**"The shapes involve `√(1+4u²)`, so they live in a quadratic extension of `F`;
your torsion argument only controls the embeddings of `F`."**
`ShapeField.sPolyAlpha_sq` / `sPolyBeta_sq` exhibit an explicit integer
polynomial `s` with `s(a)² = 1 + 4u(a)²` modulo the minimal polynomial — proved
by `linear_combination c(a)·hf` with an explicit cofactor, over an arbitrary
commutative ring.  So `√(1+4u²) ∈ F`, the shape field `E` equals `F`, and the
statement transports under every embedding at once.  This replaces the earlier
numerical `is_square()` check with a certificate a reader can verify by one
polynomial multiplication.

**"The passage `torsion of order m ⟹ Rogers value in (1/m)π²ℤ` carries
convention-dependent factors 6, 12, 24; your `Q = 2040` may be wrong."**
`regulator_quotient_eq_robust` shows the conclusion holds for **any** bound
`85 ≤ Q ≤ 10¹⁵⁰`.  The numerical certificate has error `1.1·10⁻³⁰¹` while the
separation of rationals of denominator `≤ Q` is `1/Q²`, so there are ~146 orders
of magnitude of slack.  A stray factor of 6 or 24 — or of `10¹⁴⁶` — changes
nothing.  The normalization question therefore does not have to be settled for
the conclusion to stand.

### What is *not* machine-checked

Stated as explicit hypotheses in Lean, and cited in the write-up:

1. **The Bloch–Wigner functional equations**, axiomatized as the structure
   `BlochWignerLaws` (`D(conj z) = −D z`, `D((1−z)⁻¹) = D z`, `D(x : ℝ) = 0`).
   Standard — Zagier, *The dilogarithm function*, §I.2.  Constructing `D` inside
   Lean is a separate infrastructure project.
2. **Torsion ⇒ rational with denominator dividing 2040** — Merkurjev–Suslin plus
   Zickert Thm 1.1.  Enters `regulator_value` as the hypothesis `torsion`.  See
   the robustness theorem above: only the *existence* of a bound below `10¹⁵⁰`
   matters, not its exact value.
3. **The 301-digit numerical evaluation.**  Certified evaluation of Rogers
   dilogarithms to that precision inside Lean is not attempted; the bound enters
   `regulator_value` as the hypothesis `numeric`.  It is reproduced by
   `scripts/`.

Everything else in steps 4 and 6 — the palindromic decompositions, the
totally-real root counts with exact multiplicities, the passage from a real trace
to the unit circle, the reality of the chart parameter, the shape cancellation,
and the rational reconstruction — is machine-checked.
