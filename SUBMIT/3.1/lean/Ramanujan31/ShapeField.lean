import Mathlib.Tactic
import Ramanujan31.ChartSymmetry

/-!
# The shapes are defined over the endpoint field

A referee objection to the write-up: the tetrahedron shapes of the chart are
rational functions of `u` **and `√(1 + 4u²)`**, so a priori they live in a
quadratic extension of the endpoint field `F`, not in `F` itself.  Since the
torsion argument runs over the embeddings of the field the class is defined
over, a genuine quadratic extension would break it — the argument would only
control the embeddings of `F`, not the extra complex places above them.

The objection is answered by an **algebraic certificate**, not a numerical
check: `1 + 4u²` is the square of an explicit element of `F`.

For each endpoint we give two integer polynomials, `uPoly` and `sPoly`, and
prove — over an arbitrary commutative ring, assuming only `f a = 0` — that

* `uPoly a * denominator = numerator`, i.e. `uPoly` **is** the chart parameter;
* `sPoly a ^ 2 = 1 + 4 * uPoly a ^ 2`, i.e. `√(1 + 4u²) = sPoly a ∈ F`.

Both are `linear_combination c(a) * hf` for an explicit cofactor `c`, so the
verification is one polynomial multiplication — checkable by hand, and here by
the kernel.

Because `sPoly a` lies in `F` itself, the statement transports under **every**
embedding of `F` at once; there is no place-by-place subtlety, and the shape
field `E` equals `F`.  That is exactly what the torsion argument needs.

The certificates were found in Sage (`scripts/`), but nothing here depends on
Sage: the identities are re-derived from scratch by `linear_combination`.
-/

namespace ShapeField

section Alpha

variable {K : Type*} [CommRing K]

/-- The degree-12 palindromic minimal polynomial of the alpha eigenvalue. -/
def fAlphaPoly (a : K) : K :=
  a ^ 12 - 3 * a ^ 11 + 4 * a ^ 10 - 5 * a ^ 9 + 6 * a ^ 8 - 7 * a ^ 7
    + 7 * a ^ 6 - 7 * a ^ 5 + 6 * a ^ 4 - 5 * a ^ 3 + 4 * a ^ 2 - 3 * a + 1

/-- The chart parameter `u`, written as a **polynomial** in `a` (valid modulo
`fAlphaPoly`). -/
def uPolyAlpha (a : K) : K :=
  -6 * a ^ 11 + 14 * a ^ 10 - 15 * a ^ 9 + 21 * a ^ 8 - 23 * a ^ 7 + 28 * a ^ 6
    - 25 * a ^ 5 + 27 * a ^ 4 - 20 * a ^ 3 + 18 * a ^ 2 - 13 * a + 10

/-- The certificate: `sPolyAlpha a` is a square root of `1 + 4u²` **inside `F`**. -/
def sPolyAlpha (a : K) : K :=
  12 * a ^ 11 - 30 * a ^ 10 + 32 * a ^ 9 - 42 * a ^ 8 + 50 * a ^ 7 - 58 * a ^ 6
    + 54 * a ^ 5 - 56 * a ^ 4 + 42 * a ^ 3 - 38 * a ^ 2 + 28 * a - 21

/-- **The polynomial form of the chart is correct.**

`u = (a + a¹²)/(a⁵ + a⁸)` — this is `chartUAlpha` after clearing the common
factor `a` — and modulo the minimal polynomial that rational function equals the
polynomial `uPolyAlpha`. -/
theorem uPolyAlpha_spec (a : K) (hf : fAlphaPoly a = 0) :
    uPolyAlpha a * (a ^ 5 + a ^ 8) = a + a ^ 12 := by
  unfold uPolyAlpha
  unfold fAlphaPoly at hf
  linear_combination
    (-6 * a ^ 7 - 4 * a ^ 6 - 3 * a ^ 5 - 8 * a ^ 4 - 5 * a ^ 3 - 3 * a ^ 2 - a) * hf

/-- **The square-root certificate for the alpha endpoint.**

`1 + 4u²` is a square in `F = ℚ(a)`, with explicit square root `sPolyAlpha a`.
Hence the tetrahedron shapes lie in `F`, not in a quadratic extension. -/
theorem sPolyAlpha_sq (a : K) (hf : fAlphaPoly a = 0) :
    sPolyAlpha a ^ 2 = 1 + 4 * uPolyAlpha a ^ 2 := by
  unfold sPolyAlpha uPolyAlpha
  unfold fAlphaPoly at hf
  linear_combination
    (-48 * a ^ 9 + 20 * a ^ 8 + 12 * a ^ 7 + 104 * a ^ 6 + 12 * a ^ 5
      + 104 * a ^ 4 + 12 * a ^ 3 + 56 * a ^ 2 - 16 * a + 40) * hf

end Alpha

section Beta

variable {K : Type*} [CommRing K]

/-- The degree-16 palindromic minimal polynomial of the beta eigenvalue. -/
def fBetaPoly (b : K) : K :=
  b ^ 16 - 7 * b ^ 15 + 22 * b ^ 14 - 48 * b ^ 13 + 87 * b ^ 12 - 133 * b ^ 11
    + 178 * b ^ 10 - 211 * b ^ 9 + 223 * b ^ 8 - 211 * b ^ 7 + 178 * b ^ 6
    - 133 * b ^ 5 + 87 * b ^ 4 - 48 * b ^ 3 + 22 * b ^ 2 - 7 * b + 1

/-- The beta chart parameter as a polynomial in `b`. -/
def uPolyBeta (b : K) : K :=
  -6 * b ^ 15 + 41 * b ^ 14 - 125 * b ^ 13 + 266 * b ^ 12 - 474 * b ^ 11
    + 711 * b ^ 10 - 935 * b ^ 9 + 1088 * b ^ 8 - 1127 * b ^ 7 + 1043 * b ^ 6
    - 857 * b ^ 5 + 620 * b ^ 4 - 389 * b ^ 3 + 202 * b ^ 2 - 85 * b + 21

/-- The beta square-root certificate. -/
def sPolyBeta (b : K) : K :=
  60 * b ^ 15 - 394 * b ^ 14 + 1150 * b ^ 13 - 2386 * b ^ 12 + 4196 * b ^ 11
    - 6176 * b ^ 10 + 8022 * b ^ 9 - 9202 * b ^ 8 + 9406 * b ^ 7 - 8592 * b ^ 6
    + 6958 * b ^ 5 - 4960 * b ^ 4 + 3064 * b ^ 3 - 1546 * b ^ 2 + 646 * b - 137

/-- **The polynomial form of the beta chart is correct.** -/
theorem uPolyBeta_spec (b : K) (hf : fBetaPoly b = 0) :
    uPolyBeta b * (b ^ 3 + b ^ 4) = b + b ^ 6 := by
  unfold uPolyBeta
  unfold fBetaPoly at hf
  linear_combination (-6 * b ^ 3 - 7 * b ^ 2 - b) * hf

/-- **The square-root certificate for the beta endpoint.** -/
theorem sPolyBeta_sq (b : K) (hf : fBetaPoly b = 0) :
    sPolyBeta b ^ 2 = 1 + 4 * uPolyBeta b ^ 2 := by
  unfold sPolyBeta uPolyBeta
  unfold fBetaPoly at hf
  linear_combination
    (3456 * b ^ 14 - 21120 * b ^ 13 + 56640 * b ^ 12 - 111744 * b ^ 11
      + 190968 * b ^ 10 - 268816 * b ^ 9 + 339484 * b ^ 8 - 374884 * b ^ 7
      + 370012 * b ^ 6 - 325084 * b ^ 5 + 251276 * b ^ 4 - 169884 * b ^ 3
      + 98124 * b ^ 2 - 43696 * b + 17004) * hf

end Beta

/-! ## Tying the polynomial chart to the rational chart -/

section Field

variable {K : Type*} [Field K]

/-- Over a field, the polynomial chart agrees with `chartUAlpha` wherever the
latter is defined. -/
theorem chartUAlpha_eq_uPoly (a : K) (ha : a ≠ 0) (h3 : 1 + a ^ 3 ≠ 0)
    (hf : fAlphaPoly a = 0) :
    chartUAlpha a = uPolyAlpha a := by
  have hden : a ^ 4 * (1 + a ^ 3) ≠ 0 := mul_ne_zero (pow_ne_zero _ ha) h3
  have hspec := uPolyAlpha_spec a hf
  unfold chartUAlpha
  rw [div_eq_iff hden]
  refine mul_left_cancel₀ ha ?_
  linear_combination -hspec

/-- Same for the beta chart.  Here `1 + b ≠ 0` is not an extra hypothesis: it
follows from the minimal polynomial, since `fBetaPoly (-1) = 1597 ≠ 0` (and
`1597` is precisely the prime appearing in `disc F_β = 17¹⁴ · 1597`). -/
theorem chartUBeta_eq_uPoly [CharZero K] (b : K) (hb : b ≠ 0)
    (hf : fBetaPoly b = 0) :
    chartUBeta b = uPolyBeta b := by
  have hb1 : (1 : K) + b ≠ 0 := by
    intro h
    have hb' : b = -1 := by linear_combination h
    rw [hb'] at hf
    unfold fBetaPoly at hf
    norm_num at hf
  have hspec := uPolyBeta_spec b hf
  unfold chartUBeta
  rw [div_eq_iff (pow_ne_zero 2 hb)]
  refine mul_left_cancel₀ (mul_ne_zero hb hb1) ?_
  linear_combination -hspec

end Field

end ShapeField
