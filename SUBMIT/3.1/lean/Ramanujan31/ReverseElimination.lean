import Ramanujan31.Elimination
import Ramanujan31.SubresultantSign

/-!
# Exact reverse elimination on the nonsingular A-polynomial locus

The forward elimination theorem shows that every point of the shape chart lies
on the official A-polynomial.  This file records the converse away from the
vanishing locus of the linear-subresultant coefficient.

All identities below are displayed polynomial identities checked by `ring`.
The computer algebra calculation used to discover them is not trusted by Lean.
-/

/-- Constant coefficient of the normalized linear subresultant of the
holonomy and edge equations. -/
def subresultantB (X L : ℝ) : ℝ :=
  -X * (X + L) * (X ^ 3 + L)
    * (-X ^ 4 * L + 2 * X ^ 4 + 3 * X ^ 3 * L + X ^ 2 * L ^ 2
      + 2 * X ^ 2 * L + X ^ 2 + 3 * X * L + 2 * L ^ 2 - L)

/-- The rational root selected by the linear subresultant. -/
noncomputable def reconstructedR (X L : ℝ) : ℝ :=
  -subresultantB X L / subresultantA X L

private theorem holonomy_clear_denominator
    (X L a b : ℝ) (ha : a ≠ 0) :
    a ^ 2 * holonomyEquation X L (-b / a)
      = (L + X ^ 3) * (a ^ 2 - b ^ 2) + X * b * a * (L + X) := by
  unfold holonomyEquation chartT
  field_simp [ha]
  ring

/-- Clearing the denominator in the reconstructed holonomy equation gives
exactly the official A-polynomial. -/
theorem reconstructedR_holonomy_identity
    (X L : ℝ) (ha : subresultantA X L ≠ 0) :
    subresultantA X L ^ 2
        * holonomyEquation X L (reconstructedR X L)
      = (X ^ 3 + L) ^ 4 * aPolynomialX X L := by
  calc
    subresultantA X L ^ 2
          * holonomyEquation X L (reconstructedR X L)
        =
          (L + X ^ 3)
              * (subresultantA X L ^ 2 - subresultantB X L ^ 2)
            + X * subresultantB X L * subresultantA X L * (L + X) := by
              simpa [reconstructedR] using
                holonomy_clear_denominator X L
                  (subresultantA X L) (subresultantB X L) ha
    _ = (X ^ 3 + L) ^ 4 * aPolynomialX X L := by
      unfold subresultantB subresultantA subresultantA4 subresultantA3
        subresultantA2 subresultantA1 subresultantA0 aPolynomialX
      ring

private def bezoutP19 (X : ℝ) : ℝ :=
  X ^ 19 + 3 * X ^ 18 - 23 * X ^ 16 - 71 * X ^ 15 - 109 * X ^ 14
    - 82 * X ^ 13 + 55 * X ^ 12 + 192 * X ^ 11 + 321 * X ^ 10
    + 205 * X ^ 9 + 201 * X ^ 8 - 84 * X ^ 7 + 48 * X ^ 6
    - 183 * X ^ 5 + 68 * X ^ 4 - 77 * X ^ 3 + 64 * X ^ 2
    - 20 * X + 2

private def bezoutU0 (X : ℝ) : ℝ :=
  -X ^ 4
    * (X ^ 25 + X ^ 24 - 8 * X ^ 23 - 38 * X ^ 22 - 65 * X ^ 21
      + 4 * X ^ 20 + 289 * X ^ 19 + 746 * X ^ 18 + 1013 * X ^ 17
      + 625 * X ^ 16 - 539 * X ^ 15 - 1715 * X ^ 14 - 2231 * X ^ 13
      - 1337 * X ^ 12 - 523 * X ^ 11 + 667 * X ^ 10 + 12 * X ^ 9
      + 632 * X ^ 8 - 448 * X ^ 7 + 577 * X ^ 6 - 388 * X ^ 5
      + 319 * X ^ 4 - 238 * X ^ 3 + 112 * X ^ 2 - 25 * X + 2)

private def bezoutU1 (X : ℝ) : ℝ :=
  X ^ 3
    * (X ^ 24 - X ^ 23 + X ^ 22 + 29 * X ^ 21 + 81 * X ^ 20
      - 25 * X ^ 19 - 620 * X ^ 18 - 1695 * X ^ 17 - 2475 * X ^ 16
      - 1508 * X ^ 15 + 1140 * X ^ 14 + 4678 * X ^ 13
      + 5462 * X ^ 12 + 5080 * X ^ 11 + 903 * X ^ 10 + 415 * X ^ 9
      - 3017 * X ^ 8 + 47 * X ^ 7 - 2033 * X ^ 6 + 1500 * X ^ 5
      - 734 * X ^ 4 + 818 * X ^ 3 - 485 * X ^ 2 + 112 * X - 9)

private def bezoutU2 (X : ℝ) : ℝ :=
  -(X ^ 25 - 4 * X ^ 24 - 23 * X ^ 23 - 42 * X ^ 22 + 40 * X ^ 21
    + 366 * X ^ 20 + 939 * X ^ 19 + 1449 * X ^ 18 + 1194 * X ^ 17
    + 54 * X ^ 16 - 2206 * X ^ 15 - 3368 * X ^ 14 - 4634 * X ^ 13
    - 2206 * X ^ 12 - 2319 * X ^ 11 + 1661 * X ^ 10 - 484 * X ^ 9
    + 2391 * X ^ 8 - 894 * X ^ 7 + 1046 * X ^ 6 - 848 * X ^ 5
    + 336 * X ^ 4 - 192 * X ^ 3 + 100 * X ^ 2 - 24 * X + 2)

private def bezoutU3 (X : ℝ) : ℝ :=
  (X ^ 4 + 3 * X ^ 2 + 1) * bezoutP19 X

/-- First multiplier in a Bézout certificate for the official A-polynomial
and the linear-subresultant coefficient, viewed as polynomials in `L`. -/
def subresultantBezoutU (X L : ℝ) : ℝ :=
  bezoutU3 X * L ^ 3 + bezoutU2 X * L ^ 2 + bezoutU1 X * L + bezoutU0 X

private def bezoutV0 (X : ℝ) : ℝ :=
  X ^ 6
    * (X ^ 24 + X ^ 23 - 10 * X ^ 22 - 36 * X ^ 21 - 29 * X ^ 20
      + 101 * X ^ 19 + 339 * X ^ 18 + 430 * X ^ 17 + 74 * X ^ 16
      - 584 * X ^ 15 - 955 * X ^ 14 - 530 * X ^ 13 + 175 * X ^ 12
      + 601 * X ^ 11 + 137 * X ^ 10 - 116 * X ^ 9 - 344 * X ^ 8
      + 143 * X ^ 7 + 62 * X ^ 6 + 127 * X ^ 5 - 85 * X ^ 4
      - 23 * X ^ 3 + 4 * X ^ 2 + 5 * X + 1)

private def bezoutV1 (X : ℝ) : ℝ :=
  -X ^ 2
    * (X ^ 28 - X ^ 27 - 12 * X ^ 26 - 17 * X ^ 25 + 41 * X ^ 24
      + 166 * X ^ 23 + 162 * X ^ 22 - 259 * X ^ 21 - 1016 * X ^ 20
      - 1329 * X ^ 19 - 567 * X ^ 18 + 1210 * X ^ 17
      + 2288 * X ^ 16 + 2447 * X ^ 15 + 720 * X ^ 14 - 137 * X ^ 13
      - 1504 * X ^ 12 - 666 * X ^ 11 - 428 * X ^ 10 + 724 * X ^ 9
      + 510 * X ^ 8 + 156 * X ^ 7 - 410 * X ^ 6 - 72 * X ^ 5
      - 127 * X ^ 4 + 271 * X ^ 3 - 132 * X ^ 2 + 27 * X - 2)

private def bezoutV2 (X : ℝ) : ℝ :=
  2 * X ^ 28 + X ^ 27 - 17 * X ^ 26 - 49 * X ^ 25 - 23 * X ^ 24
    + 142 * X ^ 23 + 341 * X ^ 22 + 330 * X ^ 21 - 47 * X ^ 20
    - 195 * X ^ 19 - 79 * X ^ 18 + 639 * X ^ 17 - 30 * X ^ 16
    - 619 * X ^ 15 - 2355 * X ^ 14 - 1614 * X ^ 13 - 1285 * X ^ 12
    + 383 * X ^ 11 + 543 * X ^ 10 + 893 * X ^ 9 + 207 * X ^ 8
    + 216 * X ^ 7 - 595 * X ^ 6 + 234 * X ^ 5 - 296 * X ^ 4
    + 325 * X ^ 3 - 144 * X ^ 2 + 28 * X - 2

private def bezoutV3 (X : ℝ) : ℝ :=
  -(X ^ 26 + 2 * X ^ 25 - 3 * X ^ 24 - 23 * X ^ 23 - 48 * X ^ 22
    - 35 * X ^ 21 + 40 * X ^ 20 + 148 * X ^ 19 + 67 * X ^ 18
    - 163 * X ^ 17 - 656 * X ^ 16 - 523 * X ^ 15 - 214 * X ^ 14
    + 1027 * X ^ 13 + 1315 * X ^ 12 + 1536 * X ^ 11 + 546 * X ^ 10
    + 61 * X ^ 9 - 597 * X ^ 8 - 253 * X ^ 7 - 505 * X ^ 6
    + 360 * X ^ 5 - 218 * X ^ 4 + 335 * X ^ 3 - 200 * X ^ 2
    + 48 * X - 4)

private def bezoutV4 (X : ℝ) : ℝ :=
  -bezoutP19 X

/-- Second multiplier in the Bézout certificate. -/
def subresultantBezoutV (X L : ℝ) : ℝ :=
  bezoutV4 X * L ^ 4 + bezoutV3 X * L ^ 3 + bezoutV2 X * L ^ 2
    + bezoutV1 X * L + bezoutV0 X

/-- Right-hand side of the Bézout certificate.  Its only real zeros are
`X = -1`, `X = 0`, and `X = 1`. -/
def subresultantBezoutD (X : ℝ) : ℝ :=
  (X - 1) ^ 8 * (X + 1) ^ 8 * X ^ 12 * (X ^ 2 + X + 1) ^ 5

/-- Exact Bézout identity.  It was discovered by the extended Euclidean
algorithm over `ℚ(X)[L]` and is checked here from the displayed coefficients. -/
theorem subresultant_bezout (X L : ℝ) :
    subresultantBezoutU X L * aPolynomialX X L
        + subresultantBezoutV X L * subresultantA X L
      = subresultantBezoutD X := by
  unfold subresultantBezoutU subresultantBezoutV subresultantBezoutD
    bezoutU3 bezoutU2 bezoutU1 bezoutU0 bezoutV4 bezoutV3 bezoutV2 bezoutV1
    bezoutV0 bezoutP19 aPolynomialX subresultantA subresultantA4
    subresultantA3 subresultantA2 subresultantA1 subresultantA0
  ring

/-- The linear-subresultant coefficient cannot vanish on the official
A-polynomial when `0 < X < 1`. -/
theorem subresultantA_ne_zero_of_aPolynomialX
    {X L : ℝ} (hX : 0 < X) (hX1 : X < 1)
    (hA : aPolynomialX X L = 0) :
    subresultantA X L ≠ 0 := by
  intro ha
  have hcert := subresultant_bezout X L
  rw [hA, ha] at hcert
  norm_num at hcert
  have hXm1 : X - 1 ≠ 0 := by linarith
  have hXp1 : X + 1 ≠ 0 := by linarith
  have hX0 : X ≠ 0 := ne_of_gt hX
  have hquad : X ^ 2 + X + 1 ≠ 0 := by nlinarith [sq_nonneg X]
  have hD : subresultantBezoutD X ≠ 0 := by
    unfold subresultantBezoutD
    exact mul_ne_zero
      (mul_ne_zero
        (mul_ne_zero (pow_ne_zero _ hXm1) (pow_ne_zero _ hXp1))
        (pow_ne_zero _ hX0))
      (pow_ne_zero _ hquad)
  exact hD hcert.symm

private def edgePseudoQ0 (X L : ℝ) : ℝ :=
  -X
    * (-X ^ 9 - X ^ 7 * L - X ^ 6 * L ^ 2 + 2 * X ^ 7
      + 2 * X ^ 6 * L + 4 * X ^ 5 * L ^ 2 + X ^ 4 * L ^ 3
      + 2 * X ^ 5 * L + 2 * X ^ 4 * L ^ 2 + X ^ 5 + 4 * X ^ 4 * L
      + 2 * X ^ 3 * L ^ 2 + 2 * X ^ 2 * L ^ 3 - X ^ 3 * L
      - X ^ 2 * L ^ 2 - L ^ 3)

private def edgePseudoQ1 (X L : ℝ) : ℝ :=
  (X - 1) * (X ^ 3 + L)
    * (-X ^ 7 - 2 * X ^ 6 - X ^ 5 * L - X ^ 5 + X ^ 3 * L ^ 2
      - X ^ 4 + X ^ 2 * L ^ 2 + X ^ 2 * L + 2 * X * L ^ 2 + L ^ 2)

private def edgePseudoQ2 (X L : ℝ) : ℝ :=
  -(X ^ 3 + L) ^ 2
    * (-X ^ 5 + 2 * X ^ 4 + X ^ 3 * L + X ^ 2 + 2 * X * L - L)

private def edgePseudoQ3 (X L : ℝ) : ℝ :=
  (X ^ 2 + X + 1) * (X ^ 3 + L) ^ 3

/-- Quotient in the pseudo-division of the degree-five edge equation by the
quadratic holonomy equation. -/
def edgePseudoQuotient (X L r : ℝ) : ℝ :=
  edgePseudoQ3 X L * r ^ 3 + edgePseudoQ2 X L * r ^ 2
    + edgePseudoQ1 X L * r + edgePseudoQ0 X L

/-- Exact pseudo-division certificate.  The remainder is the normalized linear
subresultant multiplied by `-X²`. -/
theorem edge_pseudo_division (X L r : ℝ) :
    (X ^ 3 + L) ^ 4 * edgeEquation X r
      =
        edgePseudoQuotient X L r * holonomyEquation X L r
          - X ^ 2 * (subresultantA X L * r + subresultantB X L) := by
  unfold edgePseudoQuotient edgePseudoQ3 edgePseudoQ2 edgePseudoQ1 edgePseudoQ0
    subresultantB subresultantA subresultantA4 subresultantA3 subresultantA2
    subresultantA1 subresultantA0 edgeEquation holonomyEquation chartT
  ring

/-- The reconstructed root satisfies the normalized linear subresultant. -/
theorem reconstructedR_linear
    (X L : ℝ) (ha : subresultantA X L ≠ 0) :
    subresultantA X L * reconstructedR X L + subresultantB X L = 0 := by
  unfold reconstructedR
  field_simp [ha]
  ring

/-- On the positive nonsingular locus of the official A-polynomial, the
reconstructed root satisfies both exact chart equations. -/
theorem chart_of_aPolynomialX
    {X L : ℝ}
    (hX : 0 < X) (hX1 : X < 1) (hL : 0 < L)
    (hA : aPolynomialX X L = 0) :
    holonomyEquation X L (reconstructedR X L) = 0
      ∧ edgeEquation X (reconstructedR X L) = 0 := by
  have ha := subresultantA_ne_zero_of_aPolynomialX hX hX1 hA
  have hHcert := reconstructedR_holonomy_identity X L ha
  rw [hA, mul_zero] at hHcert
  have ha2 : subresultantA X L ^ 2 ≠ 0 := pow_ne_zero _ ha
  have hH : holonomyEquation X L (reconstructedR X L) = 0 :=
    (mul_eq_zero.mp hHcert).resolve_left ha2
  have hlin := reconstructedR_linear X L ha
  have hEcert := edge_pseudo_division X L (reconstructedR X L)
  have hEprod :
      (X ^ 3 + L) ^ 4 * edgeEquation X (reconstructedR X L) = 0 := by
    simpa [hH, hlin] using hEcert
  have hsum : X ^ 3 + L ≠ 0 := by positivity
  have hpow : (X ^ 3 + L) ^ 4 ≠ 0 := pow_ne_zero _ hsum
  exact ⟨hH, (mul_eq_zero.mp hEprod).resolve_left hpow⟩
