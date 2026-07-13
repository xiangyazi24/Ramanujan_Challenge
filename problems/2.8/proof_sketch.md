# Problem 2.8: Proof Sketch — CMF encoding of the Chudnovsky formula

## Theorem
The 4×4 matrix recurrence of Problem 2.8 encodes the Chudnovsky series for √10005/π.

## Key identifications

### 1. The constant R
R = 151931373056001 = 1 - j(τ₁₆₃)/1728

where j is the modular j-function and τ₁₆₃ = (1+√-163)/2.

Equivalently: R - 1 = 640320³/1728 = 53360³.

### 2. The cubic fingerprint
With u = 2n+3 and k = n+1:
w = u(3u-2)(3u+2) = (2k+1)(6k+1)(6k+5)

This is the ratio of consecutive Chudnovsky hypergeometric terms:
(6k)! / ((3k)!(k!)³) has ratio involving (6k+1)(6k+2)...(6k+6)/((3k+1)(3k+2)(3k+3)(k+1)³)

### 3. The Chudnovsky series
√10005/π = (1/426880) Σ_{k≥0} (-1)^k (6k)!/((3k)!(k!)³) × (545140134k + 13591409) / 640320^{3k}

This is a PROVEN identity (Chudnovsky 1988, Borwein-Borwein 1987) via:
- CM theory of Q(√-163)
- The Ramanujan-type series classification by level
- Singular values of modular functions

### 4. The 4×4 → scalar reduction
The 4×4 matrix product A·M(0)·...·M(N-1) gives a 2×4 matrix.
Each column j gives P_{N,j}/Q_{N,j} → √10005/π.

The scalar recurrence (order ≤ 4) from minors of M(n) should reproduce
the Chudnovsky hypergeometric recurrence after gauge transformation.

## Proof structure

**Lemma 1 (CMF encoding).** The scalar recurrence extracted from M(n)
is equivalent, after rational gauge transformation, to the Clausen-type
recurrence for the Chudnovsky partial sums.

*Proof.* Compute the 4×4 minors, extract the scalar order-4 recurrence,
and show it matches the known Clausen hypergeometric recurrence
(n+1)³ u_{n+1} = ... with parameters determined by R = 1-j(τ₁₆₃)/1728.

**Lemma 2 (Initial conditions).** The initial matrix A selects the
solution corresponding to the Chudnovsky series partial sums.

*Proof.* Verify numerically to 500+ digits that the first few ratios
P_{N,j}/Q_{N,j} match √10005/π.

**Theorem (Main).** lim P_{N,j}/Q_{N,j} = √10005/π for j = 1,2,3,4.

*Proof.* By Lemma 1, the ratios equal the Chudnovsky partial sums
(or a Möbius transform thereof). By the proven Chudnovsky identity
(CM theory of discriminant -163), the limit is √10005/π.

## No-Go: Rational Gauge Does NOT Exist (ChatGPT Q4660)

A 100-digit numerical + symbolic audit proved:
**No G(n) ∈ GL_4(Q(n)) satisfies M(n)·G(n+1) = G(n)·T_{n+1}.**

Reason: det M(n) ~ -c·n^8 but det T → 1/(R-1)^2 (constant).
So det G(n+1)/det G(n) ~ -K/n^8, but for rational g(n), g(n+1)/g(n) → 1.

The gauge must be FACTORIAL or HYPERGEOMETRIC (involving (n!)^k factors).
This is consistent with the Clausen ₃F₂ structure: the natural gauge
includes (6k)!/((3k)!(k!)³) type factors.

## Corrected approach
1. Use FACTORIAL gauge: G(n) = D(n) · G_rat(n) where D(n) is diagonal
   with entries like (n!)^a · (2n)!^b / (3n)!^c
2. Or: extract scalar recurrence from M(n) minors directly (bypassing gauge)
3. Or: show M(n)'s scalar recurrence matches Chudnovsky after gauge normalization
4. The 12-digit numerical match already confirms the identity; the algebraic
   certificate needs the factorial gauge structure

## Implementation plan (revised)
1. Extract scalar order-4 recurrence from M(n) via Casorati minors
2. Apply factorial gauge to normalize
3. Compare with Chudnovsky/Clausen recurrence
4. Cite Chudnovsky CM theorem for the evaluation
