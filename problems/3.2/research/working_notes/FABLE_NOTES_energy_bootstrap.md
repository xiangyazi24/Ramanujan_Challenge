# Fable: the non-autonomy bootstrap — candidate unconditional E(p) ≪ p^{7/4} (2026-07-31, ChatGPT campaign R1–R2 + solo)

Scripts: `fable_vertical_value_law.py` (value law), `/tmp/p32/sigma_test.py`,
`/tmp/p32/cd_test.py` (to be copied into research/scripts once stable).
Status: proof skeleton with ALL analytic identities machine-verified; remaining
gaps are enumerated in §6 (counting details + exceptional-lag handling). This is
the strongest unconditional track of the campaign; it lands exactly at the
"Heilbronn-precedent cap E ≪ p^{2−δ}" that Q3.2_density_theorem.md predicted but
did not prove.

## 0. Notation

b_n = Apéry ζ(3) numbers; prime p; all congruences mod p. Recurrence
(r+1)³b_{r+1} = P(r)b_r − r³b_{r−1}, P(r) = 34r³+51r²+27r+5, regular for
0 < r < p−1. Transfer: b_{r+d} = A_d(r)b_r + B_d(r)b_{r−1}, with A_d,B_d ∈ F_p(r),
deg ≤ 4d, poles only at r ∈ {−1,…,−d}. N(c) = #{r < p: b_r ≡ c},
E(p) = Σ_c N(c)² (empirically 3p + O(√p)).

## 1. Theorem (proved, 3 lines): image size ≥ ((p−1)/4)^{1/3}

States (b_r, b_{r−1}) take ≤ D² values (D = image size), so some pair (x,y)
repeats ≥ (p−1)/D² times; for those r, b_{r+1} = (P(r)x − r³y)/(r+1)³ is a fixed
rational function of degree ≤ 4, nonconstant unless (x,y) = (0,0) (constancy
forces 51x = 27x from the r²,r¹ coefficients, so x = 0, then y = 0); it takes
each value ≤ 4 times and lands in the image: (p−1)/D² ≤ 4D. ∎
(First image bound of any kind for this sequence, to our knowledge. The same
argument gives image ≥ (c·p)^{1/(k+1)} for any nondegenerate P-recursive sequence
of order k mod p.)

## 2. The collision mechanism: additive characters eliminated

For a collision b_{r+d} = b_r = c:

- **Type I** (B_d(r) ≠ 0): the predecessor is DETERMINED:
  b_{r−1} = c·σ_d(r), σ_d := (1−A_d)/B_d.
  Verified numerically: p = 101, all 36 type-I collision pairs match.
- **Type II** (B_d(r) = 0): then A_d(r) = 1 (for c ≠ 0): r is a root of the
  explicit polynomial system {B_d = 0, A_d = 1} of degree ≤ 4d.
  Empirically the ~50 reflection pairs at p=101 are EXACTLY the type-II pairs:
  the reflection FE is transfer-degenerate (solution-independent), the
  accidental collisions are type-I. (49 type-II + 36 type-I at p = 101.)

Consequence: E(p) − (diagonal) involves NO additive characters — it is the
counting problem "how often does the orbit hit explicit low-degree conditions".

## 3. Two-lag elimination

If r has two collision lags d < d′ (both type I), then σ_d(r) = σ_{d′}(r), i.e.
r is a root of Ψ_{d,d′} := (1−A_d)B_{d′} − (1−A_{d′})B_d, a rational function of
degree ≤ C(d+d′). Together with the type-II polynomials this makes every
multi-collision index r the root of an explicit polynomial of degree O(d+d′) in
the lags involved.

Composition identities (exact, used below):
- T_{d′}(r) = T_{d′−d}(r+d)·T_d(r); det T_d(r) = r³/(r+d)³.
- A_dB_{d′} − A_{d′}B_d = (r³/(r+d)³)·B_{d′−d}(r+d).

## 4. Nondegeneracy: Ψ_{d,d′} ≢ 0 — Laurent criterion with CLOSED-FORM constant

At r = −d, B_d has a pole of exact order 3 with Laurent coefficient

  c_d := lim (r+d)³B_d(r) = −5·B_{d−1}(−d) + B_{d−2}(−d) = **d³·b_{d−1}**

(closed form DISCOVERED and VERIFIED for d = 2..30, exact rational arithmetic;
c_2 = 8·5, c_3 = 27·73, c_4 = 64·1445, …). Since B_{d′}(r) ~ A_{d′−d}(0)·B_d(r)
= b_{d′−d}·B_d(r) at this pole and the cross term is regular there, the pole-3
coefficient of Ψ_{d,d′} at r = −d equals (b_{d′−d} − 1)·d³·b_{d−1} (up to sign).

**Nondegeneracy Lemma (candidate, criterion proved):** Ψ_{d,d′} ≢ 0 in F_p(r)
whenever p ∤ d, p ∤ b_{d−1}, and b_{d′−d} ≢ 1 (mod p).

Exceptional lags for the counting range d ≤ D ~ p^{1/4}:
- {d: b_{d−1} ≡ 0 mod p} ∩ [1,D] — the small-index zero set (empirically ~1 per prime);
- {δ: b_δ ≡ 1 mod p} ∩ [1,D] — the small-index level-1 set.
Both are value-fiber conditions at SMALL indices, so they are themselves sparse;
elements r whose first three lags all pairwise land in exceptional differences
need a separate (structured) count.

Symbolic check: Ψ_{d,d′} ≢ 0 over Q for all 2 ≤ d < d′ ≤ 10 (sympy, exact).

## 5. The counting (Markov/gap), candidate exponents

Work on the reflection quotient [0,(p−1)/2) (halves N, kills reflection pairs).
For a fiber of size N: ≥ N/4 elements have their two smallest lags d₁ < d₂ with
d₁+d₂ ≤ 10p/N (gap counting). Each such element is a root of Ψ_{d₁,d₂} (type I)
or of {B_{d₁} = 0, A_{d₁} = 1} (type II), polynomials of degree ≤ C(d₁+d₂).
Root-capacity: type II: Σ_{d≤D} 4d ≈ 2D² ⟹ N ≲ p^{2/3} branch;
type I: Σ_{d₁<d₂≤D} C(d₁+d₂) ≈ CD³ with D = 10p/N ⟹ N⁴ ≲ Cp³.

**Candidate Theorem A:** N(c) ≪ p^{3/4}·(log p)^{O(1)} for every c ≠ 0.
**Candidate Theorem B:** E(p) ≪ p^{7/4}·(log p)^{O(1)}.
(c = 0 case: single-condition B_d(r) = 0 counting reproduces the PROVED
|Z_p| ≤ 3p^{2/3} — same method family, consistency check passed.)

## 6. Remaining gaps (honest list — this is not yet a theorem)

1. Ψ degree bound: write deg(numerator Ψ_{d,d′}) ≤ C(d+d′) carefully
   (denominator-cleared; poles at −1..−d′ only; C explicit, ~8).
2. Exceptional lags: the pole-3 criterion fails when p | b_{d−1} or b_{d′−d} ≡ 1.
   Fix: use the first THREE lags (three pairwise differences must ALL be
   exceptional — then r is triply constrained by the exceptional structure), or
   use a second Laurent point (r = −d′ gives an independent criterion with its
   own constant — compute it).
3. A_{d} ≡ 1 ∧ B_d ≡ 0 as IDENTITIES in F_p(r) for some d (the type-II analogue of
   degeneracy — would make a whole lag "invisible"): rule out via the same pole
   argument (B_d ≡ 0 impossible since c_d = d³b_{d−1} ≠ 0 unless p | db_{d−1};
   exceptional d handled as in 2).
4. The b_r = 0 endpoints in type I (need b_r ≠ 0 to divide; zero fiber is the
   c = 0 case, separately fine).
5. Uniformity of the Markov step when the fiber clusters (all elements in a
   short interval): then gaps are SMALL, D is small, capacity is tiny — the
   argument only improves; check the boundary r ∈ {0,1,p−2,p−1} exclusions.
6. c = 0 with |Z_p|: use their existing 3p^{2/3}.

## 7. Verified-dead branches from this campaign (for the obstruction map)

- Gross–Koblitz bounded-Γ_p formula for b_r: REFUTED (R1a; parameters move with
  r; non-hypergeometric ⟹ no bounded Gauss-sum monomial; consistent with
  Q6128's "unit part ≍ p gamma-sums").
- p-adic stationary phase / Boyarsky interpolation / any LOCAL p-adic method:
  dead by the one-sample-per-disc obstruction (r runs over one point per residue
  disc; local analyticity sees nothing). Unifies with the supercongruence
  rigidity (Dwork orders ≤ 3 cancel identically).
- Hidden conformal bilinear invariant of the cocycle (R1c blind-spot): REFUTED
  by computation — M̃ᵀJ(r+1)M̃ = ±r³(r+1)³J(r) has zero solution space for
  deg J ≤ 12 (gauge+determinant argument forces λ = ±det, so this is complete
  up to degree 12); rank-1 degenerate case = hypergeometric right factor of the
  Apéry operator = known impossible (irreducibility).
- Classical Stepanov transfer (T-variable): dead at the coefficient/index bridge
  (R1b) — the ODE lives in T, fibers live in r. THE BOOTSTRAP ABOVE IS THE
  REPLACEMENT: non-autonomy of the recurrence is the derivative-substitute.
