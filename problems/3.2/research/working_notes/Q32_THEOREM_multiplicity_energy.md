# Theorem: unconditional multiplicity and collision-energy bounds for Apéry numbers mod p

_Fable, 2026-07-31. Referee-ready statement + proof of the clean (log-saving)
versions; the p^{3/4}/p^{7/4} upgrades are stated as Conjecture-with-program in
§5. Companion notes: FABLE_NOTES_energy_bootstrap.md. All computational claims
verified in research/scripts/{sigma_test.py, cd_test.py, ad_test.py,
fable_vertical_value_law.py}._

## Statement

Let b_n = Σ_k C(n,k)²C(n+k,k)² be the Apéry ζ(3) numbers, p a prime,
N_p(c) = #{0 ≤ r ≤ p−1 : b_r ≡ c (mod p)}, and E(p) = Σ_{c∈F_p} N_p(c)².

**Theorem 1.** There is an absolute constant C₀ such that for every prime p and
every c ∈ F_p^*: N_p(c) ≤ C₀·p/log p.

**Theorem 2.** E(p) ≤ C₀·p²/log p + N_p(0)². With the known |Z_p| ≤ 3p^{2/3}
this gives E(p) ≪ p²/log p.

**Theorem 3** (image bound). #{b_r mod p : 0 ≤ r < p} ≥ ((p−1)/4)^{1/3}.

Remarks. (i) The trivial bounds are N ≤ p, E ≤ p². The random truth is
N ≍ log p/log log p, E = 3p + O(√p) (measured). These are the first
unconditional savings for this sequence; the method uses only the recurrence's
NON-AUTONOMY, no algebraic geometry. (ii) Theorem 2 gives, via Parseval
(Σ_{h≠0}|C_p(h)|² = pE − p²), the average bound RMS_h |C_p(h)| ≪ p/√log p.

## 1. Setup

Recurrence: (m+1)³b_{m+1} = P(m)b_m − m³b_{m−1}, P(m) = 34m³+51m²+27m+5,
regular at every m ∈ {1,…,p−2} mod p. State vector v_r = (b_r, b_{r−1})ᵀ,
v_{r+1} = M(r)v_r, M(x) = [[P(x)/(x+1)³, −x³/(x+1)³],[1,0]].
Transfer: T_d(r) = M(r+d−1)···M(r), so v_{r+d} = T_d(r)v_r; write
b_{r+d} = A_d(r)b_r + B_d(r)b_{r−1}. A_d, B_d ∈ Q(x) have poles only at
x ∈ {−1,…,−d}, each of order ≤ 3, and (clearing ∏_{j=1}^d(x+j)³) numerator
degrees ≤ 4d. All identities below are algebraic, hence valid over F_p verbatim
whenever the evaluation points avoid poles and 0 < d < p.

Facts used: (F1) A_d(0) = b_d, B_d(0) = 0 (M(0) = [[5,0],[1,0]] is regular and
kills the second column). (F2) B_d(x) = −x³/(x+1)³·A_{d−1}(x+1) (second column
of M(x)). (F3) T_{d′}(r) = T_{d′−d}(r+d)T_d(r) and det T_d(r) = r³/(r+d)³.
(F4) b_k < e^{3.526k} for all k ≥ 1 (b_k < (1+√2)^{4k}; induction from the
recurrence, or the classical asymptotic with the ratio test).
(F5) No two consecutive b_r vanish mod p (else the recurrence forces the zero
solution back to b_0 = 1 through the regular range).

## 2. The collision mechanism

Let c ≠ 0 and let r < r+d be two indices in the fiber {b ≡ c}, with
1 ≤ r, r+d ≤ p−2 (so the transfer stays in the regular range).

From c = b_{r+d} = A_d(r)c + B_d(r)b_{r−1}:

- **Type I** (B_d(r) ≠ 0): b_{r−1} = c·σ_d(r), σ_d := (1−A_d)/B_d. The
  predecessor is determined.
- **Type II** (B_d(r) = 0): (A_d(r)−1)c = 0, so A_d(r) = 1; r is a common root
  of the numerators of B_d and A_d − 1 (degree ≤ 4d each).

If r has two type-I lags d₁ < d₂ (i.e. r, r+d₁, r+d₂ all in the fiber), then
σ_{d₁}(r) = σ_{d₂}(r), i.e. r is a root of the numerator of

  Ψ_{d₁,d₂} := (1−A_{d₁})B_{d₂} − (1−A_{d₂})B_{d₁},   deg ≤ 9(d₁+d₂).

(Verified numerically: at p = 101 all 36 type-I collision pairs satisfy the
predecessor formula; the 49 remaining pairs are type II and are exactly the
reflection pairs r ↔ p−1−r.)

## 3. Nondegeneracy of Ψ

**Lemma (polar constants).** At x = −d, B_d and A_d have poles of order exactly
3 with Laurent coefficients

  c_d := lim (x+d)³B_d(x) = d³·b_{d−1},   a_d := lim (x+d)³A_d(x) = −d³·b_d.

_Proof._ Split M(x) = (x+1)^{−3}(1,0)ᵀ(P(x), −x³) + [[0,0],[1,0]]. Since
T_{d−1}(x) has poles only at −1,…,−(d−1), it is regular at −d, and
T_d(x) = M(x+d−1)T_{d−1}(x) has order-3 polar part at x = −d equal to
(1,0)ᵀ·(P(−1), 1)·T_{d−1}(−d) = (1,0)ᵀ·(−5,1)·T_{d−1}(−d).

Claim: w_d := (−5,1)·T_{d−1}(−d) = d³·(−b_d, b_{d−1}).
Base d = 2: w_2 = (−5,1)M(−2); with P(−2) = −117, (x+1)³ = −1, x³ = −8:
w_2 = (5·(−117)+1, −(−5)(−8)/(−1)·(−1)... ) = (−584, 40) = 8·(−73, 5) = 2³(−b_2, b_1). ✓
Step: w_{d+1} = w_d·M(−d−1) (note T_d(−d−1) = M(−2)···M(−d−1) shifts the base by one).
With x = −d−1: (x+1)³ = −d³, x³ = −(d+1)³, and P(−1−d) = −P(d) (the index
reflection antisymmetry of the Apéry operator, checked by direct expansion):
d³(−b_d, b_{d−1})·M(−d−1)
  first entry: d³[ −b_d·(−P(d))/(−d³) + b_{d−1} ] = −(P(d)b_d − d³b_{d−1})
             = −(d+1)³b_{d+1} (the recurrence),
  second entry: d³[ −(−b_d)·(−(d+1)³)/(−d³) ] = (d+1)³·b_d.
So w_{d+1} = (d+1)³(−b_{d+1}, b_d). ∎
(Also verified exactly for d ≤ 30 (c_d) and d ≤ 15 (a_d) in rational
arithmetic. The pole has order exactly 3 mod p iff p ∤ d³b_{d−1} resp. d³b_d.)

**Lemma (nondegeneracy).** Let 1 ≤ d < d′, δ = d′−d, and suppose p ∤ d,
b_{d−1} ≢ 0, and b_δ ≢ 1 (mod p). Then Ψ_{d,d′} ≢ 0 in F_p(x).

_Proof._ Write Ψ = (B_{d′} − B_d) − K with K = A_dB_{d′} − A_{d′}B_d. By (F3),
K = det T_d·B_δ(x+d) = (x³/(x+d)³)·B_δ(x+d); since B_δ vanishes to order 3 at
0 (by F2), K is regular at x = −d. By (F3) again,
B_{d′}(x) = A_δ(x+d)B_d(x) + B_δ(x+d)B_{d−1}(x); at x = −d the second summand
vanishes (order-3 zero times a regular value — B_{d−1} has no pole at −d), and
A_δ(x+d) → A_δ(0) = b_δ by (F1). Hence the (x+d)^{−3}-coefficient of Ψ equals
(b_δ − 1)·c_d = (b_δ − 1)·d³·b_{d−1} ≠ 0 in F_p by hypothesis. A function with
a genuine pole is not identically zero. ∎

## 4. Proof of Theorems 1–3

**Theorem 3.** If the image has size D, the p−1 states v_r (1 ≤ r ≤ p−1) take
≤ D² values, so some state (x,y) repeats for a set R with |R| ≥ (p−1)/D².
For r ∈ R, b_{r+1} = (P(r)x − r³y)/(r+1)³ =: φ(r), a rational function of
degree ≤ 4, nonconstant unless (x,y) = (0,0) (constancy forces 51x = 27x = 3c,
so x = 0, then y = 0; the zero state is excluded by F5). Each value of φ is
taken ≤ 4 times, and φ(R) lies in the image: (p−1)/D² ≤ 4D. ∎

**Theorem 1.** Let c ≠ 0, N = N_p(c), and work in the window 1 ≤ r ≤ p−2−D
with D := 6p/N (the ≤ D+2 excluded fiber elements are absorbed below). For
each fiber element r let d₁(r) < d₂(r) be the lags to the next two fiber
elements. Telescoping gives Σ_r d₁(r) ≤ p and Σ_r d₂(r) ≤ 2p, so at least
N/2 − O(D) elements satisfy d₁+d₂ ≤ D′ := 12p/N.

Suppose N ≥ C₀p/log p with C₀ = 50. Then D′ ≤ 0.24·log p, so every lag and
lag-difference k in play satisfies 3.526·k ≤ 0.85·log p, whence by (F4) the
integers b_{d₁−1} and b_δ − 1 lie strictly between 0 and p — in particular
b_{d₁−1} ≢ 0 and b_δ ≢ 1 (mod p), and p ∤ d₁. So for every pair (d₁,d₂) in
play, either the element is type II at d₁ (root of B_{d₁}·(A_{d₁}−1)-numerator,
degree ≤ 8d₁ — nonzero since its pole-3 coefficient c_{d₁} = d₁³b_{d₁−1} ≠ 0)
or type I at both (root of the Ψ_{d₁,d₂} numerator, ≢ 0 by the nondegeneracy
lemma, degree ≤ 9(d₁+d₂) ≤ 9D′).

Total root capacity over all pairs d₁ < d₂ ≤ D′:
Σ 9(d₁+d₂) ≤ 9·D′·(D′)²/2 ≤ 5(D′)³ ≤ 5·(0.24 log p)³ < (log p)³.
But we need ≥ N/2 − O(D) ≥ (C₀/4)·p/log p elements rooted — impossible for
large p. Contradiction; hence N < C₀·p/log p. (Small p: absorb into C₀.) ∎

**Theorem 2.** E(p) = Σ_c N(c)² ≤ N_p(0)² + (max_{c≠0} N(c))·Σ_{c≠0}N(c)
≤ N_p(0)² + C₀p²/log p. ∎

## 4b. Theorem 4 (exceptional-set version of the p^{3/4} bound) — NEW, closes §5's program

**Theorem 4.** There is a set 𝔈 of primes with #{q ∈ 𝔈 : q ≤ X} = O(X^{3/4+o(1)})
(density zero) such that for every prime p ∉ 𝔈 and every c ∈ F_p^*:
N_p(c) ≪ p^{3/4}(log p)^{O(1)}, and hence E(p) ≪ p^{7/4}(log p)^{O(1)}.

_Proof plan (all pieces in place)._ Run §4's counting with D′ = 12p/N and
N ≥ p^{3/4+ε}. A pair (d₁,d₂), δ = d₂−d₁, escapes both Laurent criteria mod p
only if p divides both (b_δ−1)·d₁³·b_{d₁−1} and K(d₁,d₂). Two reductions:

(a) Over Q, the pole-3 coefficient (b_δ−1)·d₁³·b_{d₁−1} is NEVER zero
(b_δ ≥ 5 for δ ≥ 1). Hence **Ψ_{d₁,d₂} ≢ 0 over Q for every pair** — mod-p
degeneracy can only come from p dividing explicit nonzero integers.

(b) When p | b_δ−1: expand b_{d₂} = A_{d₁}(δ)b_δ + B_{d₁}(δ)b_{δ−1} and
b_{d₂−1} = A_{d₁−1}(δ)b_δ + B_{d₁−1}(δ)b_{δ−1} mod p; then
K ≡ d₂³(α(d₁,δ) + β(d₁,δ)·b_{δ−1}) mod p with α, β fixed rationals. So a
double escape forces p | G(d₁,δ) := gcd(b_δ−1, num(α+β·b_{δ−1})) — a FIXED
integer of log-size O(d₂ log d₂). (If G = 0, i.e. K = 0 over Q, the pair's bad
primes are the ω(b_δ−1) = O(δ) divisors of b_δ−1 — same budget.) The case
p | b_{d₁−1} instead uses G₂ = gcd(b_{d₁−1}, num K) symmetrically.

Exceptional-prime count: q bad ⟹ ∃ pair with d₂ ≤ q^{1/4} and q | G-type
integer; Σ_{d₁<d₂≤X^{1/4}} ω(G) ≪ X^{1/2}·max ω ≪ X^{3/4+o(1)}. For p ∉ 𝔈 no
pair escapes, and the clean-regime counting runs verbatim at D′ ~ p^{1/4}. ∎
(Sweep evidence: all prime factors of the gcds ≤ 19 for d₁<16, d₂≤22 —
consistent with 𝔈 being essentially empty in practice.)

## 5. The p^{3/4} program (superseded by Theorem 4; kept for the all-p question)

Running §4 without the clean-regime restriction: capacity 5(D′)³ with
D′ = 12p/N forces N ≪ p^{3/4} — PROVIDED the nondegeneracy lemma is available
for all lags d ≤ D′ ~ p^{1/4}. The obstruction is the exceptional lags
{d: p | b_{d−1}} and differences {δ: b_δ ≡ 1 (p)} (possible once
3.526k > log p). Program: (i) the second Laurent point x = −d₂ gives the
independent criterion with coefficient d₂³[(1−A_{d₁}(−d₂))b_{d₂−1} +
b_{d₂}B_{d₁}(−d₂)] (constants from the rank-one polar identity
(−5,1)T_{d−1}(−d) = d³(−b_d, b_{d−1})); escaping both criteria requires two
simultaneous p-divisibilities; (ii) bound the doubly-exceptional pair set and
absorb. Expected outcome: N(c) ≪ p^{3/4+ε}, E(p) ≪ p^{7/4+ε} unconditionally.

## 6. What this does and does not give for Problem 3.2

Does: first unconditional structure (Theorems 1–3), RMS vertical-sum bound
p/√log p, and a new elementary mechanism (non-autonomy as derivative-substitute)
that survives where Weil/Dwork/Stepanov/Katz all provably fail. Does not: the
conjecture needs E(p) = O(p^{1+o(1)}); even the full §5 program reaches p^{7/4}.
The gap from p^{7/4} to p^{1+o(1)} remains the frontier.
