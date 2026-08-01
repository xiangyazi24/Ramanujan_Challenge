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

Patches from self-audit (2026-07-31, later same session):
- A_k is REGULAR at 0 with A_k(0) = b_k exactly (M(0) = [[5,0],[1,0]] is regular;
  T_k(0)e₁ = M(k−1)···M(1)(5,1)ᵀ = (b_k, b_{k−1})ᵀ). No limit subtlety.
- CLOSED FORM: B_k(x) = −x³/(x+1)³ · A_{k−1}(x+1) (second column of M(x) is
  (−x³/(x+1)³, 0)ᵀ). Hence ord₀B_k = 3 exactly, and the cross term
  B_δ(r+d)B_{d−1}(r) VANISHES to order 3 at r = −d — it contributes nothing to
  the pole-3 coefficient; no extra nonvanishing condition needed.
- Order-6 scare resolved: expanding (1−A_d)B_{d′} alone shows an order-6 pole,
  but in Ψ = (B_{d′} − B_d) − K with K = A_dB_{d′} − A_{d′}B_d =
  (r³/(r+d)³)B_δ(r+d), K is REGULAR at −d (pole³ × zero³), so the order-6
  pieces cancel exactly and the pole-3 coefficient (b_δ − 1)·d³·b_{d−1} stands.
- Boundary: transfers must not cross the singular index m = p−1; restrict the
  counting to r + d₂ ≤ p−2 and absorb the ≤ D tail elements into the O(1)/N-4
  slack (harmless in both regimes).
- SECOND closed form (verified d = 2..15): a_d := lim (r+d)³A_d(r) = **−d³·b_d**.
  Together with c_d = d³·b_{d−1}: the order-3 polar part of the transfer matrix
  T_d(r) at r = −d is RANK ONE: d³·(1,0)ᵀ⊗(−b_d, b_{d−1}) — equivalently the
  identity (−5, 1)·T_{d−1}(−d) = d³·(−b_d, b_{d−1}), provable by induction
  (the pole sits in the last factor M(r+d−1) ~ (r+d)^{−3}(1,0)ᵀ(−5,1)).
  This makes the SECOND Laurent criterion (at r = −d₂, for pairs whose δ is
  exceptional) fully explicit:
  pole-3 coeff of Ψ at −d₂ = d₂³·[(1−A_{d₁}(−d₂))·b_{d₂−1} − b_{d₂}·B_{d₁}(−d₂)]
  = d₂³·(b_{d₂−1} − b_{δ−1})  [closed form via the reflected scalar identity;
  sign fixed and verified 145/145 in §11] — an independent nonvanishing
  condition, so escaping BOTH criteria needs two simultaneous p-divisibilities.

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

## 5b. CLEAN COROLLARY — no exceptional-lag caveats (fully rigorous core)

If N := N(c) ≥ 25p/ln p (c ≠ 0), then D := 6p/N ≤ 0.24·ln p, so every lag in
play satisfies 3.53·d ≤ 0.85·ln p (margin: at the 21 constant this was 1.009·ln p
— OVER the line; 25 gives 0.85), hence the INTEGERS b_{d−1} and b_δ − 1 lie
strictly in (0, p) (b_k < (1+√2)^{4k} = e^{3.526k} for all k ≥ 1) — so mod p they cannot be 0 or 1, and the
Laurent criterion (§4) applies to EVERY pair: all Ψ_{d₁,d₂} ≢ 0, degrees ≤ 9D.
Root capacity ≤ (D²/4)·9D ≪ (ln p)³ < N/4. Contradiction. Therefore:

**Theorem 1 (unconditional):** max_{c≠0} N(c) ≪ p/log p.
**Theorem 2 (unconditional):** E(p) ≪ p²/log p.

To our knowledge these are the FIRST unconditional savings over trivial for the
Apéry value multiplicity and collision energy. The p^{3/4}/p^{7/4} upgrades need
only the exceptional-lag handling of §6.2 — the analytic core is identical.

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

## 8. Second-criterion data (K-test, 2026-07-31 late)

K(d₁,d₂) := pole-3 coefficient of Ψ at x = −d₂ = d₂³[(1−A_{d₁}(−d₂))b_{d₂−1} −
b_{d₂}B_{d₁}(−d₂)] — formula CONFIRMED against symbolic Laurent limits (4/4).
K ≠ 0 in all tested pairs. Crucially, for d₁ < 13, d₂ < 16:
max gcd(K, b_δ−1) = 37268 = 2²·7·11³ and max gcd(K, b_{d₁−1}) = 625 = 5⁴ —
ALL prime factors tiny. So a doubly-exceptional pair (both criteria ≡ 0 mod p)
requires p ≤ (small bound): for large p the escape set is EMPTY in the tested
range. If "p | gcd ⟹ p ≤ poly(d₂)" holds in general (algebraic reason hunted in
R4a), the p^{3/4} theorem follows with NO exceptional-element combinatorics at
all. Machine check of the whole collision mechanism: 60/60 fiber elements with
two forward lags root Ψ or the type-II polynomials (p = 101, 199; psi_e2e.py).

## 9. Sheltering analysis + the final gcd conjecture (solo, late session)

Counting doubly-exceptional pairs via integer sizes gives #bad pairs ≤ CD⁴/log²p
(≤ Cp/log²p at D ~ p^{1/4}), but one bad pair can shelter many fiber elements;
re-pairing sheltered elements against each other (mutual lags + Markov) caps a
single bad pair's shelter at ~p^{2/3} — not enough by itself. The REAL closing
statement is the empirically-supported arithmetic conjecture:

**GCD Conjecture (explicit, checkable):** every prime divisor of
gcd(b_δ − 1, K(d₁, d₁+δ)) is ≤ C·d₂^{O(1)}. (Data d₁<13, d₂<16: primes seen
only 2,3,5,7,11; gcd maxima 37268 = 2²·7·11³, 625 = 5⁴; the cube pattern 11³,
7³, 5³ tracks the d³-type constants, suggesting a closed form for the gcd like
the c_d = d³b_{d−1} discovery.)

Under this conjecture: no doubly-exceptional pairs for p > poly(D'), and the
p^{3/4}/p^{7/4} theorems follow with the SAME counting as the clean regime.
Publication structure: Theorems 1–3 unconditional; Theorem A/B conditional on
one explicit finite-checkable arithmetic statement (far cleaner than the
previous "missing sheaf" conditionality).

## 9b. GCD conjecture refined + sweep (d₁<16, d₂≤22)

K = d₂³·bracket. Sweep of gcd(b_δ−1, num(bracket)): max 5400 = 2³3³5² at
(2,22); ALL prime factors observed ≤ 19 (vs d₂ ≤ 22). Operative conjecture,
final form: **every prime divisor of gcd(b_δ−1, K(d₁,d₂)) is O(d₂)** — which
is exactly enough (counting uses d₂ ≤ D' ~ p^{1/4} ≪ p, so no doubly-
exceptional pair exists for large p, and Theorems A/B (p^{3/4}, p^{7/4})
follow by the clean-regime counting verbatim). Distribution mostly tiny
(4×145, 72×18, 24×11); script bracket_test.py.

## 9c. Wide sweep (d₁<20, d₂≤50): Theorem 4 exceptional set empirically EMPTY

Largest prime factor of any gcd(K, b_δ−1) or gcd(K, b_{d₁−1}) over the full
sweep: 89, at (d₁,d₂) = (3,37) — i.e. ≈ 2.4·d₂; no large residual factors at
all (script G_wide_sweep.py). The "p | gcd ⟹ p ≤ C·d₂" law holds on all ~850
pairs. Practical consequence: for every prime p > ~10³ the p^{3/4}/p^{7/4}
bounds of Theorem 4 hold OUTRIGHT in this lag range — the density-zero
exceptional set contains no prime visible in the data.

## 10. Theorem 5 (short-lag energy, c-uniform counting) — solo, parallel to R5

The Ψ polynomials do not depend on c. So for any pair (d₁,d₂) (nondegenerate
mod p), #{r: r, r+d₁, r+d₂ all in one fiber} ≤ deg ≈ 9(d₁+d₂), summing over ALL
values simultaneously. With L_Δ(r) = #same-fiber forward lags ≤ Δ of r:
Σ_r C(L_Δ(r),2) ≤ Σ_{d₁<d₂≤Δ} 9(d₁+d₂) ≈ 9Δ³, and Cauchy–Schwarz on
T_Δ := Σ_{d≤Δ} Y_d = Σ_r L_Δ(r) gives  T_Δ ≤ 2p + 6√(p·Δ³).

**Theorem 5a (clean, all p):** #collision pairs at lag ≤ c₀·log p is ≤ 2p + o(p).
**Theorem 5b (outside the Theorem-4 exceptional set):** for Δ ≤ p^{1/3},
T_Δ = O(p) — the random scale for short-lag collision mass.

Global E stays at p^{7/4}: the framework is intrinsically short-lag (transfer
degree grows with lag; no lag-reduction identity exists for single-digit r).
The amplification question (multiplicity on the Ψ family) is with R5.

## 11. R4a harvest (Q6306) — corrected sign, K closed form, (UN), and a FABRICATION flag

- ⚠️ FABRICATION: the consultant cited "Q32_paper_2026-08.tex Theorem E(1), commit
  831e1473" as repository-proved. NO SUCH FILE/COMMIT/THEOREM EXISTS. The
  universal-nonvanishing statement is NOT proved; its Plücker/descent proof
  sketch is the consultant's own construction. Treated as CONJECTURE (UN) until
  a full proof is delivered (R4b demanded it) and machine-verified.
- VERIFIED (their snippet + our cross-check, 145/145 pairs at p=101, including
  direct symbolic Laurent limits): second-pole coefficient of Ψ at x=−d₂ is
  K₀(d₁,d₂) = d₂³·(b_{d₂−1} − b_{δ−1}), δ = d₂−d₁ (sign: MINUS b_{d₂}B_{d₁};
  my §4b/§8 formulas carried a sign slip, now fixed). So the two endpoint
  criteria are: (b_δ − 1)·d₁³·b_{d₁−1} and d₂³·(b_{d₂−1} − b_{δ−1}) — double
  escape ⟺ p | gcd(b_δ − 1, b_{d₂−1} − b_{δ−1}): far cleaner integers for
  Theorem 4's exceptional-set count (heights O(d₂), same ω-summation).
- Independence confirmed by example: (1,4): p=19 kills witness 1 not 2; p=7
  kills 2 not 1. Also an INTERIOR Laurent family at x=−j, d₁<j<d₂ supplies
  further witnesses (formulas banked in Q6306.md §2).
- (UN) For all p ≥ 5, 1 ≤ h < k < p: Ψ_{h,k} ≢ 0 in F_p(x). IF (UN) holds:
  their §6 block/Hölder triple counting gives **N_p(c) ≤ 8p^{3/4} for ALL p and
  ALL c (including c = 0)** — no exceptional set, constant explicit. Sweep
  running (p ≤ 20000, h<k ≤ 30) to stress-test (UN).
- Their §7 critique of the integer-size ALL-p patch (run/occurrence hole):
  ACCEPTED — matches our own §9 sheltering analysis; the exceptional-set
  Theorem 4 is unaffected (bad primes excluded wholesale), the all-p claim
  needs (UN).

## 12. (UN) empirical status + stress tests (final numerics of the session)

- R3b's decisive stress test: primes ≤ 30000, pairs d₁<d₂≤10: ZERO flags
  (no case of Ψ ≡ 0 with nonzero predicted coefficient; u₂≡0 mod 13 retained
  and harmless). The nondegeneracy lemma survives its designated kill-shot.
- (UN) sweep: primes ≤ 20000 × all pairs h<k≤30 (~10⁶ combinations): ZERO
  identically-vanishing Ψ. Universal nonvanishing holds everywhere tested.
- Status ladder now: Theorems 1–3 (all p, audited, repaired); Theorem 4
  (exceptional-set p^{3/4}, G-integers now tiny closed forms); (UN) ⟹ all-p
  N(c) ≤ 8p^{3/4} incl. c=0 via the block/Hölder count — awaiting R4b's
  self-contained proof of (UN) for machine verification.

## 13. R8 harvest (Q6309) — downstream cashout audit + SECOND fabrication flag

- ⚠️ FABRICATION #2: R8 cited commits eda609dd / 732d85a1 / 2336ac53 as repo
  states — ALL THREE NONEXISTENT. Same consultant pattern as §11. All future
  prompts must state "no repository citations — we check". Its MATH, however,
  verifies independently:
- E(p) enters the gcd chain ONLY through |Z_p| ≤ √E(p). So E ≤ p^{7/4} gives
  |Z_p| ≤ p^{7/8} — WEAKER than the known zero-specific |Z_p| ≤ 3p^{2/3}.
  **The session's energy theorems buy zero downstream improvement for P3.2.**
  Their role: independent local value-distribution package.
- Even E = p^{1+o(1)} (the conjectural truth) gives only an N^{1/2+o(1)}
  exceptional set, NOT pointwise: the singleton-alignment countermodel (which
  IS genuine repo content — HANDOFF_P32: "Z_p = {N−p} ... H(N) = P_N") shows
  no per-prime vertical bound can cross the pointwise wall. Cross-prime k-th
  moment budget: with |Z_p| ≤ p^α, need M_k with k > 1/(1−α): α=1/2 ⟹ 3rd
  moment; α=2/3 ⟹ 4th; α=7/8 ⟹ 9th. The wall is horizontal, full stop.
- ⚠️ INTERNAL DISCREPANCY (for dm/DS to adjudicate): Q3.2_density_theorem.md
  §"E(p)=O(p^{1+o(1)}) suffices for mass X^{3/2} = URE = pointwise 3.2"
  CONTRADICTS the countermodel unless URE carries extra non-vertical input.
  One of the two claims must be corrected in the corpus.
- Bookkeeping gems kept: weighted master sum removes a log (avg log G_n ≪
  N^{2/3}, not N^{2/3}log N); d_n normalization (lcm³ vs lcm) needs unifying;
  crude pointwise constant 6 vs 3 depends on which reduction is authoritative.
- Recommended submission skeleton (adopted): Thm A (master sandwich) → Thm B
  (density-one, X^{2/3+o(1)} exceptional) → Thm C (this session's local
  package, explicitly labeled downstream-neutral) → Prop D (vertical-energy
  limitation warning) → Interface E (cross-prime moment criterion).

## 14. R7 harvest (Q6310) — the Mellin/Kummer lag reduction + closures

- ⚠️ FABRICATION #3: same phantom "commit 831e147 / Q32_paper_2026-08.tex" cited
  again (persistent confabulated repo state; ALL its repo citations = null;
  its EXTERNAL citation CFVZ Bull.Aust.Math.Soc 114 (2026) is real).
- **THE ONE GENUINE LAG REDUCTION**: b_{r+d} ≡ b_r (p) ⟺ the Mellin transform
  of the fixed K3 point-count trace, tensored by the Kummer twist χ_d − 1,
  vanishes at χ_r. χ_d has conductor O(1) INDEPENDENT of d: transfer degree
  O(d) → sheaf complexity O(1). Obstruction relocated to a NEW theorem: uniform
  two-point mod-p local limit for Kummer translates of the Apéry Mellin family.
  Per-translation √p-error ⟹ E(p) ≪ p^{3/2}; pay threshold needs AVERAGED
  error p^{1+o(1)} over all translates.
- Folding route: closed with proofs (θ does not descend past T^d−1 — the ideal
  is not differential; jet closure = Vandermonde = zero compression; 4-state
  block variant gives N⁵ ≪ d⁴ per digit ⟹ global p^{9/5} — worse than 7/4).
- Square-root convolution route: KILLED by counterexample, WE VERIFIED in one
  line: h_n = 4^{−n}C(2n,n) has (h*h)_M = 4^{−M}Σ C(2i,i)C(2M−2i,M−i) ≡ 1 for
  all M (classical identity) — E(h) small, E(h*h) maximal. No abstract
  energy-transfer inequality g → g*g can exist.
- Reflection/Lucas group: infinite dihedral, lag-invariant — PROVED closed.
  KEEP the one-end reflection tail lemma: collisions with d ≥ p−H reduce to
  collisions among the first H indices (extreme-tail removal at cost O(H²)).
- Exponent verdict: 7/4 stands this month; E ≪ p^{3/2} is the credible next
  target CONDITIONAL on the Mellin two-point local limit; p^{1+o(1)} needs the
  averaged version. Discriminating experiment (launching): stratify collision
  lags by Kummer order of d in (Z/(p−1))^× — spikes at bounded order ⟹ the
  3/2 route reduces to classifying finitely many self-twists; spikes at huge
  order ⟹ the wall is structural.

## 15. R5 harvest (Q6308) — E(p) ≪ p^{5/3} + the methodological ceiling

- ⚠️ FABRICATION #4: same phantom (commit 831e147 / Q32_paper_2026-08.tex /
  "Theorem E(1)") — its "unconditional" claims stay conditional-on-(UN) here.
- **Theorem 6 (verified by our own re-derivation): E(p) ≪ p^{5/3}**, conditional
  on nondegeneracy for lags ≤ p^{1/3} (i.e. (UN), or off the Theorem-4-type
  exceptional set). Mechanism: c-uniform TRIPLE counting — blocks of length H,
  Σ_c Σ_j C(n_{c,j},3) ≤ capacity(H) ≈ H³ ⟹ Σ_c m_c³ ≲ p²H; Cauchy–Schwarz
  E ≤ √(p·Σm³) = p^{3/2}√H; small fibers 6p²/H; H = p^{1/3} balances at 5/3.
  At H = c·log p this degenerates exactly to our clean Theorem 2 (p²/log p) —
  consistent. Parseval: RMS_h vertical ≪ p^{5/6}.
- Amplification on the Ψ/Δ family: DEAD with exact local algebra — the source
  is one-dimensional, all pair equations generate the same maximal ideal (t);
  Wronskian doubles the zero AND the degree (ratio unimproved); the universal
  product/discriminant captures all multiplicity at zero exponent gain.
- **ρ-ladder**: family-average root count R_H ≪ H^ρ ⟹ E ≪ p^{1+2/(1+ρ)}·-ish:
  ρ=3 (pointwise degrees, current) → 5/3; ρ=2 (O(1) avg roots/pair) → 3/2;
  ρ→1 → 1+o(1). Below 5/3 REQUIRES filtering candidate kernels to the
  distinguished Apéry orbit (two-base compatibility / anchored backward
  transfer / spectral input) — proved ceiling: an abstract block-permutation
  adversary realizes all one-base certificate consequences with E ≍ p^{5/3}.
- Kummer experiment (kummer_spikes.py, ~40 primes ≤ 1.8×10⁵): NO small-order
  concentration (top spikes at order ~p−1, count ≤ 10; order ≤ 20 strata avg
  2.26 vs global 2.35). Per R7's own decision rule: bounded-conductor Kummer
  translation sees no structure — the Mellin 3/2 route has no self-twist
  shortcut; the sub-5/3 wall is STRUCTURAL for this campaign.

## 16. Phase 3/2 opened (2026-08-01): rho = 2 MEASURED EXACT; germ reformulation

- ρ-measurement (rho_measure.py; p = 1009..40009, H = 4..32): R_H/H² constant
  ≈ 0.70–0.73, p-independent, H-independent; R_H/H³ ~ 0.7/H. **The ρ=2
  average-root law R_H ≤ C·H², C ≈ 0.72, is empirically exact.** Actual triple
  wastage T_H/R_H ≈ 0 — no candidate filtering needed: proving ρ=2 gives
  E(p) ≪ p^{3/2+o(1)} outright via the R5 ladder.
- Second-level determination CONFIRMED algebraically (solo): with v = (1,s),
  T_h(r)(1,s) = (1, A_{h−1}(r)+B_{h−1}(r)s); a second coincident lag h+δ forces
  σ_δ(r+h) = A_{h−1}(r) + B_{h−1}(r)·σ_h(r) — the shifted-base σ-value is an
  affine image; collisions telescope along the dual orbit.
- GERM REFORMULATION: parametrize projective solution germs by (base r, ratio s):
  σ_d(r) = s ⟺ the germ (r,s) returns at lag d. Then R_H = Σ_{germs}
  C(#returns in [1,H], 2): the ρ=2 law = "solution germs have Poisson-scale
  double-return statistics" (random model gives exactly H²/2). Per fixed δ:
  Σ_h Z(Ψ_{h,h+δ}) = total backward-return count of the one-parameter family
  γ_δ(x) (the forward-δ-return germs) — self-similar to the original problem;
  the circle must be broken by an actual counting theorem (R11 grinding the
  telescoping-product and recursion routes; R9 two-base; R10 spectral).

## 16b. Extended rho measurement + 10^7 scan cross-reference

- p = 100003, H up to 128: R_H/H² = 0.67, 0.71, 0.75, 0.75 — ρ=2 flat across
  two more octaves (mild upward drift at H ≥ 64: possibly a log factor or the
  reflection tail; still unambiguously ρ=2). Note: random-polynomial root
  statistics give EXACTLY average 1 root per polynomial — ρ=2 is precisely
  "the Ψ family has random-model root statistics" (an equidistribution-type
  statement about an explicit polynomial family).
- dm's 10⁷ scan (bn_bigscan_10M_report.md, d89add8): K(n) ≤ 4 all n ≤ 10⁷;
  mean |Z_p| = 1.0014 over 664579 primes, max |Z_p| = 12 — consistent with the
  Poisson vertical picture underlying all session theorems.
