# ORACLE B (arithmetic): why the actual Apéry zero-sets cannot align

## The complementary target

Oracle A needs one horizontal (cross-prime) arithmetic fact that the
adversarial anchored-star / no-go family CANNOT fake. Your job is to SUPPLY
that fact from the arithmetic of the ACTUAL Apéry numbers b_n — turning the
empirical "V°/S → 1, Apéry zeros are pairwise uncorrelated" into a theorem, or
isolating the exact modular/crystalline input it needs.

The concrete pseudorandomness that must be established: as p varies, the first
zero r_p = min Z_p (equivalently the doublet gap h_p = p-1-2r_p, equivalently
the whole set Z_p) is EQUIDISTRIBUTED / DECORRELATED enough that
  S(θ) = Σ_{p∈D} [e(θ m_1(p)) + e(θ m_2(p))]  has √-cancellation on minor arcs,
i.e. Σ_p over the actual zero locations does not concentrate. The adversary's
Z_p* = {m_0 mod p, p-1-m_0 mod p} DOES concentrate (all hit one m_0); prove
the real sequence cannot.

## Arithmetic handles on Z_p (all in proof.tex — read the cited items)

- Lucas/Gessel: b_n ≡ Π b_{n_i} mod p over base-p digits. So Z_p is
  determined by the first-block zeros; the whole mod-p behavior is a tiny
  automaton on the digit b_{n_i}.
- Beukers congruence: b_{(p-1)/2} ≡ a_p(f) mod p, f = η(2z)^4 η(4z)^4 ∈
  S_4(Γ_0(8)), the weight-4 level-8 newform (LMFDB 8.4.a.a). So the CENTER
  value is a modular Frobenius trace. Non-ordinary primes (p | a_p(f)) are
  exactly the Z(p)-odd primes (rem:open).
- Reflection: b_{p-1-j} ≡ b_j (rem:orbit), so Z_p is palindromic; the
  projective orbit π_n = [B_n:D_n] ∈ P^1(F_p) is a palindrome.
- Sym² squareness (Caruso–Fürnsinn–Vargas-Montoya–Zudilin 2026): the
  Hasse–Witt polynomial H_p(t) = Σ b_j t^j factors as Δ^{ε_p} B_p(t)², with
  B_p anti-palindromic of degree ~(p-3)/2. Z_p = zeros of H_p = zeros of B_p²
  (plus discriminant). So r_p is a root of a SQUARE of an explicit
  Hasse-type polynomial tied to the K3/elliptic pencil
  E_t: y² = x(x-1)(x-t(1-t)).
- Gap polynomials N_h, PGL₂ generation (T(1),T(2),T(3) generate PGL₂(F_p)),
  exact random models: uniform g∈V_p^- gives f=g² → K~Poisson(1/2) rigorously
  (prop:random-sq). The GAP: H_p is deterministic, not uniform-random in V_p^-.

## The exact missing statement

Make ONE of these precise and prove it (or reduce AMTD to it):
- (B1) Zero-location equidistribution: for fixed q and residue c, the density
  of primes p with r_p ≡ c (mod q) is the "expected" one; more generally the
  h_p are equidistributed in short progressions enough to kill S(θ) on minor
  arcs. Route: r_p is the smallest root of B_p; relate r_p mod q to Frobenius
  of the elliptic/K3 pencil via the Hasse polynomial, then Chebotarev/
  Sato–Tate over the modular tower. WATCH the quantifier trap that killed the
  gap-polynomial Chebotarev route (Q5291): Chebotarev controls roots of a
  FIXED polynomial as p varies — here B_p itself moves with p. You must use a
  FIXED geometric family (the pencil E_t) whose Frobenius gives r_p, not a
  moving polynomial. This is the crux; get it right.
- (B2) Sym²-to-decorrelation: use H_p = Δ B_p² to show the character sum
  Σ_r e_p(a r) over Z_p (= F_p(a)) inherits square-root cancellation from a
  bound on B_p, via a Weil/Bombieri bound for the curve cut out by B_p.
  Determine whether B_p, as a specialization of a FIXED bivariate/relative
  Hasse invariant of the pencil, has bounded geometric complexity (genus,
  #components) uniformly in p — the property FFK-type theorems need. If YES,
  Oracle A's route opens; if NO, prove the obstruction.
- (B3) Frobenius zero-fiber selection: prove the anchored-star family is
  arithmetically IMPOSSIBLE for b_n. Concretely: the map p ↦ (which reflection
  orbit Z_p occupies) is governed by Frobenius on a fixed motive, so it cannot
  be a single fixed m_0 mod p for a positive-density set of p (that would
  force a_p-type traces into a fixed residue too often, contradicting
  Sato–Tate/Deligne). Make this a real theorem: bound
  #{p ≤ X : m_0 ∈ Z_p} for fixed m_0 — is it o(π(X))? (Compare: fixed target
  a=0 gives the non-ordinary primes, density 0. A fixed nonzero target m_0
  moving with p as m_0 mod p is the real question.)

## DEAD ends (do not repeat — from prior audits)

- Chebotarev on a MOVING polynomial (gap poly / B_p): quantifiers reversed
  (Q5291). Must anchor to a fixed geometric family.
- Sato–Tate for the center value alone: controls only b_{(p-1)/2}, density 0.
- Katz equidistribution needing FIXED rank/conductor: b_m recurrence has
  growing complexity (rem:linnik). Only usable if B_p is a fixed-family
  specialization (that is exactly (B2)'s question).
- Pretentious / multiplicative: B(m) not multiplicative.
- Stewart / Corvaja–Zannier: constant-coefficient recurrences only.

## Ranked goals

(G1) Prove (B1) or (B2) or (B3) strongly enough to give AMTD / eq:short-arc
     (⟹ fully unconditional). Polylog loss OK.
(G2) Prove #{p≤X : m_0 ∈ Z_p} = o(π(X)) for every fixed m_0 (the minimal
     anti-concentration killing anchored star), or the fixed-family bounded-
     complexity property of B_p (opening Oracle A's route).
(G3) Reduce AMTD to a clean Frobenius/Sato–Tate statement over the fixed
     pencil E_t with all quantifiers correct.
(G4) Sharp obstruction: precisely why the b_n arithmetic does not yet yield
     cross-prime decorrelation, with the exact missing geometric input named.

## Computational duty (FIRST)

For p ≤ 20000: compute r_p, h_p, Z_p; test (a) equidistribution of r_p/p and
h_p/p (KS test vs the reflection-constrained uniform); (b) for several fixed
m_0, count #{p≤X : m_0 ∈ Z_p} and its growth vs π(X); (c) correlation of h_p
mod q with p mod q, small q; (d) verify H_p = Δ^{ε} B_p² and extract B_p's
degree/factorization type (irreducible? genus of its curve?) across p — is the
complexity bounded or growing? Write to problems/3.2/oracleB_exploration.md,
code oracleB_explore.py.

## Deliverables

problems/3.2/oracleB_result.tex (theorem/partial/obstruction, notation of
proof.tex — cite rem:orbit, rem:open, rem:squareness/Sym², Beukers, the
pencil E_t), oracleB_explore.py + oracleB_exploration.md, oracleB_verify.py
(PASS/FAIL, nonzero exit). STALL REPORT convention. Do NOT modify existing
files. No numerics-as-proof; verify identities symbolically; honest validity
domains and quantifier order (the moving-vs-fixed distinction is THE trap here
— be explicit about it in every claim).
