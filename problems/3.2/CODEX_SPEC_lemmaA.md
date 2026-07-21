# CODEX SPEC — Lemma A: Jacobi-sum skeleton + Stickelberger container for Z_p

ASSUME this target is achievable. Your job is to find the path and the explicit
formulas, not to judge feasibility. Deliver the algebra + numerics; report a
precise decision-gate verdict.

## Strategic context (from the Fable oracle — this is the route, follow it)

The full Apéry-GCD conjecture reduces (airtight, in proof.tex) to a dispersion
estimate (mh2). The oracle collapsed it:

LEMMA 0 (Lucas collapse, PROVEN + verified). For m∈(N,2N] and prime N<p≤m, the
base-p digits of m are (1, m−p), so Gessel–Lucas gives b_m ≡ 5·b_{m−p} (mod p).
Hence the detector X_p(m)=1_{b_{m−p}≡0} equals 1_{p|b_m}, and
  L_N(m) = #{p∈(N,m] : p | b_m}  (window prime-divisor count of ONE integer b_m).
So (mh2) is NOT cross-prime; it follows from the pointwise anti-concentration
  (♦)  max_{m∈(N,2N]} #{p∈(N,m] : p|b_m} ≪ N^{o(1)}.
(Verified: b_m≡5b_{m−p} and the collapse hold for all tested m,p. Size
calibration: log b_{2N}≈7·N vs Σ_{N<p≤2N}log p≈N, so SIZE PERMITS the
adversary — any proof needs non-archimedean/crystalline input on the prime
factorization of individual Apéry numbers.)

THE CRYSTALLINE DICHOTOMY (oracle Q1 — the key new mechanism). j∈Z_p means
𝔭 | S_j where S_j is a fixed-complexity toric exponential sum; this is INVISIBLE
to ℓ-adic topology. By Gross–Koblitz/Stickelberger, v_𝔭(Gauss/Jacobi sums) is
governed by DIGIT combinatorics: v_𝔭(g(ω^{−j})) = s_p(j)/(p−1) (sum of base-p
digits). So "character aspect mod p" = STICKELBERGER/digit aspect, NOT monodromy.
Concretely: b_j mod p is an alternating sum of ≤B fixed Frobenius eigenvalues
(B = uniform Betti bound of the fixed toric Kummer family); mod 𝔭 only slope-0
(unit) eigenvalues survive; if exactly ONE unit eigenvalue survives at character
j, then b_j ≡ ∓unit ≢ 0, so j∉Z_p. Hence
  Z_p ⊆ { j : (#unit eigenvalues mod 𝔭 at character j) ≠ 1 }.
This matches the empirical Poisson(1/2) law and explains the j=(p−1)/2 landmark
(b_{(p−1)/2} ≡ a_p(f), f=η(2z)^4η(4z)^4 the weight-4 level-8 newform).

DEAD ENDS (do NOT attempt — oracle-confirmed): ℓ-adic Mellin/FKM/FFK in the
j-aspect (character vanishing is Stickelberger, not geometric); two-Frobenii /
two-characteristic cohomology (no such formalism; and Lemma 0 makes p≠q moot);
Fourier detection Σ_r e_p(ar+s·b_r) (Heilbronn-type, outside all machinery);
pure dispersion/van der Corput / Hooley-Δ (blind to arithmetic, killed by the
size-permitted adversary).

## The marked coordinate you build on (verified)

Fixed Laurent polynomial (independent of p,j):
  Λ(x,y,z) = (1+x)(1+y)(1+z)·((1+y)(1+z)+xyz) / (xyz).
CT_{x,y,z} Λ^n = b_n (constant term = Apéry number, VERIFIED n≤4). Hence for
0≤j≤p−2:
  c_{p,j} := −Σ_{x,y,z∈F_p^×} Λ(x,y,z)^j ≡ b_j (mod p)   (VERIFIED mod 13, all j),
vanishing ⟺ j∈Z_p. (j=p−1: define b_{p−1}≡b_0=1.) Lift Λ(x,y,z)^j to the
character ω^j(Λ) (ω = Teichmüller) to work in characteristic 0; ω^j(t)≡t^j mod 𝔭.

## LEMMA A — three parts

### (A1) Jacobi-sum skeleton — the exact character-sum identity
Expand c_{p,j} = −Σ_{x,y,z∈F_p^×} Λ(x,y,z)^j by multiplicative orthogonality
(Σ_{t∈F_p^×} t^a = −1 if (p−1)|a, else 0; and the multinomial/Gauss-sum unfolding
of each factor of Λ^j). The numerator (1+x)^j(1+y)^j(1+z)^j((1+y)(1+z)+xyz)^j over
(xyz)^j should collapse to a SINGLE character-variable sum
  b_j ≡ −Σ_{k mod (p−1)} c(j,k) · Π_i J_i(j,k)  (mod 𝔭),
where each J_i is a Jacobi or Gauss sum with arguments AFFINE in (j,k), and there
are ≤ 5 factors (the three linear factors 1+x,1+y,1+z are standard binomial/Gauss
unfoldings; the quartic-in-monomials factor Q̄=(1+y)(1+z)+xyz is the nonstandard
one — unfold it carefully, it may need one extra summation variable). 
- The three linear factors: (1+x)^j = Σ_a C(j,a) x^a, and Σ_x x^{a-j+...} picks a
  residue; expressed via Gauss sums g(χ) and the Gross–Koblitz-friendly form.
- Q̄^j: expand ((1+y)(1+z)+xyz)^j = Σ_ℓ C(j,ℓ) ((1+y)(1+z))^{j-ℓ} (xyz)^ℓ, then
  merge with the x,y,z sums. Track how the y,z powers couple to the (1+y),(1+z)
  factors.
SUCCESS CRITERION: implement the derived skeleton and verify it reproduces
c_{p,j}=b_j mod p EXACTLY for ALL primes 5≤p≤200 and all 0≤j≤p−2 (PARI/gp or
Sage or pure Python with exact F_p arithmetic; hours of compute OK). 
CORRECTNESS ANCHOR: at j=(p−1)/2 the skeleton must reproduce the known
b_{(p−1)/2} ≡ a_p(f) mod p (Kilbourn/Ahlgren–Ono ₄F₃(1) evaluation) — verify a_p(f)
via the eta-product q-expansion for p≤200.

### (A2) Stickelberger valuation container
Apply Gross–Koblitz to each Jacobi/Gauss factor J_i(j,k): its 𝔭-adic valuation is
  v(k;j) = Σ_i (fractional-part expression ⟨(α_i j + β_i k + γ_i)/(p−1)⟩)  − (offset),
an EXPLICIT sum of base-p digit / fractional-part terms, computable in O(1) per
(j,k) once the α_i,β_i,γ_i are read off from (A1). Define
  U(j) = #{ k mod (p−1) : v(k;j) = min_k v(k;j) = 0 }   (unit-valuation terms).
Prove the theorem: if U(j)=1 then the single unit term is a 𝔭-adic unit, so
b_j ≢ 0 mod 𝔭, hence j∉Z_p. Therefore
  Z_p ⊆ { j : U(j) ≠ 1 }.
This is the FIRST structural container for Z_p beyond the size bound. State it
cleanly with the explicit U(j) formula.

### (A3) Decision gate — numerics that determine the whole program
For primes up to p≤2000 (stretch 5000): compute Z_p from the recurrence and
U(j) from the digit formula (A2) for all j. Classify every j∈Z_p as
  U=0 / U≥2-with-mod-𝔭-collision / non-ordinary (j=(p−1)/2 type),
and measure the container size |{j : U(j)≠1}| as a function of p. Report which:
  (S1) |{U≠1}| = p^{o(1)} (polylog): route FULLY ALIVE — (mh2) reduces to
       "sparse digit-defined set vs affine target j=m−p" (Mauduit–Rivat) +
       Poisson collision count.
  (S3) |{U≠1}| ≍ p^θ but digit-explicit: route ALIVE in Mauduit–Rivat form.
  (S2) generic j has U(j)≥2 and Z_p dominated by unstructured mod-𝔭 collisions:
       crystalline route alone does not give anti-concentration — report this
       honestly (it is itself a publishable structural finding about Z_p).
SECONDARY: from (A2)+per-term collision counting, attempt any unconditional
|Z_p| ≪ p^θ with θ<1 — a new theorem about Apéry numbers and a down payment on (♦).

## Deliverables
- problems/3.2/lemmaA_result.tex — the exact (A1) skeleton identity with proof,
  the (A2) container theorem Z_p⊆{U≠1} with explicit U(j), and the (A3) data +
  verdict (S1/S2/S3). Notation matching proof.tex + oracleC_result.tex (cite the
  marked coordinate Λ, prop:oracleC-*, eq:oracleB-mh2, the Lemma-0 collapse).
- problems/3.2/lemmaA_explore.py — (A1) skeleton verification p≤200, (A2) U(j)
  computation, (A3) classification + container-size table.
- problems/3.2/lemmaA_verify.py — PASS/FAIL per claim, nonzero exit on failure;
  MUST include: skeleton = b_j exact for p≤200 all j; γ_p anchor at (p−1)/2;
  U(j)=1 ⟹ j∉Z_p (no counterexample p≤2000); container Z_p⊆{U≠1}.
- If (A1)'s exact skeleton cannot be closed (the Q̄ factor resists), deliver the
  partial unfolding + a precise STALL REPORT naming exactly which character sum
  blocks, plus the (A3) data using a DIRECT computation of the unit-eigenvalue
  count (diagonalize the Frobenius on the toric cohomology numerically) so the
  decision gate still fires.

## Hard rules / negative acceptance (does NOT count as done)
- No numerics-as-proof: (A1) skeleton and (A2) container must be PROVEN
  identities/implications; numerics only verify. (A3) is explicitly empirical —
  label it so.
- A skeleton that only matches b_j "on average" or "for most j" is NOT the exact
  identity — it must be EXACT for all j, p≤200.
- Do not silently replace ω^j(Λ) by a wrong lift; ω^j(t)≡t^j mod 𝔭 must hold.
- Do not claim Z_p⊆{U≠1} without the U(j)=1⟹unit proof.
- A decision-gate verdict of (S2) is a VALID, valuable deliverable — do not fake
  (S1) to look successful.
- Do not modify existing files; deliver the three new files only.
- Run lemmaA_verify.py yourself; include its output tail in your final message.
