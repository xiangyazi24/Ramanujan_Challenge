# ORACLE A (analytic): short-arc Fourier orthogonality / bilinear dispersion

## The one target that closes fully-unconditional

Prove the short-arc Fourier orthogonality bound (proof.tex eq:short-arc,
subsection ssec:amtd) — with a POLYLOG loss allowed:
  ∫_0^1 | Σ_{P<p≤2P} Σ_{x<a/p≤x+1/N} (log p) F_p(a)/p |^2 dx
     ≪ N^{o(1)}/N · Σ_{P<p≤2P} Σ_{a=1}^{p-1} |(log p)F_p(a)/p|^2,
uniformly for √N < P ≤ N, where F_p(a) = Σ_{r∈Z_p} e_p(ar) is the Fourier
transform of the Apéry zero set Z_p = {r<p : p|b_r}.
By Gallagher's lemma this gives AMTD (hyp:amtd): V°(P,N) ≪ N^{o(1)} S(P,N),
which by prop:ap-bdh-suffices + lem:l2-uniform gives W(n)=o(n) for ALL n,
i.e. G_n = e^{o(n)} fully unconditionally. THIS IS THE GOAL.

## Precisely why it is open (do NOT re-derive — build past it)

The classical additive large sieve gives the RHS times (1 + 4P²/N). For
P>√N this factor exceeds a constant and reaches O(N) at P~N. The P²
barrier is intrinsic to the GENERIC large sieve (well-spacing of a/p alone).
Breaking it MUST use the special structure of F_p(a):
- palindromic phase: e_p(a/2)F_p(a) ∈ ℝ (from Z_p reflection r↔p-1-r);
- doublet dominance: most p have |Z_p|=2, Z_p={r_p, p-1-r_p}, so
  F_p(a) = 2 cos(π a h_p/p) e_p(a(3p-1)/2)... [verify the exact phase],
  h_p = p-1-2r_p the doublet gap;
- Sym² K3 Hasse structure: H_p(t) = Δ^{ε_p} B_p(t)² (Caruso et al.), i.e.
  the generating polynomial of (b_j mod p) is a discriminant times a square.

## DEAD ends (from W1-W5 + prior audits — do not repeat)

- Elliott–Halberstam / gen-EH / Barban–Davenport–Halberstam: WRONG object
  (fixed function in residue classes; here the tested function 1_{Z_p}
  MOVES with the modulus). Structural, not technical (ssec:amtd).
- Generic large sieve: P² barrier above.
- Dubhashi–Ranjan negative lattice, Janson: need nonneg intersection params;
  partial-period CRT is signed (rem:negative-crt).
- Martingale/Doob: cascade Δ_k ≡ dispersion itself (rem:martingale).
- Cauchy–Schwarz on E°: loose by 3000–6000× (tab:covariance).
- FFK2023 (Forey–Fresán–Kowalski) arithmetic Fourier transform: 4 mismatches
  (growing complexity, horizontal vs extension-degree, crystalline vs ℓ-adic,
  one- vs two-char dispersion) — rem:linnik. Do not claim it applies as-is;
  if you can REMOVE one mismatch, that is real progress.
- No-go model / anchored star: an ADVERSARIAL zero-set family satisfying all
  vertical bounds with V° = Ω(N²). So any proof MUST use a horizontal
  (cross-prime) arithmetic fact the adversary cannot fake.

## Attack vectors (pick what bites; combine freely)

1. DOUBLET MINOR-ARC CANCELLATION (rem:palindromic-fourier). Restrict to the
   ~60% doublet primes. The hit-position exponential sum is
   S(θ) = 2 Σ_{p∈D} e(θ(3p-1)/2) cos(π θ h_p). Diagonal AMTD ⟺
   ∫|S(θ)|² = O(T), i.e. √-cancellation on minor arcs. The prime-phase
   e(3pθ/2) is a linear phase in p (Vinogradov/Vaughan territory); the
   cosine modulation carries h_p. KEY: is h_p (mod small q) equidistributed
   / uncorrelated with the linear phase as p varies? h_p = p-1-2r_p, and r_p
   is the first zero of b_n mod p. Try Vaughan's identity on the p-sum with
   h_p as a "coefficient", reducing to bilinear sums Σ_m Σ_n α_m β_n over the
   r_p-structure. This is the most concrete route — push it hard.
2. LINNIK DISPERSION + FKM bilinear (rem:linnik). Expand V° directly, isolate
   diagonal S, seek off-diagonal cancellation. The off-diagonal is a bilinear
   Kloosterman-type form Σ_{p≠q} S_p(a)S_q(b) K_N(a/p+b/q). Reduce to a
   bounded-complexity trace-function bound for either a↦F_p(a) or m↦e_p(u b_m).
   The b_m phase comes from the NONAUTONOMOUS transfer matrix product
   M(m-1)...M(0), M(n)=[[P(n),-n^6],[1,0]] — its "complexity" grows with p,
   which is the obstruction. Can van der Corput / a q-analogue / a completion
   in the SHORT variable (the recurrence step) give a fixed-complexity
   bilinear gain? The Apéry shift-correlation target (prior Q5385) was
   |Σ_r e_p(a(b_{r+s}-b_r))| ≪ p^{1-η}; connect to it.
3. AMPLIFICATION / DFI. Bettin–Chandee / Duke–Friedlander–Iwaniec power
   saving for bilinear forms with Kloosterman fractions is the closest
   published framework. Determine EXACTLY which hypothesis on F_p(a) it needs
   and whether the palindromic/Sym² structure supplies it.
4. GALLAGHER + BLOCK REDUCTION (rem:block-decomp): cross-block covariance is
   0.03% of V°, so it suffices to prove within-block V°_P ≪ N^{o(1)} S_P for
   ONE dyadic block P~N (the top block is the whole difficulty). Focus all
   force on P ∈ (N/2, N].

## Ranked goals

(G1) Prove eq:short-arc (⟹ fully unconditional). Polylog loss OK.
(G2) Prove it for the doublet sub-sum only, or for P in a sub-range √N < P ≤
     N^{1-δ}, or conditional on a CLEAN, NAMED, strictly-weaker-than-AMTD
     Fourier/bilinear input (state it precisely; prove the reduction).
(G3) Reduce eq:short-arc to a bounded-complexity trace-function statement
     with ALL four FFK mismatches made explicit and AT LEAST ONE removed.
(G4) Sharp obstruction: the exact analytic lemma that is missing, in the
     style of the W1-W5 G4 deliverables (with the impossibility content: what
     abstract F_p(a)-family defeats the route).

## Computational duty (FIRST)

Extend scripts/variance_computation.py-style code: for the top block
P∈(N/2,N], N up to 8192+, compute (a) the doublet exponential sum S(θ) on a
fine θ-grid, its L²-norm vs T, and its large peaks' locations (major vs minor
arcs); (b) correlation of h_p with linear phases e(3pθ/2) at rational θ=a/q,
small q; (c) the true off-diagonal E° vs the bilinear-form prediction. Report
whether the doublet √-cancellation is empirically clean. Write to
problems/3.2/oracleA_exploration.md, code oracleA_explore.py.

## Deliverables

problems/3.2/oracleA_result.tex (theorem/partial/obstruction, notation of
proof.tex — cite eq:short-arc, hyp:amtd, rem:linnik, rem:palindromic-fourier,
prop:ap-bdh-suffices), oracleA_explore.py + oracleA_exploration.md,
oracleA_verify.py for any new identity (PASS/FAIL, nonzero exit). STALL REPORT
convention. Do NOT modify existing files. No numerics-as-proof; verify every
identity symbolically. Honest validity domains (recall the C4 and prop:column
lessons: an identity used in a special case must be stated for that case, and
checked against actual polynomials/data before use).
