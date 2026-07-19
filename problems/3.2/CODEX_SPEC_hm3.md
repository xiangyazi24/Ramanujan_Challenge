# CODEX SPEC W1: the third factorial moment (HM)_3

## Context

For the Apéry numbers b_n (recurrence (n+1)^3 b_{n+1} = P(n)b_n - n^3 b_{n-1},
P(n)=34n^3+51n^2+27n+5, b_0=1, b_1=5), let Z_p = {r in [0,p): p | b_r},
K_X(m) = #{p in (X,2X] : m mod p in Z_p}, lambda_X = sum_{X<p<=2X} Z(p)/p.
Falling factorial (K)_k = K(K-1)...(K-k+1).

PROVED (paper, ssec:high-moment): sum_{m<X^2}(K_X(m))_2 <= 5 X^2 lambda_X^2
(trivial: each pair of primes has <= X^2/pq + 1 <= 5X^2/pq CRT representatives
below X^2).

PROVED (thm:hm-pointwise): (HM)_k for ONE fixed k — namely
sum_{m<X^2}(K_X(m))_k << X^{2+o(1)} lambda_X^k — implies
log G_n << n^{2/3+2/k+o(1)} for EVERY n; any k > 6 closes the full
conjecture G_n = e^{o(n)}.

## Target

Prove (HM)_3, i.e. sum_{m<X^2}(K_X(m))_3 << X^{2+o(1)} lambda_X^3 — or any
nontrivial improvement over the trivial bound. Precisely, the trivial bound is
  sum_{m<X^2}(K_X(m))_3 <= #{zero-triples with CRT rep < X^2}
  <= (sum_{X<p<=2X} Z(p))^3 ~ (2X lambda_X)^3,
since for k=3 distinct primes p1p2p3 > X^3 > X^2, each residue triple has AT
MOST ONE representative below X^2. So the trivial/target ratio is X^3/X^2 = X.
Ranked goals:
  (G1) full (HM)_3;
  (G2) sum (K)_3 << X^{2+theta} lambda_X^3 for some fixed theta < 1;
  (G3) (HM)_3 under a clean auxiliary hypothesis STRICTLY WEAKER than
       Poisson/AP-BDH (name it precisely);
  (G4) a sharp obstruction theorem: what minimal arithmetic input is
       equivalent to (HM)_3.

## Tools known DEAD for this (do not rediscover; from prior audits)

Bombieri-Vinogradov (wrong quantifiers), classical large sieve (Q^2 barrier,
off by one log at k=2 and worse at k=3), Turan-Kubilius (needs pq < X^2),
Sato-Tate/non-ordinary primes (density 0), pretentious theory (not
multiplicative), Stewart & Corvaja-Zannier (constant-coefficient recurrences
only), abc/smooth numbers, hypergeometric motives / Katz equidistribution
(fixed rank needed), Fourier/bilinear Kloosterman with |S_p(a)| <= Z(p)
(no cancellation available for O(1)-size zero sets).

## Weapons available (all proved, in problems/3.2/proof.tex + nv_theorem.tex)

- Z(p) <= (3^{4/3}/2) p^{2/3}; sub-interval zero count O(|I|^{2/3}).
- Lucas: p | b_{m mod p} ==> p | b_m. So m with K_X(m) >= 3 has THREE prime
  factors of the single integer b_m in (X,2X], and log b_m ~ 3.5255 m.
- Gap polynomials N_h, nonvanishing mod p (h <= p); NEW: full-range bordered
  nonvanishing thm:nv-range (Delta_{h,k} != 0 mod p, h<k<p); uniform fiber
  bound max_a N_p(a) = O(p^{3/4}); collision energy E(p) = O(p^{7/4}).
- Codegree lemma: integers m != m' in (N,2N] sharing a bad prime p > P_0
  force p | N_{m'-m}(m), an integer of height (CN^3)^{m'-m} — at most
  O(h log N/log P_0) shared primes per pair at gap h.
- Reflection r <-> p-1-r; no consecutive zeros; block system p^3 a_n mod p.
- Wronskian a_n b_{n-1} - a_{n-1} b_n = 6/n^3.

## Structural observations to build on

(1) The k=3 sum decomposes by the TRIPLE of quotients. For a zero-triple
    ((p_i, r_i)) with CRT rep m < X^2, write m = q_i p_i + r_i, q_i < 2X.
    The three linear equations m = q_i p_i + r_i couple the primes.
(2) Fix m: its bad primes p | b_m (Lucas). So sum_m (K)_3 counts ordered
    triples of DISTINCT primes in (X,2X] all dividing b_m. Consider the
    integers c_m = gcd(b_m, prod_{X<p<=2X} p). Then sum_m (K)_3 =
    sum_m omega_X(c_m)(omega_X-1)(omega_X-2) where omega_X = # prime factors
    of b_m in the window. Height: omega_X(c_m) <= 3.53 m/log X pointwise.
    The AVERAGE of omega_X is X^2 lambda_X / X^2 = lambda_X (each prime
    divides b_m for Z(p) residues of m mod p — exact count per prime:
    #{m < X^2 : p | b_m} = Z(p)(X^2/p + O(1))). The k=3 question is whether
    the third moment of omega_X over m < X^2 is Poisson-like.
(3) Pair-correlation refinement: for two primes p, q, the joint count
    #{m < X^2 : p|b_m, q|b_m} = Z(p)Z(q) X^2/(pq) + E(p,q) with
    |E(p,q)| <= Z(p)Z(q) (CRT). The k=2 bound wastes nothing; for k=3 the
    error terms E over triples are what must cancel. Investigate whether
    the codegree lemma controls SIGNED sums of E(p,q) over prime pairs.
(4) Try a dispersion/Linnik expansion in the shortest variable, using
    Cauchy-Schwarz to reduce (HM)_3 to a bilinear form in (one prime) x
    (pairs of primes), then bound the inner pair-count by the PROVED (HM)_2
    machinery plus the codegree lemma. Track exactly where it fails if it
    fails; the failure point is the (G4) deliverable.

## Computational duty (do this FIRST, it informs the proof attempt)

Extend scripts/p32_hm_check.py methodology: for X in {256, 512, 1024, 2048,
4096 if feasible}: enumerate ALL zero-triples with CRT rep < X^2 (not via the
m-scatter but directly over triples, so you can INSPECT them). Report:
- the exact count vs X^2 lambda^3 (the (HM)_3 ratio, cross-check my R3 table);
- the distribution of the triples: quotient patterns, residue patterns
  (reflection pairs?), whether contributing triples cluster on special m
  (e.g. m near multiples of primes, m with small b_m-smooth structure);
- the pair-error statistics: histogram of E(p,q) over pairs, and the
  conditional third-moment given pair counts.
Write findings to problems/3.2/hm3_exploration.md; scripts to
problems/3.2/hm3_explore.py.

## Deliverables

1. problems/3.2/hm3_result.tex — theorem/partial theorem + full proof, OR a
   precise obstruction theorem (G4) with proof, notation matching proof.tex.
2. problems/3.2/hm3_explore.py + hm3_exploration.md (computational duty).
3. If a proof of G1/G2 is found: a verification script hm3_verify.py testing
   every new identity/lemma numerically (PASS/FAIL, nonzero exit on failure).
4. STALL REPORT convention as before if all goals fail.

## Hard constraints

No numerics-as-proof. Do not modify existing files. Verify any new identity
symbolically before use. Honesty about conditional steps.
