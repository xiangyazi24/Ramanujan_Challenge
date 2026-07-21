# CODEX SPEC — π3: unconditional pair cross-correlation #{p: p|b_m and p|b_{m+g}} ≪ 1+g

ASSUME achievable; short run. Deliver a clean unconditional theorem + proof.

## Context (Fable pivot π3, after the crystalline route was killed at G1)

The crystalline eigenvalue route is a Hasse-invariant tautology (Lemma B′, G1):
Frobenius factor (X−p)(X²−a_{p,j}X+p³), slope-0 count ∈{0,1}, so "p|b_j ⟺
non-ordinary" restates rather than solves. Program pivots. This spec harvests the
one NEW unconditional structural constraint Fable flagged as worth a short run.

## Target

Prove, unconditionally, for m∈(N,2N] and 0<g≤N:
   D(m,g) := #{p prime ∈ (N, 2N] : p | b_m AND p | b_{m+g}} ≪ 1 + g,
with an absolute implied constant (no log factor). This SHARPENS the codegree
lemma (lem:codegree in proof.tex, which gives O(g·log N/log P_0) = O(g) only for
P_0=N^δ, and carries a log). The improvement to a clean O(1+g) with no log and
uniform in N is the deliverable.

## Route (Fable)

Use the Casoratian identity a_{n+1}b_n − a_n b_{n+1} = 6/(n+1)³ (Lemma
lem:wronsk in proof.tex) together with transfer-matrix height bounds, NOT just
the gap-polynomial degree count. Mechanism to develop:
- For a prime p∈(N,2N] dividing both b_m and b_{m+g}: since p>N≥m/2 and
  p>N≥(m+g)/2 requires m+g<2p, both m,m+g have base-p digits (1, ·), so by Lucas
  p|b_m ⟺ p|b_{m mod p} and the two residues m mod p, (m+g) mod p differ by g
  (no wrap, as g≤N<p). So p | b_r and p | b_{r+g} with r=m mod p, r+g<p.
- Two zeros at gap g in the FIRST block force p | N_g(r) (gap polynomial,
  deg 3(g−1)) — that is the degree route giving O(g). To remove the log and get
  a clean absolute constant, bound the number of primes p | N_g(m) directly: the
  integer N_g(m) with m∈(N,2N] has an explicit height, and the Casoratian /
  Wronskian gives a SHARPER height or a factorization constraint than the naive
  |N_g(m)| ≤ (CN³)^g. Work out the exact height of N_g(m) from the transfer
  matrix product M(m+g−1)···M(m), M(n)=[[P(n),−n^6],[1,0]], and whether the
  Casoratian 6/(n+1)³ forces a common-factor cancellation that caps the number
  of primes p>N dividing it by O(1+g) rather than O(g log N/log N)=O(g).
- Alternative clean route: the resultant / Wronskian gives that b_r, b_{r+g}
  simultaneously zero mod p implies p divides a specific integer of size
  exp(O(g log N)) with at most O(1+g) prime factors ABOVE N (since each prime
  >N contributes ≥ log N to the log-height, and the height is O(g log N)).
  Make this precise: log|target| = O(g log N) ⟹ #{p>N : p|target} ≤
  log|target|/log N = O(g). Then use the Casoratian to remove one factor of the
  implied constant / the +1 structure. State exactly what O(1+g) needs.

## Honest scope (state in the result)

This does NOT kill the anchored star: the star's second point is forced by the
reflection congruence r ↦ p−1−r at gap g=p−1−2r, which is costless (the pair is
one reflected orbit, not an independent coincidence). So π3 is a CONSTRAINT
INVENTORY result — a new unconditional bound on genuine gap-g coincidences — not
a route to (♦). File it as such. Its value: it is the first unconditional
cross-correlation bound on the actual Apéry zero set beyond no-consecutive-zeros
(Lemma lem:no-consec), which it generalizes (g=1 case: D(m,1)=0).

## Deliverables
- problems/3.2/pi3_result.tex — the theorem D(m,g)≪1+g with full proof (cite
  lem:wronsk, lem:codegree, lem:no-consec, lem:gap-poly, the transfer matrix
  M(n)), and the honest anchored-star caveat.
- problems/3.2/pi3_verify.py — PASS/FAIL: compute D(m,g) directly for N up to a
  few thousand, all g up to say 50, confirm D(m,g)≤C(1+g) with a small explicit C;
  confirm D(m,1)=0 (no-consecutive-zeros generalization).

## Hard rules
- Unconditional. No numerics-as-proof (numerics only verify the constant).
- If the clean O(1+g) resists and only O(g) or O(g log N/log N) is provable,
  deliver that with a precise note on what blocks removing the log — do NOT
  overclaim O(1+g).
- Do not modify existing files. Run pi3_verify.py; include its tail.
