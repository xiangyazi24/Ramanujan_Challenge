# CODEX SPEC — SECOND DEEP STRIKE: the (4.8)/(4.9) escape hatches (the campaign's chosen breakthrough)

Read CODEX_LASTWALL_report.md FIRST (your predecessor's complete survey; §4.3-4.4
defines the hatches). Report to CODEX_STRIKE2_report.md. First line in capitals:
either the breakthrough statement or NO UNBOUNDED L PROVED.

## The chosen target
N = p-2, D = sqrt(N)*L, L = L(p) -> infinity arbitrary slow. P_D = # bases r in I_p
with at least one collision return within (0, D] (= primitive/first-return support).
Q_D = same-base pair count. PROVE EITHER:
  (4.8)  P_D * Q_D << N^2        [=> capacity => S_D << N => E << N^{3/2}/L => 3/2 BROKEN]
  (4.9)  P_D << N/(L^2 log D)    [=> (4.8) via banked Q_D << N L^2 log D]
or any statement implying them.

## Why this is the thinnest wall segment (verified analysis, use it)
1. (4.9) asks ONLY a factor L^2 log D below the trivial P_D <= N. Empirical truth:
   P_D ~ S_D ~ 1.4 D = 1.4 sqrt(N) L — power-level room below the target.
2. The predecessor's reflection-symmetric no-go word (§4.4) does NOT block (4.8):
   in that word Q concentrates on ~q^2 bases, P*Q ~ q^8 << N^2 ~ q^10. VERIFY this
   claim first (compute P,Q for the word exactly, q=5,7,11); if it holds, (4.8) is
   compatible with every known abstract obstruction — the proof must use word-level
   inputs plus at most mild arithmetic, and no known counterexample class exists.
3. Chain check (verify exactly): P_D <= S_D and union-of-zero-loci structure:
   P_D = #{r : exists d <= D with N_d(r) = 0 mod p} = |union of root sets|.
   Since union <= sum, P_D <= S_1(D); but the useful direction is DISTINCTNESS:
   P counts each base once however many returns it has.

## Attack vectors (drive each to a terminal verdict)
V1. FIRST-RETURN RENEWAL: P_D = sum_{d<=D} P^prim_d (first-return exactly d; the
    predecessor PROVED primitive = first-return + exact split/renewal identities).
    The banked restart machinery IS a first-return theory for the fiber [0:1]:
    the paper's R_rel telescoping bounds spacings between consecutive zeros. Build
    the analogous FIRST-RETURN spacing theory for arbitrary bases: a base r with
    first return d means N_d(r)=0 and N_e(r)!=0 for all e<d. The set of such r =
    Z(N_d) minus union of earlier zero sets. Question: does the family have a
    SIEVE structure — the addition law N_{d}(r)=0 & N_{e}(r)=0 => N_{d-e}(r+e)=0
    (cascade) means zero sets are NOT independent: overlaps propagate. Derive the
    exact inclusion-exclusion/cascade constraints and determine whether they force
    |union| << sum with the needed factor, at least for the range d in (D/2, D]
    (dyadic shells: P_D <= sum over shells of |union over shell|; within a shell
    the degrees are comparable ~ 3D).
V2. AVERAGED CHEBOTAREV OVER THE FAMILY AT FIXED p VIA THE MOMENT THEOREMS:
    sum_{d<=D} R_d = S_1(D) is unknown, but the banked Sp-full theorems (h<=32)
    and the (4h-1)sqrt(p) bounds give EXPONENTIAL-SUM control per d. P_D via
    detector: P_D = sum_r [exists d: N_d(r)=0] — use a mollified/sieve detector
    (e.g. count r weighted by (number of returns +1)^{-1}: exactly P_D = sum_r
    k_r/(k_r) trivial... instead Cauchy: P_D >= S_D^2/(S_D+2Q_D) is the capacity
    LOWER bound; an UPPER bound on P needs the zero-detector first moment with
    multiplicity REMOVED — write P_D = S_1(D) - (multiplicity excess) and the
    excess = sum_r (k_r - 1)_+ >= S - P: so P <= S trivially and nothing better
    without first moments — CONCLUDE honestly whether V2 is circular; if yes say
    so fast and move on.
V3. (4.8) DIRECTLY VIA CLIQUE GEOMETRY: Q_D large requires triangle-rich collision
    graph (verified: Q = # triangles exactly, edges = S_1(D)); triangle-rich on
    few vertices (P small) or many vertices (P large)? Kruskal-Katona: T <= E^{3/2}
    and T-rich forces edge concentration on ~T^{1/3} vertices-ish. Derive the exact
    tradeoff: P_D * Q_D <= P * c S^{3/2} and S <= ... hmm S can be up to 1.5NL^2;
    get the best unconditional P*Q bound from {Q <= min(22 S_1-weighted form,
    (sqrt2/3) S^{3/2}), P <= min(N, S), S <= 1.5 D^2} and check how far from N^2
    it lands; then identify the minimal extra input closing the gap and hunt for
    it among banked arithmetic (apparition, -51 stratification, reflection parity,
    adjacent resultants A_m = +-prod((j!)^3 b_j)^6, the h<=32 certificates).
V4. NUMERICAL CALIBRATION (do early): exact P_D, Q_D, P*Q/N^2, first-return
    spacing distribution, shell unions vs sums, for p in {997,1999,4001,7919},
    D = sqrt(p)*log p and p^{0.6}: measure WHERE the union saves over the sum
    (cascade overlap rate) — this tells V1 exactly what to prove.

## Rules
No effort ceiling. Machine-verify every intermediate identity before relying on it.
Multiple vectors to terminal verdicts (PROVED/DEAD+reason/CONDITIONAL+exact residual).
If you prove (4.8) or (4.9): first line in capitals + full proof + numerical check.
