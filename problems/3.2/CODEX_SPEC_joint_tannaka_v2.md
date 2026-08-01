# CODEX SPEC — joint Tannakian moments v2: integral traces via Q6457 recipe (max)

## Mission

Round 2 of the joint-moment fingerprint. Round 1 (CODEX_JOINT_TANNAKA.md)
correctly stalled: the integral eigen-traces were underdetermined. That gap
is now CLOSED by the Q6457 recipe (chatgpt-answers/Q6457-2da55176.md at repo
root, also problems/3.2/chatgpt-answers/Q6457.md — read §5-§6). Implement it,
verify its internal checksums, then compute the Mellin moment fingerprints
with the CORRECT normalization.

## The recipe (from Q6457; implement exactly)

For good p > 3, eps = Legendre(-3,p)   # CORRECTED 08-01: five-prime empirical certification (29,31,37,41,43); Q6457's (2|p) fails at 37,41,43 (missing the degree-six isogeny (-3)-twist in the gauge), A_p(t) = Apéry truncation
sum_{k=0}^{p-1} b_k t^k mod p (b_k via recurrence mod p), for t in F_p^*:
d = t^2-34t+1, s = Legendre(d,p).

- s = +1 (split): u1,u2 = (1-t±sqrt(d))/16 in F_p;
  a_i = p+1-#E_{u_i}(F_p) for the model E_u: y^2+(1-2u)xy+u^2y=x^3;
  f_i = a_i^2 - p; GATE: f1 == f2. T_+ = T_- = f1; P = 2 f1.
  SPEEDUP allowed: for p > 16, a_p(E_u) is determined by its centered
  residue mod p (Hasse |a| <= 2 sqrt(p) < p/2): compute the Hasse-Witt
  residue via one O(p) sum (Legendre sum over x of the model's discriminant
  quadratic, i.e. a ≡ -sum_x chi(Delta_x) mod p style) OR just do the O(p)
  exact count per u — either way document which you used and cross-check
  the two methods at p=29 for all split u.
- s = -1 (inert): build F_{p^2} = F_p[z]/(z^2-d), u = (1-t+z)/16;
  a2 = p^2+1-#E_u(F_{p^2}) (exact count over F_{p^2}, O(p^2) per t);
  T_+ = eps*a2 - p; T_- = -T_+; P = 0.   # eps = (-3|p)
- s = 0 (rational branch point, exists iff Legendre(2,p)=1): u = (1-t)/16;
  a = p+1-#E_u(F_p); f = a^2-p; T_- = p; T_+ = f-p; P = f.

MANDATORY GATES per prime (abort prime on failure, report):
 (g1) T_+ + T_- == P pointwise (integers);
 (g2) T_+(t) mod p == A_p(t) for all t in F_p^*;
 (g3) T_-(t) mod p == s*A_p(t) mod p;
 (g4) |T_±(t)| <= 3p;
 (g5) Mellin inversion: -sum_t T_+(t) t^{-r} ≡ b_r mod p for all
      1 <= r <= p-2 (and the r=0 endpoint sees b_0+b_{p-1});
 (g6) at p=29: cross-check split traces against the raw counts already in
      CRON_pushforward_check.py conventions (t=2 -> f=7).

## Moments (the fingerprint; only after all gates pass)

Primes: all good p in [29, 149] (extend to 199 if total runtime < ~25 min;
report the exact set). Characters: chi_k(g)=e^{2 pi i k/(p-1)} for a
generator g; k = 1..p-2 (exclude trivial k=0); compute
S_±(chi_k) = sum_t chi_k(t) T_±(t) via one FFT per sign per prime
(store T_± indexed by discrete log).

Normalization: s_±(chi) = S_±(chi) / p^{3/2}.

Compute per prime and report trend + average over the largest 5 primes:
1. mu_2, mu_4, mu_6 of |s_+|: predictions printed alongside —
   Mellin group of A_+ has Tannakian dim 2; if G_+ = SL_2 (std):
   E|s|^2 = 1, E|s|^4 = 2, E|s|^6 = 5 (Catalan); alternatives:
   O_2-normalizer pattern mu_4 = 3; finite groups give rational
   spectra — print observed vs each.
   NOTE traces may be complex on non-self-dual characters; use absolute
   moments |s|^{2m} as primary (phase-robust), and also report real
   moments restricted to real characters (chi = chi^{-1}, i.e. k = (p-1)/2)
   separately.
2. mu_2, mu_4, mu_6 of |s_-|: dim 4; predictions: Sp_4 std: E|s|^4 = 3;
   Sym^3(SL_2): E|s|^4 = 4. This single number is the G_- dichotomy.
3. Joint: C22 = avg |s_+|^2 |s_-|^2; product prediction = mu_2(+)*mu_2(-);
   graph (Sym^3 coupling) prediction = higher (print both, from the
   observed single-object moments); naive covariance avg s_+ conj(s_-)
   (expect ~0 in BOTH cases — print to confirm uninformative).
4. Twisted correlations: C22 shifted by chi -> chi*eta for eta of order
   2, 3, 4, 6 (quadratic eta = character with k=(p-1)/2 etc.):
   avg |s_+(chi)|^2 |s_-(chi eta)|^2 — a twisted-graph detector.
5. Bad characters: report for each prime any chi where |S_±| deviates
   above the Deligne ceilings 2p^{3/2}/4p^{3/2} (should be none off the
   trivial character; the branch/endpoint corrections are already inside
   the trace table).

## Deliverables

- Script: research/scripts/q32_joint_tannaka_v2.py (stdlib + numpy ok).
- Report: CODEX_JOINT_TANNAKA_V2.md (problems/3.2 root): gate table,
  moment tables with predictions side by side, the G_- dichotomy verdict
  (Sp_4 vs Sym^3SL_2 by mu_4), the product-vs-graph verdict via C22, a
  LIMITATIONS section (finite-p fingerprint != theorem; arithmetic vs
  geometric group caveats), and an explicit statement of which Q6457
  checksums were verified.

## Hard constraints

- Do NOT touch: proof.tex, ERRATA.md, CAMPAIGN_MAP*, lean/, any existing
  script/report including the v1 joint-tannaka files. New files only.
- No external LLM/bridge dispatch. Exact arithmetic for all F_p/F_{p^2}
  work; floats only in the final FFT/moment stage (document precision).
- If any gate fails at any prime: STOP, write the stall report with the
  witness (that would falsify the Q6457 recipe — critical information
  either way).

## Acceptance

python3 research/scripts/q32_joint_tannaka_v2.py exits 0, prints
GATES VERIFIED per prime then the tables; report matches output.
