# DOCTRINE — Ramanujan Challenge Automode (Session 4)

## Goal
All 10 Ramanujan Challenge problems addressed (9 unconditional + P3.2 conditional).
Deadline: Aug 1 2026. Current: **10/10 addressed**.
Focus: strengthen P3.2 (make unconditional if possible), polish all proofs.

## Current state
- **2.1** ✅ SOLVED (Q(√5) gauge, Poincaré roots)
- **2.2** ✅ SOLVED (Aptekarev recurrence)
- **2.3** ✅ SOLVED (Lambert × derangement)
- **2.4** ✅ SOLVED (weight-4 HPL)
- **2.5** ✅ SOLVED (Rivoal-Zudilin Catalan CMF)
- **2.6** ✅ SOLVED (ζ(2)+ζ(3) GF ODE)
- **2.7** ✅ SOLVED (adjoint certificate)
- **2.8** ✅ SOLVED (Chudnovsky CMF)
- **3.1** ✅ SOLVED (A-polynomial / Mahler measure)
- **3.2** ✅ ADDRESSED (conditional on Hypothesis Z; equivalence proved)

## P3.2 — the remaining gap

The proof is **structurally complete**: unconditional reduction to Hypothesis Z
(Z(p) = o(p) on average), plus unconditional O(√n) bound for small primes.

Hypothesis Z evidence is strong (p ≤ 10^4): mean Z(p) = 0.96, max = 8,
Poisson model matches. But Z(p) = O(1) is genuinely open — it contains
the non-ordinary prime density problem (Jin-Ma-Ono 2016) as a special case.

### Avenues to strengthen

(a) **Extend Z(p) to p ≤ 10^5 on uisai2** — more data, tighter power-law fit
(b) **Density-1 result** — prove gcd = e^{o(n)} for density-1 of n without Hypothesis Z
    (needs avg-Z over n, not over p — weaker and likely provable from Mertens)
(c) **Hasse-Weil / character sum bound** — if b_j mod p can be expressed as a
    character sum of bounded complexity, Weil gives Z(p) = O(√p), done
(d) **Direct finite-field hypergeometric** — Kilbourn (2006), Greene (1987)

Terminal: P3.2 is a Section 3 open problem. The conditional proof + computational
evidence + equivalence result is a strong submission even without (c)/(d).

## Resources
- mini: Python/mpmath, LaTeX, ChatGPT bridge (dm1-dm4)
- uisai2: available for heavy computation
- ChatGPT Pro: 4 functional tabs (dm1-dm4), SOL model
- Fable agents: dispatched for Z(p) bound literature
