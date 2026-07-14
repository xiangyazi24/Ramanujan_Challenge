# Ramanujan Challenge — Status (2026-07-14)

## Score: 9/10 solved

**Deadline:** August 1, 2026

## SOLVED (9 problems)

| Problem | Topic | Method | Status |
|---------|-------|--------|--------|
| 2.1 | PCF → π | Q(√5) gauge, Poincaré roots | ✅ proof.tex |
| 2.2 | γ Apéry | Aptekarev recurrence (index shift m=n+3) | ✅ proof.tex |
| 2.3 | π+e | Series identification | ✅ proof.tex |
| 2.4 | harmonic+polylog | Weight-4 HPLs symbolic summation | ✅ proof.tex |
| 2.5 | Catalan CMF | Rivoal-Zudilin connection | ✅ proof.tex |
| 2.6 | ζ(2)+ζ(3) | GF ODE connection formula | ✅ proof.tex |
| 2.7 | 4-term ζ(2)+ζ(3) | Adjoint certificate (Lagrange bracket) | ✅ proof.tex |
| 2.8 | √10005/π | Chudnovsky formula in CMF disguise | ✅ proof.tex |
| 3.1 | knot π² | A-polynomial / Mahler measure | ✅ proof.tex |

## REMAINING (1 problem)

### P3.2 — gcd(d_n a_n, d_n b_n) = e^{o(n)} for Apéry sequences

**Status:** CONJECTURE (Section 3 = open problems). Computational evidence + partial proof in progress.

**What we have:**
- Computational evidence to n=200: log(gcd)/n decreasing from 0.36 to 0.06
- Supercongruence tower rigidity (established)
- p-adic analysis for small primes

**What we need:**
- Prime-counting estimate: show #{bad primes p <= n : p | gcd} = o(n/log n)
- Or equivalently: for most primes p, a_n uses its full d_n denominator budget

**ChatGPT questions dispatched (3 tabs running):**
- dm2: p-adic machinery / Lucas theory approach
- dm3: modular form / Chebotarev connection
- dm4: proof strategy / literature review

**Key scripts:**
- `scripts/p32_gcd_analysis.py` — computational evidence
- `problems/3.2/proof.tex` — partial proof
- `problems/3.2/notes.md` — notes and connections
