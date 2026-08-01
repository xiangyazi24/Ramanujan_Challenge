# Ramanujan Challenge

Attacking all 10 problems from [The Ramanujan Challenge for AI](https://www.ramanujanmachine.com/ramanujan-challenge/) (July 2026).

**Deadline: August 1, 2026, 23:59 UTC**

## Problems

### Section 2: Proven (proofs known, not yet public)

| # | Name | Target | Status |
|---|------|--------|--------|
| 2.1 | Polynomial continued fraction for π | PCF = 6/(3-π) | Open |
| 2.2 | Euler's γ as Apéry limit | p_n/q_n → γ (4-term recurrence) | Open |
| 2.3 | π + e as Apéry limit | p_n/q_n → π+e (5-term recurrence) | Lean formalized (unconditional) |
| 2.4 | Harmonic/polylogarithm + zeta | Double sum = Li₄ + log + ζ combo | Open |
| 2.5 | Catalan's G rational approx | 3×3 matrix recurrence → G | Open |
| 2.6 | Series for ζ(2)+ζ(3) | 2077/720 + Σu_j = ζ(2)+ζ(3) | Lean formalized (unconditional) |
| 2.7 | 4-term recurrence for ζ(2)+ζ(3) | p_n/q_n → ζ(2)+ζ(3) | Open |
| 2.8 | Fast approx of √10005/π | 4×4 matrix recurrence → √10005/π | Open |

### Section 3: Conjectures (open problems)

| # | Name | Target | Status |
|---|------|--------|--------|
| 3.1 | Knot integral for π² | ∫ log-form over A-poly curve = 4π²/85 | Open |
| 3.2 | Apéry irrationality-measure optimality | gcd(d_n a_n, d_n b_n) = e^{o(n)} | Open |

## Approach

- **Fable oracle** for strategic route-finding on each problem
- **ChatGPT Pro** for exploratory computation and CAS derivation
- **Lean 4** formal proofs (leveraging Ripple + Q-series projects)
- **PDF** writeups for each solution

## Connections to Existing Work

- **Ripple**: zeta values, Catalan's constant, GPAC framework
- **Q-series/Chan**: q-series, modular forms, Ramanujan-style identities
- **HolonomicCRN**: holonomic recurrences, periods, ζ(3)

## Directory Structure

```
problems/          -- one subdir per problem (2.1/, 2.2/, ..., 3.2/)
  2.1/
    notes.md       -- working notes
    verify.py      -- numerical verification
    proof.tex      -- human-readable proof
  ...
lean/              -- Lean 4 project for formal proofs
papers/            -- compiled PDF proofs
scripts/           -- shared utilities (CAS, verification)
chatgpt-answers/   -- cross-problem oracle answers
```

The sole authoritative checkout is this repository.  Problem 3.2's current
frontier, failed-route log, exact experiments, and dedicated ChatGPT archive
are under [`problems/3.2/research/`](problems/3.2/research/) and
[`problems/3.2/chatgpt-answers/`](problems/3.2/chatgpt-answers/).
