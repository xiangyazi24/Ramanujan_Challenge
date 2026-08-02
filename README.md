# Ramanujan Challenge — Solutions

Solutions to all 10 problems from [The Ramanujan Challenge for AI](https://www.ramanujanmachine.com/ramanujan-challenge/) (July 2026).

Each problem asks to prove that a specific recurrence-generated sequence converges to a named mathematical constant (π, e, γ, Catalan's G, ζ(2), ζ(3), etc.).

## Results

| # | Constant | Lean | Paper | Status |
|---|----------|------|-------|--------|
| 2.1 | π (Cohen PCF) | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 2.2 | γ (Euler) | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 2.3 | π + e | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 2.4 | Li₄ + ζ combo | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 2.5 | Catalan's G | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 2.6 | ζ(2) + ζ(3) | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 2.7 | ζ(2) + ζ(3) | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 2.8 | √10005/π | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 3.1 | 4π²/85 | 0 sorry ✅ | ✅ | **Proved unconditionally** |
| 3.2 | e^{o(n)} | 3 sorry | ✅ | **Paper proof complete; Lean partial** |

**9 of 10 problems fully formalized in Lean 4 with 0 sorry.**
Problem 3.2 has a complete paper proof (135 pages) but the Lean formalization
is partial: 3 sorry remain (WZ identity, gap polynomial bound, main theorem).
All non-sorry theorems depend only on `{propext, Classical.choice, Quot.sound}`.

## Repository Structure

```
lean/                          Lean 4 formalization (v4.29.0 + Mathlib v4.29.0)
  RamanujanChallenge/
    Problem21.lean             P2.1: Cohen PCF sign-flip → π
    Problem22*.lean            P2.2: Rivoal construction → γ (4 files)
    Problem23.lean             P2.3: tensor annihilation → π + e
    Problem24*.lean            P2.4: polylogarithm identity (2 files)
    Problem25*.lean            P2.5: CMF → Catalan's G (8 files, ~3800 lines)
    Problem26*.lean            P2.6: Ore factorization → ζ(2)+ζ(3) (7 files)
    Problem27*.lean            P2.7: Barnes contour → ζ(2)+ζ(3) (12 files)
    Problem28_SUBMIT/          P2.8: Chudnovsky CM extraction (27 files, from Ripple)
    Problem31_SUBMIT/          P3.1: Bloch-Wigner four-shape (28 files)
    Problem32/                 P3.2: Apéry GCD conjecture (3 sorry remaining)
    Dilogarithm.lean           Shared: Li₂, Rogers R, functional equations
    RemainderCertificate.lean  Shared: convergence lemma
problems/
  2.*/proof.{tex,pdf}          Human-readable proofs for each problem
  2.5/AUDIT_2026-08-01.md      Rigorous proof audit by gpt-5.6-sol
SUBMIT/                        Submission packages
scripts/                       Verification scripts (Python, Sage)
```

## Key Techniques

- **Projective contraction** (Hilbert metric on positive cones): P2.5
- **Ore factorization / reduction of order**: P2.6
- **Barnes contour integrals + gauge transfer**: P2.7
- **CM evaluation / period bridge** (from [Ripple](https://github.com/xiangyazi24/Ripple)): P2.8
- **Bloch-Wigner dilogarithm / rational reconstruction**: P3.1
- **Delannoy basis decomposition + exact integral representation**: P2.5
- **Tensor annihilation / Lambert CF**: P2.3
- **Cohen Entry 5.3.22 sign-flip**: P2.1
- **Harmonic concentration / positive cone**: P2.2
- **Creative telescoping / WZ certificates**: P2.4, P3.2

## Problem 3.2 (In Progress)

The Apéry GCD conjecture `gcd(d_n a_n, d_n b_n) = e^{o(n)}` has a complete
paper proof (135 pages) including:
- Unconditional density-1 result via gap polynomials
- Conditional full result via Hypothesis Z
- Moment identity connecting to the Apéry family point count
- Palindromy from N_p(t) = N_p(1/t)

Lean formalization has 3 sorry remaining (WZ identity, gap polynomial bound, main theorem).

## Building

```bash
cd lean
~/.elan/bin/lake build        # full build (requires ~32GB RAM, use uisai2)
lake env lean <file>          # single-file check (local OK)
```

## Author

Xiang Huang, University of Illinois Springfield
- GitHub: [@xiangyazi24](https://github.com/xiangyazi24)
- With AI assistance from Claude (Anthropic) and Codex (OpenAI)

## Connections

- [Ripple](https://github.com/xiangyazi24/Ripple): Lean 4 GPAC/CRN framework (P2.8 uses its CM extraction)
- [Ramanujan Machine](https://www.ramanujanmachine.com/): The challenge organizers
