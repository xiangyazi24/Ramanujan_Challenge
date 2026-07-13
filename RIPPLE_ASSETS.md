# Ripple + Q-Series Assets for Ramanujan Challenge

## Ripple (~/repos/Ripple) — 0 sorry throughout

### Direct hits

| Problem | Ripple file | What's available |
|---------|------------|-----------------|
| **2.8** | `Chudnovsky1989.lean` (882 lines) | Coefficient recurrence, ₃F₂ bridge (`a_eq_3F2_coeff`), summability, Clausen connection |
| **2.8** | `Modular/CMEvaluation163.lean` (5422 lines) | `class_number_neg_163_eq_one`, Klein J infrastructure, CM reduction |
| **2.8** | `Modular/KleinJ.lean` (172 lines) | Klein J definition |
| **2.8** | `Modular/SingularModuli.lean` | Singular moduli machinery |
| **2.5** | `CatalanCertified.lean` (748 lines) | Complete Catalan constant via PIVP (ODE system) |
| **3.2** | `Number/Apery.lean` + 15+ Frobenius files | Apéry sequences, certificates, generating functions, Poincaré-Perron |
| **2.1** | `Hypergeometric/ThreeFtwo.lean` (474 lines) | ₃F₂ hypergeometric definition |
| **2.1** | `Hypergeometric/Clausen.lean` (710 lines) | Clausen identity |
| **2.8** | `Ramanujan1914.lean` | Original Ramanujan 1/π series |

### Gap (the one hole in Chudnovsky)
```lean
-- The CM period-derivative evaluation is a HYPOTHESIS, not proven:
-- Hypergeometric.chudnovskyCM163GaussDerivativeCombination = 640320^(3/2) / (12π)
```
This is the final CM/modular step. Everything else chains from it.

## Q-Series (~/repos/Q-series-and-Chan-s-work) — 172k lines

### Direct hits

| Problem | Q-series file | What's available |
|---------|--------------|-----------------|
| **2.1** | `Chapter07.lean` | Rogers-Ramanujan CF (level-5 structure matches 2.1's φ^{-10}) |
| **2.1** | `Pending/RamanujanQuintic.lean` | Quintic cyclotomic infrastructure, 5th root of unity products |
| **2.1** | `Pending/RamanujanQuinticJTP.lean` | Jacobi triple product for quintic |
| **2.8** | `Chapter20.lean` | Discriminant modular form (eta, τ(n)) |

### Shared infrastructure
- `QCalculus.lean` (3170 lines) — q-calculus, q-Pochhammer, q-binomial
- `Chapter01_GenFun.lean` — generating function framework
- Power series operations, formal q-expansion tools
- Partition function computations and congruences

## invitation-to-qseries-lean (public, 255k lines, 26.5k theorems)
This is the published formalization of Chan's book. Full infrastructure for:
- q-series manipulations
- Jacobi triple product
- Ramanujan theta functions
- Partition identities
- Rogers-Ramanujan type identities
