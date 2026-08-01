# [CRIT-2H] critical-value audit

Verdict: **PASS through h=30**. For every `2 <= h <= 30`, the audit gives `s_h=4h-4`; the stronger full-Morse condition holds throughout the certified range. This is a finite-range certificate, not an all-h proof.

Exact integer resultants were computed through `h=12`. For larger h, degree/content-preserving reductions at auxiliary primes [1009, 65537, 1000003] provide characteristic-zero certificates.

## Conventions and certificate logic

The computation uses exactly

```text
P(X)=34X^3+51X^2+27X+5
N_1=1, N_2=P(X+1)
N_(h+1)=P(X+h)N_h-(X+h)^6 N_(h-1)
A_h=prod_(j=1)^h (X+j)^3
C_h=(A_h N_h'-A_h' N_h)/prod_(j=1)^h (X+j)^2
V_h(T)=primitive_part Res_X(C_h,N_h-T A_h)
```

For a modular row to count as a characteristic-zero certificate, the script checks `ell>3h`, preservation of the degrees and leading coefficients of `N_h` and `C_h`, `gcd(N_h,q_h)=1`, squarefreeness of `q_h`, nonzero reduction of the raw resultant (content survives), and `deg(V_h mod ell)=4h-4`. Only then is the modular gcd/simple-root count used.

## Gates

- Gate 0 — **PASS**. The h=1 cofactor equals `-((X+1)^2+(X+1)(Y+1)+(Y+1)^2)` and splits into two linear factors over `Q(omega)`, `omega^2+omega+1=0`.
- Gate 1 — **PASS**. All 78 pole-value identities for `1 <= j <= h <= 12` agree exactly.
- Gate 2 — **PASS**. Exact quotient, degree `4h-4`, and leading coefficient `-3 lc(N_h)` hold for all `1 <= h <= 30`.
- Gate 3 — **PASS**. `C_h(-h-1-X)=C_h(X)` holds exactly through h=30; the even `V_h`/`W_h(T^2)` law holds exactly through h=12 and at every eligible modular specialization thereafter.
- Gate 4 — **PASS**. For h=2,...,6, direct Q factorization, Singular absolute factorization over Qbar, direct modular factorization, and the modular [CRIT-2H] certificate all agree.

### Gate 4 details

| h | bidegree G_h | Q factors | absolute factors | modular prime | mod factors | [CRIT-2H] |
|---:|:---:|---:|---:|---:|---:|:---:|
| 2 | (5,5) | 1 | 1 | 1009 | 1 | PASS |
| 3 | (8,8) | 1 | 1 | 1009 | 1 | PASS |
| 4 | (11,11) | 1 | 1 | 1009 | 1 | PASS |
| 5 | (14,14) | 1 | 1 | 1009 | 1 | PASS |
| 6 | (17,17) | 1 | 1 | 1009 | 1 | PASS |

## Per-h audit

| h | source | deg V | deg gcd | v_T(V) | s_h | 2h-1 | 4h-4 | full Morse | [CRIT-2H] |
|---:|:---|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 2 | exact ZZ | 4 | 0 | 0 | 4 | 3 | 4 | yes | PASS |
| 3 | exact ZZ | 8 | 0 | 0 | 8 | 5 | 8 | yes | PASS |
| 4 | exact ZZ | 12 | 0 | 0 | 12 | 7 | 12 | yes | PASS |
| 5 | exact ZZ | 16 | 0 | 0 | 16 | 9 | 16 | yes | PASS |
| 6 | exact ZZ | 20 | 0 | 0 | 20 | 11 | 20 | yes | PASS |
| 7 | exact ZZ | 24 | 0 | 0 | 24 | 13 | 24 | yes | PASS |
| 8 | exact ZZ | 28 | 0 | 0 | 28 | 15 | 28 | yes | PASS |
| 9 | exact ZZ | 32 | 0 | 0 | 32 | 17 | 32 | yes | PASS |
| 10 | exact ZZ | 36 | 0 | 0 | 36 | 19 | 36 | yes | PASS |
| 11 | exact ZZ | 40 | 0 | 0 | 40 | 21 | 40 | yes | PASS |
| 12 | exact ZZ | 44 | 0 | 0 | 44 | 23 | 44 | yes | PASS |
| 13 | mod 1009,65537,1000003 | 48 | 0 | 0 | 48 | 25 | 48 | yes | PASS |
| 14 | mod 1009,65537,1000003 | 52 | 0 | 0 | 52 | 27 | 52 | yes | PASS |
| 15 | mod 1009,65537,1000003 | 56 | 0 | 0 | 56 | 29 | 56 | yes | PASS |
| 16 | mod 1009,65537,1000003 | 60 | 0 | 0 | 60 | 31 | 60 | yes | PASS |
| 17 | mod 1009,65537,1000003 | 64 | 0 | 0 | 64 | 33 | 64 | yes | PASS |
| 18 | mod 65537,1000003 | 68 | 0 | 0 | 68 | 35 | 68 | yes | PASS |
| 19 | mod 1009,65537,1000003 | 72 | 0 | 0 | 72 | 37 | 72 | yes | PASS |
| 20 | mod 1009,65537,1000003 | 76 | 0 | 0 | 76 | 39 | 76 | yes | PASS |
| 21 | mod 1009,65537,1000003 | 80 | 0 | 0 | 80 | 41 | 80 | yes | PASS |
| 22 | mod 1009,65537,1000003 | 84 | 0 | 0 | 84 | 43 | 84 | yes | PASS |
| 23 | mod 1009,65537,1000003 | 88 | 0 | 0 | 88 | 45 | 88 | yes | PASS |
| 24 | mod 1009,65537,1000003 | 92 | 0 | 0 | 92 | 47 | 92 | yes | PASS |
| 25 | mod 1009,65537,1000003 | 96 | 0 | 0 | 96 | 49 | 96 | yes | PASS |
| 26 | mod 1009,65537,1000003 | 100 | 0 | 0 | 100 | 51 | 100 | yes | PASS |
| 27 | mod 1009,65537,1000003 | 104 | 0 | 0 | 104 | 53 | 104 | yes | PASS |
| 28 | mod 1009,65537,1000003 | 108 | 0 | 0 | 108 | 55 | 108 | yes | PASS |
| 29 | mod 1009,65537,1000003 | 112 | 0 | 0 | 112 | 57 | 112 | yes | PASS |
| 30 | mod 1009,65537,1000003 | 116 | 0 | 0 | 116 | 59 | 116 | yes | PASS |

## Degree/content preservation audit

The following attempted reductions were rejected and were not used as characteristic-zero certificates:

| h | ell | tame | content survives | deg V found | deg V expected | reason |
|---:|---:|:---:|:---:|---:|---:|:---|
| 18 | 1009 | yes | yes | 66 | 68 | V degree drop, N degree drop, C degree drop |

Every h still has at least two independent eligible primes in the default run. A rejected reduction is a bad-reduction finding, not a failure of the characteristic-zero polynomial.

## Independent implementation checks

- **PASS**: the exact h=2 payload is `625 T^4 + 541064 T^2 + 22717712`, matching the banked baseline.
- **PASS**: SymPy independently recomputed the integer scalar resultant at `h=12, T=1`; its absolute value equals `content(raw V_12) * |primitive V_12(1)|` (1989 decimal digits).
- The modular reduction of every exact `V_h`, h<=12, was separately recomputed at all three primes and compared up to a nonzero scalar.

## Sequence and symmetry analysis

For every `2 <= h <= 30`, `s_h=4h-4`. Thus the offset from the Morse maximum is identically 0, while the surplus above the [CRIT-2H] threshold is `2h-3`. No repeated or zero critical value occurs in the certified range.

Every `V_h` is even. Since its constant term is nonzero and it is squarefree, its `4h-4` roots form exactly `2h-2` disjoint mirror pairs `{a,-a}`. Beyond this forced pairing, the audit finds no collision: `W_h(U)` is squarefree and has nonzero constant term in every certified row.

## Exact coefficient payload

`crit2h_results.json` stores every coefficient of the primitive integer `V_h` for `1 <= h <= 12` in ascending degree order, as decimal strings (to avoid JSON integer-width loss), together with contents and SHA-256 digests.

## Reproduction

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -u CRON_crit2h.py
```

Environment used: Sage 10.9; Python 3.14.3.
