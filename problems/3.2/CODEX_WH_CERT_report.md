# W_h irreducibility certificates

## Verdict

**PASS for every `2 <= h <= 32`.**  In this range `W_h(0) != 0`, `W_h` is irreducible over `Q`, and `gcd(W_h,W_h')=1`.

No reducible `W_h` was found.

This is an exact finite-range certificate, not a numerical-root test and not an all-h proof.

## Definitions and symbolic calibration

The computation uses the recurrence from `campaign3_questions/CTX.txt`:

```text
P(X)=34X^3+51X^2+27X+5
N_1=1,  N_2=P(X+1)
N_(h+1)=P(X+h)N_h-(X+h)^6 N_(h-1)
q_h=product_(a=1)^h (X+a),  D_h=q_h^3
C_h=q_h N_h' - 3q_h' N_h
V_h(T)=Res_X(C_h,N_h-TD_h)
V_h(T)=content_h W_h(T^2), with W_h primitive
```

SymPy recomputed the whole chain over `ZZ` for `h=2,...,6`.  It checked `D_h N_h'-D_h'N_h=q_h^2 C_h`, all expected degrees, exact resultant formation, vanishing of every odd coefficient of `V_h`, and the exact descent identity `W_h(T^2)=primitive_part(V_h(T))`.

| h | deg N | deg D | deg C | deg V | deg W | V even | W(0) != 0 |
|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 2 | 3 | 6 | 4 | 4 | 2 | yes | yes |
| 3 | 6 | 9 | 8 | 8 | 4 | yes | yes |
| 4 | 9 | 12 | 12 | 12 | 6 | yes | yes |
| 5 | 12 | 15 | 16 | 16 | 8 | yes | yes |
| 6 | 15 | 18 | 20 | 20 | 10 | yes | yes |

The h=2,...,6 primitive `V_h` coefficient vectors agree coefficient-for-coefficient with the independently banked exact payload in `crit2h_results.json`.

## Certificate logic

For `h<=12`, the exact primitive integer `W_h` was factored directly over `Q`; the exact gcd with its derivative was also computed.

For `h>12`, every resultant was formed directly in `F_ell[T]`, not by reducing a large characteristic-zero resultant.  A row was admitted only when `ell>3h`, the degrees of `N_h`, `C_h`, `V_h`, and `W_h` were preserved, `gcd(N_h,q_h)=1`, `q_h` was squarefree, `V_h` was even, and `W_h(0)` was nonzero. Thus the modular polynomial is a nonzero scalar multiple of the reduction of the primitive characteristic-zero `W_h`.

An irreducible good reduction proves irreducibility over `Q`.  Otherwise, a factorization with degrees `d_1,...,d_r` restricts the degree of any rational factor to a subset sum of those degrees.  Intersecting these subset-sum sets over good squarefree reductions proves irreducibility when only `0` and `deg W_h` remain.  Eight primes were tested for every `h>12`, as required in the all-reducible case.  Squarefreeness over `Q` is independently certified by one squarefree, degree-preserving reduction.

## Results

| h | deg W | W(0) != 0 | irreducible over Q | certificate | squarefree |
|---:|---:|:---:|:---:|---|:---:|
| 2 | 2 | yes | yes | direct Q factorization | yes |
| 3 | 4 | yes | yes | direct Q factorization | yes |
| 4 | 6 | yes | yes | direct Q factorization | yes |
| 5 | 8 | yes | yes | direct Q factorization | yes |
| 6 | 10 | yes | yes | direct Q factorization | yes |
| 7 | 12 | yes | yes | direct Q factorization | yes |
| 8 | 14 | yes | yes | direct Q factorization | yes |
| 9 | 16 | yes | yes | direct Q factorization | yes |
| 10 | 18 | yes | yes | direct Q factorization | yes |
| 11 | 20 | yes | yes | direct Q factorization | yes |
| 12 | 22 | yes | yes | direct Q factorization | yes |
| 13 | 24 | yes | yes | degree patterns mod 101,103,107 | yes |
| 14 | 26 | yes | yes | degree patterns mod 101,103,107 | yes |
| 15 | 28 | yes | yes | degree patterns mod 101,103 | yes |
| 16 | 30 | yes | yes | degree patterns mod 101,103,107,109,113 | yes |
| 17 | 32 | yes | yes | degree patterns mod 101,107,109,113,127 | yes |
| 18 | 34 | yes | yes | degree patterns mod 101,103,107,109,113,127,131,137 | yes |
| 19 | 36 | yes | yes | degree patterns mod 101,103 | yes |
| 20 | 38 | yes | yes | irreducible mod 103 | yes |
| 21 | 40 | yes | yes | degree patterns mod 101,103,107,109,127 | yes |
| 22 | 42 | yes | yes | degree patterns mod 101,103,107,109 | yes |
| 23 | 44 | yes | yes | degree patterns mod 101,103,107,109 | yes |
| 24 | 46 | yes | yes | degree patterns mod 101,103 | yes |
| 25 | 48 | yes | yes | degree patterns mod 101,103,107 | yes |
| 26 | 50 | yes | yes | degree patterns mod 101,103,107,113,127 | yes |
| 27 | 52 | yes | yes | degree patterns mod 101,103,109,113,127,131 | yes |
| 28 | 54 | yes | yes | degree patterns mod 101,103,107,109,127 | yes |
| 29 | 56 | yes | yes | degree patterns mod 101,103,107,109 | yes |
| 30 | 58 | yes | yes | degree patterns mod 101,103 | yes |
| 31 | 60 | yes | yes | degree patterns mod 101,103,107 | yes |
| 32 | 62 | yes | yes | degree patterns mod 101,103 | yes |

## Modular degree-pattern ledger

- `h=13`: 101: 6+18; 103: 1+2+6+15; 107: 3+5+16; 109: 1+23; 113: 1+2+3+18; 127: 1^2+2+20; 131: 1+23; 137: 1^2+2^2+8+10.  Intersection = `[0, 24]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=14`: 101: 1+2+6+17; 103: 2+24; 107: 1+6+19; 127: 1+25; 131: 1+11+14; 137: 1^2+3+4+17.  Intersection = `[0, 26]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=15`: 101: 1^2+26; 103: 5+6+17; 107: 1^2+26; 109: 1^3+2^2+5+6+10; 113: 2+26; 127: 8+20; 131: 2^2+3+6+15; 137: 1^2+3+23.  Intersection = `[0, 28]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=16`: 101: 1^2+3+25; 103: 4^2+22; 107: 1+3+26; 109: 3+4+5+7+11; 113: 1^2+28; 127: 1+2+6+21; 131: 9+21; 137: 1+29.  Intersection = `[0, 30]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=17`: 101: 4+5+10+13; 107: 1^2+2+6+8+14; 109: 1+2+6^2+17; 113: 2+13+17; 127: 5+27; 131: 1^2+4+5+21.  Intersection = `[0, 32]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=18`: 101: 1+3+5^2+20; 103: 1+6+27; 107: 1+14+19; 109: 1+2^2+5+24; 113: 1+5+6+7+15; 127: 1^2+2+3+6+7+14; 131: 1^2+4^2+10+14; 137: 4+6+24.  Intersection = `[0, 34]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=19`: 101: 14+22; 103: 5+31; 107: 3+16+17; 109: 2+3^2+4+5+7+12; 113: 10+26; 127: 1+9+26; 131: 1^2+2+4+5+6+7+10; 137: 1+4+12+19.  Intersection = `[0, 36]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=20`: 101: 2^2+6+12+16; 103: 38; 107: 1+10+27; 109: 2+5+7+8+16; 113: 1+3+34; 127: 1^2+2+14+20; 131: 1^3+3^2+29; 137: 7+31.  Intersection = `[0, 38]`; method = `irreducible reduction`; squarefree prime = `101`.
- `h=21`: 101: 1^2+2+14+22; 103: 1^3+2^2+13+20; 107: 1^2+2+5+6+11+14; 109: 14+26; 127: 2+16+22; 131: 2+3+4+9+22; 137: 1^2+2+5+6+8+17.  Intersection = `[0, 40]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=22`: 101: 1^2+4+36; 103: 1+2+3+36; 107: 1^2+2+4^2+5^2+20; 109: 3+19+20; 113: 1^2+6+34; 127: 1^2+3^2+9+25; 131: 1^3+7+8+24; 137: 1+13+28.  Intersection = `[0, 42]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=23`: 101: 1^3+41; 103: 1+5+8+30; 107: 1^2+42; 109: 8+11+25; 113: 1^3+3+8+30; 127: 4+7+9+24; 131: 3+4^2+33; 137: 1+3^2+7^2+11+12.  Intersection = `[0, 44]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=24`: 101: 2+13+31; 103: 1+45; 107: 2+5+39; 109: 4+5+9+28; 113: 3+4+6+15+18; 127: 1+5+6+15+19; 131: 1^2+6+38; 137: 1+2+9+34.  Intersection = `[0, 46]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=25`: 101: 1+3+8^2+10+18; 103: 4+5+39; 107: 1+6^3+29; 109: 1+2+4+5+9+11+16; 113: 1+2+7+38; 127: 5+6+13+24; 131: 1+19+28; 137: 3+4+7+34.  Intersection = `[0, 48]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=26`: 101: 1+2+3+5+39; 103: 1+2^2+3+42; 107: 2^2+5+41; 113: 3+5+7+35; 127: 25^2; 131: 1^2+3+5+40; 137: 1^2+2+12+34.  Intersection = `[0, 50]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=27`: 101: 2+3+8+18+21; 103: 1+2+3+46; 109: 2^2+10+38; 113: 1+2+20+29; 127: 1+2+5+7+8+9+20; 131: 3+49; 137: 1+2+14+35.  Intersection = `[0, 52]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=28`: 101: 1+3^4+4+6+31; 103: 1+2+3^2+5+40; 107: 2^3+3+6+39; 109: 2+7+15+30; 127: 1^3+2+18+31; 131: 1^2+52; 137: 1+2+9+42.  Intersection = `[0, 54]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=29`: 101: 1^3+2+3+19+29; 103: 1+3+52; 107: 1+2+4+49; 109: 2+7+47; 113: 2^2+52; 127: 4+19+33; 131: 1+2^2+51; 137: 1^2+5+6+10+33.  Intersection = `[0, 56]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=30`: 101: 1+28+29; 103: 21+37; 107: 22+36; 109: 10+48; 113: 1+6+7+20+24; 127: 1+2^2+4+5+8+14+22; 131: 1+2^2+3+50; 137: 2+12+20+24.  Intersection = `[0, 58]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=31`: 101: 1+12+47; 103: 2+12+20+26; 107: 2+16+42; 109: 1+5+54; 113: 1+3+56; 127: 1+6+8+45; 131: 1+59; 137: 1+13+46.  Intersection = `[0, 60]`; method = `degree-pattern intersection`; squarefree prime = `101`.
- `h=32`: 101: 11+51; 103: 1+3+17+41; 107: 1^2+2+7+15+36; 109: 8+9+19+26; 113: 3+6+53; 127: 2+3+11+46; 131: 1+26+35; 137: 1^2+17+43.  Intersection = `[0, 62]`; method = `degree-pattern intersection`; squarefree prime = `101`.

## Reproducibility

Run:

```bash
sage -python CODEX_wh_cert.py
```

CAS versions for this run: SymPy 1.14.0; Sage 10.9.

The full machine-readable record, including structural good-reduction gates, factor patterns, intersection traces, coefficient hashes, and timings, is in `CODEX_WH_CERT_results.json`.
