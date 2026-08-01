# Marked two-gap decoupling numerical diagnostics

## Coverage and exact conventions

- Main scan: 300 log-spaced primes in `[10^3, 10^6]` (requested sample size 300), including `1069, 1193, 1223, 1231, 1499`.
- At each prime, `H=floor(sqrt(p))`; triples use `0 <= s < s+h < s+k <= p-2` and `0<h<k<=H`.
- `T_refl` counts a triple once if any of its three pairs is `(r,p-1-r)`. All counts, recurrence values, polynomial coefficients, degrees, and root counts are exact.
- Diagnostic subsample: 20 primes, 50 stratified-random `(h,k)` pairs per prime; seed `3333824`.
- Wall time: 941.5 seconds. Incremental ledger: `research/scripts/marked_dec_progress.txt`.

## Marked triples by dyadic range

| range | sampled primes | mean (T-T_refl)/H^2 | max | argmax | max T | max T_refl |
|---|---:|---:|---:|---:|---:|---:|
| (512, 1024] | 1 | 0.0020811655 | 0.0020811655 | 1021 | 2 (1021) | 0 (1021) |
| (1024, 2048] | 30 | 0.00014255599 | 0.001953125 | 1049 | 4 (1699) | 4 (1069) |
| (2048, 4096] | 31 | 0.00033781385 | 0.0018903592 | 2143 | 8 (2351) | 8 (2351) |
| (4096, 8192] | 30 | 0.00012046631 | 0.00089106705 | 4597 | 6 (6353) | 4 (6353) |
| (8192, 16384] | 30 | 5.6044144e-05 | 0.00043402778 | 9403 | 4 (13001) | 4 (13001) |
| (16384, 32768] | 30 | 9.606159e-06 | 0.00012018509 | 16747 | 4 (32749) | 4 (32749) |
| (32768, 65536] | 30 | 1.0728954e-05 | 7.0616482e-05 | 56999 | 10 (59693) | 8 (59693) |
| (65536, 131072] | 30 | 6.3272146e-06 | 3.0046271e-05 | 67003 | 12 (103919) | 12 (103919) |
| (131072, 262144] | 30 | 1.4008912e-06 | 1.3353898e-05 | 150401 | 4 (217667) | 4 (217667) |
| (262144, 524288] | 30 | 2.1480443e-06 | 8.7791495e-06 | 455899 | 6 (387799) | 4 (387799) |
| (524288, 1048576] | 28 | 8.1511829e-07 | 4.8197359e-06 | 831253 | 4 (999983) | 4 (793717) |

OLS on `(log p, log max_block((T-T_refl)/H^2))` gives slope **-0.966477** from 11 nonzero dyadic maxima. A persistent positive slope would be the danger signal; slope near zero is the prediction up to `p^epsilon`.

## Row convention and verification

Writing `C_h(s)=prod_{j=1}^h(s+j)^3`, direct iteration of the Apéry recurrence gives

`C_h(s)b_{s+h}=U_h(s)b_s-s^3 U_{h-1}(s+1)b_{s-1}`.

Since `(s+h)^3D_h(s)=s^3C_h(s)` for `D_h=prod_{j=0}^{h-1}(s+j)^3`, the implemented row is exactly

`R_h=(s^3 U_h-(s+h)^3D_h, -s^9 U_{h-1}(s+1))`,

acting on `(s^3 b_s,b_{s-1})`. Its dot product is `s^6 C_h(s)(b_{s+h}-b_s)`. Thus the zero test is equivalent at admissible nonsingular indices. Random verification passed **400/400** points (20 per diagnostic prime).

## Phantom ratio

For every sampled pair the script constructs both continuant rows as dense polynomials, forms their cross determinant, and counts its roots by exhaustive exact evaluation at every `s in F_p` (a compiled C loop using the same continuant recurrence). `actual` is the exact marked-triple count for that same `(h,k)`.

| h+k bin | anomaly class | pairs | pairs with roots | actual sum | root sum | actual/root (pooled) | median conditional ratio | max ratio |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| small (h+k <= H/2) | ordinary | 340 | 340 | 0 | 10386 | 0 | 0 | 0 |
| medium (H/2 < h+k <= H) | ordinary | 340 | 340 | 0 | 23877 | 0 | 0 | 0 |
| large (H < h+k < 2H) | ordinary | 320 | 320 | 0 | 59618 | 0 | 0 | 0 |

## Degree and leading coefficient check

- Sampled determinants: **1000**.
- Exact degree `3(h+k)+9`: **1000/1000**.
- Leading coefficient equal to the stated `-c_(k-h-1)`: **0/1000** (counting only nominal-degree cases).
- Leading coefficient equal to the coefficient derived from the implemented row, `-c_(k-h-1)+c_(k-1)-c_(h-1)`: **1000/1000**.
- Apparition events `p | c_(k-h-1)`: **0**; nominal-degree drops among them: **0**.

The discrepancy is structural, not numerical: `U_m` has leading coefficient `c_m`, while both terms of `s^3U_m-(s+m)^3D_m` have degree `3m+3`, so the first row entry has leading coefficient `c_m-1`. Taking the cross determinant therefore adds `c_(k-1)-c_(h-1)` to the claimed coefficient. For example `(h,k)=(1,2)` gives leading coefficient `32`, not `-1`, over characteristic zero. Thus item 5's stated coefficient is incompatible with item 4's stated row normalization.

## Verdict

The marked-triple scan has global maximum `(T-T_refl)/H^2=0.0020811655` at `p=1021`, and dyadic-max slope `-0.966477`. The phantom scan found 0 anomalous root counts among 1000 pairs.

On the sampled range, the triple-count and phantom-root data support the proposed `H^2 p^epsilon` scale and reveal no structured exceptional family. However, the advertised leading-coefficient/apparition mechanism is not supported under the required row convention: its formula is algebraically incompatible with that normalization. The numerical evidence therefore supports `[GAP-MARKED-DEC]` as a counting conjecture, but not the stated leading-coefficient rationale without a normalization correction.
