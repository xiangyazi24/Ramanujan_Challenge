# Large-scale top-half Apéry bad-prime scan

Date: 2026-07-19

## Headline

The scan reached the stretch target

\[
P_{\max}=N_{\max}=2{,}000{,}000.
\]

For every `1 <= n <= 2,000,000`, the top-half bad-prime count is at most
three:

\[
\boxed{\max B(n)=3.}
\]

There are 53 maximizers and no value with `B(n) >= 4`.  The complete argmax
list is

```text
321, 11576, 18444, 22101, 26164, 47066, 47859, 63887, 64555,
105433, 182715, 184862, 218132, 305455, 318441, 321882, 328313,
338471, 403315, 406011, 408460, 460045, 496157, 516294, 553239,
653968, 660487, 663216, 751905, 809290, 828985, 837983, 934904,
984120, 1006403, 1049082, 1098586, 1100091, 1133730, 1139584,
1150146, 1173575, 1254122, 1456509, 1468196, 1477328, 1578917,
1580628, 1637703, 1689988, 1698214, 1816821, 1824940
```

The global histogram is:

| `B(n)` | Count | Fraction |
|---:|---:|---:|
| 0 | 1,896,672 | 94.833600% |
| 1 | 100,670 | 5.033500% |
| 2 | 2,605 | 0.130250% |
| 3 | 53 | 0.002650% |

The global mean is `0.053019500000`, the population variance is
`0.052972432620`, and `Var(B)/E(B) = 0.999112262842`.

Following Step 2 of the specification and the established 200,000 reference
computation, this report uses primes `p >= 7`.  There is a minor discrepancy
between that mandated algorithm and the literal all-prime definition: for
`p=5`, `Z_5={1,3}`, which would add one hit at each of `n=6,8`.  Including
those two exceptional small hits would change only the global `B=0,1` counts
to `1,896,670` and `100,672`.  They do not affect the headline or any dyadic
or window table, all of which concern `n >= 1,024`; the mandatory prefix
table, like the established reference, uses `p >= 7`.

## Implementation and run

The scanner uses the division-free state

\[
A_m=(m!)^3b_m,
\qquad
A_{m+1}=(34m^3+51m^2+27m+5)A_m-m^6A_{m-1}\pmod p.
\]

For `m<p`, `(m!)^3` is nonzero modulo `p`, so `A_m` and `b_m` have exactly
the same zero positions.  Forward differences advance the cubic coefficient
and `m^6`; consequently each recurrence step uses only two logical modular
multiplications.  These are 32-bit Montgomery products, and the hot loop has
no division.  Primes are assigned dynamically, largest first, to 12 pthreads
(the detected online CPU count was 12).  The run used macOS 15.7.1 on arm64
(`Darwin 24.6.0`) and Apple clang 17.0.0.

Each recorded pair `(p,r)` contributes once to `n=p+r`.  Since `0<=r<p`,
this is exactly the condition `n/2<p<=n` and `r=n-p`, so one scatter produces
all top-half counts.

Compilation and the full command were:

```sh
cc -O2 -o problems/3.2/bn_bigscan problems/3.2/bn_bigscan.c -lpthread
problems/3.2/bn_bigscan 200000 200000 12 /tmp/data_zp_pairs_200000.bin
problems/3.2/bn_bigscan 2000000 2000000 12 problems/3.2/data_zp_pairs.bin
```

The middle command was the mandatory preflight; the final command was not
started until it passed.

The sieve found 148,930 primes in `[7,2,000,000]`.  The exact workload was
142,913,531,052 recurrence steps.  The mandatory 200,000 pre-run used
1,709,564,841 steps and took 4.503 seconds.  Scaling by the exact work ratio
predicted 376.4 seconds (6.27 minutes); a separate 500,000 benchmark predicted
about 383 seconds (6.4 minutes).  This was far below the three-hour cutoff, so
the 2,000,000 stretch target was used rather than dropping to 1,000,000.

The actual kernel time was 345.695 seconds at 0.413 billion recurrence steps
per second.  Sieve, sorting, binary output, scattering, validation, and all
statistics brought total wall time to 345.780 seconds.

## Correctness checks

Before the large scan, the Montgomery/finite-difference kernel was compared
position by position with a simple independent modular recurrence for seven
selected primes through 3137.  A separate inverse-based implementation of the
original `b_m` recurrence also agreed on all 606 pairs for every prime through
5000.

The required 200,000 gate passed exactly:

| `B(n)` | Count | Fraction |
|---:|---:|---:|
| 0 | 187,494 | 93.747% |
| 1 | 12,094 | 6.047% |
| 2 | 400 | 0.200% |
| 3 | 12 | 0.006% |

Thus `max B(n)=3` on this prefix.  Its mean, population variance, and
variance/mean ratio are respectively `0.064650000000`, `0.064830377500`, and
`1.002790061872`.

The full run also rechecked the saved million-prime reference: for primes
`7 <= p <= 1,000,000`, there are 78,495 primes and 78,462 zero pairs, with
maximum `Z(p)=12` at `p=159,977`.  The matched histogram was
`Z=0: 47,632`, `Z=1: 2`, `Z=2: 23,729`, `Z=4: 6,045`, `Z=6: 951`,
`Z=8: 123`, `Z=10: 12`, and `Z=12: 1`.

Finally, an independent parser checked the completed binary file: every `p`
is prime, every `0 <= r < p`, records are strictly sorted with no duplicates,
and an independent scatter reproduced all headline, dyadic, weighted, and
window statistics in this report.

## Per-prime zero statistics

The full file contains 149,112 zero pairs, so the mean zero count over the
148,930 scanned primes is

\[
\operatorname{mean}_p Z(p)=1.001222050628.
\]

The maximum remains 12, now attained at both `p=159,977` and `p=1,823,963`.
The only odd-zero exceptions are `p=11` and `p=3137`, each with `Z(p)=1`.

| `Z(p)` | Prime count | Fraction |
|---:|---:|---:|
| 0 | 90,375 | 60.682871% |
| 1 | 2 | 0.001343% |
| 2 | 45,020 | 30.228967% |
| 4 | 11,375 | 7.637816% |
| 6 | 1,875 | 1.258981% |
| 8 | 257 | 0.172564% |
| 10 | 24 | 0.016115% |
| 12 | 2 | 0.001343% |

After excluding the two odd-zero exceptions, the pair-count distribution
`K=Z/2` is close to `Poisson(1/2)`:

| `K` | `Z` | Observed fraction | `Poisson(1/2)` |
|---:|---:|---:|---:|
| 0 | 0 | 0.606836861 | 0.606530660 |
| 1 | 2 | 0.302293726 | 0.303265330 |
| 2 | 4 | 0.076379190 | 0.075816332 |
| 3 | 6 | 0.012589976 | 0.012636055 |
| 4 | 8 | 0.001725666 | 0.001579507 |
| 5 | 10 | 0.000161152 | 0.000157951 |
| 6 | 12 | 0.000013429 | 0.000013163 |

## Dyadic `B(n)` histogram and Poisson prediction

For each integer dyadic shell `(N,min(2N,N_MAX)]`, the prediction below is
the literal mixture specified by the prompt: for every integer `n` in the
shell, set `lambda_n=log(2)/log(n)` and sum its Poisson probabilities.  Entries
are fractions and are shown as `observed / predicted`.  The final shell is
partial because `2*1,048,576 > 2,000,000`.

| `N` | Upper end | Mean `lambda_n` | `B=0` | `B=1` | `B=2` | `B=3` | `B>=4` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,024 | 2,048 | 0.094786401 | 0.908203125 / 0.909570221 | 0.086914062 / 0.086208840 | 0.004882812 / 0.004088449 | 0 / 0.000129359 | 0 / 0.000003131 |
| 2,048 | 4,096 | 0.086576596 | 0.925781250 / 0.917067424 | 0.072265625 / 0.079392335 | 0.001953125 / 0.003438695 | 0 / 0.000099354 | 0 / 0.000002192 |
| 4,096 | 8,192 | 0.079675497 | 0.924560547 / 0.923417480 | 0.073486328 / 0.073570686 | 0.001953125 / 0.002932298 | 0 / 0.000077956 | 0 / 0.000001580 |
| 8,192 | 16,384 | 0.073793484 | 0.925292969 / 0.928864632 | 0.071411133 / 0.068541893 | 0.003173828 / 0.002530021 | 0.000122070 / 0.000062287 | 0 / 0.000001168 |
| 16,384 | 32,768 | 0.068720408 | 0.928283691 / 0.933588522 | 0.069030762 / 0.064154874 | 0.002502441 / 0.002205171 | 0.000183105 / 0.000050551 | 0 / 0.000000882 |
| 32,768 | 65,536 | 0.064300141 | 0.934295654 / 0.937724166 | 0.063598633 / 0.060294480 | 0.001983643 / 0.001939087 | 0.000122070 / 0.000041589 | 0 / 0.000000678 |
| 65,536 | 131,072 | 0.060414285 | 0.941238403 / 0.941374970 | 0.056976318 / 0.056871467 | 0.001770020 / 0.001718408 | 0.000015259 / 0.000034626 | 0 / 0.000000530 |
| 131,072 | 262,144 | 0.056971438 | 0.942214966 / 0.944621458 | 0.056068420 / 0.053815626 | 0.001693726 / 0.001533362 | 0.000022888 / 0.000029134 | 0 / 0.000000420 |
| 262,144 | 524,288 | 0.053899915 | 0.946109772 / 0.947527263 | 0.052558899 / 0.051070983 | 0.001289368 / 0.001376671 | 0.000041962 / 0.000024746 | 0 / 0.000000337 |
| 524,288 | 1,048,576 | 0.051142707 | 0.948938370 / 0.950143335 | 0.049713135 / 0.048592370 | 0.001327515 / 0.001242825 | 0.000020981 / 0.000021196 | 0 / 0.000000274 |
| 1,048,576 | 2,000,000 (partial) | 0.048751703 | 0.951230997 / 0.952417774 | 0.047614943 / 0.046431605 | 0.001135141 / 0.001131992 | 0.000018919 / 0.000018402 | 0 / 0.000000227 |

## Dyadic variance and logarithmic weight

Here `W_top(n)=sum log(p)` over the scattered bad primes, using natural
logarithms.  `Var(B)` is the population variance over the shell.  The table
reports one maximizer of `W_top(n)/n` in each shell.

| `N` | Upper end | Mean `B` | `Var(B)` | `Var/Mean` | Max `W_top/n` | Maximizing `n` | `W_top(n)` | `B(n)` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,024 | 2,048 | 0.096679687500 | 0.097098350525 | 1.004330413510 | 0.012490810025 | 1,041 | 13.002933236 | 2 |
| 2,048 | 4,096 | 0.076171875000 | 0.074275970459 | 0.975110176282 | 0.00592120343634 | 2,539 | 15.0339355249 | 2 |
| 4,096 | 8,192 | 0.077392578125 | 0.075309216976 | 0.973080607995 | 0.00372042539872 | 4,398 | 16.3624309036 | 2 |
| 8,192 | 16,384 | 0.078125000000 | 0.079101562500 | 1.012500000000 | 0.00238229434416 | 11,576 | 27.577439328 | 3 |
| 16,384 | 32,768 | 0.074584960938 | 0.075125560164 | 1.007248099619 | 0.00154681988272 | 18,444 | 28.5295459169 | 3 |
| 32,768 | 65,536 | 0.067932128906 | 0.068017061800 | 1.001250261031 | 0.000654434681486 | 47,066 | 30.8016227188 | 3 |
| 65,536 | 131,072 | 0.060562133789 | 0.060525953537 | 0.999402592842 | 0.000321157382280 | 68,310 | 21.9382607836 | 2 |
| 131,072 | 262,144 | 0.059524536133 | 0.059506146004 | 0.999691049614 | 0.000192954144159 | 182,715 | 35.2556164499 | 3 |
| 262,144 | 524,288 | 0.055263519287 | 0.055039968094 | 0.995954814391 | 0.000120528461233 | 305,455 | 36.8160211259 | 3 |
| 524,288 | 1,048,576 | 0.052431106567 | 0.052462999938 | 1.000608291010 | 0.000070679067029 | 553,239 | 39.1024163640 | 3 |
| 1,048,576 | 2,000,000 (partial) | 0.049941981703 | 0.049831575404 | 0.997789308809 | 0.0000381257535341 | 1,049,082 | 39.9970417691 | 3 |

The shell variance/mean ratios stay close to one, while the maximum
`W_top(n)/n` decreases from `1.249e-2` in the first reported shell to
`3.813e-5` in the final partial shell.

## Windowed localized dispersion

For every complete shell `(N,2N]` with `N >= 2^17`, the shell was split into
64 consecutive windows of length `N/64`.  If `S_j` is the sum of `B(n)` in
window `j`, the reported variance is the sample variance of the 64 values
`S_j` (denominator 63).  Under the iid, locally constant-mean Poisson
approximation, this variance is predicted by the window mean.  `Model mean`
sums `log(2)/log(n)` over each window and averages over the 64 windows.

| `N` | Window length | Observed mean | Model mean | Sample variance | Var/observed mean | Var/model mean | Maximum sum | Maximizing window |
|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| 131,072 | 2,048 | 121.906250000 | 116.677504216 | 115.610119048 | 0.948352681 | 0.990851834 | 153 | `(161792,163840]` |
| 262,144 | 4,096 | 226.359375000 | 220.774050485 | 281.376736111 | 1.243053159 | 1.274500946 | 280 | `(495616,499712]` |
| 524,288 | 8,192 | 429.515625000 | 418.961055275 | 567.777529762 | 1.321901921 | 1.355203599 | 481 | `(753664,761856]` |

Only complete shells are included: `N=1,048,576` would require data through
2,097,152, beyond this run's `N_MAX=2,000,000`.  The last two complete shells
show moderate local overdispersion relative to the independent model, but no
large localized pile-up: their largest window sums are 280 and 481 versus
observed means 226.36 and 429.52.

## Raw pair file

`data_zp_pairs.bin` is headerless.  Each record is exactly eight bytes:

```text
offset +0: uint32 p, little-endian
offset +4: uint32 r, little-endian
```

Records are sorted lexicographically by `(p,r)`.  The file contains every
zero pair for all scanned primes, including pairs with `p+r>N_MAX`; this is
why it contains 149,112 records while only 106,039 records enter the scatter.

- Size: 1,192,896 bytes
- Records: 149,112
- SHA-256: `8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d`
- Scanner source SHA-256: `2306b07d7fd07a8c987763d2a6daf156edccba6724c4b601ef38edaab73bd4c0`
