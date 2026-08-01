# Large scan of Apéry zero sets \(Z_p\)

Date: 2026-08-01

## Scope and headline

The threaded scan covered every prime \(7\le p<1000000\): **78495 primes** and **78462 zeros** in total. The mean was **0.999579591057**, and the largest zero set had size **12**. The complete per-prime data are in `CRON_zp_bigscan.csv`.

The scanner uses \(c_n=(n!)^3b_n\), whose zero set equals that of \(b_n\) for \(n<p\). Its division-free recurrence is evaluated with 32-bit Montgomery products and finite differences. This is algebraically the original Apéry recurrence, not a probabilistic shortcut.

## Mandatory sanity gates

All gates passed before the large scan:

- The full \(p<30000\) scan reproduced 3,242 primes, mean `1.0185` (exact sum 3,302), distribution `{0:1933, 1:2, 2:1008, 4:260, 6:34, 8:5}`, odd primes `{11,3137}`, and records `(7,0),(11,1),(17,2),(181,4),(379,6),(3727,8)`.
- Every scanned prime passed the parity/midpoint identity: `zp_size` is odd exactly when the midpoint state is zero.
- The original inverse-table recurrence and the independent factorial-table binomial sum agreed at all required spots:

| \(p\) | recurrence \(b_{(p-1)/2}\bmod p\) | binomial sum |
|---:|---:|---:|
| 13 | 9 | 9 |
| 101 | 51 | 51 |
| 3137 | 0 | 0 |

## Distribution of \(\lvert Z_p\rvert\)

| \(\lvert Z_p\rvert\) | prime count | fraction |
|---:|---:|---:|
| 0 | 47632 | 0.606815721 |
| 1 | 2 | 0.000025479 |
| 2 | 23729 | 0.302299510 |
| 4 | 6045 | 0.077011275 |
| 6 | 951 | 0.012115421 |
| 8 | 123 | 0.001566979 |
| 10 | 12 | 0.000152876 |
| 12 | 1 | 0.000012740 |

## Parity and midpoint primes

There are **2** odd cases. The complete list is:

```text
11, 3137
```

New midpoint primes above the previous \(p<30000\) scan:

```text
none
```

## Record breakers

| \(p\) | \(\lvert Z_p\rvert\) |
|---:|---:|
| 7 | 0 |
| 11 | 1 |
| 17 | 2 |
| 181 | 4 |
| 379 | 6 |
| 3727 | 8 |
| 88609 | 10 |
| 159977 | 12 |

## Running mean

The row at cutoff \(X\) uses primes \(7\le p<X\).

| cutoff \(X\) | primes | zeros | mean \(\lvert Z_p\rvert\) |
|---:|---:|---:|---:|
| 100000 | 9589 | 9490 | 0.989675670039 |
| 200000 | 17981 | 18124 | 1.007952839108 |
| 300000 | 25994 | 26086 | 1.003539278295 |
| 400000 | 33857 | 33876 | 1.000561183802 |
| 500000 | 41535 | 41414 | 0.997086794270 |
| 600000 | 49095 | 49100 | 1.000101843365 |
| 700000 | 56540 | 56514 | 0.999540148567 |
| 800000 | 63948 | 63928 | 0.999687245887 |
| 900000 | 71271 | 71186 | 0.998807369056 |
| 1000000 | 78495 | 78462 | 0.999579591057 |

## Poisson-pair fit

Restricting to the **78493 even rows**, put \(K=\lvert Z_p\rvert/2\). The maximum-likelihood fit is \(\widehat\lambda=0.499789790172\). Predicted counts are \(N e^{-\lambda}\lambda^k/k!\).

| \(k\) | \(\lvert Z_p\rvert=2k\) | observed | predicted | Pearson contribution |
|---:|---:|---:|---:|---:|
| 0 | 0 | 47632 | 47618.419881 | 0.003873 |
| 1 | 2 | 23729 | 23799.200080 | 0.207068 |
| 2 | 4 | 6045 | 5947.298607 | 1.605025 |
| 3 | 6 | 951 | 990.799708 | 1.598725 |
| 4 | 8 | 123 | 123.797894 | 0.005143 |
| 5 | 10 | 12 | 12.374585 | 0.011339 |
| 6 | 12 | 1 | 1.030782 | 0.000919 |

Pearson \(\chi^2=3.432092\) on **5 df**, using bins \(k=0,\ldots,6\). This convention reproduces the mandatory prefix value \(\chi^2=2.47\) on 3 df.

## Performance and reproducibility

- Threads: 8
- Exact recurrence steps: 37550245023
- Large-scan kernel time: 362.321 seconds
- Kernel throughput: 0.104 billion steps/second
- Sanity-gate time: 1.229 seconds
- End-to-end time for this invocation: 363.564 seconds

Compile and run from the repository root / output directory:

```sh
cc -O3 -march=native -Wall -Wextra -Wpedantic \
  -o /tmp/CRON_zp_bigscan problems/3.2/CRON_zp_bigscan.c \
  -lpthread -lm
cd problems/3.2
/tmp/CRON_zp_bigscan 1000000 8
```

## Final summary

```text
range: 7 <= p < 1000000
primes: 78495
mean |Z_p|: 0.999579591057
max |Z_p|: 12
odd/midpoint primes: 11, 3137
new midpoint primes above 30000: none
Poisson-pair lambda-hat: 0.499789790172
Poisson-pair chi^2 (df=5): 3.432092
sanity gates: PASS
```
