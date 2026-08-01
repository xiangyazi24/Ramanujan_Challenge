# Empirical factorial moments of the top-half Apéry bad-prime count

Date: 2026-08-01

## Scope and implementation

For every prime \(7\le p<N\), the scanner computes

\[
Z_p=\{0\le z<p:p\mid b_z\}
\]

from \(c_z=(z!)^3b_z\bmod p\) and scatters each zero pair \((p,z)\) to
\(n=p+z\). Consequently the scatter count is exactly

\[
H(n)=\#\{p\in(n/2,n]:n-p\in Z_p\}.
\]

The implementation is [CRON_moments_scan.c](CRON_moments_scan.c). Its hot
loop advances the cubic coefficient and \(r^6\) by finite differences, uses
two 32-bit Montgomery products per recurrence step, and distributes primes
largest-first over pthread workers. The analysis range below is
\(10^4\le n<2{,}000{,}000\), so omitting the exceptional primes below 7 has
no effect.

Compilation and runs:

    cc -O3 -march=native -Wall -Wextra -Wpedantic \
      -o /tmp/CRON_moments_scan problems/3.2/CRON_moments_scan.c -lpthread -lm
    /tmp/CRON_moments_scan 30000 1000 12
    /tmp/CRON_moments_scan 200000 10000 12
    /tmp/CRON_moments_scan 2000000 10000 12

The machine was an arm64 Mac running macOS 15.7.1, with Apple clang 17.0.0
and 12 online CPUs. At \(N=200{,}000\), 1,709,564,841 recurrence steps took
5.959 seconds in the scan kernel. The full \(N=2{,}000{,}000\) run completed
142,913,531,052 steps in 626.306 seconds (10.44 minutes), with total wall time
626.483 seconds. It found 148,930 primes and 149,112 zero pairs. The source
SHA-256 at the time of the run was
c2d551dde4c81745366c8b1cb7f42116ec68cfa5ae7fcdd9383b2b001501df1f.

## Verification

Three independent checks passed.

1. At startup, the optimized Montgomery/finite-difference kernel was compared
   position-by-position with a direct modular implementation for seven primes
   through 3137.
2. At \(N=30{,}000\), CRON_fresh_scan.py and the C scanner agreed exactly on
   the range \(1000\le n<30{,}000\):

   | \(H\) | Count |
   |---:|---:|
   | 0 | 26,847 |
   | 1 | 2,073 |
   | 2 | 76 |
   | 3 | 4 |

   Both gave \(\max H=3\).
3. The existing binary zero-pair file data_zp_pairs.bin (SHA-256
   8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d)
   was parsed independently. Its 149,112 pairs reproduced the full histogram,
   all six moments, all predictions and the top-20 list below.

The threaded run also checked that the completed recurrence-step count equaled
the precomputed total and that the merged pair list was strictly sorted, in
range, and duplicate-free.

## Histogram and maximum

There are 1,990,000 sampled integers in
\([10{,}000,2{,}000{,}000)\). Exactly 106,039 zero pairs scatter to an
integer below the cutoff.

| \(H(n)\) | Count | Fraction |
|---:|---:|---:|
| 0 | 1,887,467 | 94.847587939698% |
| 1 | 99,902 | 5.020201005025% |
| 2 | 2,579 | 0.129597989950% |
| 3 | 52 | 0.002613065327% |

Thus

\[
\boxed{\max_{10^4\le n<2\cdot10^6}H(n)=3}.
\]

There are 52 maximizers. The requested top 20, ordered by decreasing \(H\)
and then increasing \(n\), are:

| Rank | \(n\) | \(H(n)\) |
|---:|---:|---:|
| 1 | 11,576 | 3 |
| 2 | 18,444 | 3 |
| 3 | 22,101 | 3 |
| 4 | 26,164 | 3 |
| 5 | 47,066 | 3 |
| 6 | 47,859 | 3 |
| 7 | 63,887 | 3 |
| 8 | 64,555 | 3 |
| 9 | 105,433 | 3 |
| 10 | 182,715 | 3 |
| 11 | 184,862 | 3 |
| 12 | 218,132 | 3 |
| 13 | 305,455 | 3 |
| 14 | 318,441 | 3 |
| 15 | 321,882 | 3 |
| 16 | 328,313 | 3 |
| 17 | 338,471 | 3 |
| 18 | 403,315 | 3 |
| 19 | 406,011 | 3 |
| 20 | 408,460 | 3 |

## Factorial moments and Poisson test

For \(k=1,\ldots,6\), the observed factorial moment is

\[
M_k=\frac1{1{,}990{,}000}
 \sum_{n=10^4}^{1{,}999{,}999}(H(n))_k.
\]

The comparison column is the literal sample average of
\(\lambda_n^k\), where \(\lambda_n=\log 2/\log n\), rather than the power of
an averaged parameter.

| \(k\) | Observed \(M_k\) | Mean of \(\lambda_n^k\) | Ratio |
|---:|---:|---:|---:|
| 1 | \(5.287236180904523\mathsf{e}{-2}\) | \(5.148401739052273\mathsf{e}{-2}\) | 1.026966512888679 |
| 2 | \(2.748743718592965\mathsf{e}{-3}\) | \(2.666806039322719\mathsf{e}{-3}\) | 1.030725023890772 |
| 3 | \(1.567839195979900\mathsf{e}{-4}\) | \(1.391049525242686\mathsf{e}{-4}\) | 1.127090853006380 |
| 4 | 0 | \(7.314232070573247\mathsf{e}{-6}\) | 0 |
| 5 | 0 | \(3.881244518331225\mathsf{e}{-7}\) | 0 |
| 6 | 0 | \(2.081159032475822\mathsf{e}{-8}\) | 0 |

The zero fourth through sixth factorial moments are forced by the observed
maximum \(H=3\). For scale, the same Poisson mixture predicts only
0.581293 sampled integers with \(H\ge4\); equivalently, its fourth factorial
moment corresponds to 0.606472 expected four-subsets of hits. Seeing no
\(H\ge4\) value at this cutoff therefore has low resolving power for the
fourth and higher moments.

## Record structure probe

There is no \(n\) in the scan range with \(H(n)\ge4\), so the requested list
of bad primes, residues \(z=n-p\), ratios \(z/p\), classes \(p\bmod24\), and
orders

\[
\frac{p-1}{\gcd(p-1,n-1)}
\]

is empty. The C program prints all of these fields automatically for every
such \(n\) if a larger run finds one.

## Five-line interpretation

1. For \(k=1,2,3\), the ratios \(1.027,1.031,1.127\) stay close to one and show no growing-moment inflation.
2. For \(k=4,5,6\), the ratios are zero because the entire sample has \(H\le3\), not because the predicted moments are numerically large.
3. The Poisson mixture predicts only \(0.581\) observations with \(H\ge4\), so the missing upper tail is unsurprising at this sample size.
4. The data therefore support the program's moment-health heuristic through the orders that the sample can resolve, and show no empirical obstruction through \(k=6\).
5. Since \(\log(2\cdot10^6)\approx14.5\) while positive observed moments stop at \(k=3\), this run does not by itself verify the stronger \(k\asymp\log X\) hypothesis.
