# Test A report: exact zero-detector complexity

## Run and exact checks

`CRON_testA_detector.c` was compiled and run with:

```text
cc -O2 CRON_testA_detector.c -o CRON_testA_detector
./CRON_testA_detector > results.csv
```

The run covered all 455 primes in `[500,4000]`. It completed without an
assertion failure and printed progress after every 50 primes. For every prime,
the program used integer arithmetic to check

```text
|S_p| E_p >= (p-1)^2,
```

which is exactly the rational Cauchy--Schwarz inequality
`|S_p| >= (p-1)^2/E_p`, without integer division or rounding. It constructed

```text
Q_p(x) = product_{v in S_p, v != 0} (1 - x/v)
```

when `0 in S_p`, checked its degree and normalization, evaluated it by Horner's
rule at every value in the image, and then checked
`Q_p(T_p(r)) = 1_{T_p(r)=0}` for every row index. When `0` was absent, it used
the specified zero-function convention and recorded degree `|S_p|`.

The required hook printed:

```text
# spot p=13 Z={}
# spot p=17 Z={3,13}
# spot p=29 Z={}
```

## Summary

The statistics below are over all 455 tested primes.

| Quantity | Mean | Minimum | Maximum |
|---|---:|---:|---:|
| `|S_p|/p` | 0.393596 | 0.374781 (`p=571`) | 0.414538 (`p=509`) |
| `E_p/p^(5/3)` | 0.020947 | 0.011700 (`p=3947`) | 0.048847 (`p=503`) |
| `E_p/p` | 2.995865 | 2.796282 (`p=1291`) | 3.220453 (`p=1193`) |
| `deg(Q_min)/p` | 0.393369 | 0.373030 (`p=571`) | 0.414538 (`p=509`) |

There were 291 empty-zero rows, recorded as zero functions, and 164 rows with
nonempty zero sets. On the nonempty-zero rows alone, the mean detector-degree
ratio was 0.392634. The observed zero-count distribution was:

| `|Z_p|` | 0 | 1 | 2 | 4 | 6 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| Number of primes | 291 | 1 | 129 | 30 | 3 | 1 |

## Which birthday model fits?

For a genuinely uniform random map with `m=p-1` independent inputs and `p`
possible values, a fixed value is missed with probability `(1-1/p)^m`.
Therefore

```text
E[|S_p|] = p (1 - (1-1/p)^(p-1))
          -> (1-e^(-1))p
          = 0.632121 p.
```

That prediction does **not** fit these rows: the measured mean is 0.393596.

The program also exactly checked, for every tested prime and every applicable
index, the row pairing

```text
T_p(r) = T_p(p-1-r),    1 <= r <= p-2.
```

Thus the row has only about `p/2` effective representatives: index zero, the
central fixed index, and `(p-3)/2` paired indices. Treating those representatives
as uniform samples gives

```text
E[|S_p|] approximately p (1 - (1-1/p)^((p+1)/2))
                         -> (1-e^(-1/2))p
                         = 0.393469 p.
```

This is the normalization behind the earlier `1-e^(-1/2)` figure. It matches
the data to three decimals: the measured limiting constant is **0.394**. The
mean over `p>=3500` is 0.393596, only 0.000127 above `1-e^(-1/2)`.

## Collision-energy growth

For `m` independent uniform samples,

```text
E[sum_v N(v)^2] = m + m(m-1)/p,
```

because the diagonal sample pairs always collide and each distinct ordered
pair collides with probability `1/p`. With `m=p-1`, this is asymptotic to
`2p`. For the observed doubled-pair row, the dominant contribution is instead

```text
4 [p/2 + (p/2)^2/p] = 3p + O(1),
```

which agrees with the measured mean `E_p/p = 2.995865`.

An ordinary least-squares regression of `log E_p` on `log p` over all 455
primes gives

```text
E_p = 2.992 p^1.00012,
R^2 = 0.99825,
95% slope interval = [0.99626, 1.00398].
```

The measured growth exponent is therefore **1.000**, not `5/3`. Consequently
`E_p/p^(5/3)` decays like `3p^(-2/3)`.

| Prime bin | Primes | Mean `|S_p|/p` | Mean `E_p/p` | Mean `E_p/p^(5/3)` | Mean `deg/p` |
|---|---:|---:|---:|---:|---:|
| 500--999 | 73 | 0.393667 | 2.988260 | 0.037037 | 0.393110 |
| 1000--1499 | 71 | 0.393587 | 2.996117 | 0.026003 | 0.393316 |
| 1500--1999 | 64 | 0.392179 | 3.015637 | 0.020841 | 0.391975 |
| 2000--2499 | 64 | 0.393336 | 3.000835 | 0.017539 | 0.393186 |
| 2500--2999 | 63 | 0.393888 | 2.990470 | 0.015270 | 0.393761 |
| 3000--3499 | 59 | 0.395024 | 2.984204 | 0.013606 | 0.394914 |
| 3500--4000 | 61 | 0.393596 | 2.995568 | 0.012457 | 0.393500 |

The detector-degree profile is therefore linear and stable near `0.394p`.
The `1/p` difference between `|S_p|/p` and `deg/p` on nonempty-zero rows is
negligible at this scale.

## Verdict

`4000` itself is composite. The scale prediction at the endpoint is
`(1-e^(-1/2))*4000 = 1573.88`; the nearest tested prime with a nontrivial
indicator is `p=3943`, where `deg Q_min=1571`. **Fixed-degree detectors are off
by factor 1574 at p=4000** when normalized against degree one (more generally,
a proposed fixed degree `d` is short by a factor about `1574/d`). This is
linear-in-`p` behavior, not a fixed-degree phenomenon.
