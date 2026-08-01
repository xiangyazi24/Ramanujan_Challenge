# Centered Apéry pair dispersion by shift scale

Date: 2026-08-01

## Verdict

The total centered pair sum is small at every scanned scale:

\[
  -0.0731 \le \frac{D(X)}{X\lambda_X^2} \le 0.1633,
  \qquad X=2^{13},\ldots,2^{19}.
\]

That total conceals a stable shell pattern.  The shell
`X/2 <= |p-q| < X` is negative at all seven scales, and its deviation grows
from `-7.94` at `X=8192` to `-177.00` at `X=524288`.  At the three largest
scales its size is `6.79`, `12.13`, and `14.36` times the square root of its
unsigned intersection mass.  Positive deviations in intermediate shells
cancel most of this negative edge-shell contribution.  Thus the raw centered
sum is **not shellwise square-root-cancelling** at large `X`, even though its
total Poisson normalization stays small.

There is no small-shift concentration and no separate reflection diagonal in
this range.  The persistent large-shift deficit is consistent with the
deterministic support geometry of the truncated rows `A_{p,X}` rather than a
new cross-prime algebraic correlation.  A dispersion argument using this
centering would need to subtract or otherwise absorb that boundary profile.

## Definitions and conventions

For ordered distinct prime pairs, let

\[
\begin{aligned}
 C_j(X)&=\sum_{\substack{p\ne q\\2^j\le |p-q|<2^{j+1}}}
              |A_{p,X}\cap A_{q,X}|,\\
 W_j(X)&=\sum_{\substack{p\ne q\\2^j\le |p-q|<2^{j+1}}}
              |A_{p,X}|\,|A_{q,X}|,\\
 E_j(X)&=W_j(X)/X,\qquad D_j(X)=C_j(X)-E_j(X).
\end{aligned}
\]

Thus `C_j` is both the unsigned intersection mass and the count of ordered
contributing `(p,q,n)` triples.  The diagnostic `z_j` below is
`|D_j|/sqrt(C_j)` and is left blank when `C_j=0`.  `R_j` counts the ordered
triples in which both residues lie in non-fixed reflection orbits.  `F_j`
counts the stricter folded-coordinate diagonal

\[
  \min(r,p-1-r)=\min(s,q-1-s).
\]

The implementation also counts the fixed-fixed channel
`2r=p-1` and `2s=q-1`; it is zero throughout.  Since the only midpoint-zero
primes below one million are `11` and `3137`, both outside every reported
prime window, all reported intersections have `R_j=C_j`.  Every `F_j` is
zero.  Consequently reflection symmetry duplicates the row positions but
does not isolate the observed large-shift bias.

## Totals and cancellation

Here `S_X=sum_p |A_{p,X}|`, `lambda_X=S_X/X`, and
`L1=sum_j |D_j|`.  The cancellation factor is `L1/|D|`.

| X | S_X | lambda_X | M_2=C | E | D | D/(X lambda_X^2) | L1 | cancellation factor |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8,192 | 640 | 0.078125000 | 58 | 49.836914 | +8.163086 | +0.163262 | 29.015625 | 3.554 |
| 16,384 | 1,222 | 0.074584961 | 100 | 90.981323 | +9.018677 | +0.098951 | 25.923706 | 2.874 |
| 32,768 | 2,226 | 0.067932129 | 154 | 151.067566 | +2.932434 | +0.019392 | 33.960754 | 11.581 |
| 65,536 | 3,969 | 0.060562134 | 238 | 240.246033 | -2.246033 | -0.009344 | 53.092285 | 23.638 |
| 131,072 | 7,802 | 0.059524536 | 462 | 464.279999 | -2.279999 | -0.004909 | 132.632690 | 58.172 |
| 262,144 | 14,487 | 0.055263519 | 742 | 800.484581 | -58.484581 | -0.073051 | 237.878899 | 4.067 |
| 524,288 | 27,489 | 0.052431107 | 1,458 | 1,441.165592 | +16.834408 | +0.011680 | 424.194279 | 25.198 |

## Per-scale verdicts

- `X=8192`: the `[X/2,X)` shell is negative at `3.97 sqrt(C_j)`, but
  `C_j=4`; this is a finite-sample warning rather than a stable verdict.
- `X=16384`: every nonempty shell is within `1.71 sqrt(C_j)`; consistent with
  square-root fluctuation at this scale.
- `X=32768`: every nonempty shell is within `1.52 sqrt(C_j)`; consistent with
  square-root fluctuation.
- `X=65536`: the `[X/2,X)` shell is `-17.34`, or `2.74 sqrt(C_j)`; a mild
  large-shift deficit appears, with no small-shift bias.
- `X=131072`: the `[X/2,X)` shell is `-49.91` (`6.79 sqrt(C_j)`), partly
  cancelled by `+26.12` in `[X/32,X/16)`.  This is structured shell bias.
- `X=262144`: the `[X/2,X)` shell is `-105.74` (`12.13 sqrt(C_j)`), while
  `[X/16,X/8)` contributes `+43.57`.  Shellwise square-root cancellation
  fails strongly.
- `X=524288`: the `[X/2,X)` shell is `-177.00` (`14.36 sqrt(C_j)`), offset by
  several positive intermediate shells, most notably `+57.10` in
  `[X/32,X/16)`.  The small total `D` is cancellation between structured
  shells, not uniform shellwise randomness.

## Complete shell tables

### X = 8,192

| j | shift shell | C_j | W_j | E_j | D_j | z_j | R_j | F_j |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | [1, 2) | 0 | 0 | 0.000000 | +0.000000 | -- | 0 | 0 |
| 1 | [2, 4) | 0 | 196 | 0.023926 | -0.023926 | -- | 0 | 0 |
| 2 | [4, 8) | 0 | 388 | 0.047363 | -0.047363 | -- | 0 | 0 |
| 3 | [8, 16) | 0 | 744 | 0.090820 | -0.090820 | -- | 0 | 0 |
| 4 | [16, 32) | 0 | 1,450 | 0.177002 | -0.177002 | -- | 0 | 0 |
| 5 | [32, 64) | 0 | 3,012 | 0.367676 | -0.367676 | -- | 0 | 0 |
| 6 | [64, 128) | 0 | 5,720 | 0.698242 | -0.698242 | -- | 0 | 0 |
| 7 | [128, 256) | 2 | 11,636 | 1.420410 | +0.579590 | 0.410 | 2 | 0 |
| 8 | [256, 512) | 10 | 22,614 | 2.760498 | +7.239502 | 2.289 | 10 | 0 |
| 9 | [512, 1,024) | 12 | 45,748 | 5.584473 | +6.415527 | 1.852 | 12 | 0 |
| 10 | [1,024, 2,048) | 12 | 85,196 | 10.399902 | +1.600098 | 0.462 | 12 | 0 |
| 11 | [2,048, 4,096) | 18 | 124,890 | 15.245361 | +2.754639 | 0.649 | 18 | 0 |
| 12 | [4,096, 8,192) | 4 | 97,774 | 11.935303 | -7.935303 | 3.968 | 4 | 0 |
| 13 | [8,192, 16,384) | 0 | 8,896 | 1.085938 | -1.085938 | -- | 0 | 0 |

### X = 16,384

| j | shift shell | C_j | W_j | E_j | D_j | z_j | R_j | F_j |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | [1, 2) | 0 | 0 | 0.000000 | +0.000000 | -- | 0 | 0 |
| 1 | [2, 4) | 0 | 264 | 0.016113 | -0.016113 | -- | 0 | 0 |
| 2 | [4, 8) | 0 | 822 | 0.050171 | -0.050171 | -- | 0 | 0 |
| 3 | [8, 16) | 0 | 1,520 | 0.092773 | -0.092773 | -- | 0 | 0 |
| 4 | [16, 32) | 0 | 3,082 | 0.188110 | -0.188110 | -- | 0 | 0 |
| 5 | [32, 64) | 0 | 5,326 | 0.325073 | -0.325073 | -- | 0 | 0 |
| 6 | [64, 128) | 4 | 11,464 | 0.699707 | +3.300293 | 1.650 | 4 | 0 |
| 7 | [128, 256) | 0 | 23,382 | 1.427124 | -1.427124 | -- | 0 | 0 |
| 8 | [256, 512) | 2 | 43,758 | 2.670776 | -0.670776 | 0.474 | 2 | 0 |
| 9 | [512, 1,024) | 6 | 88,472 | 5.399902 | +0.600098 | 0.245 | 6 | 0 |
| 10 | [1,024, 2,048) | 18 | 176,170 | 10.752563 | +7.247437 | 1.708 | 18 | 0 |
| 11 | [2,048, 4,096) | 20 | 320,260 | 19.547119 | +0.452881 | 0.101 | 20 | 0 |
| 12 | [4,096, 8,192) | 34 | 460,874 | 28.129517 | +5.870483 | 1.007 | 34 | 0 |
| 13 | [8,192, 16,384) | 16 | 330,846 | 20.193237 | -4.193237 | 1.048 | 16 | 0 |
| 14 | [16,384, 32,768) | 0 | 24,398 | 1.489136 | -1.489136 | -- | 0 | 0 |

### X = 32,768

| j | shift shell | C_j | W_j | E_j | D_j | z_j | R_j | F_j |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | [1, 2) | 0 | 0 | 0.000000 | +0.000000 | -- | 0 | 0 |
| 1 | [2, 4) | 0 | 446 | 0.013611 | -0.013611 | -- | 0 | 0 |
| 2 | [4, 8) | 0 | 1,232 | 0.037598 | -0.037598 | -- | 0 | 0 |
| 3 | [8, 16) | 0 | 2,140 | 0.065308 | -0.065308 | -- | 0 | 0 |
| 4 | [16, 32) | 0 | 4,444 | 0.135620 | -0.135620 | -- | 0 | 0 |
| 5 | [32, 64) | 0 | 8,394 | 0.256165 | -0.256165 | -- | 0 | 0 |
| 6 | [64, 128) | 2 | 18,760 | 0.572510 | +1.427490 | 1.009 | 2 | 0 |
| 7 | [128, 256) | 2 | 36,150 | 1.103210 | +0.896790 | 0.634 | 2 | 0 |
| 8 | [256, 512) | 2 | 75,066 | 2.290833 | -0.290833 | 0.206 | 2 | 0 |
| 9 | [512, 1,024) | 4 | 147,714 | 4.507874 | -0.507874 | 0.254 | 4 | 0 |
| 10 | [1,024, 2,048) | 6 | 291,848 | 8.906494 | -2.906494 | 1.187 | 6 | 0 |
| 11 | [2,048, 4,096) | 22 | 570,740 | 17.417603 | +4.582397 | 0.977 | 22 | 0 |
| 12 | [4,096, 8,192) | 28 | 1,047,814 | 31.976746 | -3.976746 | 0.752 | 28 | 0 |
| 13 | [8,192, 16,384) | 58 | 1,522,404 | 46.460083 | +11.539917 | 1.515 | 58 | 0 |
| 14 | [16,384, 32,768) | 30 | 1,131,076 | 34.517700 | -4.517700 | 0.825 | 30 | 0 |
| 15 | [32,768, 65,536) | 0 | 91,954 | 2.806213 | -2.806213 | -- | 0 | 0 |

### X = 65,536

| j | shift shell | C_j | W_j | E_j | D_j | z_j | R_j | F_j |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | [1, 2) | 0 | 0 | 0.000000 | +0.000000 | -- | 0 | 0 |
| 1 | [2, 4) | 0 | 646 | 0.009857 | -0.009857 | -- | 0 | 0 |
| 2 | [4, 8) | 0 | 1,836 | 0.028015 | -0.028015 | -- | 0 | 0 |
| 3 | [8, 16) | 0 | 3,208 | 0.048950 | -0.048950 | -- | 0 | 0 |
| 4 | [16, 32) | 4 | 7,470 | 0.113983 | +3.886017 | 1.943 | 4 | 0 |
| 5 | [32, 64) | 0 | 13,350 | 0.203705 | -0.203705 | -- | 0 | 0 |
| 6 | [64, 128) | 0 | 27,996 | 0.427185 | -0.427185 | -- | 0 | 0 |
| 7 | [128, 256) | 2 | 55,136 | 0.841309 | +1.158691 | 0.819 | 2 | 0 |
| 8 | [256, 512) | 0 | 113,468 | 1.731384 | -1.731384 | -- | 0 | 0 |
| 9 | [512, 1,024) | 2 | 226,376 | 3.454224 | -1.454224 | 1.028 | 2 | 0 |
| 10 | [1,024, 2,048) | 12 | 454,074 | 6.928619 | +5.071381 | 1.464 | 12 | 0 |
| 11 | [2,048, 4,096) | 18 | 901,510 | 13.755951 | +4.244049 | 1.000 | 18 | 0 |
| 12 | [4,096, 8,192) | 38 | 1,765,344 | 26.937012 | +11.062988 | 1.795 | 38 | 0 |
| 13 | [8,192, 16,384) | 48 | 3,253,554 | 49.645294 | -1.645294 | 0.237 | 48 | 0 |
| 14 | [16,384, 32,768) | 74 | 4,870,450 | 74.317169 | -0.317169 | 0.037 | 74 | 0 |
| 15 | [32,768, 65,536) | 40 | 3,757,518 | 57.335175 | -17.335175 | 2.741 | 40 | 0 |
| 16 | [65,536, 131,072) | 0 | 292,828 | 4.468201 | -4.468201 | -- | 0 | 0 |

### X = 131,072

| j | shift shell | C_j | W_j | E_j | D_j | z_j | R_j | F_j |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | [1, 2) | 0 | 0 | 0.000000 | +0.000000 | -- | 0 | 0 |
| 1 | [2, 4) | 0 | 1,334 | 0.010178 | -0.010178 | -- | 0 | 0 |
| 2 | [4, 8) | 0 | 3,534 | 0.026962 | -0.026962 | -- | 0 | 0 |
| 3 | [8, 16) | 0 | 6,350 | 0.048447 | -0.048447 | -- | 0 | 0 |
| 4 | [16, 32) | 0 | 14,942 | 0.113998 | -0.113998 | -- | 0 | 0 |
| 5 | [32, 64) | 0 | 27,038 | 0.206284 | -0.206284 | -- | 0 | 0 |
| 6 | [64, 128) | 4 | 57,282 | 0.437027 | +3.562973 | 1.781 | 4 | 0 |
| 7 | [128, 256) | 8 | 115,576 | 0.881775 | +7.118225 | 2.517 | 8 | 0 |
| 8 | [256, 512) | 2 | 229,608 | 1.751770 | +0.248230 | 0.176 | 2 | 0 |
| 9 | [512, 1,024) | 6 | 461,942 | 3.524338 | +2.475662 | 1.011 | 6 | 0 |
| 10 | [1,024, 2,048) | 14 | 919,280 | 7.013550 | +6.986450 | 1.867 | 14 | 0 |
| 11 | [2,048, 4,096) | 10 | 1,839,924 | 14.037506 | -4.037506 | 1.277 | 10 | 0 |
| 12 | [4,096, 8,192) | 54 | 3,654,558 | 27.882065 | +26.117935 | 3.554 | 54 | 0 |
| 13 | [8,192, 16,384) | 62 | 7,127,872 | 54.381348 | +7.618652 | 0.968 | 62 | 0 |
| 14 | [16,384, 32,768) | 110 | 12,969,808 | 98.951782 | +11.048218 | 1.053 | 110 | 0 |
| 15 | [32,768, 65,536) | 138 | 18,772,294 | 143.221237 | -5.221237 | 0.444 | 138 | 0 |
| 16 | [65,536, 131,072) | 54 | 13,619,086 | 103.905380 | -49.905380 | 6.791 | 54 | 0 |
| 17 | [131,072, 262,144) | 0 | 1,033,680 | 7.886353 | -7.886353 | -- | 0 | 0 |

### X = 262,144

| j | shift shell | C_j | W_j | E_j | D_j | z_j | R_j | F_j |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | [1, 2) | 0 | 0 | 0.000000 | +0.000000 | -- | 0 | 0 |
| 1 | [2, 4) | 0 | 2,024 | 0.007721 | -0.007721 | -- | 0 | 0 |
| 2 | [4, 8) | 0 | 6,170 | 0.023537 | -0.023537 | -- | 0 | 0 |
| 3 | [8, 16) | 0 | 11,262 | 0.042961 | -0.042961 | -- | 0 | 0 |
| 4 | [16, 32) | 0 | 24,778 | 0.094521 | -0.094521 | -- | 0 | 0 |
| 5 | [32, 64) | 0 | 46,110 | 0.175896 | -0.175896 | -- | 0 | 0 |
| 6 | [64, 128) | 0 | 99,114 | 0.378090 | -0.378090 | -- | 0 | 0 |
| 7 | [128, 256) | 4 | 193,408 | 0.737793 | +3.262207 | 1.631 | 4 | 0 |
| 8 | [256, 512) | 2 | 395,000 | 1.506805 | +0.493195 | 0.349 | 2 | 0 |
| 9 | [512, 1,024) | 6 | 781,620 | 2.981644 | +3.018356 | 1.232 | 6 | 0 |
| 10 | [1,024, 2,048) | 2 | 1,565,400 | 5.971527 | -3.971527 | 2.808 | 2 | 0 |
| 11 | [2,048, 4,096) | 16 | 3,144,650 | 11.995888 | +4.004112 | 1.001 | 16 | 0 |
| 12 | [4,096, 8,192) | 34 | 6,244,872 | 23.822296 | +10.177704 | 1.745 | 34 | 0 |
| 13 | [8,192, 16,384) | 58 | 12,398,898 | 47.298042 | +10.701958 | 1.405 | 58 | 0 |
| 14 | [16,384, 32,768) | 136 | 24,229,826 | 92.429451 | +43.570549 | 3.736 | 136 | 0 |
| 15 | [32,768, 65,536) | 184 | 44,441,514 | 169.530922 | +14.469078 | 1.067 | 184 | 0 |
| 16 | [65,536, 131,072) | 224 | 64,895,798 | 247.557823 | -23.557823 | 1.574 | 224 | 0 |
| 17 | [131,072, 262,144) | 76 | 47,641,486 | 181.737846 | -105.737846 | 12.129 | 76 | 0 |
| 18 | [262,144, 524,288) | 0 | 3,720,300 | 14.191818 | -14.191818 | -- | 0 | 0 |

### X = 524,288

| j | shift shell | C_j | W_j | E_j | D_j | z_j | R_j | F_j |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 0 | [1, 2) | 0 | 0 | 0.000000 | +0.000000 | -- | 0 | 0 |
| 1 | [2, 4) | 0 | 3,434 | 0.006550 | -0.006550 | -- | 0 | 0 |
| 2 | [4, 8) | 0 | 11,102 | 0.021175 | -0.021175 | -- | 0 | 0 |
| 3 | [8, 16) | 2 | 19,560 | 0.037308 | +1.962692 | 1.388 | 2 | 0 |
| 4 | [16, 32) | 0 | 44,886 | 0.085613 | -0.085613 | -- | 0 | 0 |
| 5 | [32, 64) | 0 | 84,096 | 0.160400 | -0.160400 | -- | 0 | 0 |
| 6 | [64, 128) | 0 | 176,198 | 0.336071 | -0.336071 | -- | 0 | 0 |
| 7 | [128, 256) | 0 | 347,326 | 0.662472 | -0.662472 | -- | 0 | 0 |
| 8 | [256, 512) | 4 | 705,016 | 1.344711 | +2.655289 | 1.328 | 4 | 0 |
| 9 | [512, 1,024) | 12 | 1,396,576 | 2.663757 | +9.336243 | 2.695 | 12 | 0 |
| 10 | [1,024, 2,048) | 8 | 2,808,054 | 5.355938 | +2.644062 | 0.935 | 8 | 0 |
| 11 | [2,048, 4,096) | 20 | 5,615,048 | 10.709854 | +9.290146 | 2.077 | 20 | 0 |
| 12 | [4,096, 8,192) | 22 | 11,217,816 | 21.396286 | +0.603714 | 0.129 | 22 | 0 |
| 13 | [8,192, 16,384) | 74 | 22,398,412 | 42.721581 | +31.278419 | 3.636 | 74 | 0 |
| 14 | [16,384, 32,768) | 142 | 44,513,764 | 84.903267 | +57.096733 | 4.791 | 142 | 0 |
| 15 | [32,768, 65,536) | 210 | 86,934,226 | 165.813877 | +44.186123 | 3.049 | 210 | 0 |
| 16 | [65,536, 131,072) | 338 | 159,484,154 | 304.191883 | +33.808117 | 1.839 | 338 | 0 |
| 17 | [131,072, 262,144) | 474 | 234,014,478 | 446.347195 | +27.652805 | 1.270 | 474 | 0 |
| 18 | [262,144, 524,288) | 152 | 172,491,308 | 329.001060 | -177.001060 | 14.357 | 152 | 0 |
| 19 | [524,288, 1,048,576) | 0 | 13,320,372 | 25.406593 | -25.406593 | -- | 0 | 0 |

## Sanity gates

All required gates passed and abort the program on failure.

1. The transformed Montgomery recurrence agreed exactly, zero by zero, with
   the original Apéry recurrence using an inverse table at
   `p=7,11,17,31,181,379,3137`.
2. The hit lists reproduced `M_2(4000)=18`, `M_2(8000)=54`,
   `M_2(8192)=58`, and `M_2(524288)=1458`.  The five intervening dyadic
   `M_2` values also agreed with `CRON_garqi_moments_report.md`.
3. For every `X`, the exact integer identities

   \[
     \sum_j C_j=M_2,
     \quad \sum_j W_j=W,
     \quad \sum_j(C_jX-W_j)+W=M_2X
   \]

   closed exactly.
4. Five deterministic nonempty prime pairs at `X=8192` were counted both
   from the hit lists and directly from the definition; all matched:

   | p | q | witness n | hit-list intersection | direct intersection |
   |---:|---:|---:|---:|---:|
   | 7,559 | 7,937 | 8,768 | 1 | 1 |
   | 6,827 | 8,293 | 8,959 | 1 | 1 |
   | 6,781 | 7,433 | 9,587 | 1 | 1 |
   | 6,379 | 7,333 | 9,810 | 1 | 1 |
   | 5,573 | 6,151 | 9,841 | 1 | 1 |

5. The 81,950 zero records were sorted, duplicate-free, in range, and closed
   under `r -> p-1-r` before any dispersion statistic was accepted.

## Implementation and reproduction

`CRON_pair_dispersion.c` scans every prime `7 <= p <= 2^20` with eight
pthread workers.  It retains the 81,950 zero pairs, builds the actual prime
list for each hit `n`, and enumerates only the resulting collisions.  The
benchmark uses prefix sums of `|A_{p,X}|` over sorted primes for each dyadic
shift shell, costing `O(pi(2X) log X)` rather than looping over roughly three
billion prime pairs.

The final arm64 Mac mini run processed 41,162,092,072 recurrence steps in
220.329 seconds; all analysis and gates brought the total to 222.152 seconds.
The scan throughput was `1.86821e8` recurrence steps per second.  The exact
source SHA-256 for that run was

```text
3c651ed921aded1d942f40c274f993d61178228bbb3eac51c595a6b50c3b5544
```

From the repository root:

```sh
cc -O3 -march=native -Wall -Wextra -Wpedantic \
  -o /tmp/CRON_pair_dispersion \
  problems/3.2/CRON_pair_dispersion.c -lpthread -lm
/tmp/CRON_pair_dispersion 8
```
