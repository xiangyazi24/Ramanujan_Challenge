# COINC second-moment numerics report

All finite-field values and all counts below were computed with Python 3
standard-library integer arithmetic. No floating-point arithmetic was used for
a finite-field value or a count. Decimal powers, displayed ratios, and square
roots are the only floating-point quantities. The gap rows were generated from
the recurrence for `N_h(r)` and then divided by
`prod_{j=1}^h (r+j)^3` modulo `p`. As an independent check, every generated
Delta value used at the largest requested windows was compared with the direct
determinant from the `b,c` recurrence; all 594,395 checked values passed.

## T1 - D1 Delta identity

For every requested regular, nonwrapping pair:

```
p= 97: tested= 2385, passed= 2385, failed=0
p=199: tested= 5445, passed= 5445, failed=0
p=499: tested=14445, passed=14445, failed=0
total : tested=22275, passed=22275, failed=0
```

There are no failing `(p,r,h)` triples to list.

The singularity condition is redundant inside the stipulated domain:
`1 <= r < r+h <= p-2` implies `2 <= r+j <= p-2`, so the number of excluded
pairs there is zero.

For completeness, extending `r` to all residues creates exactly
`sum_{h=1}^{30} h = 465` singular `(r,h)` pairs per prime, namely
`r = -j (mod p)` for some `1 <= j <= h`. The denominator vanishes and D1 is
undefined; the regular recurrence cannot be continued through that singular
step. The numerator `N_h(r)` was nonzero at all 465 such pairs for `p=199,499`.
For `p=97`, it was nonzero at 455 pairs and zero (an indeterminate `0/0`) at

```
(h,j) = (26,1), (26,26), (27,2), (27,26), (28,3), (28,26),
        (29,4), (29,26), (30,5), (30,26),  with r=-j mod 97.
```

T1 VERDICT: D1 passes all 22,275 applicable cases. The excluded cases are
genuine poles or occasional `0/0` values, not counterexamples to the regular
identity.

## T2 - exact Parseval bookkeeping

`B_t` was represented exactly as a coefficient vector in the group ring of
`zeta_p`, not as a complex approximation. For each `t=1,...,p-1`, its
coefficient vector had sum `#S_H`. The cyclic autocorrelation represented
`|B_t|^2`. After summing over nonzero `t`, all nonconstant coefficients were
identical, so reduction by `1+zeta+...+zeta^(p-1)=0` gave the exact integer on
the left below.

```
p= 97 H= 5: #S= 460, N_coinc= 2668, LHS= 47196, RHS= 47196, discrepancy=0
p= 97 H=10: #S= 895, N_coinc= 9285, LHS= 99620, RHS= 99620, discrepancy=0
p= 97 H=20: #S=1690, N_coinc=31132, LHS=163704, RHS=163704, discrepancy=0
p=199 H= 5: #S= 970, N_coinc= 5970, LHS=247130, RHS=247130, discrepancy=0
p=199 H=10: #S=1915, N_coinc=20771, LHS=466204, RHS=466204, discrepancy=0
p=199 H=20: #S=3730, N_coinc=74736, LHS=959564, RHS=959564, discrepancy=0
```

Here `RHS = p*N_coinc-(#S_H)^2`. All six exact identities pass.

T2 VERDICT: D3 has correct signs, exclusion of `t=0`, and normalization in all
requested cells.

## T3 - N_coinc structure

The components form the following disjoint partition of ordered coincidence
pairs:

```
diag = identical indices
zz   = distinct zero-zero indices, exactly Z*(Z-1), Z=sum_{h<=H} C_h
same = same-h, off-diagonal, nonzero coincidences, excluding forced pairs
refl = nonzero equality pairs forced by a tested exact reflection
gen  = remaining cross-h, nonzero coincidences
```

For the per-component `K` bookkeeping, the five displayed contributions are

```
diag/#S, zz/#S, same/#S, refl/#S, (gen-(#S)^2/p)/#S.
```

This convention assigns the single uniform baseline to `gen`, and makes the
five numbers add exactly to
`K=(N_coinc-(#S)^2/p)/#S`. Exact component counts and the exact fraction for K
precede the decimal component contributions.

```
p=199 H= 4 S= 778 N=   4158 Z=  4 | diag=   778 zz=   12 same=  952 refl=0 gen=   2416
  K=1.434925 = 222158/154822 | Ki=1.000000 +0.015424 +1.223650 +0.000000 -0.804149
p=199 H= 6 S=1161 N=   8253 Z=  5 | diag=  1161 zz=   20 same= 1344 refl=0 gen=   5728
  K=1.274356 = 294426/231039 | Ki=1.000000 +0.017227 +1.157623 +0.000000 -0.900493
p=199 H=15 S=2835 N=  44027 Z= 17 | diag=  2835 zz=  272 same= 2796 refl=0 gen=  38124
  K=1.283575 = 724148/564165 | Ki=1.000000 +0.095944 +0.986243 +0.000000 -0.798612
p=199 H=33 S=5940 N= 186854 Z= 50 | diag=  5940 zz= 2450 same= 5332 refl=0 gen= 173132
  K=1.607656 = 1900346/1182060 | Ki=1.000000 +0.412458 +0.897643 +0.000000 -0.702445

p=499 H= 5 S= 2470 N=   15124 Z=  6 | diag=  2470 zz=  30 same= 2956 refl=0 gen=    9668
  K=1.173177 = 1445976/1232530 | Ki=1.000000 +0.012146 +1.196761 +0.000000 -1.035730
p=499 H= 8 S= 3940 N=   34648 Z= 12 | diag=  3940 zz= 132 same= 4316 refl=0 gen=   26260
  K=0.898117 = 1765752/1966060 | Ki=1.000000 +0.033503 +1.095431 +0.000000 -1.230817
p=499 H=23 S=11155 N=  260787 Z= 41 | diag= 11155 zz=1640 same=10920 refl=0 gen=  237072
  K=1.023776 = 5698688/5566345 | Ki=1.000000 +0.147019 +0.978933 +0.000000 -1.102177
p=499 H=61 S=28426 N= 1654442 Z= 92 | diag= 28426 zz=8372 same=26568 refl=0 gen= 1591076
  K=1.235785 = 17529082/14184574 | Ki=1.000000 +0.294519 +0.934637 +0.000000 -0.993372

p=997 H= 6 S= 5949 N=   41823 Z=  3 | diag=  5949 zz=    6 same= 6592 refl=0 gen=   29276
  K=1.063356 = 6306930/5931153 | Ki=1.000000 +0.001009 +1.108085 +0.000000 -1.045737
p=997 H=10 S= 9895 N=  107691 Z=  9 | diag=  9895 zz=   72 same=10532 refl=0 gen=   87192
  K=0.958601 = 9456902/9865315 | Ki=1.000000 +0.007276 +1.064376 +0.000000 -1.113051
p=997 H=32 S=31312 N= 1016402 Z= 42 | diag= 31312 zz= 1722 same=31160 refl=0 gen=  952208
  K=1.054244 = 32911450/31218064 | Ki=1.000000 +0.054995 +0.995146 +0.000000 -0.995897
p=997 H=96 S=90864 N= 8380388 Z=136 | diag= 90864 zz=18360 same=86196 refl=0 gen= 8184968
  K=1.092602 = 98980340/90591408 | Ki=1.000000 +0.202060 +0.948627 +0.000000 -1.058085

p=1499 H=  7 S= 10451 N=    81807 Z=  9 | diag= 10451 zz=   72 same=  8860 refl=0 gen=    62424
  K=0.855691 = 13405292/15666049 | Ki=1.000000 +0.006889 +0.847766 +0.000000 -0.998964
p=1499 H= 12 S= 17886 N=   230150 Z= 20 | diag= 17886 zz=  380 same= 16008 refl=0 gen=   195876
  K=0.935651 = 25085854/26811114 | Ki=1.000000 +0.021246 +0.895002 +0.000000 -0.980596
p=1499 H= 39 S= 57603 N=  2276693 Z= 55 | diag= 57603 zz= 2970 same= 55292 refl=0 gen=  2160828
  K=1.096243 = 94657198/86346897 | Ki=1.000000 +0.051560 +0.959881 +0.000000 -0.915197
p=1499 H=125 S=179250 N= 21636754 Z=164 | diag=179250 zz=26732 same=169680 refl=0 gen= 21261092
  K=1.127415 = 302931746/268695750 | Ki=1.000000 +0.149132 +0.946611 +0.000000 -0.968328
```

Reflection scan: affine reversal candidates `r*=A-r-h` were tested for
`A=p-3,...,p+3`. The unique exact candidate was `A=p-1`, and it was exact with
the opposite sign:

```
Delta_{p-1-r-h,h} = -Delta_{r,h}.
p= 199, H= 33: 5940/5940 anti-equalities; 50 equalities, all zero
p= 499, H= 61: 28426/28426 anti-equalities; 92 equalities, all zero
p= 997, H= 96: 90864/90864 anti-equalities; 136 equalities, all zero
p=1499, H=125: 179250/179250 anti-equalities; 164 equalities, all zero
```

Thus reflection creates no nonzero forced equality component. Its zero
equalities are already in `diag` when the reflected index is fixed and in `zz`
otherwise. None of the neighboring affine reversals was an exact equality or
anti-equality.

T3 VERDICT: [COINC] LOOKS TRUE in this range. Across all 16 cells,
`0.855691 <= K <= 1.607656`, with no growth from `p=199` to `p=1499`, including
the `H=ceil(p^0.66)` cells. The persistent positive pieces are the diagonal
(exactly 1 unit of K) and same-h coincidences (about 0.85 to 1.22 units). The
generic cross-h term is correspondingly below its globally assigned uniform
baseline. Reflection is not a threat. The circular zero-zero term is the most
H-sensitive piece and reaches 0.412458 at the smallest prime, but at the
largest windows its contribution falls through 0.412458, 0.294519, 0.202060,
0.149132 as p grows. It is the component to monitor, but these data do not show
it threatening bounded K.

## T4 - difference-curve anomaly scan

The following is the requested upper-triangular matrix. Row `h` lists
`P(h,h')/997` for `h'=h+1,...,20`, rounded to three decimals. All `P` counts
used for thresholding were exact integers.

```
 1: 0.945 0.995 0.973 0.867 1.023 0.911 0.923 0.899 0.979 1.005 1.065 0.899 1.021 1.033 1.085 0.923 1.015 0.949 0.925
 2: 0.945 1.036 0.951 1.016 1.001 1.008 0.939 0.996 0.929 0.956 0.933 0.960 0.989 0.962 1.049 0.982 0.977 0.960
 3: 0.993 1.005 0.923 0.895 0.947 0.999 0.925 1.009 0.905 1.001 1.007 0.933 1.001 1.069 0.987 0.963 1.019
 4: 1.041 0.960 0.955 0.998 1.025 0.926 1.031 0.932 1.015 1.040 0.927 0.948 0.925 1.032 1.061 0.952
 5: 1.013 0.917 0.971 0.987 0.997 1.009 0.931 0.945 1.033 0.951 0.943 0.965 0.937 0.981 0.939
 6: 0.935 0.912 0.953 0.988 0.953 0.950 0.965 0.968 0.949 1.018 0.967 1.004 0.959 0.954
 7: 1.005 1.043 1.047 1.005 0.969 1.011 0.945 0.955 1.013 0.993 0.941 0.963 0.971
 8: 0.973 1.048 1.013 0.966 0.941 0.918 0.955 0.986 0.967 1.044 1.035 0.962
 9: 0.975 0.995 0.961 1.015 0.997 1.029 1.005 0.949 0.911 0.977 0.895
10: 1.031 0.910 0.993 0.930 0.999 1.000 0.979 0.940 0.967 0.954
11: 0.955 0.961 0.911 0.935 0.983 0.967 0.935 1.035 0.959
12: 1.025 0.972 0.977 1.032 0.939 0.952 0.965 0.996
13: 1.025 0.889 0.957 0.973 0.979 0.975 1.005
14: 1.039 0.986 1.007 1.006 0.929 1.000
15: 1.063 0.887 0.935 1.007 1.003
16: 0.991 0.972 0.913 1.070
17: 1.021 0.951 0.995
18: 0.949 0.994
19: 0.951
```

There are no flagged cross-gap pairs: none has `P>2.5p` or `P<0.3p`. The
minimum is `P(1,5)=864=0.866600p`; the maximum is
`P(1,16)=1082=1.085256p`. Consequently there is no flagged parity,
divisibility, or `h+h'` explanation to make, and no visible grouping by those
features at this scale.

Same-h off-diagonal counts, after subtracting the diagonal `r=r'`, are:

```
h= 1: P=1984, P/p=1.989970    h=11: P= 942, P/p=0.944835
h= 2: P= 920, P/p=0.922768    h=12: P= 984, P/p=0.986961
h= 3: P= 932, P/p=0.934804    h=13: P= 996, P/p=0.998997
h= 4: P= 884, P/p=0.886660    h=14: P= 912, P/p=0.914744
h= 5: P= 956, P/p=0.958877    h=15: P= 992, P/p=0.994985
h= 6: P= 916, P/p=0.918756    h=16: P=1052, P/p=1.055165
h= 7: P= 934, P/p=0.936810    h=17: P= 970, P/p=0.972919
h= 8: P=1016, P/p=1.019057    h=18: P= 950, P/p=0.952859
h= 9: P=1004, P/p=1.007021    h=19: P= 924, P/p=0.926780
h=10: P= 994, P/p=0.996991    h=20: P= 978, P/p=0.980943
```

The same-h `h=1` value is the sole conspicuous entry. It is explained by
`Delta_{r,1}=(r+1)^(-3)` and `997-1` being divisible by 3, so the cube map on
`F_997^*` has three-element fibers (up to the two window exclusions).

T4 VERDICT: the cross-gap family is strikingly generic for `h,h'<=20`; all 190
ratios lie in `[0.866600,1.085256]`, far inside the requested anomaly band.

## T5 - cumulative C_h profile

At each decile, the fields are
`H : sumC, sumC/H, sumC/(H+sqrt(pH)), H+sqrt(pH)`. The square-root curve is
shown to three decimals; `sumC` is exact.

```
p=997, Hmax=floor(p^0.66)=95, total sumC=135
 10:   9, 0.900, 0.082, 109.850
 19:  23, 1.211, 0.147, 156.634
 29:  38, 1.310, 0.191, 199.038
 38:  49, 1.289, 0.211, 232.643
 48:  62, 1.292, 0.232, 266.760
 57:  72, 1.263, 0.244, 295.388
 67:  91, 1.358, 0.280, 325.455
 76: 108, 1.421, 0.307, 351.267
 86: 123, 1.430, 0.325, 378.817
 95: 135, 1.421, 0.335, 402.758

p=1499, Hmax=floor(p^0.66)=124, total sumC=164
 13:  20, 1.538, 0.131, 152.596
 25:  44, 1.760, 0.201, 218.585
 38:  55, 1.447, 0.199, 276.667
 50:  67, 1.340, 0.207, 323.770
 62:  79, 1.274, 0.215, 366.857
 75:  93, 1.240, 0.227, 410.298
 87: 111, 1.276, 0.248, 448.127
100: 124, 1.240, 0.255, 487.169
112: 144, 1.286, 0.276, 521.741
124: 164, 1.323, 0.295, 555.133

p=1999, Hmax=floor(p^0.66)=150, total sumC=213
 15:  25, 1.667, 0.133, 188.162
 30:  57, 1.900, 0.207, 274.888
 45:  72, 1.600, 0.209, 344.925
 60:  94, 1.567, 0.231, 406.324
 75: 115, 1.533, 0.249, 462.201
 90: 133, 1.478, 0.259, 514.158
105: 154, 1.467, 0.273, 563.143
120: 184, 1.533, 0.302, 609.775
135: 199, 1.474, 0.304, 654.485
150: 213, 1.420, 0.305, 697.586
```

As a scale-free tracking diagnostic, the coefficient of variation of the ten
decile ratios was:

```
p= 997: CV(sumC/H)=0.115, CV(sumC/(H+sqrt(pH)))=0.327
p=1499: CV(sumC/H)=0.115, CV(sumC/(H+sqrt(pH)))=0.196
p=1999: CV(sumC/H)=0.084, CV(sumC/(H+sqrt(pH)))=0.209
```

T5 VERDICT: `sum_{h<=H} C_h` tracks `H`, not `H+sqrt(pH)`. The ratio to H is
roughly flat, whereas the ratio to the square-root curve rises systematically
through the deciles. At the largest windows the exact totals are only 135,
164, and 213, versus square-root-curve values 402.758, 555.133, and 697.586.

## Overall verdict

D1 and D3 pass exactly. In the measured range, [COINC] is strongly consistent
with a bounded absolute K near 1. The proposed reflection does not add nonzero
coincidences because it is an anti-symmetry. The only structurally sensitive
term is zero-zero, but its normalized contribution decreases with p at the
largest windows. Both the difference-curve scan and the direct `C_h` profile
support generic cross-gap behavior and cumulative root counts of order H.
