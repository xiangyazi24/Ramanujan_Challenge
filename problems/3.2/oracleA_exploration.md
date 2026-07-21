# Oracle A: top-block Fourier exploration

This is a reproducible diagnostic for the analytic oracle. It is evidence, not a proof of `eq:short-arc` or AMTD.

## Definitions and reproducibility

The data source is `problems/3.2/data_zp_pairs.bin` (149,112 records, SHA-256 `8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d`). Its 2,565 stored nonempty zero sets passed sortedness, reflection, and no-consecutive-zero checks. The analyzed prime block is always `N/2 < p <= N`; the default scales are `8192, 16384, 32768, 65536`.

For a doublet `Z_p={r_p,p-1-r_p}`, with `r_p<(p-1)/2`, set

```text
h_p  = p-1-2r_p
m_1  = p+r_p              = (3p-1-h_p)/2
m_2  = 2p-1-r_p           = (3p-1+h_p)/2
S(theta) = sum_p [e(theta m_1)+e(theta m_2)]
         = 2 sum_p e(theta(3p-1)/2) cos(pi theta h_p).
```

The FFT grid has size the next power of two at least `32N`. Because it is longer than the diameter of the integer hit set, sampled Parseval is alias-free and equals the exact integral up to floating-point roundoff. A sampled point is called a major-arc point when `|theta-a/q| <= 16/(qN)` for some reduced `a/q` with `q <= 16`; this cutoff is only a declared diagnostic convention.

The true covariance experiment is separate: on `I_N=(N,2N]`, it uses `Omega_p={m in I_N: m mod p in Z_p}` and both unit and `log p` weights. Every lift `r+kp` in the shell is included; because `p>N/2`, each residue has one or two lifts. Its exact finite Fourier comparison is

```text
E^o = sum_{p!=q} w_p w_q (|Omega_p cap Omega_q|-A_p A_q/N)
    = (1/N) sum_{k=1}^{N-1} (|sum_p w_p F_p(k)|^2
                                  -sum_p w_p^2 |F_p(k)|^2),
F_p(k)=sum_{m in Omega_p} e(-k(m-N-1)/N).
```

Run:

```sh
python3 problems/3.2/oracleA_explore.py
```

Runtime dependencies: Python `3.9.6`, NumPy `1.26.4`. The report used `12` visible CPUs, although NumPy decides its own FFT threading.

## Doublet Fourier norm and minor-arc maxima

Here `T=2K`, `F=sum_m lambda(m)(lambda(m)-1)`, and the exact identity is `integral |S|^2 = T+F`. The minor maximum is over sampled local maxima outside the declared major arcs.

| N | K | complete T | complete positions in I_N | all doublet shell lifts | F | (T+F)/T | FFT grid | Parseval rel. error | max minor theta | max minor / sqrt(T) | nearest a/q | N offset |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|
| 8192 | 134 | 268 | 162 | 374 | 6 | 1.022388 | 262,144 | 0.000e+00 | 0.183731079 | 3.3760 | 2/11 | +15.670 |
| 16384 | 282 | 564 | 339 | 775 | 18 | 1.031915 | 524,288 | 1.953e-15 | 0.061416626 | 3.4017 | 1/16 | -17.750 |
| 32768 | 509 | 1018 | 625 | 1408 | 32 | 1.031434 | 1,048,576 | 1.126e-14 | 0.170152664 | 3.4349 | 1/6 | +114.229 |
| 65536 | 877 | 1754 | 1072 | 2430 | 48 | 1.027366 | 2,097,152 | 2.271e-15 | 0.102297306 | 3.6885 | 1/10 | +150.556 |

Top `8` sampled local maxima on `0 <= theta <= 1/2` at each scale:

### Peaks for N=8192

| theta | |S| | |S|/sqrt(T) | arc | nearest a/q | N offset |
|---:|---:|---:|:---:|:---:|---:|
| 0.000000000 | 268.0000 | 16.3707 | major | 0/1 | +0.000 |
| 0.183731079 | 55.2673 | 3.3760 | minor | 2/11 | +15.670 |
| 0.134159088 | 51.9505 | 3.1734 | minor | 2/15 | +6.765 |
| 0.024024963 | 49.4102 | 3.0182 | minor | 0/1 | +196.812 |
| 0.198791504 | 47.0597 | 2.8746 | minor | 1/5 | -9.900 |
| 0.316894531 | 46.4583 | 2.8379 | minor | 5/16 | +36.000 |
| 0.106693268 | 46.1395 | 2.8184 | minor | 1/9 | -36.191 |
| 0.234798431 | 45.4684 | 2.7774 | minor | 3/13 | +33.007 |

### Peaks for N=16384

| theta | |S| | |S|/sqrt(T) | arc | nearest a/q | N offset |
|---:|---:|---:|:---:|:---:|---:|
| 0.000000000 | 564.0000 | 23.7487 | major | 0/1 | +0.000 |
| 0.061416626 | 80.7856 | 3.4017 | minor | 1/16 | -17.750 |
| 0.304628372 | 76.8433 | 3.2357 | minor | 4/13 | -50.200 |
| 0.285606384 | 68.5556 | 2.8867 | major | 2/7 | -1.768 |
| 0.331329346 | 66.5339 | 2.8016 | minor | 1/3 | -32.833 |
| 0.192657471 | 66.1307 | 2.7846 | minor | 3/16 | +84.500 |
| 0.132814407 | 66.0366 | 2.7806 | minor | 2/15 | -8.502 |
| 0.296766281 | 65.9226 | 2.7758 | minor | 3/10 | -52.981 |

### Peaks for N=32768

| theta | |S| | |S|/sqrt(T) | arc | nearest a/q | N offset |
|---:|---:|---:|:---:|:---:|---:|
| 0.000000000 | 1018.0000 | 31.9061 | major | 0/1 | +0.000 |
| 0.170152664 | 109.5946 | 3.4349 | minor | 1/6 | +114.229 |
| 0.060119629 | 108.6450 | 3.4051 | minor | 1/16 | -78.000 |
| 0.271627426 | 107.0572 | 3.3554 | minor | 3/11 | -36.040 |
| 0.207726479 | 106.5766 | 3.3403 | minor | 3/14 | -214.933 |
| 0.198372841 | 103.8443 | 3.2547 | minor | 1/5 | -53.319 |
| 0.143852234 | 101.5070 | 3.1814 | minor | 1/7 | +32.607 |
| 0.212293625 | 101.3322 | 3.1759 | minor | 3/14 | -65.277 |

### Peaks for N=65536

| theta | |S| | |S|/sqrt(T) | arc | nearest a/q | N offset |
|---:|---:|---:|:---:|:---:|---:|
| 0.000000000 | 1754.0000 | 41.8808 | major | 0/1 | +0.000 |
| 0.102297306 | 154.4787 | 3.6885 | minor | 1/10 | +150.556 |
| 0.294124603 | 134.3722 | 3.2084 | minor | 3/10 | -385.050 |
| 0.140850544 | 134.3383 | 3.2076 | minor | 1/7 | -131.504 |
| 0.483538628 | 134.2832 | 3.2063 | minor | 1/2 | -1078.812 |
| 0.492629051 | 133.1225 | 3.1786 | minor | 1/2 | -483.062 |
| 0.396017551 | 132.1390 | 3.1551 | minor | 2/5 | -260.994 |
| 0.246010303 | 131.9947 | 3.1517 | minor | 1/4 | -261.469 |

## Small-denominator gap/phase correlations

At each reduced `theta=a/q`, `0<a<q<=16`, the relevant statistic is the connected complex correlation `Corr(e(3p theta/2), cos(pi theta h_p))`. The raw-gap column instead uses `h_p/p`. Cases where the prime phase or cosine is constant are marked degenerate and excluded from maxima. `|S|/sqrt(T)` is also recorded because a large value can come from a nonzero marginal even when connected correlation is small.

| N | max modulation corr. | theta | median modulation corr. | max raw-gap corr. | theta | max rational |S|/sqrt(T) | theta |
|---:|---:|:---:|---:|---:|:---:|---:|:---:|
| 8192 | 0.2204 | 11/16 | 0.0666 | 0.1494 | 8/9 | 2.4331 | 11/16 |
| 16384 | 0.1412 | 3/4 | 0.0422 | 0.1229 | 7/9 | 2.4086 | 3/4 |
| 32768 | 0.0905 | 1/4 | 0.0378 | 0.0843 | 9/14 | 1.9793 | 1/4 |
| 65536 | 0.0657 | 13/14 | 0.0290 | 0.0539 | 3/13 | 1.9361 | 13/14 |

For the largest scale `N=65536`, maxima over reduced numerators at each denominator are:

| q | max modulation corr. (a/q) | max raw-gap corr. (a/q) | max |S|/sqrt(T) (a/q) |
|---:|:---|:---|:---|
| 2 | 0.0034 (1/2) | 0.0526 (1/2) | 0.1433 (1/2) |
| 3 | 0.0338 (1/3) | 0.0087 (1/3) | 1.0267 (2/3) |
| 4 | 0.0392 (1/4) | 0.0436 (1/4) | 1.1550 (1/4) |
| 5 | 0.0432 (4/5) | 0.0266 (4/5) | 1.2684 (4/5) |
| 6 | 0.0613 (5/6) | 0.0526 (1/6) | 1.7908 (5/6) |
| 7 | 0.0423 (5/7) | 0.0501 (3/7) | 1.1943 (2/7) |
| 8 | 0.0447 (3/8) | 0.0530 (5/8) | 1.3342 (3/8) |
| 9 | 0.0353 (5/9) | 0.0453 (7/9) | 1.2684 (5/9) |
| 10 | 0.0357 (1/10) | 0.0459 (3/10) | 1.0656 (1/10) |
| 11 | 0.0510 (10/11) | 0.0281 (5/11) | 1.5548 (10/11) |
| 12 | 0.0291 (5/12) | 0.0436 (11/12) | 0.9009 (5/12) |
| 13 | 0.0290 (7/13) | 0.0539 (3/13) | 0.8499 (1/13) |
| 14 | 0.0657 (13/14) | 0.0468 (9/14) | 1.9361 (13/14) |
| 15 | 0.0398 (14/15) | 0.0266 (13/15) | 1.1534 (14/15) |
| 16 | 0.0357 (13/16) | 0.0395 (15/16) | 1.0467 (13/16) |

## Collision bilinear model for the complete doublets

For each observed pair of prime labels `p<q`, the bilinear model independently replaces the first zeros by uniforms `r_p in [1,(p-3)/2]`, `r_q in [1,(q-3)/2]`. It exactly counts the four linear collision equations for `m_i(p)=m_j(q)` and doubles the unordered sum. Thus

```text
F_bil = 2 sum_{p<q} (C_A(p,q)+C_B(p,q)+C_C(p,q)+C_D(p,q))
                         / (((p-3)/2)((q-3)/2)).
```

The iid fixed-margin comparator is `F_iid=4K(K-1)/N`. The centered entries below subtract that same comparator; `F_bil` is a geometric independent-gap prediction, while `F` is the observed ordered collision count.

| N | observed F | F_bil | F_iid | F-F_iid | F_bil-F_iid | observed/bilinear |
|---:|---:|---:|---:|---:|---:|---:|
| 8192 | 6 | 7.9412 | 8.7021 | -2.7021 | -0.7609 | 0.7556 |
| 16384 | 18 | 17.2380 | 19.3462 | -1.3462 | -2.1082 | 1.0442 |
| 32768 | 32 | 28.8432 | 31.5640 | +0.4360 | -2.7207 | 1.1094 |
| 65536 | 48 | 42.5446 | 46.8904 | +1.1096 | -4.3458 | 1.1282 |

## True shell off-diagonal and Fourier reconstruction

The random-CRT prediction has mean zero. Its displayed standard deviation is the exact fixed-margin formula from `prop:random-crt-moment`, with the observed `A_p` and weights. Agreement of `E_direct` and `E_Fourier` is an identity check, not statistical evidence by itself.

| N | weight | active columns | hits | diagonal S_c | E_direct | E_Fourier | difference | E/S_c | random sd | E/sd | V/S_c |
|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8192 | 1 | 177 | 625 | 624.6696 | +2.6467 | +2.6467 | +0.00e+00 | +0.00424 | 9.7272 | +0.272 | 1.00424 |
| 8192 | log p | 177 | 625 | 46987.4351 | +233.3274 | +233.3274 | -9.09e-12 | +0.00497 | 731.7053 | +0.319 | 1.00497 |
| 16384 | 1 | 363 | 1241 | 1240.6791 | +2.3218 | +2.3218 | +0.00e+00 | +0.00187 | 13.6847 | +0.170 | 1.00187 |
| 16384 | log p | 363 | 1241 | 108781.1898 | +201.7107 | +201.7107 | +7.82e-11 | +0.00185 | 1199.8857 | +0.168 | 1.00185 |
| 32768 | 1 | 669 | 2381 | 2380.6747 | -16.6838 | -16.6838 | +0.00e+00 | -0.00701 | 18.5818 | -0.898 | 0.99299 |
| 32768 | log p | 669 | 2381 | 240379.4314 | -1655.3913 | -1655.3913 | -2.66e-10 | -0.00689 | 1876.2537 | -0.882 | 0.99311 |
| 65536 | 1 | 1145 | 4073 | 4072.7218 | +13.1451 | +13.1451 | +0.00e+00 | +0.00323 | 22.4866 | +0.585 | 1.00323 |
| 65536 | log p | 1145 | 4073 | 469724.4382 | +1597.1849 | +1597.1849 | +4.40e-10 | +0.00340 | 2593.4983 | +0.616 | 1.00340 |

## What the computation says (and does not say)

- Complete-doublet Parseval ratios `(T+F)/T` lie in `1.0224..1.0319`. Thus the integrated norm is cleanly at the diagonal scale on every tested block.
- With the explicitly declared `q<=16` major arcs removed, the largest sampled minor-arc values are `3.38..3.69` times `sqrt(T)`. This is compatible with square-root-scale cancellation plus an extreme-value factor; it is not a uniform minor-arc bound.
- The independent-uniform-gap bilinear model predicts the observed collision energy only to constant-factor/fluctuation accuracy: observed/predicted ranges from `0.756` to `1.128`. It is a useful null geometry, not an arithmetic theorem about `r_p`.
- In the true shell with logarithmic weights, `|E^o|/S_c` ranges from `0.00185` to `0.00689`, the ratios are nonmonotone, and the largest fixed-margin null-model magnitude is `0.88` standard deviations. The direct/Fourier computations agree to the displayed precision.
- Small-denominator correlations must be read denominator by denominator: exact or near degeneracies of the prime phase are major-arc structure, while the connected nondegenerate correlations are finite-sample diagnostics. None supplies the horizontal cross-prime estimate required by `eq:short-arc`.
- The two shell-count columns quantify a boundary/lifting issue: the complete reflected-pair sum contains only the designated `k=1` lifts `p+r`, while a true shell column can omit these below `N` and include higher lifts (`2p+r`, and sometimes `3p+r`). Thus the complete sum is not the exact Fourier transform of the shell columns. Any analytic reduction from this `S(theta)` to shell AMTD must insert the interval-dependent lift sums (or average/translate the interval); the complete-pair Parseval identity alone does not make that reduction.
- Consequently these computations support the doublet-cancellation heuristic but do not remove the generic `P^2` large-sieve barrier and do not prove `hyp:amtd`.
