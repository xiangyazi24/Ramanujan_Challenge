# Test F: dyadic gap profile for exact and fixed-degree detectors

## Scope and reproducibility

The full requested window was used: 208 primes in [500, 2000], hence 21528 pairs. No range reduction was needed. The C run took 4.25 seconds.

Reproduce with <code>python3 CRON_testF_driver.py</code> from this directory.

## Reconstruction of the missing definitions

Put $P=p-1$, $Q=q-1$, $D=q-p>0$, and $d=p-q=-D$. For an observable $F$, the centered cyclic row is

$$h_{p,F}(r)=F(T_p(r))-P^{-1}\sum_{s=0}^{P-1}F(T_p(s)),\qquad r\in\mathbb Z/P\mathbb Z.$$

The displayed equations missing from Q6420 can be reconstructed from Sections 2.3--2.5 as follows. With

$$\widehat h_{p,F}(u)=\sum_{r=0}^{P-1}h_{p,F}(r)e^{-2\pi iur/P},\qquad K_W(\theta)=\sum_r W(r)e^{2\pi ir\theta},$$

the shifted pair has the exact Fourier reconstruction

$$C_{p,q}(F;d,W)=\frac1{PQ}\sum_{u\bmod P}\sum_{v\bmod Q}\widehat h_{p,F}(u)\widehat h_{q,F}(v)e^{2\pi ivd/Q}K_W(u/P+v/Q).$$

For the exact indicator, additive orthogonality gives, for nonzero $u$, $\widehat h_p(u)=p^{-1}\sum_{a=1}^{p-1}\mathfrak N_{p,a}(u)$, where $\mathfrak N_{p,a}(u)=\sum_r e_p(aT_p(r))e^{-2\pi iur/P}$. Substitution for both primes is the fourfold $(a,b,u,v)$ formula. Centering sets the $u=0$ and $v=0$ coefficients to zero.

The cyclic centering also recovers Q6420 Section 2.1's endpoint conversion. If $z_p=|Z_p|$ and $g_p(r)=1_{Z_p}(r)-z_p/p$ on $\mathbb F_p$, then for $0\le r\le p-2$ one has $g_p(r)=h_p(r)+z_p/[p(p-1)]$, while the separate endpoint $r=p-1$ contributes $-z_p/p$. The fallback statistic below is defined directly with $h_p$, so it needs no hidden endpoint correction.

**Flagged fallback.** Q6420's export omits the displayed definition of $W$ and the ambient integer interval. Therefore the primary tables use the specification's explicit fallback

$$C^{\rm align}_{p,q}(F)=\sum_{r=0}^{p-2}h_{p,F}(r)h_{q,F}(r).$$

I also computed the linkage recovered from the prose: $d=p-q$ and $W=1_{[D,P)}$, so that

$$C^{\rm shift}_{p,q}(F)=\sum_{r=D}^{P-1}h_{p,F}(r)h_{q,F}(r-D).$$

This is the maximal hard window on which both original (non-wrapped) indices are admissible. It is empty when $D\ge P$. The unreported smooth/ambient $W$ in Q6420 is not guessed.

## Independence benchmark

Let $x_1,\dots,x_n$ and $y_1,\dots,y_m$ be centered rows, $Q_x=\sum x_i^2$, $Q_y=\sum y_j^2$, and let $\pi$ be a uniform random injection of the $n$ first-row positions into the $m$ second-row positions. Then

$$\mathbb E[y_{\pi(i)}^2]=Q_y/m,\qquad \mathbb E[y_{\pi(i)}y_{\pi(j)}]=-Q_y/[m(m-1)]\quad(i\ne j),$$

$$\operatorname{Var}\!\left(\sum_i x_i y_{\pi(i)}\right)=Q_xQ_y/(m-1).$$

This is the exact finite-population benchmark used for aligned and reflection-reduced rows. For a fixed subset $S$ of first-row positions in the shifted-overlap statistic, the exact variant used is

$$\operatorname{Var}=\frac{Q_y}{m-1}\left(\sum_{i\in S}x_i^2-\frac{(\sum_{i\in S}x_i)^2}{m}\right).$$

In every bin, <code>RMS ratio</code> means $\sqrt{\sum C_{p,q}^2/\sum V_{p,q}}$; <code>sum z</code> means $\sum C_{p,q}/\sqrt{\sum V_{p,q}}$. Powers use the canonical integer representative $T_p(r)\in[0,p)$ before centering.

## Primary aligned statistic

| D bin | pairs | F | sum C | RMS | benchmark RMS | RMS ratio | sum z |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [2,4) | 37 | indicator | -0.01176 | 0.00104 | 0.01780 | 0.059 | -0.11 |
| [2,4) | 37 | T | 2.6878e+07 | 7.5730e+06 | 7.1373e+06 | 1.061 | +0.62 |
| [2,4) | 37 | T^2 | -8.3917e+13 | 2.6387e+13 | 2.4050e+13 | 1.097 | -0.57 |
| [2,4) | 37 | T^3 | -4.4264e+20 | 8.2473e+19 | 7.5385e+19 | 1.094 | -0.97 |
| [4,8) | 119 | indicator | -0.07133 | 0.00215 | 0.02434 | 0.088 | -0.27 |
| [4,8) | 119 | T | 1.5606e+08 | 7.7275e+06 | 6.7577e+06 | 1.144 | +2.12 |
| [4,8) | 119 | T^2 | 3.3641e+14 | 2.5460e+13 | 2.1746e+13 | 1.171 | +1.42 |
| [4,8) | 119 | T^3 | 7.3962e+20 | 8.3074e+19 | 6.6119e+19 | 1.256 | +1.03 |
| [8,16) | 210 | indicator | -0.10170 | 0.00178 | 0.02183 | 0.082 | -0.32 |
| [8,16) | 210 | T | 3.0249e+08 | 6.2719e+06 | 6.3943e+06 | 0.981 | +3.26 |
| [8,16) | 210 | T^2 | 4.0360e+14 | 1.8543e+13 | 1.9913e+13 | 0.931 | +1.40 |
| [8,16) | 210 | T^3 | 8.0651e+20 | 5.3961e+19 | 5.8856e+19 | 0.917 | +0.95 |
| [16,32) | 465 | indicator | 0.60153 | 0.04613 | 0.02910 | 1.585 | +0.96 |
| [16,32) | 465 | T | 3.9992e+08 | 7.0088e+06 | 6.3794e+06 | 1.099 | +2.91 |
| [16,32) | 465 | T^2 | 1.2209e+14 | 2.2116e+13 | 2.0150e+13 | 1.098 | +0.28 |
| [16,32) | 465 | T^3 | -7.8007e+20 | 6.3832e+19 | 6.0069e+19 | 1.063 | -0.60 |
| [32,64) | 867 | indicator | 0.30827 | 0.03383 | 0.02789 | 1.213 | +0.38 |
| [32,64) | 867 | T | 1.1732e+09 | 6.6655e+06 | 6.3813e+06 | 1.045 | +6.24 |
| [32,64) | 867 | T^2 | 1.5448e+15 | 2.0360e+13 | 1.9867e+13 | 1.025 | +2.64 |
| [32,64) | 867 | T^3 | 2.2359e+21 | 5.8695e+19 | 5.8156e+19 | 1.009 | +1.31 |
| [64,128) | 1746 | indicator | 0.74462 | 0.03368 | 0.02666 | 1.263 | +0.67 |
| [64,128) | 1746 | T | 2.1501e+09 | 6.4285e+06 | 6.3477e+06 | 1.013 | +8.11 |
| [64,128) | 1746 | T^2 | 2.3902e+15 | 1.9504e+13 | 1.9680e+13 | 0.991 | +2.91 |
| [64,128) | 1746 | T^3 | 3.0526e+21 | 5.6985e+19 | 5.7247e+19 | 0.995 | +1.28 |
| [128,256) | 3197 | indicator | 3.65704 | 0.04296 | 0.02668 | 1.610 | +2.42 |
| [128,256) | 3197 | T | 4.3280e+09 | 6.3177e+06 | 5.9756e+06 | 1.057 | +12.81 |
| [128,256) | 3197 | T^2 | 5.9979e+15 | 1.8303e+13 | 1.7447e+13 | 1.049 | +6.08 |
| [128,256) | 3197 | T^3 | 1.0117e+22 | 5.0204e+19 | 4.7747e+19 | 1.051 | +3.75 |
| [256,512) | 5561 | indicator | -0.50993 | 0.02317 | 0.02431 | 0.953 | -0.28 |
| [256,512) | 5561 | T | 5.1629e+09 | 5.4475e+06 | 5.4242e+06 | 1.004 | +12.76 |
| [256,512) | 5561 | T^2 | 4.3538e+15 | 1.4474e+13 | 1.4561e+13 | 0.994 | +4.01 |
| [256,512) | 5561 | T^3 | 5.9207e+21 | 3.6618e+19 | 3.6593e+19 | 1.001 | +2.17 |
| [512,1024) | 7192 | indicator | 4.44406 | 0.03325 | 0.02232 | 1.490 | +2.35 |
| [512,1024) | 7192 | T | 7.3878e+09 | 4.3180e+06 | 4.1424e+06 | 1.042 | +21.03 |
| [512,1024) | 7192 | T^2 | 6.6745e+15 | 9.0044e+12 | 8.9135e+12 | 1.010 | +8.83 |
| [512,1024) | 7192 | T^3 | 8.3346e+21 | 1.8094e+19 | 1.8139e+19 | 0.998 | +5.42 |
| [1024,2048) | 2134 | indicator | -0.02340 | 0.02151 | 0.02202 | 0.977 | -0.02 |
| [1024,2048) | 2134 | T | 1.6661e+09 | 2.8268e+06 | 2.7758e+06 | 1.018 | +12.99 |
| [1024,2048) | 2134 | T^2 | 1.0245e+15 | 4.0739e+12 | 4.1205e+12 | 0.989 | +5.38 |
| [1024,2048) | 2134 | T^3 | 7.7907e+20 | 5.4970e+18 | 5.5758e+18 | 0.986 | +3.02 |

## Midpoint and reflection removal (indicator)

<code>no-mid</code> replaces the midpoint indicator by zero and re-centers on the full cycle. <code>reflection-reduced</code> then projects to the reflection-even part, discards the midpoint orbit, keeps representatives $0\le r<(p-1)/2$, and re-centers the representative row. The program also checked $T_p(r)=T_p(-r\bmod p-1)$ exactly for every row.

| D bin | full R | no-mid R | reduced R | midpoint share of excess | total forced share | sum z: full → reduced |
| --- | --- | --- | --- | --- | --- | --- |
| [2,4) | 0.059 | 0.059 | 0.041 | n/a | n/a | -0.11 → -0.08 |
| [4,8) | 0.088 | 0.088 | 0.062 | n/a | n/a | -0.27 → -0.19 |
| [8,16) | 0.082 | 0.082 | 0.057 | n/a | n/a | -0.32 → -0.22 |
| [16,32) | 1.585 | 1.585 | 0.070 | +0.0% | +165.7% | +0.96 → -0.45 |
| [32,64) | 1.213 | 1.213 | 1.716 | +0.0% | -313.3% | +0.38 → +1.16 |
| [64,128) | 1.263 | 1.263 | 1.268 | +0.0% | -1.8% | +0.67 → +0.55 |
| [128,256) | 1.610 | 1.610 | 1.613 | +0.0% | -0.6% | +2.42 → +1.81 |
| [256,512) | 0.953 | 0.953 | 0.779 | n/a | n/a | -0.28 → -0.37 |
| [512,1024) | 1.490 | 1.490 | 1.489 | +0.0% | +0.4% | +2.35 → +1.68 |
| [1024,2048) | 0.977 | 0.977 | 1.379 | n/a | n/a | -0.02 → +0.84 |

The percentages compare normalized excess energy $R^2-1$; negative percentages mean that removal increased, rather than explained, the observed excess. <code>n/a</code> means the original bin had no positive excess.

There were 0 midpoint hits in this prime window. Consequently midpoint deletion changes no bin (0% of every positive excess); all changes in the reduced column come from passing to one reflection-orbit representative and re-centering.

## Reconstructed shifted linkage

| D bin | nonempty/all | F | sum C | RMS | benchmark RMS | RMS ratio | sum z |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [2,4) | 37/37 | indicator | -0.01178 | 0.00104 | 0.01780 | 0.059 | -0.11 |
| [2,4) | 37/37 | T | -4.4757e+06 | 6.3301e+06 | 7.1247e+06 | 0.888 | -0.10 |
| [2,4) | 37/37 | T^2 | 3.4663e+13 | 2.3030e+13 | 2.4033e+13 | 0.958 | +0.24 |
| [2,4) | 37/37 | T^3 | 2.6926e+20 | 7.5199e+19 | 7.5354e+19 | 0.998 | +0.59 |
| [4,8) | 119/119 | indicator | -0.07181 | 0.00217 | 0.02434 | 0.089 | -0.27 |
| [4,8) | 119/119 | T | 1.4273e+08 | 6.6857e+06 | 6.7350e+06 | 0.993 | +1.94 |
| [4,8) | 119/119 | T^2 | 5.1541e+14 | 2.1403e+13 | 2.1708e+13 | 0.986 | +2.18 |
| [4,8) | 119/119 | T^3 | 1.4310e+21 | 6.4131e+19 | 6.6040e+19 | 0.971 | +1.99 |
| [8,16) | 210/210 | indicator | -0.10284 | 0.00181 | 0.02183 | 0.083 | -0.33 |
| [8,16) | 210/210 | T | 6.0890e+07 | 6.1738e+06 | 6.3610e+06 | 0.971 | +0.66 |
| [8,16) | 210/210 | T^2 | 1.1367e+14 | 1.7976e+13 | 1.9843e+13 | 0.906 | +0.40 |
| [8,16) | 210/210 | T^3 | 1.5412e+19 | 4.9983e+19 | 5.8682e+19 | 0.852 | +0.02 |
| [16,32) | 465/465 | indicator | 0.60731 | 0.04605 | 0.02868 | 1.606 | +0.98 |
| [16,32) | 465/465 | T | -5.0611e+07 | 6.4085e+06 | 6.3220e+06 | 1.014 | -0.37 |
| [16,32) | 465/465 | T^2 | 1.3068e+14 | 2.1926e+13 | 2.0011e+13 | 1.096 | +0.30 |
| [16,32) | 465/465 | T^3 | 7.0628e+20 | 6.9120e+19 | 5.9702e+19 | 1.158 | +0.55 |
| [32,64) | 867/867 | indicator | 0.34276 | 0.03346 | 0.02737 | 1.222 | +0.43 |
| [32,64) | 867/867 | T | -7.9478e+06 | 6.3377e+06 | 6.2737e+06 | 1.010 | -0.04 |
| [32,64) | 867/867 | T^2 | -3.1770e+14 | 1.9513e+13 | 1.9585e+13 | 0.996 | -0.55 |
| [32,64) | 867/867 | T^3 | -1.6350e+21 | 5.6712e+19 | 5.7390e+19 | 0.988 | -0.97 |
| [64,128) | 1746/1746 | indicator | -1.14801 | 0.00230 | 0.02517 | 0.091 | -1.09 |
| [64,128) | 1746/1746 | T | -5.5514e+07 | 5.9804e+06 | 6.1394e+06 | 0.974 | -0.22 |
| [64,128) | 1746/1746 | T^2 | -2.6794e+13 | 1.8360e+13 | 1.9112e+13 | 0.961 | -0.03 |
| [64,128) | 1746/1746 | T^3 | -1.5175e+20 | 5.2951e+19 | 5.5679e+19 | 0.951 | -0.07 |
| [128,256) | 3197/3197 | indicator | 2.09625 | 0.03534 | 0.02374 | 1.489 | +1.56 |
| [128,256) | 3197/3197 | T | -2.6100e+08 | 5.5001e+06 | 5.5689e+06 | 0.988 | -0.83 |
| [128,256) | 3197/3197 | T^2 | -1.1180e+15 | 1.6362e+13 | 1.6383e+13 | 0.999 | -1.21 |
| [128,256) | 3197/3197 | T^3 | -4.1135e+21 | 4.5418e+19 | 4.4967e+19 | 1.010 | -1.62 |
| [256,512) | 5558/5561 | indicator | -0.71484 | 0.01341 | 0.01824 | 0.735 | -0.53 |
| [256,512) | 5558/5561 | T | 4.5702e+07 | 4.6344e+06 | 4.5905e+06 | 1.010 | +0.13 |
| [256,512) | 5558/5561 | T^2 | 4.5530e+14 | 1.2641e+13 | 1.2580e+13 | 1.005 | +0.49 |
| [256,512) | 5558/5561 | T^3 | 1.3975e+21 | 3.1979e+19 | 3.1920e+19 | 1.002 | +0.59 |
| [512,1024) | 4466/7192 | indicator | -0.48991 | 5.2046e-04 | 0.00857 | 0.061 | -0.67 |
| [512,1024) | 4466/7192 | T | 5.1948e+07 | 2.2834e+06 | 2.3719e+06 | 0.963 | +0.26 |
| [512,1024) | 4466/7192 | T^2 | 2.4662e+14 | 5.3936e+12 | 5.6533e+12 | 0.954 | +0.51 |
| [512,1024) | 4466/7192 | T^3 | 4.5313e+20 | 1.1528e+19 | 1.2217e+19 | 0.944 | +0.44 |
| [1024,2048) | 0/2134 | indicator | 0 | 0 | 0 | 0.000 | +0.00 |
| [1024,2048) | 0/2134 | T | 0 | 0 | 0 | 0.000 | +0.00 |
| [1024,2048) | 0/2134 | T^2 | 0 | 0 | 0 | 0.000 | +0.00 |
| [1024,2048) | 0/2134 | T^3 | 0 | 0 | 0 | 0.000 | +0.00 |

## Fourier self-check and near-resonant bookkeeping

| window | D bin | sample pair | direct | Fourier | absolute error |
| --- | --- | --- | --- | --- | --- |
| aligned | [2,4) | (1487,1489) | -0.00269 | -0.00269 | 3.4968e-15 |
| shifted_overlap | [2,4) | (1487,1489) | -0.00270 | -0.00270 | 5.7103e-15 |
| aligned | [4,8) | (571,577) | -0.00702 | -0.00702 | 6.5919e-17 |
| shifted_overlap | [4,8) | (571,577) | -0.00709 | -0.00709 | 5.3004e-15 |
| aligned | [8,16) | (617,631) | -0.00649 | -0.00649 | 2.1225e-14 |
| shifted_overlap | [8,16) | (617,631) | -0.00664 | -0.00664 | 7.6241e-16 |
| aligned | [16,32) | (541,571) | -0.00741 | -0.00741 | 7.7091e-15 |
| shifted_overlap | [16,32) | (541,571) | -0.00429 | -0.00429 | 4.5216e-15 |
| aligned | [32,64) | (503,541) | -0.00398 | -0.00398 | 3.7453e-15 |
| shifted_overlap | [32,64) | (503,541) | -0.00454 | -0.00454 | 1.4956e-14 |
| aligned | [64,128) | (503,571) | -0.00797 | -0.00797 | 1.5491e-15 |
| shifted_overlap | [64,128) | (503,571) | -0.00541 | -0.00541 | 1.4466e-14 |
| aligned | [128,256) | (503,631) | -0.00398 | -0.00398 | 1.9863e-16 |
| shifted_overlap | [128,256) | (503,631) | -0.00243 | -0.00243 | 4.2110e-16 |
| aligned | [256,512) | (503,761) | 0.98406 | 0.98406 | 3.1197e-14 |
| shifted_overlap | [256,512) | (503,761) | -1.4678e-04 | -1.4678e-04 | 3.2440e-15 |
| aligned | [512,1024) | (503,1049) | -0.00398 | -0.00398 | 6.3404e-15 |
| shifted_overlap | [512,1024) | (541,1069) | -0.00537 | -0.00537 | 2.1927e-15 |
| aligned | [1024,2048) | (503,1531) | -0.00398 | -0.00398 | 1.9125e-14 |

The maximum direct-versus-Fourier absolute error was 3.120e-14. Thus the reconstruction succeeded numerically.

For each pair, the near block is $v=-u$, $1\le |u|\le\lfloor((p+q)/2)/D\rfloor$, truncated to the signed frequency ranges. Both signs are included. The last two columns profile the exact indicator over every cyclic shift $|\delta|\le D$.

| D bin | near nonempty/all | sum near | RMS near | near/full shifted RMS | cyclic max RMS ratio | endpoint ±D attains max |
| --- | --- | --- | --- | --- | --- | --- |
| [2,4) | 37/37 | -0.01170 | 0.00104 | 0.995 | 0.059 | 4/4 |
| [4,8) | 119/119 | -0.06586 | 0.00193 | 0.893 | 0.088 | 12/12 |
| [8,16) | 210/210 | 0.07288 | 0.00676 | 3.747 | 5.452 | 17/20 |
| [16,32) | 465/465 | 0.11825 | 0.00650 | 0.141 | 6.316 | 57/71 |
| [32,64) | 867/867 | 0.11582 | 0.00453 | 0.135 | 8.207 | 80/121 |
| [64,128) | 1746/1746 | 0.03416 | 0.00264 | 1.146 | 9.767 | 115/230 |
| [128,256) | 3197/3197 | -0.08255 | 0.00161 | 0.046 | 12.555 | 90/436 |
| [256,512) | 5558/5561 | -0.05090 | 9.9104e-04 | 0.074 | 14.783 | 29/741 |
| [512,1024) | 4466/7192 | 0.04900 | 2.4128e-04 | 0.464 | 16.617 | 10/979 |
| [1024,2048) | 0/2134 | 0 | 0 | n/a | 17.539 | 1/312 |

## Ground-truth and structural checks

| check | scope | result |
| --- | --- | --- |
| $p=17$ | $Z_{17}=\{3,13\}$ | PASS |
| Parity-law sample | first 20 window primes | PASS (0 failures) |
| Parity law, stronger run | all 208 window primes | PASS (0 failures) |
| Reflection | all 208 complete rows | PASS (0 failures) |
| Fourier reconstruction | 19 selected bin/window checks | PASS (max error 3.120e-14) |

## Verdict

| D bin | pairs | indicator R / z | T R | T^2 R | T^3 R | predicted split? |
| --- | --- | --- | --- | --- | --- | --- |
| [2,4) | 37 | 0.059 / -0.11 (DEFICIT) | 1.061 (~1) | 1.097 (~1) | 1.094 (~1) | NO |
| [4,8) | 119 | 0.088 / -0.27 (DEFICIT) | 1.144 (~1) | 1.171 (~1) | 1.256 (EXCESS) | NO |
| [8,16) | 210 | 0.082 / -0.32 (DEFICIT) | 0.981 (~1) | 0.931 (~1) | 0.917 (~1) | NO |
| [16,32) | 465 | 1.585 / +0.96 (EXCESS) | 1.099 (~1) | 1.098 (~1) | 1.063 (~1) | YES |
| [32,64) | 867 | 1.213 / +0.38 (~1) | 1.045 (~1) | 1.025 (~1) | 1.009 (~1) | NO |
| [64,128) | 1746 | 1.263 / +0.67 (EXCESS) | 1.013 (~1) | 0.991 (~1) | 0.995 (~1) | YES |
| [128,256) | 3197 | 1.610 / +2.42 (EXCESS) | 1.057 (~1) | 1.049 (~1) | 1.051 (~1) | YES |
| [256,512) | 5561 | 0.953 / -0.28 (~1) | 1.004 (~1) | 0.994 (~1) | 1.001 (~1) | NO |
| [512,1024) | 7192 | 1.490 / +2.35 (EXCESS) | 1.042 (~1) | 1.010 (~1) | 0.998 (~1) | YES |
| [1024,2048) | 2134 | 0.977 / -0.02 (~1) | 1.018 (~1) | 0.989 (~1) | 0.986 (~1) | NO |

**Signature verdict: MIXED.**

All three fixed-degree surrogates are simultaneously in the ~1 band in 9/10 bins. The exact indicator has excess RMS in 4/10 bins, a deficit in 3/10, and is in-band in 3/10. Thus the predicted surrogate/indicator split is clear in 4/10 adequately populated bins, but it is not uniform across the dyadic profile.

The table uses a declared descriptive rule: <code>~1</code> means $0.80\le R\le1.25$; <code>EXCESS</code> means $R>1.25$; and <code>DEFICIT</code> means $R<0.80$. The separately reported sum z records signed drift, but pair statistics sharing a prime are not independent, so it is not used as a formal z-test. Bins with fewer than 30 pairs are not used for the overall verdict. These thresholds are reporting conventions, not formal significance tests.
