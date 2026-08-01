# Joint Tannakian moments v2: exact integral-trace fingerprint

## Verdict

All corrected Q6457-recipe trace gates passed before any FFT was computed.  The exact prime set was `29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199`.

The exact trace stage took 17.852 seconds and the complete run took 18.033 seconds on this machine; after covering [29,149], the computation therefore extended through 199, well below the 25-minute cutoff.

## Exact gates

`g1` is the integral sum identity, `g2`/`g3` are the two Apéry residues, `g4` is the pointwise `3p` bound, `g5` is full Mellin inversion (with the endpoint alias checked separately), and `g6` is the p=29 raw-count checksum.

| p | eps=(-3/p) | split | inert | branch | g1 | g2 | g3 | g4 | g5 | g6 | max |T+|/p | max |T-|/p | trace SHA |
|---:|---:|---:|---:|---:|:---:|:---:|:---:|:---:|:---:|:---:|---:|---:|:---|
| 29 | -1 | 13 | 15 | 0 | PASS | PASS | PASS | PASS | PASS | PASS | 2.310 | 2.310 | `c08527d6379075e1` |
| 31 | +1 | 13 | 15 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.097 | 2.097 | `87e1c113335d5a94` |
| 37 | +1 | 17 | 19 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 1.703 | 1.703 | `2444db48553bfd14` |
| 41 | -1 | 18 | 20 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.659 | 2.659 | `f5bdbd78b68d3dfe` |
| 43 | +1 | 20 | 22 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.488 | 2.488 | `a8605d74acf4f314` |
| 47 | -1 | 21 | 23 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.191 | 2.191 | `249c521e081125eb` |
| 53 | -1 | 25 | 27 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 1.830 | 1.830 | `776b4289b7903f2d` |
| 59 | -1 | 28 | 30 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.661 | 2.661 | `cf417de9b99d55ef` |
| 61 | +1 | 29 | 31 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.541 | 2.541 | `6f1f667fd82911a5` |
| 67 | +1 | 32 | 34 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.821 | 2.821 | `f0590c40376cbf31` |
| 71 | -1 | 33 | 35 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.042 | 2.042 | `b61976784a7bd9a6` |
| 73 | +1 | 34 | 36 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.507 | 2.507 | `3029231a3fa9e6a1` |
| 79 | +1 | 37 | 39 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.722 | 2.722 | `53f6cdd30c33b198` |
| 83 | -1 | 40 | 42 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.904 | 2.904 | `d678bdb9eff5ba0a` |
| 89 | -1 | 42 | 44 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.640 | 2.640 | `704029fc0bce74e6` |
| 97 | +1 | 46 | 48 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.031 | 2.031 | `c77268da561afee9` |
| 101 | -1 | 49 | 51 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.802 | 2.802 | `46270c3944571ab6` |
| 103 | +1 | 49 | 51 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.883 | 2.883 | `d3670ef92fcd8bb4` |
| 107 | -1 | 52 | 54 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.589 | 2.589 | `0c6306ef0ef07c74` |
| 109 | +1 | 53 | 55 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.670 | 2.670 | `bb7a4915bbd6f89c` |
| 113 | -1 | 54 | 56 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.398 | 2.398 | `f5726e29efeaba2e` |
| 127 | +1 | 61 | 63 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.827 | 2.827 | `188f4b6f4bb46d80` |
| 131 | -1 | 64 | 66 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.710 | 2.710 | `d3f8ebdb7f939701` |
| 137 | -1 | 66 | 68 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.547 | 2.547 | `0f6ba53b44f5cdb5` |
| 139 | +1 | 68 | 70 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.496 | 2.496 | `2275920b4065b8a9` |
| 149 | -1 | 73 | 75 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.866 | 2.866 | `3df79fa12be12f4f` |
| 151 | +1 | 73 | 75 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.219 | 2.219 | `2dcf205bce066dc4` |
| 157 | +1 | 77 | 79 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.822 | 2.822 | `56fe02f3f59f2c85` |
| 163 | +1 | 80 | 82 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.681 | 2.681 | `3f62fff6c51e1c36` |
| 167 | -1 | 81 | 83 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.593 | 2.593 | `66f7903325565140` |
| 173 | -1 | 85 | 87 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.468 | 2.468 | `1974d7ef92cbfebe` |
| 179 | -1 | 88 | 90 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.352 | 2.352 | `25e2f4d4dc8b43a5` |
| 181 | +1 | 89 | 91 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.735 | 2.735 | `3e8f7a35319b8ca7` |
| 191 | -1 | 93 | 95 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.801 | 2.801 | `03d1d37558dee996` |
| 193 | +1 | 94 | 96 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.762 | 2.762 | `cf71ec081c7c0d93` |
| 197 | -1 | 97 | 99 | 0 | PASS | PASS | PASS | PASS | PASS | -- | 2.685 | 2.685 | `c862664f8316eb71` |
| 199 | +1 | 97 | 99 | 2 | PASS | PASS | PASS | PASS | PASS | -- | 2.940 | 2.940 | `5e945da2c82cdfeb` |

Q6457 checksums verified explicitly:

- The recurrence produced all 29 Apéry coefficients modulo 29 identically to the direct binomial formula and to `CRON_pushforward_check.py`.
- At p=29, all 13 split fibres (26 source parameters) agreed between exact point counts, centered Franel Hasse--Witt residues, and the raw CRON count convention.
- All 15 inert p=29 counts were unchanged after replacing the chosen square root by its conjugate.
- The named checksum is `t=2`, sources (8, 10), elliptic traces (-6, 6), and `f=7`.
- For every prime and every nonzero t, the integral sum, both residual congruences, and both `3p` bounds passed; the branch count was two exactly when `(2/p)=+1`.
- The inert normalization used the corrected certified sign `eps=(-3/p)` at every prime; `(2/p)` was used only for the independent branch-count checksum.
- For every prime, Mellin inversion passed for every `1 <= r <= p-2`; at `r=0` the raw residue was exactly `b_0+b_{p-1}`.

The F_p traces used direct O(p) point counts.  The p=29 comparison against the centered Hasse--Witt polynomial is the required two-method cross-check.  Inert traces used exact O(p^2) point counts in `F_p[z]/(z^2-d)`, with the extension quadratic character evaluated exactly through the norm to F_p.

## Absolute Mellin moments

For a primitive root g, the script stores `T(g^j)` and computes `S(chi_k)=sum_j exp(2*pi*i*k*j/(p-1)) T(g^j)`.  The normalization is `s=S/p^(3/2)`.  The table excludes only k=0.

| p | plus mu2 | plus mu4 | plus mu6 | minus mu2 | minus mu4 | minus mu6 |
|---:|---:|---:|---:|---:|---:|---:|
| 29 | 0.764516 | 1.154706 | 2.342900 | 0.762967 | 1.452346 | 3.539476 |
| 31 | 0.874085 | 1.775177 | 4.781647 | 0.831985 | 2.031254 | 7.936213 |
| 37 | 0.857327 | 1.826119 | 5.437566 | 0.820363 | 1.568280 | 4.132319 |
| 41 | 1.215827 | 3.271419 | 11.016267 | 1.070628 | 2.495028 | 7.645865 |
| 43 | 0.903985 | 1.826532 | 4.361443 | 0.903220 | 2.616656 | 10.882980 |
| 47 | 1.350563 | 3.547862 | 12.295696 | 1.210772 | 3.532593 | 12.288257 |
| 53 | 0.696618 | 1.049346 | 2.156504 | 0.691350 | 1.543777 | 5.277537 |
| 59 | 1.282226 | 2.683625 | 6.495851 | 1.245242 | 4.027622 | 19.785341 |
| 61 | 0.968629 | 1.918524 | 5.133887 | 0.947176 | 2.322719 | 7.241769 |
| 67 | 0.877693 | 1.492874 | 3.148822 | 0.837539 | 1.614491 | 4.575545 |
| 71 | 0.934868 | 1.938154 | 5.864202 | 0.841916 | 1.514308 | 3.770685 |
| 73 | 1.037704 | 1.590049 | 2.902083 | 1.050701 | 3.335406 | 13.730482 |
| 79 | 1.013556 | 2.085963 | 5.226653 | 1.018310 | 2.259642 | 6.583052 |
| 83 | 0.926170 | 1.532374 | 3.600716 | 0.920986 | 3.001550 | 16.272720 |
| 89 | 1.097532 | 2.856030 | 9.726282 | 1.028041 | 2.636115 | 10.414989 |
| 97 | 0.800373 | 1.344732 | 2.952822 | 0.725855 | 1.374624 | 3.751318 |
| 101 | 1.201324 | 2.522269 | 6.444756 | 1.197604 | 3.337873 | 12.668840 |
| 103 | 0.918994 | 1.748423 | 4.257434 | 0.938265 | 2.401354 | 8.150916 |
| 107 | 1.055210 | 2.190991 | 5.854113 | 1.038319 | 2.630109 | 10.439263 |
| 109 | 0.946881 | 1.794486 | 4.310784 | 0.936627 | 2.172290 | 6.842686 |
| 113 | 1.081386 | 2.314200 | 6.177380 | 1.023904 | 2.511161 | 9.149677 |
| 127 | 1.226745 | 2.774106 | 7.674216 | 1.165477 | 4.229175 | 20.124210 |
| 131 | 0.749483 | 1.312310 | 2.800281 | 0.731337 | 1.740672 | 6.095030 |
| 137 | 0.964536 | 2.095235 | 6.091818 | 0.916626 | 2.086849 | 6.200338 |
| 139 | 1.111473 | 2.210356 | 5.337126 | 1.111076 | 4.446822 | 27.154684 |
| 149 | 0.959885 | 1.672278 | 3.664609 | 0.954673 | 2.412208 | 8.441279 |
| 151 | 0.763340 | 1.233577 | 2.508314 | 0.713193 | 1.295896 | 3.258107 |
| 157 | 1.010712 | 2.008012 | 5.065482 | 1.010080 | 3.154188 | 15.732265 |
| 163 | 0.988159 | 1.981522 | 4.881655 | 0.983250 | 2.630354 | 9.079752 |
| 167 | 1.159574 | 2.758343 | 8.042099 | 1.113796 | 3.579343 | 17.899563 |
| 173 | 0.972684 | 1.923931 | 4.675291 | 0.969186 | 3.004585 | 14.717594 |
| 179 | 0.954015 | 1.758272 | 4.025045 | 0.953730 | 2.652906 | 9.690478 |
| 181 | 1.025133 | 2.164499 | 5.750139 | 1.021425 | 3.240112 | 15.302011 |
| 191 | 1.198724 | 2.669414 | 7.362944 | 1.158552 | 4.271710 | 23.341583 |
| 193 | 1.055572 | 2.307702 | 6.272125 | 1.038618 | 2.868585 | 12.012375 |
| 197 | 0.806357 | 1.608065 | 3.965822 | 0.793512 | 1.731849 | 5.200305 |
| 199 | 0.918371 | 1.722189 | 4.161145 | 0.916420 | 2.308810 | 8.683214 |

### Trend and predictions

The largest-five-prime average uses p = 181, 191, 193, 197, 199.

| object/model | mu2 | mu4 | mu6 |
|:---|---:|---:|---:|
| plus, first five observed | 0.923148 | 1.970791 | 5.587965 |
| plus, largest five observed | 1.000832 | 2.094374 | 5.502435 |
| plus: SL2 standard | 1 | 2 | 5 |
| plus: O2 normalizer | 1 | 3 | 10 |
| plus: finite subgroup | group-dependent rational/discrete | group-dependent rational/discrete | group-dependent rational/discrete |
| minus, first five observed | 0.877832 | 2.032713 | 6.827371 |
| minus, largest five observed | 0.985705 | 2.884213 | 12.907898 |
| minus: Sp4 standard | 1 | 3 | 14 |
| minus: Sym^3(SL2) | 1 | 4 | 34 |

The plus fourth moment is numerically closer to **SL2**.
The G_- dichotomy by the requested fourth-moment test favors **Sp4**.

## Joint moments

The finite-p product prediction is the product of the two observed second moments.  For a Sym^3 graph, SU2 representation theory gives `E(|x|^2 |Sym^3 x|^2)=2`; the calibrated graph column is therefore twice the observed product prediction.

| p | C22 | product | Sym3 graph | avg s+ conj(s-) |
|---:|---:|---:|---:|:---|
| 29 | 1.109340 | 0.583300 | 1.166601 | -0.508389+0.000000i |
| 31 | 0.427891 | 0.727225 | 1.454450 | -0.178633+0.000000i |
| 37 | 1.081147 | 0.703319 | 1.406638 | 0.163258+0.000000i |
| 41 | 0.990392 | 1.301698 | 2.603395 | -0.327331+0.000000i |
| 43 | 1.222302 | 0.816498 | 1.632995 | -0.344572+0.000000i |
| 47 | 1.621386 | 1.635223 | 3.270446 | 0.257633+0.000000i |
| 53 | 0.490968 | 0.481607 | 0.963213 | -0.056280+0.000000i |
| 59 | 1.058354 | 1.596682 | 3.193363 | -0.075879+0.000000i |
| 61 | 0.724431 | 0.917462 | 1.834923 | -0.072654+0.000000i |
| 67 | 0.609350 | 0.735102 | 1.470205 | 0.103749+0.000000i |
| 71 | 0.828845 | 0.787080 | 1.574160 | -0.097352+0.000000i |
| 73 | 0.902795 | 1.090318 | 2.180635 | 0.126575+0.000000i |
| 79 | 1.218964 | 1.032115 | 2.064229 | -0.205032+0.000000i |
| 83 | 0.712129 | 0.852990 | 1.705980 | -0.159199+0.000000i |
| 89 | 0.980534 | 1.128307 | 2.256615 | 0.110863+0.000000i |
| 97 | 0.536772 | 0.580955 | 1.161910 | -0.218899+0.000000i |
| 101 | 1.291684 | 1.438710 | 2.877420 | 0.141610+0.000000i |
| 103 | 0.888618 | 0.862260 | 1.724520 | 0.006240+0.000000i |
| 107 | 0.916162 | 1.095644 | 2.191287 | -0.281752+0.000000i |
| 109 | 1.138459 | 0.886874 | 1.773748 | 0.036002+0.000000i |
| 113 | 1.092156 | 1.107236 | 2.214473 | 0.078189+0.000000i |
| 127 | 1.642179 | 1.429743 | 2.859486 | -0.045669+0.000000i |
| 131 | 0.667026 | 0.548125 | 1.096250 | -0.199188+0.000000i |
| 137 | 0.689360 | 0.884119 | 1.768238 | -0.109513+0.000000i |
| 139 | 1.251971 | 1.234931 | 2.469862 | 0.058901+0.000000i |
| 149 | 0.836657 | 0.916376 | 1.832752 | 0.009460+0.000000i |
| 151 | 0.509500 | 0.544408 | 1.088817 | -0.043410+0.000000i |
| 157 | 1.034278 | 1.020900 | 2.041800 | -0.012488+0.000000i |
| 163 | 1.142801 | 0.971608 | 1.943215 | -0.054599+0.000000i |
| 167 | 1.187697 | 1.291529 | 2.583058 | 0.159286+0.000000i |
| 173 | 1.055769 | 0.942712 | 1.885424 | -0.005767+0.000000i |
| 179 | 0.721946 | 0.909872 | 1.819744 | 0.047381+0.000000i |
| 181 | 0.736801 | 1.047097 | 2.094193 | 0.010861+0.000000i |
| 191 | 1.641064 | 1.388784 | 2.777569 | 0.102200+0.000000i |
| 193 | 1.172748 | 1.096336 | 2.192672 | -0.063519+0.000000i |
| 197 | 0.774831 | 0.639854 | 1.279707 | -0.117060+0.000000i |
| 199 | 0.862810 | 0.841613 | 1.683227 | -0.195219+0.000000i |

Largest-five averages: `C22=1.037651`, `product=1.002737`, `Sym3 graph=2.005474`.  The unshifted C22 fingerprint is closer to **product**.
The corresponding naive covariance is `-0.052547+0.000000i` (absolute value 0.052547); as predicted, this statistic is not the product/graph discriminator.

## Twisted C22 correlations

Each entry is `specified/generic`: the first average follows the literal k=1,...,p-2 range in the spec, while the second removes the single k for which chi*eta becomes the exceptional trivial character.

| p | eta order 2 | eta order 3 | eta order 4 | eta order 6 |
|---:|:---:|:---:|:---:|:---:|
| 29 | 0.5861/0.6061 | -- | 0.9005/0.9351 | -- |
| 31 | 0.6397/0.6549 | 0.5780/0.5956 | -- | 0.7839/0.8061 |
| 37 | 0.9918/1.0013 | 1.0145/0.9054 | 0.8996/0.8539 | 0.7476/0.7468 |
| 41 | 1.4969/1.5362 | -- | 1.1659/1.1703 | -- |
| 43 | 0.7578/0.7767 | 0.6555/0.6718 | -- | 0.9867/1.0108 |
| 47 | 1.9571/1.9156 | -- | -- | -- |
| 53 | 0.6845/0.6961 | -- | 0.4192/0.4215 | -- |
| 59 | 1.5538/1.4998 | -- | -- | -- |
| 61 | 1.0644/1.0536 | 0.7373/0.7364 | 0.8119/0.8258 | 0.8938/0.8988 |
| 67 | 0.5728/0.5769 | 0.5791/0.5797 | -- | 0.7276/0.6879 |
| 71 | 0.9072/0.8900 | -- | -- | -- |
| 73 | 1.0258/1.0403 | 1.0240/1.0359 | 1.1903/1.2054 | 1.1213/1.1364 |
| 79 | 1.0365/1.0439 | 0.8765/0.8846 | -- | 1.1246/1.1394 |
| 83 | 0.7768/0.7860 | -- | -- | -- |
| 89 | 1.0658/1.0726 | -- | 1.0856/1.0982 | -- |
| 97 | 0.6579/0.6610 | 0.6136/0.5993 | 0.6615/0.6683 | 0.5520/0.5248 |
| 101 | 1.4971/1.5035 | -- | 1.4260/1.4361 | -- |
| 103 | 0.9463/0.9556 | 1.0637/1.0742 | -- | 0.8063/0.8141 |
| 107 | 0.8223/0.8217 | -- | -- | -- |
| 109 | 1.1976/1.1771 | 0.7109/0.7116 | 0.9478/0.9557 | 1.0643/1.0714 |
| 113 | 1.4428/1.4550 | -- | 0.8162/0.8100 | -- |
| 127 | 1.2628/1.2485 | 1.3828/1.3511 | -- | 1.7039/1.7145 |
| 131 | 0.6360/0.5820 | -- | -- | -- |
| 137 | 0.7670/0.7641 | -- | 0.7661/0.7696 | -- |
| 139 | 1.6103/1.6220 | 1.3149/1.3240 | -- | 1.1373/1.1456 |
| 149 | 0.8908/0.8968 | -- | 1.2260/1.2280 | -- |
| 151 | 0.4761/0.4732 | 0.4399/0.4410 | -- | 0.6305/0.6102 |
| 157 | 0.9176/0.9219 | 1.1901/1.1961 | 0.9266/0.9321 | 1.3758/1.3841 |
| 163 | 0.8936/0.8887 | 1.1143/1.1139 | -- | 1.0942/1.0934 |
| 167 | 1.3996/1.4077 | -- | -- | -- |
| 173 | 0.8162/0.8156 | -- | 0.9019/0.9068 | -- |
| 179 | 1.2171/1.2235 | -- | -- | -- |
| 181 | 1.0818/1.0842 | 1.3298/1.3364 | 1.0057/1.0113 | 1.0528/1.0574 |
| 191 | 1.9329/1.9431 | -- | -- | -- |
| 193 | 1.0367/1.0376 | 0.9784/0.9830 | 1.2286/1.2212 | 1.3071/1.3072 |
| 197 | 0.7565/0.7584 | -- | 0.6789/0.6763 | -- |
| 199 | 0.6988/0.7023 | 0.7175/0.7211 | -- | 0.7940/0.7976 |

Largest-five twisted averages:

| eta order | primes supporting eta | specified | generic |
|---:|:---|---:|---:|
| 2 | 181, 191, 193, 197, 199 | 1.101315 | 1.105132 |
| 3 | 181, 193, 199 | 1.008533 | 1.013504 |
| 4 | 181, 193, 197 | 0.971089 | 0.969584 |
| 6 | 181, 193, 199 | 1.051302 | 1.054039 |

Against the same largest-five product/graph baselines, the generic shifted detectors favor `order 2: product, order 3: product, order 4: product, order 6: product`.

## Real-character restriction

The only nontrivial real character is quadratic, k=(p-1)/2.  Its transform was recomputed as an exact signed integer sum before division by p^(3/2).

| p | s+(quadratic) | |s+|^2 | |s+|^4 | |s+|^6 | s-(quadratic) | |s-|^2 | |s-|^4 | |s-|^6 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29 | -1.267851 | 1.607446 | 2.583883 | 4.153452 | -0.384197 | 0.147608 | 0.021788 | 0.003216 |
| 31 | 0.926995 | 0.859320 | 0.738431 | 0.634548 | 0.463498 | 0.214830 | 0.046152 | 0.009915 |
| 37 | 0.719801 | 0.518113 | 0.268442 | 0.139083 | 0.586504 | 0.343988 | 0.118327 | 0.040703 |
| 41 | 0.129510 | 0.016773 | 0.000281 | 0.000005 | 0.304729 | 0.092860 | 0.008623 | 0.000801 |
| 43 | -0.184417 | 0.034010 | 0.001157 | 0.000039 | -0.198603 | 0.039443 | 0.001556 | 0.000061 |
| 47 | -2.222113 | 4.937788 | 24.381754 | 120.391940 | -1.800036 | 3.240130 | 10.498441 | 34.016313 |
| 53 | 0.627194 | 0.393372 | 0.154741 | 0.060871 | 0.964116 | 0.929519 | 0.864006 | 0.803110 |
| 59 | 1.474003 | 2.172686 | 4.720563 | 10.256299 | 0.158875 | 0.025241 | 0.000637 | 0.000016 |
| 61 | -1.154431 | 1.332711 | 1.776118 | 2.367051 | -0.696856 | 0.485609 | 0.235816 | 0.114514 |
| 67 | -0.342804 | 0.117514 | 0.013810 | 0.001623 | 0.452209 | 0.204493 | 0.041817 | 0.008551 |
| 71 | -1.691582 | 2.861449 | 8.187889 | 23.429226 | -2.172981 | 4.721844 | 22.295813 | 105.277357 |
| 73 | -0.246909 | 0.060964 | 0.003717 | 0.000227 | 1.353188 | 1.831118 | 3.352993 | 6.139725 |
| 79 | 0.934250 | 0.872823 | 0.761821 | 0.664935 | 0.979823 | 0.960054 | 0.921703 | 0.884885 |
| 83 | -0.312101 | 0.097407 | 0.009488 | 0.000924 | 1.682170 | 2.829697 | 8.007186 | 22.657909 |
| 89 | -1.274379 | 1.624043 | 2.637516 | 4.283439 | -0.133393 | 0.017794 | 0.000317 | 0.000006 |
| 97 | 0.500346 | 0.250346 | 0.062673 | 0.015690 | -0.406138 | 0.164948 | 0.027208 | 0.004488 |
| 101 | -1.542800 | 2.380233 | 5.665507 | 13.485225 | -1.808800 | 3.271758 | 10.704403 | 35.022222 |
| 103 | 0.926018 | 0.857510 | 0.735323 | 0.630547 | 2.357137 | 5.556096 | 30.870202 | 171.517800 |
| 107 | 0.704724 | 0.496636 | 0.246647 | 0.122494 | 2.233432 | 4.988220 | 24.882339 | 124.118580 |
| 109 | 1.752207 | 3.070229 | 9.426308 | 28.940927 | 1.395439 | 1.947249 | 3.791779 | 7.383538 |
| 113 | 0.407923 | 0.166401 | 0.027689 | 0.004608 | 0.506158 | 0.256196 | 0.065636 | 0.016816 |
| 127 | -0.983778 | 0.967819 | 0.936674 | 0.906531 | -2.202321 | 4.850218 | 23.524614 | 114.099506 |
| 131 | 1.795428 | 3.223563 | 10.391361 | 33.497213 | 1.360577 | 1.851171 | 3.426833 | 6.343652 |
| 137 | -1.355747 | 1.838050 | 3.378429 | 6.209722 | -1.626398 | 2.645169 | 6.996921 | 18.508040 |
| 139 | 0.417382 | 0.174208 | 0.030348 | 0.005287 | 0.102515 | 0.010509 | 0.000110 | 0.000001 |
| 149 | -0.166046 | 0.027571 | 0.000760 | 0.000021 | -0.903904 | 0.817043 | 0.667559 | 0.545424 |
| 151 | -0.728637 | 0.530912 | 0.281868 | 0.149647 | -0.500130 | 0.250130 | 0.062565 | 0.015649 |
| 157 | -1.597190 | 2.551017 | 6.507689 | 16.601228 | -1.730374 | 2.994195 | 8.965206 | 26.843580 |
| 163 | -1.458883 | 2.128339 | 4.529825 | 9.641002 | -0.357513 | 0.127815 | 0.016337 | 0.002088 |
| 167 | -0.187200 | 0.035044 | 0.001228 | 0.000043 | 1.310402 | 1.717155 | 2.948620 | 5.063236 |
| 173 | 1.241947 | 1.542432 | 2.379096 | 3.669594 | 2.504988 | 6.274966 | 39.375203 | 247.078073 |
| 179 | -1.287760 | 1.658325 | 2.750042 | 4.560463 | -0.791697 | 0.626784 | 0.392858 | 0.246237 |
| 181 | 0.992975 | 0.986000 | 0.972196 | 0.958585 | 2.419607 | 5.854499 | 34.275161 | 200.663903 |
| 191 | 0.074252 | 0.005513 | 0.000030 | 0.000000 | 0.580375 | 0.336835 | 0.113458 | 0.038217 |
| 193 | -1.074875 | 1.155357 | 1.334849 | 1.542226 | -2.865836 | 8.213018 | 67.453668 | 553.998203 |
| 197 | -0.392763 | 0.154263 | 0.023797 | 0.003671 | -0.941763 | 0.886918 | 0.786623 | 0.697670 |
| 199 | -0.031348 | 0.000983 | 0.000001 | 0.000000 | -1.858052 | 3.452359 | 11.918782 | 41.147913 |

| largest-five real-character average | mu2 | mu4 | mu6 |
|:---|---:|---:|---:|
| plus quadratic samples | 0.460423 | 0.466175 | 0.500896 |
| minus quadratic samples | 3.748726 | 22.909538 | 159.309181 |

## Deligne-ceiling audit

| p | plus violations above 2 p^(3/2) | minus violations above 4 p^(3/2) |
|---:|:---|:---|
| 29 | none | none |
| 31 | none | none |
| 37 | none | none |
| 41 | k=10, k=30 | none |
| 43 | none | none |
| 47 | k=5, k=7, k=23, k=39, k=41 | none |
| 53 | none | none |
| 59 | none | none |
| 61 | none | none |
| 67 | none | none |
| 71 | k=2, k=68 | none |
| 73 | none | none |
| 79 | none | none |
| 83 | none | none |
| 89 | k=18, k=24, k=33, k=55, k=64, k=70 | none |
| 97 | none | none |
| 101 | none | none |
| 103 | none | none |
| 107 | none | none |
| 109 | none | none |
| 113 | k=11, k=101 | none |
| 127 | none | none |
| 131 | none | none |
| 137 | k=2, k=47, k=89, k=134 | none |
| 139 | none | none |
| 149 | none | none |
| 151 | none | none |
| 157 | none | none |
| 163 | none | none |
| 167 | k=60, k=106 | none |
| 173 | none | none |
| 179 | none | none |
| 181 | none | none |
| 191 | k=37, k=94, k=96, k=153 | none |
| 193 | none | none |
| 197 | none | none |
| 199 | none | none |

Contrary to the parenthetical expectation in the spec, the nontrivial-character audit found plus-side violations at p = 41, 47, 71, 89, 113, 137, 167, 191; minus-side violation primes: none.  These ceilings are a requested diagnostic, not one of gates g1--g6, so the characters are reported rather than silently discarded from the moments.

As an exact spot-check, at p=41 the order-four characters k=10 and k=30 both have `S_+=-574=-14p`, obtained by grouping the four integer residue-class sums in discrete-log order; thus their violation is not FFT roundoff.

## Numerical precision

All F_p and F_{p^2} operations, point counts, traces, congruences, and inversion gates were exact integers.  The FFT used NumPy complex128 after those gates passed.  Exact trivial- and quadratic-character sums plus Parseval were checked against the FFT at every prime; table entries are rounded to six decimals.

## Limitations

This is a finite-p fingerprint, not a theorem identifying a Tannakian group.  Moment proximity does not prove connectedness, Zariski density, or exclusion of finite/imprimitive subgroups.  It also does not identify arithmetic and geometric groups: finite arithmetic component groups, determinant phases, and exceptional characters can survive while leaving absolute moments nearly unchanged.  The product-versus-graph comparison is therefore evidence about the sampled Frobenius distributions, not a compatible-system or Goursat theorem.
