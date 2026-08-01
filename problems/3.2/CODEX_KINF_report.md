H0_EMPIRICAL <= 40: YES; CERTIFIED YES RANGE: h=2..60

# K_infinity branch-table scan

## Verdict

The 120-decimal-digit mpmath pass returned `YES` for `h=2..60`. The proof-grade pass then returned `YES` for exactly `h=2..60`: Arb isolated all quotient critical roots, every explicit interval-Newton inclusion succeeded, every squared critical-value ball excluded zero, and all nonmirror value balls were disjoint.

Thus the empirical failure cutoff is `none <= 60`.  In the whole tested range there are `2h-2` certified mirror-orbits (equivalently `4h-4` critical points), stronger than the required `h` good orbits.

This is consistent with every banked full-Morse certificate through `h=30` (indeed through `h=32`): no tested height is refuted.

## Exact reconstruction and counting convention

The exact recurrence used by the program is

```text
N_0=0, N_1=1,
N_(h+1)(x)=P(x+h)N_h(x)-(x+h)^6N_(h-1)(x),
delta_h(x)=N_h(x)/q_h(x)^3,
A_h(x)=q_h N_h'-3q_h' N_h.
```

With `s=2x+h+1` and `u=s^2`, reflection makes `A_h` an exact degree-`2h-2` polynomial `J_h(u)`.  Both `N_h(x)^2` and `q_h(x)^6` are also exact polynomials in `u`, so the program evaluates

```text
Psi_h(u)=delta_h(x)^2
```

without choosing a square root of `u`.  The roots of `J_h` are exactly the mirror-orbits of critical points, with multiplicity; their `Psi_h` values are exactly the scaled `U_(r,s)` of Q6723 (3.5).

There is a coordinate ambiguity in Q6723 that matters in code.  If the standard cell is `-1 < Re(z) < 0`, reflection in that same orientation is

```text
(r,s,z) -> (s-1,r+1,-1-z),
```

not `(s,r,-z)`.  The latter is correct only after switching to the opposite pole orientation, whose local coordinate lies in `0 < Re(z) < 1`.  The quotient variable above is orientation-free and therefore prevents an off-by-one double count.  The finite-cell identity was independently checked from the Q6723 `F/G` recurrence:

| h | j | (r,s) | relative error in (3.3) |
|---:|---:|:---:|---:|
| 2 | 1 | (0,1) | `4.751681e-119` |
| 7 | 2 | (1,5) | `2.4491103e-111` |
| 12 | 5 | (4,7) | `2.6254668e-99` |

The transfer anchors were also checked symbolically in the implementation: `lambda_+- = 17 +- 12 sqrt(2)` and the two normalized diagonal drifts used by the limiting object are both `-3/2`; consequently no `n^d` ratio remains in `K_infinity`.

## K_infinity bulk template (uncertified)

The code implements `K_infinity(z)=phi(-z)phi(z)+z^6 gamma(-z)gamma(z)` through the normalized `F_n/b_n,G_n/b_n` recurrence and Richardson extrapolation.  Arb-accelerated Newton found the following finite numerical roots of `H_infinity=zK_infinity'-3K_infinity`.  As an independent asymptotic diagnostic, exact central `K_(m,m)` roots were also computed for `m=10,15,20,25,30`.  Neither template calculation is used as a proof of any fixed-height row.

| Re(z) | Im(z) | residual (fit) | shorter-fit check | Newton correction |
|---:|---:|---:|---:|---:|
| `-0.5000000000000000000008359435742539470351` | `-0.1103178000763257967025342734633929459559` | `9.1606938e-102` | `8.9465807e-18` | `1.8443096e-52` |
| `-0.5000000000000000000008359435742539470351` | `0.1103178000763257967025342734633929459559` | `9.1606938e-102` | `8.9465807e-18` | `1.8443096e-52` |

The outer seed attempts did not converge to additional finite roots; their final diagnostics are retained in `CRON_kinf_results.json`.  The finite central-cell outer branches move outward with `m`, so they are not silently counted as finite `K_infinity` branches.  Without a certified outer-contour count this is evidence of escape, not an exhaustive theorem about all finite `H_infinity` zeros.

## Per-height machine table

`min sep` is the minimum distance between midpoint values.  `cert margin` is the outward-rounded lower bound on the same pairwise differences (also including distance from zero in the global minimum).  `rel margin` divides the pairwise lower bound by the larger value magnitude.  All value balls use the exact Apéry pole scales, via the algebraically equivalent `Psi_h`.

| h | critical zeros / quotient orbits | mpmath | Arb | min sep | cert margin | rel margin | bits |
|---:|:---:|:---:|:---:|---:|---:|---:|---:|
| 2 | 4 / 2 | YES | YES | `7.772047919e+02` | `4.424880405e+01` | `9.461335317e-01` | 642 |
| 3 | 8 / 4 | YES | YES | `1.989299294e+03` | `1.989299294e+03` | `4.037854270e-01` | 647 |
| 4 | 12 / 6 | YES | YES | `8.027905125e+05` | `5.525010233e+05` | `5.271339668e-01` | 652 |
| 5 | 16 / 8 | YES | YES | `2.174635554e+07` | `2.174635554e+07` | `1.270226445e-01` | 657 |
| 6 | 20 / 10 | YES | YES | `1.490291627e+10` | `1.490291627e+10` | `1.882936466e-01` | 662 |
| 7 | 24 / 12 | YES | YES | `1.509463653e+12` | `1.509463653e+12` | `6.147598305e-02` | 667 |
| 8 | 28 / 14 | YES | YES | `1.371550871e+15` | `1.371550871e+15` | `9.756935266e-02` | 672 |
| 9 | 32 / 16 | YES | YES | `2.349164873e+17` | `2.349164873e+17` | `3.612363725e-02` | 677 |
| 10 | 36 / 18 | YES | YES | `2.544947509e+20` | `2.544947509e+20` | `5.985735237e-02` | 682 |
| 11 | 40 / 20 | YES | YES | `5.780064638e+22` | `5.780064638e+22` | `2.374089424e-02` | 687 |
| 12 | 44 / 22 | YES | YES | `7.060912918e+25` | `7.060912918e+25` | `4.052216935e-02` | 692 |
| 13 | 48 / 24 | YES | YES | `1.912183116e+28` | `1.912183116e+28` | `1.678353469e-02` | 697 |
| 14 | 52 / 26 | YES | YES | `2.548810354e+31` | `2.548810354e+31` | `2.927254896e-02` | 702 |
| 15 | 56 / 28 | YES | YES | `7.783705977e+33` | `7.783705977e+33` | `1.248987336e-02` | 707 |
| 16 | 60 / 30 | YES | YES | `1.108529968e+37` | `1.108529968e+37` | `2.214487256e-02` | 712 |
| 17 | 64 / 32 | YES | YES | `3.693895068e+39` | `3.693895068e+39` | `9.655373181e-03` | 717 |
| 18 | 68 / 34 | YES | YES | `5.541326850e+42` | `5.541326850e+42` | `1.734186106e-02` | 722 |
| 19 | 72 / 36 | YES | YES | `1.972914054e+45` | `1.972914054e+45` | `7.686668154e-03` | 727 |
| 20 | 76 / 38 | YES | YES | `3.086182719e+48` | `3.086182719e+48` | `1.395054524e-02` | 732 |
| 21 | 80 / 40 | YES | YES | `1.157408383e+51` | `1.157408383e+51` | `6.263971179e-03` | 737 |
| 22 | 84 / 42 | YES | YES | `1.873976857e+54` | `1.873976857e+54` | `1.146649345e-02` | 742 |
| 23 | 88 / 44 | YES | YES | `7.328496339e+56` | `7.328496339e+56` | `5.202574621e-03` | 747 |
| 24 | 92 / 46 | YES | YES | `1.221298925e+60` | `1.221298925e+60` | `9.592291556e-03` | 752 |
| 25 | 96 / 48 | YES | YES | `4.943535169e+62` | `4.943535169e+62` | `4.389747628e-03` | 757 |
| 26 | 100 / 50 | YES | YES | `8.442846553e+65` | `8.442846553e+65` | `8.143255184e-03` | 762 |
| 27 | 104 / 52 | YES | YES | `3.517504080e+68` | `3.517504080e+68` | `3.753526357e-03` | 767 |
| 28 | 108 / 54 | YES | YES | `6.135364508e+71` | `6.135364508e+71` | `6.999737003e-03` | 772 |
| 29 | 112 / 56 | YES | YES | `2.619573724e+74` | `2.619573724e+74` | `3.246222597e-03` | 777 |
| 30 | 116 / 58 | YES | YES | `4.653646691e+77` | `4.653646691e+77` | `6.081436086e-03` | 782 |
| 31 | 120 / 60 | YES | YES | `2.029263354e+80` | `2.029263354e+80` | `2.835221333e-03` | 787 |
| 32 | 124 / 62 | YES | YES | `3.663373817e+83` | `3.663373817e+83` | `5.332821206e-03` | 792 |
| 33 | 128 / 64 | YES | YES | `1.626984611e+86` | `1.626984611e+86` | `2.497605165e-03` | 797 |
| 34 | 132 / 66 | YES | YES | `2.979194292e+89` | `2.979194292e+89` | `4.714490788e-03` | 802 |
| 35 | 136 / 68 | YES | YES | `1.344563160e+92` | `1.344563160e+92` | `2.216889281e-03` | 807 |
| 36 | 140 / 70 | YES | YES | `2.493412321e+95` | `2.493412321e+95` | `4.197857177e-03` | 812 |
| 37 | 144 / 72 | YES | YES | `1.141429612e+98` | `1.141429612e+98` | `1.980970480e-03` | 817 |
| 38 | 148 / 74 | YES | YES | `2.140865951e+101` | `2.140865951e+101` | `3.761769091e-03` | 822 |
| 39 | 152 / 76 | YES | YES | `9.925264573e+103` | `9.925264573e+103` | `1.780800471e-03` | 827 |
| 40 | 156 / 78 | YES | YES | `1.880719820e+107` | `1.880719820e+107` | `3.390304335e-03` | 832 |
| 41 | 160 / 80 | YES | YES | `8.818692910e+109` | `8.818692910e+109` | `1.609505225e-03` | 837 |
| 42 | 164 / 82 | YES | YES | `1.686594091e+113` | `1.686594091e+113` | `3.071291085e-03` | 842 |
| 43 | 168 / 84 | YES | YES | `7.989735536e+115` | `7.989735536e+115` | `1.461786873e-03` | 847 |
| 44 | 172 / 86 | YES | YES | `1.540997564e+119` | `1.540997564e+119` | `2.795295221e-03` | 852 |
| 45 | 176 / 88 | YES | YES | `7.367975291e+121` | `7.367975291e+121` | `1.333509112e-03` | 857 |
| 46 | 180 / 90 | YES | YES | `1.432076680e+125` | `1.432076680e+125` | `2.554913008e-03` | 862 |
| 47 | 184 / 92 | YES | YES | `6.905191072e+127` | `6.905191072e+127` | `1.221404479e-03` | 867 |
| 48 | 188 / 94 | YES | YES | `1.351657347e+131` | `1.351657347e+131` | `2.344268236e-03` | 872 |
| 49 | 192 / 96 | YES | YES | `6.567871113e+133` | `6.567871113e+133` | `1.122864207e-03` | 877 |
| 50 | 196 / 98 | YES | YES | `1.294032777e+137` | `1.294032777e+137` | `2.158648929e-03` | 882 |
| 51 | 200 / 100 | YES | YES | `6.332513190e+139` | `6.332513190e+139` | `1.035785074e-03` | 887 |
| 52 | 204 / 102 | YES | YES | `1.255194973e+143` | `1.255194973e+143` | `1.994241039e-03` | 892 |
| 53 | 208 / 104 | YES | YES | `6.182608938e+145` | `6.182608938e+145` | `9.584562429e-04` | 897 |
| 54 | 212 / 106 | YES | YES | `1.232337483e+149` | `1.232337483e+149` | `1.847930606e-03` | 902 |
| 55 | 216 / 108 | YES | YES | `6.106673331e+151` | `6.106673331e+151` | `8.894745841e-04` | 907 |
| 56 | 220 / 110 | YES | YES | `1.223528376e+155` | `1.223528376e+155` | `1.717154983e-03` | 912 |
| 57 | 224 / 112 | YES | YES | `6.096939394e+157` | `6.096939394e+157` | `8.276805765e-04` | 917 |
| 58 | 228 / 114 | YES | YES | `1.227492776e+161` | `1.227492776e+161` | `1.599789695e-03` | 922 |
| 59 | 232 / 116 | YES | YES | `6.148487033e+163` | `6.148487033e+163` | `7.721092616e-04` | 927 |
| 60 | 236 / 118 | YES | YES | `1.243467684e+167` | `1.243467684e+167` | `1.494061488e-03` | 932 |

## Inline machine output

```text
h=02 UNCERTIFIED=YES        CERTIFIED=YES        zeros=4 orbits=2 good=2 rel_margin=9.461335317e-01
h=03 UNCERTIFIED=YES        CERTIFIED=YES        zeros=8 orbits=4 good=4 rel_margin=4.037854270e-01
h=04 UNCERTIFIED=YES        CERTIFIED=YES        zeros=12 orbits=6 good=6 rel_margin=5.271339668e-01
h=05 UNCERTIFIED=YES        CERTIFIED=YES        zeros=16 orbits=8 good=8 rel_margin=1.270226445e-01
h=06 UNCERTIFIED=YES        CERTIFIED=YES        zeros=20 orbits=10 good=10 rel_margin=1.882936466e-01
h=07 UNCERTIFIED=YES        CERTIFIED=YES        zeros=24 orbits=12 good=12 rel_margin=6.147598305e-02
h=08 UNCERTIFIED=YES        CERTIFIED=YES        zeros=28 orbits=14 good=14 rel_margin=9.756935266e-02
h=09 UNCERTIFIED=YES        CERTIFIED=YES        zeros=32 orbits=16 good=16 rel_margin=3.612363725e-02
h=10 UNCERTIFIED=YES        CERTIFIED=YES        zeros=36 orbits=18 good=18 rel_margin=5.985735237e-02
h=11 UNCERTIFIED=YES        CERTIFIED=YES        zeros=40 orbits=20 good=20 rel_margin=2.374089424e-02
h=12 UNCERTIFIED=YES        CERTIFIED=YES        zeros=44 orbits=22 good=22 rel_margin=4.052216935e-02
h=13 UNCERTIFIED=YES        CERTIFIED=YES        zeros=48 orbits=24 good=24 rel_margin=1.678353469e-02
h=14 UNCERTIFIED=YES        CERTIFIED=YES        zeros=52 orbits=26 good=26 rel_margin=2.927254896e-02
h=15 UNCERTIFIED=YES        CERTIFIED=YES        zeros=56 orbits=28 good=28 rel_margin=1.248987336e-02
h=16 UNCERTIFIED=YES        CERTIFIED=YES        zeros=60 orbits=30 good=30 rel_margin=2.214487256e-02
h=17 UNCERTIFIED=YES        CERTIFIED=YES        zeros=64 orbits=32 good=32 rel_margin=9.655373181e-03
h=18 UNCERTIFIED=YES        CERTIFIED=YES        zeros=68 orbits=34 good=34 rel_margin=1.734186106e-02
h=19 UNCERTIFIED=YES        CERTIFIED=YES        zeros=72 orbits=36 good=36 rel_margin=7.686668154e-03
h=20 UNCERTIFIED=YES        CERTIFIED=YES        zeros=76 orbits=38 good=38 rel_margin=1.395054524e-02
h=21 UNCERTIFIED=YES        CERTIFIED=YES        zeros=80 orbits=40 good=40 rel_margin=6.263971179e-03
h=22 UNCERTIFIED=YES        CERTIFIED=YES        zeros=84 orbits=42 good=42 rel_margin=1.146649345e-02
h=23 UNCERTIFIED=YES        CERTIFIED=YES        zeros=88 orbits=44 good=44 rel_margin=5.202574621e-03
h=24 UNCERTIFIED=YES        CERTIFIED=YES        zeros=92 orbits=46 good=46 rel_margin=9.592291556e-03
h=25 UNCERTIFIED=YES        CERTIFIED=YES        zeros=96 orbits=48 good=48 rel_margin=4.389747628e-03
h=26 UNCERTIFIED=YES        CERTIFIED=YES        zeros=100 orbits=50 good=50 rel_margin=8.143255184e-03
h=27 UNCERTIFIED=YES        CERTIFIED=YES        zeros=104 orbits=52 good=52 rel_margin=3.753526357e-03
h=28 UNCERTIFIED=YES        CERTIFIED=YES        zeros=108 orbits=54 good=54 rel_margin=6.999737003e-03
h=29 UNCERTIFIED=YES        CERTIFIED=YES        zeros=112 orbits=56 good=56 rel_margin=3.246222597e-03
h=30 UNCERTIFIED=YES        CERTIFIED=YES        zeros=116 orbits=58 good=58 rel_margin=6.081436086e-03
h=31 UNCERTIFIED=YES        CERTIFIED=YES        zeros=120 orbits=60 good=60 rel_margin=2.835221333e-03
h=32 UNCERTIFIED=YES        CERTIFIED=YES        zeros=124 orbits=62 good=62 rel_margin=5.332821206e-03
h=33 UNCERTIFIED=YES        CERTIFIED=YES        zeros=128 orbits=64 good=64 rel_margin=2.497605165e-03
h=34 UNCERTIFIED=YES        CERTIFIED=YES        zeros=132 orbits=66 good=66 rel_margin=4.714490788e-03
h=35 UNCERTIFIED=YES        CERTIFIED=YES        zeros=136 orbits=68 good=68 rel_margin=2.216889281e-03
h=36 UNCERTIFIED=YES        CERTIFIED=YES        zeros=140 orbits=70 good=70 rel_margin=4.197857177e-03
h=37 UNCERTIFIED=YES        CERTIFIED=YES        zeros=144 orbits=72 good=72 rel_margin=1.980970480e-03
h=38 UNCERTIFIED=YES        CERTIFIED=YES        zeros=148 orbits=74 good=74 rel_margin=3.761769091e-03
h=39 UNCERTIFIED=YES        CERTIFIED=YES        zeros=152 orbits=76 good=76 rel_margin=1.780800471e-03
h=40 UNCERTIFIED=YES        CERTIFIED=YES        zeros=156 orbits=78 good=78 rel_margin=3.390304335e-03
h=41 UNCERTIFIED=YES        CERTIFIED=YES        zeros=160 orbits=80 good=80 rel_margin=1.609505225e-03
h=42 UNCERTIFIED=YES        CERTIFIED=YES        zeros=164 orbits=82 good=82 rel_margin=3.071291085e-03
h=43 UNCERTIFIED=YES        CERTIFIED=YES        zeros=168 orbits=84 good=84 rel_margin=1.461786873e-03
h=44 UNCERTIFIED=YES        CERTIFIED=YES        zeros=172 orbits=86 good=86 rel_margin=2.795295221e-03
h=45 UNCERTIFIED=YES        CERTIFIED=YES        zeros=176 orbits=88 good=88 rel_margin=1.333509112e-03
h=46 UNCERTIFIED=YES        CERTIFIED=YES        zeros=180 orbits=90 good=90 rel_margin=2.554913008e-03
h=47 UNCERTIFIED=YES        CERTIFIED=YES        zeros=184 orbits=92 good=92 rel_margin=1.221404479e-03
h=48 UNCERTIFIED=YES        CERTIFIED=YES        zeros=188 orbits=94 good=94 rel_margin=2.344268236e-03
h=49 UNCERTIFIED=YES        CERTIFIED=YES        zeros=192 orbits=96 good=96 rel_margin=1.122864207e-03
h=50 UNCERTIFIED=YES        CERTIFIED=YES        zeros=196 orbits=98 good=98 rel_margin=2.158648929e-03
h=51 UNCERTIFIED=YES        CERTIFIED=YES        zeros=200 orbits=100 good=100 rel_margin=1.035785074e-03
h=52 UNCERTIFIED=YES        CERTIFIED=YES        zeros=204 orbits=102 good=102 rel_margin=1.994241039e-03
h=53 UNCERTIFIED=YES        CERTIFIED=YES        zeros=208 orbits=104 good=104 rel_margin=9.584562429e-04
h=54 UNCERTIFIED=YES        CERTIFIED=YES        zeros=212 orbits=106 good=106 rel_margin=1.847930606e-03
h=55 UNCERTIFIED=YES        CERTIFIED=YES        zeros=216 orbits=108 good=108 rel_margin=8.894745841e-04
h=56 UNCERTIFIED=YES        CERTIFIED=YES        zeros=220 orbits=110 good=110 rel_margin=1.717154983e-03
h=57 UNCERTIFIED=YES        CERTIFIED=YES        zeros=224 orbits=112 good=112 rel_margin=8.276805765e-04
h=58 UNCERTIFIED=YES        CERTIFIED=YES        zeros=228 orbits=114 good=114 rel_margin=1.599789695e-03
h=59 UNCERTIFIED=YES        CERTIFIED=YES        zeros=232 orbits=116 good=116 rel_margin=7.721092616e-04
h=60 UNCERTIFIED=YES        CERTIFIED=YES        zeros=236 orbits=118 good=118 rel_margin=1.494061488e-03
```

## Caveats

1. `H0_EMPIRICAL` is the finite-scan convention requested in `CODEX_SPEC_CRON_kinf.md`: the largest failed height in `2..60`.  It is not the universal tail threshold `h0` of Q6723 Section 9.
2. The fixed-height Arb rows are genuine finite certificates, but this run does not certify the all-`h` predicates `PERSIST`, `EXHAUST`, `SAME`, `CROSS`, and `COUNT` on inverse-length boxes.  Therefore the finite result alone does not close the all-height campaign prize.
3. The displayed `K_infinity` roots use accelerated normalized recurrences and central-branch extrapolation without the explicit Q6708 conjugation-tail ball or an outer-contour count.  They are intentionally labelled uncertified; every fixed-height verdict is recomputed from exact integer polynomials and does not depend on them.
4. `python-flint` was found in the pinned uv cache rather than the system Python import path.  The script re-executes through uv and records the exact versions in the JSON payload.
5. The uncertified pass takes the midpoints of FLINT root balls as Newton seeds, then discards every radius and recomputes roots and values in mpmath.  It is a 100+-digit floating-point regression pass, not an implementation-independent root finder; the Arb pass is the certificate.

Reproduction command:

```bash
python3 -u CRON_kinf_branch.py --min-h 2 --max-h 60 --digits 120 --cert-digits 110
```
