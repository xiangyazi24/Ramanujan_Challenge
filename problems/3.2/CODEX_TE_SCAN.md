# TE_{5/4} Apéry collision-energy scan

## Coverage and conventions

- Exact range tested: every prime `p` in `[1000, 100000]` (9424 primes; first `1009`, last `99991`).
- Run started: `2026-08-01 08:42:52 CDT`; elapsed: `2339.5` seconds.
- The histogram uses exactly `1 <= r <= p-2`. `E_off` counts ordered pairs `(r,s)` with `r != s`; all energy computations are integer-exact.
- The involution `r -> p-1-r` has one fixed point and `(p-3)/2` two-cycles, so its ordered forced contribution is `E_refl=p-3`.
- Gap counts are directed and modular: an unordered pair `r<s` adds one count at `h=s-r` and one at `p-h`. Both are removed as forced exactly when `r+s=p-1`.

The incremental per-prime ledger, including all six requested normalizations, is `research/scripts/te_scan_progress.txt`.

## Dyadic block maxima: total energy

| block | coverage | primes | max E/p^(5/4) (p) | max E/p (p) | max E/(p log^2 p) (p) |
|---|---:|---:|---:|---:|---:|
| (512, 1024] | 1009-1021 | 4 | 0.36906038 (1021) | 2.08619 (1021) | 0.043458094 (1021) |
| (1024, 2048] | 1031-2039 | 137 | 0.37738904 (1193) | 2.217938 (1193) | 0.044194124 (1193) |
| (2048, 4096] | 2053-4093 | 255 | 0.313703 (2137) | 2.1469008 (2791) | 0.036282812 (2137) |
| (4096, 8192] | 4099-8191 | 464 | 0.26126776 (4133) | 2.1386436 (5323) | 0.030213442 (4133) |
| (8192, 16384] | 8209-16381 | 872 | 0.21748882 (8719) | 2.1016172 (8719) | 0.025528595 (8719) |
| (16384, 32768] | 16411-32749 | 1612 | 0.17966135 (16427) | 2.0600367 (17989) | 0.021587515 (16427) |
| (32768, 65536] | 32771-65521 | 3030 | 0.15105502 (32831) | 2.0453206 (47837) | 0.018802358 (32831) |
| (65536, 131072] | 65537-99991 | 3050 | 0.12633333 (65551) | 2.0365058 (70619) | 0.016434404 (65551) |

The first and last rows are partial dyadic blocks when the stated scan endpoints cut them. The regression excludes those partial rows. Each complete block contributes the prime attaining the largest raw `E_off` in that block; ordinary least squares on `(log p, log E_off)` gives

- fitted slope: **0.985580** from 6 complete-block maxima;
- danger threshold `slope > 1.25`: **not triggered**.

## Reflection-separated block maxima

| block | max (E-E_refl)/p^(5/4) (p) | max (E-E_refl)/p (p) | max (E-E_refl)/(p log^2 p) (p) |
|---|---:|---:|---:|
| (512, 1024] | 0.19267378 (1021) | 1.0891283 (1021) | 0.022687982 (1021) |
| (1024, 2048] | 0.20766381 (1193) | 1.2204526 (1193) | 0.02431846 (1193) |
| (2048, 4096] | 0.16683108 (2137) | 1.1479756 (2791) | 0.019295642 (2137) |
| (4096, 8192] | 0.136639 (4133) | 1.1392072 (5323) | 0.015801163 (4133) |
| (8192, 16384] | 0.11403801 (8719) | 1.1019612 (8719) | 0.013385655 (8719) |
| (16384, 32768] | 0.09196327 (17359) | 1.0602035 (17989) | 0.011077197 (17359) |
| (32768, 65536] | 0.076772047 (32831) | 1.0453833 (47837) | 0.0095560907 (32831) |
| (65536, 131072] | 0.063839767 (65551) | 1.0365482 (70619) | 0.0083151066 (70619) |

## Gap spectrum

Inspected primes are the largest 50 tested primes together with the final five spike primes by `E_off/p`. The danger test is exact: `C_p(h)^4 > p`, equivalent to `C_p(h) > p^(1/4)`.

| p | top five non-reflection (h, C_p(h)) | max > p^(1/4)? |
|---:|---|:---:|
| 1069 | (528, 8), (541, 8), (51, 6), (111, 6), (115, 6) | YES |
| 1193 | (234, 10), (959, 10), (276, 8), (917, 8), (34, 6) | YES |
| 1223 | (463, 10), (760, 10), (298, 8), (439, 8), (526, 8) | YES |
| 1231 | (363, 8), (868, 8), (61, 6), (127, 6), (184, 6) | YES |
| 1499 | (184, 8), (335, 8), (1164, 8), (1315, 8), (26, 6) | YES |
| 99431 | (5836, 12), (11062, 12), (35793, 12), (46259, 12), (47422, 12) | no |
| 99439 | (15373, 12), (34965, 12), (64474, 12), (84066, 12), (4966, 10) | no |
| 99469 | (38472, 14), (60997, 14), (1115, 10), (21600, 10), (27048, 10) | no |
| 99487 | (1092, 10), (3004, 10), (7508, 10), (14735, 10), (27423, 10) | no |
| 99497 | (24541, 12), (74956, 12), (4456, 10), (5334, 10), (14827, 10) | no |
| 99523 | (48326, 12), (51197, 12), (18392, 10), (20637, 10), (49383, 10) | no |
| 99527 | (712, 10), (8810, 10), (19763, 10), (21299, 10), (30202, 10) | no |
| 99529 | (30152, 12), (69377, 12), (6869, 10), (7915, 10), (29850, 10) | no |
| 99551 | (5555, 14), (93996, 14), (3220, 10), (16979, 10), (18151, 10) | no |
| 99559 | (35465, 12), (64094, 12), (8610, 10), (8856, 10), (10902, 10) | no |
| 99563 | (2700, 10), (10062, 10), (11460, 10), (15320, 10), (19040, 10) | no |
| 99571 | (9874, 10), (15801, 10), (20102, 10), (42556, 10), (57015, 10) | no |
| 99577 | (11533, 10), (15569, 10), (25202, 10), (36978, 10), (37408, 10) | no |
| 99581 | (33958, 12), (65623, 12), (14893, 10), (21323, 10), (25473, 10) | no |
| 99607 | (1875, 12), (97732, 12), (3165, 10), (5964, 10), (7815, 10) | no |
| 99611 | (1497, 10), (11750, 10), (12928, 10), (13540, 10), (18580, 10) | no |
| 99623 | (1110, 10), (18774, 10), (24210, 10), (75413, 10), (80849, 10) | no |
| 99643 | (6837, 12), (92806, 12), (5056, 10), (17594, 10), (20346, 10) | no |
| 99661 | (22435, 10), (23604, 10), (25158, 10), (32632, 10), (47548, 10) | no |
| 99667 | (3049, 10), (6215, 10), (31946, 10), (33825, 10), (34474, 10) | no |
| 99679 | (14514, 12), (85165, 12), (13960, 10), (15858, 10), (16845, 10) | no |
| 99689 | (193, 10), (3489, 10), (8433, 10), (19991, 10), (21319, 10) | no |
| 99707 | (7048, 10), (9264, 10), (10341, 10), (29231, 10), (31461, 10) | no |
| 99709 | (21602, 12), (78107, 12), (6138, 10), (9754, 10), (22964, 10) | no |
| 99713 | (431, 10), (8858, 10), (10432, 10), (14298, 10), (32358, 10) | no |
| 99719 | (7384, 10), (8886, 10), (20101, 10), (26253, 10), (28249, 10) | no |
| 99721 | (6858, 10), (25906, 10), (28454, 10), (33211, 10), (41399, 10) | no |
| 99733 | (14276, 12), (85457, 12), (11542, 10), (21267, 10), (23799, 10) | no |
| 99761 | (35877, 12), (63884, 12), (4277, 10), (6970, 10), (27729, 10) | no |
| 99767 | (6438, 10), (7859, 10), (9140, 10), (11602, 10), (16286, 10) | no |
| 99787 | (5455, 10), (12223, 10), (15184, 10), (17337, 10), (17795, 10) | no |
| 99793 | (27397, 12), (72396, 12), (2109, 10), (3714, 10), (3927, 10) | no |
| 99809 | (41646, 12), (58163, 12), (1164, 10), (3708, 10), (4511, 10) | no |
| 99817 | (41812, 14), (58005, 14), (39153, 12), (60664, 12), (8045, 10) | no |
| 99823 | (1964, 12), (23844, 12), (35021, 12), (42149, 12), (57674, 12) | no |
| 99829 | (47315, 12), (52514, 12), (14558, 10), (18345, 10), (19694, 10) | no |
| 99833 | (745, 10), (13884, 10), (14599, 10), (22729, 10), (23481, 10) | no |
| 99839 | (158, 10), (10714, 10), (13408, 10), (16541, 10), (21723, 10) | no |
| 99859 | (18559, 10), (27070, 10), (27213, 10), (27364, 10), (41359, 10) | no |
| 99871 | (795, 12), (99076, 12), (4446, 10), (16395, 10), (18285, 10) | no |
| 99877 | (277, 12), (33948, 12), (65929, 12), (99600, 12), (497, 10) | no |
| 99881 | (7761, 12), (92120, 12), (288, 10), (3017, 10), (28556, 10) | no |
| 99901 | (27782, 12), (40402, 12), (59499, 12), (72119, 12), (10920, 10) | no |
| 99907 | (46353, 12), (53554, 12), (10641, 10), (27430, 10), (29237, 10) | no |
| 99923 | (5332, 10), (8704, 10), (16406, 10), (16960, 10), (23197, 10) | no |
| 99929 | (15873, 10), (23518, 10), (24584, 10), (27556, 10), (27985, 10) | no |
| 99961 | (6641, 10), (7451, 10), (17806, 10), (20029, 10), (21879, 10) | no |
| 99971 | (1832, 10), (6470, 10), (27632, 10), (28702, 10), (39213, 10) | no |
| 99989 | (38964, 14), (61025, 14), (4205, 10), (13767, 10), (15308, 10) | no |
| 99991 | (5572, 10), (7346, 10), (12885, 10), (15294, 10), (19378, 10) | no |

Danger-signal primes: 1069, 1193, 1223, 1231, 1499.
Among the largest 50 tested primes: none.

## Spike autopsies

A spike is labeled *concentrated* when one non-reflection gap accounts for at least 25% of all nonforced unordered pairs, and *diffuse* otherwise. An AP is reported as structural only at length at least 3; length 2 is tautological.

| rank | p | E/p | Z_p | max N_p(a), a | gap diagnosis | popular-fiber structure |
|---:|---:|---:|---:|---|---|---|
| 1 | 1193 | 2.217938 | 2 | 10, a=136 | diffuse; h=234, count=10, share=1.374% | positions [155, 156, 348, 355, 564, 628, 837, 844, 1036, 1037]; reflection pairs=5, center=no; longest AP=2 (step 1); adjacent-gap mode=1 x2 |
| 2 | 1223 | 2.1978741 | 0 | 8, a=11 | diffuse; h=463, count=10, share=1.362% | positions [302, 309, 508, 595, 627, 714, 913, 920]; reflection pairs=4, center=no; longest AP=4 (step 206); adjacent-gap mode=7 x2 |
| 3 | 1231 | 2.1575955 | 0 | 12, a=569 | diffuse; h=363, count=8, share=1.120% | positions [41, 321, 421, 523, 591, 603, 627, 639, 707, 809, 909, 1189]; reflection pairs=6, center=no; longest AP=2 (step 12); adjacent-gap mode=12 x2 |
| 4 | 1499 | 2.1534356 | 2 | 8, a=100 | diffuse; h=184, count=8, share=0.924% | positions [107, 202, 405, 613, 885, 1093, 1296, 1391]; reflection pairs=4, center=no; longest AP=2 (step 95); adjacent-gap mode=95 x2 |
| 5 | 1069 | 2.1496726 | 6 | 10, a=851 | diffuse; h=528, count=8, share=1.299% | positions [29, 155, 261, 369, 521, 547, 699, 807, 913, 1039]; reflection pairs=5, center=no; longest AP=2 (step 26); adjacent-gap mode=106 x2 |

## Envelope test and verdict

Violations of `E_off <= 100 p log^2 p`: 0.
Violating primes: none.

Globally, `max E_off/p^(5/4) = 0.37738904` at `p=1193`, while `max E_off/p = 2.217938` at `p=1193`. After subtracting reflection, `max (E_off-E_refl)/p^(5/4) = 0.20766381` at `p=1193`.

**Verdict.** No genuine resonance was found: the fitted block-max slope stays below 1.25, the 100 p log^2 p envelope has no violation, and no p^(1/4) gap signal survives at the top of the tested range. The 5 gap-threshold crossings are confined to low-p spike primes (largest p=1499), are all diffuse under the stated 25% test, and do not persist among the largest 50 primes; they are not a growing resonance. The data therefore support attacking the stronger p^(1+eps) collision-energy theorem; the 5/4 target is comfortably supported as the fallback.
