# Last-wall empirical ground truth

All finite-field values and counts below are exact integers produced by pure Python 3 standard-library code in `CODEX_lastwall_empirics.py`. Here `N=p-2`, `r` ranges over `1,...,N`, and `d_D(r)` counts admissible nonwrapping collisions with `r+d<=N`. The cutoffs `p^0.6` and `p^0.66` were ceiled by exact integer-power comparisons; `sqrt(p) log(p)` (natural logarithm) was evaluated with 80-digit Decimal arithmetic. Ratios only are rounded.

For each prime, the orbit was generated from both Apery recurrences. Every regular value was checked against

`Delta_(r,d) prod_(j=1)^d (r+j)^3 = N_d(r) (mod p)`,

and all collision pairs were independently reconstructed by grouping equal normalized projective orbit values. `R_d` in Section 4 means the full number of roots of `N_d` in `F_p`, including residues outside the nonwrapping window; it is therefore distinct from `C_d`.

## 1. Exact window statistics and distributions

Histogram notation is `multiplicity:number of bases`. Deciles use nearest rank: `q_j` is the entry of rank `ceil(jN/10)` in the sorted list (with `q_0` the minimum). Thus the histogram plus deciles and top ten give the requested full distribution summary.

### p=997

`N=995`, `|Z_p intersect [1,N]|=0`, Delta-identity checks `193815` (all passed).

#### ceil(sqrt(p) log(p)): D=219

`S_D=313`, `Q_D=58`, `max_r d_D(r)=4`.

- Full histogram: `0:729, 1:228, 2:31, 3:5, 4:2`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=1, q90=1, q100=4`.
- Top ten (ties ordered by increasing `r`):

  1. r=347 (m=4, r/p=0.3480, Z_p=no, small=no, gaps=[100, 126, 176, 202], exact-mirror gaps=none, min |2r+d-(p-1)|=100)
  2. r=447 (m=4, r/p=0.4483, Z_p=no, small=no, gaps=[26, 76, 102, 202], exact-mirror gaps=[102], min |2r+d-(p-1)|=0)
  3. r=394 (m=3, r/p=0.3952, Z_p=no, small=no, gaps=[61, 147, 208], exact-mirror gaps=[208], min |2r+d-(p-1)|=0)
  4. r=399 (m=3, r/p=0.4002, Z_p=no, small=no, gaps=[87, 111, 198], exact-mirror gaps=[198], min |2r+d-(p-1)|=0)
  5. r=406 (m=3, r/p=0.4072, Z_p=no, small=no, gaps=[10, 174, 184], exact-mirror gaps=[184], min |2r+d-(p-1)|=0)
  6. r=445 (m=3, r/p=0.4463, Z_p=no, small=no, gaps=[22, 84, 106], exact-mirror gaps=[106], min |2r+d-(p-1)|=0)
  7. r=473 (m=3, r/p=0.4744, Z_p=no, small=no, gaps=[50, 76, 176], exact-mirror gaps=[50], min |2r+d-(p-1)|=0)
  8. r=7 (m=2, r/p=0.0070, Z_p=no, small=yes, gaps=[62, 73], exact-mirror gaps=none, min |2r+d-(p-1)|=909)
  9. r=26 (m=2, r/p=0.0261, Z_p=no, small=yes, gaps=[102, 194], exact-mirror gaps=none, min |2r+d-(p-1)|=750)
  10. r=79 (m=2, r/p=0.0792, Z_p=no, small=no, gaps=[18, 80], exact-mirror gaps=none, min |2r+d-(p-1)|=758)

#### ceil(p^0.6): D=63

`S_D=87`, `Q_D=0`, `max_r d_D(r)=1`.

- Full histogram: `0:908, 1:87`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=0, q100=1`.
- Top ten (ties ordered by increasing `r`):

  1. r=7 (m=1, r/p=0.0070, Z_p=no, small=yes, gaps=[62], exact-mirror gaps=none, min |2r+d-(p-1)|=920)
  2. r=18 (m=1, r/p=0.0181, Z_p=no, small=yes, gaps=[44], exact-mirror gaps=none, min |2r+d-(p-1)|=916)
  3. r=69 (m=1, r/p=0.0692, Z_p=no, small=no, gaps=[11], exact-mirror gaps=none, min |2r+d-(p-1)|=847)
  4. r=79 (m=1, r/p=0.0792, Z_p=no, small=no, gaps=[18], exact-mirror gaps=none, min |2r+d-(p-1)|=820)
  5. r=86 (m=1, r/p=0.0863, Z_p=no, small=no, gaps=[20], exact-mirror gaps=none, min |2r+d-(p-1)|=804)
  6. r=97 (m=1, r/p=0.0973, Z_p=no, small=no, gaps=[62], exact-mirror gaps=none, min |2r+d-(p-1)|=740)
  7. r=115 (m=1, r/p=0.1153, Z_p=no, small=no, gaps=[52], exact-mirror gaps=none, min |2r+d-(p-1)|=714)
  8. r=144 (m=1, r/p=0.1444, Z_p=no, small=no, gaps=[39], exact-mirror gaps=none, min |2r+d-(p-1)|=669)
  9. r=180 (m=1, r/p=0.1805, Z_p=no, small=no, gaps=[33], exact-mirror gaps=none, min |2r+d-(p-1)|=603)
  10. r=198 (m=1, r/p=0.1986, Z_p=no, small=no, gaps=[56], exact-mirror gaps=none, min |2r+d-(p-1)|=544)

#### ceil(p^0.66): D=96

`S_D=136`, `Q_D=8`, `max_r d_D(r)=2`.

- Full histogram: `0:867, 1:120, 2:8`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=1, q100=2`.
- Top ten (ties ordered by increasing `r`):

  1. r=7 (m=2, r/p=0.0070, Z_p=no, small=yes, gaps=[62, 73], exact-mirror gaps=none, min |2r+d-(p-1)|=909)
  2. r=79 (m=2, r/p=0.0792, Z_p=no, small=no, gaps=[18, 80], exact-mirror gaps=none, min |2r+d-(p-1)|=758)
  3. r=445 (m=2, r/p=0.4463, Z_p=no, small=no, gaps=[22, 84], exact-mirror gaps=none, min |2r+d-(p-1)|=22)
  4. r=447 (m=2, r/p=0.4483, Z_p=no, small=no, gaps=[26, 76], exact-mirror gaps=none, min |2r+d-(p-1)|=26)
  5. r=467 (m=2, r/p=0.4684, Z_p=no, small=no, gaps=[62, 84], exact-mirror gaps=[62], min |2r+d-(p-1)|=0)
  6. r=473 (m=2, r/p=0.4744, Z_p=no, small=no, gaps=[50, 76], exact-mirror gaps=[50], min |2r+d-(p-1)|=0)
  7. r=837 (m=2, r/p=0.8395, Z_p=no, small=no, gaps=[62, 80], exact-mirror gaps=none, min |2r+d-(p-1)|=740)
  8. r=916 (m=2, r/p=0.9188, Z_p=no, small=no, gaps=[11, 73], exact-mirror gaps=none, min |2r+d-(p-1)|=847)
  9. r=17 (m=1, r/p=0.0171, Z_p=no, small=yes, gaps=[82], exact-mirror gaps=none, min |2r+d-(p-1)|=880)
  10. r=18 (m=1, r/p=0.0181, Z_p=no, small=yes, gaps=[44], exact-mirror gaps=none, min |2r+d-(p-1)|=916)

### p=1999

`N=1997`, `|Z_p intersect [1,N]|=0`, Delta-identity checks `621010` (all passed).

#### ceil(sqrt(p) log(p)): D=340

`S_D=498`, `Q_D=84`, `max_r d_D(r)=4`.

- Full histogram: `0:1566, 1:379, 2:39, 3:11, 4:2`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=1, q90=1, q100=4`.
- Top ten (ties ordered by increasing `r`):

  1. r=810 (m=4, r/p=0.4052, Z_p=no, small=no, gaps=[169, 186, 192, 209], exact-mirror gaps=none, min |2r+d-(p-1)|=169)
  2. r=979 (m=4, r/p=0.4897, Z_p=no, small=no, gaps=[17, 23, 40, 209], exact-mirror gaps=[40], min |2r+d-(p-1)|=0)
  3. r=714 (m=3, r/p=0.3572, Z_p=no, small=no, gaps=[213, 230, 340], exact-mirror gaps=none, min |2r+d-(p-1)|=230)
  4. r=843 (m=3, r/p=0.4217, Z_p=no, small=no, gaps=[66, 246, 312], exact-mirror gaps=[312], min |2r+d-(p-1)|=0)
  5. r=857 (m=3, r/p=0.4287, Z_p=no, small=no, gaps=[59, 225, 284], exact-mirror gaps=[284], min |2r+d-(p-1)|=0)
  6. r=864 (m=3, r/p=0.4322, Z_p=no, small=no, gaps=[67, 203, 270], exact-mirror gaps=[270], min |2r+d-(p-1)|=0)
  7. r=884 (m=3, r/p=0.4422, Z_p=no, small=no, gaps=[113, 117, 230], exact-mirror gaps=[230], min |2r+d-(p-1)|=0)
  8. r=894 (m=3, r/p=0.4472, Z_p=no, small=no, gaps=[19, 191, 210], exact-mirror gaps=[210], min |2r+d-(p-1)|=0)
  9. r=912 (m=3, r/p=0.4562, Z_p=no, small=no, gaps=[16, 158, 174], exact-mirror gaps=[174], min |2r+d-(p-1)|=0)
  10. r=920 (m=3, r/p=0.4602, Z_p=no, small=no, gaps=[65, 93, 158], exact-mirror gaps=[158], min |2r+d-(p-1)|=0)

#### ceil(p^0.6): D=96

`S_D=144`, `Q_D=8`, `max_r d_D(r)=3`.

- Full histogram: `0:1860, 1:131, 2:5, 3:1`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=0, q100=3`.
- Top ten (ties ordered by increasing `r`):

  1. r=979 (m=3, r/p=0.4897, Z_p=no, small=no, gaps=[17, 23, 40], exact-mirror gaps=[40], min |2r+d-(p-1)|=0)
  2. r=548 (m=2, r/p=0.2741, Z_p=no, small=no, gaps=[17, 92], exact-mirror gaps=none, min |2r+d-(p-1)|=810)
  3. r=920 (m=2, r/p=0.4602, Z_p=no, small=no, gaps=[65, 93], exact-mirror gaps=none, min |2r+d-(p-1)|=65)
  4. r=985 (m=2, r/p=0.4927, Z_p=no, small=no, gaps=[28, 93], exact-mirror gaps=[28], min |2r+d-(p-1)|=0)
  5. r=996 (m=2, r/p=0.4982, Z_p=no, small=no, gaps=[6, 23], exact-mirror gaps=[6], min |2r+d-(p-1)|=0)
  6. r=1358 (m=2, r/p=0.6793, Z_p=no, small=no, gaps=[75, 92], exact-mirror gaps=none, min |2r+d-(p-1)|=793)
  7. r=23 (m=1, r/p=0.0115, Z_p=no, small=yes, gaps=[13], exact-mirror gaps=none, min |2r+d-(p-1)|=1939)
  8. r=30 (m=1, r/p=0.0150, Z_p=no, small=yes, gaps=[31], exact-mirror gaps=none, min |2r+d-(p-1)|=1907)
  9. r=42 (m=1, r/p=0.0210, Z_p=no, small=yes, gaps=[90], exact-mirror gaps=none, min |2r+d-(p-1)|=1824)
  10. r=53 (m=1, r/p=0.0265, Z_p=no, small=yes, gaps=[52], exact-mirror gaps=none, min |2r+d-(p-1)|=1840)

#### ceil(p^0.66): D=151

`S_D=213`, `Q_D=16`, `max_r d_D(r)=3`.

- Full histogram: `0:1798, 1:187, 2:10, 3:2`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=0, q100=3`.
- Top ten (ties ordered by increasing `r`):

  1. r=927 (m=3, r/p=0.4637, Z_p=no, small=no, gaps=[17, 127, 144], exact-mirror gaps=[144], min |2r+d-(p-1)|=0)
  2. r=979 (m=3, r/p=0.4897, Z_p=no, small=no, gaps=[17, 23, 40], exact-mirror gaps=[40], min |2r+d-(p-1)|=0)
  3. r=488 (m=2, r/p=0.2441, Z_p=no, small=no, gaps=[10, 127], exact-mirror gaps=none, min |2r+d-(p-1)|=895)
  4. r=548 (m=2, r/p=0.2741, Z_p=no, small=no, gaps=[17, 92], exact-mirror gaps=none, min |2r+d-(p-1)|=810)
  5. r=884 (m=2, r/p=0.4422, Z_p=no, small=no, gaps=[113, 117], exact-mirror gaps=none, min |2r+d-(p-1)|=113)
  6. r=920 (m=2, r/p=0.4602, Z_p=no, small=no, gaps=[65, 93], exact-mirror gaps=none, min |2r+d-(p-1)|=65)
  7. r=944 (m=2, r/p=0.4722, Z_p=no, small=no, gaps=[110, 127], exact-mirror gaps=[110], min |2r+d-(p-1)|=0)
  8. r=985 (m=2, r/p=0.4927, Z_p=no, small=no, gaps=[28, 93], exact-mirror gaps=[28], min |2r+d-(p-1)|=0)
  9. r=996 (m=2, r/p=0.4982, Z_p=no, small=no, gaps=[6, 23], exact-mirror gaps=[6], min |2r+d-(p-1)|=0)
  10. r=997 (m=2, r/p=0.4987, Z_p=no, small=no, gaps=[4, 117], exact-mirror gaps=[4], min |2r+d-(p-1)|=0)

### p=4001

`N=3999`, `|Z_p intersect [1,N]|=2`, Delta-identity checks `1961400` (all passed).

#### ceil(sqrt(p) log(p)): D=525

`S_D=754`, `Q_D=78`, `max_r d_D(r)=5`.

- Full histogram: `0:3304, 1:651, 2:32, 3:10, 4:1, 5:1`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=1, q100=5`.
- Top ten (ties ordered by increasing `r`):

  1. r=1792 (m=5, r/p=0.4479, Z_p=no, small=no, gaps=[3, 58, 358, 413, 416], exact-mirror gaps=[416], min |2r+d-(p-1)|=0)
  2. r=1795 (m=4, r/p=0.4486, Z_p=no, small=no, gaps=[55, 355, 410, 413], exact-mirror gaps=[410], min |2r+d-(p-1)|=0)
  3. r=1631 (m=3, r/p=0.4076, Z_p=no, small=no, gaps=[175, 304, 434], exact-mirror gaps=none, min |2r+d-(p-1)|=304)
  4. r=1705 (m=3, r/p=0.4261, Z_p=no, small=no, gaps=[62, 222, 368], exact-mirror gaps=none, min |2r+d-(p-1)|=222)
  5. r=1743 (m=3, r/p=0.4356, Z_p=no, small=no, gaps=[75, 439, 514], exact-mirror gaps=[514], min |2r+d-(p-1)|=0)
  6. r=1767 (m=3, r/p=0.4416, Z_p=no, small=no, gaps=[160, 306, 466], exact-mirror gaps=[466], min |2r+d-(p-1)|=0)
  7. r=1802 (m=3, r/p=0.4504, Z_p=no, small=no, gaps=[52, 344, 396], exact-mirror gaps=[396], min |2r+d-(p-1)|=0)
  8. r=1806 (m=3, r/p=0.4514, Z_p=no, small=no, gaps=[129, 259, 388], exact-mirror gaps=[388], min |2r+d-(p-1)|=0)
  9. r=1850 (m=3, r/p=0.4624, Z_p=no, small=no, gaps=[300, 355, 358], exact-mirror gaps=[300], min |2r+d-(p-1)|=0)
  10. r=1876 (m=3, r/p=0.4689, Z_p=no, small=no, gaps=[70, 178, 248], exact-mirror gaps=[248], min |2r+d-(p-1)|=0)

#### ceil(p^0.6): D=145

`S_D=198`, `Q_D=2`, `max_r d_D(r)=2`.

- Full histogram: `0:3803, 1:194, 2:2`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=0, q100=2`.
- Top ten (ties ordered by increasing `r`):

  1. r=1792 (m=2, r/p=0.4479, Z_p=no, small=no, gaps=[3, 58], exact-mirror gaps=none, min |2r+d-(p-1)|=358)
  2. r=2150 (m=2, r/p=0.5374, Z_p=no, small=no, gaps=[55, 58], exact-mirror gaps=none, min |2r+d-(p-1)|=355)
  3. r=9 (m=1, r/p=0.0022, Z_p=no, small=yes, gaps=[90], exact-mirror gaps=none, min |2r+d-(p-1)|=3892)
  4. r=12 (m=1, r/p=0.0030, Z_p=no, small=yes, gaps=[7], exact-mirror gaps=none, min |2r+d-(p-1)|=3969)
  5. r=38 (m=1, r/p=0.0095, Z_p=no, small=yes, gaps=[93], exact-mirror gaps=none, min |2r+d-(p-1)|=3831)
  6. r=62 (m=1, r/p=0.0155, Z_p=no, small=yes, gaps=[24], exact-mirror gaps=none, min |2r+d-(p-1)|=3852)
  7. r=79 (m=1, r/p=0.0197, Z_p=no, small=no, gaps=[58], exact-mirror gaps=none, min |2r+d-(p-1)|=3784)
  8. r=154 (m=1, r/p=0.0385, Z_p=no, small=no, gaps=[88], exact-mirror gaps=none, min |2r+d-(p-1)|=3604)
  9. r=179 (m=1, r/p=0.0447, Z_p=no, small=no, gaps=[79], exact-mirror gaps=none, min |2r+d-(p-1)|=3563)
  10. r=182 (m=1, r/p=0.0455, Z_p=no, small=no, gaps=[133], exact-mirror gaps=none, min |2r+d-(p-1)|=3503)

#### ceil(p^0.66): D=239

`S_D=345`, `Q_D=6`, `max_r d_D(r)=2`.

- Full histogram: `0:3660, 1:333, 2:6`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=0, q100=2`.
- Top ten (ties ordered by increasing `r`):

  1. r=1705 (m=2, r/p=0.4261, Z_p=no, small=no, gaps=[62, 222], exact-mirror gaps=none, min |2r+d-(p-1)|=368)
  2. r=1792 (m=2, r/p=0.4479, Z_p=no, small=no, gaps=[3, 58], exact-mirror gaps=none, min |2r+d-(p-1)|=358)
  3. r=1876 (m=2, r/p=0.4689, Z_p=no, small=no, gaps=[70, 178], exact-mirror gaps=none, min |2r+d-(p-1)|=70)
  4. r=1946 (m=2, r/p=0.4864, Z_p=no, small=no, gaps=[108, 178], exact-mirror gaps=[108], min |2r+d-(p-1)|=0)
  5. r=2073 (m=2, r/p=0.5181, Z_p=no, small=no, gaps=[160, 222], exact-mirror gaps=none, min |2r+d-(p-1)|=306)
  6. r=2150 (m=2, r/p=0.5374, Z_p=no, small=no, gaps=[55, 58], exact-mirror gaps=none, min |2r+d-(p-1)|=355)
  7. r=9 (m=1, r/p=0.0022, Z_p=no, small=yes, gaps=[90], exact-mirror gaps=none, min |2r+d-(p-1)|=3892)
  8. r=12 (m=1, r/p=0.0030, Z_p=no, small=yes, gaps=[7], exact-mirror gaps=none, min |2r+d-(p-1)|=3969)
  9. r=38 (m=1, r/p=0.0095, Z_p=no, small=yes, gaps=[93], exact-mirror gaps=none, min |2r+d-(p-1)|=3831)
  10. r=62 (m=1, r/p=0.0155, Z_p=no, small=yes, gaps=[24], exact-mirror gaps=none, min |2r+d-(p-1)|=3852)

### p=7919

`N=7917`, `|Z_p intersect [1,N]|=0`, Delta-identity checks `6006083` (all passed).

#### ceil(sqrt(p) log(p)): D=799

`S_D=1229`, `Q_D=142`, `max_r d_D(r)=5`.

- Full histogram: `0:6803, 1:1022, 2:72, 3:18, 4:1, 5:1`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=1, q100=5`.
- Top ten (ties ordered by increasing `r`):

  1. r=3879 (m=5, r/p=0.4898, Z_p=no, small=no, gaps=[66, 76, 84, 94, 160], exact-mirror gaps=[160], min |2r+d-(p-1)|=0)
  2. r=3945 (m=4, r/p=0.4982, Z_p=no, small=no, gaps=[10, 18, 28, 94], exact-mirror gaps=[28], min |2r+d-(p-1)|=0)
  3. r=3293 (m=3, r/p=0.4158, Z_p=no, small=no, gaps=[363, 539, 793], exact-mirror gaps=none, min |2r+d-(p-1)|=539)
  4. r=3571 (m=3, r/p=0.4509, Z_p=no, small=no, gaps=[291, 485, 776], exact-mirror gaps=[776], min |2r+d-(p-1)|=0)
  5. r=3586 (m=3, r/p=0.4528, Z_p=no, small=no, gaps=[113, 633, 746], exact-mirror gaps=[746], min |2r+d-(p-1)|=0)
  6. r=3593 (m=3, r/p=0.4537, Z_p=no, small=no, gaps=[363, 369, 732], exact-mirror gaps=[732], min |2r+d-(p-1)|=0)
  7. r=3611 (m=3, r/p=0.4560, Z_p=no, small=no, gaps=[161, 535, 696], exact-mirror gaps=[696], min |2r+d-(p-1)|=0)
  8. r=3653 (m=3, r/p=0.4613, Z_p=no, small=no, gaps=[31, 581, 612], exact-mirror gaps=[612], min |2r+d-(p-1)|=0)
  9. r=3656 (m=3, r/p=0.4617, Z_p=no, small=no, gaps=[176, 430, 606], exact-mirror gaps=[606], min |2r+d-(p-1)|=0)
  10. r=3678 (m=3, r/p=0.4645, Z_p=no, small=no, gaps=[250, 312, 562], exact-mirror gaps=[562], min |2r+d-(p-1)|=0)

#### ceil(p^0.6): D=219

`S_D=349`, `Q_D=32`, `max_r d_D(r)=5`.

- Full histogram: `0:7588, 1:317, 2:7, 3:3, 4:1, 5:1`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=0, q100=5`.
- Top ten (ties ordered by increasing `r`):

  1. r=3879 (m=5, r/p=0.4898, Z_p=no, small=no, gaps=[66, 76, 84, 94, 160], exact-mirror gaps=[160], min |2r+d-(p-1)|=0)
  2. r=3945 (m=4, r/p=0.4982, Z_p=no, small=no, gaps=[10, 18, 28, 94], exact-mirror gaps=[28], min |2r+d-(p-1)|=0)
  3. r=3860 (m=3, r/p=0.4874, Z_p=no, small=no, gaps=[24, 174, 198], exact-mirror gaps=[198], min |2r+d-(p-1)|=0)
  4. r=3893 (m=3, r/p=0.4916, Z_p=no, small=no, gaps=[17, 115, 132], exact-mirror gaps=[132], min |2r+d-(p-1)|=0)
  5. r=3955 (m=3, r/p=0.4994, Z_p=no, small=no, gaps=[8, 18, 84], exact-mirror gaps=[8], min |2r+d-(p-1)|=0)
  6. r=962 (m=2, r/p=0.1215, Z_p=no, small=no, gaps=[95, 148], exact-mirror gaps=none, min |2r+d-(p-1)|=5846)
  7. r=2675 (m=2, r/p=0.3378, Z_p=no, small=no, gaps=[19, 196], exact-mirror gaps=none, min |2r+d-(p-1)|=2372)
  8. r=3884 (m=2, r/p=0.4905, Z_p=no, small=no, gaps=[150, 174], exact-mirror gaps=[150], min |2r+d-(p-1)|=0)
  9. r=3910 (m=2, r/p=0.4937, Z_p=no, small=no, gaps=[98, 115], exact-mirror gaps=[98], min |2r+d-(p-1)|=0)
  10. r=3963 (m=2, r/p=0.5004, Z_p=no, small=no, gaps=[10, 76], exact-mirror gaps=none, min |2r+d-(p-1)|=18)

#### ceil(p^0.66): D=375

`S_D=575`, `Q_D=50`, `max_r d_D(r)=5`.

- Full histogram: `0:7379, 1:510, 2:22, 3:4, 4:1, 5:1`.
- Deciles: `q0=0, q10=0, q20=0, q30=0, q40=0, q50=0, q60=0, q70=0, q80=0, q90=0, q100=5`.
- Top ten (ties ordered by increasing `r`):

  1. r=3879 (m=5, r/p=0.4898, Z_p=no, small=no, gaps=[66, 76, 84, 94, 160], exact-mirror gaps=[160], min |2r+d-(p-1)|=0)
  2. r=3945 (m=4, r/p=0.4982, Z_p=no, small=no, gaps=[10, 18, 28, 94], exact-mirror gaps=[28], min |2r+d-(p-1)|=0)
  3. r=3802 (m=3, r/p=0.4801, Z_p=no, small=no, gaps=[3, 311, 314], exact-mirror gaps=[314], min |2r+d-(p-1)|=0)
  4. r=3860 (m=3, r/p=0.4874, Z_p=no, small=no, gaps=[24, 174, 198], exact-mirror gaps=[198], min |2r+d-(p-1)|=0)
  5. r=3893 (m=3, r/p=0.4916, Z_p=no, small=no, gaps=[17, 115, 132], exact-mirror gaps=[132], min |2r+d-(p-1)|=0)
  6. r=3955 (m=3, r/p=0.4994, Z_p=no, small=no, gaps=[8, 18, 84], exact-mirror gaps=[8], min |2r+d-(p-1)|=0)
  7. r=962 (m=2, r/p=0.1215, Z_p=no, small=no, gaps=[95, 148], exact-mirror gaps=none, min |2r+d-(p-1)|=5846)
  8. r=2675 (m=2, r/p=0.3378, Z_p=no, small=no, gaps=[19, 196], exact-mirror gaps=none, min |2r+d-(p-1)|=2372)
  9. r=3593 (m=2, r/p=0.4537, Z_p=no, small=no, gaps=[363, 369], exact-mirror gaps=none, min |2r+d-(p-1)|=363)
  10. r=3678 (m=2, r/p=0.4645, Z_p=no, small=no, gaps=[250, 312], exact-mirror gaps=none, min |2r+d-(p-1)|=250)

## 2. Vector-7 premise: maximum return multiplicity

The comparison columns are `max/D^(2/3)` and `max/log(p)`. The final column lists every maximizing base, not merely the first ten. `mirror` records whether at least one of that base's collisions is the exact forced mirror collision `2r+d=p-1`; `near2` is the minimum doubled distance `|2r+d-(p-1)|`. `small` means `r<=ceil((log p)^2)`.

Each cell is written as a code block, followed by the complete maximizing-base list.

`p=997, D=219, max=4, D^(2/3)=36.332566, log(p)=6.904751, max/D^(2/3)=0.110094, max/log(p)=0.579311`

Maximizing bases: 347(Z=N,small=N,mirror=N,near2=100,gaps=[100, 126, 176, 202]); 447(Z=N,small=N,mirror=Y,near2=0,gaps=[26, 76, 102, 202]).

`p=997, D=63, max=1, D^(2/3)=15.832896, log(p)=6.904751, max/D^(2/3)=0.063160, max/log(p)=0.144828`

Maximizing bases: 7(Z=N,small=Y,mirror=N,near2=920,gaps=[62]); 18(Z=N,small=Y,mirror=N,near2=916,gaps=[44]); 69(Z=N,small=N,mirror=N,near2=847,gaps=[11]); 79(Z=N,small=N,mirror=N,near2=820,gaps=[18]); 86(Z=N,small=N,mirror=N,near2=804,gaps=[20]); 97(Z=N,small=N,mirror=N,near2=740,gaps=[62]); 115(Z=N,small=N,mirror=N,near2=714,gaps=[52]); 144(Z=N,small=N,mirror=N,near2=669,gaps=[39]); 180(Z=N,small=N,mirror=N,near2=603,gaps=[33]); 198(Z=N,small=N,mirror=N,near2=544,gaps=[56]); 206(Z=N,small=N,mirror=N,near2=566,gaps=[18]); 211(Z=N,small=N,mirror=N,near2=567,gaps=[7]); 242(Z=N,small=N,mirror=N,near2=469,gaps=[43]); 246(Z=N,small=N,mirror=N,near2=458,gaps=[46]); 249(Z=N,small=N,mirror=N,near2=471,gaps=[27]); 253(Z=N,small=N,mirror=N,near2=430,gaps=[60]); 265(Z=N,small=N,mirror=N,near2=436,gaps=[30]); 269(Z=N,small=N,mirror=N,near2=421,gaps=[37]); 274(Z=N,small=N,mirror=N,near2=425,gaps=[23]); 275(Z=N,small=N,mirror=N,near2=383,gaps=[63]); 325(Z=N,small=N,mirror=N,near2=328,gaps=[18]); 346(Z=N,small=N,mirror=N,near2=248,gaps=[56]); 376(Z=N,small=N,mirror=N,near2=227,gaps=[17]); 380(Z=N,small=N,mirror=N,near2=174,gaps=[62]); 394(Z=N,small=N,mirror=N,near2=147,gaps=[61]); 406(Z=N,small=N,mirror=N,near2=174,gaps=[10]); 445(Z=N,small=N,mirror=N,near2=84,gaps=[22]); 447(Z=N,small=N,mirror=N,near2=76,gaps=[26]); 467(Z=N,small=N,mirror=Y,near2=0,gaps=[62]); 468(Z=N,small=N,mirror=Y,near2=0,gaps=[60]); 469(Z=N,small=N,mirror=Y,near2=0,gaps=[58]); 470(Z=N,small=N,mirror=Y,near2=0,gaps=[56]); 471(Z=N,small=N,mirror=Y,near2=0,gaps=[54]); 472(Z=N,small=N,mirror=Y,near2=0,gaps=[52]); 473(Z=N,small=N,mirror=Y,near2=0,gaps=[50]); 474(Z=N,small=N,mirror=Y,near2=0,gaps=[48]); 475(Z=N,small=N,mirror=Y,near2=0,gaps=[46]); 476(Z=N,small=N,mirror=Y,near2=0,gaps=[44]); 477(Z=N,small=N,mirror=Y,near2=0,gaps=[42]); 478(Z=N,small=N,mirror=Y,near2=0,gaps=[40]); 479(Z=N,small=N,mirror=Y,near2=0,gaps=[38]); 480(Z=N,small=N,mirror=Y,near2=0,gaps=[36]); 481(Z=N,small=N,mirror=Y,near2=0,gaps=[34]); 482(Z=N,small=N,mirror=Y,near2=0,gaps=[32]); 483(Z=N,small=N,mirror=Y,near2=0,gaps=[30]); 484(Z=N,small=N,mirror=Y,near2=0,gaps=[28]); 485(Z=N,small=N,mirror=Y,near2=0,gaps=[26]); 486(Z=N,small=N,mirror=Y,near2=0,gaps=[24]); 487(Z=N,small=N,mirror=Y,near2=0,gaps=[22]); 488(Z=N,small=N,mirror=Y,near2=0,gaps=[20]); 489(Z=N,small=N,mirror=Y,near2=0,gaps=[18]); 490(Z=N,small=N,mirror=Y,near2=0,gaps=[16]); 491(Z=N,small=N,mirror=Y,near2=0,gaps=[14]); 492(Z=N,small=N,mirror=Y,near2=0,gaps=[12]); 493(Z=N,small=N,mirror=Y,near2=0,gaps=[10]); 494(Z=N,small=N,mirror=Y,near2=0,gaps=[8]); 495(Z=N,small=N,mirror=Y,near2=0,gaps=[6]); 496(Z=N,small=N,mirror=Y,near2=0,gaps=[4]); 497(Z=N,small=N,mirror=Y,near2=0,gaps=[2]); 523(Z=N,small=N,mirror=N,near2=76,gaps=[26]); 529(Z=N,small=N,mirror=N,near2=84,gaps=[22]); 541(Z=N,small=N,mirror=N,near2=147,gaps=[61]); 554(Z=N,small=N,mirror=N,near2=174,gaps=[62]); 580(Z=N,small=N,mirror=N,near2=174,gaps=[10]); 594(Z=N,small=N,mirror=N,near2=248,gaps=[56]); 603(Z=N,small=N,mirror=N,near2=227,gaps=[17]); 653(Z=N,small=N,mirror=N,near2=328,gaps=[18]); 658(Z=N,small=N,mirror=N,near2=383,gaps=[63]); 683(Z=N,small=N,mirror=N,near2=430,gaps=[60]); 690(Z=N,small=N,mirror=N,near2=421,gaps=[37]); 699(Z=N,small=N,mirror=N,near2=425,gaps=[23]); 701(Z=N,small=N,mirror=N,near2=436,gaps=[30]); 704(Z=N,small=N,mirror=N,near2=458,gaps=[46]); 711(Z=N,small=N,mirror=N,near2=469,gaps=[43]); 720(Z=N,small=N,mirror=N,near2=471,gaps=[27]); 742(Z=N,small=N,mirror=N,near2=544,gaps=[56]); 772(Z=N,small=N,mirror=N,near2=566,gaps=[18]); 778(Z=N,small=N,mirror=N,near2=567,gaps=[7]); 783(Z=N,small=N,mirror=N,near2=603,gaps=[33]); 813(Z=N,small=N,mirror=N,near2=669,gaps=[39]); 829(Z=N,small=N,mirror=N,near2=714,gaps=[52]); 837(Z=N,small=N,mirror=N,near2=740,gaps=[62]); 890(Z=N,small=N,mirror=N,near2=804,gaps=[20]); 899(Z=N,small=N,mirror=N,near2=820,gaps=[18]); 916(Z=N,small=N,mirror=N,near2=847,gaps=[11]); 927(Z=N,small=N,mirror=N,near2=920,gaps=[62]); 934(Z=N,small=N,mirror=N,near2=916,gaps=[44]).

`p=997, D=96, max=2, D^(2/3)=20.965931, log(p)=6.904751, max/D^(2/3)=0.095393, max/log(p)=0.289656`

Maximizing bases: 7(Z=N,small=Y,mirror=N,near2=909,gaps=[62, 73]); 79(Z=N,small=N,mirror=N,near2=758,gaps=[18, 80]); 445(Z=N,small=N,mirror=N,near2=22,gaps=[22, 84]); 447(Z=N,small=N,mirror=N,near2=26,gaps=[26, 76]); 467(Z=N,small=N,mirror=Y,near2=0,gaps=[62, 84]); 473(Z=N,small=N,mirror=Y,near2=0,gaps=[50, 76]); 837(Z=N,small=N,mirror=N,near2=740,gaps=[62, 80]); 916(Z=N,small=N,mirror=N,near2=847,gaps=[11, 73]).

`p=1999, D=340, max=4, D^(2/3)=48.713868, log(p)=7.600402, max/D^(2/3)=0.082112, max/log(p)=0.526288`

Maximizing bases: 810(Z=N,small=N,mirror=N,near2=169,gaps=[169, 186, 192, 209]); 979(Z=N,small=N,mirror=Y,near2=0,gaps=[17, 23, 40, 209]).

`p=1999, D=96, max=3, D^(2/3)=20.965931, log(p)=7.600402, max/D^(2/3)=0.143089, max/log(p)=0.394716`

Maximizing bases: 979(Z=N,small=N,mirror=Y,near2=0,gaps=[17, 23, 40]).

`p=1999, D=151, max=3, D^(2/3)=28.356413, log(p)=7.600402, max/D^(2/3)=0.105796, max/log(p)=0.394716`

Maximizing bases: 927(Z=N,small=N,mirror=Y,near2=0,gaps=[17, 127, 144]); 979(Z=N,small=N,mirror=Y,near2=0,gaps=[17, 23, 40]).

`p=4001, D=525, max=5, D^(2/3)=65.078800, log(p)=8.294300, max/D^(2/3)=0.076830, max/log(p)=0.602824`

Maximizing bases: 1792(Z=N,small=N,mirror=Y,near2=0,gaps=[3, 58, 358, 413, 416]).

`p=4001, D=145, max=2, D^(2/3)=27.600186, log(p)=8.294300, max/D^(2/3)=0.072463, max/log(p)=0.241129`

Maximizing bases: 1792(Z=N,small=N,mirror=N,near2=358,gaps=[3, 58]); 2150(Z=N,small=N,mirror=N,near2=355,gaps=[55, 58]).

`p=4001, D=239, max=2, D^(2/3)=38.512224, log(p)=8.294300, max/D^(2/3)=0.051932, max/log(p)=0.241129`

Maximizing bases: 1705(Z=N,small=N,mirror=N,near2=368,gaps=[62, 222]); 1792(Z=N,small=N,mirror=N,near2=358,gaps=[3, 58]); 1876(Z=N,small=N,mirror=N,near2=70,gaps=[70, 178]); 1946(Z=N,small=N,mirror=Y,near2=0,gaps=[108, 178]); 2073(Z=N,small=N,mirror=N,near2=306,gaps=[160, 222]); 2150(Z=N,small=N,mirror=N,near2=355,gaps=[55, 58]).

`p=7919, D=799, max=5, D^(2/3)=86.105558, log(p)=8.977020, max/D^(2/3)=0.058068, max/log(p)=0.556978`

Maximizing bases: 3879(Z=N,small=N,mirror=Y,near2=0,gaps=[66, 76, 84, 94, 160]).

`p=7919, D=219, max=5, D^(2/3)=36.332566, log(p)=8.977020, max/D^(2/3)=0.137618, max/log(p)=0.556978`

Maximizing bases: 3879(Z=N,small=N,mirror=Y,near2=0,gaps=[66, 76, 84, 94, 160]).

`p=7919, D=375, max=5, D^(2/3)=52.002096, log(p)=8.977020, max/D^(2/3)=0.096150, max/log(p)=0.556978`

Maximizing bases: 3879(Z=N,small=N,mirror=Y,near2=0,gaps=[66, 76, 84, 94, 160]).


Across all 12 cells the maximum lies in `1..5` while `D^(2/3)` grows by a much larger factor. These data strongly reject `c D^(2/3)` at the tested scale. They are compatible with an absolute bound and also with logarithmic growth too slow to resolve from four primes; the bounded description is the sharper empirical fit, not a proof.

## 3. Vector-8 premise: primitive decomposition

For an endpoint collision `(r,d)`, its split multiplicity is the number of `d'` with `0<d'<d` for which both `(r,d')` and `(r+d',d-d')` collide. Because the determinant and projective-fiber computations agreed pair by pair, this is also the exact number of intermediate occurrences of the same projective value. Histogram notation is `split multiplicity:number of endpoint collisions`. In every cell, `P_D` is the zero bin, every positive-bin collision splits, and the exact renewal checksum is `sum_k k*n_k=Q_D`.

The per-cell records below are `(p,D,P_D,S_D,P_D/S_D,nonprimitive,split_histogram,split_witnesses,Q_D,failures)`.

- `(997,219,266,313,0.849840,47,[0:266, 1:38, 2:7, 3:2],58,58,0)`
- `(997,63,87,87,1.000000,0,[0:87],0,0,0)`
- `(997,96,128,136,0.941176,8,[0:128, 1:8],8,8,0)`
- `(1999,340,431,498,0.865462,67,[0:431, 1:52, 2:13, 3:2],84,84,0)`
- `(1999,96,137,144,0.951389,7,[0:137, 1:6, 2:1],8,8,0)`
- `(1999,151,199,213,0.934272,14,[0:199, 1:12, 2:2],16,16,0)`
- `(4001,525,695,754,0.921751,59,[0:695, 1:44, 2:12, 3:2, 4:1],78,78,0)`
- `(4001,145,196,198,0.989899,2,[0:196, 1:2],2,2,0)`
- `(4001,239,339,345,0.982609,6,[0:339, 1:6],6,6,0)`
- `(7919,799,1114,1229,0.906428,115,[0:1114, 1:92, 2:20, 3:2, 4:1],142,142,0)`
- `(7919,219,329,349,0.942693,20,[0:329, 1:12, 2:5, 3:2, 4:1],32,32,0)`
- `(7919,375,538,575,0.935652,37,[0:538, 1:28, 2:6, 3:2, 4:1],50,50,0)`

The renewal claim therefore passes every endpoint collision in all 12 cells, with no unsplit nonprimitive collision. The primitive share remains high, so the data validate decomposition but do not by themselves supply a small primitive bound.

## 4. Exact small-d full root counts

Each line gives every exact `d:R_d` value through `K=ceil((log p)^2)`. The mean is also shown as the exact fraction `sum R_d/K`.

### p=997, K=48

`max R_d=7` at `d=[18]`; mean `=62/48=1.291667`.

`1:0, 2:1, 3:0, 4:1, 5:0, 6:1, 7:2, 8:1, 9:0, 10:3, 11:2, 12:1, 13:0, 14:1, 15:0, 16:1, 17:2, 18:7, 19:0, 20:3, 21:0, 22:3, 23:2, 24:1, 25:0, 26:3, 27:2, 28:1, 29:0, 30:3, 31:0, 32:1, 33:2, 34:1, 35:0, 36:1, 37:2, 38:1, 39:2, 40:1, 41:0, 42:1, 43:2, 44:3, 45:0, 46:3, 47:0, 48:1`

### p=1999, K=58

`max R_d=6` at `d=[17]`; mean `=91/58=1.568966`.

`1:0, 2:1, 3:0, 4:1, 5:4, 6:1, 7:0, 8:3, 9:0, 10:5, 11:0, 12:1, 13:4, 14:1, 15:4, 16:3, 17:6, 18:1, 19:4, 20:1, 21:2, 22:1, 23:2, 24:1, 25:2, 26:1, 27:0, 28:5, 29:2, 30:1, 31:4, 32:1, 33:0, 34:3, 35:0, 36:1, 37:0, 38:1, 39:0, 40:1, 41:0, 42:3, 43:0, 44:1, 45:0, 46:3, 47:2, 48:1, 49:0, 50:1, 51:0, 52:5, 53:0, 54:1, 55:0, 56:3, 57:0, 58:3`

### p=4001, K=69

`max R_d=8` at `d=[7]`; mean `=92/69=1.333333`.

`1:0, 2:3, 3:2, 4:1, 5:0, 6:1, 7:8, 8:1, 9:2, 10:1, 11:0, 12:1, 13:0, 14:1, 15:4, 16:1, 17:0, 18:1, 19:0, 20:5, 21:2, 22:1, 23:0, 24:3, 25:2, 26:1, 27:2, 28:1, 29:0, 30:3, 31:0, 32:1, 33:0, 34:5, 35:0, 36:1, 37:0, 38:1, 39:0, 40:3, 41:0, 42:1, 43:0, 44:1, 45:2, 46:1, 47:2, 48:1, 49:2, 50:3, 51:0, 52:3, 53:0, 54:1, 55:2, 56:1, 57:2, 58:5, 59:0, 60:1, 61:0, 62:3, 63:0, 64:1, 65:0, 66:1, 67:0, 68:1, 69:0`

### p=7919, K=81

`max R_d=7` at `d=[6]`; mean `=134/81=1.654321`.

`1:0, 2:3, 3:4, 4:1, 5:0, 6:7, 7:0, 8:1, 9:2, 10:3, 11:0, 12:3, 13:2, 14:1, 15:4, 16:5, 17:4, 18:3, 19:4, 20:1, 21:0, 22:1, 23:0, 24:3, 25:0, 26:1, 27:0, 28:3, 29:0, 30:3, 31:2, 32:1, 33:4, 34:3, 35:2, 36:3, 37:0, 38:1, 39:0, 40:3, 41:2, 42:1, 43:0, 44:3, 45:0, 46:1, 47:2, 48:1, 49:0, 50:1, 51:0, 52:1, 53:6, 54:1, 55:0, 56:1, 57:0, 58:1, 59:0, 60:1, 61:0, 62:1, 63:2, 64:5, 65:6, 66:5, 67:0, 68:1, 69:2, 70:1, 71:0, 72:1, 73:0, 74:3, 75:0, 76:3, 77:2, 78:1, 79:0, 80:1, 81:0`

The largest observed small-segment root count is `8`. Means stay of constant size. Thus `R_d=O(1)` on this polylogarithmic segment is empirically available for these four primes (with observed constant `8`), but the finite scan is not a uniform theorem in `p`.

## 5. Verdicts for the three deep-strike premises

- **Vector 7 -- HOLDS EMPIRICALLY IN THE BOUNDED/LOG-SIZED FORM, NOT THE `D^(2/3)` FORM.** The maxima remain tiny compared with both `D^(2/3)` and the available window length. Maximizers are not systematically Apery zeros, small bases, or exact mirror centers; the per-cell annotations give the full exceptions. This supports anti-concentration, but does not prove the uniform bound needed by the strike.

- **Vector 8 -- RENEWAL CLAIM VERIFIED EXACTLY; THE NEEDED PRIMITIVE BOUND REMAINS OPEN.** Every nonprimitive endpoint splits and the weighted split histogram equals `Q_D` in every cell. However, most endpoint collisions are primitive at these scales, so decomposition alone has not reduced `S_D` to a demonstrably smaller quantity.

- **Small-d input -- HOLDS EMPIRICALLY.** Complete root counts through `ceil((log p)^2)` have constant-sized means and a small global maximum. This makes the proposed `R_d=O(1)` input numerically plausible on the entire polylog segment, while supplying no proof beyond the four tested primes.
