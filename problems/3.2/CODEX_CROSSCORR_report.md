# Exact cross-gap correlation experiment

This report executes `CODEX_SPEC_crosscorr.md` on the exact windows I_h={1,...,M-h}, M=p-2. All delta values and coincidence counts use integer arithmetic modulo p. Every E_(h,k) is retained exactly as (p*count-|I_h||I_k|)/p; decimal output is only presentation.

The scales are integer-exact ceilings: H^5>=p^2 for exponent 0.4, H^2>=p for exponent 0.5, and H^5>=p^3 for exponent 0.6.

For same gaps, E_(h,h) has its original definition, including r=r'. The raw r!=r' coincidence count is also reported separately. Thus E_(h,h)=P_off(h)+|I_h|-|I_h|^2/p, which is the normalization in which the off-diagonal H_h component contributes about +p.

## p=997 (p mod 3 = 1)

H values: 16, 32, 63.

### H=16

All 120 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 4499372/997 = 4512.910732.
- sum_(h<k) |E_(h,k)|/(p H) = 0.282906.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.558302.
- max_d |A_d|/p = 0.276895, attained at d=1.
- V(H) = 29677093487/1976107784162 = 0.015018.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 1: 275236/997; +0.276895; 7: -192366/997; -0.193525; 12: -175774/997; -0.176833; 3: -143160/997; -0.144023; 8: -142370/997; -0.143228.
- Signed A_d/p distribution deciles: q0=-0.1935 q10=-0.1637 q20=-0.1434 q30=-0.1084 q40=-0.0501 q50=-0.0302 q60=+0.0025 q70=+0.0593 q80=+0.0786 q90=+0.1025 q100=+0.2769.
- Ordered-d profile samples (d: A_d/p): 2:-0.0790; 3:-0.1440; 5:+0.0649; 6:-0.0302; 8:-0.1432; 9:-0.0308; 11:-0.0205; 12:-0.1768; 14:+0.0370; 15:+0.1063.
- Mean |E_(h,k)|/p diagnostics: even d=0.03836, odd d=0.03716, d<=10=0.03675, d>10=0.04448, k=2h=0.03259, k!=2h=0.03809.
- h=1 cross row: mean |E_(1,k)|/p=0.05650; maximum |E_(1,k)|/p=0.12339 at k=5.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.8867 q10=+0.9168 q20=+0.9228 q30=+0.9358 q40=+0.9448 q50=+0.9729 q60=+0.9950 q70=+0.9980 q80=+1.0070 q90=+1.0371 q100=+1.9900.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-1.1562 q10=-0.7406 q20=-0.3857 q30=-0.2302 q40=-0.1567 q50=-0.0596 q60=+0.0247 q70=+0.0308 q80=+0.0629 q90=+0.1291 q100=+31.3533.
- The same normalized deciles restricted to generic h>=2: q0=-1.1562 q10=-0.7619 q20=-0.4353 q30=-0.2389 q40=-0.1803 q50=-0.1214 q60=+0.0112 q70=+0.0274 q80=+0.0394 q90=+0.0938 q100=+0.1439.
- h=1 self row: off-diagonal count=1984 (1.989970 p); E_(1,1)/p=+1.992970.

### H=32

All 496 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 17847944/997 = 17901.648947.
- sum_(h<k) |E_(h,k)|/(p H) = 0.561110.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.553663.
- max_d |A_d|/p = 0.407147, attained at d=3.
- V(H) = 52355204537/1976107784162 = 0.026494.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 3: -404708/997; -0.407147; 6: -356086/997; -0.358232; 15: 335704/997; +0.337727; 14: 262430/997; +0.264012; 16: 260410/997; +0.261980.
- Signed A_d/p distribution deciles: q0=-0.4071 q10=-0.1618 q20=-0.0745 q30=-0.0428 q40=+0.0039 q50=+0.0091 q60=+0.0247 q70=+0.0855 q80=+0.1597 q90=+0.2579 q100=+0.3377.
- Ordered-d profile samples (d: A_d/p): 4:+0.1571; 7:-0.1841; 10:+0.2579; 13:-0.0745; 16:+0.2620; 19:-0.0837; 22:-0.1618; 25:-0.0455; 28:+0.0134; 31:+0.0039.
- Mean |E_(h,k)|/p diagnostics: even d=0.03762, odd d=0.03487, d<=10=0.03693, d>10=0.03536, k=2h=0.03124, k!=2h=0.03637.
- h=1 cross row: mean |E_(1,k)|/p=0.05477; maximum |E_(1,k)|/p=0.12339 at k=5.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.8245 q10=+0.8794 q20=+0.9151 q30=+0.9240 q40=+0.9400 q50=+0.9559 q60=+0.9705 q70=+0.9771 q80=+0.9934 q90=+1.0062 q100=+1.9900.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-1.1562 q10=-0.3716 q20=-0.2039 q30=-0.1282 q40=-0.0810 q50=-0.0441 q60=-0.0092 q70=+0.0016 q80=+0.0220 q90=+0.0600 q100=+31.3533.
- The same normalized deciles restricted to generic h>=2: q0=-1.1562 q10=-0.3857 q20=-0.2157 q30=-0.1311 q40=-0.0874 q50=-0.0482 q60=-0.0101 q70=+0.0003 q80=+0.0111 q90=+0.0336 q100=+0.1439.
- h=1 self row: off-diagonal count=1984 (1.989970 p); E_(1,1)/p=+1.992970.

### H=63

All 1953 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 68056864/997 = 68261.648947.
- sum_(h<k) |E_(h,k)|/(p H) = 1.086779.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.544688.
- max_d |A_d|/p = 0.640733, attained at d=5.
- V(H) = 1007261686112/20749131733701 = 0.048545.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 5: 636894/997; +0.640733; 22: -611480/997; -0.615165; 15: 514266/997; +0.517366; 16: 448498/997; +0.451201; 24: 382308/997; +0.384612.
- Signed A_d/p distribution deciles: q0=-0.6152 q10=-0.2809 q20=-0.1442 q30=-0.0687 q40=-0.0147 q50=+0.0045 q60=+0.0341 q70=+0.0697 q80=+0.1896 q90=+0.2735 q100=+0.6407.
- Ordered-d profile samples (d: A_d/p): 7:-0.2958; 13:+0.1341; 19:-0.0543; 25:-0.2394; 31:+0.1945; 38:+0.1703; 44:-0.1478; 50:+0.3189; 56:-0.0568; 62:-0.0834.
- Mean |E_(h,k)|/p diagnostics: even d=0.03530, odd d=0.03482, d<=10=0.03555, d>10=0.03485, k=2h=0.02788, k!=2h=0.03517.
- h=1 cross row: mean |E_(1,k)|/p=0.05708; maximum |E_(1,k)|/p=0.17434 at k=58.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.7803 q10=+0.8229 q20=+0.8818 q30=+0.9027 q40=+0.9204 q50=+0.9268 q60=+0.9416 q70=+0.9657 q80=+0.9741 q90=+0.9950 q100=+1.9900.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-1.1562 q10=-0.2039 q20=-0.1110 q30=-0.0801 q40=-0.0460 q50=-0.0272 q60=-0.0179 q70=-0.0080 q80=+0.0075 q90=+0.0271 q100=+31.3533.
- The same normalized deciles restricted to generic h>=2: q0=-1.1562 q10=-0.2098 q20=-0.1126 q30=-0.0805 q40=-0.0470 q50=-0.0280 q60=-0.0196 q70=-0.0090 q80=+0.0036 q90=+0.0244 q100=+0.1439.
- h=1 self row: off-diagonal count=1984 (1.989970 p); E_(1,1)/p=+1.992970.

## p=1999 (p mod 3 = 1)

H values: 21, 45, 96.

### H=21

All 210 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 22557502/1999 = 11284.393197.
- sum_(h<k) |E_(h,k)|/(p H) = 0.268810.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.572312.
- max_d |A_d|/p = 0.177571, attained at d=9.
- V(H) = 2478660870916/335328503832021 = 0.007392.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 9: -709574/1999; -0.177571; 13: 650618/1999; +0.162817; 2: 627922/1999; +0.157138; 3: 460282/1999; +0.115186; 1: 423746/1999; +0.106043.
- Signed A_d/p distribution deciles: q0=-0.1776 q10=-0.0934 q20=-0.0627 q30=-0.0401 q40=-0.0071 q50=+0.0018 q60=+0.0233 q70=+0.0638 q80=+0.0755 q90=+0.1194 q100=+0.1628.
- Ordered-d profile samples (d: A_d/p): 2:+0.1571; 4:-0.0923; 6:-0.0017; 8:-0.0371; 10:+0.0646; 12:-0.0869; 14:+0.0678; 16:+0.0634; 18:-0.0061; 20:-0.0085.
- Mean |E_(h,k)|/p diagnostics: even d=0.02632, odd d=0.02739, d<=10=0.02552, d>10=0.03072, k=2h=0.01570, k!=2h=0.02744.
- h=1 cross row: mean |E_(1,k)|/p=0.04717; maximum |E_(1,k)|/p=0.13856 at k=15.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.8944 q10=+0.9485 q20=+0.9585 q30=+0.9665 q40=+0.9725 q50=+0.9825 q60=+0.9965 q70=+1.0215 q80=+1.0385 q90=+1.0725 q100=+1.9950.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-0.6152 q10=-0.4452 q20=-0.1471 q30=-0.0871 q40=-0.0586 q50=-0.0322 q60=+0.0142 q70=+0.1399 q80=+0.2185 q90=+0.8890 q100=+44.5535.
- The same normalized deciles restricted to generic h>=2: q0=-0.6152 q10=-0.4526 q20=-0.1480 q30=-0.0973 q40=-0.0626 q50=-0.0345 q60=-0.0044 q70=+0.1015 q80=+0.1778 q90=+0.3630 q100=+1.1182.
- h=1 self row: off-diagonal count=3988 (1.994997 p); E_(1,1)/p=+1.996496.

### H=45

All 990 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 102665894/1999 = 51358.626313.
- sum_(h<k) |E_(h,k)|/(p H) = 0.570937.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.567260.
- max_d |A_d|/p = 0.402124, attained at d=16.
- V(H) = 10559599219052/718561079640045 = 0.014695.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 16: 1606888/1999; +0.402124; 13: 1545840/1999; +0.386847; 25: 1213348/1999; +0.303641; 29: -878540/1999; -0.219855; 14: 607350/1999; +0.151989.
- Signed A_d/p distribution deciles: q0=-0.2199 q10=-0.0872 q20=-0.0392 q30=-0.0243 q40=-0.0134 q50=+0.0038 q60=+0.0240 q70=+0.0518 q80=+0.0950 q90=+0.1366 q100=+0.4021.
- Ordered-d profile samples (d: A_d/p): 5:-0.0285; 9:-0.0685; 14:+0.1520; 18:-0.0061; 22:+0.1369; 27:+0.0109; 31:+0.0156; 36:+0.1359; 40:-0.0568; 44:-0.0185.
- Mean |E_(h,k)|/p diagnostics: even d=0.02557, odd d=0.02632, d<=10=0.02576, d>10=0.02608, k=2h=0.01801, k!=2h=0.02613.
- h=1 cross row: mean |E_(1,k)|/p=0.04182; maximum |E_(1,k)|/p=0.14811 at k=42.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.8664 q10=+0.9139 q20=+0.9441 q30=+0.9577 q40=+0.9665 q50=+0.9735 q60=+0.9881 q70=+1.0013 q80=+1.0133 q90=+1.0385 q100=+1.9950.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-0.6152 q10=-0.1454 q20=-0.0915 q30=-0.0684 q40=-0.0440 q50=-0.0196 q60=+0.0050 q70=+0.0233 q80=+0.0418 q90=+0.1982 q100=+44.5535.
- The same normalized deciles restricted to generic h>=2: q0=-0.6152 q10=-0.1458 q20=-0.0943 q30=-0.0687 q40=-0.0462 q50=-0.0199 q60=+0.0003 q70=+0.0163 q80=+0.0413 q90=+0.1593 q100=+1.1182.
- h=1 self row: off-diagonal count=3988 (1.994997 p); E_(1,1)/p=+1.996496.

### H=96

All 4560 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 449022126/1999 = 224623.374687.
- sum_(h<k) |E_(h,k)|/(p H) = 1.170499.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.545138.
- max_d |A_d|/p = 0.549291, attained at d=13.
- V(H) = 15280842403861/383232575808024 = 0.039874.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 13: 2194968/1999; +0.549291; 29: -2148854/1999; -0.537751; 20: -2101854/1999; -0.525989; 56: 1966998/1999; +0.492242; 58: 1863512/1999; +0.466344.
- Signed A_d/p distribution deciles: q0=-0.5378 q10=-0.1929 q20=-0.0984 q30=-0.0701 q40=-0.0329 q50=+0.0228 q60=+0.0570 q70=+0.1134 q80=+0.1826 q90=+0.2589 q100=+0.5493.
- Ordered-d profile samples (d: A_d/p): 10:+0.1490; 19:-0.2540; 29:-0.5378; 38:-0.0323; 48:+0.2401; 57:+0.2509; 67:-0.1419; 76:+0.1734; 86:-0.0703; 95:-0.0351.
- Mean |E_(h,k)|/p diagnostics: even d=0.02473, odd d=0.02455, d<=10=0.02440, d>10=0.02470, k=2h=0.02284, k!=2h=0.02466.
- h=1 cross row: mean |E_(1,k)|/p=0.03388; maximum |E_(1,k)|/p=0.14811 at k=42.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.8304 q10=+0.8789 q20=+0.8944 q30=+0.9175 q40=+0.9425 q50=+0.9515 q60=+0.9665 q70=+0.9740 q80=+0.9905 q90=+1.0140 q100=+1.9950.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-0.6152 q10=-0.0957 q20=-0.0625 q30=-0.0476 q40=-0.0368 q50=-0.0222 q60=-0.0119 q70=+0.0013 q80=+0.0132 q90=+0.0414 q100=+44.5535.
- The same normalized deciles restricted to generic h>=2: q0=-0.6152 q10=-0.0970 q20=-0.0629 q30=-0.0478 q40=-0.0372 q50=-0.0226 q60=-0.0136 q70=+0.0006 q80=+0.0125 q90=+0.0407 q100=+1.1182.
- h=1 self row: off-diagonal count=3988 (1.994997 p); E_(1,1)/p=+1.996496.

## p=4001 (p mod 3 = 2)

H values: 28, 64, 145.

### H=28

All 378 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 97953496/4001 = 24482.253437.
- sum_(h<k) |E_(h,k)|/(p H) = 0.218537.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.493686.
- max_d |A_d|/p = 0.217571, attained at d=18.
- V(H) = 12533910067972/1793792672112007 = 0.006987.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 18: 3482876/4001; +0.217571; 2: 3220292/4001; +0.201168; 5: 2889594/4001; +0.180509; 6: -2003176/4001; -0.125136; 16: -1715070/4001; -0.107138.
- Signed A_d/p distribution deciles: q0=-0.1251 q10=-0.0952 q20=-0.0375 q30=-0.0110 q40=-0.0006 q50=+0.0085 q60=+0.0245 q70=+0.0372 q80=+0.0609 q90=+0.1234 q100=+0.2176.
- Ordered-d profile samples (d: A_d/p): 3:+0.0085; 6:-0.1251; 9:+0.0327; 11:+0.0402; 14:+0.0040; 17:+0.0628; 19:-0.0091; 22:+0.0337; 25:+0.0365; 27:+0.0005.
- Mean |E_(h,k)|/p diagnostics: even d=0.01661, odd d=0.01579, d<=10=0.01601, d>10=0.01645, k=2h=0.01729, k!=2h=0.01615.
- h=1 cross row: mean |E_(1,k)|/p=0.00045; maximum |E_(1,k)|/p=0.00175 at k=7.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.0000 q10=+0.9580 q20=+0.9657 q30=+0.9799 q40=+0.9877 q50=+0.9943 q60=+0.9980 q70=+1.0025 q80=+1.0152 q90=+1.0306 q100=+1.0482.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-63.2061 q10=-0.1683 q20=-0.1078 q30=-0.0457 q40=-0.0335 q50=-0.0039 q60=+0.0092 q70=+0.0491 q80=+0.0847 q90=+0.1484 q100=+0.6429.
- The same normalized deciles restricted to generic h>=2: q0=-1.1620 q10=-0.1423 q20=-0.0952 q30=-0.0385 q40=-0.0301 q50=+0.0027 q60=+0.0099 q70=+0.0548 q80=+0.0889 q90=+0.1576 q100=+0.6429.
- h=1 self row: off-diagonal count=0 (0.000000 p); E_(1,1)/p=+0.000749.

### H=64

All 2016 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 554762932/4001 = 138656.068983.
- sum_(h<k) |E_(h,k)|/(p H) = 0.541490.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.535174.
- max_d |A_d|/p = 0.316976, attained at d=18.
- V(H) = 33747774141315/2050048768128008 = 0.016462.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 18: 5074150/4001; +0.316976; 20: 4482068/4001; +0.279989; 3: 4325726/4001; +0.270223; 25: 4171788/4001; +0.260606; 34: -4066294/4001; -0.254016.
- Signed A_d/p distribution deciles: q0=-0.2540 q10=-0.1358 q20=-0.1074 q30=-0.0654 q40=-0.0154 q50=+0.0051 q60=+0.0220 q70=+0.0434 q80=+0.0879 q90=+0.1931 q100=+0.3170.
- Ordered-d profile samples (d: A_d/p): 7:-0.0691; 13:+0.1091; 19:+0.2439; 26:+0.0220; 32:+0.0175; 38:-0.1222; 45:+0.0191; 51:+0.0390; 57:+0.0383; 63:-0.0000.
- Mean |E_(h,k)|/p diagnostics: even d=0.01678, odd d=0.01759, d<=10=0.01728, d>10=0.01715, k=2h=0.02063, k!=2h=0.01713.
- h=1 cross row: mean |E_(1,k)|/p=0.00039; maximum |E_(1,k)|/p=0.00175 at k=7.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.0000 q10=+0.9489 q20=+0.9598 q30=+0.9707 q40=+0.9800 q50=+0.9903 q60=+0.9948 q70=+0.9978 q80=+1.0074 q90=+1.0186 q100=+1.0482.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-63.2061 q10=-0.1094 q20=-0.0582 q30=-0.0345 q40=-0.0248 q50=-0.0064 q60=+0.0084 q70=+0.0131 q80=+0.0314 q90=+0.0848 q100=+0.6429.
- The same normalized deciles restricted to generic h>=2: q0=-1.1620 q10=-0.0970 q20=-0.0558 q30=-0.0337 q40=-0.0231 q50=-0.0025 q60=+0.0086 q70=+0.0131 q80=+0.0340 q90=+0.0866 q100=+0.6429.
- h=1 self row: off-diagonal count=0 (0.000000 p); E_(1,1)/p=+0.000749.

### H=145

All 10440 cross-gap pairs were evaluated; no k-h restriction was used.

- sum_(h<k) |E_(h,k)| = 2892875464/4001 = 723038.106473.
- sum_(h<k) |E_(h,k)|/(p H) = 1.246306.
- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = 0.543677.
- max_d |A_d|/p = 0.797503, attained at d=3.
- V(H) = 1466354706293576/37157133922320145 = 0.039464.
- Top five shells by |A_d| (d: A_d exact; A_d/p): 3: 12766424/4001; +0.797503; 98: 8049876/4001; +0.502866; 6: -7966154/4001; -0.497636; 110: 7395302/4001; +0.461975; 69: 7188274/4001; +0.449043.
- Signed A_d/p distribution deciles: q0=-0.4976 q10=-0.2258 q20=-0.1266 q30=-0.0829 q40=-0.0355 q50=-0.0025 q60=+0.0457 q70=+0.1289 q80=+0.1861 q90=+0.2550 q100=+0.7975.
- Ordered-d profile samples (d: A_d/p): 15:+0.2863; 29:-0.1443; 44:-0.0625; 58:-0.0644; 72:-0.0135; 87:-0.0160; 101:-0.1319; 116:+0.2095; 130:-0.0768; 144:-0.0008.
- Mean |E_(h,k)|/p diagnostics: even d=0.01745, odd d=0.01718, d<=10=0.01788, d>10=0.01722, k=2h=0.01846, k!=2h=0.01730.
- h=1 cross row: mean |E_(1,k)|/p=0.00044; maximum |E_(1,k)|/p=0.00202 at k=120.
- Same-gap off-diagonal raw-count/p distribution deciles: q0=+0.0000 q10=+0.9127 q20=+0.9274 q30=+0.9414 q40=+0.9524 q50=+0.9623 q60=+0.9758 q70=+0.9837 q80=+0.9948 q90=+1.0051 q100=+1.0482.
- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: q0=-63.2061 q10=-0.0533 q20=-0.0412 q30=-0.0315 q40=-0.0254 q50=-0.0158 q60=-0.0093 q70=+0.0024 q80=+0.0086 q90=+0.0255 q100=+0.6429.
- The same normalized deciles restricted to generic h>=2: q0=-1.1620 q10=-0.0524 q20=-0.0399 q30=-0.0310 q40=-0.0253 q50=-0.0157 q60=-0.0083 q70=+0.0028 q80=+0.0087 q90=+0.0255 q100=+0.6429.
- h=1 self row: off-diagonal count=0 (0.000000 p); E_(1,1)/p=+0.000749.

## Five verdicts

VERDICT 1 (coverage and arithmetic): PASS. All nine (p,H) cases and all h<k pairs were computed, without the optional k-h<=200 truncation. The h=1 Casoratian identity delta_1(r)=(r+1)^(-3) was independently asserted at every admissible r.

VERDICT 2 (cross-gap L1 mass): NOT SUPPORTED as a uniform asymptotic bound sum|E| << pH. The observed ratio sum|E|/(pH) ranges from 0.2185 to 1.2463, while its rescaling by sqrt(p)/H is strikingly stable in [0.4937,0.5723]. Thus the data fit sum|E| about c H^2 sqrt(p), the generic absolute-noise scale, which exceeds pH once H grows past sqrt(p). At the nine finite test points the unscaled ratio is still at most the displayed 1.2463. This L1 behavior does not contradict signed shell cancellation or the square-mean target.

VERDICT 3 (shell variance and structured d): No stable anomalous shell was detected. Across all cases max_d |A_d|/p <= 0.7975 and V(H) <= 0.0485. The largest within-case ratios between the paired mean-|E| diagnostics were even/odd=1.079, small/large d=1.210, and k=2h/non-doubling=1.748; top-five d values move with p and H rather than identifying a persistent parity, small-d, or k=2h mechanism.

VERDICT 4 (same-gap component census): PASS for h>=2. The centered self rows have the predicted +p main term after the diagonal is included and the random baseline is subtracted; the displayed (E_(h,h)-p)/(h sqrt(p)) deciles are bounded at all nine scales. The separately displayed raw off-diagonal counts are correspondingly near p.

VERDICT 5 (h=1 cube-root exception): PASS. For p=997 and 1999, both 1 mod 3, E_(1,1)/p lies near 2 (observed range 1.9930 to 1.9965); for p=4001, 2 mod 3, it lies near 0 (observed range 0.0007 to 0.0007). The cross row also detects the cube map: mean |E_(1,k)|/p is 0.0339--0.0571 for p=1 mod 3, but only 0.0004--0.0005 for p=2 mod 3, where cubing is bijective and only window-boundary errors remain.

The empirical conclusion is therefore mixed in the intended useful way: the stronger cross-gap L1 statement is not supported, while the signed fixed-difference shell statistics, same-gap component census, and the cube-root exception behave as predicted.
