# Q4313 independent fixed-gcd census

Implementation: independent Java standard-library exact arithmetic.

```text
min q                   15001
max q                   20000
primes through qmax     2262
selected states         560
sign records            1120
largest Apery index     258743
largest bit length      1315993
elapsed                  00:44:40
plus ledger differences 0
plus extra q powers     0
minus extra q powers    0
```

## Prefix q <= 20000

Selected states: 560; sign records: 1120; minus raw-window states: 70; M5-capable states: 40.

### plus

```text
records                 560
G=1                     351
median G                1
median nontrivial G     19
maximum G               21299 = 19^2*59
maximum state           (15683,14400)
maximum omega           3
largest prime factor    2543 at (15439,2185)
rho(log log G~log q)  0.470141377  R2=0.0136320943
slope(log G~q)        0.000112388009  R2=0.0125022365
size certified          560 / 560
max logG/(2log pmin)   0.516431788742 at (15683,14400)
```

### minus

```text
records                 560
G=1                     343
median G                1
median nontrivial G     19
maximum G               72539 = 17^2*251
maximum state           (16087,639)
maximum omega           3
largest prime factor    1973 at (19273,9815)
rho(log log G~log q)  0.423390234  R2=0.00844525611
slope(log G~q)        0.000106610091  R2=0.00745539298
size certified          560 / 560
max logG/(2log pmin)   0.586673642894 at (16087,639)
```

### minus raw-window

```text
records                 70
G=1                     40
median G                1
median nontrivial G     19
maximum G               72539 = 17^2*251
maximum state           (16087,639)
maximum omega           2
largest prime factor    251 at (16087,639)
rho(log log G~log q)  0.895732256  R2=0.0378323286
slope(log G~q)        0.000191407075  R2=0.0203159171
size certified          70 / 70
max logG/(2log pmin)   0.586673642894 at (16087,639)
```

### minus M5-capable

```text
records                 40
G=1                     24
median G                1
median nontrivial G     17
maximum G               72539 = 17^2*251
maximum state           (16087,639)
maximum omega           2
largest prime factor    251 at (16087,639)
rho(log log G~log q)  -0.520436244  R2=0.0133492040
slope(log G~q)        -0.000228378275  R2=0.0260318573
size certified          40 / 40
max logG/(2log pmin)   0.586673642894 at (16087,639)
```

## Direct raw-leaf audit

```text
raw plus leaves        0
raw minus leaves       0
assigned plus leaves   0
assigned minus leaves  0
raw minus M5 leaves    0
```

