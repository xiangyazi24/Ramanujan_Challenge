# Q4316 extension analysis

Exact input: the independently generated `records.tsv`.

## Cumulative through q <= 20000

### plus

```text
records                 560
G=1                     351
median G                1
median nontrivial G     19
maximum G               21299 = 19^2*59 at (15683,14400)
maximum omega           3
largest prime factor    2543 at (15439,2185)
rho(log log G~log q)  0.470141377  R2=0.0136320943  n=209
size certified          560 / 560
raw leaves              0
assigned leaves         0
```

### minus all

```text
records                 560
G=1                     343
median G                1
median nontrivial G     19
maximum G               72539 = 17^2*251 at (16087,639)
maximum omega           3
largest prime factor    1973 at (19273,9815)
rho(log log G~log q)  0.423390234  R2=0.00844525611  n=217
size certified          560 / 560
raw leaves              0
assigned leaves         0
```

### minus raw-window

```text
records                 70
G=1                     40
median G                1
median nontrivial G     19
maximum G               72539 = 17^2*251 at (16087,639)
maximum omega           2
largest prime factor    251 at (16087,639)
rho(log log G~log q)  0.895732256  R2=0.0378323286  n=30
size certified          70 / 70
raw leaves              0
assigned leaves         0
```

### minus M5-capable

```text
records                 40
G=1                     24
median G                1
median nontrivial G     17
maximum G               72539 = 17^2*251 at (16087,639)
maximum omega           2
largest prime factor    251 at (16087,639)
rho(log log G~log q)  -0.520436244  R2=0.0133492040  n=16
size certified          40 / 40
raw leaves              0
assigned leaves         0
```

## New band 15000 < q <= 20000

### plus

```text
records                 560
G=1                     351
median G                1
median nontrivial G     19
maximum G               21299 = 19^2*59 at (15683,14400)
maximum omega           3
largest prime factor    2543 at (15439,2185)
rho(log log G~log q)  0.470141377  R2=0.0136320943  n=209
size certified          560 / 560
raw leaves              0
assigned leaves         0
```

### minus all

```text
records                 560
G=1                     343
median G                1
median nontrivial G     19
maximum G               72539 = 17^2*251 at (16087,639)
maximum omega           3
largest prime factor    1973 at (19273,9815)
rho(log log G~log q)  0.423390234  R2=0.00844525611  n=217
size certified          560 / 560
raw leaves              0
assigned leaves         0
```

### minus raw-window

```text
records                 70
G=1                     40
median G                1
median nontrivial G     19
maximum G               72539 = 17^2*251 at (16087,639)
maximum omega           2
largest prime factor    251 at (16087,639)
rho(log log G~log q)  0.895732256  R2=0.0378323286  n=30
size certified          70 / 70
raw leaves              0
assigned leaves         0
```

### minus M5-capable

```text
records                 40
G=1                     24
median G                1
median nontrivial G     17
maximum G               72539 = 17^2*251 at (16087,639)
maximum omega           2
largest prime factor    251 at (16087,639)
rho(log log G~log q)  -0.520436244  R2=0.0133492040  n=16
size certified          40 / 40
raw leaves              0
assigned leaves         0
```

## Required point checks

```text
minus state (5647,4553)   MISSING
minus maximum through 15000 G=72539 = 17^2*251 at (16087,639)
plus largest factor      2543 at (15439,2185)
minus all largest factor 1973 at (19273,9815)
minus raw largest factor 251 at (16087,639)
minus M5 largest factor  251 at (16087,639)
```

## Largest retained prime factor by dyadic q block

| q block | plus | minus all | minus raw | minus M5 |
|---|---:|---:|---:|---:|
|[16,32)|0|0|0|0|
|[32,64)|0|0|0|0|
|[64,128)|0|0|0|0|
|[128,256)|0|0|0|0|
|[256,512)|0|0|0|0|
|[512,1024)|0|0|0|0|
|[1024,2048)|0|0|0|0|
|[2048,4096)|0|0|0|0|
|[4096,8192)|0|0|0|0|
|[8192,16384)|2543|251|251|251|
|[16384,32768)|499|1973|193|31|

## Dyadic-maximum trend diagnostics

- plus: log(block max factor) slope -2.34942, R^2 1.00000, blocks 2.
- minus all: log(block max factor) slope 2.97463, R^2 1.00000, blocks 2.
- minus raw-window: log(block max factor) slope -0.379087, R^2 1.00000, blocks 2.
- minus M5-capable: log(block max factor) slope -3.01735, R^2 1.00000, blocks 2.
