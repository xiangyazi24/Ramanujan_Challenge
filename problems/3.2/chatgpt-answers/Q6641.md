ANSWER Q6641 ab3cdae8

Reprint of Q6597 null-model prediction table for the supercritical Q_D experiment.

```text
MODEL SETUP

Let N = p - 1 be the orbit length. Treat the projective orbit colors as
independent random colors in F_p, with the mirror involution

    r <-> p-1-r

removed as a deterministic constraint.

Define

    C_p(a,g) = #{r : pi(r)=pi(r+a)=pi(r+a+g)}

and

    Q_D = sum_{a>=1,g>=1,a+g<=D} C_p(a,g).

The number of lag pairs is

    M_D = D(D-1)/2.

Pure random model:
    all positions independent colors with probability 1/p.

Mirror corrected model:
    each color class has a forced reflected partner, giving an even-gap
    deterministic collision contribution.
```

## 1. Expectation and variance of Q_D

For a fixed triple of positions, pure random coloring gives

```text
P[three positions have same color] = 1/p^2.
```

There are approximately

```text
N * M_D
```

ordered starting triples, hence

```text
E_random[Q_D]

    = (p-1) D(D-1)/(2 p^2)
    = D(D-1)/(2p) + O(D^2/p^2).
```

For the variance, write Q_D as a sum of indicator variables.

The diagonal contribution is

```text
N*M_D*(1/p^2)*(1-1/p^2)
```

and the overlap terms are dominated by pairs of triples sharing one or two
positions.

The leading scale is therefore

```text
Var_random(Q_D)

    = D^2/(2p) + O(D^3/p^2).
```

Equivalently, for

```text
D = sqrt(p)*L
```

we obtain

```text
E_random[Q_D] ~ L^2/2

Var_random(Q_D) ~ L^2/2.
```

### Mirror correction

The mirror involution forces

```text
pi(r)=pi(p-1-r).
```

For an even lag h:

```text
h = 2t
```

a pair separated by h already has one deterministic collision channel.

The correction changes the effective triple probability:

```text
P_mirror[triple collision]

 = 1/p^2 + 1/(2p) * 1_{all relevant gaps even}
 + lower order terms.
```

Therefore the corrected expectation is

```text
E_mirror[Q_D]

 = E_random[Q_D]
   + (1/4)*D^2/p
   + O(D/p).
```

The 1/4 comes from:

```text
half of the lag pairs are even,
each even lag has one forced partner,
and the third equality costs probability ~1/p.
```

Thus:

```text
E_mirror[Q_D]

 = 3D^2/(4p) + O(D^2/p^2+D/p).
```

At

```text
D=sqrt(p)*L
```

this predicts

```text
E_mirror[Q_D] ~ 3L^2/4.
```

The variance remains Poisson scale:

```text
Var_mirror(Q_D)

 = 3D^2/(4p) + O(D^3/p^2).
```

---

## 2. Distribution of individual C_p(a,g)

For a fixed pair (a,g), the number of possible starting positions is

```text
N-a-g.
```

Pure random model:

```text
C_p(a,g) ~ Poisson(lambda)

lambda = (p-1-a-g)/p^2

      = 1/p + O(D/p^2).
```

Hence individual triple collisions are extremely sparse.

With mirror correction:

```text
lambda(a,g)

 = 1/p
   + 1/(2p)*1_{a,g parity compatible}
   + O(1/p^2).
```

For the M_D approximately independent strata, the maximum satisfies a
Gumbel extreme-value law.

If X_i are Poisson(lambda), the maximum m solves

```text
M_D * P[X>=m] ~ 1.
```

Using the Poisson tail approximation:

```text
m log(m/(e lambda)) = log M_D.
```

Therefore

```text
m_max

 ~ log(M_D)/log(log(M_D)/lambda).
```

For

```text
D=sqrt(p) log p
```

we have

```text
M_D ~ p(log p)^2/2.
```

Therefore

```text
log M_D ~ log p
```

and

```text
m_max ~ log p / log(log p).
```

A 3-sigma anomaly is:

```text
observed max > mu_max + 3 sigma_max
```

where the Gumbel centering is

```text
b_M = solution of

M_D*P(Poisson(lambda)>=b_M)=1
```

and scale

```text
a_M = 1/log(b_M/lambda).
```

---

# 3. Per-stratum predictions

## (i) Mirror forced lags

Even gaps:

```text
a+g even
```

have deterministic contribution.

Mean:

```text
E[C_p(a,g)]

 = 1/(2p)+1/p^2.
```

Maximum:

```text
max C_p(a,g)

 ~ 1 + Poisson extreme tail.
```

A value:

```text
C>=4
```

is already unusual.

---

## (ii) Axis strips

Define

```text
min(a,g)<=G
```

with

```text
G=N^(1/3).
```

Number of pairs:

```text
M_axis ~ 2GD.
```

Expected mass:

```text
E[Q_axis]

 ~ 2GD/p.
```

For D=sqrt(p)L:

```text
E[Q_axis]

 ~ 2p^(-1/6)*L.
```

The axis strip is null-small.

---

## (iii) a or g dividing p-1

Number of divisor lags:

```text
tau(p-1).
```

Expected contribution:

```text
E_div(p-1)

 ~ tau(p-1)*D/p^2.
```

Since

```text
tau(n)=n^(o(1)),
```

this is negligible.

---

## (iv) a or g dividing p+1

Identical:

```text
E_div(p+1)

 ~ tau(p+1)*D/p^2.
```

The p+1 family is a resonance candidate because of the
Q(sqrt(-51)) splitting layer.

---

## (v) Parity classes

Four classes:

```text
(odd,odd)
(odd,even)
(even,odd)
(even,even)
```

Pure random:

```text
each receives 1/4 of Q_D.
```

Mirror corrected:

```text
(even,even):

    mean = 1/4 random + forced term

other three:

    mean = 1/4 random.
```

Thus:

```text
Q_even_even/Q_total

 -> approximately 1/2.
```

rather than 1/4.

---

## (vi) Dyadic ratio classes

For

```text
2^j <= a/g < 2^(j+1)
```

the number of pairs is proportional to the area.

Hence:

```text
E[Q_j]

 = area_j/p.
```

All dyadic ratio bands should have the same normalized density:

```text
Q_j/area_j

 approximately constant.
```

A deviation larger than sqrt(area) fluctuations indicates a resonance.

---

# 4. Detecting moving resonance families

A thin family h=h(p) contributes only O(1) lags.

Therefore pooled averages dilute it.

The best statistic is:

```text
S1 = number of primes where the stratum maximum exceeds its
     null 99% threshold.
```

For each prime:

```text
flag_p = 1 if max_stratum(p)>T_99(p).
```

Under the null:

```text
flag_p ~ Bernoulli(0.01).
```

For P primes:

```text
sum flag_p ~ Binomial(P,0.01).
```

For approximately 150 primes:

```text
mean = 1.5
sigma = sqrt(1.485)=1.22.
```

Therefore:

```text
10 or more flagged primes
```

is overwhelmingly non-null.

Second best statistic:

```text
pooled stratum z-score

Z=(observed-null mean)/sqrt(null variance).
```

Per-prime maximum is preferred because a moving resonance changes location.

---

# 5. Raw resultant contamination

Suppose:

```text
p | b_j
```

for some j.

Then every pair:

```text
j <= h < k <= D
```

is polluted.

Triangle size:

```text
T_j=(D-j)(D-j-1)/2.
```

Under the null:

```text
P[p|b_j]~1/p.
```

Expected number of polluted triangles:

```text
E[number]

 = D/p.
```

For

```text
D=sqrt(p)log p
```

we get:

```text
E[number]

 = log(p)/sqrt(p).
```

This tends to zero.

However, conditional on one zero, the contamination is huge:

```text
polluted entries

 ~ D^2/2.
```

Hence raw resultants have a heavy-tailed false signal:

```text
rare event probability ~1/sqrt(p)
large amplitude ~p log^2(p).
```

The regular-root count removes this amplification.

Therefore:

```text
raw resultant statistics are non-diagnostic
at D=sqrt(p)log p.
```

---

# Final audit answers

```text
(a) Least-confident derivation:

The weakest step is the mirror-corrected variance formula.

The expectation correction is robust because it is a first-moment
pair-counting argument.

The variance requires controlling dependencies between overlapping
triples after quotienting by the mirror involution. The stated formula

    Var(Q_D) ~ 3D^2/(4p)

is a null-model prediction, not a theorem.

The computation should explicitly compare the empirical covariance matrix
of C_p(a,g) strata with the predicted overlap covariance.
```

```text
(b) Likely shared blind spot:

The largest danger is not the mean but hidden low-dimensional strata.

A thin family of O(1) lags can carry the whole anomaly while disappearing
inside band averages.

Therefore the experiment should always report:

1. maximum per stratum,
2. number of primes exceeding null quantiles,
3. location of exceptional lags,

not only Q_D averages.

The raw resultant experiment is especially vulnerable because one Apery
zero creates a quadratic-size fake signal.
```