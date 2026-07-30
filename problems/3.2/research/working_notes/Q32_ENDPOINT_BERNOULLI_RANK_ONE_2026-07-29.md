# Rank-one Bernoulli defect at all Apéry block endpoints

## 0. Result

Let

```text
b_0=1, b_1=5,
(n+1)^3 b_(n+1)
 =(34n^3+51n^2+27n+5)b_n-n^3b_(n-1)
```

be the Apéry numbers for `zeta(3)`.  For every prime `p>=5` and every
integer `m>=1`,

```text
b_(mp)
 ==b_m
   -m^3(17b_m-b_(m-1))p^3 B_(p-3)/18              (mod p^4),  (D)

b_(mp-1)
 ==b_(m-1)
   +m^3(17b_(m-1)-b_m)p^3 B_(p-3)/18              (mod p^4).  (R)
```

All rational congruences are in `Z_(p)`; the displayed denominators
are units because `p>=5`.

Put

```text
e_p=(b_(p-1)-1)/p^3                                (mod p),

E_m=m^3(b_(m-1)-17b_m)/12,
F_m=m^3(17b_(m-1)-b_m)/12.
```

The case `m=1` of `(R)` gives

```text
e_p==2B_(p-3)/3                                    (mod p),
```

and `(D)`--`(R)` become the rank-one endpoint laws

```text
b_(mp)  ==b_m     +p^3 e_p E_m                    (mod p^4),
b_(mp-1)==b_(m-1) +p^3 e_p F_m                    (mod p^4).  (0.1)
```

Thus the two endpoint defects contain only one prime-dependent scalar.
The coefficient vector has the closed form

```text
[E_m]    m^3 [-17   1] [b_m    ]
[F_m]  = --- [ -1  17] [b_(m-1)],
           12
```

and the inner matrix squares to `288 I`.

There is a stronger one-order lift.  For `p>=7`, put

```text
Delta_p=b_(p-1)-1.
```

Then, for every `m>=1`,

```text
b_(mp)-b_m     ==E_m Delta_p                       (mod p^5),
b_(mp-1)-b_(m-1)==F_m Delta_p                      (mod p^5). (0.2)
```

Thus the endpoint defect is rank one through two consecutive
`p`-adic digits.  The prime `p=5` satisfies `(D)`--`(R)` but is an
exception to `(0.2)`.

The special cases

```text
b_p    ==5-(14/3)p^3B_(p-3),
b_(2p) ==73-(1648/3)p^3B_(p-3),
b_(2p-1)==5+(16/3)p^3B_(p-3),
b_(3p) ==1445-36738p^3B_(p-3)                    (mod p^4)
```

recover the formulas proved separately in the recent
Apéry-supercongruence literature.

The direct congruence `(D)` also has a precise published precursor.
The `(r,s,t)=(2,2,0)` specialization of Ji-Cai Liu, *An extension of
Gauss congruences for Apéry numbers*, arXiv:2404.16636, Theorem 1.1,
gives

```text
b_(mp^q)-b_(mp^(q-1))
 ==(2/3)p^(3q)B_(p-3)E_m                         (mod p^(3q+1)).
```

At `q=1` this is `(D)`, after the finite sum in Liu's coefficient is
reduced by `(2.5)`.  The reflected companion `(R)`, the two-entry
matrix form, and especially the strong law `(0.2)` are not supplied
by that prime-power-index theorem; Sections 2--4 give a self-contained
derivation of all of them.

## 1. Three standard harmonic congruences

Write `H=H_(p-1)`, `H2=H_(p-1)^(2)`, and
`H3=H_(p-1)^(3)`.  For `p>=5`,

```text
H  ==-(p^2/3)B_(p-3)                              (mod p^3),
H2 == (2p/3)B_(p-3)                               (mod p^2),
H3 ==0                                            (mod p).    (1.1)
```

We will also use the following refinement of Jacobsthal's binomial
congruence.  If `a=b+c` are nonnegative integers, then

```text
binom(ap,bp)
 ==binom(a,b)[1-abc p^3 B_(p-3)/3]                (mod p^4).  (1.2)
```

For completeness, remove the multiples of `p` from the three
factorials.  The logarithm of the unit block indexed by `t` is

```text
sum_(s=1)^(p-1) log(1+tp/s)
 ==tpH-(t^2p^2/2)H2
 ==-t(t+1)p^3B_(p-3)/3                            (mod p^4)
```

by `(1.1)`; the cubic-log term vanishes by `H3==0 (mod p)`.
Furthermore

```text
sum_(t=0)^(a-1)t(t+1)
-sum_(t=0)^(b-1)t(t+1)
-sum_(t=0)^(c-1)t(t+1)
=abc.
```

Exponentiating proves `(1.2)`.  This argument does not require
`a,b,c<p`.

## 2. Direct blocks

Set

```text
T(n,k)=binom(n,k)^2 binom(n+k,k)^2,
t_(m,j)=T(m,j).
```

Split `k=jp+s`.  For the anchor `s=0`, applying `(1.2)` to both
binomial coefficients gives

```text
T(mp,jp)
 ==t_(m,j)[1-(4/3)m^2j p^3B_(p-3)]                (mod p^4).  (2.1)
```

For `0<=j<=m-1` and `1<=s<=p-1`, elementary product expansion gives

```text
(1/p)binom(mp,jp+s)
 ==(-1)^(s-1)(m-j)binom(m,j)/s
   *[1-p((m-j)H_(s-1)+jH_s)]                      (mod p^2),

binom((m+j)p+s,jp+s)
 ==binom(m+j,j)[1+mpH_s]                          (mod p^2).
```

Consequently

```text
T(mp,jp+s)
 ==p^2(m-j)^2t_(m,j)
   [1/s^2+2p(m-j)/s^3]                            (mod p^4).  (2.2)
```

Summing `(2.2)` over `s` and using `(1.1)` leaves

```text
sum_(s=1)^(p-1)T(mp,jp+s)
 ==(2/3)p^3B_(p-3)(m-j)^2t_(m,j)                 (mod p^4).  (2.3)
```

Thus

```text
b_(mp)-b_m
 ==(p^3B_(p-3)/3)
   [-4m^2 sum_j j t_(m,j)
     +2 sum_j (m-j)^2t_(m,j)]                    (mod p^4).  (2.4)
```

The required finite-sum reduction is

```text
24m^2 sum_j j t_(m,j)
-12 sum_j (m-j)^2t_(m,j)
=m^3(17b_m-b_(m-1)).                              (2.5)
```

Here is a one-line hypergeometric certificate.  Extend
`t_(m,m+1)=0`, put

```text
g_(m,j)=-4j^4(4m+3j)/(m+j)^2,

W_(m,j)=24m^2j-12(m-j)^2-17m^3
        +m^3(m-j)^2/(m+j)^2.
```

Since

```text
t_(m,j+1)/t_(m,j)
 =[(m-j)(m+j+1)/(j+1)^2]^2,
```

direct simplification gives

```text
g_(m,j+1)t_(m,j+1)-g_(m,j)t_(m,j)
 =W_(m,j)t_(m,j).                                 (2.6)
```

The boundary terms vanish.  Also

```text
((m-j)/(m+j))^2 t_(m,j)=t_(m-1,j).
```

Summing `(2.6)` is exactly `(2.5)`.  Substitution in `(2.4)` proves
`(D)`.

## 3. Reflected blocks

Put `u_(m,j)=t_(m-1,j)`.  The anchors are

```text
T(mp-1,jp)
 ==u_(m,j)[1-(4/3)m^2j p^3B_(p-3)]                (mod p^4).  (3.1)
```

For `0<=j<=m-1` and `1<=s<=p-1`,

```text
binom(mp-1,jp+s)
 ==(-1)^s binom(m-1,j)[1-mpH_s]                   (mod p^2),

(1/p)binom((m+j)p+s-1,jp+s)
 ==(m+j)binom(m+j-1,j)/s
   *[1+p((m+j)H_(s-1)-jH_s)]                     (mod p^2).
```

Hence

```text
T(mp-1,jp+s)
 ==p^2(m+j)^2u_(m,j)
   [1/s^2-2p(m+j)/s^3]                            (mod p^4),  (3.2)
```

and therefore

```text
b_(mp-1)-b_(m-1)
 ==(p^3B_(p-3)/3)
   [-4m^2 sum_j j u_(m,j)
     +2 sum_j (m+j)^2u_(m,j)]                    (mod p^4).  (3.3)
```

The second square sum here is the same as in `(2.4)`:

```text
(m+j)^2u_(m,j)=(m-j)^2t_(m,j).
```

It remains only to compare the two first moments.  The exact identity

```text
sum_j j[t_(m,j)-u_(m,j)]
 =(2m/3)(b_m+b_(m-1))                             (3.4)
```

has the equally short certificate

```text
h_(m,j)=-4j^4/[3m(m+j)^2].
```

Indeed, with `rho=((m-j)/(m+j))^2`,

```text
h_(m,j+1)t_(m,j+1)-h_(m,j)t_(m,j)
 =[j(1-rho)-(2m/3)(1+rho)]t_(m,j).
```

Summing proves `(3.4)`.  Combining `(2.4)`, `(2.5)`, `(3.3)`, and
`(3.4)` yields

```text
-4m^2 sum_j j u_(m,j)+2 sum_j(m+j)^2u_(m,j)
 =m^3(17b_(m-1)-b_m)/6,
```

which proves `(R)`.

## 4. One more order: proof of the strong rank-one law

Assume `p>=7`.  The elementary pairings `s` with `p-s` give

```text
2H_(p-1)+pH_(p-1)^(2)==0                           (mod p^4),
H_(p-1)^(3)==0                                     (mod p^2),
H_(p-1)^(4)==0                                     (mod p).   (4.1)
```

For example, the paired summand in the first line is
`p^3/s^4 (mod p^4)`, and its half-range sum vanishes because the
fourth-power sum does.  Also

```text
sum_(s=1)^(p-1) H_s^(2)/s^2
 =[(H_(p-1)^(2))^2+H_(p-1)^(4)]/2
 ==0                                               (mod p).   (4.2)
```

Keeping one more term in the unit-block logarithm from Section 1
gives the refined anchor formula

```text
binom(ap,bp)
 ==binom(a,b){
   1+bc[pH_(p-1)-(a-1)p^2H_(p-1)^(2)/2]
 }                                                 (mod p^5), (4.3)
```

where `a=b+c`.  The cubic and quartic logarithm terms vanish by
`(4.1)`, and the logarithm already has valuation at least three, so
its square does not enter.

### 4.1 Direct side

Extending the two product expansions in Section 2 by one term gives,
with `a=m-j`,

```text
T(mp,jp+s)
 ==p^2a^2t_(m,j)/s^2
   *{1+2pa/s
      +p^2[3a^2/s^2-2m^2H_s^(2)]}                 (mod p^5).
                                                                    (4.4)
```

After summing over `s`, the terms on the second line vanish by
`(4.1)`--`(4.2)`, as does the cubic-harmonic term.  Thus the entire
non-anchor block is simply

```text
p^2H_(p-1)^(2)(m-j)^2t_(m,j)                       (mod p^5). (4.5)
```

Applying `(4.3)` to the two anchor binomials gives the anchor
correction

```text
t_(m,j){
 2j(2m-j)pH_(p-1)
 -j(2m^2-2m+j)p^2H_(p-1)^(2)
}.                                                 (4.6)
```

Use the first congruence in `(4.1)` to replace
`pH_(p-1)` by `-p^2H_(p-1)^(2)/2`.  Adding `(4.5)` and `(4.6)` leaves

```text
p^2H_(p-1)^(2)
 [sum_j(m-j)^2t_(m,j)-2m^2sum_j jt_(m,j)]
 =p^2H_(p-1)^(2)E_m,                              (4.7)
```

where the last equality is exactly the telescoping identity `(2.5)`.

### 4.2 Reflected side

The analogous expansion is

```text
T(mp-1,jp+s)
 ==p^2(m+j)^2u_(m,j)/s^2
   *{1-2p(m+j)/s
      +p^2[3(m+j)^2/s^2-2m^2H_s^(2)]}             (mod p^5).
                                                                    (4.8)
```

Again, every new summed term vanishes by `(4.1)`--`(4.2)`.  The
anchors have the same correction `(4.6)`, with `t_(m,j)` replaced by
`u_(m,j)`.  Hence the reflected defect is

```text
p^2H_(p-1)^(2)
 [sum_j(m+j)^2u_(m,j)-2m^2sum_j ju_(m,j)]
 =p^2H_(p-1)^(2)F_m.                              (4.9)
```

At `m=1`, `(4.9)` says

```text
Delta_p=b_(p-1)-1
 ==p^2H_(p-1)^(2)                                 (mod p^5). (4.10)
```

Substitution of `(4.10)` into `(4.7)` and `(4.9)` proves the strong
rank-one laws `(0.2)`.

## 5. Independent computation

The script

```text
../scripts/q32_endpoint_bernoulli_rank_one_audit.py
```

checks the equivalent rank-one laws `(0.1)`, the half-harmonic
Bernoulli carrier, the strong laws `(0.2)`, and integrality of the
displayed coefficient sequences.  Its default run covers all primes
`5<=p<=997` and `1<=m<=12`:

```text
166 Bernoulli-carrier checks,
1992 direct endpoint checks,
1992 reflected endpoint checks,
1980 direct mod-p^5 rank-one checks,
1980 reflected mod-p^5 rank-one checks,
24 coefficient-integrality checks,
zero failures.
```

A second adversarial run over all primes `7<=p<=97` and every
`1<=m<=100` gives 4400 strong endpoint checks and zero failures,
including many cases with `m>=p`.

## 6. What this changes, and what it does not

The obstruction recorded in Section 5 of
`Q32_TARGET_CUBIC_BLOCK_LAW_2026-07-28.md` was too pessimistic at the
first omitted order.  The two endpoint quotients are not two free
prime-dependent digits; modulo `p` they are one Bernoulli scalar
times the fixed vectors `(E_m,F_m)`.

The mod-`p^5` target lift is now unconditional.  The re-entering
endpoint term is one scalar `Delta_p`, and direct/reflected samples
at the same prime eliminate it; see
`Q32_QUARTIC_TARGET_BERNOULLI_ELIMINATION_2026-07-29.md`.
The stronger law `(0.2)` carries the same cancellation through one
additional digit and yields the common cofactor modulo `R^5` after
the degree-five jet is included.

This does **not** yet bound the target radical.  Even an arbitrarily
high family of congruence identities can merely determine the
already-existing cofactor `b_n/R` modulo higher powers of `R`; an
independent small-height or Archimedean constraint is still required.
