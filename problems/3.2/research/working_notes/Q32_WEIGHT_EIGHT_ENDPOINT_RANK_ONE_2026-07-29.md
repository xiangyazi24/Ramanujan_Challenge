# Computational weight-eight endpoint rank-one law

Date: 2026-07-29. Owner: Codex.

## 0. Status

This note records an exact computational discovery.  The displayed
all-`m` law is **not yet a symbolic theorem**.  It was reconstructed
from exact modular data through `p<=400`, then independently checked
through `p<=1000` and `m<=20` with zero failures.

The result is strong evidence that the endpoint tower remains rank one
after the proved precision-seven layer, even though the new endpoint
digit is not a constant multiple of the complete-block scalar
`S_6/p`.

## 1. Coordinates

For `p>=11`, `p!=769`, put

```text
Delta=b_(p-1)-1,
H=b_p-5+7Delta,

w=[b_(2p)-73+824Delta-(752/5)H]/p^6,

v=[b_(2p-1)-5-8Delta-(336/5)H
   +(103/769)p^6w]/p^7.
```

The proved endpoint laws imply `w,v in Z_(p)`.  Let the already-proved
lower-grade coefficients be denoted by

```text
E_m=m^3(b_(m-1)-17b_m)/12,
F_m=m^3(17b_(m-1)-b_m)/12,

P_m=-L_m/24,                  Q_m=-M_m/24,

R_m,S_m
```

with `L_m,M_m,R_m,S_m` as in
`Q32_WEIGHT_SEVEN_ENDPOINT_RANK_ONE_2026-07-29.md`.

## 2. Reconstructed law

Define

```text
P0(m)=3845m^4-29268m^3+36974m^2-9112,
Q0(m)=45371m^4-58102m^2+536,

N=305911296=2^9*3^3*22129,

C_m=-m^3[P0(m)b_m+Q0(m)b_(m-1)]/N,

D_m= m^3[Q0(m)b_m+P0(-m)b_(m-1)]/N.
```

The conjectural endpoint law is

```text
b_(mp)-b_m
 ==E_m Delta+P_m H+p^6R_m w+p^7C_m v              (mod p^8), (2.1)

b_(mp-1)-b_(m-1)
 ==F_m Delta+Q_m H+p^6S_m w+p^7D_m v              (mod p^8). (2.2)
```

The normalization is exact:

```text
C_1=D_1=C_2=0,                    D_2=1.
```

Thus `v` is precisely the reflected `m=2` residual, while every
tested endpoint residual is its fixed rational multiple.

The new fixed denominator prime is

```text
22129.
```

It is an exception to this normalization, just as `769` was at the
preceding endpoint grade.  It is not needed for the separate
precision-eight target Casoratian law, which uses `v` directly.

## 3. Reconstruction and hostile checks

For each prime, the script retains the full `p`-adic digits of `w`
rather than replacing it by its residue modulo `p`.  It forms

```text
d_m=[
 b_(mp)-b_m-E_mDelta-P_mH-p^6R_mw
]/p^7                                                   (mod p),

f_m=[
 b_(mp-1)-b_(m-1)-F_mDelta-Q_mH-p^6S_mw
]/p^7                                                   (mod p).
```

The anchor is `f_2=v`.  Chinese-remainder rational reconstruction
over the initial prime set gives the coefficients `C_m,D_m`.
Solving the symmetry ansatz

```text
C_m/m^3=A(m)b_m+B(m)b_(m-1),

D_m/m^3=-B(m)b_m-A(-m)b_(m-1)
```

finds degree four, with

```text
A(m)=-P0(m)/N,                 B(m)=-Q0(m)/N.
```

The formula then passes disjoint larger primes and two additional
quotient rows.  The final exact audit through `p<=1000,m<=20` gives

```text
direct rank-one checks                  3260
reflected rank-one checks               3260
residual divisibility checks            6520
zero anchor primes                         0
failures                                   0
```

This excludes a numerical accident of a few primes, but it is still
finite evidence.

## 4. A proved reduction of the primitive MHS quotient

The lifted-reversal gap in the first version of this note can be closed.
If

```text
xi=S_6/p, eta=H(2,4)/p, A=H(2,2,3), B=H(2,5)  (mod p),
```

then, for every `p>=11`,

```text
3 eta=2 xi,             3 A=14 xi,             2 B=-7 xi.
```

The proof uses the depth-two Bernoulli formula at weight seven,
first-order lifted reversal modulo `p^2`, and one stuffle identity.
In particular,

```text
H(2,4)-H(4,2)
 ==p[2H(4,3)+4H(5,2)]                    (mod p^2),
```

so the quotient correction omitted by ordinary reversal is explicit.
The full proof and an exact audit through `p<=1000` are

```text
Q32_WEIGHT_SEVEN_MHS_RANK_ONE_2026-07-29.md
../scripts/q32_weight_seven_mhs_rank_one_audit.py
```

Thus no new finite-MHS dimension obstructs the conjectural law.  This
does not yet prove it: changing from the primitive MHS basis to
`(Delta,H,w,v)` also contains lifted lower-coordinate terms.  A naive
scalar projection of the primitive block vector omits those terms.

## 5. Remaining exact proof obligation

A proof must extend the block calculation by one effective grade:

1. expand the anchor and nonanchor blocks through weight seven;
2. retain the lifted lower-grade coordinates together with the now
   one-dimensional weight-seven finite-MHS quotient;
3. change from that block basis to `(Delta,H,w,v)`;
4. prove the resulting termwise direct/reflected identities;
5. telescope their two sums to `C_m,D_m`.

The likely final step is a termwise plane identity followed by a
Gosper certificate, as at the preceding grade.  Until those identities
are printed and checked symbolically, `(2.1)--(2.2)` remains a
high-confidence conjectural law.

## 6. Reproducer

The exact script is

```text
../scripts/q32_weight_eight_endpoint_rank_audit.py
```

Run

```text
python3 q32_weight_eight_endpoint_rank_audit.py \
  --prime-limit 1000 --quotient-limit 20
```

It performs no floating-point calculation.
