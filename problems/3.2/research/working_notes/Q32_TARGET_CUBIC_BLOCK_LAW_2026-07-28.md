# Target-cubic block law and the common quotient modulo `R^3`

Date: 2026-07-28. Owner: Codex.

## 0. Verdict

The universal quadratic jet has one further target-selective lift.
The characteristic-dependent endpoint error at order `p^3` propagates
in the two-dimensional solution space.  On a target `p|b_r`, its
component along the first Apéry solution vanishes modulo `p`; the
remaining component is an explicit multiple of the second Apéry
solution `a_r`.

Let

```text
(r+1)^3u_(r+1)=B(r)u_r-r^3u_(r-1),
B(r)=34r^3+51r^2+27r+5,

b_0=1, b_1=5,
a_0=0, a_1=6.
```

Define the universal cubic jet

```text
T_r(X)=b_r+XG_r+X^2H_r+X^3K_r             (mod X^4)          (0.1)
```

by the shifted recurrence in Section 1.  If `p>=5` is prime,
`0<=r<=p-1`, and `p|b_r`, then for every `m>=1`,

```text
b_(mp+r)
 ==b_m T_r(mp)-m^3p^3 b_(m-1)a_r/6         (mod p^4),        (D4)
```

and for every `m>=0`,

```text
b_((m+1)p-1-r)
 ==b_m T_r(-(m+1)p)
   +(m+1)^3p^3 b_(m+1)a_r/6                (mod p^4).        (R4)
```

In the top block these become

```text
b_(p+r)
 ==5T_r(p)-p^3a_r/6                         (mod p^4),        (0.2)

b_(p-1-r)
 ==T_r(-p)+5p^3a_r/6                        (mod p^4),        (0.3)

b_(2p-1-r)
 ==5T_r(-2p)+(292/3)p^3a_r                  (mod p^4).        (0.4)
```

The coefficient in (0.4) is `2^3b_2/6=8*73/6=292/3`.
Every rational number displayed is `p`-integral because `r<p`.

This closes the "second-layer collapse" question in a stronger form:
the first endpoint layer not already contained in the universal jet
introduces no new sequence.  It is exactly the classical second Apéry
solution.

For a fixed upper index `n`, let `S` be any set of top targets, put

```text
R=product_(p in S)p,       A=b_n/R.
```

Then the local formulas determine the same integer `A` modulo `p^3`
for every `p in S`, and hence determine `A modulo R^3`.  This is a
one-digit upgrade of the preceding prime-cube quadratic-jet law,
which determined `A modulo R^2`; Gessel's original prime-square
first-jet law determined `A modulo R`.

It still does not prove `log R=o(n)`.  Even at the maximal possible
top primorial,

```text
log A>=3.02549...n+o(n),       3log R<=1.5n+o(n).
```

Thus `R^3` is exponentially smaller than the actual quotient `A`.
For a forbidden but arbitrarily small spike `log R>=epsilon n`, any
CRT reconstruction argument would require an order tending past
`3.52549/epsilon`, not the fixed order three proved here.

## 1. Two exact shifted fundamental solutions

Let `mathcal U_r(X),mathcal V_r(X)` be the rational functions defined
by

```text
mathcal U_0(X)=1,
mathcal U_1(X)=B(X)/(X+1)^3,

mathcal V_0(X)=0,
mathcal V_1(X)=1,                                      (1.1)
```

and, for either `mathcal W=mathcal U` or `mathcal V`,

```text
(X+r+1)^3 mathcal W_(r+1)(X)
 =B(X+r)mathcal W_r(X)-(X+r)^3mathcal W_(r-1)(X).
                                                               (1.2)
```

For `0<=r<=p-1` and `X` divisible by `p`, every denominator in
(1.2) is a `p`-unit.  Define `T_r` to be the Taylor polynomial

```text
T_r(X)=mathcal U_r(X) mod X^4.                       (1.3)
```

Since

```text
B(X)/(X+1)^3=5+12X-7X^3+O(X^4),
```

the initial cubic jets are

```text
T_0(X)=1,       T_1(X)=5+12X-7X^3.
```

The constant, linear, and quadratic coefficients are the already
proved sequences `b_r,G_r,H_r`.  At `X=0`, the second fundamental
solution has initial values `(0,1)`, whereas `a_r` has initial values
`(0,6)`.  Therefore

```text
mathcal V_r(0)=a_r/6.                                (1.4)
```

## 2. Exact direct block decomposition

Fix an integer `M>=1`.  As a sequence in `r`, `b_(M+r)` obeys
(1.2) at `X=M`.  Its values at `r=0,1` give the exact identity

```text
b_(M+r)
 =b_M mathcal U_r(M)
  -M^3 b_(M-1)/(M+1)^3 mathcal V_r(M).              (2.1)
```

Indeed, (2.1) is clear at `r=0`; at `r=1` it is exactly the Apéry
recurrence at index `M`; thereafter the second-order recurrence
propagates it.

Take `M=mp`.  The endpoint supercongruences give

```text
b_(mp)==b_m,       b_(mp-1)==b_(m-1)                (mod p^3).
                                                               (2.2)
```

Subtract `b_m mathcal U_r(mp)` from (2.1).  The first endpoint error
is

```text
(b_(mp)-b_m)mathcal U_r(mp).
```

It is divisible by `p^3`; on a target,

```text
mathcal U_r(mp)==b_r==0                              (mod p),
```

so it is divisible by `p^4`.  In the second term of (2.1), reduction
after division by `p^3` uses (1.4) and (2.2):

```text
-m^3p^3 b_(mp-1)/(mp+1)^3 mathcal V_r(mp)
 ==-m^3p^3 b_(m-1)a_r/6                            (mod p^4).
```

Finally `mathcal U_r(mp)==T_r(mp) mod p^4`.  This proves `(D4)`.

## 3. Exact reflected block decomposition

Put

```text
M=(m+1)p,       X=-M,       c_r=b_(M-1-r).
```

Using `B(-x-1)=-B(x)`, the reflected sequence `c_r` again obeys
(1.2), now at `X=-M`.  Its two endpoint values give

```text
b_(M-1-r)
 =b_(M-1)mathcal U_r(-M)
  -(-M)^3 b_M/(1-M)^3 mathcal V_r(-M).              (3.1)
```

The reflected endpoint supercongruences are

```text
b_((m+1)p-1)==b_m,
b_((m+1)p)==b_(m+1)                                 (mod p^3).
                                                               (3.2)
```

As on the direct side, the first endpoint error in (3.1) acquires a
fourth power of `p` from the target condition.  Since

```text
-(-M)^3=(m+1)^3p^3,       (1-M)^3==1 (mod p),
```

the second term reduces to

```text
(m+1)^3p^3 b_(m+1)a_r/6                            (mod p^4).
```

Replacing `mathcal U_r(-M)` by `T_r(-M)` proves `(R4)`.

## 4. The common quotient modulo `R^3`

Fix `n`, let `S` be a set of primes `p` with

```text
n/2<p<=n,       p|b_n,
```

and put `r_p=n-p`.  Lucas gives `p|b_(r_p)`.  Write

```text
R_p=R/p.
```

Dividing (0.2) by `p` gives

```text
A R_p
 ==5[b_(r_p)/p+G_(r_p)+pH_(r_p)+p^2K_(r_p)]
   -p^2a_(r_p)/6                                  (mod p^3).  (4.1)
```

Every `R_p` is a unit modulo `p^3`.  Hence

```text
A
 ==R_p^(-1){
      5[b_(r_p)/p+G_(r_p)+pH_(r_p)+p^2K_(r_p)]
      -p^2a_(r_p)/6
   }                                               (mod p^3).  (4.2)
```

The right sides for all targets are therefore the reductions of one
integer modulo pairwise coprime moduli `p^3`.  CRT identifies one
class modulo

```text
product_(p in S)p^3=R^3.                           (4.3)
```

This is an actual cross-characteristic relation among the divided
digits.  It is not enough for reconstruction because its canonical
CRT representative can have height `3log R`, while the known integer
`A` has height about `3.52549n-log R`.

## 5. Exact ceiling and next gate

If `log R>=epsilon n`, then

```text
log A=(3.52549...n)-log R+O(log n),
log R^3>=3epsilon n.
```

For small fixed `epsilon`, the modulus still does not identify `A`.
To prove `log R=o(n)` by this mechanism one needs at least one of:

```text
1. target-selective block laws of arbitrarily high order, with only
   finitely many globally controlled solution coordinates;

2. a representative of the class (4.2) whose height is o(n), rather
   than the generic CRT height 3log R;

3. a nonzero low-height polynomial satisfied by A modulo R^3.
```

The first two omitted endpoint digits are now controlled exactly.
The theorem proved in
`Q32_ENDPOINT_BERNOULLI_RANK_ONE_2026-07-29.md` gives, first,

```text
b_(mp)-b_m
 ==-m^3(17b_m-b_(m-1))p^3B_(p-3)/18          (mod p^4),

b_(mp-1)-b_(m-1)
 == m^3(17b_(m-1)-b_m)p^3B_(p-3)/18          (mod p^4).
```

Thus these are not two free prime-dependent digits: they are one
Bernoulli scalar times a fixed two-vector in `m`.  More strongly, for
`p>=7`, with `Delta_p=b_(p-1)-1`,

```text
b_(mp)-b_m
 ==m^3(b_(m-1)-17b_m)Delta_p/12                (mod p^5),

b_(mp-1)-b_(m-1)
 ==m^3(17b_(m-1)-b_m)Delta_p/12                (mod p^5).
```

The same-index direct and reflected formulas eliminate `Delta_p`
with the fixed weights `8,7` through one further digit.  The complete
derivation in

```text
Q32_QUARTIC_TARGET_BERNOULLI_ELIMINATION_2026-07-29.md
```

therefore determines the common cofactor `A` modulo `R^5`, not merely
modulo `R^3`.

The global warning remains: even a completely controlled congruence
tower only identifies the already-existing integer `A` modulo larger
powers of `R`.  To bound `R`, the tower still needs an independent
small-height, non-stabilization, positivity, or Archimedean input.

## 6. Reproducer

The script

```text
../scripts/q32_target_cubic_block_law_audit.py
```

checks:

1. the exact recurrence construction of `b_r,a_r,G_r,H_r,K_r`;
2. both endpoint supercongruences;
3. `(D4)` for all requested quotients and every target digit;
4. `(R4)` for all requested quotients and every target digit;
5. all three top-block specializations (0.2)--(0.4);
6. the common-quotient congruence (4.2).

All congruences involving rational numbers are evaluated only after
asserting that their denominators are units modulo the relevant
prime power.
