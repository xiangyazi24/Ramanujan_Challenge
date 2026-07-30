# Bernoulli-free quartic and quintic target laws for the Apéry cofactor

## 0. The fourth common digit

Let `p>=7` be prime, let `0<=r<p`, and assume `p|b_r`.  Put

```text
n=p+r,                  s=p-1-r,
x=b_n/p,
e_p=(b_(p-1)-1)/p^3     (mod p).
```

Then `p|b_s` as well.  The exact block decompositions and the
rank-one endpoint theorem give two residues `D_(p,r),Z_(p,s)` modulo
`p^4`, defined below, such that

```text
D_(p,r)==x(1+(7/5)p^3e_p)                       (mod p^4),
Z_(p,s)==x(1-(8/5)p^3e_p)                       (mod p^4).    (0.1)
```

Therefore

```text
15x==8D_(p,r)+7Z_(p,s)                          (mod p^4).    (0.2)
```

The Bernoulli scalar has disappeared.

For one fixed upper index `n`, let `S` be any set of target primes in
`(n/2,n]`, and put

```text
R=product_(p in S)p,              A=b_n/R.
```

Writing `R_p=R/p`, `(0.2)` gives

```text
A
 ==[15R_p]^(-1)
   [8D_(p,n-p)+7Z_(p,2p-1-n)]                    (mod p^4).   (0.3)
```

All right sides are reductions of the same integer `A`.  CRT thus
determines

```text
A modulo R^4.                                                   (0.4)
```

This is the fourth consecutive common quotient digit:

```text
Gessel's linear jet              A modulo R,
the quadratic prime-cube jet     A modulo R^2,
the target-cubic law             A modulo R^3,
Bernoulli elimination            A modulo R^4.
```

The stronger endpoint theorem modulo `p^5` gives one more digit with
the **same weights**.  If

```text
Delta_p=b_(p-1)-1,
```

then the degree-five residues `D5_(p,r),Z5_(p,s)` defined in Section
4 satisfy

```text
D5_(p,r)==x(1+(7/5)Delta_p)                         (mod p^5),
Z5_(p,s)==x(1-(8/5)Delta_p)                         (mod p^5),

15x==8D5_(p,r)+7Z5_(p,s)                            (mod p^5). (0.5)
```

Consequently the same-index construction determines

```text
A modulo R^5.                                                   (0.6)
```

## 1. Shifted fundamental solutions

Let `mathcal U_r(X),mathcal V_r(X)` be defined by

```text
U_0=1, U_1=B(X)/(1+X)^3,
V_0=0, V_1=1,

(X+r+1)^3W_(r+1)
 =B(X+r)W_r-(X+r)^3W_(r-1),
```

where `B(X)=34X^3+51X^2+27X+5`.  Define

```text
T_r^[4](X)=mathcal U_r(X) mod X^5
          =b_r+XG_r+X^2H_r+X^3K_r+X^4L_r,

mathcal J_r(X)=mathcal V_r(X)/(1+X)^3,
A_r=mathcal J_r(0)=a_r/6,
C_r=mathcal J_r'(0).
```

For `r<p`, all these Taylor coefficients are `p`-integral.  The exact
identities from the cubic-law note are

```text
b_(M+r)
 =b_M mathcal U_r(M)-M^3b_(M-1)mathcal J_r(M),      (1.1)

b_(M-1-r)
 =b_(M-1)mathcal U_r(-M)
  -(-M)^3b_M mathcal J_r(-M).                       (1.2)
```

## 2. Direct expansion

The endpoint theorem gives

```text
b_p    ==5-7p^3e_p                                 (mod p^4),
b_(p-1)==1+p^3e_p                                  (mod p^4). (2.1)
```

Apply `(1.1)` with `M=p`.  Because the second term already has a
factor `p^3`, the correction in `b_(p-1)` first occurs at `p^6`.
Taylor expansion gives, modulo `p^5`,

```text
b_(p+r)
 ==5T_r^[4](p)
   -p^3(A_r+pC_r)
   -7p^3e_p(b_r+pG_r).                              (2.2)
```

Define the integral `p`-adic residue

```text
D_(p,r)
 =[5T_r^[4](p)-p^3(A_r+pC_r)]/p                    (mod p^4).
                                                                    (2.3)
```

Writing `b_r=p beta`, the prime-square direct jet says

```text
x==5(beta+G_r)                                     (mod p).
```

Divide `(2.2)` by `p` and substitute this last congruence only inside
the term already multiplied by `p^3`.  This proves the first formula
in `(0.1)`.

## 3. Reflected expansion of the same integer

Since `n=2p-1-s`, apply `(1.2)` with `M=2p`.  The endpoint theorem
gives

```text
b_(2p-1)==5+8p^3e_p                                (mod p^4),
b_(2p)  ==73-824p^3e_p                             (mod p^4). (3.1)
```

The defect in `b_(2p)` is again multiplied by an exterior `p^3` and
does not enter modulo `p^5`.  Thus

```text
b_(2p-1-s)
 ==5T_s^[4](-2p)
   +584p^3(A_s-2pC_s)
   +8p^3e_p(b_s-2pG_s)                             (mod p^5). (3.2)
```

Define

```text
Z_(p,s)
 =[5T_s^[4](-2p)+584p^3(A_s-2pC_s)]/p              (mod p^4).
                                                                    (3.3)
```

Writing `b_s=p gamma`, the reflected prime-square jet gives

```text
x==5(gamma-2G_s)                                   (mod p).
```

Division of `(3.2)` by `p` proves the second formula in `(0.1)`.
Multiplying the direct equation by `8`, the reflected equation by
`7`, and adding cancels

```text
8*7+7*(-8)=0,
```

which proves `(0.2)`--`(0.4)`.

## 4. The fifth common digit

Let

```text
T_r^[5](X)=mathcal U_r(X) mod X^6,
Q_r=[X^2]mathcal J_r(X).
```

The endpoint rank-one theorem one order deeper states, for `p>=7`,

```text
b_(mp)-b_m==E_m Delta_p                         (mod p^5),
b_(mp-1)-b_(m-1)==F_m Delta_p                  (mod p^5),

E_m=m^3(b_(m-1)-17b_m)/12,
F_m=m^3(17b_(m-1)-b_m)/12.                        (4.1)
```

In particular,

```text
b_p==5-7Delta_p,
b_(2p-1)==5+8Delta_p,
b_(2p)==73-824Delta_p                              (mod p^5). (4.2)
```

The exact direct decomposition, reduced modulo `p^6`, is

```text
b_(p+r)
 ==5T_r^[5](p)-p^3(A_r+pC_r+p^2Q_r)
   -7Delta_p mathcal U_r(p)                        (mod p^6). (4.3)
```

The product `p^3 Delta_p mathcal J_r(p)` has valuation at least six
and is absent.  Define

```text
D5_(p,r)
 =[5T_r^[5](p)-p^3(A_r+pC_r+p^2Q_r)]/p            (mod p^5).
                                                                    (4.4)
```

Since `p|mathcal U_r(p)` on a target, division of `(4.3)` is valid.
Moreover the prime-square jet gives

```text
mathcal U_r(p)/p==x/5                              (mod p^2).
```

As `p^3|Delta_p`, this proves

```text
D5_(p,r)==x(1+(7/5)Delta_p)                        (mod p^5). (4.5)
```

The reflected decomposition gives, modulo `p^6`,

```text
b_(2p-1-s)
 ==5T_s^[5](-2p)
   +584p^3(A_s-2pC_s+4p^2Q_s)
   +8Delta_p mathcal U_s(-2p).                     (4.6)
```

Here the `-824Delta_p` correction in `b_(2p)` is multiplied by
`8p^3`, so it first contributes modulo `p^6`.  Define

```text
Z5_(p,s)
 =[5T_s^[5](-2p)
   +584p^3(A_s-2pC_s+4p^2Q_s)]/p                  (mod p^5).
                                                                    (4.7)
```

The reflected prime-square jet gives

```text
mathcal U_s(-2p)/p==x/5                            (mod p^2),
```

and hence

```text
Z5_(p,s)==x(1-(8/5)Delta_p)                        (mod p^5). (4.8)
```

The same integer weights as before cancel the full endpoint defect:

```text
8*7+7*(-8)=0.
```

Thus `(0.5)` follows.  Replacing `x` by `A(R/p)` and applying CRT
over the target primes proves `(0.6)`.

## 5. Audit

The independent reproducer is

```text
../scripts/q32_quartic_target_bernoulli_elimination_audit.py.
```

It constructs both shifted fundamental solutions directly in
`Q[X]/(X^6)` and checks:

1. the direct mod-`p^5` formula `(2.2)`;
2. the reflected formula `(3.2)` for the same integer;
3. both normalized equations `(0.1)`;
4. the Bernoulli-free combination `(0.2)`;
5. recovery of the actual common cofactor in `(0.3)`;
6. the two mod-`p^6` block formulas `(4.3)` and `(4.6)`;
7. the rank-one elimination `(0.5)`;
8. recovery of the actual common cofactor modulo `p^5`.

For all primes `7<=p<=2000` and every target digit:

```text
281 direct quartic and 281 direct quintic checks,
281 reflected quartic and 281 reflected quintic checks,
281 eliminations at each of the two precisions,
281 common-cofactor checks at each of mod-p^4 and mod-p^5,
zero failures.
```

## 6. Exact ceiling

This is positive structure, but fixed order five is still below the
global height threshold.  At the maximal top primorial,

```text
5log R<=2.5n+o(n),
log A>=3.02549...n+o(n).
```

For a hypothetical spike `log R>=epsilon n`, reconstruction of the
real integer `A` would require order exceeding
`(3.52549-epsilon)/epsilon`; no fixed order suffices uniformly as
`epsilon` tends to zero.

The important change is local: the first prime-dependent Frobenius
scalar can be removed through two consecutive quotient digits by two
descriptions of the **same** integer.  The next question is now
exact.  To obtain the next block congruence modulo `p^7`, determine
the rank of the endpoint defect modulo `p^6` and whether sufficiently
many reflected block descriptions of the same `b_n` annihilate it
without importing linear-height coefficients.
