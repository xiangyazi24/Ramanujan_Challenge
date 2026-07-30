# The seventh endpoint grade: symbolic rank-one reduction

Date: 2026-07-29. Owner: Codex.

## 0. Statement

Let

```text
b_n=sum_(j=0)^n binom(n,j)^2 binom(n+j,j)^2,
Delta_p=b_(p-1)-1,
H_p=b_p-5+7Delta_p.
```

For a prime `p>=11`, put

```text
X=p^2 S_2,
Y=p^4 S_4,
Z=p^5 sum_(k=1)^(p-1) H_k^(2)/k^3,
```

where `S_a=sum_(k=1)^(p-1)k^(-a)`.  The endpoint expansion is
controlled by the finite-MHS congruence

```text
X^2-5Y-2Z==0                                      (mod p^7). (H6)
```

Equivalently,

```text
S_2^2-5S_4-2p sum_k H_k^(2)/k^3==0               (mod p^3).
```

Define the weight-five moments `L_m,M_m` as in
`Q32_WEIGHT_FIVE_ENDPOINT_AND_SEXTIC_TARGET_2026-07-29.md`, and put

```text
P_m=-L_m/24,             Q_m=-M_m/24,

E_m=m^3(b_(m-1)-17b_m)/12,
F_m=m^3(17b_(m-1)-b_m)/12.
```

Then, modulo `p^7`,

```text
b_(mp)-b_m
 ==E_m Delta_p+P_m H_p+p^6 R_m w_p,

b_(mp-1)-b_(m-1)
 ==F_m Delta_p+Q_m H_p+p^6 S_m w_p,              (0.1)
```

where

```text
w_p=[b_(2p)-73+824Delta_p-(752/5)H_p]/p^6,

R_m=
 m^3[(60m^3-14m^2-51)b_m+(22m^2+3)b_(m-1)]
 ------------------------------------------------,
                    288*769

S_m=
 m^3[-(22m^2+3)b_m+(60m^3+14m^2+51)b_(m-1)]
 ------------------------------------------------,
                    288*769.                      (0.2)
```

Here `R_1=S_1=0`, `R_2=1`, and `S_2=-103/769`.  Identity `(H6)` is
proved in `Q32_H6_STAR_PROOF_2026-07-29.md`, so `(0.1)` is
unconditional for `p>=11`, `p!=769`.  The primes `7`
and `769` are genuine fixed exceptions to this normalization and may
be omitted at a global logarithmic cost `O(1)`.

Moreover, H6 gives the nonlinear relation

```text
w_p == (24*769/5)(Delta_p/p^3)^2 mod p.
```

This relation restores a fixed order-seven target carrier; see
`Q32_ORDER_SEVEN_TARGET_FIXED_LAW_2026-07-29.md`.

The purpose of this note is to give the symbolic block proof of
`(0.1)`.  Numerical verification of `(H6)` is not used as a proof:
the two-layer Faulhaber proof of `(H6)` remains a separate local
obligation.

## 1. Raw block expansion

Write

```text
t_(m,j)=binom(m,j)^2 binom(m+j,j)^2,
a=m-j,                    c=m+j.
```

For the anchor `T(mp,jp)`, expansion of the two unit factorial
ratios through `p^6` gives

```text
T(mp,jp)
 ==t_(m,j){
      1-2jm^2X
       -jm^2(2j^2+m^2-2)Y
       +2j^2m^4X^2
    }                                               (mod p^7). (1.1)
```

For the direct nonanchor block, summing `T(mp,jp+k)` over
`1<=k<p` gives

```text
a^2t_(m,j){
 X+(3a^2-3a-m^2)Y
  -m^2(1+2j)X^2-4m^3Z
}                                                  (mod p^7). (1.2)
```

For the reflected nonanchor block,

```text
c^2t_(m-1,j){
 X+(3c^2+3c-m^2)Y
  -m^2(1+2j)X^2+4m^3Z
}                                                  (mod p^7). (1.3)
```

For reference, the direct relative logarithm before summation is

```text
l_1=2a/k,
l_2=a^2/k^2-2m^2H_k^(2),
l_3=2a^3/(3k^3)+4m^2jH_k^(3),
l_4=a^4/(2k^4)-m^2(m^2+6j^2)H_k^(4).
```

After exponentiation, the relative `p^4` coefficient is

```text
5a^4/k^4
-6a^2m^2H_k^(2)/k^2
+8am^2jH_k^(3)/k
+2m^4(H_k^(2))^2
-m^2(m^2+6j^2)H_k^(4).
```

The reflected formula changes `a` to `c` and reverses the sign of
the `H_k^(3)` term.  The weight-six sums in this last line vanish
modulo `p` for `p>=11`; the surviving lower-weight terms reduce to
`(1.2)` and `(1.3)`.  In the anchor, the raw `X` coefficient is
`-2jm^2`; the often-used extra constant `8` appears only after
replacing `X` by `Delta_p`, and must not be inserted into `(1.1)`.

## 2. The residual plane is termwise

For the direct block define its four raw coefficients by

```text
x^+ =a^2-2jm^2,
y^+ =-jm^2(2j^2+m^2-2)+a^2(3a^2-3a-m^2),
u^+ =2j^2m^4-a^2m^2(1+2j),
z^+ =-4m^3a^2.                                    (2.1)
```

Thus its contribution beyond `t_(m,j)` is

```text
t_(m,j)(x^+X+y^+Y+u^+X^2+z^+Z).
```

Let

```text
K_(m,j)=-4jm^2(2j^2+m^2+8)/5,

C^+_(m,j)
 =4-12a/5+12a^2/5-4m^2/5+8m^3,

ell^+_(m,j)=K_(m,j)+a^2C^+_(m,j).                 (2.2)
```

This is exactly the summand defining `L_m`.  Direct expansion gives
the following polynomial identity, before summing in `j`:

```text
10z^+-4y^+-20x^++5ell^+ =0.                      (2.3)
```

The reflected coefficients are

```text
x^- =c^2-2jm^2,
y^- =-jm^2(2j^2+m^2-2)+c^2(3c^2+3c-m^2),
u^- =2j^2m^4-c^2m^2(1+2j),
z^- =4m^3c^2,

C^-_(m,j)
 =4+12c/5+12c^2/5-4m^2/5-8m^3,

ell^-_(m,j)=K_(m,j)+c^2C^-_(m,j),
```

and satisfy the identical relation

```text
10z^--4y^--20x^-+5ell^- =0.                      (2.4)
```

Now sum the direct coefficients and denote them by
`A_X,A_Y,A_U,A_Z`.  The already-proved cubic and weight-five laws
say

```text
A_X=E_m,                    sum_j t_(m,j)ell^+=L_m.
```

Subtracting `E_m Delta_p+P_mH_p`, with

```text
Delta_p==X+5Y-X^2+4Z,
H_p==30Y+24Z                                      (mod p^7), (2.5)
```

leaves

```text
alpha_m=A_Y-5E_m-30P_m,
beta_m =A_U+E_m,
gamma_m=A_Z-4E_m-24P_m.                           (2.6)
```

Summing `(2.3)` and using `P_m=-L_m/24` gives

```text
5gamma_m=2alpha_m.                                (2.7)
```

The same argument using `(2.4)` proves the reflected analogue.
Hence every endpoint residual lies in the fixed plane

```text
alpha Y+beta X^2+(2alpha/5)Z.                     (2.8)
```

This is the structural reason the next grade is still rank one:
`(H6)` reduces `(2.8)` to

```text
(alpha+5beta)(Y+2Z/5).                            (2.9)
```

## 3. Closed carrier and a telescoping certificate

It remains to evaluate

```text
C_m=alpha_m+5beta_m
   =sum_(j=0)^m t_(m,j)
      [y^+_(m,j)+5u^+_(m,j)+(5/4)ell^+_(m,j)].
                                                               (3.1)
```

Put

```text
q^+_(m,j)=y^+_(m,j)+5u^+_(m,j)+(5/4)ell^+_(m,j),

rho_(m,j)=
 m^3[
   60m^3-14m^2-51
   +(22m^2+3)(m-j)^2/(m+j)^2
 ]/12.
```

Let

```text
Q^+_(m,j)=
 6j^3-15j^2m^2-18j^2m-30jm^3+12jm^2+9j
 -15m^4+58m^3+12m,

G^+_(m,j)=t_(m,j) j^4Q^+_(m,j)/(3(m+j)^2).
```

Using

```text
t_(m,j+1)/t_(m,j)
 =[(m-j)(m+j+1)/(j+1)^2]^2,
```

one checks the exact Gosper certificate

```text
t_(m,j)[q^+_(m,j)-rho_(m,j)]
 =G^+_(m,j+1)-G^+_(m,j).                         (3.2)
```

Both boundary terms vanish.  Moreover,

```text
t_(m,j)(m-j)^2/(m+j)^2=t_(m-1,j).
```

Summing `(3.2)` proves

```text
C_m=
 m^3[
  (60m^3-14m^2-51)b_m+(22m^2+3)b_(m-1)
 ]/12.                                            (3.3)
```

For completeness, the reflected sum has its own certificate.  Put

```text
q^-_(m,j)=y^-_(m,j)+5u^-_(m,j)+(5/4)ell^-_(m,j),

rho^-_(m,j)=
 m^3[
   60m^3+14m^2+51
   -(22m^2+3)(m+j)^2/(m-j)^2
 ]/12,

Q^-_(m,j)=
 6j^3-15j^2m^2+18j^2m+30jm^3+12jm^2+9j
 -15m^4-58m^3-12m,

G^-_(m,j)=
 t_(m-1,j)j^4Q^-_(m,j)/(3(m-j)^2).
```

For `0<=j<m`,

```text
t_(m-1,j)[q^-_(m,j)-rho^-_(m,j)]
 =G^-_(m,j+1)-G^-_(m,j).                         (3.4)
```

The terminal value is interpreted after cancelling the zero of
`t_(m-1,m)`.  Explicitly,

```text
G^-_(m,m)
 =-m(2m-1)^2(22m^2+3)t_(m-1,m-1)/3
 =-m^3(22m^2+3)t_(m,m)/12.
```

This is exactly the missing `j=m` term when `b_m` is rewritten in
terms of `t_(m-1,j)`.  Therefore

```text
C^-_m=
 m^3[
  -(22m^2+3)b_m+(60m^3+14m^2+51)b_(m-1)
 ]/12.                                            (3.5)
```

## 4. Normalization at `m=2`

The four raw endpoint formulas at `m=1,2` are

```text
b_(p-1)-1 = X+5Y-X^2+4Z,
b_p-5     =-7X-5Y+7X^2-4Z,

b_(2p)-73 =-824X-3592Y+5312X^2-1280Z,
b_(2p-1)-5=   8X+1144Y-320X^2+1280Z              (mod p^7).
```

Consequently,

```text
D_2=b_(2p)-73+824Delta_p-(752/5)H_p
    =(24/5)(935X^2-830Y-332Z),

F_2=b_(2p-1)-5-8Delta_p-(336/5)H_p
    =-(24/5)(65X^2+190Y+76Z).
```

The exact elimination is

```text
769F_2+103D_2
 =222336(X^2-5Y-2Z).                              (4.1)
```

Also,

```text
alpha_2+5beta_2=18456=24*769.
```

Combining `(2.9)`, `(3.3)`, and `(3.5)` with this normalization
gives precisely `(0.1)`--`(0.2)`.

## 5. Proof boundary

The endpoint block calculation and both all-`m` carrier evaluations
above are exact polynomial identities.  The remaining local proof
boundary is only `(H6)`.  Direct computation currently gives:

```text
120/120 primes 11<=p<=683: pass,
p=7: residue 294 mod 343.
```

The reversal-plus-stuffle equations are linearly dependent modulo
`p^2` and cannot prove `(H6)`.  The viable proof is the direct
two-layer Faulhaber expansion of the inclusive double sum, as
recorded in Fable mailbox r31.
