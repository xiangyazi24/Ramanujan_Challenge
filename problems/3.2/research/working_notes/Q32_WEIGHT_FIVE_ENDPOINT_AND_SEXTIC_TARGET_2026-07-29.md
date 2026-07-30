# Weight-five endpoint rank one and the sextic same-index target law

Date: 2026-07-29. Owner: Codex.

## 0. Results

Let `b_n` be the Apéry numbers for `zeta(3)`.  Put

```text
Delta_p=b_(p-1)-1.
```

For every prime `p>=7` and every `m>=1`, the endpoint defects one
order beyond the `Delta_p` law still have rank one:

```text
b_(mp)-b_m
 ==E_m Delta_p+p^5B_(p-5)L_m                    (mod p^6),

b_(mp-1)-b_(m-1)
 ==F_m Delta_p+p^5B_(p-5)M_m                   (mod p^6),    (0.1)
```

where

```text
E_m=m^3(b_(m-1)-17b_m)/12,
F_m=m^3(17b_(m-1)-b_m)/12,
```

and the prime-independent coefficients `L_m,M_m` are the explicit
finite sums in Section 3.  In particular,

```text
L_1=-24,          M_1=0,
L_2=-18048/5,     M_2=-8064/5.                  (0.2)
```

Thus, if

```text
h_p=(b_p-5+7Delta_p)/p^5                         (mod p),
```

then

```text
h_p==-24B_(p-5),

(b_(2p-1)-5-8Delta_p)/p^5
 ==(336/5)h_p                                    (mod p).    (0.3)
```

Now let `0<=r<p`, assume `p|b_r`, and put

```text
n=p+r,        s=p-1-r,        x=b_n/p.
```

The direct and reflected descriptions of the same `b_n` give
explicit residues `D6_(p,r),Z6_(p,s)` modulo `p^6` such that

```text
D6_(p,r)==x(1-p^5h_p/5),
Z6_(p,s)==x(1-336p^5h_p/25)                      (mod p^6).
```

Consequently

```text
331x==336D6_(p,r)-5Z6_(p,s)                      (mod p^6).  (0.4)
```

For every target prime other than `331`, this determines `x mod
p^6`.  For a fixed upper index, CRT therefore determines the common
cofactor modulo the sixth power of the product of all target primes
other than `331`.  Omitting that one fixed prime changes the relevant
logarithm by at most `log 331`.

## 1. Harmonic lemma at weight five

Write

```text
S_a=sum_(k=1)^(p-1) 1/k^a,
H_k^(a)=sum_(u=1)^k 1/u^a,
B=B_(p-5).
```

For `p>=7`, the following congruences hold:

```text
S_3==-(6/5)p^2B                                  (mod p^3),
S_4== (4/5)pB                                    (mod p^2), (1.1)

sum_k H_k^(2)/k^2 ==(2/5)pB                      (mod p^2),
sum_k H_k^(2)/k^3 ==-2B                          (mod p),
sum_k H_k^(3)/k^2 == 2B                          (mod p).    (1.2)
```

The first line is the standard one-step refinement of the finite
power-sum congruence:

```text
S_a==a p B_(p-1-a)/(a+1)                         (mod p^2)
```

for even `a`, and

```text
S_a==-a(a+1)p^2B_(p-2-a)/(2(a+2))                (mod p^3)
```

for odd `a`.  Only `a=4` and `a=3` are used here.

For the first identity in `(1.2)`, stuffle gives the exact equality

```text
sum_k H_k^(2)/k^2=(S_2^2+S_4)/2,
```

and `p|S_2`.  For the last two, write strict double harmonic sums.
The elementary Faulhaber evaluation

```text
sum_(1<=u<k<=p-1) 1/(u^a k^b)
 ==(-1)^b binom(a+b,a)B_(p-a-b)/(a+b)            (mod p)
```

gives `-2B` at `(a,b)=(2,3)` and `2B` at `(3,2)`;
the diagonal `S_5` vanishes modulo `p`.

The exact product for the reflected endpoint is

```text
T(p-1,k)
 =p^2/k^2 (1-p/k)^2
  product_(u<k)(1-p^2/u^2)^2,
```

where `T(N,k)=binom(N,k)^2binom(N+k,k)^2`.  Expanding through
relative order `p^3` gives

```text
Delta_p
 ==p^2S_2-2p^3S_3
   +p^4[S_4-2sum_k H_(k-1)^(2)/k^2]
   +4p^5sum_k H_(k-1)^(2)/k^3                  (mod p^6).
```

Use `(1.1)`--`(1.2)` and

```text
sum_k H_(k-1)^(2)/k^2=(S_2^2-S_4)/2
```

to obtain the key scalar identity

```text
p^2S_2==Delta_p+4p^5B                            (mod p^6). (1.3)
```

Finally, pairing `k` with `p-k` yields

```text
pS_1==-Delta_p/2-(9/5)p^5B                       (mod p^6). (1.4)
```

Indeed,

```text
2pS_1==-p^2S_2-p^3S_3-p^4S_4                    (mod p^6),
```

and `(1.1)`, `(1.3)` finish the calculation.

## 2. Anchor and nonanchor blocks modulo `p^6`

Put

```text
t_(m,j)=binom(m,j)^2binom(m+j,j)^2,
a=m-j,       c=m+j.
```

For `A=B+C`, remove the multiples of `p` from the factorials and set

```text
D_q(A,B)=sum_(u=0)^(A-1)u^q
          -sum_(u=0)^(B-1)u^q
          -sum_(u=0)^(C-1)u^q.
```

The unit-block logarithm through `p^5` is

```text
pD_1S_1-(p^2/2)D_2S_2+(p^3/3)D_3S_3
 -(p^4/4)D_4S_4.                                  (2.1)
```

The square of `(2.1)` has valuation at least six.  Substitution of
Section 1 shows that the part beyond the `Delta_p` correction is

```text
-(p^5B/5)(9D_1+10D_2+2D_3+D_4).                  (2.2)
```

Apply `(2.2)` to the two binomial coefficients in the anchor
`T(mp,jp)`, square them, and multiply.  The result is

```text
T(mp,jp)
 ==t_(m,j){
    1-2m^2j Delta_p+p^5B K_(m,j)
   }                                               (mod p^6), (2.3)

K_(m,j)=-4jm^2(2j^2+m^2+8)/5.
```

For `1<=k<=p-1`, the direct nonanchor logarithm from the exact
factorial ratios gives

```text
T(mp,jp+k)
 ==p^2a^2t_(m,j)/k^2
   {1+2pa/k
     +p^2[3a^2/k^2-2m^2H_k^(2)]
     +p^3[4a^3/k^3
           -4am^2H_k^(2)/k+4m^2jH_k^(3)]}
                                                    (mod p^6). (2.4)
```

After summing over `k` and applying Section 1, `(2.4)` becomes

```text
a^2t_(m,j){
 Delta_p+p^5B C^+_(m,j)
},

C^+_(m,j)
 =4-12a/5+12a^2/5-4m^2/5+8m^3.                   (2.5)
```

For the reflected endpoint use the exact factor

```text
T(mp-1,jp+k)
 =[(ap-k)/(cp+k)]^2 T(mp,jp+k).
```

It changes the bracket in `(2.4)` to

```text
1-2pc/k
 +p^2[3c^2/k^2-2m^2H_k^(2)]
 +p^3[-4c^3/k^3
       +4cm^2H_k^(2)/k+4m^2jH_k^(3)].
```

Thus the reflected nonanchor sum is

```text
c^2t_(m-1,j){
 Delta_p+p^5B C^-_(m,j)
},

C^-_(m,j)
 =4+12c/5+12c^2/5-4m^2/5-8m^3.                  (2.6)
```

## 3. The endpoint coefficient sequences

Define

```text
P^+_(m,j)=K_(m,j)+(m-j)^2 C^+_(m,j),
P^-_(m,j)=K_(m,j)+(m+j)^2 C^-_(m,j),

L_m=sum_(j=0)^m t_(m,j)P^+_(m,j),
M_m=sum_(j=0)^(m-1) t_(m-1,j)P^-_(m,j).           (3.1)
```

Adding anchors and nonanchors in `(2.3)`--`(2.6)` gives

```text
b_(mp)-b_m
 ==Delta_p sum_j[(m-j)^2-2m^2j]t_(m,j)
   +p^5B L_m                                      (mod p^6),

b_(mp-1)-b_(m-1)
 ==Delta_p{
    sum_j(m+j)^2t_(m-1,j)-2m^2sum_j jt_(m-1,j)
   }
   +p^5B M_m                                      (mod p^6).
```

The two `Delta_p` sums are exactly `E_m,F_m`, by the telescoping
certificates in the rank-one endpoint note.  This proves `(0.1)`.
Direct substitution in `(3.1)` gives `(0.2)` and hence `(0.3)`.

## 4. Same-index elimination through the sixth quotient digit

Let `mathcal U_r(X),mathcal V_r(X)` be the shifted fundamental
solutions from the preceding target-law notes, and put

```text
mathcal J_r(X)=mathcal V_r(X)/(1+X)^3,
T_r^[6](X)=mathcal U_r(X) mod X^7,
J_r^[3](X)=mathcal J_r(X) mod X^4.
```

The exact direct decomposition is

```text
b_(p+r)=b_p mathcal U_r(p)-p^3b_(p-1)mathcal J_r(p).
```

Using

```text
b_p==5-7Delta_p+p^5h_p                            (mod p^6)
```

and the target divisibility `p|mathcal U_r(p)`, define

```text
D6_(p,r)
 =[(5-7Delta_p)T_r^[6](p)
   -p^3(1+Delta_p)J_r^[3](p)]/p                   (mod p^6).
                                                                    (4.1)
```

Since `mathcal U_r(p)/p==x/5 (mod p)`,

```text
D6_(p,r)==x(1-p^5h_p/5)                           (mod p^6). (4.2)
```

On the reflected side,

```text
b_(2p-1-s)
 =b_(2p-1)mathcal U_s(-2p)
  +8p^3b_(2p)mathcal J_s(-2p).
```

Use

```text
b_(2p-1)==5+8Delta_p+(336/5)p^5h_p                (mod p^6),
b_(2p)  ==73-824Delta_p                           (mod p^5),
```

and define

```text
Z6_(p,s)
 =[(5+8Delta_p)T_s^[6](-2p)
   +8p^3(73-824Delta_p)J_s^[3](-2p)]/p            (mod p^6).
                                                                    (4.3)
```

Then

```text
Z6_(p,s)==x(1-336p^5h_p/25)                       (mod p^6). (4.4)
```

The combination `336*(4.2)-5*(4.4)` proves `(0.4)`.

For a set `S` of top target primes, put

```text
R=product_(p in S)p,       A=b_n/R,
R_*=product_(p in S, p!=331)p.
```

For each `p|R_*`, `R/p` and `331` are units modulo `p^6`, and

```text
A
 ==[331(R/p)]^(-1)
   [336D6_(p,n-p)-5Z6_(p,2p-1-n)]                (mod p^6).
```

CRT determines `A modulo R_*^6`.

## 5. Independent audit

The reproducer

```text
../scripts/q32_weight_five_endpoint_sextic_target_audit.py
```

checks independently:

1. the endpoint laws `(0.1)` from their finite sums;
2. `h_p==-24B_(p-5)` and the `336/5` reflected ratio;
3. the exact degree-six direct and reflected recurrence jets;
4. `(4.2)`, `(4.4)`, and the `331` elimination;
5. recovery of the actual common cofactor modulo `p^6`.

The initial run over all primes `7<=p<=1000` gives 165 endpoint-prime
checks and 163 target-digit checks on each side, with zero failures.

## 6. Global status

This is the sixth consecutive common quotient digit.  It is a
genuine continuation of the local rank-one phenomenon, but order six
still sits just below the first information-theoretic threshold:

```text
7*(n/2)=3.5n < log b_n=3.525494...n+o(n).
```

Moreover, a residue class alone reconstructs the existing integer
and is not a contradiction.  The next useful targets are therefore:

```text
1. extend the endpoint filtration one more grade and test whether
   the weight-seven defect is again rank one;
2. produce a non-stabilization or small-lift theorem for the common
   quotient tower.
```
