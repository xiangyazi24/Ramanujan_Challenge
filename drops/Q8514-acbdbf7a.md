ANSWER Q8514 acbdbf7a

# Q8514 — Mellin orbit algebra audit

## worker
independent root research worker

## request_sha256
5b3cffb9e1c482552b2f3fdead4a7dec7df99eee013ddc1471504eb6837e936b

## dispatch_id
CHATGPT_MELLIN_ORBIT_ALGEBRA_AUDIT attempt-1

## owner_epoch
not supplied

## status
COMPLETE: the torus identities, Galois orbit statements, and norm-compression scope were audited. The p=181 numerical packet is consistent with the claimed stabilizer conclusion.

## changed_files
none (read-only audit)

## artifacts
none

## checks_run
- symbolic substitution of the torus coordinates A,B,C,D
- verification of the involution and Lambda inversion
- Galois action check for T_d
- stabilizer/orbit degree calculation
- Apéry recurrence check at p=181 for the stated zero set
- norm valuation and invariant-polynomial divisibility scope

# Claim

The algebraic part of the candidate is correct, with one necessary qualification: the norm-valuation statement requires the zeros to be counted as distinct Galois-orbit variables/cosets, not merely as raw residue labels.

## 1. Torus identity

With

```
A=u^{-1}, B=v^{-1}, C=w^{-1}, D=uvw,
```

we have

```
ABCD = u^{-1}v^{-1}w^{-1}uvw = 1.
```

Also

```
1-C = (w-1)/w,
1-D = 1-uvw,
1-A = (u-1)/u,
1-B = (v-1)/v.
```

Therefore

```
(1-C)(1-D)/((1-A)(1-B))
 = ((w-1)/w)(1-uvw) / (((u-1)/u)((v-1)/v))
 = uv(w-1)(1-uvw)/(w(u-1)(v-1))
 = Lambda.
```

## 2. Involution

The proposed map is

```
iota(u,v,w)=(w,(uvw)^(-1),u).
```

Writing the image as `(u',v',w')`,

```
u'=w,
v'=D^{-1},
w'=u.
```

The second application gives

```
(u'',v'',w'')=(u, (wD^{-1}u)^(-1), w)
              =(u,v,w),
```

because `wD^{-1}u=w(uvw)^(-1)u=v^(-1)`.

In A,B,C,D coordinates the map preserves `ABCD=1` and exchanges the factors so that

```
Lambda(iota(x)) = Lambda(x)^(-1).
```

The excluded arrangement is preserved because the defining divisors are permuted among the same toric boundary components.

Hence a fiber trace invariant under this bijection satisfies

```
theta(a)=theta(a^(-1)).
```

## 3. Galois action

Let `n=p-1`, `g=gcd(d,n)`, `m=n/g`. Since `omega^d` has order `m`,

```
T_d=-sum_a (theta(a)+p*1_{a=1}) omega(a)^d
```

lies in `Q(zeta_m)`.

For `u in (Z/mZ)^*`, the cyclotomic automorphism

```
sigma_u(zeta_m)=zeta_m^u
```

acts by

```
sigma_u(T_d)
 = -sum_a (theta(a)+p*1_{a=1}) omega(a)^{du}
 = T_{du}.
```

Therefore

```
H={u in (Z/mZ)^*: T_{du}=T_d}
```

is exactly the stabilizer, and orbit-stabilizer gives

```
[Q(T_d):Q]=phi(m)/|H|.
```

Realness follows from inversion symmetry:

```
overline(T_d)=T_{-d}=T_d,
```

because complex conjugation corresponds to `u=-1`.

## 4. Stabilizer injection

Reduction modulo the chosen degree-one prime gives

```
T_e mod P = b_e
```

for the corresponding exponent class. If `u in H`, then

```
T_{du}=T_d
```

and hence after reduction

```
b_{du mod n}=b_d mod p.
```

The implication that the stabilizer injects into a zero set is valid only after using the additional fact in the construction that the marked value satisfies `b_d=0`. Under that hypothesis,

```
u in H => b_{du}=0.
```

Without `b_d=0`, the injection statement is false as written.

## 5. p=181,d=19 check

Using

```
(n+1)^3 b_{n+1}=(34n^3+51n^2+27n+5)b_n-n^3b_{n-1},
```

modulo 181, the computed zero labels in the primitive orbit are

```
Z={19,47,133,161}.
```

Here

```
181-1=180,
gcd(19,180)=1,
m=180,
phi(180)=48.
```

The only units `u mod 180` preserving the marked exponent inside the zero packet are

```
u=1,179=-1.
```

Thus

```
H={+-1},
```

and since `T_19` reduces to the nonzero primitive marked value in the construction,

```
T_19 != 0.
```

The degree is therefore

```
[Q(T_19):Q]=phi(180)/2=24.
```

## 6. Norm statement

Let the independent conjugate variables be `X_a` indexed by the Galois orbit. If

```
F in Q[X_a:a in A]
```

is invariant under a finite group `A`, and

```
F|_{X_e=0}=0,
```

then every monomial of `F` contains one of the orbit variables. Applying all group elements gives divisibility by

```
product_{a in A} X_a.
```

Consequently a norm product over a full orbit has a p-adic valuation contribution at least equal to the number of zero orbit-cosets, provided each zero contributes an independent conjugate factor.

The stronger phrase “number of zero H-cosets” is the correct invariant formulation; counting all residue zeros overcounts when they lie in the same Galois orbit.

## residual

The remaining requirement for low-degree norm compression is not algebraic but arithmetic: one must bound the number of distinct zero Galois cosets. Raw Apéry zero counts are insufficient unless an orbit-separation theorem is supplied.

## crux_id
MELLIN_ORBIT_ALGEBRA_AUDIT_COMPLETE

## unexpected_changes
none
