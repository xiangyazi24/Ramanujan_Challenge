ANSWER Q1066 4dfbede2

# Audit of the exact Padé/Newton-content ratio

## Executive verdict

The proposed formula

```
D_n = 2^{e2(n)} 3^{e3(n)} Gamma_n
```

is not proved by the supplied finite computation. The gcd identity (T) is correct. The remaining statement (C) is a much stronger local statement about the Padé minors. It requires exact control of the 2-adic and 3-adic parts and a proof that every prime p>=5 contributes identically to D_n and Gamma_n. The finite scan does not establish this.

The narrow missing theorem is a primewise content comparison lemma:

```
v_p(gcd_j L_j(n)(b_n-b_j)) = v_p(gcd_j L_j(n)b_j)
```

for every p>=5.

No proof of this lemma follows from (T) alone.

---

# 1. Verification of (T)

Let

```
D = gcd_j L_j(b_n-b_j)
G = gcd_j L_j b_j.
```

Since

```
sum_j L_j = 1,
```
we have the integer identity

```
sum_j L_j(b_n-b_j) = b_n - sum_j L_j b_j.
```

Therefore every common divisor of D and G divides b_n. Hence

```
gcd(D,G) | gcd(D,b_n).
```

Conversely, if a number q divides D and b_n, then

```
L_j b_j = L_j b_n - L_j(b_n-b_j),
```

and both terms on the right are divisible by q. Hence q divides every generator of G, so

```
gcd(D,b_n) | G.
```

The same argument gives

```
gcd(G,b_n)=gcd(D,b_n).
```

This is an equality of ideals in Z, hence includes all prime powers, not only prime support.

Thus (T) is proved.

---

# 2. Primewise analysis of (C)

## Primes p>=5

For p>=5 the desired statement is

```
v_p(D_n)=v_p(Gamma_n).
```

The only information supplied by (T) is

```
min(v_p(D_n),v_p(Gamma_n))
 = min(v_p(D_n),v_p(b_n))
 = min(v_p(Gamma_n),v_p(b_n)).
```

This does not imply equality of v_p(D_n) and v_p(Gamma_n).

A missing step is required:

```
if p^a | Gamma_n then p^a | D_n,
```

or equivalently

```
L_j b_j == 0 (mod p^a) for all j
    ==> L_j(b_n-b_j) == 0 (mod p^a) for all j.
```

For p>=5 this would follow if one could prove a uniform congruence

```
b_n L_j == 0 (mod p^a)
```

whenever the Gamma generators have p^a-content. The current Apéry-Lucas/Dwork input does not provide this lifting statement.

Conclusion: the p>=5 part of (C) remains open.

---

## The 2-adic excess

The claimed formula

```
e2(n)=2 if n=2^a-1,
otherwise min(5,2+r(n))
```

is an empirical pattern, not a consequence of the displayed identities.

A proof would require computing

```
min_j v_2(L_j(n)(b_n-b_j))
-
min_j v_2(L_j(n)b_j)
```

and showing the difference is exactly the claimed value.

The Apéry recurrence and parity congruences are relevant, but no supplied lemma controls the minimising index j or excludes cancellation among the generators.

Conclusion: unproved.

---

## The 3-adic excess

Similarly, the statement

```
e3(n)=1 iff n=2*3^a or 2*3^a+1
```

requires a complete valuation analysis of the same gcd generators.

The congruence pattern may suggest where the extra factor comes from, but a recurrence congruence alone does not establish the exact gcd valuation.

Conclusion: unproved.

---

# 3. Relation with the old Dwork/content quantity

If the old Q888 quantity was defined by

```
c_n = gcd_j L_j(n)b_j,
```

then it is literally Gamma_n.

Under (C),

```
D_n = 2^e2 3^e3 c_n.
```

Using (T),

```
gcd(D_n,c_n)=gcd(c_n,b_n).
```

But the left side becomes

```
gcd(2^e2 3^e3 c_n,c_n)=c_n.
```

Therefore (C) would imply

```
c_n | b_n.
```

This is exactly the formerly open full divisibility statement.

Hence proving (C) would solve the old content divisibility problem.

---

# 4. Computation status

A claimed scan through n=1200 is evidence only. An extended executable verification requires the repository code or an independently executed exact script.

The necessary verifier is:

```python
for n in range(4,N):
    M=n//2
    L=[]
    for j in range(M):
        L.append((-1)**(M-1-j)
          * C(n,j)
          * C(n-j-1,M-1-j))
    Gamma=gcd(abs(L[j]*b[j]) for j in range(M))
    D=gcd(abs(L[j]*(b[n]-b[j])) for j in range(M))
    assert D == 2**e2(n)*3**e3(n)*Gamma
```

No new numerical range is certified here without actually running this code against the repository data.

---

# 5. Height consequence

Even if (C) is proved, the correct consequence is only

```
log D_n = log Gamma_n + O(1),
```

because the ratio is supported only at 2 and 3 with bounded exponent.

This does not imply

```
log D_n=o(n).
```

An independent estimate on Gamma_n is still required.

---

# Final classification

| Statement | Status |
|---|---|
| Pairwise gcd identity (T) | Proven |
| Identification Gamma_n=c_n (if Q888 used same definition) | Yes |
| p>=5 equality of valuations | Missing p-adic lifting lemma |
| Exact 2-adic excess | Open |
| Exact 3-adic excess | Open |
| Formula (C) | Not proved |
| Height conclusion | Only log D=log Gamma+O(1) |

The narrowest sufficient theorem is the p-adic comparison lemma for every prime, together with exact 2-adic and 3-adic generator valuation lemmas. The finite scan strongly motivates these statements but does not replace them.