ANSWER Q791 bb0f9fec

# Q791: additive-Sidon zero-set work versus LGTQ-2

## 1. Audit of the available zero-set statements

The additive-Sidon/difference-sparse computations separate into two classes.

### Proven uniformly for a fixed prime characteristic

The valid statements are of the form:

```
Z_p = {u mod p : p | b_u}
```

has controlled additive structure, for example bounds on

```
E_+(Z_p)=#{(a,b,c,d) in Z_p^4:a+b=c+d}
```

or reflection-normalized difference multiplicities:

```
#{u in Z_p: u and u+h both belong to Z_p}
```

after the reflection identification. These statements constrain pairs living inside one field F_p.

### Only computationally scanned

Statements such as:

```
max_h #{fixed gap zero pairs}
```

for many numerical primes are evidence only. They do not imply an all-prime theorem unless a proof was supplied.

In particular, a finite scan of Apéry zero sets cannot be promoted to an asymptotic LGTQ-2 estimate.

## 2. Why one-characteristic difference sparsity does not control LGTQ-2

The LGTQ-2 pair is

```
p, ell=p-d,
 u=n-q p,
 v=u+q d.
```

The two divisibility conditions are

```
p | b_u,
ell | b_v.
```

The first condition lives in F_p:

```
u in Z_p.
```

The second lives in F_ell:

```
v in Z_ell.
```

Although the integer representatives differ by the deterministic shift qd, the characteristics are different. A theorem controlling

```
Z_p intersect (Z_p+h)
```

says nothing directly about

```
Z_p and Z_{p-d}.
```

There is no valid projection from a cross-characteristic pair to a one-characteristic additive configuration.

## 3. Reflection-preserving countermodel

Take abstract zero sets satisfying the strongest one-prime difference condition:

```
|Z_p intersect (Z_p+h)| <= M(h)
```

with the required reflection symmetry

```
u in Z_p iff -1-u in Z_p.
```

Choose independently for each neighboring prime pair (p,p-d):

```
u_p in Z_p,

u_p+qd in Z_{p-d}.
```

Because the second set is a different characteristic, this does not increase any additive energy of either individual Z_p.

For a fixed row q,h one can arrange

```
#{p: p and p-d prime, u_p in Z_p, u_p+qd in Z_{p-d}}
~ n/log^2(n)
```

while keeping

```
E_+(Z_p)
```

and every reflection-minimal difference multiplicity bounded as required.

Therefore one-prime Sidon information alone cannot imply LGTQ-2.

## 4. Consequence for LGTQ-2

No unconditional saving follows from the existing additive-energy theorem.

The desired estimate is a mixed-characteristic correlation:

```
#{p: p|b_{n-qp}, p-d|b_{n-qp+qd}} = o(n/log^2 n).
```

This requires new information coupling neighboring characteristics.

Possible sufficient statements would be:

1. a two-modulus correlation estimate;
2. a resultant estimate involving both

```
b_u mod p
```

and

```
b_{u+qd} mod (p-d);
```

3. a character sum with both moduli present.

The existing additive-Sidon theorem gives none of these.

## 5. Exact surviving implication

The only valid deduction is local:

For each fixed p,

```
#{h: u,u+h in Z_p}
```

is bounded by the one-prime difference theorem.

It may reduce repeated gaps inside a single characteristic, but it does not bound the number of affine neighboring-prime matches.

## 6. Conclusion

The additive-Sidon/difference-sparse work is useful as a diagnostic and may control aligned structures after both endpoints are forced into one characteristic. It does not close LGTQ-2, because LGTQ-2 is intrinsically a cross-characteristic affine correlation problem.

Any proof of LGTQ-2 must introduce a genuinely two-prime invariant; promoting scanned zero-set sparsity to a theorem would be an invalid step.
