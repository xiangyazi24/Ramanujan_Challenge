ANSWER Q2699 daacf09f

# P3.2 — selected-state repair of the AH-M1 carrier

## Verdict

There **is** an exact characteristic-zero repair of the singular `OR` defect that is smaller and cleaner than multiplying back `b_r` or `b_s`.

The right object is not the raw determinant and not its first jet.  It is the determinant **together with the two endpoint principal cofactors** (equivalently, the rank-2 / second-determinantal-divisor mark) at the same half-integer center.

Write, as in `problems/3.2/meso_result.tex`,

```text
T_H^(d) = 2^(3(H-1)) N_H((d-1)/2)
```

for even `d`, with `N_0=0`, `N_1=1` and

```text
N_(m+1)(x) = P(x+m) N_m(x) - (x+m)^6 N_(m-1)(x).
```

Then the two endpoint principal cofactors of the dyadically cleared tridiagonal matrix are exactly

```text
T_(H-1)^(d),      T_(H-1)^(d+2).
```

For an actual pair of distinct lower reflection orbits

```text
0 < r < s < (p-1)/2,
h = s-r,
d = p-1-2r,
H = r+s+1,
```

one has

```text
p = d+H-h.
```

Define the atomic marked carrier

```text
p0(H,h,d) = d+H-h,

C_atom(H;h,d)
 = gcd(
     p0(H,h,d),
     T_H^(d),
     T_(H-1)^(d),
     T_(H-1)^(d+2)
   ).                                                   (A)
```

Assume

```text
d even,
h>=1,
d>2h,
H>=h+3,
H == h+1 (mod 2),
p0=d+H-h is an odd prime.
```

Put

```text
r=(H-h-1)/2,
s=(H+h-1)/2.
```

Then `0<r<s<(p0-1)/2` and the following is exact:

```text
boxed:
C_atom(H;h,d) = p0
    <=>
p0 | b_r  and  p0 | b_s.                            (B)
```

Otherwise `C_atom(H;h,d)=1`.

Thus a single lower zero cannot create an `H`-tail: at the crossing it gives rank defect one, whereas (A) requires rank defect at least two.  No Apéry value occurs in the definition.

If the requested indexing must literally be only `(h,d)`, use the dyadic projection

```text
C_X(h,d)
 = product over H satisfying
       H>=h+3,
       H == h+1 (mod 2),
       X < d+H-h <= 2X
   of C_atom(H;h,d).                                  (C)
```

No primality test is needed in the definition.  For a prime `X<p<=2X`, a factor `p` of `C_X(h,d)` can only come from the unique `H=p-d+h`, because the linear factor `d+H-h` itself lies in `(X,2X]`.

Consequently

```text
boxed:
p | C_X(h,d)
<=>
there are distinct lower zero representatives r<s with
s-r=h,
p-1-2r=d.                                           (D)
```

This is the clean dyadic two-parameter carrier I recommend.  It satisfies the requested logical conditions (i), (ii), and (iv) exactly.  Its dyadic radical aggregate has **no spectral plateau at all**; it is exactly the genuine lower-orbit-pair mass.  What is *not* presently proved is the final `O(X polylog X)` estimate for that genuine pair mass.  The remaining input is the vertical second moment of `Z(p)`, not another determinant identity.

More precisely, if `z_p` is the number of lower representatives of noncentral reflection orbits, then summing over all admissible `(h,d)` gives the exact identity

```text
sum_(h,d) log rad_(X,2X]( C_X(h,d) )
  = sum_(X<p<=2X) binom(z_p,2) log p.                  (E)
```

Hence the desired height bound would follow from the genuinely arithmetic estimate

```text
sum_(X<p<=2X) z_p^2 log p << X (log X)^A             (Z2)
```

for some fixed `A`.  The present repository does not prove `(Z2)`.  The known pointwise `Z(p) << p^(2/3)` only gives a superlinear bound of order at best `X^(7/3+o(1))` here.  So the repair removes the false spectral mass; it does not manufacture the still-missing average zero-count theorem.

That is the main classification.

---

# 1. Source audit

I audited the connector-visible current `main` at

```text
47fa0e653f52c4a71e9a8c26b31ca9f66f6bbe86
```

and the previous hostile audit `drops/Q2693-84db657a.md` on `chatgpt-drop`.

The relevant established pieces are:

1. `problems/3.2/meso_result.tex`:
   - the continuant recurrence for `N_h`;
   - the symmetric tridiagonal realization `det A_h(x)=N_h(x)`;
   - reflection
     ```text
     N_h(-x-h-1)=(-1)^(h-1) N_h(x);
     ```
   - the integral half-center values `T_h^(d)`;
   - the root strip `-h < Re(alpha) < -1`.

2. `problems/3.2/fiber_bound.tex`:
   - the factorial-gauged solution
     ```text
     Y_m = Pi_m(r) b_(r+m);
     ```
   - the exact fundamental-solution formula
     ```text
     Y_m = N_m(r) Y_1 + B_m(r) Y_0;
     ```
   - the return-row / bordered determinant machinery and the exact Casoratian factorization.

3. `problems/3.2/nv_theorem.tex`:
   - nondegeneration of the bordered equal-fiber certificate in the nonsingular range;
   - the fact that one must distinguish genuine rank conditions from structural factors.

4. `drops/Q2693-84db657a.md`:
   - on a singular crossing `p=d+2r+1`, with `t=H-1-r`, the raw AH-M1 carrier satisfies
     ```text
     T_H^(d)
       == unit * b_r * b_t                 (mod p),
     ```
     hence its zero condition is `b_r=0 OR b_t=0`;
   - at `X=240,p=271`, 470 of the 474 singular-fan hits are XOR false positives.

I found no definition or theorem named `Racah K` in the current `Ramanujan_Challenge` repository, and no established theorem identifying a first-jet sequence with a Racah polynomial and proving the required coprimality with the endpoint Apéry value.  I therefore do not import such a theorem from another project or silently assume it.

---

# 2. The exact two-border rank theorem

Let

```text
x=(d-1)/2,
n=H-1.
```

Let `A_H(x)` be the `(H-1)x(H-1)` symmetric tridiagonal matrix from `meso_result.tex`:

```text
diagonal i:       P(x+i),
off-diagonal i:   (x+i+1)^3.
```

For even `d`, clear the half denominators by

```text
M_(H,d) = 8 A_H((d-1)/2)  in Mat_(H-1)(Z).             (2.1)
```

Then

```text
det M_(H,d) = T_H^(d).                                 (2.2)
```

The last and first principal cofactors are

```text
cof_nn M_(H,d) = T_(H-1)^(d),
cof_11 M_(H,d) = T_(H-1)^(d+2).                        (2.3)
```

The reason is simply that deleting the last row/column leaves the leading tridiagonal block, while deleting the first shifts `x` to `x+1`.

Now take an odd prime `p` and suppose

```text
p=d+2r+1,
1<=r<=H-2,
t=H-1-r,
r,t<p.                                                  (2.4)
```

The scaled off-diagonal in position `i` is

```text
8(x+i+1)^3 = (d+2i+1)^3.                               (2.5)
```

Hence modulo `p` exactly the `i=r` off-diagonal vanishes.  The matrix splits as

```text
M_(H,d) mod p = L_r direct_sum R_t.                    (2.6)
```

All internal off-diagonals in both blocks are nonzero.  Therefore each block has nullity at most one.

The two block determinants are the factorial-gauged Apéry values:

```text
det L_r = unit * (r!)^3 b_r,
det R_t = unit * (t!)^3 b_t.                           (2.7)
```

For the right block this is literally the recurrence

```text
X_0=1,
X_1=5,
X_(j+1)=P(j)X_j-j^6 X_(j-1),
X_j=(j!)^3 b_j.
```

For the left block it is the reflected version, exactly as in Q2693.

Thus

```text
nullity M_(H,d) mod p
 = 1_(p|b_r) + 1_(p|b_t).                              (2.8)
```

This immediately gives the theorem.

## Theorem 2.1 — two-border AND detector

Under (2.4),

```text
p divides all three of
  T_H^(d), T_(H-1)^(d), T_(H-1)^(d+2)

<=>

p|b_r and p|b_t.                                       (2.9)
```

### Proof

If both blocks are singular, the full matrix has nullity two, so its determinant and every `(n-1)x(n-1)` cofactor vanish.

If exactly one block is singular, the full matrix has nullity one.  Its kernel is supported on that singular block.  Because an irreducible tridiagonal block has nonzero endpoint coordinates in every nonzero kernel vector, the principal cofactor at the corresponding outer endpoint is nonzero.  Therefore the two endpoint cofactors cannot both vanish.

If neither block is singular, the determinant is nonzero.

That proves (2.9).

This is the precise rank-theoretic replacement for the raw determinant's `OR`.

---

# 3. The canonical Smith version and the smallest practical version

There are two equivalent ways to package the same rank condition.

## 3.1 Canonical Smith carrier

Let

```text
Delta_(n-1)(M_(H,d))
```

be the gcd of all `(n-1)x(n-1)` minors of `M_(H,d)`.  This is the next-to-top determinantal divisor (the product of the first `n-1` Smith invariant factors).

A prime divides this determinantal divisor iff

```text
rank(M mod p) <= n-2.
```

Therefore on the crossing (2.4),

```text
p | Delta_(n-1)(M_(H,d))
<=> p|b_r and p|b_t.                                   (3.1)
```

This is the most canonical selected-state repair of AH-M1.

It also has an immediate structural height divisor.  Delete the last row and the first column.  The resulting minor is triangular, with determinant

```text
product_(i=1..H-2) (d+2i+1)^3.                        (3.2)
```

Hence

```text
Delta_(n-1)(M_(H,d))
  | product_(i=1..H-2) (d+2i+1)^3.                    (3.3)
```

So its prime support is forced onto the explicit crossing progression.  This is exactly what the raw determinant lacked.

## 3.2 Smaller practical carrier

Computing every cofactor is unnecessary for the large-prime AH range.  Define

```text
Qcirc_(H,d) = product_(i=1..H-2) (d+2i+1),

G_(H,d) = gcd(
   T_H^(d),
   T_(H-1)^(d),
   T_(H-1)^(d+2),
   Qcirc_(H,d)
).                                                     (3.4)
```

Then `G_(H,d) | Qcirc_(H,d)`, and for every prime in the AH dyadic range its support is exactly the same AND support as the full Smith carrier.

If one wants to delete the diagonal self-pair `r=t`, note that on a crossing

```text
r=t <=> p=d+H.                                         (3.5)
```

Thus the distinct-orbit version is

```text
Gneq_(H,d) = sat_(d+H)( G_(H,d) ),                    (3.6)
```

where `sat_m` removes all prime-primary factors supported on the primes dividing `m`.

For a genuine distinct pair `r!=t`, the crossing prime does not divide `d+H`:

```text
d+H = p+(t-r),
```

and `0<|t-r|<p`, so the saturation does not remove it.

This is an exact characteristic-zero replacement for the explicit marked product in Q2693, but it never mentions `b_r` or `b_t` in its definition.

---

# 4. Gap coordinates and the requested `C_(h,d)`

The only subtlety is indexing.

For the actual lower pair

```text
r<s,
h=s-r,
d=p-1-2r,
```

the old singular AH height is **not** the gap `h`.  It is

```text
H=r+s+1.                                              (4.1)
```

Equivalently,

```text
p=d+H-h.                                              (4.2)
```

This is why a rank-two repair of the old AH matrix naturally has the extra sum-coordinate `H`.

The atomic carrier (A) removes this nuisance by inserting the linear candidate prime `p0=d+H-h`.  For the dyadic scale, define (C):

```text
C_X(h,d)
 = product_H gcd(
       d+H-h,
       T_H^(d),
       T_(H-1)^(d),
       T_(H-1)^(d+2)
   ),                                                  (4.3)
```

where

```text
H>=h+3,
H == h+1 (mod 2),
X < d+H-h <= 2X,
d even,
d>2h.                                                 (4.4)
```

The parity and inequality give

```text
r=(H-h-1)/2,
s=(H+h-1)/2,
0<r<s<(p0-1)/2.                                       (4.5)
```

## Theorem 4.1 — exact dyadic gap support

For every prime `X<p<=2X`,

```text
p | C_X(h,d)
```

if and only if

```text
r=(p-d-1)/2,
s=r+h
```

are distinct lower representatives and

```text
p|b_r,
p|b_s.                                               (4.6)
```

### Proof

If the pair exists, take `H=p-d+h`.  Then (4.5) gives the original `r,s`, and Theorem 2.1 makes all three continuant factors divisible by `p`; the linear factor is literally `p`.

Conversely, suppose `p` divides one factor in (4.3).  It divides the linear integer `p0=d+H-h`, which lies in `(X,2X]`.  Since `p>X`, one has `p0<2p`; hence `p|p0` forces `p0=p`.  Theorem 2.1 then forces both zero conditions in (4.6).

No raw determinant root away from the crossing survives, and a single zero block cannot survive the two endpoint cofactors.

This is the requested dyadic selected-state carrier in literal `(h,d)` coordinates.

---

# 5. Exact radical aggregate: the plateau is gone

Let `z_p` be the number of lower representatives of **noncentral** reflection orbits:

```text
z_p = #{ 0<r<(p-1)/2 : p|b_r }.
```

Each unordered pair `r<s` determines exactly one

```text
h=s-r,
d=p-1-2r.
```

Conversely Theorem 4.1 recovers that pair from `(p,h,d)`.

Therefore

```text
boxed:
Sum_(admissible h,d)
  log rad_(X,2X]( C_X(h,d) )
=
Sum_(X<p<=2X) binom(z_p,2) log p.                      (5.1)
```

This is the strongest reason to prefer the two-border/Smith mark over the edge determinant alone: the dyadic radical is now **exactly genuine pair mass**.  There is no spectral plateau and no regular spectral false-positive term hidden in (5.1).

Equation (5.1) also states the remaining analytic problem without renaming SDC:

```text
need:  Sum z_p^2 log p << X polylog X.                 (5.2)
```

This is a vertical second moment of local zero counts.  It is not a cross-prime dispersion statement and not a restatement of the earlier TCUS/SDC hypotheses.

Current status:

- finite data (`zp_million_output.txt` and the HM3 computations) make an `X polylog X` scale plausible;
- the current theorem `Z(p)<<p^(2/3)` is far too weak by itself;
- no proof of (5.2) is present in the repository.

Thus I would call (A)/(C) a **successful carrier repair but not yet a completed global height theorem**.

---

# 6. Concrete fixture: `p=271`

Q2693 records, exactly,

```text
Z_271 intersect [0,134] = {17,41}.
```

Take the two distinct lower orbits

```text
r=17,
s=41,
h=24,
d=271-1-34=236,
H=r+s+1=59.
```

Then

```text
p0=d+H-h=236+59-24=271.
```

The scaled `58x58` tridiagonal matrix `M_(59,236)` has its unique zero off-diagonal modulo `271` at `i=17`; it splits into blocks of sizes `17` and `41`.

Both block determinants vanish because

```text
271|b_17,
271|b_41.
```

Hence

```text
271 divides
T_59^(236),
T_58^(236),
T_58^(238),
```

and therefore

```text
C_atom(59;24,236)=271.                                (6.1)
```

By contrast, in the 474-hit singular fan from Q2693, the raw determinant retained every OR hit.  The rank-two mark retains only the AND intersections.  At this prime the two low zeros give four ordered `(r,t)` intersections if the diagonal is allowed; after restricting to distinct lower orbits there are the two orientations, or one unordered pair.  The 470 XOR hits disappear identically, without using an empirical threshold.

This is the requested concrete regression fixture.

---

# 7. Audit of the direct gap-edge carrier

There is a very small two-parameter necessary carrier worth recording.  Define

```text
E_(h,d)
 = 2^(3(h-1)) N_h(-(d+1)/2).                           (7.1)
```

Using the centered integer polynomial from `meso_result.tex`,

```text
M_h(Y)=2^(3(h-1)) N_h((Y-h-1)/2),
```

this is simply

```text
E_(h,d)=M_h(h-d).                                      (7.2)
```

Reflection also gives

```text
E_(h,d)=(-1)^(h-1) T_h^(d-2h).                        (7.3)
```

If `p=d+2r+1`, `p|b_r`, and `r+h<(p-1)/2`, then `-(d+1)/2 == r (mod p)`.  The exact fundamental-solution identity from `fiber_bound.tex` says

```text
Pi_h(r)b_(r+h)
 = N_h(r)(r+1)^3 b_(r+1) + B_h(r)b_r.                 (7.4)
```

Since `b_r=0`, no consecutive zeros gives `b_(r+1)!=0`, and every factorial factor is a `p`-unit.  Therefore

```text
boxed:
p | E_(h,d)  <=>  p | b_(r+h),
provided p|b_r and r+h<(p-1)/2.                       (7.5)
```

So `E_(h,d)` by itself already kills the **specific** single-zero `h`-tail for the target orbit.

Also `d>2h` for two distinct lower orbits, so

```text
-(d+1)/2 < -h.
```

The root-strip theorem in `meso_result.tex` shows `E_(h,d)` is a nonzero characteristic-zero integer.

For the `p=271` fixture,

```text
E_(24,236)=(-1)^23 T_24^(188),
```

and with the mark `271|b_17`, equation (7.5) says

```text
271|E_(h,236)
<=>
271|b_(17+h)
```

through the physical lower range.  Since the only lower zeros are `17,41`, the target prime occurs only at `h=24`, not on a tail.

### Why I do not stop at `E_(h,d)`

`E` is a transition-entry condition.  Without the initial-state hypothesis `p|b_r`, an ordinary root `N_h(r)=0` can still occur.  In other words, `E` selects the zero-to-zero **direction** of the local recurrence, not the fact that the Apéry solution is the solution occupying that direction.

The two-border construction in Sections 2--4 supplies exactly that missing state mark without inserting `b_r`.

---

# 8. Gcd with the reflected one-sided border: exact no-gain theorem

One might try to intersect (7.1) with the reflected long edge.  Let

```text
x=-(d+1)/2,
1<=h<d/2,

E1 = 2^(3(h-1))       N_h(x),
E2 = 2^(3(d-h-1))     N_(d-h)(x).                     (8.1)
```

This gives no new information.

## Proposition 8.1

One has the exact absolute-value identity

```text
|E2|
 = |E1| * ((d-2h-1)!!)^6.                             (8.2)
```

Consequently

```text
gcd(|E1|,|E2|)=|E1|.                                  (8.3)
```

### Proof

At the central point `x=-(d+1)/2`, the `(d-1)x(d-1)` tridiagonal matrix satisfies

```text
A_d(x) = -J A_d(x) J,
```

so, because `d-1` is odd, it has a zero mode.  Every off-diagonal is nonzero in characteristic zero, hence the kernel is one-dimensional.

For a kernel vector normalized at the left endpoint, the standard continuant formula is

```text
v_i = (-1)^(i-1) N_i(x)
      / product_(j=1..i-1) (x+j+1)^3.                 (8.4)
```

The reflected vector `Jv` spans the same one-dimensional kernel, so the components at `i=h` and `i=d-h` agree up to sign.  Clearing the products gives

```text
N_(d-h)(x)
 = +/- N_h(x)
     product_(j=h..d-h-1) (x+j+1)^3.                  (8.5)
```

After multiplying by the powers of two in (8.1), the linear factors are

```text
2j+1-d,
```

which run symmetrically through

```text
-(d-2h-1), ...,-1,1,...,(d-2h-1).
```

Their product, cubed, has absolute value `((d-2h-1)!!)^6`.  This proves (8.2).

So the reflected-border gcd is literally the original edge carrier again.  It cannot improve the radical.

---

# 9. First-jet audit

There are two different first-jet ideas, and they should not be conflated.

## 9.1 First jet of the structurally central orbit block — fails

For even `d`, `N_d(-(d+1)/2)=0` identically by reflection.  Its first jet does **not** mark whether the Apéry orbit is a zero orbit.

The repository gives the concrete actual zero set

```text
Z_37={17,19}.
```

This is one reflection orbit, with `r=17,d=2`.  Since

```text
N_2(x)=P(x+1),
```

at the characteristic-zero center `x=-3/2`,

```text
N_2'(-3/2)=P'(-1/2)=3/2.
```

For the integral centered polynomial

```text
M_2(Y)=8 N_2((Y-3)/2),
```

one gets

```text
M_2'(0)=6,
```

so

```text
37 does not divide M_2'(0),
```

even though `37|b_17`.  Thus the central first jet is not a zero-orbit mark.

## 9.2 First jet of the **old singular AH determinant** — promising but conditional

This is a different object.  Put

```text
W_n(z)=N_(n+1)(-1+z),
X_n=W_n(0)=(n!)^3 b_n,
K_n=W_n'(0).                                          (9.1)
```

At a raw singular crossing, let

```text
H=r+t+1,
x0=-r-1.
```

Because the coupling coefficient across the crossing is `z^6`, the restart factorization actually holds through the first five jets:

```text
N_H(x0+z)
 == (-1)^r W_r(-z) W_t(z)        (mod z^6).            (9.2)
```

In particular

```text
N_H(x0)=(-1)^r X_r X_t,
N_H'(x0)=(-1)^r ( X_r K_t - K_r X_t ).                (9.3)
```

Therefore the value/first-jet gcd has exact singular support

```text
(pair: X_r=X_t=0)
OR
(left singleton: X_r=K_r=0)
OR
(right singleton: X_t=K_t=0).                         (9.4)
```

So a theorem of the form

```text
p|X_n and n<p/2  ==>  p does not divide K_n            (K-simple)
```

would indeed make the first-jet gcd convert `OR` to `AND`.

But `(K-simple)` is **not proved in the current repository**, and the desired actual-pair raw height is

```text
H=r+s+1=p-d+h,
```

which still depends on `p`, not just `(h,d)`.

Thus the first-jet route has two precise missing pieces:

1. endpoint simplicity / coprimality `gcd(X_n,K_n)` in the required characteristic range;
2. elimination of the `p`-dependent sum-coordinate `H` if one insists on a single non-dyadic `C_(h,d)`.

This is the exact point where a genuine Racah-`K` theorem could help.  Since no such theorem is defined or established in the present repo, I do not use it.

The two-border rank test avoids `(K-simple)` entirely.

---

# 10. Racah `K` audit

There is no connector-visible `Racah K` definition in the current project state, so only a scoped statement is possible.

There are two natural meanings a `K` could have here:

1. **Endpoint jet** `K_n=W_n'(0)` from (9.1).  Then the exact missing theorem is `(K-simple)` above.  A product/resultant formula proving that every prime divisor of `gcd(X_n,K_n)` lies below the physical range would make the first-jet carrier valid.

2. **Casoratian / Christoffel-Darboux kernel** of the two local fundamental solutions.  This does not mark the Apéry state.  In the nonsingular physical interval the fundamental transfer matrix has determinant
   ```text
   product_(j=1..h) (r+j)^6,
   ```
   which is a `p`-unit.  Such a kernel is therefore a nonvanishing normalization, not a second zero condition.

So the first missing hypothesis is not “some Kloosterman cancellation” or another name for SDC.  It is a concrete local algebraic statement: **which primes can divide the endpoint value and its selected jet simultaneously?**

Until a `Racah K` is explicitly identified with one of these objects and its coprimality theorem is proved in this repository, it cannot be used as the carrier mark.

---

# 11. `2x2` marked determinant audit

For the direct gap transfer, the two fundamental solutions form

```text
F_h(x) = [[N_(h+1)(x), B_(h+1)(x)],
          [N_h(x),     B_h(x)]].                        (11.1)
```

The recurrence gives the exact Casoratian

```text
det F_h(x)
 = product_(j=1..h) (x+j)^6.                           (11.2)
```

At an actual lower pair `r,r+h<p`, every factor in (11.2) is a `p`-unit.  Hence the local `2x2` determinant is **nonzero**, even when the selected entry `N_h(r)` vanishes.

This is the linear-algebra reason a bare `2x2` transfer determinant cannot be the desired pair carrier: zero-to-zero for a specified second-order solution is a vanishing **matrix entry / boundary functional**, not a loss of invertibility of the transfer.

If one instead forms two symmetric boundary rows, reflection makes them dependent in exactly the way quantified by Proposition 8.1.  One gets the same edge carrier times an explicit sixth power, not an independent gcd.

The old singular AH matrix is different: the crossing already makes one off-diagonal zero and splits the matrix.  There, two selected zero blocks really do mean rank defect two, and the endpoint-cofactor/Smith construction is exactly the right rank test.  That is why Sections 2--4 work whereas a direct-gap `2x2` determinant does not.

---

# 12. Height audit

There are three different height statements, and only one is currently justified at the desired structural level.

## 12.1 Raw characteristic height of the edge integer

For `E_(h,d)` in (7.1), Hadamard gives

```text
log |E_(h,d)| << h log(d+h).                            (12.1)
```

Summed naively over an `O(X^2)` rectangle this is cubic-scale.  The root-strip theorem even shows that on subregions with the evaluation point a fixed positive distance to the left of the strip, these values are exponentially large in `h`; so ordinary absolute size is not the mechanism that will give `O(X polylog X)`.

## 12.2 Fixed-gap dyadic large-prime radical of the edge integer

For a fixed `h` and prime `p`, the condition

```text
p|E_(h,d), p>d
```

corresponds to a root `r=(p-d-1)/2` of `N_h(r)` modulo `p`.  Since

```text
deg N_h=3(h-1),
```

there are at most `3(h-1)` such `d` for each `(p,h)`.  Hence

```text
sum_d log rad_(X,2X](E_(h,d)) << h X.                  (12.2)
```

A dyadic block `h~H` gives only

```text
<< X H^2,                                               (12.3)
```

which is `X polylog X` only when `H` itself is polylogarithmic.  Thus the edge carrier alone does not close the all-height problem.

## 12.3 Marked dyadic carrier

For `C_X(h,d)` the exact identity (5.1) replaces the `H^2` spectral root count by the actual pair count.  This is the correct structural endpoint:

```text
spectral plateau cost = 0,
remaining cost = genuine lower-orbit pairs.             (12.4)
```

An `O(X polylog X)` theorem is now equivalent to a vertical second-moment estimate like `(Z2)`.  That is plausible and sharply scoped, but not proved by the present local algebra.

---

# 13. Final classification

### What works

The smallest robust **atomic** repair is

```text
C_atom(H;h,d)
 = gcd(
     d+H-h,
     T_H^(d),
     T_(H-1)^(d),
     T_(H-1)^(d+2)
   ).
```

For a prime candidate `p0=d+H-h` it is exactly `p0` iff the two distinct lower representatives

```text
r=(H-h-1)/2,
s=(H+h-1)/2
```

are both Apéry zeros modulo `p0`.

The literal two-parameter dyadic carrier is the product `C_X(h,d)` over admissible `H`.  Its target-window radical support is exactly the actual distinct lower-orbit-pair support.

### Why it fixes Q2693

The raw determinant detects nullity at least one, hence `OR`.  The two endpoint cofactors detect nullity at least two, hence `AND`.  A single zero block cannot generate a tail.

### Why it is preferable to the first jet

The first jet detects multiplicity only after an extra endpoint-simplicity theorem.  The cofactor rank test has no such exception.

### Why the direct reflected-border gcd fails

The two edges are exactly proportional by `((d-2h-1)!!)^6`, so the gcd is unchanged.

### Why a bare local `2x2` determinant fails

The direct-gap transfer is invertible; its determinant is a factorial unit.  Pair information is a selected boundary state, not transfer singularity.

### Remaining theorem

The carrier problem is repaired.  The remaining quantitative gate is

```text
sum_(X<p<=2X) z_p^2 log p << X polylog X,
```

or an equivalent average bound for genuine lower reflection-orbit pairs.  Current `Z(p)<<p^(2/3)` does not prove it.

So the answer is **not** to reuse the failed unmarked `A_(h,d)`.  Use its dyadically cleared tridiagonal matrix only as one entry of the rank test; the actual carrier is the two-border/Smith gcd above.