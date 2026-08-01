# Irreducibility and branch-disjointness report

## 1. Verdict

Status labels in this report have the following meanings.

- PROVED: a symbolic argument valid for every indicated index.
- VERIFIED-16: an exact computer-assisted proof for every indicated index up
  to 16. This is not numerical sampling.
- OPEN: not proved for unbounded indices.

The monodromy route can be completed from substantially less than full Morse
genericity. The missing primitivity step does not follow merely from
"transitive + a 3-cycle + a transposition"; that group-theoretic assertion is
false in general. For the present covers, however, `2h-1` simple nonzero
critical values rule out every decomposition of the rational function. This
supplies primitivity and gives the following uniform conditional theorem.

Let

```text
q_h(x) = product_(j=1)^h (x+j),
D_h(x) = q_h(x)^3,
delta_h(x) = N_h(x)/D_h(x),
A_h(x) = N_h'(x) q_h(x) - 3 N_h(x) q_h'(x),
C_h(T) = Res_x(N_h(x)-T q_h(x)^3, A_h(x)).
```

Factor the primitive part of `C_h` over Qbar as

```text
C_h(T) = product_a (T-a)^m(a).
```

Let

```text
s_h = #{a != 0 : m(a)=1},
CRIT(h): s_h >= 2h-1.
```

The stronger full-Morse condition used by the finite certificates is

```text
BG(h): C_h(0) != 0 and gcd(C_h,C_h')=1.
```

It gives `s_h=4h-4`, hence implies CRIT(h).

The main conclusions are:

1. PROVED (all-h conditional monodromy). For every h >= 2, CRIT(h) implies

```text
GeomMon(delta_h) = S_(3h).
```

2. PROVED (T3 conditional). For every h >= 2, CRIT(h) implies that

```text
H_h(x,y) = F_(h,h)(x,y)/(x-y)
```

is irreducible over Qbar.

3. PROVED (T1 conditional). If 2 <= h < k and CRIT(h), CRIT(k) hold, then
`F_(h,k)` is irreducible over Qbar. If h=1<k, CRIT(k) alone suffices.

4. PROVED (exact characterization of T2). Put

```text
C_h^*(T) = C_h(T)/T^ord_T(C_h).
```

The nonzero finite branch loci are disjoint exactly when

```text
gcd(C_h^*,C_k^*) = 1.
```

A failure is exactly a common nonzero critical value. Equivalently, it is an
affine singularity of `F_(h,k)=0` with both coordinates away from the poles.
Such a failure need not make the fiber product reducible.

5. VERIFIED-16 (unconditional finite range). Exact certificates at two good
primes prove BG(h) for every `2 <= h <= 16` and prove
`gcd(C_h,C_k)=1` for all 105 pairs `2 <= h < k <= 16`. Consequently:

```text
T1 holds for every 1 <= h < k <= 16.
T2 holds for every 1 <= h < k <= 16.
T3 holds for every 2 <= h <= 16.
GeomMon(delta_h) = S_(3h) for every 2 <= h <= 16.
```

6. OPEN (unbounded indices). I did not prove CRIT(h) for every h, nor the
required pairwise coprimality for every distinct pair. Thus T1-T3 are not
claimed unconditionally for all h. The weakest explicit sufficient statements
isolated by this route are

```text
[CRIT-2H] s_h >= 2h-1 for every h >= 2.
[BD-ALL]  Res_T(C_h^*,C_k^*) != 0 for every 2 <= h < k.
```

CRIT-2H closes T1 and T3 by the proof below. BD-ALL closes T2. The stronger
full-Morse statement BG(h) for every h would imply CRIT-2H.

## 2. Elementary structure of the covers

The recurrence gives

```text
deg N_h = 3h-3,                 deg q_h^3 = 3h,
gcd(N_h,q_h) = 1.
```

The last assertion also follows directly from the supplied polar-value
identity. For `1 <= j <= h`, put

```text
r_(h,j) = product_(1 <= m <= h, m != j) (m-j)
        = (-1)^(j-1) (j-1)! (h-j)!.
```

The supplied Apery-product evaluation becomes

```text
N_h(-j) = r_(h,j)^3 b_(j-1) b_(h-j).
```

Hence `N_h(-j) != 0` and the polar part is especially simple:

```text
delta_h(x)
  = b_(j-1)b_(h-j)/(x+j)^3 + O((x+j)^(-2)).       (2.1)
```

If `ell_h=lc(N_h)`, then `ell_1=1` and

```text
ell_(h+1) = 34 ell_h - ell_(h-1),
```

so `ell_h != 0`. At infinity,

```text
delta_h(x) = ell_h/x^3 + O(x^(-4)).                (2.2)
```

The reflection identity is

```text
N_h(-x-h-1) = (-1)^(h-1) N_h(x),
q_h(-x-h-1) = (-1)^h q_h(x),
delta_h(-x-h-1) = -delta_h(x).                     (2.3)
```

In particular, finite critical points occur in reflected pairs and their
critical values are negatives of one another. Thus `C_h(T)` is even, up to a
nonzero scalar. This symmetry does not itself prove squarefreeness.

Direct differentiation gives

```text
delta_h'(x) = A_h(x)/q_h(x)^4.                     (2.4)
```

The leading term of `A_h` is `-3 ell_h x^(4h-4)`, so

```text
deg A_h = 4h-4.
```

At a zero `-j` of `q_h`,

```text
A_h(-j) = -3 N_h(-j) q_h'(-j) != 0.
```

Therefore `gcd(A_h,q_h)=1` and

```text
deg_T C_h = 4h-4.                                  (2.5)
```

This degree statement is unconditional.

## 3. The critical-value criterion

PROVED. Assume BG(h). Then the branch profile of

```text
delta_h : P1_x -> P1_T
```

is exactly

```text
over T=0:        (3,1^(3h-3));
over T=infinity: (3^h);
over 4h-4 other, pairwise distinct values:
                  (2,1^(3h-2)).                    (3.1)
```

Proof. The roots of `C_h` are the values `delta_h(alpha)` with
`A_h(alpha)=0`. Because `A_h` is coprime to `q_h`, no pole contributes to
this resultant. A repeated root of `A_h`, or two distinct roots of `A_h`
having the same critical value, would give a repeated factor of `C_h`.
Thus `gcd(C_h,C_h')=1` says that the `4h-4` finite critical points are
simple and have pairwise distinct values.

Also

```text
C_h(0) = 0  iff  gcd(N_h,A_h) != 1
         iff  gcd(N_h,N_h') != 1,
```

where the last equivalence uses `gcd(N_h,q_h)=1` and
`A_h=N_h' q_h` at a root of `N_h`. Hence `C_h(0) != 0` says that all
`3h-3` finite numerator roots are simple. Equation (2.2) supplies the one
index-3 zero at infinity, and (2.1) supplies the h index-3 poles. The simple
zeros of `A_h` give the remaining simple ramification points. Finally,

```text
2 + 2h + (4h-4) = 6h-2 = 2 deg(delta_h)-2,
```

so Riemann-Hurwitz leaves no unlisted ramification. This proves (3.1).

Full BG(h) is not needed to compute CRIT(h). In characteristic zero, put

```text
g_h      = gcd(C_h,C_h'),
rad_h    = C_h/g_h,
repeat_h = gcd(rad_h,g_h),
simple_h = rad_h/repeat_h.
```

Remove one factor T from `simple_h` if present. The degree of the result is
exactly `s_h`: a factor appears in `simple_h` precisely when its multiplicity
in `C_h` is one. This is the one-variable [CRIT-2H] test; it remains useful
even if full Morse squarefreeness first fails.

## 4. Primitivity and full symmetric monodromy

### 4.1 The missing primitivity lemma

PROVED. If h >= 2 and CRIT(h), then `delta_h` is indecomposable over Qbar.

Suppose otherwise that

```text
delta_h = g o u,
deg u = a > 1,
deg g = b > 1,
ab = 3h.                                           (4.1)
```

For a value t, let `defect_f(t)` be the sum of `e_f(x)-1` over its fiber.
For each `z in g^(-1)(t)`, composition gives the contribution

```text
(e_g(z)-1)a + defect_u(z).                         (4.2)
```

A simple critical value of `delta_h` has total defect 1. Since `a>1`, (4.2)
forces g to be unramified over that value and forces exactly one simple
critical point of u above it. Distinct simple critical values give distinct
critical points of u. Hence

```text
s_h <= 2a-2.                                       (4.3)
```

If `b>=3`, then `a<=h`, so (4.3) gives

```text
s_h <= 2h-2,                                       (4.4)
```

contrary to CRIT(h).

It remains to exclude `b=2`. Here `a=3h/2`. A degree-2 outer map cannot have
a double pole: composition through it would give even local degree, whereas
every pole of `delta_h` has local degree 3 by (2.1). Thus g has two simple
poles. Above either pole, every point of the u-fiber has u-index 3. Each pole
fiber therefore consumes

```text
a-a/3 = 2a/3 = h
```

units of ramification of u, for `2h` units in total. In addition,
`x=infinity` is an index-3 zero of `delta_h`. The outer map cannot supply an
index 2 there, so u has index 3 and consumes two further units. Since the
total ramification of u is

```text
2a-2 = 3h-2,
```

at most `h-4` units remain for simple critical points. This is again less
than `2h-1`; when `h<4`, the forced ramification already makes the
decomposition impossible.

Therefore `delta_h` is indecomposable for every h>=2 under CRIT(h). The
threshold `2h-1` is exactly one above the general `b>=3` upper bound.

### 4.2 From indecomposable to S_(3h)

By the intermediate-field theorem and Luroth's theorem, indecomposability of
a rational function is equivalent to primitivity of its geometric monodromy
action. CRIT(h) supplies at least one literal transposition as local inertia
around a simple nonzero critical value.

A primitive permutation group containing a transposition is the full
symmetric group. For completeness, take the graph whose edges are the
supports of all conjugates of that transposition. Its connected components
form a block system. Primitivity makes the graph connected, and edge
transpositions of a connected graph generate the full symmetric group.

Consequently

```text
GeomMon(delta_h) = S_(3h).                          (4.5)
```

Notice what was and was not used. The index-3 zero and the index-3 poles are
needed to control the outer-degree-2 case. But the bare assertion
"transitive + a 3-cycle + a transposition" would not prove primitivity: an
imprimitive wreath product can contain all three kinds of elements.

## 5. Fiber products

### 5.1 Same gap: T3

Let `n=3h`. The irreducible components of the normalized generic fiber product

```text
delta_h(x) = delta_h(y)
```

correspond to the orbits of geometric monodromy on ordered pairs of sheets.
The diagonal action of `S_n` has exactly two orbits:

```text
{(i,i)}
and
{(i,j): i != j}.
```

Thus the diagonal is one component and the off-diagonal locus is one
geometrically irreducible component. Since `gcd(N_h,D_h)=1`, clearing
denominators adds no vertical or horizontal component. Moreover `x-y` occurs
with multiplicity one, as can be checked at any unramified diagonal point.
Therefore

```text
F_(h,h)(x,y) = (x-y) H_h(x,y)
```

with `H_h` irreducible over Qbar. This proves T3 under CRIT(h).

The excluded h=1 case behaves exactly as previously found: `delta_1` is a
degree-3 power map and the off-diagonal quadratic splits over Qbar into the
two cube-root-of-unity lines.

### 5.2 Distinct gaps: T1

First take `2 <= h < k`, and put `n=3h`, `m=3k`. Under CRIT(h), CRIT(k), the
two geometric Galois closures have groups `S_n` and `S_m`. Their intersection
over `Qbar(T)` gives a common quotient. Since `n,m>=6` and `n!=m`, the only
possibilities are

```text
the trivial quotient,
or the common C2 sign quotient.                     (5.1)
```

In the first case the compositum group is `S_n x S_m`, which is transitive on
the `n*m` ordered pairs of sheets.

In the second case it is

```text
{(sigma,tau) in S_n x S_m : sign(sigma)=sign(tau)}. (5.2)
```

This subgroup is also transitive on ordered pairs. To send `(i,j)` to
`(i',j')`, choose permutations doing so in both coordinates. Their parities
can be changed independently by multiplying by odd permutations in the two
point stabilizers; these exist because `n,m>=3`. Thus the parities can be
made equal.

The normalized fiber product is therefore connected in either case. Clearing
denominators again adds no vertical or horizontal component, so `F_(h,k)` is
irreducible over Qbar.

If h=1, the first cover is cyclic of degree 3. It has no nontrivial common
quotient with `S_(3k)` for k>=2, since `S_(3k)` has no quotient `C3`. The
product action is transitive, giving T1 in this boundary case as well.

An important consequence is that branch-disjointness is not needed for T1.
Even a common sign field does not split the distinct-degree fiber product.

## 6. Branch-disjointness and its exact failure mode

Without any simplicity assumption, the nonzero finite branch values of
`delta_h` are exactly the roots of `C_h^*`. Therefore, for h != k,

```text
Branch(delta_h) intersect Branch(delta_k)
  is contained in {0,infinity}

iff

gcd(C_h^*,C_k^*)=1.                                (6.1)
```

If (6.1) fails, a common root t gives critical points alpha, beta with

```text
A_h(alpha)=0,
A_k(beta)=0,
delta_h(alpha)=delta_k(beta)=t.
```

After clearing denominators, this is exactly

```text
F_(h,k)(alpha,beta)
= partial_x F_(h,k)(alpha,beta)
= partial_y F_(h,k)(alpha,beta)
= 0,                                                (6.2)
```

with neither coordinate a pole. Thus failures of T2 are precisely the
non-pole affine singularities (6.2). Under BG(h), BG(k), the two critical
points are unique and simple, and the singularity is an ordinary node. Such
singularities do not by themselves create extra components; T1 is still
controlled by the monodromy argument in Section 5.

For h=1, `delta_1=1/(x+1)^3` has no nonzero finite branch value, so T2 is
automatic against every k>=2.

## 7. Local pole calculation

Put `x=-j+u`, `y=-i+v`. The first nonzero homogeneous part of
`F_(h,k)(x,y)` has degree 3 and is

```text
r_(h,j)^3 r_(k,i)^3
  * (b_(j-1)b_(h-j) v^3 - b_(i-1)b_(k-i) u^3).     (7.1)
```

This follows immediately from (2.1), or directly from the cleared equation.
The verification script checked every coefficient below total degree 3 and
every degree-3 coefficient at all pole pairs for

```text
(h,k)=(2,3):  2*3 = 6 checks,
(h,k)=(3,5):  3*5 = 15 checks.
```

All 21 checks were exact over ZZ.

Over Qbar the cubic tangent cone (7.1) splits into three lines. Consequently,
the pole-local Newton polygon identifies the three local branches but does not
by itself forbid a global factor from grouping them. The leading constants
are fully controlled, but a compatibility theorem transporting a proposed
grouping between all pole pairs would still be required. The monodromy proof
avoids this unresolved global transport step.

## 8. Exact certificates through h=16

The verifier is

```text
CODEX_irred_verify.py
```

and runs with SymPy over ZZ and finite prime fields. The command is

```text
PYTHONDONTWRITEBYTECODE=1 python3 -u CODEX_irred_verify.py
```

### 8.1 Why the modular test proves a characteristic-zero statement

For a fixed h, `C_h` has the unconditional degree `4h-4`. The script evaluates
the resultant at `4h-3` nonzero values of T and interpolates it modulo p. The
point T=0 is intentionally not used for interpolation: there the degree in x
of `N_h-Tq_h^3` drops by three, and specializing before taking the resultant
would omit a power of `lc(A_h)`. Three additional nonzero T-values check the
interpolation, and h=2,3,4 are independently compared with direct resultants
over ZZ.

For each good prime, the script verifies

```text
deg(C_h mod p) = 4h-4,
C_h(0) mod p != 0,
gcd(C_h mod p, C_h' mod p) = 1,
C_h(T) mod p is even,
simple-nonzero degree = 4h-4 >= 2h-1,
gcd(C_h mod p,C_k mod p) = 1 for all h<k.
```

Degree preservation is crucial. If an integer polynomial had a repeated
factor over Q, its degree-preserving reduction would still have a repeated
nonconstant factor. More generally, reduction cannot increase the total
degree of multiplicity-one factors, so a modular lower bound for `s_h` is also
a characteristic-zero lower bound. Likewise, a common Q-factor of
`C_h,C_k` would remain a common nonconstant factor after a degree-preserving
reduction. Therefore one good prime already proves the Q-statements; the
second prime is an independent implementation audit.

### 8.2 Certificate results

At p=1009:

```text
h=2,...,16: expected degree, nonzero C_h(0), squarefree C_h, and
            s_h=4h-4 all pass.
105/105 pairs: gcd(C_h,C_k)=1.
coefficient SHA256:
9b3e572e3656cd2b63543738808f967dd6bbedf159181198bac6fa0170d68d7e
```

At p=65537:

```text
h=2,...,16: expected degree, nonzero C_h(0), squarefree C_h, and
            s_h=4h-4 all pass.
105/105 pairs: gcd(C_h,C_k)=1.
coefficient SHA256:
4d7fdc3ec199c4b10d83d2f33d0d302d83b4833b31b380e528bab1cc51eb7dd3
```

The script also checks all 136 polar coefficients for h<=16, the reflection
identity for every h<=16, and the 21 local tangent cones in Section 7.

As a small direct audit, the primitive critical-value polynomial at h=2 is

```text
625 T^4 + 541064 T^2 + 22717712.
```

It is squarefree and has no zero root, in agreement with both modular runs.

### 8.3 Independent Frobenius certificates for h<=6

For `R_h(x)=N_h(x)-q_h(x)^3`, exact factorization over Q proves irreducibility
for h=2,...,6. Thus each specialized Galois group is transitive. The first
listed good prime gives cycle type `[1,3h-1]`; a transitive group containing
such a cycle is 2-transitive. The second listed prime gives exactly one
2-cycle and only odd remaining cycles. Raising that Frobenius element to the
lcm of its odd cycle lengths leaves a transposition. A 2-transitive group
containing a transposition is `S_(3h)`.

The exact factor-degree data are

```text
h=2: p=29  [1,5];   p=37  [1,2,3].
h=3: p=43  [1,8];   p=29  [1,2,3,3].
h=4: p=41  [1,11];  p=89  [1,2,9].
h=5: p=79  [1,14];  p=89  [1,1,2,11].
h=6: p=151 [1,17];  p=73  [1,1,1,2,13].
```

All modular factorizations are squarefree, so the displayed degrees are valid
Frobenius cycle types. The specialized `S_(3h)` embeds in the generic
arithmetic group, forcing that group to be `S_(3h)`. The geometric group is a
normal subgroup and contains the branch transposition, hence is also
`S_(3h)`. These calculations independently confirm the uniform branch-profile
proof for h<=6.

## 9. Why the all-h proof stops at CRIT-2H and BD-ALL

Differentiating the exact addition law is valid, but it does not give a closed
induction for critical points. In schematic form it produces

```text
delta'_(h+d)
 = 3(x+h+1)^2 * (bilinear delta terms)
 + (x+h+1)^3 * (four mixed delta/delta' terms).
```

A zero of the new derivative is therefore not inherited from a zero of one
lower derivative, and equality of two new critical values does not reduce to
adjacent coprimality of the N-polynomials. The polar product formula controls
the special fiber over infinity, but critical values are global data at the
roots of `A_h`.

The computations support the full-Morse statement BG(h), which is stronger
than CRIT-2H, as well as BD-ALL. A finite census proves neither unbounded
statement. Nor does full symmetric monodromy imply branch-disjointness: two
full-symmetric covers may share isolated branch values while remaining
linearly disjoint and while their fiber product stays irreducible.

Thus the exact frontier of the present route is:

```text
To finish T1 and T3 for all h by this route: prove CRIT-2H.
To finish T2 for all h: prove BD-ALL.
```

No axiom, genericity substitution, or unproved local-to-global grouping is
used in the conclusions labeled PROVED or VERIFIED-16.
