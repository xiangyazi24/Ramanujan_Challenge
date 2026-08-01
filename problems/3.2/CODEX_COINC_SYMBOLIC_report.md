# Difference-curve symbolic factorization report

All computations below were exact, over `QQ`, with SymPy 1.14.0.  I used

```text
D_h(x) = prod_{j=1}^h (x+j)^3,
delta_h(x) = N_h(x)/D_h(x),
F_{h,k}(x,y) = N_h(x)D_k(y)-N_k(y)D_h(x).
```

Every call to `factor_list` in S1 and S2 completed over `QQ`; none of the
answers below is merely mod-p evidence.  Total degree is denoted by `tdeg`,
and bidegree in `(x,y)` by `bideg`.

## S1 - Same-gap curves

The expected diagonal factor is always present:

```text
F_{h,h}(x,y) = u_h (x-y) H_h(x,y),
```

where `H_h` is the primitive nonlinear factor returned by SymPy.  Both
`x-y` and `H_h` are irreducible over `QQ`.  The complete factor data are:

- `h=2`: `u_h=-1`; factor degrees `1` and `8`; nonlinear bidegree `(5,5)`;
  `H_h` has 33 terms.
- `h=3`: `u_h=-1`; factor degrees `1` and `14`; nonlinear bidegree `(8,8)`;
  `H_h` has 78 terms.
- `h=4`: `u_h=-2`; factor degrees `1` and `20`; nonlinear bidegree `(11,11)`;
  `H_h` has 141 terms.
- `h=5`: `u_h=-1`; factor degrees `1` and `26`; nonlinear bidegree `(14,14)`;
  `H_h` has 222 terms.
- `h=6`: `u_h=-1`; factor degrees `1` and `32`; nonlinear bidegree `(17,17)`;
  `H_h` has 321 terms.
- `h=7`: `u_h=-1`; factor degrees `1` and `38`; nonlinear bidegree `(20,20)`;
  `H_h` has 438 terms.
- `h=8`: `u_h=-4`; factor degrees `1` and `44`; nonlinear bidegree `(23,23)`;
  `H_h` has 573 terms.

Thus, uniformly in the tested range,

```text
tdeg(H_h) = 6h-4,       bideg(H_h) = (3h-1,3h-1),
c(h,h) = 2.
```

The exact nonlinear factor, without a many-thousand-character coefficient
dump, is unambiguously displayed by

```text
H_h = primitive_part_QQ((N_h(x)D_h(y)-N_h(y)D_h(x))/(x-y)),
```

with the signs and contents `u_h` listed above.  For a fully expanded sanity
check, the smallest nonlinear factor is

```text
H_2 =
  34*x^5*y^3 + 153*x^5*y^2 + 231*x^5*y + 117*x^5
+ 34*x^4*y^4 + 459*x^4*y^3 + 1608*x^4*y^2 + 2196*x^4*y + 1053*x^4
+ 34*x^3*y^5 + 459*x^3*y^4 + 2730*x^3*y^3 + 7245*x^3*y^2
+ 8676*x^3*y + 3861*x^3
+ 153*x^2*y^5 + 1608*x^2*y^4 + 7245*x^2*y^3 + 16071*x^2*y^2
+ 17190*x^2*y + 7099*x^2
+ 231*x*y^5 + 2196*x*y^4 + 8676*x*y^3 + 17190*x*y^2
+ 16837*x*y + 6498*x
+ 117*y^5 + 1053*y^4 + 3861*y^3 + 7099*y^2 + 6498*y + 2364.
```

There is no mirror-line factor.  The reflection identities give

```text
N_h(-x-h-1) = (-1)^(h-1) N_h(x),
D_h(-x-h-1) = (-1)^h D_h(x),
delta_h(-x-h-1) = -delta_h(x).
```

Consequently, on the candidate mirror line `y=-x-h-1`,

```text
F_{h,h}(x,-x-h-1) = 2*(-1)^h N_h(x)D_h(x),
```

which is not the zero polynomial.  Hence there is no constant `c_h` for
which `x+y+c_h` is a component.  The geometric reflection line would have
`c_h=h+1`, but it reverses the sign of `delta_h`; it only meets the
same-value curve in finitely many zero/pole points.

For completeness, the unrequested boundary case `h=1` is exceptional:

```text
F_{1,1} = (y-x)*((y+1)^2+(x+1)*(y+1)+(x+1)^2).
```

The quadratic factor is irreducible over `QQ` but splits into two lines over
`Qbar` by adjoining the nontrivial cube roots of unity.  This is why the
nontrivial S1 scan naturally starts at `h=2`.  Over `F_p` (away from `p=3`),
this boundary curve has three `p`-sized lines when `p=1 mod 3`, but only the
diagonal is `F_p`-defined when `p=2 mod 3`; the conjugate pair then has only
its common pole point over `F_p` and supplies no `p` main term.

## S2 - Distinct-gap curves

All 28 requested distinct-gap polynomials are irreducible in `QQ[x,y]`.
There are no exceptional factors.  In every case

```text
tdeg(F_{h,k}) = 3(h+k-1),       bideg(F_{h,k}) = (3h,3k).
```

The complete tested list, grouped by the smaller gap, is:

- `h=1`: `(1,2), (1,3), (1,4), (1,5), (1,6), (1,7), (1,8)` are irreducible,
  of total degrees `6,9,12,15,18,21,24`.
- `h=2`: `(2,3), (2,4), (2,5), (2,6), (2,7), (2,8)` are irreducible,
  of total degrees `12,15,18,21,24,27`.
- `h=3`: `(3,4), (3,5), (3,6), (3,7), (3,8)` are irreducible,
  of total degrees `18,21,24,27,30`.
- `h=4`: `(4,5), (4,6), (4,7), (4,8)` are irreducible,
  of total degrees `24,27,30,33`.
- `h=5`: `(5,6), (5,7), (5,8)` are irreducible,
  of total degrees `30,33,36`.
- `h=6`: `(6,7), (6,8)` are irreducible, of total degrees `36,39`.
- `h=7`: `(7,8)` is irreducible, of total degree `42`.

The parity factor of an even-gap numerator does not induce a curve
component.  If `h` is even, write

```text
N_h(x) = (2x+h+1) M_h(x).
```

At `x=-(h+1)/2`, the first term of `F_{h,k}` vanishes, but the remaining
term is

```text
-N_k(y) D_h(-(h+1)/2),
```

which is nonzero as a polynomial in `y`.  Thus `2x+h+1` is not a factor.
The analogous statement holds in the `y` coordinate.  When both gaps are
even, the two forced roots give one point on `X_{h,k}`, not a one-dimensional
component.  Exact factorization also rules out every mixed rational line of
the suggested form involving `2x+h+1` and `2y+k+1`.

## S3 - Genus and singularity budget

The required plane-curve arithmetic-genus bound for total degree `e` is

```text
g_a = (e-1)(e-2)/2.
```

In this range one can cheaply determine more than that.  Put

```text
q_h(x) = prod_{j=1}^h (x+j),
A_h(x) = N_h'(x)q_h(x)-3N_h(x)q_h'(x),
C_h(T) = Res_x(N_h(x)-T*q_h(x)^3, A_h(x)).
```

Thus `A_h/q_h^4` is the derivative of `delta_h`.  Exact `QQ` calculations
for every `1 <= h <= 8` gave

```text
gcd(N_h,N_h') = gcd(N_h,q_h) = 1,
gcd(A_h,A_h') = 1,
deg(C_h) = 4h-4,
gcd(C_h,C_h') = 1,
gcd(C_h,C_k) = 1 for every h != k.
```

For `h=1`, `A_1=-3` and `C_1` is constant, as the displayed degree formula
indicates.  These certificates say that all nonzero finite branch values of
the maps `delta_h : P1 -> P1` are simple and that two distinct tested maps
share no such branch value.  Their only shared branch values are `0` and
`infinity`.

This gives the complete singular-point count of the ordinary plane
projective closure.  The affine singularities are exactly the `hk` pole
pairs

```text
(x,y)=(-i,-j),       1 <= i <= h, 1 <= j <= k.
```

There are no regular affine singularities, because such a point would give
a common root of `C_h` and `C_k`.  If `a_h=lc(N_h)`, the equation at the line
at infinity is

```text
X^(3h-3) Y^(3k-3) (a_h*Y^3-a_k*X^3) = 0.
```

The three noncoordinate roots of the cubic factor are smooth.  The point
`[1:0:0]` is singular for every requested pair, and `[0:1:0]` is also
singular when `h >= 2`.  Therefore the number of distinct singular points
over `Qbar` is

```text
hk+1  if h=1,
hk+2  if h>=2.
```

The same branch certificates prove absolute irreducibility of every
distinct-gap curve in the scan.  Indeed, a nontrivial common intermediate
cover of `delta_h` and `delta_k` could branch only over their common branch
set `{0,infinity}`.  Over `Qbar`, such a cover is a power map.  This would
force both rational functions to be nontrivial perfect powers, which their
simple numerator zeros exclude.  For `h=1`, `delta_1` is a cube, but
`delta_k` for `k>=2` is not, so the conclusion is unchanged.

One can then also compute the exact geometric genus of the normalization.
The degree over the common value line is `(3h)(3k)=9hk`.  The ramification
profiles of `delta_h` are

```text
over 0:         one index-3 point and 3h-3 unramified points;
over infinity: h index-3 points;
elsewhere:     4h-4 distinct simple branch values.
```

The total ramification of the normalized fiber product is

```text
R = 3k(4h-4) + 3h(4k-4) + (6h+6k-6) + 6hk
  = 30hk-6h-6k-6.
```

Riemann-Hurwitz therefore gives, for every requested distinct pair,

```text
g = 6hk-3h-3k-2.
```

The requested per-pair budgets are listed next as
`(h,k): degree, g_a upper bound, singular-point count, exact g`:

- `(1,2): 6, 10, 3, 1`; `(1,3): 9, 28, 4, 4`;
  `(1,4): 12, 55, 5, 7`; `(1,5): 15, 91, 6, 10`;
  `(1,6): 18, 136, 7, 13`; `(1,7): 21, 190, 8, 16`;
  `(1,8): 24, 253, 9, 19`.
- `(2,3): 12, 55, 8, 19`; `(2,4): 15, 91, 10, 28`;
  `(2,5): 18, 136, 12, 37`; `(2,6): 21, 190, 14, 46`;
  `(2,7): 24, 253, 16, 55`; `(2,8): 27, 325, 18, 64`.
- `(3,4): 18, 136, 14, 49`; `(3,5): 21, 190, 17, 64`;
  `(3,6): 24, 253, 20, 79`; `(3,7): 27, 325, 23, 94`;
  `(3,8): 30, 406, 26, 109`.
- `(4,5): 24, 253, 22, 91`; `(4,6): 27, 325, 26, 112`;
  `(4,7): 30, 406, 30, 133`; `(4,8): 33, 496, 34, 154`.
- `(5,6): 30, 406, 32, 145`; `(5,7): 33, 496, 37, 172`;
  `(5,8): 36, 595, 42, 199`.
- `(6,7): 36, 595, 44, 211`; `(6,8): 39, 703, 50, 244`.
- `(7,8): 42, 820, 58, 289`.

Thus the elementary Hasse-Weil constant `2g_a` is always safe, while the
branch calculation supplies the sharper constant `2g` for the smooth
projective normalization (apart from the usual `O(degree)` adjustment when
returning to the affine curve).

For the same-gap nonlinear factors from S1, the result `H_h` irreducible over
`QQ` means the arithmetic monodromy of `delta_h` is 2-transitive.  Since
`C_h` is squarefree and nonconstant for `h>=2`, the monodromy contains a
transposition.  Hence it is the full symmetric group `S_(3h)`; the geometric
monodromy also contains that transposition and is normal, so it is again
`S_(3h)`.  Consequently each `H_h` in S1 is absolutely irreducible, not just
irreducible over `QQ`.

## S4 - Exact update identity and its leverage

Extend the numerator sequence by `N_0=0`.  Start with

```text
G_{m+g}(s) = G_g(s+m) G_m(s)
```

and set

```text
s=x+1,       m=h,       g=d-1.
```

The second row of `G_{d-1}(x+h+1)` is

```text
[-N_{d-1}(x+h+1), N_d(x+h)],
```

and the second column of `G_h(x+1)` is

```text
[(x+h+1)^6 N_h(x), N_{h+1}(x)]^T.
```

Their scalar product is the bottom-right entry of `G_{h+d-1}(x+1)`, namely
`N_{h+d}(x)`.  This proves the exact addition formula

```text
N_{h+d}(x)
  = N_d(x+h) N_{h+1}(x)
    - (x+h+1)^6 N_{d-1}(x+h+1) N_h(x).            (U1)
```

Using instead `m=h-1` and `g=d` gives the equally useful companion identity

```text
N_{h+d}(x)
  = N_{d+1}(x+h-1) N_h(x)
    - (x+h)^6 N_d(x+h) N_{h-1}(x).                (U2)
```

Since

```text
D_{h+d}(y) = D_h(y) D_d(y+h),
```

substitution of (U1) into the definition of `F` gives the requested exact
lower-index identity:

```text
F_{h,h+d}(x,y)
  = N_h(x) D_h(y) D_d(y+h)
    - D_h(x) N_d(y+h) N_{h+1}(y)
    + D_h(x) (y+h+1)^6 N_{d-1}(y+h+1) N_h(y).     (F-update)
```

Every `N` on the right has index in `{h,h+1,d,d-1}`, rather than `h+d`.
After normalizing by denominators, the cleanest cocycle form is

```text
delta_{h+d}(x)
  = (x+h+1)^3 *
    (delta_{h+1}(x) delta_d(x+h)
     - delta_h(x) delta_{d-1}(x+h+1)).             (D-update)
```

Thus the curve itself can be written as

```text
delta_h(x)
  = (y+h+1)^3 *
    (delta_{h+1}(y) delta_d(y+h)
     - delta_h(y) delta_{d-1}(y+h+1)).
```

For the required check `(h,d)=(3,2)`, (U1) specializes to

```text
N_5(x) = N_2(x+3) N_4(x) - (x+4)^6 N_3(x),
N_2(x+3) = 34*x^3+459*x^2+2067*x+3105 = P(x+4).
```

An exact SymPy expansion returned

```text
update_h3_d2_zero = True
F_h3_h5_three_term_zero = True
alternate_update_h3_d2_zero = True
```

As a wider index audit, exact expansion also verified (U1), (U2), and
(F-update) in all 36 cases with `h>=1`, `d>=1`, and `h+d<=9`; all residuals
were zero.  The reflection-line substitution was independently checked for
all `1 <= h <= 8`.

The identity is genuine self-similarity, but at the rank-2 state level.  A
long gap `h+d` is obtained from the boundary state `(N_h,N_{h+1})` and a
shifted length-`d` transfer state `(N_{d-1},N_d)`.  It therefore offers a
block/induction mechanism if one propagates the two-component state (or its
projective ratio).  It does not by itself identify `X_{h,h+d}` with a single
lower-gap difference curve, and it does not force a factorization: the
shift, the boundary state at `h`, and the bilinear determinant remain
coupled.  The complete irreducibility found in S2 is consistent with this
more limited but still useful cocycle self-similarity.

## S5 - Verdict

For distinct gaps, the evidence is uniform and in the tested range is stronger than heuristic: all 28 curves with `1 <= h < h' <= 8` are irreducible over `QQ`, their nontrivial critical-value sets are pairwise disjoint, and the common-cover argument makes them absolutely irreducible; this strongly suggests that generic distinct-gap curves are absolutely irreducible whenever the same squarefreeness and branch-disjointness persist.  The complete systematic component list found here is: for every same gap `2 <= h <= 8`, the diagonal `x=y` and one absolutely irreducible off-diagonal component `H_h=0`, each contributing one `p+O_h(sqrt(p))` main term and hence about `2p` together; for distinct gaps, one absolutely irreducible component, contributing about `p`; and, only in the boundary case `(1,1)`, the diagonal plus a `QQ`-irreducible quadratic that geometrically splits into the two cube-root-of-unity lines, giving about `3p` when `p=1 mod 3` and only the diagonal `p` main term when `p=2 mod 3`.  There is no mirror component: reflection sends `delta_h` to `-delta_h`.  The forced linear factor of even `N_h` creates isolated zero-fiber points, not curve components, so it contributes no additional `p` main term.  These component main terms concern full `F_p` curves; restriction to the collision windows still requires the separate incomplete-point-count analysis.
