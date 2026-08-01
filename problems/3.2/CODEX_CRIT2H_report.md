# [CRIT-2H] follow-up report

## 1. Verdict

The labels used below are:

- PROVED: a symbolic argument valid for every stated h.
- VERIFIED-32: an exact finite-field certificate which proves the corresponding
  characteristic-zero statement for every `2 <= h <= 32`.
- HEURISTIC: structural evidence, not a proof.
- OPEN: not proved for unbounded h.

The main outcome is a complete all-h proof of Task A.1.

1. PROVED: `N_h` is squarefree over Q for every `h >= 2`.
2. PROVED: `C_h(0) != 0`, hence `ord_(T=0) C_h=0`, for every `h >= 2`.
   For the raw resultant there is the exact identity

   ```text
   C_h^raw(0)
     = -27 lambda_h^4 S_h^6 B_h^2 Disc(N_h),

   S_h = product_(i=0)^(h-1) i!,
   B_h = product_(i=0)^(h-1) b_i,
   lambda_h = lc(N_h).
   ```

3. VERIFIED-32: `C_h` is squarefree and nonzero at zero for all
   `2 <= h <= 32`. Thus `s_h=4h-4` and `[CRIT-2H]` holds throughout this
   range.
4. VERIFIED-32: `gcd(C_h,C_k)=1` for all 465 pairs
   `2 <= h < k <= 32`. Of these, 360 are new relative to the old range
   `h,k <= 16`.
5. OPEN: `[CRIT-2H]` for unbounded h. The reflection quotient gives a clean
   half-degree formulation, but neither Bezout nor the fiber-degree bound gives
   the required linear collision control.

There is also a correction to the proposed multiplicity-cap route. The two
conditions

```text
deg gcd(C_h,C_h') <= 2h-3,
ord_(T=0) C_h = 0
```

do not imply the required number of multiplicity-one roots. A correct simple
sufficient bound is

```text
deg gcd(C_h,C_h') <= h-2.
```

The exact quotient criterion is given in Section 4.

## 2. Normalization

This report follows the notation of `CODEX_IRRED_THEOREM_report.md`:

```text
P(x)       = 34x^3+51x^2+27x+5,
N_0(x)     = 0,
N_1(x)     = 1,
N_(h+1)(x) = P(x+h)N_h(x)-(x+h)^6N_(h-1)(x),
q_h(x)     = product_(j=1)^h (x+j),
delta_h(x) = N_h(x)/q_h(x)^3,
A_h(x)     = q_h(x)N_h'(x)-3q_h'(x)N_h(x),
C_h^raw(T) = Res_x(N_h(x)-Tq_h(x)^3,A_h(x)).
```

Here

```text
deg N_h = 3h-3,
deg A_h = 4h-4,
lc(A_h) = -3 lambda_h.
```

Dividing `C_h^raw` by nonzero integer content does not change its zero set,
its order at zero, its squarefreeness, or any of the gcd conclusions.

The banked pole identity is

```text
N_h(-j)
 = (-1)^(j-1) ((j-1)!(h-j)!)^3 b_(j-1)b_(h-j) != 0,
1 <= j <= h.                                           (2.1)
```

All `b_i` are positive Apery numbers.

## 3. Task A.1: an all-h squarefreeness proof

### 3.1 The tridiagonal pencil

Put `m=h-1`. For `0 <= a <= 1`, let `J_(h,a)(x)` be the real symmetric
`m` by `m` tridiagonal matrix with

```text
diagonal entry i:       P(x+i),
off-diagonal entry i:   a (x+i+1)^3,
1 <= i <= m-1.
```

Expansion along the last row gives

```text
det J_(h,1)(x)=N_h(x),
det J_(h,0)(x)=product_(i=1)^(h-1) P(x+i).              (3.1)
```

The leading coefficient of `det J_(h,a)` is the determinant of the constant
tridiagonal matrix with diagonal 34 and off-diagonal a. This matrix is
strictly diagonally dominant for `0 <= a <= 1`, so the leading coefficient is
positive. Hence the degree stays equal to `3h-3` throughout the deformation.

### 3.2 A scalar localization lemma

PROVED. If

```text
|P(z)| <= |z|^3+|z+1|^3,                               (3.2)
```

then

```text
-1 < Re(z) < 0.                                        (3.3)
```

It is enough to exclude `Re(z)>=0`, since

```text
P(-1-z)=-P(z)
```

and the right side of (3.2) is invariant under `z -> -1-z`.

Write `z=u+iv`, `r=|z|`, and `s=|z+1|`. For `u>=0`, one has `s>=r` and

```text
|2z+1| >= s.
```

Use the exact factorization

```text
P(z)=(2z+1)(17z(z+1)+5).
```

The following difference is strictly positive for `u>=0`:

```text
|17z(z+1)+5|^2-(|z|^2+|z+1|^2)^2
 = 285u^4+570u^3+570u^2v^2+451u^2
   +570uv^2+166u+285v^4+115v^2+24.                    (3.4)
```

Therefore

```text
|P(z)|
 > s(r^2+s^2)
 >= r^3+s^3,
```

which contradicts (3.2). This proves the lemma.

If `det J_(h,a)(z)=0`, take a nonzero kernel vector and choose a coordinate
of maximal absolute value. Its row equation gives

```text
|P(z+i)| <= a(|z+i|^3+|z+i+1|^3)
           <= |z+i|^3+|z+i+1|^3
```

for some `1 <= i <= h-1`. The lemma therefore confines every root to one of
the disjoint strips

```text
S_i = {z : -i-1 < Re(z) < -i},
1 <= i <= h-1.                                         (3.5)
```

No root can cross a boundary of these strips during `0 <= a <= 1`. The
leading coefficient stays nonzero, and the coefficients vary continuously,
so roots cannot escape to infinity either. At `a=0`, the three roots in
`S_i` are exactly

```text
-i-1/2,
-i-1/2 + i sqrt(51)/34,
-i-1/2 - i sqrt(51)/34.
```

It follows by the argument principle, or equivalently continuity of the
multiset of polynomial roots, that

```text
N_h has exactly three roots counted with multiplicity in every S_i. (3.6)
```

### 3.3 Every real root is simple

The signs in (2.1) alternate at consecutive negative integers. Hence every
strip `S_i` contains at least one real root.

At a real root `r`, none of the off-diagonal entries of `J_(h,1)(r)` vanishes:
the only possible zeros are negative integers, and (2.1) says those are not
roots of `N_h`. Thus the unreduced tridiagonal matrix `J_(h,1)(r)` has a
one-dimensional kernel.

Its derivative has diagonal `P'(r+i)` and off-diagonal
`3(r+i+1)^2`. It is positive definite. Indeed, for an interior row, with
`t=r+i`, the diagonal-dominance margin is

```text
P'(t)-3t^2-3(t+1)^2 = 24(2t+1)^2 >= 0.                (3.7)
```

The endpoint margins are strictly positive:

```text
P'(t)-3(t+1)^2 = 3(33t^2+32t+8) > 0,
P'(t)-3t^2     = 3(33t^2+34t+9) > 0.                  (3.8)
```

The graph is connected, so irreducible diagonal dominance makes the symmetric
matrix positive definite. If `v` spans the real kernel of `J_(h,1)(r)`, then
rank `J_(h,1)(r)=m-1` and

```text
adj(J_(h,1)(r))=c vv^T
```

for some real `c!=0`. Jacobi's determinant formula now gives

```text
N_h'(r)
 = tr(adj(J_(h,1)(r)) J_(h,1)'(r)) != 0.              (3.9)
```

Thus every real root is simple.

### 3.4 Nonreal repeated roots are impossible

Suppose a nonreal root `alpha` had multiplicity at least two. Since `N_h` has
real coefficients, `conj(alpha)` would be a distinct root of the same
multiplicity. Both have the same real part, hence lie in the same strip. They
would contribute at least four roots counted with multiplicity to that strip,
contradicting (3.6).

Therefore:

```text
gcd(N_h,N_h')=1 for every h>=2.                         (3.10)
```

This proof includes the forced central root for even h and shows that it is
simple. It does not use, or assume, all-h irreducibility of the symmetrized
polynomial `M_h`.

### 3.5 Exact constant-term identity

Let `n=3h-3`, `k=4h-4`, and `lambda_h=lc(N_h)`. At `T=0`, the first argument
of the resultant drops in degree from `3h` to `3h-3`. Using the root-product
definition of the resultant with the roots of `A_h` gives the exact
specialization factor

```text
C_h^raw(0)=(-3lambda_h)^3 Res(N_h,A_h).                 (3.11)
```

At a root of `N_h`,

```text
A_h=q_h N_h'.
```

The degrees match, so multiplicativity of the resultant gives

```text
Res(N_h,A_h)=Res(N_h,q_h)Res(N_h,N_h').                 (3.12)
```

Put

```text
S_h=product_(i=0)^(h-1) i!,
B_h=product_(i=0)^(h-1) b_i.
```

Equation (2.1) yields

```text
Res(N_h,q_h)=(-1)^(h(h-1)/2) S_h^6 B_h^2.              (3.13)
```

With the standard discriminant convention,

```text
Res(N_h,N_h')
 = (-1)^(n(n-1)/2) lambda_h Disc(N_h).                 (3.14)
```

The sum of the two sign exponents is

```text
h(h-1)/2+n(n-1)/2=(h-1)(5h-6),
```

which is always even. Combining (3.11)--(3.14) gives

```text
C_h^raw(0)
 = -27 lambda_h^4 S_h^6 B_h^2 Disc(N_h).               (3.15)
```

Every displayed factor is nonzero: `S_h,B_h>0`, `lambda_h>0` from the
leading tridiagonal matrix, and `Disc(N_h)!=0` by (3.10). Hence

```text
C_h(0)!=0 and ord_(T=0) C_h=0 for every h>=2.           (3.16)
```

The raw identity (3.15) was also checked directly over ZZ for `h=2,3,4`.

## 4. Task A.2: reflection quotient and multiplicity bounds

### 4.1 Correction to the proposed gcd threshold

Let the roots of a degree `4h-4` polynomial C have multiplicities `r_a`.
Then

```text
deg gcd(C,C')=sum_a (r_a-1),
```

whereas the degree carried by non-simple roots is

```text
sum_(r_a>=2) r_a.
```

These are not the same quantity.

For a concrete counterprofile, take distinct nonzero numbers
`a_1,...,a_(h-2),b_1,b_2` and put

```text
D(Z)=product_(i=1)^(h-2)(Z-a_i)^2 (Z-b_1)(Z-b_2),
C(T)=D(T^2).
```

Then `C(0)!=0` and

```text
deg gcd(C,C')=2h-4 <= 2h-3,
```

but C has only four multiplicity-one roots. This fails `[CRIT-2H]` for every
`h>=3`. Thus the proposed cap does not prove the stated target.

### 4.2 The exact quotient

Put

```text
s=2x+h+1,
U=s^2.
```

The alternative variable `x(x+h+1)` differs from U by an affine change.
Reflection sends `s` to `-s`. Since

```text
delta_h(-h-1-x)=-delta_h(x),
A_h(-h-1-x)=A_h(x),
```

there are rational functions `R_h(U), Phi_h(U)` and a polynomial `J_h(U)`
such that

```text
delta_h(x)=s R_h(U),
Phi_h(U)=delta_h(x)^2=U R_h(U)^2,
A_h(x)=J_h(U),
deg J_h=2h-2,
deg Phi_h=3h.                                          (4.1)
```

Section 3 proves that no critical point lies over `delta_h=0`. The center is
either a noncritical simple zero (even h) or a pole (odd h). Thus the
`4h-4` finite critical points form `2h-2` free reflection pairs, exactly the
roots of `J_h` counted with multiplicity.

Write

```text
C_h(T)=c_h D_h(T^2),
deg D_h=2h-2.
```

Because `D_h(0)!=0`,

```text
deg gcd(C_h,C_h')=2 deg gcd(D_h,D_h').                 (4.2)
```

The multiplicity of a root z of `D_h` is the total critical defect among the
quotient critical points with squared value `Phi_h(U)=z`. In particular,
reflection creates the pair of values `t,-t`; it does not create a
same-value collision because `t!=0`.

Let

```text
g_D    = gcd(D_h,D_h'),
rad_D  = D_h/g_D,
rep_D  = gcd(rad_D,g_D),
simp_D = rad_D/rep_D.
```

Then the exact relation is

```text
s_h=2 deg(simp_D).                                      (4.3)
```

Since the left side is even, `[CRIT-2H]` is equivalent to

```text
deg(simp_D) >= h.                                       (4.4)
```

If `e=deg g_D`, every non-simple root of multiplicity r uses r degree while
contributing `r-1` to e, and `r <= 2(r-1)`. Therefore

```text
deg(simp_D) >= (2h-2)-2e.                              (4.5)
```

A correct convenient sufficient cap is consequently

```text
e <= (h-2)/2,
```

or, by (4.2),

```text
deg gcd(C_h,C_h') <= h-2.                              (4.6)
```

### 4.3 Best bounds obtained from the quotient alone

Assume first that `J_h` is squarefree. If `r_z` quotient critical points have
the same squared critical value z, each contributes local degree two in the
fiber of the degree-`3h` map `Phi_h`. Hence

```text
r_z <= floor(3h/2).                                    (4.7)
```

There are `m=2h-2` quotient critical points in total. The number of unordered
colliding pairs is

```text
E_h=sum_z binom(r_z,2).
```

Without (4.7), the quotient gives only

```text
E_h <= binom(2h-2,2).
```

With (4.7), convexity improves this to

```text
R=min(2h-2,floor(3h/2)),
E_h <= binom(R,2)+binom(2h-2-R,2),                     (4.8)
```

where `binom(a,2)=0` for `a<=1`.

This is still quadratic in h. The off-diagonal quotient curve

```text
[Phi_h(U)-Phi_h(V)]/(U-V)=0
```

has bidegree at most `(3h-1,3h-1)`. Intersecting it with
`J_h(U)=J_h(V)=0` cannot beat the `m(m-1)` ordered critical grid by plain
Bezout. Thus reflection removes the forced sign pairing and improves constants,
but it does not supply the linear cap (4.6).

Without all-h squarefreeness of `J_h`, even (4.7) must be replaced by the
weaker ramification-defect accounting in a degree-`3h` fiber. No useful
all-h gcd bound follows. Therefore the multiplicity-cap part of Task A.2 is
OPEN beyond the corrected reduction (4.4)--(4.6).

## 5. Tasks B and C: exact certificates through h=32

### 5.1 Method

The existing verifier `CODEX_irred_verify.py` was run with height 32. For a
prime ell and each h it:

1. constructs `N_h,q_h,A_h` over `F_ell` from the exact recurrence;
2. obtains `C_h(T)` by interpolation from `4h-3` nonzero T-values;
3. checks three additional resultant values not used in the interpolation;
4. checks the expected degree `4h-4`, evenness, nonzero constant term, and
   `gcd(C_h,C_h')=1`;
5. compares `h=2,3,4` with direct resultants over ZZ;
6. checks every pairwise gcd among the resulting critical-value polynomials.

Avoiding `T=0` during interpolation is essential because the x-degree of
`N_h-Tq_h^3` drops by three there. Degree preservation at the end is the gate
which makes reduction a proof over Q, rather than sampling.

The commands were

```text
PYTHONDONTWRITEBYTECODE=1 python3 -u CODEX_irred_verify.py \
  --height 32 --primes 65537

PYTHONDONTWRITEBYTECODE=1 python3 -u CODEX_irred_verify.py \
  --height 32 --primes 104729
```

Environment: Python 3.9.6, SymPy 1.14.0.

The old prime 1009 is not good for the extended range: it passes through
`h=17`, but at `h=18` one has

```text
lambda_18 = 0 mod 1009,
```

and the expected resultant degree is not preserved. This is a modular
degeneration, not a characteristic-zero failure, and it was rejected by the
degree gate.

### 5.2 Per-h profile for the new range

Both good primes give the following profile.

| h | deg C_h | threshold 2h-1 | s_h | ord_0 C_h | deg gcd(C_h,C_h') | verdict |
|---:|---:|---:|---:|---:|---:|:---|
| 17 | 64  | 33 | 64  | 0 | 0 | PASS |
| 18 | 68  | 35 | 68  | 0 | 0 | PASS |
| 19 | 72  | 37 | 72  | 0 | 0 | PASS |
| 20 | 76  | 39 | 76  | 0 | 0 | PASS |
| 21 | 80  | 41 | 80  | 0 | 0 | PASS |
| 22 | 84  | 43 | 84  | 0 | 0 | PASS |
| 23 | 88  | 45 | 88  | 0 | 0 | PASS |
| 24 | 92  | 47 | 92  | 0 | 0 | PASS |
| 25 | 96  | 49 | 96  | 0 | 0 | PASS |
| 26 | 100 | 51 | 100 | 0 | 0 | PASS |
| 27 | 104 | 53 | 104 | 0 | 0 | PASS |
| 28 | 108 | 55 | 108 | 0 | 0 | PASS |
| 29 | 112 | 57 | 112 | 0 | 0 | PASS |
| 30 | 116 | 59 | 116 | 0 | 0 | PASS |
| 31 | 120 | 61 | 120 | 0 | 0 | PASS |
| 32 | 124 | 63 | 124 | 0 | 0 | PASS |

Together with the old range, this proves exactly:

```text
BG(h) for every 2 <= h <= 32,
s_h=4h-4 for every 2 <= h <= 32,
gcd(C_h,C_k)=1 for all 465 pairs 2 <= h < k <= 32.
```

There are 105 old pairs with `k<=16`, so the extension certifies 360 new
pairs. Each prime independently checks all 465 pairs.

Certificate digests (coefficients in ascending T-degree, canonically reduced
modulo the indicated prime) are:

```text
p=65537:  2ab63d4472cd37cfd7666ea86b906f229750e80620dda81dad99b948bb769f1a
p=104729: fc112eb324b5d0e2aeed8e064eb5d7a99f8e150a116a8067cd495769711aa444
```

### 5.3 Multiplicity and zero-column failure search

For every `2 <= h <= 32`, at each of the two good primes:

```text
deg gcd(C_h,C_h') = 0,
ord_(T=0) C_h     = 0,
all 4h-4 roots have multiplicity one.
```

Because the reductions preserve degree, these are characteristic-zero
certificates. No collision, higher critical point, or zero critical value was
found. In the quotient notation, every `D_h` has `2h-2` simple nonzero roots.

The first failure encountered while extending the old computation was instead
the bad reduction `(ell,h)=(1009,18)`, explained completely by the vanishing
of the leading coefficient `lambda_18 mod 1009`. A second good prime restores
the full degree and squarefree profile, so this is not a structural candidate
for failure of `[CRIT-2H]` over Q.

## 6. Asymptotic heuristic and exact frontier

HEURISTIC. After quotienting reflection, the constrained object is a degree
`2h-2` polynomial `D_h` of squared critical values. Reflection explains only
the passage `t -> -t`; it does not force two different quotient critical
points to have the same squared value. In a generic one-parameter-free
algebraic configuration, `Disc(D_h)` should therefore be nonzero. The full
Morse certificates through h=32 support this heuristic with no observed loss.

An actual failure of `[CRIT-2H]` would be substantial. Since `s_h` is even,
failure means `s_h <= 2h-2`. With `deg C_h=4h-4`, at least `2h-2` degrees of C
would then lie on non-simple roots. Consequently

```text
deg gcd(C_h,C_h') >= h-1.                              (6.1)
```

Equivalently, on the quotient,

```text
deg gcd(D_h,D_h') >= ceil((h-1)/2).                    (6.2)
```

Thus a failure cannot be a one-degree perturbation of the observed profile;
it requires a linear-scale discriminant degeneration, whether concentrated
in one high-multiplicity value or spread among many double values.

PROVED all h: the zero column is absent by Section 3. The remaining all-h
problem is exactly to rule out the linear-scale degeneration (6.2). Neither
the pole formula, reflection, the quotient Bezout count, nor the degree of the
fiber currently does this. Hence the final status is

```text
[C_h(0) != 0, all h]                 PROVED
[N_h squarefree, all h]              PROVED
[BG(h), 2 <= h <= 32]                VERIFIED-32
[CRIT-2H, all h]                     OPEN
[pairwise gcd(C_h,C_k), h,k <= 32]   VERIFIED-32
[pairwise gcd(C_h,C_k), all h != k]  OPEN
```
