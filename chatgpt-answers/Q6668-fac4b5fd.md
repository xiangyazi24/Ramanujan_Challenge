ANSWER Q6668 fac4b5fd

# Verdict

The redirected Montes route on the point-side quotient A_h does not prove all-h irreducibility, nor does it isolate an infinite irreducible subfamily at p=3 or p=17.

It does produce a clean negative theorem and a complete low-height calibration:

- At p=3, the first residual type already has at least two distinct branches for every h >= 7 with 3 not dividing h. Therefore no single higher-order side and no one-leaf Montes certificate can occur on either infinite class h = 1 or 2 mod 3.

- For h divisible by 3, the residual continuant has a genuine base-3 hierarchy. In particular R_(9k) is identically zero mod 3. The first-correction recurrence is then insufficient to collapse the tree to one uniform type.

- At p=17, pole collisions are absent only for h <= 17. For h > 17 the ordinary model is unsaturated. Moreover every even h with 17 not dividing h has a Q_17-linear branch coming from the top two centered coefficients of A_h. Thus p=17 also cannot yield a uniform irreducibility theorem.

- The small cases h=3,4,5 can be worked completely at p=3. Their Q_3 factor-degree multisets are respectively [2,2], [2,4], and [4,4].

- Independently of Montes, the newly landed exact certificate file proves that A_h is irreducible over Q, squarefree, nonzero at zero, and noncollapsed for every 2 <= h <= 32. This is a finite certificate theorem, not an all-h mechanism.

Thus the deliverable is alternative (ii): an exact account of where the higher-order analyses lose, together with explicit h=3,4,5 polynomials, local Montes data, finite irreducibility data, and an honest Galois-group status.

# 1. Setup and normalization

Write

```plain text
q_h(X) = product_(j=1)^h (X+j),
K_h(X) = q_h(X) N_h'(X) - 3 q_h'(X) N_h(X).
```

This K_h is the critical-point polynomial called C_h(X) in the updated campaign notation. It satisfies

```plain text
K_h(-h-1-X) = K_h(X),
deg K_h = 4h-4.
```

Put

```plain text
z = X + (h+1)/2,
U = z^2.
```

Let A_h(U) be the primitive associate with positive leading coefficient such that K_h(X) is a nonzero rational multiple of A_h(z^2). Then

```plain text
deg A_h = 2h-2.
```

For local calculations at 3 it is convenient to use the integral unit-scaled variable

```plain text
V = (2X+h+1)^2 = 4U.
```

Since 4 is a unit in Z_3, replacing U by V changes neither Q_3 factor degrees nor the Montes tree.

# 2. The exact mod-3 continuant and its three-step monodromy

Modulo 3,

```plain text
P(u) = (u-1)^3.
```

Define R_h in F_3[X] by

```plain text
R_1 = 1,
R_2 = X,
R_(h+1) = (X+h-1) R_h - (X+h)^2 R_(h-1).
```

Then the Frobenius identity gives

```plain text
N_h = R_h^3 mod 3.
```

Put

```plain text
q = X^3-X.
```

Multiplying the three periodic transfer matrices gives a matrix B with

```plain text
trace(B) = q,
det(B) = q^2,
B = -q I + J,
J^2 = 0.
```

Consequently, for m >= 0, with the terms containing m(-q)^(m-1) interpreted as zero when m=0,

```plain text
R_(3m+1) = (-q)^m + m(-q)^(m-1)(X+1),
R_(3m+2) = X(-q)^m + m(-q)^(m-1)(2X^2+1),
R_(3m+3) = (2-m)(-q)^m.
```

In particular,

```plain text
R_(3m) = -m(-q)^(m-1),
R_(9k) = 0 in F_3[X].
```

This is the first decisive structural fact. The residual type is controlled by the base-3 digits of h; it is not governed by one periodic nonzero polynomial.

In the centered coordinate T=2X+h+1 one has

```plain text
X^3-X = T-T^3 mod 3.
```

Thus the same nilpotent three-step structure survives the reflection descent to A_h.

# 3. First correction and the precise residual object

Choose integer lifts and write

```plain text
N_h = R_h^3 + 3 S_h mod 9.
```

Set

```plain text
E(u) = (P(u)-(u-1)^3)/3
     = 11u^3+18u^2+8u+2.
```

For a=X+h, the first-correction recurrence from Q6579 is

```plain text
S_(h+1)
 = E(a)R_h^3 + (a-1)^3 S_h - a^6 S_(h-1)
   + (a-1)a^2 R_h R_(h-1) R_(h+1) mod 3,
S_1 = 0,
S_2 = E(X+1) mod 3.
```

When the exact 3-content of K_h is one,

```plain text
K_h/3
 = q_h(R_h^2 R_h' + S_h') - q_h' R_h^3 mod 3.
```

This formula remains the correct first Montes residual. The obstruction is not that the formula is wrong. The obstruction is that its factorization is forced to branch for two infinite congruence classes of h, and it degenerates to deeper levels on a third class.

# 4. An infinite negative theorem at p=3

Let lambda_h be the leading coefficient of N_h. It satisfies

```plain text
lambda_(h+1) = 34 lambda_h - lambda_(h-1),
lambda_1 = lambda_2 = 1 mod 3.
```

Hence

```plain text
lambda_h mod 3 = 1,1,0,2,2,0,...,
```

with period six. Therefore lambda_h is a 3-adic unit exactly when 3 does not divide h.

The leading coefficient of K_h is

```plain text
lc(K_h) = -3 lambda_h.
```

Thus, when 3 does not divide h, the exact coefficient-content valuation is

```plain text
v_3(content(K_h)) = 1.
```

Now use the exact pole evaluation

```plain text
K_h(-j)/3
 = -((j-1)!(h-j)!)^4 b_(j-1)b_(h-j),
1 <= j <= h.
```

If h >= 7, then for every j at least one of j-1 and h-j is at least 3. Therefore

```plain text
K_h(-j)/3 = 0 mod 3
```

at every pole. The pole residues cover all three elements of F_3. In the centered coordinate T=2X+h+1, the first residual polynomial is even and vanishes at

```plain text
T = 0, 1, -1.
```

Evenness makes the root at T=0 have multiplicity at least two. Hence its V=T^2 descent is divisible by

```plain text
V(V-1).
```

Therefore:

```plain text
THEOREM 3-RESIDUAL-SPLIT.
For every h >= 7 with 3 not dividing h,
the primitive reduction A_h(V) mod 3 is divisible by V(V-1).
In particular the first Montes level has at least two distinct branches,
and A_h is reducible over Q_3.
```

This excludes a one-leaf p=3 proof on both infinite classes h=1 mod 3 and h=2 mod 3.

For h divisible by 3, lambda_h is no longer a unit and the degree drops at the first residual level. At h divisible by 9 the stronger identity R_h=0 mod 3 forces a further jump in type order. The P6 recurrence then depends on S_h and does not reduce to a scalar periodic sequence. This is exactly where the remaining h=0 mod 3 branch stalls: it is a growing MacLane tree indexed by deeper 3-adic data, not one coprime-slope side.

# 5. Complete p=3 calibration at h=3,4,5

The repository normalization U=z^2 gives the following primitive polynomials.

```plain text
A_3(U)
 = 1155 U^4 - 120 U^3 + 256 U^2 - 148 U + 25.
```

```plain text
A_4(U)
 = 40177664 U^6 - 91658240 U^5 + 105719552 U^4
   - 60570880 U^3 - 19055888 U^2 + 29693160 U
   + 361017.
```

```plain text
A_5(U)
 = 1332869 U^8 - 9608320 U^7 + 29298796 U^6
   - 45582524 U^5 + 31760187 U^4 - 11908448 U^3
   + 22394112 U^2 - 9687296 U + 1364224.
```

Each is irreducible over Q by the exact characteristic-zero factorization now stored in CODEX_AH_CERT_report.md.

## 5.1 Height h=3

Use V=4U and multiply by a 3-adic unit. The polynomial becomes

```plain text
F_3(V)
 = 1155V^4 - 480V^3 + 4096V^2 - 9472V + 6400.
```

Its coefficient valuations in increasing degree are

```plain text
0, 0, 0, 1, 1.
```

The horizontal residual is

```plain text
(V+1)^2.
```

Take the MacLane representative phi=V+1. The first three phi-adic Taylor coefficients have valuations

```plain text
v_3(c_0), v_3(c_1), v_3(c_2) = 1, 2, 0.
```

The lower side has length two and slope -1/2. Its residual polynomial is linear, so this repeated branch is one irreducible ramified quadratic. The other ordinary side also has length two with linear residual.

Therefore

```plain text
A_3 over Q_3 has factor degrees [2,2].
```

## 5.2 Height h=4

Modulo 3,

```plain text
A_4(V) = V^2 Q_4(V),
Q_4(V) = 2V^4 + V^3 + 2V^2 + 2V + 1.
```

The quartic Q_4 is irreducible over F_3. The nonhorizontal side has coefficient valuations

```plain text
4, 2, 0
```

and residual polynomial

```plain text
Y^2 + 2Y + 2,
```

which is irreducible over F_3. Hence

```plain text
A_4 over Q_3 has factor degrees [2,4].
```

## 5.3 Height h=5

Modulo 3,

```plain text
A_5(V) = (V-1)^4 Q_4(V),
```

with the same irreducible quartic Q_4. For phi=V-1, the relevant translated coefficient valuations are

```plain text
4, 4, 2, 3, 0.
```

The lower side has length four and slope -1. Its residual polynomial is

```plain text
2Z^4 + Z^2 + 1.
```

After making it monic this is

```plain text
Z^4 + 2Z^2 + 2.
```

It is irreducible over F_3: the quadratic in Z^2 is irreducible over F_3, and its roots have nonsquare norm in F_9. Thus the repeated phi-branch lifts as one irreducible quartic. Consequently

```plain text
A_5 over Q_3 has factor degrees [4,4].
```

These three examples show the actual phenomenon. Higher-order analysis resolves inseparable first residual factors, but it resolves them into several local factors rather than one all-degree factor.

# 6. The p=17 route

The prime 17 has two separate obstructions.

## 6.1 Pole collisions

For h <= 17, the poles -1,...,-h remain distinct mod 17. For h > 17, some j and j+17 are both poles. The endpoint law gives

```plain text
N_h(-j)
 = sign * ((j-1)!(h-j)!)^3 b_(j-1)b_(h-j).
```

Whenever j+17 <= h, the factor (h-j)! is divisible by 17. Hence q_h and N_h already share roots after reduction. The ordinary critical-point model is unsaturated, and the pole-cluster factors must be removed before any Newton polygon is interpreted.

Therefore there is no infinite collision-free h-range at the fixed prime 17. The collision-free range is exactly the finite range h <= 17.

## 6.2 A forced linear local branch for even h

The leading coefficients satisfy the proved valuation law

```plain text
v_17(lambda_h) = 0                         for h odd,
v_17(lambda_h) = 1+v_17(h/2)              for h even.
```

Center N_h and q_h as

```plain text
N_h = lambda_h z^(3h-3) + u_h z^(3h-5) + ...,
q_h = z^h + q_2 z^(h-2) + ....
```

The top two coefficients of K_h are

```plain text
[z^(4h-4)] K_h = -3 lambda_h,
[z^(4h-6)] K_h = -5u_h + 3q_2 lambda_h.
```

The known closed formula is

```plain text
u_h
 = (-(h-1)(32h^2-64h-11)lambda_h
    +5h lambda_(h-1))/256.
```

If h is even and 17 does not divide h, then lambda_h is divisible by 17 while lambda_(h-1) is a unit. Hence u_h is a unit, and the second displayed coefficient is a unit. In the descended polynomial A_h, the leading coefficient is divisible by 17 while the next coefficient is a unit.

Thus the 17-adic Newton polygon has a final side of horizontal length one. Its residual polynomial is linear, so A_h has a Q_17-linear factor.

Therefore:

```plain text
THEOREM 17-LINEAR-BRANCH.
For every even h with 17 not dividing h,
A_h has a linear factor over Q_17.
```

This is another infinite negative theorem for local irreducibility. The prime 17 can still supply useful small cycles after global transitivity is known, but it cannot prove transitivity uniformly.

# 7. Factorization and Galois data for h=3,...,8

The exact factorization over Q is now stronger than the original requested calibration:

- A_3 is irreducible of degree 4.

- A_4 is irreducible of degree 6.

- A_5 is irreducible of degree 8.

- A_6 is irreducible of degree 10.

- A_7 is irreducible of degree 12.

- A_8 is irreducible of degree 14.

The same exact certificate file proves irreducibility for every height through h=32.

The current repository output does not compute the full global Galois groups of A_3,...,A_8. Irreducibility proves only transitivity. It would be incorrect to promote the expected groups S_(2h-2) to theorems without additional Frobenius-cycle and discriminant certificates.

One additional exact global datum is available at h=4. Modulo 5,

```plain text
A_4(U)
 = (U-1)(U+1)(U^4+4U^2+2).
```

The quartic factor is irreducible over F_5. Hence Gal(A_4/Q) contains a permutation of cycle type [1,1,4], which is odd. Therefore the discriminant of A_4 is not a square and the Galois group is not contained in A_6. This still does not determine the exact transitive subgroup of S_6.

For h=3,5,6,7,8, the exact group and discriminant square class are not recorded in CODEX_AH_CERT_report.md. Closing that finite subtask requires a standard cycle-certificate run: obtain an irreducible good-prime reduction, a long prime cycle excluding blocks, and an odd cycle type or exact discriminant. That computation is finite and cheap, but it is separate from the all-h Montes mechanism and was not performed by the current certificate script.

# 8. The Capell non-square condition for h=3,4,5

There is a different non-square condition relevant to lifting from A_h(U) to the full even critical polynomial A_h(z^2).

Let alpha be a root of irreducible A_h. If alpha were a square in Q(alpha), then its norm to Q would be a rational square. Since deg A_h is even,

```plain text
Norm(alpha) = A_h(0)/lc(A_h).
```

Thus a nonsquare value of this ratio is a sufficient Capell certificate for irreducibility of A_h(z^2).

For h=3,

```plain text
A_3(0)/lc(A_3) = 25/1155 = 5/231,
```

which is not a square in Q.

For h=4,

```plain text
A_4(0)/lc(A_4) = 361017/40177664.
```

The numerator and denominator are coprime. Moreover

```plain text
361017 = 3^4 * 4457,
40177664 = 2^12 * 9809,
66^2 < 4457 < 67^2,
99^2 < 9809 < 100^2.
```

Hence the reduced numerator and denominator are not both squares, so the ratio is not a rational square.

For h=5,

```plain text
A_5(0)/lc(A_5) = 1364224/1332869
               = 1168^2/1332869.
```

The numerator and denominator are coprime, and

```plain text
1154^2 < 1332869 < 1155^2.
```

Thus this ratio is also not a rational square.

Consequently the full even critical-point polynomials A_h(z^2), equivalently K_h(X) up to a scalar and a linear change of variable, are irreducible over Q for h=3,4,5.

This Capell condition is not the same as the discriminant-square test used to distinguish S_n from A_n. The former is proved here for h=3,4,5; the latter remains part of the unrun finite Galois certification except for the odd Frobenius datum at h=4.

# 9. Final route verdict

The p=3 Montes route is not merely unfinished. It is structurally incapable of proving A_h irreducible for all h by one local type:

```plain text
h >= 7, 3 not dividing h:
first residual splits by V(V-1);

h divisible by 3:
degree drop and deeper base-3 type tree;

h divisible by 9:
R_h is identically zero mod 3.
```

The p=17 route is likewise structurally nonuniform:

```plain text
h > 17:
pole collisions and unsaturated reduction;

even h, 17 not dividing h:
forced Q_17-linear branch.
```

Therefore no infinite-family irreducibility theorem emerges from these two structural primes.

The positive result is finite but substantial:

```plain text
A_h irreducible over Q,
A_h squarefree,
A_h(0) != 0,
NONCOLLAPSE_h,
```

for every 2 <= h <= 32, with exact reproducible certificates.

For an all-h theorem, the next viable route must use primes depending on h, or a global recurrence/Galois argument. Fixed-prime Montes data can still provide cycle types and block exclusions after transitivity is obtained, but cannot itself supply uniform transitivity.

# Repository sources used

- problems/3.2/CODEX_AH_CERT_report.md, commit 346b43992bb6a52a8d67d57694c92416b6a8e63e.

- problems/3.2/CODEX_ah_cert.py, same commit.

- drops/Q6579-38379a8c.md on chatgpt-drop, containing the P5/P6/P7 derivations.

- problems/3.2/CODEX_MONTES_report.md, the older M_h local-factor pipeline, used only as a comparison and not as evidence about A_h.