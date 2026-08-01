ANSWER Q6708 d43acede

# Uniform C2 conjugation with polynomial drift: complete proof

## 0. Verdict and one normalization correction

The desired statement is true, with the following precise normalization.

If, in an eigenbasis of the limiting matrix,

```plain text
B_n(z) = Lambda(z) + C(z)/n + O_C2(n^-2),
Lambda = diag(lambda_+,lambda_-),
```

then the polynomial exponents are not the raw diagonal entries C_ii. They are

```plain text
d_+(z) = C_++(z)/lambda_+(z),
d_-(z) = C_--(z)/lambda_-(z).
```

Indeed,

```plain text
lambda_i + C_ii/n = lambda_i (1 + d_i/n).
```

The proof has four exact stages.

1. Diagonalize the limiting matrix analytically on a parameter neighborhood.

1. Remove the off-diagonal 1/n term by the explicit homological conjugacy I+H/n. The remainder is then genuinely O_C2(n^-2) and summable.

1. Divide out the diagonal products

```plain text
prod_{k=n0}^{n-1} lambda_i (1+d_i/k).
```

The remaining column equations are solved by a mixed forward/backward Green operator. That operator is a contraction in an explicitly weighted C^2 Banach norm.

1. Compare the finite products with lambda_i^n n^d_i by an absolutely convergent logarithmic product.

The construction actually gives the exact cocycle identity

```plain text
P_{n+1}(z)^(-1) A_n(z) P_n(z)
 = diag(lambda_+(z)(1+d_+(z)/n),
        lambda_-(z)(1+d_-(z)/n))
```

for every n>=n0, not merely an asymptotic diagonalization.

---

# 1. Parameter domain and the weighted C2 algebra

Let Omega be an open subset of C, let K be compact, and choose a second compact set K1 such that

```plain text
K subset int(K1) subset K1 subset Omega.
```

All estimates below are made on K1; the final conclusion is restricted to K. This small enlargement is what makes analyticity automatic: all fixed-point iterates are holomorphic on int(K1) and converge uniformly there.

For 0<eps<=1, define

```plain text
||f||_(2,eps)
 := sup_K1 |f|
  + eps sup_K1 |f'|
  + (eps^2/2) sup_K1 |f''|.
```

For vectors use the maximum of the component norms. For matrices use the induced row-sum norm.

This is a Banach algebra norm. Leibniz gives

```plain text
||fg||_(2,eps) <= ||f||_(2,eps) ||g||_(2,eps),
||MN||_(2,eps) <= ||M||_(2,eps) ||N||_(2,eps).
```

It is equivalent to the ordinary C^2 norm. In particular,

```plain text
max_{0<=j<=2} sup_K |f^(j)|
 <= 2 eps^(-2) ||f||_(2,eps).                 (1.1)
```

The derivative weights are essential. Multiplication by

```plain text
rho(z)=lambda_-(z)/lambda_+(z)
```

need not be a contraction in the unweighted C^2 norm even when sup |rho|<1, because derivatives of rho contribute. The equivalent weighted norm makes the contraction quantitative.

Put

```plain text
q0  = sup_K1 |rho| < 1,
rho1 = sup_K1 |rho'|,
rho2 = sup_K1 |rho''|,
```

and choose

```plain text
eps0 = min(
  1,
  (1-q0)/(4 max(1,rho1)),
  sqrt((1-q0)/(2 max(1,rho2)))
).
```

Use eps=eps0 from now on. Then

```plain text
qbar := ||rho||_(2,eps)
 <= q0 + eps rho1 + (eps^2/2) rho2
 <= (1+q0)/2 < 1.                            (1.2)
```

All norms below mean ||.||_(2,eps) on K1.

---

# 2. Precise hypotheses

Assume that for every n>=Nexp, A_n(z) is holomorphic near K1 and

```plain text
A_n = A_0 + A_1/n + R_n,
||R_n|| <= M_2/n^2.                           (2.1)
```

Assume that on a neighborhood of K1 there are holomorphic nonzero eigenvalues and a holomorphic invertible eigenvector matrix

```plain text
S^(-1) A_0 S = Lambda
              = diag(lambda_+,lambda_-),
```

with

```plain text
sup_K1 |lambda_-/lambda_+| = q0 < 1.          (2.2)
```

Define

```plain text
S0 = ||S||,
T0 = ||S^(-1)||,
L0 = max(||lambda_+||,||lambda_-||),
I0 = max(||lambda_+^(-1)||,||lambda_-^(-1)||),
G0 = ||(lambda_+-lambda_-)^(-1)||.
```

In the eigenbasis put

```plain text
B_n := S^(-1) A_n S
     = Lambda + C/n + F_n,
C   := S^(-1) A_1 S.
```

Then one may take

```plain text
C1 = ||C|| <= T0 ||A_1|| S0,
E2 = sup_{n>=Nexp} n^2 ||F_n|| <= T0 M_2 S0.  (2.3)
```

Let

```plain text
Cdiag = diag(C_++,C_--),
d_+ = C_++/lambda_+,
d_- = C_--/lambda_-,
d0  = max(||d_+||,||d_-||) <= I0 C1.          (2.4)
```

These d_+,- are the required polynomial correction exponents.

---

# 3. Remove the off-diagonal 1/n term

Define the off-diagonal matrix H by

```plain text
H_++ = H_-- = 0,
H_+- = -C_+-/(lambda_+-lambda_-),
H_-+ = -C_-+/(lambda_--lambda_+).
```

Equivalently,

```plain text
C + Lambda H - H Lambda = Cdiag.              (3.1)
```

Set

```plain text
H0 := ||H|| <= G0 C1,
T_n := I + H/n.
```

For n>=2H0, T_n is invertible and the Neumann series gives

```plain text
||T_n^(-1)|| <= 2.                            (3.2)
```

Define the diagonal drift matrix

```plain text
D_n := Lambda + Cdiag/n
     = diag(a_+,n,a_-,n),

a_i,n := lambda_i(1+d_i/n).                  (3.3)
```

Make the change X_n=T_n Z_n. Then

```plain text
Z_{n+1} = (D_n+E_n) Z_n,                      (3.4)
```

where

```plain text
E_n = T_{n+1}^(-1)
      [B_n T_n - T_{n+1}D_n].
```

The bracket can be computed exactly. First,

```plain text
B_n T_n
 = Lambda
 + (Lambda H+C)/n
 + C H/n^2
 + F_n
 + F_n H/n,
```

while

```plain text
T_{n+1}D_n
 = Lambda + Cdiag/n
 + H Lambda/(n+1)
 + H Cdiag/[n(n+1)].
```

Using (3.1),

```plain text
B_nT_n-T_{n+1}D_n
 = F_n
 + C H/n^2
 + F_nH/n
 + H Lambda/[n(n+1)]
 - H Cdiag/[n(n+1)].                          (3.5)
```

Consequently

```plain text
||B_nT_n-T_{n+1}D_n|| <= B0/n^2,
```

with the explicit constant

```plain text
B0 := E2(1+H0) + 2 C1 H0 + L0 H0.            (3.6)
```

By (3.2),

```plain text
||E_n|| <= L/n^2,
L := 2B0.                                     (3.7)
```

This is the summable perturbation promised in the question:

```plain text
sum_{n>=n0} ||E_n|| <= L/(n0-1).
```

No 1/n off-diagonal term remains.

---

# 4. Bounds for the diagonal ratios

For n>=2d0, each a_i,n is invertible and

```plain text
||a_i,n^(-1)||
 <= 2 ||lambda_i^(-1)||
 <= Ainv,
Ainv := 2I0.                                  (4.1)
```

Put

```plain text
r_n := a_-,n/a_+,n
     = rho (1+d_-/n)(1+d_+/n)^(-1).
```

Because ||(1+d_+/n)^(-1)||<=2,

```plain text
||r_n-rho||
 <= ||rho|| * ||d_--d_+||/n * 2
 <= Cr/n,
Cr := 4 qbar d0.                              (4.2)
```

Hence, if

```plain text
n >= 2Cr/(1-qbar),
```

then

```plain text
||r_n|| <= q1,
q1 := (1+qbar)/2 < 1.                         (4.3)
```

Every product of consecutive r_n therefore satisfies

```plain text
||prod_{ell=u}^v r_ell|| <= q1^(v-u+1).       (4.4)
```

This estimate includes two parameter derivatives because it is an estimate in the weighted C^2 Banach algebra.

---

# 5. The Green operators after division by the drift products

We now construct an exact matrix Q_n=I+U_n such that

```plain text
(D_n+E_n)Q_n = Q_{n+1}D_n.                   (5.1)
```

For a column j in {+,-}, write

```plain text
Q_n e_j = e_j+u_j,n.
```

The jth column of (5.1) is equivalent to

```plain text
u_j,n+1 - (D_n/a_j,n)u_j,n
 = a_j,n^(-1) E_n(e_j+u_j,n).                (5.2)
```

The left side has one neutral component and one hyperbolic component. We invert it with a mixed forward/backward Green operator.

## 5.1 The plus column

For a forcing sequence g_n=(g_n^+,g_n^-), define

```plain text
(G_+ g)_n^+
 := -sum_{k=n}^infinity g_k^+,

(G_+ g)_n^-
 := sum_{k=n0}^{n-1}
      (prod_{ell=k+1}^{n-1} r_ell) g_k^-,
```

where the second sum is zero at n=n0 and an empty product equals one.

Then directly

```plain text
x_n=G_+g
```

satisfies

```plain text
x_{n+1}^+ = x_n^+ + g_n^+,
x_{n+1}^- = r_n x_n^- + g_n^-,
x_n^- at n=n0 is 0.                          (5.3)
```

## 5.2 The minus column

The expanding ratio in the plus component is r_n^(-1). Define

```plain text
(G_- g)_n^+
 := -sum_{k=n}^infinity
      (prod_{ell=n}^k r_ell) g_k^+,

(G_- g)_n^-
 := -sum_{k=n}^infinity g_k^-.
```

A one-line index shift gives

```plain text
x_{n+1}^+ = r_n^(-1)x_n^+ + g_n^+,
x_{n+1}^- = x_n^- + g_n^-.                   (5.4)
```

Thus G_- imposes the unique boundary condition for which the expanding component tends to zero at infinity.

## 5.3 Explicit Green bounds

For a sequence x=(x_n)_{n>=n0}, use

```plain text
||x||_[s] := sup_{n>=n0} n^s ||x_n||.
```

Define

```plain text
CG := 2 + 2/(1-q1) + 1/(1-q1)^2.             (5.5)
```

Then for both signs,

```plain text
||G_+- g||_[1] <= CG ||g||_[2],              (5.6)
||G_+- g||_[1] <= (CG/n0) ||g||_[3].         (5.7)
```

Here is the complete estimate.

For a neutral tail,

```plain text
n sum_{k=n}^infinity 1/k^2 <= 2,
n sum_{k=n}^infinity 1/k^3 <= 2/n0.          (5.8)
```

For the stable forward sum, put ell=n-1-k. Since k>=n0,

```plain text
n/k^2 = (k+1+ell)/k^2
       <= 1/n0 + (1+ell)/n0^2.
```

Therefore

```plain text
n sum_{k=n0}^{n-1} q1^(n-1-k)/k^2
 <= 1/[n0(1-q1)] + 1/[n0^2(1-q1)^2].        (5.9)
```

Similarly,

```plain text
n/k^3
 <= 1/n0^2 + (1+ell)/n0^3,
```

which gives

```plain text
n sum_{k=n0}^{n-1} q1^(n-1-k)/k^3
 <= 1/[n0^2(1-q1)]
  + 1/[n0^3(1-q1)^2].                        (5.10)
```

For the backward hyperbolic tail,

```plain text
n sum_{k=n}^infinity q1^(k-n+1)/k^2
 <= q1/[n(1-q1)],                            (5.11)

n sum_{k=n}^infinity q1^(k-n+1)/k^3
 <= q1/[n^2(1-q1)].                          (5.12)
```

Equations (5.8)-(5.12) imply (5.6)-(5.7) with the stated CG.

---

# 6. Banach fixed point in C2(K)

For j in {+,-}, define

```plain text
(N_j u)_n
 := a_j,n^(-1) E_n(e_j+u_n),

F_j(u) := G_j(N_j u).                         (6.1)
```

Work in the Banach space

```plain text
X := {u : ||u||_[1] < infinity}
```

with the C^2(K1) weighted norm above.

From (3.7) and (4.1),

```plain text
||(N_j u)_n||
 <= Ainv L/n^2 * (1+||u_n||).                (6.2)
```

Let

```plain text
Rc := 2 CG Ainv L.                            (6.3)
```

Consider the closed ball

```plain text
B_Rc := {u in X : ||u||_[1] <= Rc}.
```

If n0>=Rc, then ||u_n||<=Rc/n0<=1, so (5.6) and (6.2) give

```plain text
||F_j(u)||_[1]
 <= CG Ainv L(1+Rc/n0)
 <= 2 CG Ainv L
 = Rc.                                        (6.4)
```

Thus F_j maps the ball into itself.

For two elements u,v,

```plain text
(N_j u-N_j v)_n
 = a_j,n^(-1)E_n(u_n-v_n),
```

and hence

```plain text
||N_j u-N_j v||_[3]
 <= Ainv L ||u-v||_[1].                      (6.5)
```

By (5.7),

```plain text
||F_j(u)-F_j(v)||_[1]
 <= theta ||u-v||_[1],

theta := CG Ainv L/n0.                       (6.6)
```

Therefore F_j is a contraction whenever n0>CG Ainv L.

For the single explicit threshold used below, impose

```plain text
n0 >= 4Rc = 8 CG Ainv L.                     (6.7)
```

Then

```plain text
theta <= 1/8,
Rc/n0 <= 1/4.                                (6.8)
```

Banach's theorem gives unique fixed points u_+,u_- in B_Rc. By construction, the matrix

```plain text
U_n := [u_+,n  u_-,n],
Q_n := I+U_n
```

satisfies (5.1). Moreover, because each column has norm at most Rc/n,

```plain text
||U_n|| <= 2Rc/n.                            (6.9)
```

Under (6.7), ||U_n||<=1/2, so every Q_n is invertible and

```plain text
||Q_n^(-1)|| <= 2.                           (6.10)
```

Every fixed-point iterate starting from zero is holomorphic in z; convergence occurs in the C^2(K1) norm and hence uniformly on K1. The fixed points, Q_n, and all matrices constructed from them are holomorphic on int(K1). This proves the requested analytic parameter dependence.

---

# 7. The exact conjugating sequence and its C2 rate

Return to the original basis and define

```plain text
P_n := S T_n Q_n
     = S(I+H/n)(I+U_n),
P_infinity := S.                             (7.1)
```

Combining the definitions of T_n and Q_n gives the exact identity

```plain text
A_n P_n = P_{n+1}D_n,                        (7.2)
```

or

```plain text
P_{n+1}^(-1) A_n P_n
 = diag(lambda_+(1+d_+/n),
        lambda_-(1+d_-/n)).                  (7.3)
```

This is the desired analytic diagonalization.

From ||H||<=H0 and (6.9),

```plain text
||(I+H/n)(I+U_n)-I||
 <= H0/n + 2Rc/n + 2H0Rc/n^2.
```

Therefore

```plain text
||P_n-P_infinity|| <= CP/n,                  (7.4)
```

with

```plain text
CP := S0 (H0 + 2Rc + 2H0Rc/n0).              (7.5)
```

Using (1.1), the ordinary two-derivative estimate is

```plain text
max_{0<=j<=2} sup_K
 ||partial_z^j(P_n-P_infinity)||
 <= 2 eps^(-2) CP/n.                         (7.6)
```

This proves items (i) and (ii) for the conjugating sequence.

---

# 8. Explicit threshold n0(K)

All preceding requirements are simultaneously satisfied by any integer

```plain text
n0 >= Nstar,                                 (8.1)
```

where

```plain text
Bd := d0^2 + d0/2,
Cr := 4 qbar d0,
Rc := 2 CG Ainv L,

Nstar := max of
  Nexp,
  2,
  ceil(2H0),
  ceil(2d0),
  ceil(2Cr/(1-qbar)),
  ceil(4Rc),
  ceil(4Bd).                                  (8.2)
```

Every symbol on the right was defined explicitly in (1.2), (2.3)-(2.4), (3.6)-(3.7), (4.1)-(4.3), and (5.5).

The role of the individual terms is:

```plain text
2H0                 invert I+H/n,
2d0                 invert 1+d_i/n,
2Cr/(1-qbar)        make ||r_n||<=q1<1,
4Rc                  contraction, ball invariance, and invert Q_n,
4Bd                  control the logarithmic product tail.
```

---

# 9. Polynomial factors and the uniform asymptotic formula

For i in {+,-}, put

```plain text
p_i(n,z)
 := prod_{k=n0}^{n-1} a_i,k(z)
  = lambda_i(z)^(n-n0)
    prod_{k=n0}^{n-1}(1+d_i(z)/k).            (9.1)
```

The exact solutions associated with (7.2) are

```plain text
Y_i,n^raw(z) := p_i(n,z) P_n(z)e_i.           (9.2)
```

We now compare the finite product with n^d_i.

For k>=n0, define in the Banach algebra

```plain text
theta_i,k
 := Log(1+d_i/k) - d_i log(1+1/k),           (9.3)
```

where Log is the power-series branch. It is well-defined because

```plain text
||d_i/k|| <= 1/2.
```

For ||x||<=1/2,

```plain text
||Log(1+x)-x||
 <= ||x||^2/[2(1-||x||)]
 <= ||x||^2.
```

Also,

```plain text
0 <= 1/k-log(1+1/k) <= 1/(2k^2).
```

Consequently

```plain text
||theta_i,k|| <= Bd/k^2,
Bd=d0^2+d0/2.                                 (9.4)
```

The series

```plain text
Theta_i := sum_{k=n0}^infinity theta_i,k      (9.5)
```

therefore converges absolutely in C^2(K1) and is holomorphic. Define the nonzero holomorphic scalar

```plain text
c_i(z)
 := lambda_i(z)^(-n0)
    n0^(-d_i(z))
    exp(Theta_i(z)).                          (9.6)
```

Here n0^(-d_i)=exp(-d_i log n0).

Since

```plain text
prod_{k=n0}^{n-1}(1+1/k)^d_i
 = (n/n0)^d_i,
```

we obtain the exact identity

```plain text
p_i(n)
 = lambda_i^n n^d_i c_i
   exp(-tau_i,n),                             (9.7)

where

tau_i,n := sum_{k=n}^infinity theta_i,k.
```

By (9.4),

```plain text
||tau_i,n|| <= Bd/(n-1) <= 2Bd/n.            (9.8)
```

The threshold n0>=4Bd makes ||tau_i,n||<=1/2, and therefore

```plain text
||exp(-tau_i,n)-1||
 <= e ||tau_i,n||
 <= 2e Bd/n.                                  (9.9)
```

Normalize the solutions by

```plain text
Y_i,n := c_i^(-1) Y_i,n^raw.                 (9.10)
```

Let

```plain text
v_i(z):=S(z)e_i.
```

Equations (7.4), (9.7), and (9.9) give

```plain text
lambda_i(z)^(-n) n^(-d_i(z)) Y_i,n(z)
 = exp(-tau_i,n(z)) P_n(z)e_i.               (9.11)
```

Thus

```plain text
||lambda_i^(-n)n^(-d_i)Y_i,n-v_i||
 <= Casym/n,                                  (9.12)
```

where one may take

```plain text
Casym := e (CP + 2Bd S0).                     (9.13)
```

Returning to the ordinary C^2 norm by (1.1),

```plain text
max_{0<=j<=2} sup_{z in K}
 |partial_z^j[
   lambda_i(z)^(-n)n^(-d_i(z))Y_i,n(z)-v_i(z)
 ]|
 <= 2 eps^(-2) Casym/n.                       (9.14)
```

In particular, uniformly for z in K,

```plain text
Y_i,n(z)
 = lambda_i(z)^n n^d_i(z)
   (v_i(z)+O_C2,K(1/n)).                      (9.15)
```

This is the requested Benzaid-Lutz/Levinson form with the polynomial factors included.

Notice that (9.14) controls two derivatives of the normalized error. Derivatives of the explicit factor n^d_i(z) naturally contain powers of log n; those are part of the stated main term, not part of the error.

---

# 10. Uniform dichotomy

The two columns are linearly independent because P_n is invertible. Every solution on n>=n0 has the exact representation

```plain text
Y_n
 = alpha p_+(n)P_n e_+
 + beta  p_-(n)P_n e_-.                       (10.1)
```

Moreover,

```plain text
p_-(n)/p_+(n)
 = p_-(m)/p_+(m)
   prod_{k=m}^{n-1} r_k.
```

By (4.4), for n>=m>=n0,

```plain text
|| [p_-(n)/p_+(n)]/[p_-(m)/p_+(m)] ||
 <= q1^(n-m).                                 (10.2)
```

Thus the splitting is a uniform exponential dichotomy in the parameter-dependent C^2 algebra. The explicit powers n^d_+,- do not weaken hyperbolicity; they are already incorporated into the exact diagonal products.

---

# 11. Dependence on the degeneration distance

Let

```plain text
Sigma = {z in Omega : lambda_+(z)=lambda_-(z)}
```

be the eigenvalue-collision locus, and suppose

```plain text
tau = dist(K1,Sigma) > 0.
```

The proof uses the degeneration distance only through the analytic eigenbasis constants and the inverse spectral gap

```plain text
G0=||(lambda_+-lambda_-)^(-1)||_(2,eps).
```

For any scalar holomorphic g with

```plain text
delta = inf_K1 |g| > 0,
M1g = sup_K1 |g'|,
M2g = sup_K1 |g''|,
```

one has explicitly

```plain text
||g^(-1)||_(2,eps)
 <= delta^(-1)
  + eps M1g delta^(-2)
  + (eps^2/2)
      [M2g delta^(-2)+2 M1g^2 delta^(-3)].   (11.1)
```

Apply (11.1) to

```plain text
g=lambda_+-lambda_-
```

and separately to lambda_+ and lambda_-. If the coefficients are holomorphic on a complex r-neighborhood of K1, Cauchy's inequalities give

```plain text
sup_K1 |g'| <= r^(-1) sup_{K1^r}|g|,
sup_K1 |g''| <= 2r^(-2) sup_{K1^r}|g|.       (11.2)
```

One may take r to be any fixed fraction of the distance to the degeneration and coefficient-pole loci, for example r=tau/4 when no other singularity is nearer. Hence the constants in (11.1), and therefore H0,L,n0,CP,Casym, are explicit rational expressions in

```plain text
q0,
1-q0,
delta_gap^(-1),
tau^(-1),
S0,T0,
||A_1||_C2,
M_2.
```

More precisely, the worst inverse-gap terms are at most cubic in delta_gap^(-1) and quadratic in tau^(-1) because of the second derivative in (11.1).

A distance tau alone cannot bound delta_gap without a normalization of the polynomial or analytic function defining the discriminant: multiplying that function by a tiny constant leaves its zero set unchanged. For a concrete polynomial family, both tau and delta_gap=inf_K1|lambda_+-lambda_-| are directly machine-certifiable, and formulas (11.1)-(11.2) give the required explicit constants.

A global analytic eigenbasis can also have a topological obstruction on a non-simply-connected parameter set. The theorem therefore assumes the displayed analytic S. On one archimedean cell, K1 is taken contractible, so the eigenline bundles are trivial. In the Apery application A_infinity is constant, and S can be chosen constant, eliminating this issue entirely.

---

# 12. Apery specialization

For the Apery hyperbolic limit,

```plain text
lambda_+ = (1+sqrt(2))^4 = 17+12sqrt(2),
lambda_- = (1-sqrt(2))^4 = 17-12sqrt(2),

lambda_+-lambda_- = 24sqrt(2),
lambda_+ lambda_- = 1,

q0 = |lambda_-/lambda_+|
   = (17-12sqrt(2))^2
   approximately 8.6655e-4.                  (12.1)
```

The eigenvalues and eigenvectors are constant in z, so

```plain text
rho1=rho2=0,
eps=1,
qbar=q0.
```

If

```plain text
C(z)=S^(-1)A_1(z)S,
```

then the exact exponents are

```plain text
d_+(z)=C_++(z)/(17+12sqrt(2)),
d_-(z)=C_--(z)/(17-12sqrt(2)).                (12.2)
```

The first conjugacy is especially explicit:

```plain text
H_+-(z) = -C_+-(z)/(24sqrt(2)),
H_-+(z) =  C_-+(z)/(24sqrt(2)).               (12.3)
```

Because A_1(z) is polynomial in the rescaled cell parameter, its C^2(K1) norm is a finite explicit polynomial coefficient bound. Substitution into (3.6), (5.5), (6.3), and (8.2) gives the numerical branch threshold required by the downstream K_infinity computation.

---

# 13. Clean quotable interface

## Uniform hyperbolic conjugation theorem with 1/n drift

Let K subset int(K1) subset Omega be compact parameter sets. Let A_n:Omega->GL_2(C) be holomorphic and suppose that for n>=Nexp

```plain text
||A_n-A_0-A_1/n||_C2(K1) <= M_2/n^2.
```

Suppose there is a holomorphic S:Omega->GL_2(C) near K1 such that

```plain text
S^(-1)A_0S=diag(lambda_+,lambda_-),
inf_K1 |lambda_+lambda_-|>0,
sup_K1 |lambda_-/lambda_+|<1.
```

Define

```plain text
C=S^(-1)A_1S,
d_i=C_ii/lambda_i.
```

Choose the weighted C^2 norm and the constants

```plain text
qbar,H0,L,d0,q1,CG,Ainv,Rc,Nstar
```

by equations (1.2), (3.6)-(3.7), (4.1)-(4.3), (5.5), (6.3), and (8.2). Then for every integer n0>=Nstar there exist holomorphic invertible matrices P_n(z), n>=n0, such that

```plain text
P_{n+1}^(-1)A_nP_n
 = diag(lambda_+(1+d_+/n),
        lambda_-(1+d_-/n))                   (13.1)
```

for every z in K1 and every n>=n0, and

```plain text
max_{0<=j<=2} sup_{z in K}
 ||partial_z^j(P_n(z)-S(z))||
 <= 2 eps^(-2)CP/n.                          (13.2)
```

There are two holomorphic fundamental solutions Y_n^+,- satisfying

```plain text
max_{0<=j<=2} sup_{z in K}
 |partial_z^j[
   lambda_i(z)^(-n)n^(-d_i(z))Y_n^i(z)-S(z)e_i
 ]|
 <= 2 eps^(-2)Casym/n,                       (13.3)
```

where

```plain text
CP=S0(H0+2Rc+2H0Rc/n0),
Casym=e(CP+2(d0^2+d0/2)S0).                  (13.4)
```

In particular,

```plain text
Y_n^i(z)
 = lambda_i(z)^n n^d_i(z)
   (S(z)e_i+O_C2,K(1/n))                     (13.5)
```

uniformly for z in K. The relative stable/unstable ratio contracts by q1^(n-m) for every n>=m>=n0.

This proves the complete uniform conjugation lemma required by [ARCH-CELL], including the polynomial factors, two parameter derivatives, an explicit contraction threshold, and exact dependence on the spectral and degeneration constants.