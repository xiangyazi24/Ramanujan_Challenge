ANSWER Q6723 985c02e6

# K_infinity branch-table and h0 certification procedure

# 0. Executive verdict

The downstream computation is finite and can be implemented now with SymPy plus complex ball arithmetic. The Apéry specialization is substantially simpler than the general conjugation theorem:

```plain text
A_infinity = [[34,-1],[1,0]],
A_1        = [[-51,3],[0,0]],
lambda_+   = 17+12sqrt(2),
lambda_-   = 17-12sqrt(2),
q          = lambda_-/lambda_+ = lambda_-^2.
```

In the eigenbasis of A_infinity, both normalized diagonal drifts are exactly

```plain text
d_+ = d_- = -3/2.
```

Thus the polynomial factor

```plain text
g_n = product_(k=n0)^(n-1) (1-3/(2k))
```

is common to both hyperbolic modes and independent of z. It cancels in the ratios F_n(z)/b_n and G_n(z)/b_n. Consequently the limiting cell object really is

```plain text
K_infinity(z)
 = phi(-z) phi(z) + z^6 gamma(-z) gamma(z),
```

with no hidden n^(d(z)-d(0)) factor.

There are two distinct numerical tasks.

1. Local branch certification: isolate all zeros of

```plain text
H_infinity(z) := z K_infinity'(z)-3K_infinity(z)
```

in the standard cell and certify their critical values.

1. Global cell squeeze: after multiplying by the exact pole scales b_r b_s, certify that at least h mirror-orbits of squared critical values are simple, nonzero, and globally distinct.

The first task is routine interval Newton plus the conjugation error bound. The second is the decisive test. A bare O(1/r+1/s) error is not enough near the central cells, where adjacent same-branch scales differ by O(h^-2). The certified program must therefore use either:

```plain text
(a) direct finite-r,s interval evaluation for the requested h, or
(b) a second-order inverse-length branch expansion with an O((r+1)^-3+(s+1)^-3) remainder.
```

The algorithm below includes both modes. No numerical value of h0 is claimed here; h0 is the integer emitted by the interval run.

---

# 1. Exact Apéry transfer data

The shifted cell recurrence is

```plain text
(n+z+1)^3 Y_(n+1)
 = P(n+z)Y_n-(n+z)^3Y_(n-1),
P(t)=34t^3+51t^2+27t+5.
```

For the state vector (Y_n,Y_(n-1))^T,

```plain text
A(n,z) =
[ P(n+z)/(n+z+1)^3    -(n+z)^3/(n+z+1)^3 ]
[ 1                    0                     ].
```

The exact expansion is

```plain text
A(n,z)=A_infinity+A_1/n+R_n(z),
```

with

```plain text
A_infinity = [[34,-1],[1,0]],
A_1        = [[-51,3],[0,0]].                 (1.1)
```

Put w=1/n. Then n^2 R_n(z)=E(w,z), where only the first row is nonzero:

```plain text
E_11(w,z) =
 ((51z+78)(1+zw)^2
  +(75z+124)w(1+zw)
  +(29z+51)w^2)
 / (1+(z+1)w)^3,

E_12(w,z) =
-((3z+6)(1+zw)^2
  +(3z+8)w(1+zw)
  +(z+3)w^2)
 / (1+(z+1)w)^3.                              (1.2)
```

These formulas are the preferred input for interval bounds. They avoid taking a supremum over infinitely many integers n.

Let K1 be the enlarged compact cell used by the conjugation theorem and set

```plain text
Z1 = sup_(z in K1) |z+1|,
Nexp = max(2,ceil(2Z1)).                       (1.3)
```

Then for 0<=w<=1/Nexp,

```plain text
|1+(z+1)w| >= 1/2.
```

Define the exact interval-computable remainder constant

```plain text
M2(K1,Nexp)
 := max_(a=0,1,2)
    sup_(z in K1, 0<=w<=1/Nexp)
    ||partial_z^a (S^-1 E(w,z) S)||.          (1.4)
```

The norm is the weighted C^2 row-sum norm from Q6708. Since A_1 is constant,

```plain text
||A_1||_C2(K1) = 54,
partial_z A_1 = partial_z^2 A_1 = 0.          (1.5)
```

Thus the requested per-cell C^2 norm of A_1(z) is the same for every cell. Only M2(K1,Nexp) depends on the chosen cell box.

---

# 2. Exact hyperbolic constants and the specialized n0(K)

Let

```plain text
Delta = lambda_+-lambda_- = 24sqrt(2),
S = [[lambda_+,lambda_-],[1,1]],
S^-1 = (1/Delta)[[1,-lambda_-],[-1,lambda_+]].
```

Then

```plain text
C := S^-1 A_1 S
   = (3/2)
     [[-lambda_+, lambda_-],
      [ lambda_+,-lambda_-]].                 (2.1)
```

Therefore

```plain text
d_+ = C_++/lambda_+ = -3/2,
d_- = C_--/lambda_- = -3/2.                  (2.2)
```

The off-diagonal homological conjugator is

```plain text
H =
[[0, -3lambda_-/(2Delta)],
 [3lambda_+/(2Delta), 0]].                    (2.3)
```

With the row-sum matrix norm, record

```plain text
C1    = ||C|| = 51,
H0    = ||H|| = 3lambda_+/(2Delta),
L0    = ||diag(lambda_+,lambda_-)|| = lambda_+,
Ainv  = 2lambda_+,
q     = lambda_-^2 approximately 8.6655e-4.  (2.4)
```

Because d_+=d_-, the diagonal ratio is exactly q; no Cr/n enlargement is needed. Use

```plain text
CG = 2 + 2/(1-q) + 1/(1-q)^2.                (2.5)
```

Let

```plain text
E2 = M2(K1,Nexp),
B0 = E2(1+H0)+2C1H0+L0H0,
L  = 2B0,
Rc = 2 CG Ainv L,
Bd = (3/2)^2+(3/2)/2 = 3.                    (2.6)
```

A conservative specialized threshold is

```plain text
n0(K1) = ceil max(
  Nexp,
  4,
  2H0,
  2Rc,
  2Bd
).                                            (2.7)
```

This is slightly stronger than necessary but simultaneously guarantees:

```plain text
I+H/n is invertible,
I+U_n is invertible,
the fixed-point contraction is <=1/2,
all diagonal drift factors are nonzero,
and the product-to-n^(-3/2) comparison is valid.
```

The exact common diagonal product is

```plain text
p_+ (n) = lambda_+^(n-n0)
          product_(k=n0)^(n-1)(1-3/(2k)),

p_- (n) = lambda_-^(n-n0)
          product_(k=n0)^(n-1)(1-3/(2k)).     (2.8)
```

Equivalently,

```plain text
product_(k=n0)^(n-1)(1-3/(2k))
 = Gamma(n-3/2)Gamma(n0)
   /(Gamma(n)Gamma(n0-3/2)).                  (2.9)
```

This is the normalization correction that must be used in every branch computation.

---

# 3. Exact cell functions

Define

```plain text
F_0(z)=1,
F_1(z)=P(z)/(z+1)^3,
G_0(z)=0,
G_1(z)=1/(z+1)^3,
```

and propagate both by the shifted recurrence. Then

```plain text
b_n = F_n(0),
phi_n(z)   = F_n(z)/b_n,
gamma_n(z) = G_n(z)/b_n.                      (3.1)
```

For a pole split

```plain text
r=j-1,
s=h-j,
r+s=h-1,
```

the exact normalized cell function is

```plain text
K_(r,s)(z)
 = phi_r(-z)phi_s(z)
   +z^6 gamma_r(-z)gamma_s(z).                (3.2)
```

The original rational function in that cell is

```plain text
f_h(-j+z)
 = b_r b_s z^(-3) K_(r,s)(z).                (3.3)
```

Define

```plain text
Phi_(r,s)(z) = K_(r,s)(z)/z^3,
H_(r,s)(z)   = z K_(r,s)'(z)-3K_(r,s)(z).    (3.4)
```

Away from z=0, the critical points are exactly the zeros of H_(r,s). If xi is such a zero, the corresponding squared critical value, i.e. a root candidate of the mirror quotient W_h, is

```plain text
U_(r,s)(xi)
 = (b_r b_s)^2
   (K_(r,s)(xi)/xi^3)^2.                     (3.5)
```

Reflection acts by

```plain text
(r,s,z,t) <-> (s,r,-z,-t),
```

so the two reflected critical points give the same U=t^2. The counting routine must quotient by this involution before comparing with the target h.

---

# 4. Computing phi, gamma, and K_infinity from the exact conjugator

Let P_n(z) be the exact conjugating matrix from Q6708, normalized so that

```plain text
P_(n+1)^(-1) A(n,z) P_n
 = diag(lambda_+(1-3/(2n)),
        lambda_-(1-3/(2n))).                 (4.1)
```

At n=n0, compute the state vectors

```plain text
Y_F(z) = (F_n0(z),F_(n0-1)(z))^T,
Y_G(z) = (G_n0(z),G_(n0-1)(z))^T,
Y_b    = Y_F(0).
```

Put

```plain text
a_F(z)=P_n0(z)^(-1)Y_F(z),
a_G(z)=P_n0(z)^(-1)Y_G(z),
a_b   =P_n0(0)^(-1)Y_b.                      (4.2)
```

Certify 0 notin (a_b)_+. Then

```plain text
phi(z)   = (a_F(z))_+/(a_b)_+,
gamma(z) = (a_G(z))_+/(a_b)_+.               (4.3)
```

The common n^(-3/2) factor and the first component of the limiting eigenvector cancel exactly. Hence

```plain text
K_infinity(z)
 = phi(-z)phi(z)+z^6gamma(-z)gamma(z).        (4.4)
```

To compute P_n0 rigorously, iterate the fixed-point map from Q6708 in the weighted C^2 sequence norm. If the contraction factor is kappa<1, an iterate Q^(m) has the a posteriori enclosure

```plain text
||Q-Q^(m)|| <= ||F(Q^(m))-Q^(m)||/(1-kappa). (4.5)
```

This gives proof-grade complex balls for P_n0, phi, gamma, and their first two z derivatives.

A simpler independent check is to compute F_N/b_N and G_N/b_N for a large N and enlarge their C^2 balls by the explicit conjugation tail. The amplitude method (4.2)-(4.3) is preferable because it does not rely on cancellation of huge numbers.

---

# 5. The finite branch table

Choose a standard compact cell Kcell inside the strip

```plain text
-1 < Re(z) < 0
```

and an enlarged compact K1 with positive distance from z=0, z=-1, and the shifted singular lattice. The outer contour of Kcell must contain every finite critical point assigned to that cell.

There are two kinds of templates.

## 5.1 Bulk template

```plain text
K_inf,inf(z) = K_infinity(z).
```

It controls cells with r,s>=n0.

## 5.2 Edge templates

For every fixed m=0,...,n0-1, compute

```plain text
K_(m,inf)(z)
 = phi_m(-z)phi(z)
   +z^6gamma_m(-z)gamma(z),                  (5.1)
```

where phi_m=F_m/b_m and gamma_m=G_m/b_m are exact rational functions. The reflected template K_(inf,m) is obtained by (m,inf,z)->(inf,m,-z) and does not need a separate numerical solve.

This finite list is essential. Discarding the edge cells loses critical-value orbits and does not by itself reach the [CRIT-2H] count.

## 5.3 Record stored for each template branch

For each template T, define

```plain text
H_T(z)=zK_T'(z)-3K_T(z),
Phi_T(z)=K_T(z)/z^3.
```

The branch table must store:

```plain text
1. A complex box B_(T,nu) containing one zero xi_(T,nu) of H_T.
2. A proof that the zero is unique and simple.
3. The ball Xi_(T,nu) for xi_(T,nu).
4. The critical-value ball Tau_(T,nu)=Phi_T(Xi_(T,nu)).
5. mu_(T,nu)=inf_B |H_T'(z)|.
6. eta_(T,nu)=inf_(boundary B) |H_T(z)|.
7. M2Phi_(T,nu)=sup_B |Phi_T''(z)|.
8. The distance from Tau_(T,nu)^2 to zero.
9. The pairwise separation matrix for the squared Tau balls.
10. An outer-contour zero count proving that the listed boxes exhaust H_T in the cell.
11. The selected mirror-orbit label used in the CRIT-2H count.
```

Use complex interval Newton or a Krawczyk operator. A sufficient test is

```plain text
N(B)=c-H_T(c)/H_T'(B) subset interior(B),
0 notin H_T'(B).                              (5.2)
```

For full exhaustiveness, certify the winding number of H_T on the outer contour and check that it equals the sum of the box multiplicities.

---

# 6. Second-order branch constants

The K_infinity table alone certifies local persistence and cross-branch separation at fixed scale. It does not resolve the O(h^-2) central same-branch gaps from an O(h^-1) remainder.

Extend the normal-form calculation by two orders. The rational transfer has a complete expansion in 1/n; solve the homological equations recursively to obtain

```plain text
phi_n(z)
 = phi(z)+phi_1(z)/n+phi_2(z)/n^2+O_C2(n^-3),

gamma_n(z)
 = gamma(z)+gamma_1(z)/n+gamma_2(z)/n^2
   +O_C2(n^-3).                               (6.1)
```

This is a finite symbolic computation. At each order, the off-diagonal coefficient is obtained by dividing by lambda_+-lambda_-; the diagonal coefficient is absorbed into the scalar product. The residual O(n^-4) system is handled by the same fixed-point theorem.

For the bulk template put

```plain text
u=1/r,
v=1/s.
```

Expand

```plain text
K_(r,s)(z)
 = K_00(z)+uK_10(z)+vK_01(z)
   +u^2K_20(z)+uvK_11(z)+v^2K_02(z)
   +R_K(z;u,v),                               (6.2)
```

with a certified bound

```plain text
||R_K||_C2
 <= C_K3 (u^3+u^2v+uv^2+v^3).               (6.3)
```

For an edge template m,inf, use the one-variable version in v=1/s.

For each simple limiting branch, implicit differentiation of

```plain text
H(z,u,v)=zK_z-3K=0
```

gives the branch location coefficients. For example,

```plain text
xi_u = -H_u/H_z,
xi_v = -H_v/H_z.                              (6.4)
```

Since Phi_z=0 at a critical point, the first critical-value derivatives are especially simple:

```plain text
Tau_u = Phi_u,
Tau_v = Phi_v.                                (6.5)
```

The second derivatives are

```plain text
Tau_uu = Phi_uu+2Phi_zu xi_u+Phi_zz xi_u^2,
Tau_uv = Phi_uv+Phi_zu xi_v+Phi_zv xi_u
         +Phi_zz xi_u xi_v,
Tau_vv = Phi_vv+2Phi_zv xi_v+Phi_zz xi_v^2.  (6.6)
```

Store these six branch constants and the cubic remainder in the table. They are the data that decide the central-cell ordering.

Also compute the Apéry scale-ratio expansion, or equivalently a certified bound for

```plain text
rho_n = b_(n+1)/b_n,
Mscale = sup_(n>=n0)
 n^3 |log(rho_n/lambda_+)
      +3/(2n)-beta_2/n^2|.                   (6.7)
```

The coefficient beta_2 and the remainder follow from the same normal-form expansion. SymPy can solve these coefficients by substituting a formal 1/n series into the recurrence.

---

# 7. Certified YES/NO algorithm for a fixed h

The direct fixed-h algorithm is exact and is the best regression test.

```plain text
function CERTIFY_H(h, precision):
    generate exact integers b_0,...,b_(h-1)
    generate exact rational jets F_n,G_n for every required cell box

    branch_records = []

    for each canonical split r+s=h-1 with r<=s:
        build K_(r,s), H_(r,s), Phi_(r,s)
        choose starting boxes from the matching bulk or edge template

        for every template branch nu:
            apply complex interval Newton to H_(r,s) on B_(T,nu)
            if uniqueness or simplicity fails:
                return NO("branch not certified")

            Xi = certified root ball
            Tau = Phi_(r,s)(Xi)
            U = (b_r*b_s)^2 * Tau^2
            store (r,s,nu,Xi,Tau,U)

        certify the cell root count by a contour argument
        if the boxes are not exhaustive:
            return NO("unaccounted critical point")

    quotient branch_records by
        (r,s,z,Tau) ~ (s,r,-z,-Tau)

    build the overlap graph of the U-balls
    mark a mirror orbit GOOD iff
        its critical point is simple,
        0 notin U,
        and its U-ball is disjoint from every nonmirror U-ball

    if number of GOOD mirror orbits >= h:
        return YES
    else:
        return NO("cell-squeeze margin not certified")
```

Here NO means not certified by this interval run, not that [CRIT-2H] is false. Increase precision or subdivide boxes before treating it as a structural failure.

The logical implication is

```plain text
CERTIFY_H(h)=YES
 => W_h has at least h simple nonzero roots
 => s_h>=2h
 => [CRIT-2H] at h.                          (7.1)
```

The parity convention in the project may state the intermediate bound as s_h>=2h-1; since s_h is even, the two formulations are equivalent at this step.

---

# 8. Exact interval inequalities

The program should use balls directly. The following scalar inequalities are useful for diagnostics and for the asymptotic tail proof.

## 8.1 Branch persistence

Suppose a template branch box B has

```plain text
eta = inf_(boundary B)|H_T|,
mu  = inf_B |H_T'|,
R   = sup_B |z|,
r0  = inf_B |z|.
```

If

```plain text
||K_(r,s)-K_T||_C2 <= delta_K,                (8.1)
```

then

```plain text
||H_(r,s)-H_T||_C0 <= (R+3)delta_K,
||H_(r,s)'-H_T'||_C0 <= (R+2)delta_K.         (8.2)
```

A sufficient persistence and simplicity condition is

```plain text
(R+3)delta_K < eta,
(R+2)delta_K < mu.                            (8.3)
```

Put

```plain text
m = mu-(R+2)delta_K > 0,
delta_xi = (R+3)delta_K/m.                    (8.4)
```

Also require delta_xi to be smaller than the certified distance from the limiting root ball to the boundary of B.

## 8.2 Critical-value radius

Let

```plain text
delta_Phi0 = delta_K/r0^3,
delta_Phi1 = (R+3)delta_K/r0^4.
```

A certified value radius is

```plain text
epsilon_Tau
 = delta_Phi0
   +delta_Phi1 delta_xi
   +(1/2)M2Phi delta_xi^2.                    (8.5)
```

The second-order table from Section 6 replaces the crude delta_K by the cubic remainder when central-cell resolution is needed.

## 8.3 Global squared-value separation

Let branch a have center-radius enclosure

```plain text
Tau_a subset Ball(c_a,e_a),
```

and similarly for branch b. Set

```plain text
E_a = 2|c_a|e_a+e_a^2,
E_b = 2|c_b|e_b+e_b^2.                        (8.6)
```

To avoid huge Apéry integers, divide by the first pole scale and put

```plain text
xi_ab = ((b_(r_b)b_(s_b))/(b_(r_a)b_(s_a)))^2. (8.7)
```

The exact interval separation inequality is

```plain text
lower_abs(c_a^2-xi_ab c_b^2)
   > E_a+xi_ab E_b.                           (8.8)
```

Equivalently,

```plain text
0 notin Tau_a^2-xi_ab Tau_b^2.               (8.9)
```

This is the load-bearing numerical inequality. It must be checked for every selected branch against every nonmirror branch.

Nonvanishing is certified by

```plain text
lower_abs(c_a) > e_a.                         (8.10)
```

For adjacent cells on the same branch, use the exact scale ratio

```plain text
((b_(r+1)b_(s-1))/(b_r b_s))^2
 = (rho_r/rho_(s-1))^2                       (8.11)
```

and the second-order Tau model. The signed interval

```plain text
Delta_same
 = (rho_r/rho_(s-1))^2
   Tau_nu(1/(r+1),1/(s-1))^2
   -Tau_nu(1/r,1/s)^2                        (8.12)
```

must exclude zero. This is the exact check that resolves the central O(h^-2) gap.

---

# 9. Tail certification and the definition of h0

The all-h proof uses a finite template computation plus compact interval branch-and-bound.

For each integer H>=2n0+1, define the following certified conditions.

```plain text
PERSIST(H):
  inequalities (8.3) hold for every bulk branch whenever r,s>=n0
  and r+s>=H-1; the analogous edge inequalities hold for every m<n0.

EXHAUST(H):
  the persistent branch boxes account for all critical points in every
  bulk and edge cell.

NONZERO(H):
  every selected squared-value ball excludes zero.

SAME(H):
  the adjacent same-branch intervals (8.12) exclude zero for every
  admissible r,s with r+s>=H-1.

CROSS(H):
  inequality (8.8) holds for every selected nonmirror pair, including
  bulk-bulk, bulk-edge, and edge-edge pairs.

COUNT(H):
  after mirror quotienting, the branch inventory supplies at least h
  singleton squared-value orbits for every h>=H.
```

The universal checks are finite after the following reduction.

1. Use u=1/r, v=1/s and the second-order Taylor models with certified cubic remainder.

1. For scale ratios use the formal expansion of rho_n with certified cubic remainder.

1. Subdivide the compact inverse-length domains by interval branch-and-bound.

1. For CROSS, first compute a pilot ordering/nearest-neighbor graph at high precision. Prove that every adjacent gap in that graph stays nonzero on the parameter boxes. If the graph can change, subdivide at the possible change locus. If no stable finite graph is certified, the method returns NO.

1. Retain integer constraints when they improve a box: r,s,r',s' are integers and r+s=r'+s'=h-1.

Define

```plain text
h0 := least H such that
  PERSIST(H) and EXHAUST(H) and NONZERO(H)
  and SAME(H) and CROSS(H) and COUNT(H)        (9.1)
```

are all interval-certified.

The program should emit both h0 and a certificate file containing:

```plain text
all template boxes,
all interval-Newton inclusions,
all contour counts,
all Taylor coefficients and remainders,
all pair-separation margins,
and the final universal-domain subdivisions.
```

---

# 10. Exact closure logic with the banked h<=40 certificates

The finite computation that closes the campaign is precisely:

```plain text
(A) certify the complete bulk and edge branch table;
(B) certify all predicates in (9.1);
(C) obtain h0<=40.                            (10.1)
```

If the interval program returns h0<=40, then:

```plain text
2<=h<=40:
  [CRIT-2H] follows from the banked exact full-Morse certificates.

h>=41:
  h>h0, so [CRIT-2H] follows from the certified cell-squeeze theorem.
```

Therefore [CRIT-2H] holds for every h>=2.

Running CERTIFY_H(h) for 2<=h<=40 is a valuable regression test, but it is not logically required once the banked certificates are accepted. The genuinely new finite output is the proof that h0<=40.

If the computed h0 is larger than 40, nothing is disproved. One must either:

```plain text
extend the exact machine certificates through h0,
sharpen the conjugation/remainder constants,
use higher-order inverse-length expansions,
or subdivide the cell and pair domains more finely.
```

---

# 11. SymPy/mpmath-style pseudocode

```plain text
function transfer_remainder_ball(K1, Nexp):
    W = interval(0,1/Nexp)
    evaluate E_11,E_12 and z-derivatives up to order 2
    E2 = norm(S^-1 E S)
    return E2

function conjugation_constants(K1):
    Nexp = max(2,ceil(2*sup_abs(z+1,K1)))
    E2 = transfer_remainder_ball(K1,Nexp)
    compute H0,C1,L0,Ainv,q,CG,B0,L,Rc
    n0 = ceil(max(Nexp,4,2*H0,2*Rc,6))
    return constants,n0

function fixed_point_conjugator(K1,n0,prec):
    initialize U_n=0 on n0<=n<=Ntail
    iterate the two Green fixed-point maps
    use residual/(1-kappa) for the tail enclosure
    return C2-ball function P_n0(z)

function limit_functions(K1,n0):
    propagate F,G and their z-jets to n0
    P = fixed_point_conjugator(K1,n0)
    aF = inverse(P)*state(F,n0)
    aG = inverse(P)*state(G,n0)
    ab = inverse(P(0))*state(F(0),n0)
    assert 0 notin ab.plus
    phi   = aF.plus/ab.plus
    gamma = aG.plus/ab.plus
    return phi,gamma

function certify_template(T,Kcell,K1):
    construct K_T,H_T,Phi_T as C2 ball functions
    locate floating roots of H_T
    inflate to disjoint complex boxes
    run interval Newton on every box
    certify outer-contour root count
    compute Xi,Tau,mu,eta,M2Phi
    compute second-order inverse-length branch coefficients
    return TemplateRecord

function certify_h(h,records):
    execute CERTIFY_H from Section 7

function certify_tail(records,H):
    certify PERSIST,EXHAUST,NONZERO,SAME,CROSS,COUNT
    on all compact parameter boxes
    return YES or NO with failed margin

function find_h0(records):
    for H from 2*n0+1 upward:
        if certify_tail(records,H)==YES:
            return H
```

For exploratory work, sympy plus mpmath.iv is sufficient. For the final certificate, use outward-rounded complex balls such as Arb/Acb through python-flint; ordinary floating-point root finding is not a proof.

---

# 12. Failure modes

The program must report which of the following occurred.

```plain text
F1. K1 intersects a pole or the denominator ball contains zero.
F2. The remainder constant makes the fixed-point contraction >=1.
F3. The dominant Apéry amplitude ball contains zero.
F4. Interval Newton fails to isolate a limiting branch.
F5. The outer contour count reveals an unaccounted critical point.
F6. A branch derivative ball contains zero: simplicity is unresolved.
F7. A critical-value ball contains zero.
F8. Two nonmirror squared-value balls overlap.
F9. The O(1/n) model resolves cross-branch gaps but not the central
    same-branch O(h^-2) gap; increase the expansion order.
F10. The pilot global branch ordering changes on the parameter domain,
     so no finite stable neighbor graph has been certified.
F11. Edge templates fail to supply the required mirror-orbit count.
F12. The resulting h0 is greater than 40.
F13. Precision or interval dependency blow-up is too large; subdivide.
```

Failures F4-F13 are failures of the certificate, not counterexamples to [CRIT-2H].

---

# 13. Quotable interface

```plain text
THEOREM SCHEMA [APERY-CELL-SQUEEZE].

Let Kcell subset int(K1) be compact standard cell domains. Construct the
bulk template K_infinity and the finite edge templates K_(m,infinity),
0<=m<n0(K1), using the exact Apéry transfer and the uniform C2
conjugation theorem. Suppose interval certificates establish:

1. exhaustive simple nonzero branch tables for every template;
2. the second-order inverse-length expansions with certified cubic tails;
3. persistence and exhaustiveness for every finite cell with h>h0;
4. nonvanishing and the pairwise separation inequality

   lower_abs(c_a^2-xi_ab c_b^2)
     > (2|c_a|e_a+e_a^2)
       +xi_ab(2|c_b|e_b+e_b^2)

   for every nonmirror branch pair;
5. after mirror quotienting, at least h singleton squared-value orbits
   for every h>h0.

Then W_h has at least h simple nonzero roots for every h>h0; hence
[CRIT-2H] holds for every h>h0.

If additionally h0<=40, the banked full-Morse certificates for 2<=h<=40
imply [CRIT-2H] for every h>=2.
```