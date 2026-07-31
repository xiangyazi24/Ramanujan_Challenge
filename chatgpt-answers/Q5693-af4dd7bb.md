ANSWER Q5693 af4dd7bb

# Full Cartier packets for Apéry P3.2: exact PDE defect, Smith saturation, and exposed-face range

## Executive verdict

I audited the current repository state at main commit [`256a37c42c864e0770ea7ebafe210c5729915b4f`](https://github.com/xiangyazi24/Ramanujan_Challenge/commit/256a37c42c864e0770ea7ebafe210c5729915b4f), especially the fixed-moment torus reformulation in Section 48 of `problems/3.2/research/working_notes/Q32_SEPARATION_ANALYSIS.md`, the coefficient-shell audit, and the earlier Fitting/Cartier saturation results.

The full Cartier packet is strictly stronger data than the torus shell, but it still does **not** add target rank.

Let

\[
\Lambda(x,y,z)=
\frac{(1+x)(1+y)(1+z)\bigl((1+y)(1+z)+xyz\bigr)}{xyz},
\qquad
\Lambda^m=\sum_{v\in\mathbf Z^3}c_m(v)X^v,
\]

where \(X^v=x^{v_1}y^{v_2}z^{v_3}\).  For

\[
n=aq+r,\qquad 0\le r<q,
\]

Cartier gives

\[
\boxed{
\mathcal C_q(\Lambda^n)=b_r\Lambda^a,
\qquad
c_n(q\mu)\equiv b_r c_a(\mu)\pmod q.
}
\tag{E1}
\]

The conclusions are:

1. The support polytope \(P=\operatorname{Newt}(\Lambda)\) has the integer-decomposition property, and
   \[
   \operatorname{Supp}(\Lambda^a)=aP\cap\mathbf Z^3.
   \]
   Its packet size is exactly
   \[
   \boxed{
   N_a=\#(aP\cap\mathbf Z^3)
   =\frac{38a^3+57a^2+31a+6}{6}.
   }
   \tag{E2}
   \]

2. The packet coefficient vector
   \[
   p_a=(c_a(\mu))_{\mu\in aP\cap\mathbf Z^3}
   \]
   is primitive: every vertex coefficient is \(1\).  Therefore the full congruence packet in (E1) has Smith content exactly one scalar:
   \[
   \boxed{
   \bigl(c_n(q\mu):\mu\in aP\cap\mathbf Z^3\bigr)_{\mathbf F_q}
   =(b_r)_{\mathbf F_q}.
   }
   \tag{E3}
   \]
   All \(N_a\) packet entries vanish if and only if the one scalar \(b_r\) vanishes.  They do not give \(N_a\) independent equations.

3. The exact logarithmic-derivative coefficient PDEs have a characteristic-\(q\) Frobenius defect of dimension \(N_a\).  In fact the complete Laurent-polynomial solution space with Newton support in \(nP\) is
   \[
   \boxed{
   \mathscr S_{q,n}
   =\{\Lambda^r H(x^q,y^q,z^q):
       \operatorname{Newt}(H)\subseteq aP\},
   \qquad
   \dim_{\mathbf F_q}\mathscr S_{q,n}=N_a.
   }
   \tag{E4}
   \]
   If \(q\mid b_r\), the whole \(q\)-packet vanishes for **every** \(H\) in this \(N_a\)-dimensional space.  Thus the packet imposes no condition on the Frobenius-defect coordinates.

4. For a nearby prime \(\ell=q+h\), with \(2ah<q\), the nonzero packet nodes \(q(aP\cap\mathbf Z^3)\) and \(\ell(aP\cap\mathbf Z^3)\) are disjoint.  Modulo \(q\), the \(\ell\)-packet is the value of an explicit linear map on the free \(H\)-coordinates; the \(q\)-target packet does not determine it.  The common origin is the sole shared node and yields the tall carrier \(c_n(0)=b_n\).

5. Exposed faces do not become rigid.  If \(F\) is an exposed face of \(P\), the surviving face defect has dimension
   \[
   \#(aF\cap\mathbf Z^3).
   \]
   If \(ah>r\), the corresponding exposed face of the \(\ell\)-packet lies outside \(nP\) and is identically zero over \(\mathbf Z\); this is universal support geometry, not target selectivity.  If \(ah\le r\), the face survives and keeps its full face-sized Frobenius defect.

6. There is one exact positive boundary range.  At the vertex \(u=(1,1,1)\), put \(t=r-ah\).  For \(t\ge0\),
   \[
   c_n(\ell a u)
   =V_n(t)
   :=\sum_{s=0}^{t}
      \binom ns\binom n{t-s}\binom{n+s}{t}^2,
   \qquad
   \log V_n(t)\ll(t+1)\log n.
   \tag{E5}
   \]
   Hence a close direct-branch target pair has the nonzero carrier
   \[
   V_n(r)V_n(r-ah),
   \]
   divisible by \(q\ell\), of logarithmic height \(O(r\log n)\).  This is useful only when \(r=o(H/\log n)\), a boundary range already contributing \(o(H/\log n)\) candidates.  The same packet gives no analogous gain from a small reflected index \(q-1-r\), because every exposed packet face is at normal distance \(r\), not \(q-1-r\), from \(nP\).

7. Consequently no nonzero characteristic-zero scalar is produced by eliminating the packet and the coefficient PDEs alone.  For two characteristics the residual defect is
   \[
   \boxed{2N_a}
   \]
   (or \(2N_a-2\) after fixing one scale coordinate in each characteristic), independent of \(h\) as long as \(2ah<\min(q,\ell)\).  Any successful corridor determinant must use an additional identity selecting the distinguished defect coordinate \(H=\Lambda^a\); it cannot follow from Cartier packet vanishing, the finite-stencil PDEs, root filters, or exposed-face triangularity alone.

The smallest remaining lemma is therefore still a cross-characteristic anti-clustering statement for the **single Cartier scalar** \(b_r\), not a high-rank packet theorem.  A constructive version is stated in Section 11.

---

## 1. The correct support and Newton polytope

The Laurent polynomial has the useful decomposition

\[
\Lambda
=x^{-1}y^{-1}z^{-1}(1+x)(1+y)^2(1+z)^2
 +(1+x)(1+y)(1+z).
\tag{1.1}
\]

Thus its support is the union

\[
\mathcal A_1=\{-1,0\}\times\{-1,0,1\}^2,
\qquad
\mathcal A_2=\{0,1\}^3,
\tag{1.2}
\]

with four overlaps, hence \(18+8-4=22\) support points.  All coefficients are positive.  Explicitly, the first box has coefficient

\[
\binom2{u_2+1}\binom2{u_3+1},
\]

and the second cube contributes one additional unit on its eight points.

The Newton polytope is

\[
\boxed{
P=\{(x,y,z)\in\mathbf R^3:
-1\le x,y,z\le1,\ x-y\le1,\ x-z\le1\}.
}
\tag{1.3}
\]

Its lattice points are exactly the 22 support points in (1.2).

### Lemma 1.1 — integer decomposition

For every \(a\ge0\),

\[
\operatorname{Supp}(\Lambda^a)=aP\cap\mathbf Z^3.
\tag{1.4}
\]

### Proof

Because all coefficients of \(\Lambda\) are positive, it is enough to prove that every lattice point of \(aP\) is a sum of \(a\) lattice points of \(P\).

Let \((m,y,z)\in aP\cap\mathbf Z^3\).

- If \(m\ge0\), choose \(m\) of the \(a\) first coordinates equal to \(1\), and the remaining \(a-m\) equal to \(0\).  At a coordinate with first entry \(1\), the second and third entries may be \(0\) or \(1\); at a coordinate with first entry \(0\), they may be \(-1,0,1\).  The attainable range for the second-coordinate sum is therefore
  \[
  m-a\le y\le a,
  \]
  exactly the inequality \(x-y\le a\); every integer in this interval is attained.  The third coordinate is independent and identical.

- If \(m<0\), choose \(-m\) first coordinates equal to \(-1\), and all remaining first coordinates equal to \(0\).  There is then no additional restriction on the second and third coordinates, and every pair \(y,z\in[-a,a]\cap\mathbf Z\) is attained.

This proves (1.4). \(\square\)

For a fixed first coordinate \(m\), the numbers of allowed second and third coordinates are

\[
2a+1\quad(-a\le m\le0),
\qquad
2a-m+1\quad(1\le m\le a).
\]

Therefore

\[
\begin{aligned}
N_a
&=(a+1)(2a+1)^2+
  \sum_{m=1}^{a}(2a-m+1)^2\\
&=(a+1)(2a+1)^2+
  \sum_{t=a+1}^{2a}t^2\\
&=\frac{38a^3+57a^2+31a+6}{6}.
\end{aligned}
\tag{1.5}
\]

For \(a=1\), this gives \(N_1=22\), as required.

---

## 2. The full Cartier packet and its primitive content

Define the Cartier operator

\[
\mathcal C_q\!\left(\sum_v d(v)X^v\right)
=\sum_{\mu\in\mathbf Z^3}d(q\mu)X^\mu.
\tag{2.1}
\]

In characteristic \(q\), the freshman's-dream identity gives

\[
\Lambda^{aq+r}
=\Lambda^r(\Lambda^q)^a
\equiv\Lambda^r\Lambda(x^q,y^q,z^q)^a.
\tag{2.2}
\]

Cartier satisfies

\[
\mathcal C_q(F(x^q,y^q,z^q)G)=F(x,y,z)\mathcal C_q(G).
\tag{2.3}
\]

Since \(0\le r<q\) and

\[
\operatorname{Supp}(\Lambda^r)\subseteq[-r,r]^3,
\]

the only exponent of \(\Lambda^r\) divisible coordinatewise by \(q\) is \(0\).  Hence

\[
\mathcal C_q(\Lambda^r)=\operatorname{CT}\Lambda^r=b_r.
\tag{2.4}
\]

Combining (2.2)--(2.4) proves

\[
\boxed{
\mathcal C_q(\Lambda^n)=b_r\Lambda^a.
}
\tag{2.5}
\]

Equivalently, for every \(\mu\in aP\cap\mathbf Z^3\),

\[
\boxed{
c_n(q\mu)\equiv b_r c_a(\mu)\pmod q.
}
\tag{2.6}
\]

This formula is valid also at \(r=0\) and \(r=q-1\).  Those boundary values are not targets because

\[
b_0=1,\qquad b_{q-1}\equiv1\pmod q.
\]

### Primitive packet lemma

Every vertex of \(P\) occurs in \(\Lambda\) with coefficient \(1\).  A vertex \(u\) has a unique extremal decomposition in \(\Lambda^a\), so

\[
c_a(au)=1.
\tag{2.7}
\]

Thus the integer vector

\[
p_a=(c_a(\mu))_{\mu\in aP\cap\mathbf Z^3}
\tag{2.8}
\]

is primitive.  There is a unimodular matrix \(U_a\in\mathrm{GL}_{N_a}(\mathbf Z)\) with

\[
U_ap_a=e_1.
\tag{2.9}
\]

Applying \(U_a\) to (2.6) gives

\[
U_a(c_n(q\mu))_\mu
\equiv(b_r,0,\ldots,0)^t\pmod q.
\tag{2.10}
\]

Therefore the \(N_a\) packet congruences have one nontrivial Smith invariant:

\[
\boxed{
\mathrm{SNF}_{q}(\text{target packet})
=\operatorname{diag}(1,\ldots,1,q).
}
\tag{2.11}
\]

Equivalently,

\[
q\mid c_n(q\mu)\ \forall\mu
\quad\Longleftrightarrow\quad
q\mid b_r.
\tag{2.12}
\]

This is the first exact saturation theorem: the full packet repeats one target scalar.

---

## 3. Exact logarithmic-derivative finite stencils

Write

\[
\Lambda=\sum_{u\in\mathcal A}\lambda_uX^u,
\qquad
c_n(v)=[X^v]\Lambda^n,
\]

and put \(	heta_i=X_i\partial/\partial X_i\), \(i=1,2,3\).  The identity

\[
\Lambda\theta_i(\Lambda^n)
=n(\theta_i\Lambda)\Lambda^n
\tag{3.1}
\]

gives, after extracting the coefficient of \(X^v\),

\[
\boxed{
\sum_{u\in\mathcal A}
\lambda_u\bigl(v_i-(n+1)u_i\bigr)c_n(v-u)=0,
\qquad i=1,2,3.
}
\tag{3.2}
\]

These are the exact three finite-stencil relations.  Their stencil is the 22-point set \(\mathcal A\subset[-1,1]^3\).

Together with multiplication,

\[
c_{n+1}(v)=\sum_{u\in\mathcal A}\lambda_uc_n(v-u),
\tag{3.3}
\]

they may also be written

\[
v_i c_{n+1}(v)
=(n+1)\sum_{u\in\mathcal A}u_i\lambda_uc_n(v-u).
\tag{3.4}
\]

### GKZ form

Introduce independent support coefficients and define

\[
\mathcal C_{n,v}(\lambda)
=[X^v]\left(\sum_{u\in\mathcal A}\lambda_uX^u\right)^n.
\]

Let \(\widetilde{\mathcal A}\) be the \(4\times22\) matrix whose columns are \((1,u)^t\).  Then

\[
\left(\sum_u\lambda_u\partial_{\lambda_u}-n\right)
\mathcal C_{n,v}=0,
\tag{3.5}
\]

\[
\left(\sum_u u_i\lambda_u\partial_{\lambda_u}-v_i\right)
\mathcal C_{n,v}=0,
\qquad i=1,2,3,
\tag{3.6}
\]

and, for every \(L\in\ker_{\mathbf Z}\widetilde{\mathcal A}\),

\[
\boxed{
\left(
\prod_{L_u>0}\partial_{\lambda_u}^{L_u}
-
\prod_{L_u<0}\partial_{\lambda_u}^{-L_u}
\right)\mathcal C_{n,v}=0.
}
\tag{3.7}
\]

Equations (3.2) and (3.5)--(3.7) are the exact logarithmic-derivative and \(A\)-hypergeometric systems requested in the question.

A useful closed coefficient formula, obtained from the factorization (1.1), is

\[
\boxed{
 c_n(e_1,e_2,e_3)
 =\sum_k
 \binom nk
 \binom n{k-e_1}
 \binom{2n-k}{n-e_2}
 \binom{2n-k}{n-e_3},
}
\tag{3.8}
\]

with the convention that an out-of-range binomial coefficient is zero.

---

## 4. Exact Frobenius defect of the coefficient PDE

Reduce (3.1) modulo \(q\).  Since \(n\equiv r\pmod q\), define

\[
D_iG:=\Lambda\theta_iG-r(\theta_i\Lambda)G.
\tag{4.1}
\]

### Theorem 4.1 — complete polynomial solution space

Among Laurent polynomials \(G\) with

\[
\operatorname{Newt}(G)\subseteq nP,
\]

the simultaneous solution space of

\[
D_iG=0,\qquad i=1,2,3,
\]

is exactly

\[
\boxed{
G=\Lambda^rH(x^q,y^q,z^q),
\qquad
\operatorname{Newt}(H)\subseteq aP.
}
\tag{4.2}
\]

Consequently

\[
\boxed{
\dim_{\mathbf F_q}\mathscr S_{q,n}=N_a.
}
\tag{4.3}
\]

### Proof

Every expression in (4.2) is a solution because

\[
\theta_iH(x^q,y^q,z^q)=0
\]

in characteristic \(q\).

Conversely, in the fraction field,

\[
D_iG=0
\quad\Longleftrightarrow\quad
\theta_i(G/\Lambda^r)=0.
\]

The common kernel of the three torus derivations on
\(\mathbf F_q(x,y,z)\) is

\[
\mathbf F_q(x^q,y^q,z^q).
\]

Hence

\[
G=\Lambda^rR(x^q,y^q,z^q)
\]

for a rational function \(R\).  The nonmonomial factors of \(\Lambda\) are

\[
1+x,\quad1+y,\quad1+z,\quad(1+y)(1+z)+xyz,
\]

all with multiplicity one.  Since \(r<q\), a denominator of \(R(x^q,y^q,z^q)\) would occur with multiplicity divisible by \(q\) and cannot be cancelled by \(\Lambda^r\).  Thus \(R=H\) is Laurent polynomial.

Finally, Newton polytopes add under multiplication.  From

\[
q\operatorname{Newt}(H)+rP
\subseteq nP=(aq+r)P,
\]

support functions give

\[
\operatorname{Newt}(H)\subseteq aP.
\]

There are \(N_a\) possible monomials by (1.5). \(\square\)

### Corollary 4.2 — target packet has zero PDE rank

For \(G=\Lambda^rH(X^q)\), Cartier gives

\[
\mathcal C_q(G)=b_rH.
\tag{4.4}
\]

If \(q\mid b_r\), then

\[
\mathcal C_q(G)=0
\]

for **every** \(H\) in the \(N_a\)-dimensional solution space.  Thus the target packet equations do not reduce the Frobenius defect dimension at all:

\[
\boxed{
\delta_q(a)=N_a.
}
\tag{4.5}
\]

If one fixes one exposed-vertex coefficient of \(H\) as a scale normalization, the remaining defect is \(N_a-1\).

This is stronger than the scalar-shell saturation: it proves that even the coefficient PDE regards the full target packet as a zero-rank condition on its characteristic-\(q\) Frobenius degrees of freedom.

---

## 5. The nearby packet modulo the first prime

Let

\[
\ell=q+h,
\qquad
2ah<q,
\tag{5.1}
\]

as holds uniformly for \(h\le A_0\log n\), \(a\le H=n^{1/3}\), and a mesoscopic block for sufficiently large \(n\).

For

\[
H(X)=\sum_{\kappa\in aP\cap\mathbf Z^3}H_\kappa X^\kappa,
\]

put \(G=\Lambda^rH(X^q)\).  At an \(\ell\)-packet node,

\[
\boxed{
 c_G(\ell\mu)
 =\sum_{\kappa\in aP\cap\mathbf Z^3}
 H_\kappa\,
 c_r\bigl(q(\mu-\kappa)+h\mu\bigr).
}
\tag{5.2}
\]

Thus the \(\ell\)-packet modulo \(q\) is the image of the explicit linear map

\[
E^{(q)}_{a,r,h}:H\longmapsto
\bigl(c_G(\ell\mu)\bigr)_{\mu\in aP\cap\mathbf Z^3}.
\tag{5.3}
\]

The \(q\)-target packet gives no equation on \(H\), so it gives no new equation on (5.3).  Any restriction on the nearby packet comes from the intrinsic rank of \(E^{(q)}_{a,r,h}\), not from targetness at \(q\).

There are two unavoidable exceptions:

1. At \(\mu=0\), the two packet lattices share the same coefficient \(c_n(0)=b_n\), and targetness forces its divisibility.  This is the known common carrier of logarithmic height \(\Theta(n)\).
2. Some exposed-face coefficients may be zero over \(\mathbf Z\) because the \(\ell\)-node lies outside \(nP\).  Those are universal support zeros.

Apart from these, one target packet does not determine the nearby packet through Cartier or the PDE.

### Packet lattices meet only at the origin

If

\[
q\mu=\ell\nu,
\qquad
\mu,\nu\in aP\cap\mathbf Z^3,
\]

then \(q\mid\nu_i\) for all \(i\).  Since \(|\nu_i|\le a<q\), one has \(
u=0\), hence \(\mu=0\).  Therefore

\[
\boxed{
q(aP\cap\mathbf Z^3)
\cap
\ell(aP\cap\mathbf Z^3)=\{0\}.
}
\tag{5.4}
\]

There is no nonzero coefficient coordinate shared by the two packets.

---

## 6. Exposed faces: exact triangularity and residual dimensions

Let \(F\) be a face exposed by an integral functional \(arphi\) normalized so that

\[
\varphi(x)\le1\quad(x\in P),
\qquad
\varphi=1\quad\text{on }F.
\]

For \(\mu\in aF\),

\[
\varphi(\ell\mu)=a\ell=n-r+ah.
\tag{6.1}
\]

Since the support of \(\Lambda^n\) satisfies \(arphi(v)\le n\), two cases occur.

### Outside-face case

If

\[
ah>r,
\tag{6.2}
\]

then

\[
c_n(\ell\mu)=0\quad\text{over }\mathbf Z
\qquad(\mu\in aF).
\tag{6.3}
\]

The whole exposed face of the nearby packet is universally zero.  These zeros are independent of \(q\mid b_r\) and add no selective rank.

### Surviving-face case

If

\[
ah\le r,
\tag{6.4}
\]

the face survives.  The normal-direction coefficient recurrences become triangular, but the free \(H\)-coordinates on the face remain.  Their exact number is

\[
\boxed{
\delta_F(a)=\#(aF\cap\mathbf Z^3).
}
\tag{6.5}
\]

Representative facet counts are:

| facet \(F\) | \(\#(aF\cap\mathbf Z^3)\) |
|---|---:|
| \(x=-1\) | \((2a+1)^2\) |
| \(x=1\) | \((a+1)^2\) |
| \(y=-1\), \(z=-1\) | \((a+1)(2a+1)\) |
| \(y=1\), \(z=1\) | \((7a^2+7a+2)/2\) |
| \(x-y=1\), \(x-z=1\) | \((a+1)(3a+2)/2\) |

Thus exposed-face normalization reduces a cubic defect to a quadratic one, but it does not eliminate it.

### Vertex formula

Let \(u\) be a vertex of \(P\), and set \(\mu=au\).  If \(ah\le r\), then the exposed-normal gap to every other \(\kappa\in aP\) is at least one, so

\[
q\varphi(\mu-\kappa)+ah>r
\]

unless \(\kappa=\mu\).  Formula (5.2) therefore becomes the exact diagonal relation

\[
\boxed{
 c_G(\ell au)=H_{au}\,c_r(ha u).
}
\tag{6.6}
\]

Hence, whenever \(c_r(ha u)\not\equiv0\pmod q\), the nearby vertex coefficient is a completely free Frobenius-defect coordinate.  If \(ah>r\), it is the universal zero from (6.3).

This answers Question 1 precisely: the target packet constrains no surviving nearby face coordinate; face triangularity merely exposes free defect coordinates.

---

## 7. The one positive exposed-vertex carrier

Take the vertex

\[
u=(1,1,1).
\]

The general coefficient formula (3.8) gives, for \(0\le t\le n\),

\[
\boxed{
V_n(t):=c_n(n-t,n-t,n-t)
=\sum_{s=0}^{t}
 \binom ns\binom n{t-s}\binom{n+s}{t}^2.
}
\tag{7.1}
\]

All terms are positive and

\[
\log V_n(t)\le 3t\log(2n)+O(\log(t+1)).
\tag{7.2}
\]

For the \(q\)-packet vertex,

\[
qa u=(n-r)u,
\]

so

\[
c_n(qa u)=V_n(r).
\tag{7.3}
\]

For \(\ell=q+h\), put

\[
t=r-ah.
\]

If \(t\ge0\),

\[
c_n(\ell a u)=V_n(t).
\tag{7.4}
\]

Therefore a close target pair satisfies

\[
q\mid V_n(r),
\qquad
\ell\mid V_n(r-ah),
\]

and the positive integer

\[
\boxed{
D^{\rm vertex}_{q,\ell}=V_n(r)V_n(r-ah)
}
\tag{7.5}
\]

is divisible by \(q\ell\), with

\[
\log D^{\rm vertex}_{q,\ell}\ll r\log n.
\tag{7.6}
\]

This is genuinely sub-\(H\) when

\[
r=o(H/\log n).
\tag{7.7}
\]

However, in a fixed direct quotient block \(r(q)=n-aq\) changes with slope \(-a\).  The number of candidate primes with \(0<r(q)\le R\) is at most

\[
R/a+1.
\]

Taking \(R=o(H/\log n)\) already gives

\[
\#\{q\in I:0<r(q)\le R\}=o(H/\log n).
\tag{7.8}
\]

Thus (7.5) handles only a boundary range which is already negligible by geometry.  It does not enter the central LAC problem.

For a reflected target the small folded index is \(q-1-r\), but every exposed \(q\)-packet face sits at normal distance \(r\) from the corresponding face of \(nP\).  The full packet therefore does not convert a small reflected index into a small-height face coefficient.

---

## 8. Two-characteristic defect and elimination

For a target prime \(q\), the characteristic-\(q\) PDE solution space has defect \(N_a\).  For a target prime \(\ell\), the characteristic-\(\ell\) solution space has another defect \(N_a\).  Over the product ring

\[
\mathbf F_q\times\mathbf F_\ell,
\]

they are independent.  Hence the two-target defect is

\[
\boxed{
\delta_{q,\ell}(a,h)=2N_a
=\frac{38a^3+57a^2+31a+6}{3}.
}
\tag{8.1}
\]

After fixing one exposed-vertex coefficient in each Frobenius factor, it is

\[
\boxed{2N_a-2.}
\tag{8.2}
\]

These dimensions are independent of \(h\), provided \(2ah<\min(q,\ell)\); \(h\) changes only which corridor coefficient functionals are evaluated.

### Smith count for the actual packets

For the distinguished array \(\Lambda^n\), each packet is a primitive scalar vector.  Thus the local Smith invariants of the two packet conditions are

\[
\boxed{
\operatorname{diag}
(1^{N_a-1},q)
\oplus
\operatorname{diag}
(1^{N_a-1},\ell).
}
\tag{8.3}
\]

The many packet rows contribute one factor \(q\) and one factor \(\ell\), not \(q^{N_a}\ell^{N_a}\).  If a determinant presentation displays the latter powers before saturation, \(N_a-1\) powers in each characteristic belong to the fake Frobenius-defect submodule.  Primitive Fitting saturation removes them.

### Scalar elimination ideal

Consider the universal coefficient system consisting of:

- the finite-stencil equations (3.2),
- support in \(nP\), and
- the vanishing \(q\)-packet.

The family (4.2), with arbitrary \(H\), satisfies all these equations.  Therefore eliminating all coefficient and \(H\)-variables produces no nonzero scalar condition.  Adding the \(\ell\)-system gives the product of two independent such families, so the scalar elimination ideal remains zero.

This proves the scoped no-go:

> **Cartier--PDE saturation theorem.**  No nonzero characteristic-zero determinant or resultant is a formal consequence of the two target packets and the logarithmic-derivative/GKZ coefficient equations alone.  A nonzero determinant must specialize the Frobenius defects to the distinguished coordinates \(H=\Lambda^a\) in both characteristics, or use another Apéry-specific identity not present in the packet/PDE module.

The theorem does not exclude a new nonlinear identity special to the actual coefficient array.  It proves that the proposed rank amplification cannot come from the number of packet zeros or from exposed-face triangularity.

---

## 9. Why a finite corridor does not close the defect

For \(\mu\in aP\), the two corresponding nodes differ by

\[
\ell\mu-q\mu=h\mu,
\qquad
\|h\mu\|_\infty\le ah.
\tag{9.1}
\]

Thus one can surround the two packet lattices by a corridor of thickness \(O(ah)\) and impose all three stencils (3.2) there.  This does not alter Theorem 4.1: every global solution

\[
\Lambda^rH(X^q)
\]

restricts to a corridor solution, and the \(q\)-packet boundary is zero for every \(H\) on the target locus.  The corridor therefore retains the \(N_a\) Frobenius degrees of freedom before any distinguished-state normalization.

On an exposed face the normal recurrence is triangular, but Section 6 shows exactly what remains: one free coordinate for every lattice point of \(aF\).  In the interior the defect is cubic in \(a\); on facets it is quadratic; on edges it is linear; at a vertex it is one.  No face has negative defect, and an identically zero outside face is nonselective.

Consequently:

- root-of-unity filters are linear combinations of packet rows and retain the Smith content (2.11);
- finite differences in the displacement \(h\mu\) operate on the free map (5.3), not on a determined packet;
- determinants of corridor stencils either vanish identically on the Frobenius family or retain free \(H\)-coordinates;
- resultants in two residue characteristics split over \(\mathbf F_q\times\mathbf F_\ell\) and cannot create a common-field root;
- exposed-face normalization reduces the defect dimension but supplies only the boundary carrier (7.5).

Thus there is no general \(o(H)\)-height \(q\ell\)-carrier from the full packet and its coefficient PDEs.

---

## 10. Height and conductor warning

The distinguished solution \(H=\Lambda^a\) has positive coefficients, and

\[
\sum_v c_n(v)=\Lambda(1,1,1)^n=40^n.
\]

The common origin coefficient is

\[
c_n(0)=b_n,
\qquad
\log b_n=\Theta(n).
\]

Raw packet coefficients and the evident cross-weighted determinants therefore remain on a linear exponential scale except in the exposed boundary range of Section 7.

Likewise, a geometric trace estimate cannot simply treat \(\Lambda^n\) as bounded conductor.  The fixed Laurent map \(\Lambda\) has bounded support, but the coefficient/PDE system at moment \(n\) has Frobenius defect dimension

\[
N_a\asymp a^3,
\]

and the relevant exponent \(n\) grows.  Reducing the exponent modulo \(q-1\) recovers the moving residue \(r\), while retaining the fixed exponent keeps a growing Kummer/power sheaf.  A Katz--Deligne citation would require an explicit compatible system and conductor bound uniform in \(n,a,q\); the packet identity alone provides neither.

---

## 11. The smallest remaining lemma

Because the full packet condition is equivalent to one scalar condition by (2.12), packet rank does not weaken the horizontal problem.  The smallest sufficient statement remains the localized adjacent collision estimate, now written intrinsically in Cartier language.

Let \(H=n^{1/3}\).  In a length-\(H\) prime block \(I\), fix the quotient \(a\le H\) and one fold branch.  For \(q\in I\), write

\[
n=aq+r_q.
\]

### Cartier scalar anti-clustering (CSAC)

For every fixed \(A_0>0\), uniformly in \(n,a,I\),

\[
\boxed{
\sum_{1\le h\le A_0\log n}
\#\left\{
\begin{array}{l}
q,q+h\in I\text{ prime}:\\
\mathcal C_q(\Lambda^n)=0,\\
\mathcal C_{q+h}(\Lambda^n)=0
\end{array}
\right\}
=o_{A_0}(H/\log n).
}
\tag{11.1}
\]

Here equality to zero is equality of the full Cartier polynomial packet.  By (2.5), it is equivalent to

\[
q\mid b_{r_q},
\qquad
q+h\mid b_{r_{q+h}}.
\]

Thus (11.1) is exactly LAC, not a weaker consequence.  Standard adjacent-gap counting then gives

\[
|T(n,I)|=o(H/\log n),
\]

and the mesoscopic block decomposition gives P3.2.

A stronger constructive lemma, sufficient for (11.1), would be:

> **Distinguished Cartier-corridor carrier.**  For every fixed \(A_0\), and every \(n,a,I\), there are nonzero integers \(D_{n,a,I,h}\), \(1\le h\le A_0\log n\), such that every close pure-cross target pair \(q,q+h\in I\) satisfies
> \[
> q(q+h)\mid D_{n,a,I,h},
> \]
> and
> \[
> \sum_{h\le A_0\log n}\log|D_{n,a,I,h}|=o_{A_0}(H).
> \]

The present packet/PDE module cannot prove this lemma: its scalar elimination ideal is zero and its residual Frobenius defect is (8.1).  A proof must use a new identity selecting the actual characteristic-zero defect coordinate \(H=\Lambda^a\), or a genuine cross-characteristic arithmetic theorem for the Cartier scalar \(b_r\).

---

## 12. Final answers to the four questions

### 1. Does one target packet constrain an \(\ell\)-packet modulo \(q\)?

Only in two nonselective ways: the shared origin \(b_n\), and coefficients forced to be zero because an exposed \(\ell\)-face lies outside \(nP\).  On every surviving face the packet leaves the face-sized Frobenius defect (6.5); at a vertex it leaves the free coordinate (6.6) whenever its diagonal coefficient is a \(q\)-unit.  No additional target constraint follows from the packet or coefficient PDEs.

### 2. Can a finite corridor give a nonzero \(q\ell\)-determinant of height \(o(H)\)?

Not in the central range by Cartier/PDE elimination.  The finite corridor retains an \(N_a\)-dimensional defect in each characteristic.  A positive carrier exists only in the direct exposed-vertex boundary range \(r=o(H/\log n)\), where (7.5) has height \(o(H)\); that range already contains \(o(H/\log n)\) candidates.

### 3. What is the remaining defect dimension?

For one characteristic it is

\[
N_a=\frac{38a^3+57a^2+31a+6}{6},
\]

or \(N_a-1\) after scale normalization.  For two characteristics it is \(2N_a\), or \(2N_a-2\) normalized.  An exposed face \(F\) retains exactly \(\#(aF\cap\mathbf Z^3)\) defect coordinates.  These dimensions are independent of \(h\) in the mesoscopic close-gap range \(2ah<\min(q,\ell)\).

### 4. Are there useful faces or quotients?

Yes only at the direct small-remainder boundary.  If \(r=o(H/\log n)\), the corner packet coefficient gives the \(o(H)\)-height carrier (7.5), and the entire range contributes \(o(H/\log n)\) candidates.  If \(ah>r\), the nearby exposed face is universally outside the support and gives no target information.  No corresponding reflected small-index range is obtained from this Cartier packet.

The full packet therefore confirms, rather than escapes, the scalar-shell obstruction: its many zeros form one primitive Cartier scalar multiplied by a large Frobenius defect module.