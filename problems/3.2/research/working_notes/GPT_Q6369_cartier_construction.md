ANSWER Q6369 1eb3f40f

# Executive verdict

There are three distinct statements here, and they should not be conflated.

1. **The obvious fixed curves do not work.** The double cover
   \[
   y^2=1-34t+t^2
   \]
   has genus zero, hence no nonzero holomorphic differentials and no Hasse--Witt matrix. Its Franel pullback splits because
   \[
   q(\phi(u))=\left(\frac{1-16u-8u^2}{1+u}\right)^2.
   \]
   The natural quartic reversion cover
   \[
   C:\quad w^2=s^8-34s^4+1
   \tag{0.1}
   \]
   has genus three and has exactly the right quarter-index geometry, but its **ordinary** Hasse--Witt matrix contains coefficients of
   \((s^8-34s^4+1)^{(p-1)/2}\), not coefficients of \(\tau\) or \(\sigma\). It already gives the wrong answer at \(p=5\).

2. **There is an explicit fixed geometric period family.** Start from the Hesse elliptic family for
   \({}_2F_1(1/3,2/3;1;z)\), pull it back by
   \[
   z(u)=\frac{27u^2}{(1-2u)^3},
   \]
   twist quadratically by \(1+u\), and then pull back to the fixed genus-three curve (0.1) using the local branch
   \[
   u=\frac{1-s^4-w}{16},\qquad \phi(u)=s^4,
   \tag{0.2}
   \]
   at \((s,w)=(0,1)\). The resulting elliptic family has normalized local period \(\tau(s^4)\); \(\sigma(s^4)\) is obtained from the same family by a rational change of the Hodge section.

3. **The quarter coefficient is an exact entry of a fixed-dimensional Cartier module with coefficients in that rank-two period system.** It is not, in the construction presently justified, an entry of the ordinary Hasse--Witt matrix of one fixed curve. The relevant module has dimension at most six: a two-dimensional branch space tensored with the three-dimensional principal-parts space
   \(H^1(\mathbf P^1,\mathcal O(-4))\). The desired entries are
   \[
   \tau_{(p-1)/4}\quad\text{at matrix position }(1,1),
   \]
   and
   \[
   \sigma_{(p-3)/4}\quad\text{at matrix position }(3,1).
   \]
   This construction includes the nonlinear reversion and the branch switch when \((-6|p)=-1\).

Thus the strongest proved statement is:

> **Proved local Cartier realization.** The selected quarter coefficient is a matrix entry of a bounded-rank, explicitly defined **local Cartier operator with coefficients in the Apéry rank-two period module** on the fixed quartic reversion cover.

What remains open is the global motivic upgrade:

> **[GAP-1]** Identify this local Cartier operator with a canonically extended global crystalline/overconvergent Frobenius operator of bounded conductor, with no choice of local trivialization hidden in the definition.

The new factorization
\[
\tau_{(p-1)/4}=2A\,U_p\pmod p
\]
is naturally explained by the \(\mu_4\)-Kummer projector introduced by \(t=s^4\). The factor \(2A\) is the oriented quartic Hasse--Witt factor of \(y^2=x^3-x\). It does **not** yet prove that the entire object is a tensor product with that elliptic curve, and it does not point specifically to a “degree-24 CM abelian variety.” The correct possible residual object is an order-dividing-24 Kummer/CM **representation or Hecke-character trace**, not a rank-one finite-order character.

No theorem below proves that the quarter values are non-motivic. In fact, the exact residue-class zero law and the factor \(2A\) are compatible with bounded-conductor geometry. What can be proved is a no-go theorem for the obvious curves and for a finite-order rank-one twist.

---

# 1. An explicit elliptic period family

## 1.1 The Hesse rank-two source

Put
\[
E_z:\quad y^2+xy+\frac{1-z}{27}y=x^3.
\tag{1.1}
\]
For this Weierstrass equation,
\[
a_1=1,\qquad a_3=\frac{1-z}{27},\qquad a_2=a_4=a_6=0.
\]
A direct invariant calculation gives
\[
\Delta(E_z)=\frac{z(1-z)^3}{27^3},
\qquad
j(E_z)=\frac{27(1+8z)^3}{z(1-z)^3}.
\tag{1.2}
\]
Indeed, writing \(a=(1-z)/27\), one has
\[
b_2=1,\quad b_4=a,\quad b_6=a^2,\quad b_8=0,
\]
so
\[
\Delta=-8a^3-27a^4+9a^3=a^3(1-27a)=a^3z,
\]
and
\[
c_4=1-24a=\frac{1+8z}{9}.
\]
This is a rational model of the Hesse pencil. Its normalized period at the cusp \(z=0\) is
\[
g(z)={}_2F_1\!\left(\begin{matrix}1/3,2/3\\1\end{matrix};z\right),
\tag{1.3}
\]
the unique power-series solution with constant term one of
\[
z(1-z)g''+(1-2z)g'-\frac29g=0.
\tag{1.4}
\]
The Picard--Fuchs identification of (1.1) with (1.4) is standard for the Hesse pencil; it can also be verified directly by Griffiths reduction. No rank or Frobenius conclusion below depends on more than the differential equation and the normalization at \(z=0\).

## 1.2 Franel pullback and the two Apéry square-root branches

Set
\[
z(u)=\frac{27u^2}{(1-2u)^3},
\qquad
\phi(u)=\frac{u(1-8u)}{1+u}.
\tag{1.5}
\]
The Franel period is
\[
h(u)=\frac1{1-2u}g(z(u)).
\tag{1.6}
\]
The supplied characteristic-zero identities are
\[
F(\phi(u))=(1+u)h(u)^2,
\tag{1.7}
\]
and
\[
q(\phi(u))
=\left(\frac{1-16u-8u^2}{1+u}\right)^2.
\tag{1.8}
\]
Therefore
\[
\boxed{
\tau(\phi(u))
=\frac{\sqrt{1+u}}{1-2u}\,g(z(u)),
}
\tag{1.9}
\]
and
\[
\boxed{
\sigma(\phi(u))
=\frac{(1+u)^{3/2}}
 {(1-16u-8u^2)(1-2u)}\,g(z(u)).
}
\tag{1.10}
\]

These formulas give an explicit geometric realization. Complete the square in (1.1): with
\[
Y=2y+x+\frac{1-z}{27},
\]
one obtains
\[
Y^2=4x^3+x^2+rac{2(1-z)}{27}x+rac{(1-z)^2}{27^2}.
\tag{1.11}
\]
Pull back by \(z=z(u)\), then take the quadratic twist
\[
(1+u)Y^2
=4x^3+x^2+rac{2(1-z(u))}{27}x+rac{(1-z(u))^2}{27^2}.
\tag{1.12}
\]
Over \(\mathbf Q(u)(\sqrt{1+u})\), the differential on the twist acquires the factor \(\sqrt{1+u}\). Multiplying the chosen Hodge section additionally by \((1-2u)^{-1}\) gives the period (1.9). The sigma section is the rationally gauged section
\[
\omega_\sigma
=\frac{1+u}{1-16u-8u^2}\,\omega_\tau.
\tag{1.13}
\]
Thus tau and sigma do not require unrelated curves; they are two local Hodge sections of the same rank-two elliptic Gauss--Manin system, with one quadratic twist and one rational gauge.

## 1.3 The nonlinear reversion as a fixed genus-three cover

The equation \(t=\phi(u)\) is equivalent to
\[
8u^2+(t-1)u+t=0.
\tag{1.14}
\]
Its discriminant is
\[
(t-1)^2-32t=t^2-34t+1=q(t).
\tag{1.15}
\]
Now put \(t=s^4\) and introduce
\[
C:\quad w^2=q(s^4)=s^8-34s^4+1.
\tag{1.16}
\]
The polynomial on the right has degree eight and is squarefree away from the fixed bad characteristics \(2,3\), hence \(C\) has genus three. The two roots of (1.14) are
\[
u_\pm=\frac{1-s^4\pm w}{16}.
\tag{1.17}
\]
At the point
\[
P_0=(s,w)=(0,1),
\]
the correct local branch is
\[
\boxed{
u=rac{1-s^4-w}{16},}
\tag{1.18}
\]
because \(u(P_0)=0\). It satisfies \(\phi(u)=s^4\) identically in \(\mathbf Q(C)\).

Pulling (1.12) and its tau differential to \(C\) via (1.18) produces an explicit elliptic curve over the fixed function field \(\mathbf Q(C)\) whose normalized local period at \(P_0\) is
\[
\tau(s^4)=\sum_{n\ge0}\tau_n s^{4n}.
\tag{1.19}
\]
The sigma section has period \(\sigma(s^4)\).

This construction completely incorporates the nonlinear reversion. It is independent of \(p\).

---

# 2. The two relevant involutions

There are two different involutions, and confusing them obscures the geometry.

## 2.1 The deck involution of the quadratic reversion

The second root of (1.14), expressed in terms of the first, is
\[
\iota(u)=\frac{1-8u}{8(1+u)}.
\tag{2.1}
\]
It satisfies
\[
\phi(\iota(u))=\phi(u).
\]
It exchanges the local branch at \(u=0\) with the branch at \(u=1/8\). It does not identify two copies of the same Taylor germ.

## 2.2 The lift of coefficient reflection \(t\mapsto1/t\)

The Möbius involution
\[
J(u)=-\frac1{8u}
\tag{2.2}
\]
satisfies
\[
\boxed{\phi(J(u))=\frac1{\phi(u)}.}
\tag{2.3}
\]
The quotient of \(\mathbf P^1_u\) by \(J\) has genus zero. Its fixed points satisfy
\[
u^2=-\frac18,
\tag{2.4}
\]
which is the elementary geometric source of the quadratic character \((-2|p)\) in the branch reversal law.

On the genus-three curve (1.16), \(t\mapsto1/t\) has the two lifts
\[
R_\pm(s,w)=\left(\frac1s,\ \pm\frac{w}{s^4}\right).
\tag{2.5}
\]
For \(R_+\), the fixed points have \(s=\pm1\), with two choices of \(w\) at each value; there are four fixed points. Riemann--Hurwitz gives
\[
g(C/R_+)=1.
\tag{2.6}
\]
The lift \(R_-\) is fixed-point-free, and
\[
g(C/R_-)=2.
\tag{2.7}
\]

A basis of holomorphic differentials on \(C\) is
\[
\omega_j=s^{j-1}\frac{ds}{w},
\qquad j=1,2,3.
\tag{2.8}
\]
Direct substitution gives
\[
R_+^*(\omega_j)=-\omega_{4-j}.
\tag{2.9}
\]
Hence the invariant line for \(R_+\) is spanned by
\[
\omega_1-\omega_3,
\]
while the invariant space for \(R_-\) is the two-dimensional span of
\[
\omega_1+\omega_3,\qquad\omega_2.
\]
These are the genus-one and genus-two quotient differential spaces.

The important negative conclusion is that passing to either quotient only extracts an eigenspace of the **same ordinary Cartier operator of \(C\)**. Since that ordinary operator has the wrong multiplier, quotienting does not repair the problem.

---

# 3. Ordinary Hasse--Witt matrices of the obvious curves

## Lemma 3.1: hyperelliptic Cartier formula

Let
\[
X:\quad y^2=f(x)
\]
be a hyperelliptic curve of genus \(g\) over \(\mathbf F_p\), \(p>2\), with \(f\) squarefree, and write
\[
f(x)^{(p-1)/2}=\sum_{m\ge0}c_mx^m.
\tag{3.1}
\]
In the basis
\[
x^{j-1}\frac{dx}{y},\qquad1\le j\le g,
\]
the Cartier--Manin matrix has entries \(c_{pi-j}\), up to the standard transpose convention between Cartier--Manin and Hasse--Witt matrices. Over \(\mathbf F_p\), the possible \(p\)-th roots of the coefficients do not alter the entries.

This formula follows immediately by writing
\[
\frac{x^{j-1}dx}{y}
=f(x)^{(p-1)/2}\frac{x^{j-1}dx}{y^p}
\]
and retaining the monomials whose exponent is congruent to \(p-1\) modulo \(p\). Only the convention of whether vectors act on the left or right changes the transpose.

## 3.2 Candidate (a): \(y^2=q(t)\)

The curve
\[
y^2=t^2-34t+1
\]
has genus zero. Therefore
\[
H^0(C,\Omega^1)=0,
\]
and there is no nontrivial Cartier or Hasse--Witt matrix to contain a quarter coefficient.

Moreover, after the Franel pullback it splits:
\[
y^2=q(\phi(u))
=\left(\frac{1-16u-8u^2}{1+u}\right)^2.
\]
Thus the double cover supplies the two algebraic branches, but no positive-genus Jacobian.

## 3.3 Candidate (b): quartic base change and reversion curve

The curve (1.16) is the correct fixed genus-three reversion curve. Its ordinary Hasse--Witt matrix is nevertheless wrong for the desired purpose.

Write
\[
(s^8-34s^4+1)^{(p-1)/2}
=\sum c_ms^m.
\]
Its Hasse--Witt matrix is \((c_{pi-j})_{1\le i,j\le3}\). Because the defining polynomial is a polynomial in \(s^4\), the matrix is sparse:

- for \(p\equiv1\pmod4\), it is diagonal;
- for \(p\equiv3\pmod4\), it is anti-diagonal together with the middle entry.

This is precisely the correct **index pattern**, but the entries are quarter coefficients of
\[
q(t)^{(p-1)/2},
\]
not of \(\tau\) or \(\sigma\).

There is an exact one-prime disproof. At \(p=5\),
\[
q(s^4)=s^8+s^4+1,
\]
and
\[
q(s^4)^2=s^{16}+2s^{12}+3s^8+2s^4+1.
\]
Therefore the genus-three Hasse--Witt matrix is, up to transpose,
\[
\begin{pmatrix}
2&0&0\\
0&3&0\\
0&0&2
\end{pmatrix}.
\tag{3.2}
\]
But
\[
\tau_{(5-1)/4}=\tau_1=\frac52\equiv0\pmod5.
\tag{3.3}
\]
So the ordinary Cartier operator of \(C\) cannot be the desired operator.

The genus-one and genus-two quotients from Section 2 inherit subblocks of (3.2), and hence also fail at \(p=5\).

This proves the strongest clean obstruction for the obvious candidates:

> **No-go theorem for the base curves.** Neither \(y^2=q(t)\), the quartic reversion curve \(w^2=q(s^4)\), nor either reflection quotient has ordinary Hasse--Witt matrix entry equal to the selected Apéry quarter coefficient for all good primes.

---

# 4. Exact branch Frobenius identities

The supplied Dwork relation is
\[
F(t)=A_p(t)F(t^p)\quad\text{in }\mathbf F_p[[t]].
\tag{4.1}
\]
The two factorization cases now give exact identities for the square-root branches.

## Lemma 4.1: split branch

If \(\chi_p=(-6|p)=+1\), then
\[
A_p=S_p^2,
\]
and
\[
\boxed{
\tau(t)=S_p(t)\tau(t^p).
}
\tag{4.2}
\]

### Proof

Square the right side and use (4.1):
\[
S_p(t)^2\tau(t^p)^2=A_p(t)F(t^p)=F(t)=\tau(t)^2.
\]
Both sides have constant term one, so their formal square roots agree. ∎

## Lemma 4.2: nonsplit branch and branch switching

If \(\chi_p=-1\), then
\[
A_p=qS_p^2,
\]
and
\[
\boxed{
\sigma(t)=S_p(t)\tau(t^p).
}
\tag{4.3}
\]

### Proof

Using \(\sigma^2=F/q\),
\[
\sigma(t)^2
=\frac{F(t)}{q(t)}
=S_p(t)^2F(t^p)
=S_p(t)^2\tau(t^p)^2.
\]
Again both formal square roots have constant term one. ∎

This is the correct branch bookkeeping. In the split case Frobenius preserves the tau line. In the nonsplit case the selected scalar \(S_p\) is an **off-diagonal entry** from the Frobenius pullback of the tau line to the sigma line. A one-dimensional scalar model misses this distinction.

After substituting \(t=s^4\), put
\[
H_p(s)=S_p(s^4)=\sum_{m\ge0}h_ms^m.
\tag{4.4}
\]
Then
\[
h_{4j}=s_j,
\qquad h_m=0\quad(4\nmid m).
\tag{4.5}
\]

---

# 5. The exact bounded-dimensional Cartier package

The following construction is elementary and complete. It is the surviving positive result.

## Lemma 5.1: coefficient extraction by a Cartier map

Let \(k=\mathbf F_p\), let
\[
H(s)=\sum_{m=0}^{4(p-1)}h_ms^m,
\]
viewed as a section of \(\mathcal O_{\mathbf P^1}(4(p-1))\). Multiplication by \(H\) gives a morphism
\[
F^*\mathcal O(-4)=\mathcal O(-4p)
\xrightarrow{\ \cdot H\ }
\mathcal O(-4).
\tag{5.1}
\]
It induces a Frobenius-semilinear operator
\[
\mathcal C_H:
H^1(\mathbf P^1,\mathcal O(-4))
\longrightarrow
H^1(\mathbf P^1,\mathcal O(-4)).
\tag{5.2}
\]
Use the Čech basis
\[
e_j=[s^{-j}],\qquad j=1,2,3.
\tag{5.3}
\]
Then the matrix entry in row \(i\), column \(j\) is
\[
\boxed{(\mathcal C_H)_{ij}=h_{pj-i}.}
\tag{5.4}
\]

### Proof

Frobenius sends
\[
e_j\longmapsto[s^{-pj}].
\]
Multiplication by \(H\) gives
\[
\sum_mh_m[s^{m-pj}].
\]
In \(H^1(\mathbf P^1,\mathcal O(-4))\), the only surviving Laurent classes are \(s^{-1},s^{-2},s^{-3}\). The coefficient of \(s^{-i}\) occurs precisely when
\[
m-pj=-i,
\]
i.e. \(m=pj-i\). ∎

## Corollary 5.2: the two quarter entries

Take \(H=H_p=S_p(s^4)\).

If \(p\equiv1\pmod4\), then
\[
(\mathcal C_{H_p})_{11}
=h_{p-1}
=s_{(p-1)/4}.
\tag{5.5}
\]
If \(p\equiv3\pmod4\), then
\[
(\mathcal C_{H_p})_{31}
=h_{p-3}
=s_{(p-3)/4}.
\tag{5.6}
\]
Thus, in the relevant branches,
\[
\boxed{
\tau_{(p-1)/4}=(\mathcal C_{H_p})_{11}
\quad(\chi_p=+1,\ p\equiv1\bmod4),
}
\tag{5.7}
\]
and
\[
\boxed{
\sigma_{(p-3)/4}=(\mathcal C_{H_p})_{31}
\quad(\chi_p=-1,\ p\equiv3\bmod4).
}
\tag{5.8}
\]

## 5.3 Incorporating the branch space

Let \(B\) be the two-dimensional formal branch space with basis
\[
e_\tau,e_\sigma.
\]
Let
\[
P=H^1(\mathbf P^1,\mathcal O(-4)),
\qquad\dim P=3.
\]
On the fixed six-dimensional space \(B\otimes P\), define the selected part of Frobenius by
\[
e_\tau\otimes v
\longmapsto
\begin{cases}
e_\tau\otimes\mathcal C_{H_p}(v),&\chi_p=+1,\\
e_\sigma\otimes\mathcal C_{H_p}(v),&\chi_p=-1.
\end{cases}
\tag{5.9}
\]
Equations (4.2)--(4.3) show that this is exactly the local branch transition induced by the Apéry period system.

The desired quarter values are therefore selected entries of one fixed-dimensional semilinear matrix, including the nonlinear reversion and the branch switch.

This proves a concrete version of the requested Cartier construction.

## 5.4 Relation with the genus-three curve

For the double cover
\[
\pi:C\to\mathbf P^1_s,
\]
branched at the eight roots of \(q(s^4)\), one has
\[
\pi_*\mathcal O_C
\simeq\mathcal O\oplus\mathcal O(-4).
\tag{5.10}
\]
Since \(H^1(\mathbf P^1,\mathcal O)=0\),
\[
H^1(C,\mathcal O_C)
\simeq H^1(\mathbf P^1,\mathcal O(-4)).
\tag{5.11}
\]
Thus the three-dimensional space in Lemma 5.1 is not artificial: it is the anti-invariant coherent cohomology of the explicit genus-three reversion curve.

However, the operator \(\mathcal C_{H_p}\) is **Cartier with coefficients in the Apéry branch Frobenius scalar**. Ordinary Cartier on \(C\) replaces \(H_p\) by
\[
q(s^4)^{(p-1)/2},
\]
which is why Section 3 gave the wrong matrix.

This is the precise answer to Question (1):

- no ordinary curve Hasse--Witt matrix among the obvious candidates works;
- the desired coefficient is an entry of a twisted/coefficient Cartier operator on the fixed genus-three reversion curve;
- the coefficient system is the rank-two elliptic period system constructed in Section 1.

---

# 6. What remains to make the construction globally motivic

The local construction above is exact, but “bounded conductor” normally means more than bounded matrix dimension. It requires a fixed global geometric or overconvergent object with a fixed singular divisor.

The explicit candidate is:

1. the genus-three base curve \(C\) of (1.16);
2. the pulled-back Hesse elliptic family (1.12) over \(\mathbf Q(C)\);
3. its rank-two relative crystalline cohomology;
4. the tau/sigma Hodge sections from (1.9)--(1.13);
5. the \(H^1(C,\mathcal O_C)\) anti-invariant piece from (5.11).

The formal Dwork identities prove that the local Frobenius entry at \(P_0\) is \(S_p(s^4)\). To turn the six-dimensional package (5.9) into a canonical global Frobenius matrix, two points still require proof.

**[GAP-1: geometric Frobenius identification].** Verify directly, for the explicit Hesse pullback and twist (1.12), that the formal Frobenius giving (4.2)--(4.3) is the crystalline Frobenius on the indicated relative de Rham basis, including every gauge factor in (1.9)--(1.13).

This is expected from the Picard--Fuchs/Dwork construction, but the equality of a normalized formal period ratio with a specific crystalline matrix entry is basis-sensitive and should not be asserted without the calculation.

**[GAP-2: global extension].** Extend the branch lines and their Frobenius across the fixed singular divisor on \(C\), and identify the induced map on a global cohomology group with the local matrix (5.9). Equivalently, construct the required coefficient F-crystal or overconvergent isocrystal globally and prove that its conductor is bounded independently of \(p\).

The singular set is visibly bounded: it is contained in the preimages of
\[
z=0,1,\infty,\qquad
u=-1,\qquad
1-2u=0,\qquad
1-16u-8u^2=0,
\]
and the branch points of \(C\). Thus no growing singular support is present. What is missing is the canonical extension and comparison theorem, not a plausible finite set of bad points.

Subject to [GAP-1] and [GAP-2], the quarter coefficient becomes an honest bounded-conductor Frobenius matrix entry in the cohomology of a fixed curve with coefficients in a fixed rank-two geometric crystal. No growing-genus construction is needed.

---

# 7. The factor \(2A\) and the quartic projector

## 7.1 Direct Hasse--Witt calculation for \(y^2=x^3-x\)

Let
\[
E_i:\quad y^2=x^3-x=x(x^2-1).
\]
For \(p=4m+1\), the Hasse invariant is the coefficient of \(x^{p-1}\) in
\[
(x^3-x)^{(p-1)/2}
=x^{2m}(x^2-1)^{2m}.
\]
The required term uses \(x^{2m}\) from \((x^2-1)^{2m}\), hence
\[
\operatorname{Ha}_p(E_i)
=(-1)^m\binom{2m}{m}.
\tag{7.1}
\]
After the standard Gauss orientation of the invariant differential, this is the normalized number denoted \(2A\) in the question:
\[
2A\equiv\binom{(p-1)/2}{(p-1)/4}\pmod p,
\tag{7.2}
\]
with the predictable sign absorbed in the choice of \(A\).

Thus \(2A\) is genuinely an oriented quartic Cartier datum.

## 7.2 Exact Mellin interpretation of the quarter coefficient

If \(S_p(t)=\sum_{n=0}^{D}s_nt^n\) and \(D<p-1\), then for \(0\le n<p-1\),
\[
\boxed{
 s_n=-\sum_{x\in\mathbf F_p^\times}S_p(x)x^{-n}.
}
\tag{7.3}
\]
Indeed,
\[
\sum_{x\in\mathbf F_p^\times}x^{m-n}
=
\begin{cases}-1,&m=n,\\0,&m\ne n,
\end{cases}
\]
for \(0\le m,n<p-1\).

When \(p\equiv1\pmod4\) and \(n=(p-1)/4\), the weight
\[
x^{-n}
\]
is a quartic character. Therefore the tau-quarter coefficient is exactly a quartic Mellin component of the Hasse/Frobenius scalar \(S_p\).

This supplies the conceptual source of the factor \(2A\): the quartic base change \(t=s^4\) and the corresponding \(\mu_4\)-eigenspace insert the same Kummer projector whose oriented Gauss/Jacobi factor is the Hasse invariant of \(E_i\).

What is **not** proved by (7.3) is the complete factorization
\[
s_{(p-1)/4}=2A\,U_p.
\]
To prove that identity geometrically, one must factor the quartic Mellin transform of the Apéry Frobenius scalar into the standard quartic Gauss factor and a residual bounded-conductor trace.

**[GAP-3: residual factor].** Construct the residual Frobenius object giving \(U_p\), or derive a Jacobi-sum identity that extracts \(2A\) from (7.3) and leaves a fixed compatible trace.

## 7.3 Why the sigma class is different

For \(p\equiv23\pmod{24}\), write
\[
p=24k+23,
\qquad n=\frac{p-3}{4}=6k+5.
\]
Then
\[
\gcd(n,p-1)=\gcd(6k+5,24k+22)=1,
\tag{7.4}
\]
because \((p-1)-4n=2\) and \(n\) is odd. Hence the multiplicative character \(x\mapsto x^n\) has full order \(p-1\), not order \(4\), \(8\), \(12\), or \(24\).

Therefore a single order-24 multiplicative character on \(\mathbf F_p^\times\) cannot explain the sigma-quarter coefficient directly. Its fixed mod-24 vanishing law comes instead from the branch-switching Frobenius and reversal parity. Any order-24 interpretation of the sigma class must use descent from an extension field or a different cyclotomic action, not the naive Mellin exponent.

---

# 8. What bounded object could give \(2A\,U_p\)?

The most natural candidate is not “a degree-24 CM abelian variety.” There are two arithmetic corrections.

1. \(\mathbf Q(\zeta_{24})\) has degree
   \[
   \varphi(24)=8.
   \]
   A simple abelian variety with full CM by this field would have dimension \(4\), not \(24\).

2. The maximal real subfield \(\mathbf Q(\zeta_{24})^+\) has degree \(4\) and is totally real, not a CM field. It does not have elliptic-curve \(j\)-invariants in the usual CM sense.

A correct bounded candidate would be one of the following.

- A rank-one algebraic Hecke character over \(\mathbf Q(\zeta_{24})\), induced to \(\mathbf Q\), giving a representation of rank at most eight.
- A trace of several Galois conjugates of such a Hecke character.
- A tensor constituent of the Apéry rank-two crystal with the quartic Kummer sheaf and the cubic/level-six Hesse structure.
- A motive inside the cohomology of a product or Kummer surface built from the \(j=1728\) elliptic curve and a CM factor of conductor supported at \(2\) and \(3\).

The third option is the one directly suggested by Sections 1 and 7. The relevant orders arise from
\[
\mu_4\quad\text{(quartic base projector)},
\qquad
\mu_3\text{ and level }6\quad\text{(Hesse/Franel source)}.
\]
The least common cyclotomic modulus is naturally supported on \(24\).

But a rank-one finite-order twist is impossible if \(U_p\) vanishes in the quarter-zero class. At every unramified prime, a finite-order character takes a root-of-unity value and is never zero. Thus:

> **Finite-order-twist obstruction.** The residual \(U_p\) cannot be the value of one fixed finite-order rank-one character at \(p\) if \(U_p=0\) for infinitely many unramified primes.

A trace of an induced CM representation can vanish, so this does not rule out a CM trace or a higher-rank Artin/Hecke object.

Most importantly, a congruence factorization modulo \(p\)
\[
\tau_{(p-1)/4}=2A\,U_p
\]
is basis-dependent and does not by itself imply a tensor decomposition of compatible systems. One needs compatible algebraic lifts and equality of Frobenius polynomials, not only one matrix entry modulo \(p\).

---

# 9. Why a global impossibility theorem is not currently available

The proposed negative criteria do not yet prove non-motivicity.

## 9.1 Interpolation degree is not an obstruction

A bounded-conductor trace function on \(\mathbf F_p\) may have ordinary polynomial interpolation degree nearly \(p-1\). Therefore high interpolation degree or growing Taylor-jet complexity does not imply growing conductor.

## 9.2 A residue modulo \(p\) has no canonical complex statistic

To contradict a bounded-weight compatible system via Sato--Tate or Chebotarev, one needs canonical algebraic numbers \(a_p\) in a fixed number field, with prescribed embeddings and Weil bounds. A value
\[
U_p\in\mathbf F_p
\]
alone has many inequivalent lifts. Its histogram as a least nonnegative residue is not a motivic invariant.

## 9.3 The known zero law is compatible with finite Galois data

The quarter zero sets are exact congruence classes modulo \(24\). Such sets are Chebotarev sets in the cyclotomic extension \(\mathbf Q(\zeta_{24})/\mathbf Q\). Therefore the vanishing pattern is evidence **for compatibility**, not an incompatibility theorem.

## 9.4 What a genuine obstruction would require

A proof that no bounded-conductor object exists would need at least one of:

1. a lower bound on the Swan/conductor of the coefficient-index sheaf as \(p\to\infty\);
2. a proof that no fixed compatible system can have the observed crystalline matrix entries;
3. a violation of Frobenius-polynomial integrality, Galois covariance, duality, or Weil-weight constraints after constructing canonical lifts;
4. a theorem excluding all bounded Jacobi-monomial or induced-Hecke-character formulas.

None of these is proved by the current data.

**[GAP-4: impossibility].** Prove a conductor lower bound or a compatible-system obstruction for the quarter coefficient. Failed searches among a finite list of CM forms do not establish this.

---

# 10. A precise computational fork

The construction above suggests a finite, falsifiable next test.

For primes \(p\equiv1\pmod4\), normalize
\[
U_p=(2A)^{-1}\tau_{(p-1)/4}\pmod p.
\]
Test all Jacobi/Hecke ansätze with character orders dividing \(24\), but impose structural conditions rather than only value matching:

1. **Galois covariance:** conjugating a character exponent by
   \((\mathbf Z/24\mathbf Z)^\times\) must conjugate the proposed algebraic lift.
2. **Norm identities:** Jacobi sums must satisfy their exact product and absolute-norm relations.
3. **Local conductor:** ramification may occur only at \(2\), \(3\), and the prescribed coefficient-system singularities.
4. **Frobenius polynomial:** a proposed rank \(r\) system must yield a fixed-degree polynomial with integral coefficients and the expected duality.
5. **Vanishing:** a rank-one eigenvalue ansatz is immediately excluded; only traces or matrix entries of rank at least two can vanish.

An exact violation rules out that ansatz class. Agreement would provide the missing formula for [GAP-3].

For the sigma class, the direct exponent has full order \(p-1\) by (7.4). Any order-24 model must therefore be tested in an extension-field/descent formulation, not by a quartic character on \(\mathbf F_p\).

---

# 11. Final theorem ledger

## Proved here

1. The Hesse family (1.1), Franel pullback, quadratic twist, and genus-three reversion curve give an explicit fixed geometric family whose local periods are \(\tau(s^4)\) and \(\sigma(s^4)\).
2. The obvious curve \(y^2=q(t)\) has genus zero and its Franel pullback splits.
3. The quartic reversion curve has genus three; the two reflection lifts have quotient genera one and two.
4. The ordinary Hasse--Witt matrices of the genus-three curve and its quotients do not encode the Apéry quarter coefficient; \(p=5\) is an explicit contradiction.
5. The branch Frobenius identities are exactly
   \[
   \tau=S_p\tau(t^p)\quad(\chi_p=+1),
   \]
   and
   \[
   \sigma=S_p\tau(t^p)\quad(\chi_p=-1).
   \]
6. The quarter coefficients are exact entries of the fixed-dimensional coefficient-Cartier operator (5.9), at positions \((1,1)\) and \((3,1)\).
7. The factor \(2A\) is the oriented Hasse invariant of \(y^2=x^3-x\) and is structurally the quartic Kummer factor in the Mellin projector.
8. A finite-order rank-one twist cannot be the residual \(U_p\) if \(U_p\) vanishes at infinitely many unramified primes.
9. The sigma-quarter Mellin exponent has full order \(p-1\) in the class \(p\equiv23\pmod{24}\), so it is not directly an order-24 character.

## Open

- **[GAP-1]** Compare the formal branch Frobenius with a specifically normalized crystalline Frobenius matrix of the explicit pulled-back Hesse family.
- **[GAP-2]** Globalize the coefficient Cartier module and prove bounded conductor.
- **[GAP-3]** Identify the residual \(U_p\) as a fixed Jacobi/Hecke/CM trace, or disprove every such bounded ansatz.
- **[GAP-4]** Prove any general non-motivicity or conductor-growth obstruction.

# Least-confident step

The least-confident step is not an algebraic identity above; it is the proposed geometric upgrade in [GAP-1]. The local period and formal Frobenius identities are exact, but matching their chosen tau/sigma trivializations to the canonical crystalline basis of the explicit Hesse twist is sensitive to covariant-versus-contravariant Frobenius, to the quadratic-twist differential, and to the rational sigma gauge. That comparison should be done by an explicit de Rham basis calculation before the local Cartier package is called a global Hasse--Witt block.

# Reliable reference metadata

- Frits Beukers, **“Consequences of Apéry’s work on \(\zeta(3)\)”**, Utrecht University preprint, 2003. The paper identifies the relevant second-order equation as the Picard--Fuchs equation of the modular elliptic family associated to \(\Gamma_1(6)\).
- Jeffrey D. Achter and Everett W. Howe, **“Hasse-Witt and Cartier-Manin matrices: A warning and a request,”** arXiv:1710.10726. This is the convention reference for the transpose/semilinearity distinction.
- David Harvey and Andrew V. Sutherland, **“Computing Hasse-Witt matrices of hyperelliptic curves in average polynomial time,”** arXiv:1402.3246. This contains the standard hyperelliptic coefficient formula used in Lemma 3.1.
- Jie Zhou, **“GKZ Hypergeometric Series for the Hesse Pencil, Chain Integrals and Orbifold Singularities,”** arXiv:1606.08352. This is a precise Hesse-pencil/hypergeometric reference.
- Masha Vlasenko, **“Higher Hasse-Witt matrices,”** arXiv:1605.06440. This is relevant to the proposed global Frobenius upgrade, but it does not by itself identify the specific tau/sigma local matrix entry.