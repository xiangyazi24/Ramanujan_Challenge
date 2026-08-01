# Horizontal Gross--Koblitz reconstruction: exact formulas and the surviving obstruction

## 0. Outcome ledger

The FAST-RECON outcome is mixed.

1. **[VERIFIED]** In addition to the carry-free Lagrange--Gross--Koblitz
   formula, the rank-two source has a fixed-dimensional finite-field
   compression.  Its Franel Hasse polynomial is one homogenized
   \({}_2F_1\) Jacobi polynomial, equivalently one Euler--Kummer sum, on the
   quadratic cover \(t=x(1-8x)/(1+x)\).  Finite Mellin inversion recovers every
   \(\tau_j\) or \(\sigma_j\).  This is the exact object promised by rank two;
   it is not an \(O(1)\)-term Jacobi evaluation of each reverted coefficient.
2. **[VERIFIED]** Substitution gives both an explicit quadratic Jacobi
   convolution and a shorter Mellin formula for every interior \(b_r\bmod p\).
   At \(p=13\) and \(p=29\), all relevant branch coefficients and all
   \(b_0,\ldots,b_{p-1}\) agree with direct calculation.
3. **[GAP-BGK]** The verified formula is **not** a bounded-length sum of
   Gauss/Jacobi monomials.  Each summand has bounded Gamma/Jacobi complexity,
   but the number of summands grows polynomially with \(j\).  Rank two of the
   differential local system does not by itself remove the nonlinear
   coefficient extraction caused by the algebraic reversion
   \(t=x(1-8x)/(1+x)\).  No bounded-length compression was found or proved
   impossible.
4. **[PARTIAL-HORIZONTAL]** Weil gives \(2\sqrt p\) for every generic fiber of
   the Euler--Kummer lift.  Fiberwise summation first gives \(O(p^{3/2})\) for
   the linear split-cover Mellin sum; Katz's irreducibility and local
   pseudoreflection theorem, followed through the degree-three pullback,
   upgrades this to \(O(p)\).  There is also an exact full-character-period
   two-prime correlation with a genuine power saving.  These estimates concern
   complex split/quadratically-twisted *linear* lifts, not zero indicators.
5. **[NEGATIVE-ZERO-HORIZONTAL]** No new saving follows for the fixed shell
   pair count or \(F_4\).  Additive zero detection exponentiates the value of
   the whole Jacobi/Mellin sum and cannot move the additive character inside
   it.  The surviving object is an incomplete mixed-characteristic CRT
   diagonal.  The only unconditional zero-event saving is the complete-period
   or sliding-shift average from CRT and the pre-existing vertical bound.

The three independent verification scripts are
`research/scripts/codex_horizontal_gk_verify.py`,
`research/scripts/codex_hgk_coefficients.py`, and
`research/scripts/codex_hgk_horizontal.py`.

## 1. Carry-free Gamma and Jacobi atoms

Let \(p\geq 7\), let \(\Gamma_p\) be Morita's Gamma function, and let
\([\alpha]_p\in\{0,\ldots,p-1\}\) be the residue of
\(\alpha\in\mathbf Z_p\).  If

\[
 0\leq s\leq [\alpha]_p\leq p-1,
 \qquad \alpha,\alpha-1,\ldots,\alpha-s+1\in\mathbf Z_p^\times,
\]

define

\[
 \mathcal C_p(\alpha,s)
 :=-\frac{\Gamma_p(\alpha+1)}
 {\Gamma_p(\alpha-s+1)\Gamma_p(s+1)}\pmod p.
 \tag{1.1}
\]

Repeated use of \(\Gamma_p(z+1)=-z\Gamma_p(z)\), with no exceptional
\(p\)-divisible factor, gives

\[
 \boxed{\mathcal C_p(\alpha,s)=\binom{\alpha}{s}\pmod p.}
 \tag{1.2}
\]

This is also a Jacobi sum.  For \(a=[\alpha]_p>0\), extend multiplicative
characters by zero and put

\[
 \overline J_p(s,a)
 :=\overline{\sum_{x\in\mathbf F_p}
       \omega_p^{-s}(x)\omega_p^a(1-x)}
 =\sum_{x\in\mathbf F_p}x^{p-1-s}(1-x)^a\in\mathbf F_p.
\]

Power-sum orthogonality gives

\[
 \boxed{\mathcal C_p(\alpha,s)=(-1)^{s+1}\overline J_p(s,a).}
 \tag{1.3}
\]

The convention-dependent trivial-character endpoints
\((a,s)=(0,0),(p-1,0),(p-1,p-1)\) are kept as the elementary value \(1\).
For \(0<s<a\leq p-1\), the Gauss-sum quotient and Gross--Koblitz give the more
recognizable formula

\[
 \boxed{
 \mathcal C_p(\alpha,s)=(-1)^s
 \frac{\Gamma_p\!\left(\frac{s}{p-1}\right)
       \Gamma_p\!\left(\frac{p-1-a}{p-1}\right)}
      {\Gamma_p\!\left(\frac{p-1-a+s}{p-1}\right)}\pmod p .}
 \tag{1.4}
\]

The endpoint values \(s=0,a\) are \(1\) and are best kept separate because a
trivial character occurs in the Gauss-sum quotient.  Thus (1.1), (1.3), and
(1.4) are three forms of the same carry-free atom.  All their arguments are
affine in the displayed discrete indices.

## 2. Explicit branch formulas

Put

\[
 \Phi(x)=\frac{1+x}{1-8x},\qquad x=t\Phi(x),\qquad
 f_k=\sum_{\ell=0}^k\binom{k}{\ell}^3.
\]

For an arbitrary \(G\), integration by parts in Lagrange--Buermann inversion
gives the useful derivative-free identity

\[
 [t^n]G(x(t))=[x^n]G(x)\Phi(x)^n
 \left(1-x\frac{\Phi'(x)}{\Phi(x)}\right).
 \tag{2.1}
\]

Here

\[
 1-x\frac{\Phi'}{\Phi}
 =\frac{1-16x-8x^2}{(1+x)(1-8x)}.
 \tag{2.2}
\]

This cancels the extra quadratic denominator in the \(\sigma\)-pullback.
Consequently the exact characteristic-zero formulas simplify to

\[
 \boxed{\sigma_n=[x^n]h(x)(1+x)^{n+1/2}(1-8x)^{-n-1},}
 \tag{2.3}
\]

\[
 \boxed{\tau_n=[x^n]h(x)(1-16x-8x^2)
 (1+x)^{n-1/2}(1-8x)^{-n-1}.}
 \tag{2.4}
\]

Define the residual Franel and one-layer coefficient kernels

\[
 \mathcal F_p(k)=\sum_{\ell=0}^k\mathcal C_p(k,\ell)^3,
 \tag{2.5}
\]

\[
 \mathcal K_p(n,\alpha,L)=
 \begin{cases}
 \displaystyle\sum_{a=0}^{L}
 \mathcal C_p(\alpha,a)
 \mathcal C_p(n+L-a,L-a)8^{L-a},&L\geq0,\\
 0,&L<0.
 \end{cases}
 \tag{2.6}
\]

Equations (2.3)--(2.4) now give the explicit Gamma/Jacobi formulas

\[
 \boxed{
 \sigma_n\equiv\sum_{m=0}^n\mathcal F_p(m)
 \mathcal K_p(n,n+1/2,n-m)\pmod p,}
 \tag{2.7}
\]

\[
 \boxed{
 \begin{aligned}
 \tau_n\equiv\sum_{m=0}^n\mathcal F_p(m)\bigl(&
 \mathcal K_p(n,n-1/2,n-m)\\
 &-16\mathcal K_p(n,n-1/2,n-m-1)\\
 &-8\mathcal K_p(n,n-1/2,n-m-2)\bigr)\pmod p.
 \end{aligned}}
 \tag{2.8}
\]

These formulas are valid without a carry correction on exactly the relevant
branch ranges

\[
 D_p=\begin{cases}
 (p-1)/2,&\left(\frac{-6}{p}\right)=+1\quad(\tau),\\
 (p-3)/2,&\left(\frac{-6}{p}\right)=-1\quad(\sigma),
 \end{cases}
 \qquad 0\leq n\leq D_p.
\]

Indeed every lower factorial index in (2.5)--(2.8) is less than \(p\), and
the falling factors are \(p\)-units.  Some upper residues at the branch
endpoint equal \(p-1\); the elementary trivial-character cases recorded
after (1.3) cover them.

### The fixed-dimensional rank-two compression

The growing Lagrange sums are not the only exact description.  Put

\[
 N_p=[-1/3]_p,\qquad M_p=[-2/3]_p,\qquad
 N_p+M_p=p-1,\qquad K_p=\min(N_p,M_p)=\lfloor(p-1)/3\rfloor,
\]

and abbreviate \(B_p(a,k)=\mathcal C_p(a,k)\).  Define

\[
 P_p(z)=\sum_{k=0}^{K_p}B_p(N_p,k)B_p(M_p,k)z^k.
 \tag{2.9}
\]

There is an elementary finite-field Euler formula

\[
 \boxed{P_p(z)=(-1)^{M_p+1}\sum_{y\in\mathbf F_p}
 y^{N_p}(1-y)^{M_p}(1-zy)^{N_p}.}
 \tag{2.10}
\]

Indeed, after expanding the last factor, the coefficient of \(z^k\) reduces
to the unique power-sum exponent \(N_p+k+(M_p-k)=p-1\); it is exactly
\(B_p(N_p,k)B_p(M_p,k)\).  Thus this identity needs no unproved finite-field
hypergeometric transformation.

The homogenized pullback is

\[
 \boxed{
 H_p(x)=\sum_{k=0}^{K_p}B_p(N_p,k)B_p(M_p,k)27^kx^{2k}
 (1-2x)^{p-1-3k}
 =\sum_{m=0}^{p-1}f_mx^m.}
 \tag{2.11}
\]

Away from \(x=1/2\), this is
\((1-2x)^{p-1}P_p(27x^2/(1-2x)^3)\); (2.11) supplies the value at the
missing point.  Let

\[
 \phi(x)=\frac{x(1-8x)}{1+x},\qquad R(x)=1-16x-8x^2,
 \qquad q(t)=t^2-34t+1.
\]

The inverse equation is

\[
 8x^2+(t-1)x+t=0,
 \qquad x=\frac{1-t\pm\sqrt{q(t)}}{16}.
 \tag{2.12}
\]

On either root, in \(\mathbf F_p\) or \(\mathbf F_{p^2}\), the relevant branch
polynomial descends as

\[
 \boxed{
 S_p(\phi(x))=
 \begin{cases}
 H_p(x)/(1+x)^{(p-1)/2},&(-6/p)=+1,\\[2mm]
 H_p(x)/(R(x)(1+x)^{(p-3)/2}),&(-6/p)=-1.
 \end{cases}}
 \tag{2.13}
\]

At removable singularities (2.13) is interpreted by polynomial continuation.
The two inverse roots give the same base-field value.  Since
\(\deg S_p\le p-2\), finite Mellin inversion is exact:

\[
 \boxed{s_{p,j}=-\sum_{t\in\mathbf F_p^\times}S_p(t)t^{-j},
 \qquad 0\le j\le D_p.}
 \tag{2.14}
\]

Equations (2.9)--(2.14) are the bounded-dimensional rank-two object that
survives algebraic reversion: one \(y\)-sum, a degree-two cover, and one Mellin
sum.  They coexist with, but do not imply, an \(O(1)\)-term Jacobi evaluation
for an individual coefficient.

### Why this is not yet the bounded formula requested in the spec

Formula (2.7), fully expanded into Jacobi atoms, has

\[
 L_\sigma(n)=\sum_{m=0}^n(m+1)(n-m+1)
 =\binom{n+3}{3}.
 \tag{2.15}
\]

Formula (2.8) has

\[
 L_\tau(n)=\binom{n+3}{3}+\binom{n+2}{3}+\binom{n+1}{3}.
 \tag{2.16}
\]

Thus the derivative-free reconstruction is \(\Theta(n^3)\), not a short sum.

The distinction is structural.  Gross--Koblitz controls Mellin coefficients
of a finite-field hypergeometric trace function.  Here the classical
\({}_2F_1(1/3,2/3;1)\) first undergoes rational pullback and multiplication,
and then the coefficient in the **different coordinate** \(t\) is extracted
after algebraic reversion.  The degree-two map has bounded geometric degree,
but Lagrange inversion is not a bounded linear operation on its Taylor
coefficients.  Rank two therefore does not imply bounded term count.

The existing characteristic-zero check already disproves one Pochhammer term
for either branch.  It does not exclude every accidental fixed-length sum of
hypergeometric terms.  That stronger exclusion, or a positive compression
identity, remains **[GAP-BGK]**.  Expanding the Jacobi atoms and geometrically
summing their fixed-dimensional index polytope may reorganize (2.7)--(2.8)
into a bounded-dimensional complete sum with boundary strata; that would not
be an \(O(1)\)-term Jacobi evaluation and is not needed for the verified
formulas here.

## 3. The quadratic Jacobi convolution

Let \(U_{p,j}\) be (2.8) for the \(\tau\)-classes and \(V_{p,j}\) be (2.7)
for the \(\sigma\)-classes, extended by zero outside \(0\leq j\leq D_p\).
Since every \(\mathcal C_p\) can be replaced by (1.3), the following are
literal finite sums of products of residual Jacobi sums:

\[
 \boxed{b_r\equiv
 \sum_{i=0}^r U_{p,i}U_{p,r-i}\pmod p}
 \qquad\left(\frac{-6}{p}\right)=+1,
 \tag{3.1}
\]

and, writing \(C_{p,r}=\sum_{i=0}^rV_{p,i}V_{p,r-i}\),

\[
 \boxed{b_r\equiv C_{p,r}-34C_{p,r-1}+C_{p,r-2}\pmod p}
 \qquad\left(\frac{-6}{p}\right)=-1,
 \tag{3.2}
\]

with \(C_{p,-1}=C_{p,-2}=0\).  Equations (3.1)--(3.2) are the exact
"quadratic Jacobi convolution."  They are explicit, but inherit the growing
length of (2.7)--(2.8).

There is also a literal double Mellin form.  For \(d=D_p\), put

\[
 K_{p,r}(t,u)=
 \sum_{i=\max(0,r-d)}^{\min(d,r)}t^{-i}u^{-(r-i)}.
 \tag{3.3}
\]

Then the unfiltered square coefficient is

\[
 \boxed{[T^r]S_p(T)^2=
 \sum_{t,u\in\mathbf F_p^\times}S_p(t)S_p(u)K_{p,r}(t,u).}
 \tag{3.4}
\]

For the \(\tau\)-branch this is \(b_r\); for the \(\sigma\)-branch one applies
the three-term \(1,-34,1\) filter in (3.2).  The kernel (3.3) is a finite
geometric progression, so this is a fixed-dimensional double character-sum
presentation even though its expansion into Jacobi monomials is long.

For Step 3 an even shorter identity is decisive.  Let

\[
 A_p(t)=\sum_{r=0}^{p-1}b_rt^r
 =\begin{cases}S_p(t)^2,&(-6/p)=+1,\\q(t)S_p(t)^2,&(-6/p)=-1.
 \end{cases}
\]

For every possible nontrivial bad index \(1\le r\le p-2\), cyclic and
ordinary coefficient extraction agree, hence

\[
 \boxed{b_r=-\sum_{t\in\mathbf F_p^\times}A_p(t)t^{-r}.}
 \tag{3.5}
\]

At \(r=0\), Mellin extraction aliases the two endpoint coefficients
\(b_0=b_{p-1}=1\); neither endpoint is a zero.  Squaring (2.13) removes the
branch distinction and gives the pointwise cover identity

\[
 \boxed{A_p(\phi(x))=\frac{H_p(x)^2}{(1+x)^{p-1}}=H_p(x)^2,
 \qquad x\in\mathbf F_p\setminus\{-1\}.}
 \tag{3.6}
\]

The cover fiber over \(t\) has exactly

\[
 \#\{x\in\mathbf F_p:\phi(x)=t\}=1+\eta_p(q(t)),
 \tag{3.7}
\]

where \(\eta_p(0)=0\).  Thus the \(x\)-cover computes a split-weighted Mellin
coefficient, not silently the unweighted coefficient in (3.5).

## 4. Exhaustive machine checks at \(p=13,29\)

Running all three independent checks

```text
python3 research/scripts/codex_horizontal_gk_verify.py
python3 research/scripts/codex_hgk_coefficients.py
python3 research/scripts/codex_hgk_horizontal.py
```

includes the outputs

```text
p=13: chi=-1, branch=sigma, D=5, all 6 branch coefficients and all 13 b_r verified
  branch residues: [1, 0, 10, 3, 0, 12]
  distinct carry-free Gamma/Jacobi atoms checked: 55
  fully expanded top-coefficient Jacobi monomials: 56
p=29: chi=+1, branch=tau, D=14, all 15 branch coefficients and all 29 b_r verified
  branch residues: [1, 17, 8, 21, 2, 9, 23, 0, 6, 20, 27, 8, 21, 12, 28]
  distinct carry-free Gamma/Jacobi atoms checked: 309
  fully expanded top-coefficient Jacobi monomials: 1695
quadratic Jacobi convolution agrees with the direct Apéry sum: VERIFIED
p=13: branch=sigma, all 6 coefficients; Jacobi/Gamma, 2F1 pullback, Mellin inversion, and b_r convolution VERIFIED
p=29: branch=tau, all 15 coefficients; Jacobi/Gamma, 2F1 pullback, Mellin inversion, and b_r convolution VERIFIED
CRT averages: pair and four-prime complete periods factor exactly; pair/four sliding means and additive zero detector VERIFIED
Apéry Mellin extraction and split quadratic-cover projection: p=13,29, every interior index VERIFIED
linear Mellin full-period identity: p=13, q=29, gcd(p-1,q-1)=4, period=84, VERIFIED
split-cover Euler/Euler-square two-prime correlations: Weil bounds and 4-term DFT collapses VERIFIED
nonlinear zero-detector interchange: explicit F_5 counterexample VERIFIED
zero events over lcm(p-1,q-1): gcd-stratified identity VERIFIED
```

For every atom used, the script independently compares the direct generalized
binomial product, (1.1), the finite-field power sum (1.3), and the interior
Gross--Koblitz quotient (1.4).  It compares both the original derivative
Lagrange formula and the integrated formulas (2.7)--(2.8) with the rank-two
recurrence, and finally compares (3.1) or (3.2) with the defining Apéry
binomial sum for every \(0\leq r<p\).  The second script independently checks
(2.9)--(2.14) on every point of the quadratic cover, including nonsplit fibers
in \(\mathbf F_{p^2}\).  The third checks (3.5)--(3.7), additive
orthogonality, complete CRT factorization, the full-character-period DFT
identity, and the invalid nonlinear interchange used in the obstruction
below.

## 5. Step 3: the exact horizontal object

For \(n/2<p\leq n\), put \(r_p=n-p\).  In character coordinates

\[
 \omega_p^{r_p}=\omega_p^{n-1},
 \qquad
 \operatorname{ord}(\omega_p^{n-1})
 =\frac{p-1}{\gcd(p-1,n-1)}.
\]

Let \(\mathcal B_p(r)\) denote the right side of (3.1) or (3.2).  For two
primes \(p<q\), the shell pair count can be written, with the appropriate
interval \(I_{p,q,N}\), as

\[
 C_{p,q}(N)=
 \sum_{s\in I_{p,q,N}}
 \mathbf 1_{\mathcal B_p(s+q-p)=0}
 \mathbf 1_{\mathcal B_q(s)=0}.
 \tag{5.1}
\]

Additive orthogonality makes the mixed-characteristic obstruction exact:

\[
 \boxed{
 C_{p,q}(N)=\frac1{pq}
 \sum_{a\in\mathbf F_p}\sum_{c\in\mathbf F_q}
 \sum_{s\in I_{p,q,N}}
 e_p\!\left(a\mathcal B_p(s+q-p)\right)
 e_q\!\left(c\mathcal B_q(s)\right).}
 \tag{5.2}
\]

Likewise, for the falling fourth moment,

\[
 F_4(N)=\sum_{N<n\leq2N}(H(n))_4,
\]

expansion over four distinct primes gives terms

\[
 \frac1{p_1p_2p_3p_4}
 \sum_{a_i\in\mathbf F_{p_i}}
 \sum_{n\in I(\mathbf p,N)}
 \prod_{i=1}^4
 e_{p_i}\!\left(a_i\mathcal B_{p_i}(n-p_i)\right).
 \tag{5.3}
\]

Equations (5.2)--(5.3), followed by either (3.1)--(3.4) or the shorter (3.5),
are the exact requested rewrite.  The shorter form makes the surviving
nonlinearity especially transparent:

\[
 e_p(ab_r)
 =e_p\!\left(-a\sum_{t\ne0}A_p(t)t^{-r}\right)
 =\prod_{t\ne0}e_p(-aA_p(t)t^{-r}),
 \tag{5.4}
\]

not \(\sum_{t\ne0}e_p(-aA_p(t)t^{-r})\).  The false interchange already fails
for a two-term sum over \(\mathbf F_5\), as checked by the third script.  Thus
the zero detector is an additive character of a *global trace value*.  If it
is forced into a sheaf description, its tensor/conductor complexity grows
with \(p\); it is not the bounded-conductor sheaf underlying an individual
Jacobi sum.

### A genuine one-prime Weil object

Let \(\widetilde\omega_p\) be a complex Teichmuller character, set
\(\chi_p=\widetilde\omega_p^{N_p}\), and define

\[
 \mathscr P_p(z)=(-1)^{M_p+1}\sum_{y\in\mathbf F_p}
 \chi_p(y)\chi_p^{-1}(1-y)\chi_p(1-zy).
 \tag{5.5}
\]

Its reduction at the chosen prime over \(p\) is (2.10).  For
\(z\notin\{0,1\}\), the combined Kummer function
\(y(1-zy)/(1-y)\) has four distinct simple zero/pole points.  The associated
rank-one tame sheaf has \(\dim H_c^1=2\) and no \(H_c^0,H_c^2\), so Weil gives

\[
 \boxed{|\mathscr P_p(z)|\le2\sqrt p.}
 \tag{5.6}
\]

At \(z=0,1\) the sum degenerates to a Jacobi or punctured-character sum of
size at most one.  This proves the fiber bound with no monodromy assumption.

Substitute \(z=27x^2/(1-2x)^3\), multiply by the relevant branch and Mellin
Kummer factors, and omit the finitely many singular \(x\).  The resulting
split-cover sum \(\mathscr C_{p,j}^{\tau}\) or
\(\mathscr C_{p,j}^{\sigma}\) is a complete \((x,y)\)-sum.  Fiberwise use of
(5.6) gives unconditionally

\[
 \boxed{|\mathscr C_{p,j}^{\bullet}|\ll p^{3/2}.}
 \tag{5.7}
\]

This is a square-root saving over the raw \(p^2\) terms.  Deligne improves
(5.7) to \(O(p)\) once the rank-two weight-one sheaf obtained after
pushforward in \(y\), pullback in \(x\), and Kummer twisting has no
geometrically constant constituent.  In this case that condition can be
checked, rather than assumed.

Up to a geometrically constant Gauss normalization, Greene's Euler sum (5.5)
is Katz's type-\((2,2)\) hypergeometric sheaf with numerator characters
\((\chi_p^{-1},\chi_p)\) and denominator characters \((\mathbf1,\mathbf1)\).
The lists are disjoint.  Katz's theorem therefore makes the sheaf
geometrically irreducible, and its local monodromy at \(z=1\) is a tame
pseudoreflection.  Its determinant is
\(\mathbf1/(\chi_p^{-1}\chi_p)=\mathbf1\), so this pseudoreflection is a
nontrivial unipotent.

Consequently the identity component of the geometric monodromy contains
\(\mathrm{SL}_2\): an irreducible rank-two algebraic subgroup containing a
nontrivial unipotent cannot be finite, a torus normalizer, or triangular.  On
the complement of its ramification, the nonconstant map
\(z(x)=27x^2/(1-2x)^3\) is finite étale.  Its fundamental group has finite
index in that of the \(z\)-line, so the Zariski closure after pullback still
contains the original identity component.  The pullback is therefore
irreducible; tensoring by any branch/Mellin Kummer character preserves this.
Thus \(H_c^2=0\), and Deligne gives the unconditional improvement

\[
 \boxed{|\mathscr C_{p,j}^{\bullet}|\ll p.}
 \tag{5.7a}
\]

All omitted boundary terms are finite in number and are \(O(p)\) even by the
raw one-variable bound, so they are absorbed.  This closes the
cover-monodromy check; it does not supply descent to the \(t\)-line.

The cover has another exact limitation.  Its reduction is, up to explicitly
removed boundary points,

\[
 -\sum_{t\ne0}(1+\eta_p(q(t)))S_p(t)t^{-j}
 =s_{p,j}+s_{p,j}^{(\eta)},
 \tag{5.8}
\]

where
\(s_{p,j}^{(\eta)}=-\sum_{t\ne0}\eta_p(q(t))S_p(t)t^{-j}\).  Thus it is not
\(s_{p,j}\) alone.  Isolating the unweighted lift requires a compatible
deck-equivariant descent to the \(t\)-line, or an independent construction of
the nonsplit part.  Equality of the reductions on the two sheets is weaker
than such a sheaf descent.  This is **[GAP-DESCENT]**.

If that descent is supplied, Deligne has one further precise target.  A
rank-two weight-one lift \(\mathcal S_p\) makes the lift of (3.5) a Mellin sum
of \(\mathcal S_p\otimes\mathcal S_p\), hence weight two.  Off geometrically
constant constituents its complete sum is \(O(p^{3/2})\).  The determinant
line in
\(\mathcal S_p\otimes\mathcal S_p=\operatorname{Sym}^2\mathcal S_p\oplus
\det\mathcal S_p\) gives an explicit exceptional character: if the Mellin
twist cancels it, a \(p^2\) main term can occur.  Extra exceptional lines would
require dihedral or finite projective monodromy, which the nontrivial unipotent
and \(\mathrm{SL}_2\) identity component exclude here.  Thus the determinant
is the only rank-one constituent.  This classification is about the complex
lift; it still does not detect \(b_r\equiv0\pmod p\).

### The partial two-prime power saving that survives

There is nevertheless a rigorous cross-prime estimate for the *linearized
split-cover lift*.  Puncture the one point above \(z=\infty\), whose omitted
contribution is an explicit parity mode, and put

\[
 \mathscr F_p(t)=
 \sum_{\substack{x:\ \phi(x)=t\\x\text{ nonsingular}}}
 \mathscr P_p(z(x))^2,
 \qquad |\mathscr F_p(t)|\le8p,
 \tag{5.9}
\]

\[
 \mathscr M_p(r)=-\sum_{t\in\mathbf F_p^\times}
 \mathscr F_p(t)\chi_{p,r}(t),
 \qquad
 \chi_{p,r}(g_p^e)=\exp\!\left(-\frac{2\pi i re}{p-1}\right).
 \tag{5.10}
\]

Modulo the Teichmuller prime, (5.10) is the
\((1+\eta_p(q(t)))\)-weighted version of (3.5), with the stated single parity
correction.  Let

\[
 a=p-1,\quad b=q-1,\quad g=\gcd(a,b),\quad L=\operatorname{lcm}(a,b).
\]

Finite Fourier orthogonality gives the exact cross-prime identity

\[
 \boxed{
 \sum_{r=0}^{L-1}\mathscr M_p(r)\overline{\mathscr M_q(r)}
 =L\sum_{k=0}^{g-1}
 \mathscr F_p(g_p^{ak/g})
 \overline{\mathscr F_q(g_q^{bk/g})}.}
 \tag{5.11}
\]

Combining (5.6), the two-sheeted fibers, and (5.11) yields

\[
 \boxed{
 \left|\sum_{r=0}^{L-1}\mathscr M_p(r)
 \overline{\mathscr M_q(r)}\right|
 \le64Lgpq.}
 \tag{5.12}
\]

The unsquared rank-two branch lift has the same DFT identity with fiber bound
\(4\sqrt p\), and therefore the sharper absolute bound
\(16Lg\sqrt{pq}\).  The third verification script checks the underlying
linear DFT identity coefficientwise at \(p=13,q=29\); (5.6) supplies the
archimedean fiber estimate.

The termwise bound for the left side is \(O(Lp^2q^2)\), so (5.12) saves a
factor \(pq/g\); for \(p,q\asymp N\) this is at least one power of \(N\).
This is the requested genuine two-prime power saving newly exposed by the
rank-two Euler formula.  Its scope is exact but limited: it is a complete
character-index period \(L\asymp N^2/g\), it concerns a complex split-cover
lift of linear coefficients, and it gives no information about whether the
defining-characteristic reduction \(b_r\) is zero.

### The strongest unconditional zero-event averages

For the periodic zero sets themselves, CRT gives

\[
 \frac1{pq}\sum_{n\bmod pq}
 \mathbf1_{\mathcal B_p(n-p)=0}
 \mathbf1_{\mathcal B_q(n-q)=0}
 =\frac{|Z_p|}{p}\frac{|Z_q|}{q}
 \ll (pq)^{-1/3}.
 \tag{5.13}
\]

More quantitatively, if \(C_{p,q}(A;X)\) is the length-\(X\) correlation of
the periodic extensions beginning at \(A\), then

\[
 \boxed{
 \sum_{A\bmod pq}C_{p,q}(A;X)=X|Z_p||Z_q|,
 \qquad
 \mathbb E_A C_{p,q}(A;X)\ll X(pq)^{-1/3}.}
 \tag{5.14}
\]

For \(p,q,X\asymp N\), this is \(O(N^{1/3})\), a genuine
\(N^{2/3}\)-saving for the average shift.  It uses only the known
\(|Z_p|\ll p^{2/3}\), so it predates the Gamma/Jacobi reconstruction and does
not control the arithmetically fixed shell interval.

For one fixed interval \(I\) of length \(X\), CRT gives only

\[
 C_{p,q}(I)\le |Z_p||Z_q|\left(\frac{X}{pq}+1\right).
 \tag{5.15}
\]

When \(p,q,X\asymp N\), the \(+1\) term is larger than the trivial interval
bound, so (5.15) contains no fixed-shell saving.  The four-prime sliding-shift
analogue of (5.14) has mean
\(X\prod_i|Z_{p_i}|/\prod_i p_i\ll
X(\prod_i p_i)^{-1/3}\), which is \(O(N^{-1/3})\) when all four primes and
\(X\) are of size \(N\); its averaging period is, however, \(\asymp N^4\).

There is also an exact character-period stratification.  With
\(Z_p^*=Z_p\cap[0,p-2]\), let \(z_{p,c}\) count its elements congruent to
\(c\pmod g\).  Then

\[
 \sum_{r=0}^{L-1}
 \mathbf1_{r\bmod(p-1)\in Z_p^*}
 \mathbf1_{r\bmod(q-1)\in Z_q^*}
 =\sum_{c\bmod g}z_{p,c}z_{q,c}
 \le |Z_p||Z_q|.
 \tag{5.16}
\]

Thus its density is \(\ll g(pq)^{-1/3}\), a power saving whenever
\(g=o((pq)^{1/3})\).  Again, this is a complete index average, not the
length-\(N\) diagonal in (5.2).

Completing the actual \(s\)-sum modulo \(pq\) merely factors it into one
\(p\)-sum and one \(q\)-sum.  The target interval has length \(\asymp N\)
while \(pq\asymp N^2\); for four primes the mismatch is \(N\) versus
\(\asymp N^4\).  Complex Weil bounds do not detect divisibility at the prime
of reduction: the cyclotomic field degree grows with \(p\), so a lift of size
\(O(p)\) can still lie in the prime above \(p\).  The missing input remains a
high-order defining-characteristic Mellin zero-density/dispersion theorem.

## 6. Final gap ledger

The cover-monodromy nondegeneracy is closed by the type-\((2,2)\)
hypergeometric argument leading to (5.7a).  The remaining gaps are:

- **[GAP-BGK]** Compress (2.7)--(2.8) to \(O(1)\) Gauss/Jacobi monomials, or
  prove that no such compression exists.  Rank two alone is insufficient.
- **[GAP-CARTIER]** A special bounded evaluation at a quarter character is not
  excluded.  It would require identifying the reverted local-coefficient
  functional with a controlled Cartier/Frobenius matrix entry; this is extra
  arithmetic, not a consequence of rank two.
- **[GAP-DESCENT]** Construct a deck-equivariant cyclotomic/sheaf descent from
  the quadratic \(x\)-cover to the \(t\)-line.  Without it, the rigorous trace
  bound applies to the split-weighted combination (5.8), not the full branch
  coefficient.
- **[GAP-ZERO]** Supply a defining-characteristic unit, valuation, or norm
  theorem converting a complex trace bound into nonvanishing modulo \(p\).
  Weil--Deligne size estimates alone do not make this conversion.
- **[GAP-INDEX-SHEAF]** Construct a bounded-conductor sheaf for
  \(r\mapsto e_p(a\mathcal B_p(r))\), uniformly in \(p,a\), or prove an
  equivalent complete-sum estimate.  The parameter-side rank-two sheaf does
  not supply this automatically.
- **[GAP-MIXED]** Obtain cancellation on the length-\(N\) CRT diagonal in
  (5.2) or (5.3), uniformly across changing residue characteristics.  Separate
  one-prime Deligne bounds do not address this step.

## 7. Primary references used

- B. Gross and N. Koblitz, [*Gauss sums and the \(p\)-adic
  \(\Gamma\)-function*](https://annals.math.princeton.edu/1979/109-3/p06),
  *Annals of Mathematics* 109 (1979), 569--581.
- J. Greene, [*Hypergeometric functions over finite
  fields*](https://www.ams.org/tran/1987-301-01/S0002-9947-1987-0879564-8/S0002-9947-1987-0879564-8.pdf),
  *Transactions of the AMS* 301 (1987), 77--101.
- N. Katz, [*Exponential Sums and Differential
  Equations*](https://books.google.com/books/about/Exponential_Sums_and_Differential_Equati.html?id=eRuTokPwUa4C),
  Annals of Mathematics Studies 124, Princeton University Press (1990),
  Section 8.4; see also his explicit summary of irreducibility and local
  monodromy in [*G2 and hypergeometric
  sheaves*](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf),
  pp. 3--4.
- P. Deligne, [*La conjecture de Weil:
  II*](https://numdam.org/articles/10.1007/BF02684780/), *Publications
  Mathématiques de l'IHÉS* 52 (1980), 137--252.
- X. Caruso, F. Fürnsinn, D. Vargas-Montoya, and W. Zudilin,
  [*Galois Groups of Apéry-like Series Modulo
  Primes*](https://arxiv.org/abs/2510.23298), arXiv:2510.23298.
