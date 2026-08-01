# Horizontal Gross--Koblitz reconstruction: exact formulas and the surviving obstruction

## 0. Outcome ledger

The FAST-RECON outcome is mixed.

1. **[VERIFIED]** On the relevant branch range, every binomial atom in the
   Lagrange--Buermann formula for \(\tau_j\) or \(\sigma_j\) has both an exact
   carry-free Morita-\(\Gamma_p\) formula and a residual Jacobi-sum formula.
   Gross--Koblitz gives the corresponding three-\(\Gamma_p\) quotient for every
   interior atom.  This yields completely explicit finite multiple sums for
   the branch coefficients.
2. **[VERIFIED]** Substitution gives an explicit quadratic Jacobi convolution
   for every \(b_r\bmod p\).  At \(p=13\) and \(p=29\), all relevant branch
   coefficients and all \(b_0,\ldots,b_{p-1}\) agree with direct calculation.
3. **[GAP-BGK]** The verified formula is **not** a bounded-length sum of
   Gauss/Jacobi monomials.  Each summand has bounded Gamma/Jacobi complexity,
   but the number of summands grows polynomially with \(j\).  Rank two of the
   differential local system does not by itself remove the nonlinear
   coefficient extraction caused by the algebraic reversion
   \(t=x(1-8x)/(1+x)\).  No bounded-length compression was found or proved
   impossible.
4. **[NEGATIVE-HORIZONTAL]** Consequently the formula does not create a
   bounded-conductor trace function in the index \(r\), and ordinary
   Weil--Deligne estimates give no new two-prime saving.  The exact
   cross-prime expression below is an incomplete mixed-characteristic CRT
   diagonal.  Completion over a full CRT period merely factorizes it into
   separate one-prime sums.

The verification script is
`research/scripts/codex_horizontal_gk_verify.py`.

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

### Why this is not yet the bounded formula requested in the spec

Formula (2.7), fully expanded into Jacobi atoms, has

\[
 L_\sigma(n)=\sum_{m=0}^n(m+1)(n-m+1)
 =\binom{n+3}{3}.
 \tag{2.9}
\]

Formula (2.8) has

\[
 L_\tau(n)=\binom{n+3}{3}+\binom{n+2}{3}+\binom{n+1}{3}.
 \tag{2.10}
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

## 4. Exhaustive machine checks at \(p=13,29\)

Running

```text
python3 research/scripts/codex_horizontal_gk_verify.py
```

produces

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
```

For every atom used, the script independently compares the direct generalized
binomial product, (1.1), the finite-field power sum (1.3), and the interior
Gross--Koblitz quotient (1.4).  It compares both the original derivative
Lagrange formula and the integrated formulas (2.7)--(2.8) with the rank-two
recurrence, and finally compares (3.1) or (3.2) with the defining Apéry
binomial sum for every \(0\leq r<p\).

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

Equations (5.2)--(5.3), followed by (1.3) and (3.1)--(3.2), are the exact
requested rewrite.  They also show precisely why it is not yet a usable
complete exponential sum over one bounded-complexity variety.

### What Weil--Deligne can and cannot do

- At one fixed prime it can estimate a genuinely complete trace-function sum
  of bounded conductor once such a sheaf is supplied.  Formula (1.3) also
  evaluates each individual Jacobi atom exactly.
- Expanding the atoms in (2.7)--(2.8) and geometrically summing their
  fixed-dimensional index polytope may produce a bounded-dimensional complete
  sum, with singular boundary strata treated separately.  Deligne could bound
  the complex realization of such a fixed-coefficient sum after a
  nondegeneracy audit.  This would still be neither a bounded Jacobi-product
  evaluation nor a sheaf in the moving index \(r\), and a complex size bound
  would not decide reduction to zero at the defining prime.
- The present \(r\mapsto\mathcal B_p(r)\) formula is a growing multiple sum,
  not a proved bounded-conductor trace function in \(r\).  Therefore one
  cannot apply a uniform Deligne bound to the inner \(s\)-sum in (5.2).
- Completing \(s\) modulo \(pq\) makes the sum factor, by CRT, into one
  \(p\)-sum and one \(q\)-sum.  This supplies no interaction between the two
  characteristics.  The relevant interval has length \(\asymp N\), whereas
  the complete period is \(pq\asymp N^2\); for four primes the mismatch is
  \(N\) versus \(\asymp N^4\).
- The additive zero detector itself averages over all \(a\in\mathbf F_p\).
  Pointwise square-root cancellation for isolated nonzero modes, even if it
  were available, would not automatically survive this mode sum and the short
  mixed-characteristic diagonal.

There is one complete-period identity:

\[
 \frac1{pq}\sum_{n\bmod pq}
 \mathbf1_{\mathcal B_p(n-p)=0}
 \mathbf1_{\mathcal B_q(n-q)=0}
 =\frac{|Z_p|}{p}\frac{|Z_q|}{q}
 \ll p^{-1/3}q^{-1/3}.
 \tag{5.4}
\]

The displayed power saving uses the known \(|Z_p|\ll p^{2/3}\), but it is
just CRT plus the vertical bound and predates the Gamma/Jacobi reconstruction.
No new incomplete-interval or genuine two-prime power saving follows from the
present formula.  The missing input remains a high-order,
defining-characteristic Mellin diagonal zero-density/dispersion theorem.

## 6. Final gap ledger

- **[GAP-BGK]** Compress (2.7)--(2.8) to \(O(1)\) Gauss/Jacobi monomials, or
  prove that no such compression exists.  Rank two alone is insufficient.
- **[GAP-CARTIER]** A special bounded evaluation at a quarter character is not
  excluded.  It would require identifying the reverted local-coefficient
  functional with a controlled Cartier/Frobenius matrix entry; this is extra
  arithmetic, not a consequence of rank two.
- **[GAP-INDEX-SHEAF]** Construct a bounded-conductor sheaf for
  \(r\mapsto e_p(a\mathcal B_p(r))\), uniformly in \(p,a\), or prove an
  equivalent complete-sum estimate.  The parameter-side rank-two sheaf does
  not supply this automatically.
- **[GAP-MIXED]** Obtain cancellation on the length-\(N\) CRT diagonal in
  (5.2) or (5.3), uniformly across changing residue characteristics.  Separate
  one-prime Deligne bounds do not address this step.
