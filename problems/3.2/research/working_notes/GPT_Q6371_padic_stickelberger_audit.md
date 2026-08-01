ANSWER Q6371 8e28d1ef

# The defining-characteristic Mellin route for Apéry values

## Executive verdict

There are now **two exact and complementary compressed presentations** of the Apéry coefficient \(b_rmod p\).

1. A **non-circular one-character Jacobi convolution**:
   
   
   \[
   b_r\equiv-\frac1{p-1}\sum_{\eta\in\widehat{\mathbf F_p^\times}}
   J(\varphi\eta\chi_r,\varphi\eta^{-1})^2
   J(\varphi\eta\chi_r^{-1},\varphi\eta^{-1})^2
   \pmod{\mathfrak p},
   \tag{A}
   \]
   
   where \(\chi_r=\omega^r\), \(\varphi=\omega^{(p-1)/2}\), and \(1\le r\le p-2\). This is one complete sum over the character group; \(r\) occurs only through the twists \(\chi_r^{\pm1}\).

2. An exact **Franel–Mellin two-summand formula**:
   
   \[
   b_r\equiv
   -\sum_x H_p(x)^2\,\omega^{-r}(\phi(x))
   +\sum_t\chi_2(q(t))A_p(t)\,\omega^{-r}(t)
   \pmod p,
   \tag{B}
   \]
   
   with explicitly stated domains below, \(\phi(x)=x(1-8x)/(1+x)\), and \(q(t)=t^2-34t+1\). This is a difference of two one-variable Mellin transforms. The finite-field identity is proved exactly. Promoting the second amplitude to an *independent integral bounded-conductor object*, rather than an expression involving \(A_p\) itself, is still the object-level Cartier/descent issue.

Formula (A) is the correct coordinate for Gross–Koblitz and Stickelberger. Its valuation ledger can be computed completely. The result is negative for the proposed shortcut:

> For \(\rho=\min(r,p-1-r)\), the lowest Stickelberger slope has multiplicity \(\rho+1\). Its contribution is exactly
> \[
> \sum_{k=0}^{\rho}\binom{\rho}{k}^2\binom{\rho+k}{k}^2=b_\rho\equiv b_r\pmod p.
> \]
> Every one of these \(\rho+1\) summands is a \(p\)-adic unit. Thus \(p\mid b_r\) is cancellation among many equal-slope Gamma monomials, not the disappearance of all unit terms and not a unique-minimum argument.

I find **no Ax–Katz, Adolphson–Sperber, Katz-Mellin, Heath-Brown–Patterson, or Wan theorem** that gives

\[
\#\{r<p:b_r\equiv0\pmod p\}=o(p)
\quad\text{or}\quad O(p^{1-\delta})
\]

from this presentation. The known theorems control, respectively, divisibility of point counts, Newton polygons of fixed exponential-sum \(L\)-functions, complex equidistribution of normalized Mellin transforms, prime-aspect complex distributions of low-order Gauss sums, or \(p\)-adic meromorphy of unit-root \(L\)-functions. None controls **exact reduction to zero at the coefficient prime while a tame character of order comparable to \(p\) varies**.

The clean missing statement is:

> **[GAP-DCM] Defining-characteristic Mellin zero-density.** For some \(\delta>0\), uniformly in primes \(p\), the explicit function \(\mathcal B_p(\chi)\) in (A) satisfies
> \[
> \#\bigl\{\chi\in\widehat{\mathbf F_p^\times}:
> \operatorname{ord}(\chi)>(\log p)^A,
> \ \overline{\mathcal B_p(\chi)}=0\in\mathbf F_p\bigr\}
> \ll_A p^{1-\delta}.
> \tag{DCM}
> \]
> The random-scale prediction is \(O(1)\). To improve the proved all-value bound \(N_p(c)\le8p^{3/4}\), one needs \(\delta>1/4\); to improve the zero-fibre exponent \(2/3\), one needs \(\delta>1/3\).

The \(p\)-adic route does **not** currently explain that zeros at generic \(r\) are impossible. It explains why low-order twists can have exceptional CM/Cartier zeros and why those structured twists are sparse in a fixed moving column. At high order, it leaves a flat block of units whose cancellation is expected to be rare but is not structurally prohibited. The census fact \(F_4=0\) is cross-prime anti-alignment; no one-prime \(p\)-adic theorem here explains it.

---

# 1. Exact Mellin compression in at most two variables

Let

\[
b_n=\sum_{k=0}^n\binom nk^2\binom{n+k}{k}^2,
\qquad
A_p(t)=\sum_{n=0}^{p-1}b_nt^n\in\mathbf F_p[t].
\]

Let the Franel numbers and their truncation be

\[
f_n=\sum_{k=0}^n\binom nk^3,
\qquad
H_p(x)=\sum_{n=0}^{p-1}f_nx^n.
\]

Put

\[
\phi(x)=\frac{x(1-8x)}{1+x},
\qquad
q(t)=t^2-34t+1.
\]

All characters below are extended by zero at \(0\). The Teichmüller character is denoted by \(\omega\), and \(\chi_2=\omega^{(p-1)/2}\) is quadratic.

## 1.1 The pullback identity

The characteristic-zero Apéry–Franel identity is

\[
F(\phi(x))=(1+x)h(x)^2,
\tag{1.1}
\]

where \(F=\sum b_nt^n\) and \(h=\sum f_nx^n\). The Lucas congruences give, in characteristic \(p\),

\[
F(t)=A_p(t)F(t)^p,
\qquad
h(x)=H_p(x)h(x)^p.
\tag{1.2}
\]

Substituting (1.1) into the first identity and using the second gives

\[
(1+x)H_p(x)^2h(x)^{2p}
=A_p(\phi(x))(1+x)^ph(x)^{2p}.
\]

Since \(h(0)=1\), cancellation is legitimate in the formal series ring. Therefore

\[
\boxed{
A_p(\phi(x))(1+x)^{p-1}=H_p(x)^2
}
\tag{1.3}
\]

in \(\mathbf F_p(x)\). In particular, for every \(x\in\mathbf F_p\setminus\{-1\}\), Fermat gives

\[
A_p(\phi(x))=H_p(x)^2.
\tag{1.4}
\]

## 1.2 The fibres of the quadratic pullback

The equation \(\phi(x)=t\) is

\[
8x^2+(t-1)x+t=0.
\tag{1.5}
\]

Its discriminant is

\[
(t-1)^2-32t=t^2-34t+1=q(t).
\]

For \(p>3\), \(x=-1\) is never a solution of (1.5), since substitution gives \(9\). Hence

\[
\#\{x\in\mathbf F_p:\phi(x)=t\}=1+\chi_2(q(t)).
\tag{1.6}
\]

The two points above \(t=0\) are \(x=0\) and \(x=1/8\); they are simply omitted in nontrivial Mellin coefficients.

## 1.3 The exact two-summand Mellin formula

For \(1\le r\le p-2\), define

\[
\mathcal F_p(r)=
\sum_{\substack{x\in\mathbf F_p\\x\ne-1,\ \phi(x)\ne0}}
H_p(x)^2\,\phi(x)^{-r},
\tag{1.7}
\]

and

\[
\mathcal Q_p(r)=
\sum_{t\in\mathbf F_p^\times}
\chi_2(q(t))A_p(t)t^{-r}.
\tag{1.8}
\]

Equivalently, these powers are the reductions of the character twists \(\omega^{-r}(\phi(x))\) and \(\omega^{-r}(t)\).

Using (1.4) and grouping by the fibres (1.6),

\[
\begin{aligned}
\mathcal F_p(r)
&=\sum_{t\in\mathbf F_p^\times}
(1+\chi_2(q(t)))A_p(t)t^{-r}\\
&=\mathcal M_p(r)+\mathcal Q_p(r),
\end{aligned}
\tag{1.9}
\]

where

\[
\mathcal M_p(r)=\sum_{t\in\mathbf F_p^\times}A_p(t)t^{-r}.
\]

Power-sum orthogonality gives

\[
\begin{aligned}
\mathcal M_p(r)
&=\sum_{n=0}^{p-1}b_n
\sum_{t\in\mathbf F_p^\times}t^{n-r}\\
&=(p-1)b_r=-b_r,
\end{aligned}
\tag{1.10}
\]

because \(n=r\) is the only possible congruence \(n\equiv r\pmod{p-1}\) for \(1\le r\le p-2\). Therefore

\[
\boxed{
b_r=-\mathcal F_p(r)+\mathcal Q_p(r)
}
\qquad(1\le r\le p-2).
\tag{1.11}
\]

This proves the requested compression as the difference of two complete one-variable sums.

### Circularity warning

The finite-field identity (1.11) is exact, but \(\mathcal Q_p\) still contains \(A_p\). Calling the pair in (1.11) an *independent fixed bounded-conductor pair* requires an integral descent identifying \(\chi_2(q)A_p\) with the correct Cartier/Frobenius companion of the Franel-square pushforward. The CFVZ square-root factorization makes this plausible and supplies the expected rank-two objects, but the pointwise identity alone does not remove the circular occurrence of \(A_p\).

Thus:

- (1.11) is a proved algebraic compression;
- the non-circular object-level statement is the still-relevant **[GAP-CARTIER]**;
- the Jacobi formula below is non-circular and needs no such qualification.

## 1.4 A non-circular one-character Jacobi formula

Let \(N=p-1\), let \(\psi=\omega\) generate \(\widehat{\mathbf F_p^\times}\), and let

\[
\varphi=\psi^{N/2}.
\]

For multiplicative characters \(\alpha,\beta\), write

\[
J(\alpha,\beta)=\sum_{x\in\mathbf F_p}\alpha(x)\beta(1-x).
\]

Fix the standard prime \(\mathfrak p_\psi\) above \(p\) for which Teichmüller values reduce to the corresponding elements of \(\mathbf F_p\). Then, for \(1\le r\le N-1\),

\[
\boxed{
 b_r\equiv-\frac1N\sum_{t=0}^{N-1}
 J(\varphi\psi^{t+r},\psi^{N/2-t})^2
 J(\varphi\psi^{t-r},\psi^{N/2-t})^2
 \pmod{\mathfrak p_\psi}.
}
\tag{1.12}
\]

Equivalently, with \(\eta=\psi^t\) and \(\chi_r=\psi^r\),

\[
\boxed{
\mathcal B_p(\chi_r):=
-\frac1N\sum_{\eta\in\widehat{\mathbf F_p^\times}}
J(\varphi\eta\chi_r,\varphi\eta^{-1})^2
J(\varphi\eta\chi_r^{-1},\varphi\eta^{-1})^2
\equiv b_r.
}
\tag{1.13}
\]

This is a single complete character sum. It is the cleanest unconditional answer to part (1): \(r\) enters only via \(\chi_r^{\pm1}\), and no coefficient \(b_j\) occurs on the right. The endpoint classes \(r=0,N\) involve the familiar trivial-character normalization correction and are better handled directly; \(b_0\equiv b_N\equiv1\pmod p\).

## 1.5 The square-root branch coefficients

In the CFVZ square classes let

\[
T_p(t)=\sum_{j=0}^{(p-1)/2}\tau_jt^j,
\qquad A_p=T_p^2,
\]

and in the other four classes let

\[
S_p(t)=\sum_{j=0}^{(p-3)/2}\sigma_jt^j,
\qquad A_p=qS_p^2.
\]

Ordinary Mellin inversion gives the exact one-variable formulas

\[
\boxed{
\tau_j=-\sum_{t\in\mathbf F_p^\times}T_p(t)t^{-j},
\qquad
\sigma_j=-\sum_{t\in\mathbf F_p^\times}S_p(t)t^{-j}
}
\tag{1.14}
\]

through the respective branch degrees. These identities are formally complete but are tautological unless \(T_p,S_p\) are independently realized as trace/Cartier functions of bounded conductor. The proved non-circular formulas obtained from Lagrange–Bürmann and Gross–Koblitz are explicit but have \(\Theta(j^3)\) Gamma/Jacobi atoms. The old fixed-dimensional geometric-summation compression reorganizes those atoms; it does not turn them into \(O(1)\) Jacobi monomials.

---

# 2. Stickelberger and Gross–Koblitz: the exact valuation ledger

Formula (1.12) is ideal for termwise \(p\)-adic analysis.

## 2.1 Stickelberger carries for a Jacobi sum

Choose \(\pi\) with \(\pi^{p-1}=-p\). For \(0\le a\le N-1\), Gross–Koblitz gives, up to the conventional choice of additive character,

\[
g(\omega^{-a})=-\pi^a\Gamma_p\!\left(\frac aN\right).
\tag{2.1}
\]

If neither character nor their product is trivial, then

\[
J(\omega^{-a},\omega^{-b})
=\frac{g(\omega^{-a})g(\omega^{-b})}
{g(\omega^{-\langle a+b\rangle_N})},
\]

and therefore

\[
\boxed{
v_pJ(\omega^{-a},\omega^{-b})
=\frac{a+b-\langle a+b\rangle_N}{N}
=\left\lfloor\frac{a+b}{N}\right\rfloor.
}
\tag{2.2}
\]

Thus a Jacobi sum has slope \(0\) or \(1\), determined by one base-\(p\) carry. When the product character is trivial, \(J(\alpha,\alpha^{-1})=-\alpha(-1)\), a unit; when one character is trivial, the nondegenerate Jacobi sum is also an explicit unit. These are finitely many boundary strata.

For the \(t\)-th term of (1.12), define the carry function using least residues of negative character exponents. Away from the trivial-character strata its valuation is

\[
\begin{aligned}
\nu_{r}(t)=&\ 2\,\operatorname{car}
\left(\frac N2+t+r,\frac N2-t\right)\\
&+2\,\operatorname{car}
\left(\frac N2+t-r,\frac N2-t\right),
\end{aligned}
\tag{2.3}
\]

where each \(\operatorname{car}\) is \(0\) or \(1\). Hence

\[
\nu_r(t)\in\{0,2,4\},
\tag{2.4}
\]

again with explicit unit corrections at the trivial-character endpoints. This completely determines the termwise Newton/Stickelberger filtration.

## 2.2 The slope-zero block

Let

\[
\rho=\min(r,p-1-r).
\]

The Apéry reflection gives \(b_r\equiv b_\rho\pmod p\). The slope-zero sector of (1.12), after the standard terminating-hypergeometric reindexing, is exactly

\[
\sum_{k=0}^{\rho}inom{\rho}{k}^2\binom{\rho+k}{k}^2.
\tag{2.5}
\]

All remaining terms of the complete Jacobi sum have positive even valuation and disappear modulo \(p\).

There is an especially transparent Morita-Gamma form. Since \(0\le k\le\rho\le(p-1)/2\), all factorials are below \(p\), and

\[
\binom{\rho}{k}\binom{\rho+k}{k}
=\frac{(\rho+k)!}{k!^2(\rho-k)!}.
\]

Morita’s functional equation gives, with the sign disappearing after squaring,

\[
\boxed{
U_{\rho,k}:=
\left(
\frac{\Gamma_p(\rho+k+1)}
{\Gamma_p(k+1)^2\Gamma_p(\rho-k+1)}
\right)^2
=\binom{\rho}{k}^2\binom{\rho+k}{k}^2.
}
\tag{2.6}
\]

Every factor \(\Gamma_p(\cdot)\) here is a \(p\)-adic unit. Consequently

\[
\boxed{
 p\mid b_r
 \iff
 \sum_{k=0}^{\rho}U_{\rho,k}\equiv0\pmod p,
 \qquad U_{\rho,k}\in\mathbf Z_p^\times.
}
\tag{2.7}
\]

This is the precise \(p\)-adic cancellation statement requested in part (2).

## 2.3 What the valuation computation proves—and what it rules out

It proves:

- every non-boundary term of the complete character formula has an explicit Stickelberger slope;
- the positive-slope pieces are irrelevant modulo \(p\);
- the unit-root/slope-zero piece is an explicit \(\rho+1\)-term Gamma sum;
- the divisibility event is entirely inside that unit block.

It rules out the most tempting shortcut:

> There is no unique lowest-valuation monomial whose nonzero residue forces \(b_r\not\equiv0\). The minimum is attained \(\rho+1\) times, which is \(\asymp p\) for generic \(r\).

Gross–Koblitz supplies the unit residues, but not cancellation control among them. The problem has not become an Ax–Katz “all terms acquire positive valuation” event. It has become a value-distribution problem for a long \(p\)-adic Gamma sum.

A useful comparison is the carry-free atom used in the branch reconstruction. For \(0<s<a<N\),

\[
\binom as\equiv(-1)^s
\frac{\Gamma_p(s/N)\Gamma_p((N-a)/N)}
{\Gamma_p((N-a+s)/N)}\pmod p.
\tag{2.8}
\]

Thus every branch atom is a bounded Gamma quotient, but the number of atoms grows. “Bounded Gamma complexity per atom” and “bounded number of atoms” are different assertions.

---

# 3. Why Ax–Katz and Adolphson–Sperber do not count the twist zeros

## 3.1 Ax–Katz is pointed in the wrong direction

Ax–Katz and its Newton-polytope refinements give lower bounds for the \(p\)-divisibility of

- the number of common zeros of polynomial systems, or
- exponential sums attached to fixed polynomial data.

Our parameter \(r\) is not a bounded-degree algebraic coordinate. Choosing a generator \(g\in\mathbf F_p^\times\), the character \(\omega^r(x)\) is \(g^{r\log_gx}\). The discrete logarithm prevents \(r\mapsto\omega^r(x)\) from being a fixed-degree algebraic function on \(\mathbf A^1\).

One can interpolate \(r\mapsto b_r\) by a polynomial in \(r\) of degree at most \(p-1\), but this destroys the hypothesis that makes Ax–Katz useful. A one-variable polynomial of degree \(\asymp p\) may have \(\asymp p\) zeros, and the Ax–Katz divisibility exponent is then trivial. Moreover, Ax–Katz controls divisibility of a zero count, not an upper bound for the number of parameters at which a character sum is zero.

The identity

\[
1_{\{\mathcal B_p(\chi)=0\}}
=\frac1p\sum_{u\in\mathbf F_p}e_p(u\mathcal B_p(\chi))
\tag{3.1}
\]

does not fix the issue. Exponentiating the *value of a trace function* is not a bounded-conductor tensor operation on the original sheaf. The interpolation conductor again grows with \(p\).

## 3.2 Adolphson–Sperber controls slopes of \(L\)-functions, not trace cancellation

Adolphson–Sperber’s Newton-polyhedron theory and their work on twisted exponential sums control:

- the degree of the associated \(L\)-function;
- lower bounds for its Newton polygon;
- in nondegenerate cases, cohomological concentration and sharpness of slope bounds.

Those results are powerful when the desired event is “all Frobenius eigenvalues have valuation at least \(\lambda\)” or “the Newton polygon jumps.” Here every \(r\) already has a slope-zero block. The event \(b_r\equiv0\) is cancellation among slope-zero contributions. A trace can vanish while the Newton polygon and all individual slopes remain generic.

Generic-Newton-polygon and Hasse-polynomial results, such as those of Blache–Férard, likewise identify parameters where the Newton polygon is non-generic. They do not identify parameters where one selected trace coefficient happens to be zero modulo the defining prime.

Therefore no direct Ax–Katz or Adolphson–Sperber theorem yields \(o(p)\), let alone a power saving, for the \(r\)-zero locus.

---

# 4. What the actual Gauss-sum and Mellin literature gives

The following distinctions are load-bearing.

## 4.1 Gross–Koblitz and Stickelberger

**B. Gross and N. Koblitz**, “Gauss sums and the \(p\)-adic \(\Gamma\)-function,” *Annals of Mathematics* 109 (1979), 569–581, DOI 10.2307/1971226.

This gives (2.1) and hence the complete valuation/carry ledger. It is termwise. It contains no theorem bounding the number of parameter values for which a sum of many equal-slope terms cancels.

## 4.2 Ax–Katz and Newton-polytope refinements

- **J. Ax**, “Zeroes of polynomials over finite fields,” *American Journal of Mathematics* 86 (1964), 255–261, DOI 10.2307/2373163.
- **N. Katz**, “On a theorem of Ax,” *American Journal of Mathematics* 93 (1971), 485–499, DOI 10.2307/2373389.
- **A. Adolphson and S. Sperber**, “\(p\)-adic estimates for exponential sums and the theorem of Chevalley–Warning,” *Annales scientifiques de l’École Normale Supérieure* 20 (1987), 545–556, DOI 10.24033/asens.1543.
- **A. Adolphson and S. Sperber**, “Exponential sums and Newton polyhedra: cohomology and estimates,” *Annals of Mathematics* 130 (1989), 367–406, DOI 10.2307/1971424.
- **A. Adolphson and S. Sperber**, “On twisted exponential sums,” *Mathematische Annalen* 290 (1991), 713–726, DOI 10.1007/BF01459269.

These provide point-count divisibility and Newton-polygon estimates for fixed algebraic data. They do not give a zero-density theorem for exact trace cancellation as a high-order tame character varies.

## 4.3 Heath-Brown–Patterson and Patterson

- **D. R. Heath-Brown and S. J. Patterson**, “The distribution of Kummer sums at prime arguments,” *Journal für die reine und angewandte Mathematik* 310 (1979), 111–130, DOI 10.1515/crll.1979.310.111.
- **S. J. Patterson**, “The distribution of general Gauss sums and similar arithmetic functions at prime arguments,” *Proceedings of the London Mathematical Society* 54 (1987), 193–215, DOI 10.1112/plms/s3-54.2.193.

These are prime-aspect distribution results for complex normalized Gauss/Kummer sums, especially low fixed orders such as cubic sums. They do not count, at one prime, high-order characters whose algebraic values reduce to zero at a prime above the same \(p\).

Complex equidistribution cannot by itself see divisibility. Two algebraic integers can have indistinguishable normalized complex arguments while one is divisible by \(\mathfrak p\) and the other is not.

## 4.4 Katz’s finite-field Mellin equidistribution

- **N. Katz**, *Gauss Sums, Kloosterman Sums, and Monodromy Groups*, Princeton University Press, 1988.
- **N. Katz**, *Convolution and Equidistribution: Sato–Tate Theorems for Finite-Field Mellin Transforms*, Annals of Mathematics Studies 180, Princeton University Press, 2012, DOI 10.23943/princeton/9780691153308.001.0001.
- **J. Fresán**, “Équirépartition de sommes exponentielles (travaux de Katz),” Séminaire Bourbaki, arXiv:1910.08572.
- **A. Forey, J. Fresán, and E. Kowalski**, “Arithmetic Fourier transforms over finite fields: generic vanishing, convolution, and equidistribution,” arXiv:2109.11961.

Katz’s theorem studies normalized Mellin transforms of a fixed suitable sheaf as characters vary, principally in an extension-field limit. Monodromy governs the limiting **complex/\(\ell\)-adic Frobenius conjugacy classes**.

Our question is finer in three ways:

1. the base field is \(\mathbf F_p\) itself rather than \(\mathbf F_{p^m}\) with fixed \(p\) and \(m\to\infty\);
2. the coefficient prime is the defining characteristic \(p\), not a fixed auxiliary \(\ell\ne p\);
3. the event is exact reduction \(\overline{\mathcal B_p(\chi)}=0\), a lattice-scale event of expected probability \(1/p\), not membership in a fixed positive-measure region of a compact monodromy group.

Even perfect Sato–Tate equidistribution at fixed macroscopic resolution would not imply a local limit theorem at resolution \(1/p\). Katz’s generic-vanishing theorem is cohomological vanishing of unwanted degrees for generic twists; it is not vanishing of the numerical Mellin trace modulo \(p\).

## 4.5 Wan and \(p\)-adic variation

- **D. Wan**, “Dwork’s conjecture on unit root zeta functions,” *Annals of Mathematics* 150 (1999), 867–927, DOI 10.2307/121058.
- **C. D. Haessig and S. Sperber**, “\(p\)-adic variation of unit root \(L\)-functions,” arXiv:1512.06258.

Wan proves \(p\)-adic meromorphy of unit-root \(L\)-functions coming from geometry; subsequent work studies variation in suitable \(p\)-adic families. This controls analytic continuation and slopes of reciprocal roots. It does not provide a uniform bound on the number of tame torsion characters of order dividing \(p-1\) at which one selected Mellin trace is zero modulo \(p\).

The tame characters are spread among roots of unity of order prime to \(p\); as \(p\) varies they do not form one fixed \(p\)-adic disk on which Strassmann gives a bounded zero count. A useful Strassmann argument would require a single integral analytic function with uniformly bounded Weierstrass degree, which is exactly what has not been constructed.

## 4.6 Related \(p\)-adic valuation and hypergeometric zero results

- **D. J. Katz, P. Langevin, S. Lee, and Y. Sapozhnikov**, “The \(p\)-adic valuations of Weil sums of binomials,” *Journal of Number Theory* 181 (2017), 1–26, DOI 10.1016/j.jnt.2017.05.020.
- **R. Blache and É. Férard**, “Newton polygons for twisted exponential sums and polynomials \(P(x^d)\),” arXiv:math/0702502; see also *Journal of Number Theory* 123 (2007), 456–472, DOI 10.1016/j.jnt.2006.06.009, for generic Newton strata.
- **D. McCarthy**, “Extending Gaussian hypergeometric series to the \(p\)-adic setting,” *International Journal of Number Theory* 8 (2012), 1581–1612, DOI 10.1142/S1793042112500844.
- **N. Saikia**, “Zeros of hypergeometric functions in the \(p\)-adic setting,” *Ramanujan Journal* 61 (2023), 1339–1355, DOI 10.1007/s11139-022-00646-5.

The Weil-sum paper studies minimum valuations as an additive parameter varies. Blache–Férard studies generic Newton polygons. McCarthy and Saikia prove identities and zero/nonzero classifications for special hypergeometric families and special parameters. None supplies a general \(O(p^{1-\delta})\) bound for the zero locus of the high-order-character family (1.13).

### Literature verdict

Among the sources above, **none gives a nontrivial zero-density estimate for (1.13)**. The nearest theorem shape is Katz-style Mellin equidistribution after an integral sheaf has been identified, but an additional defining-characteristic local-limit theorem is required.

---

# 5. The precise missing theorem

The cleanest formulation avoids all sheaf-identification ambiguity by using the explicit Jacobi expression itself.

For \(\chi\in\widehat{\mathbf F_p^\times}\), define

\[
\mathcal B_p(\chi)=
-\frac1{p-1}\sum_{\eta}
J(\varphi\eta\chi,\varphi\eta^{-1})^2
J(\varphi\eta\chi^{-1},\varphi\eta^{-1})^2.
\tag{5.1}
\]

Then \(b_r\equiv\mathcal B_p(\omega^r)\) for \(1\le r\le p-2\).

> **[GAP-DCM] Defining-characteristic Mellin zero-density.** There exist \(\delta>0\) and, for every fixed \(A>0\), a constant \(C_A\) such that
> \[
> \#\left\{\chi\in\widehat{\mathbf F_p^\times}:
> \operatorname{ord}(\chi)>(\log p)^A,
> \ \mathcal B_p(\chi)\equiv0\pmod{\mathfrak p_\psi}
> \right\}
> \le C_Ap^{1-\delta}
> \tag{5.2}
> \]
> for every prime \(p>3\), or at least outside a quantitatively sparse exceptional set of primes.

Natural strengthenings are:

\[
\#\{\chi:\mathcal B_p(\chi)=a\}
=1+O(p^{1-\delta})
\quad\text{uniformly in }a\in\mathbf F_p,
\tag{5.3}
\]

or the random-scale estimate \(p^{o(1)}\), expected to be \(O(1)\).

This is accurately described as

> **zero-density on the diagonal of a high-order-character Mellin family in defining characteristic.**

“Diagonal” means that the field characteristic, the coefficient prime at which divisibility is tested, and the order scale of the varying character all move together.

## 5.1 Why ordinary monodromy is not yet enough

Suppose [GAP-CARTIER] identifies the Franel–Mellin pair with an integral compatible system of bounded conductor and large geometric monodromy. Standard Mellin equidistribution would then predict continuous Sato–Tate behavior of normalized lifts. To deduce (5.2), one additionally needs one of:

1. a mod-\(p\) trace equidistribution theorem uniform in the coefficient characteristic;
2. an integral monodromy theorem giving near-surjectivity modulo \(p\);
3. a \(p\)-adic local limit theorem for the Mellin trace at tame characters;
4. a bounded-Weierstrass-degree interpolation over the entire tame character set.

None is presently available for this family.

## 5.2 Exponent thresholds

The currently proved all-value estimate is

\[
N_p(c)\le8p^{3/4}.
\]

Thus (5.2) is new for all fibres only if \(\delta>1/4\). The zero fibre has the stronger direct exponent \(2/3\), so improving it requires \(\delta>1/3\). The random model predicts \(\delta=1-o(1)\) in the sense that only \(O(1)\) twists vanish.

---

# 6. What the \(p\)-adic route explains structurally

## 6.1 It explains the special low-order zeros

At low character order, the Mellin family can acquire CM, self-twist, or Cartier structure. The quarter-point laws are exactly of this type: the character order is \(4\), and a class-field condition forces cancellation.

The Stickelberger ledger is compatible with this: low-order structure controls the residues of the unit terms, not their valuations.

## 6.2 Low-order character structure cannot build a bad moving column

This part is rigorous and useful. Fix an integer \(n\) and a window prime \(n/2<p\le n\). Put \(r=n-p\). Since \(p\equiv1\pmod{p-1}\),

\[
r=n-p\equiv n-1\pmod{p-1}.
\]

Therefore

\[
\operatorname{ord}(\omega_p^r)
=\frac{p-1}{\gcd(p-1,n-1)}.
\tag{6.1}
\]

If this order is at most \(T\), then

\[
p-1=ed,
\qquad e\le T,
\qquad d\mid n-1.
\]

For each pair \((e,d)\) there is at most one candidate \(p=ed+1\). Hence

\[
\boxed{
\#\{p\in(n/2,n]:\operatorname{ord}(\omega_p^{n-p})\le T\}
\le T\,\tau(n-1).
}
\tag{6.2}
\]

For polylogarithmic \(T\), this is \(n^{o(1)}\). Thus every bounded-order apparition law—including all fixed-order CM/Jacobsthal laws—is pointwise harmless.

## 6.3 It does not forbid generic high-order zeros

For high-order twists the minimum Stickelberger slope still occurs many times. There is no valuation gap. The condition is a congruence among \(\asymp p\) units.

Thus the strongest honest structural statement is:

- special low-order zeros are explainable and sparse in columns;
- generic high-order zeros are expected to be rare by random cancellation;
- the present \(p\)-adic theory does not make them impossible.

Indeed, the existence of nonempty \(Z_p\) for many primes, and the proved quarter-point zeros in the rank-two branches, demonstrate that exact unit cancellation is a real phenomenon rather than a forbidden one.

---

# 7. Comparison with the census fact \(F_4=0\)

The census statement \(F_4=0\) says that in the tested range no integer column is hit by four relevant primes. This is a statement about **simultaneous alignment across four different characteristics**.

The \(p\)-adic Mellin route is vertical: for one \(p\), it studies

\[
Z_p=\{r:\mathcal B_p(\omega^r)=0\}.
\]

Even the hypothetical bound \(|Z_p|=1\) for every prime does not prevent choosing those single residues to align at one integer \(n\). Therefore neither (5.2) nor the conjectural \(O(1)\) vertical bound explains \(F_4=0\) by itself.

What the \(p\)-adic route contributes is a potential **inverse theorem interface**:

> If many high-order Mellin zeros align in one column, prove that the corresponding integral Mellin objects acquire a bounded-order self-twist or another low-complexity degeneration.

Then (6.2) would kill the structured alternative. That is the “bad-diagonal inverse theorem” shape. It is strictly stronger than single-prime zero-density and is not supplied by Gross–Koblitz, Katz, or Wan.

So the census is evidence for high-order cross-prime pseudorandomness, not a consequence of the current \(p\)-adic analysis.

---

# 8. The strongest proved package from this route

The following statements are complete.

### Theorem 8.1 — exact one-character presentation

For \(1\le r\le p-2\), formula (1.12), equivalently (1.13), gives \(b_r\bmod p\) as one complete character sum with \(r\) entering only through \(\omega^{\pm r}\).

### Theorem 8.2 — exact Franel–Mellin presentation

For the same range, formula (1.11) gives \(b_r\) as the difference of two one-variable Mellin transforms. Its proof is only the pullback identity, the quadratic fibre count, and multiplicative orthogonality.

### Theorem 8.3 — exact Stickelberger filtration

Every Jacobi factor in (1.12) has valuation given by one carry; the product terms have slopes \(0,2,4\), apart from explicit trivial-character units.

### Theorem 8.4 — unit-cancellation criterion

With \(\rho=\min(r,p-1-r)\),

\[
p\mid b_r
\iff
\sum_{k=0}^{\rho}
\left(
\frac{\Gamma_p(\rho+k+1)}
{\Gamma_p(k+1)^2\Gamma_p(\rho-k+1)}
\right)^2
\equiv0\pmod p,
\]

and every summand is a unit.

### Theorem 8.5 — low-order twists are column-sparse

Equation (6.2) bounds the number of window primes carrying a twist of order at most \(T\) by \(T\tau(n-1)\).

What remains unproved is exactly [GAP-DCM], or the stronger cross-prime inverse theorem needed to explain \(F_4=0\).

---

# 9. A decisive next computation

The most informative experiment is not another complex absolute-value histogram. It is a **Stickelberger-stratified residue test**.

For each prime \(p\):

1. compute every term of (1.12) as a \(p\)-adic Jacobi/Gamma value;
2. record its carry slope \(0,2,4\);
3. sum only the slope-zero block and verify it equals \(b_r\bmod p\);
4. stratify zeros by
   \[
   \operatorname{ord}(\omega^r),
   \quad
   \gcd(r,p-1),
   \quad
   p\bmod24;
   \]
5. on the high-order subset, measure the full value histogram of \(\mathcal B_p(\omega^r)\in\mathbf F_p\), not merely the zero count.

The discriminating statistic is

\[
D_p^{\mathrm{hi}}=
\max_{a\in\mathbf F_p}
\left|\#\{r\in R_p^{\mathrm{hi}}:b_r=a\}
-\frac{|R_p^{\mathrm{hi}}|}{p}\right|.
\tag{9.1}
\]

A square-root-sized \(D_p^{\mathrm{hi}}\) would support the local-limit form of [GAP-DCM]; a persistent residue class or order stratum would reveal the missing self-twist. This test directly targets defining-characteristic distribution and cannot be replaced by complex Sato–Tate moments.

---

# Least-confident step

The least-confident substantive step is **not** the exact formulas or the Stickelberger ledger; those are elementary once the normalizations are fixed. It is the object-level assertion that the two amplitudes in the Franel–Mellin decomposition arise from one absolute bounded-conductor integral compatible system whose reduction at the coefficient prime is exactly the displayed \(\mathbf F_p\)-valued function.

The finite-field identity (1.11) is proved, and the Jacobi family (1.13) makes [GAP-DCM] unconditional as a statement. But applying Katz/Forey–Fresán–Kowalski machinery requires a precise integral sheaf/Cartier descent with controlled local monodromy. Until that is written down, the sheaf-theoretic route has two distinct gaps:

1. **[GAP-CARTIER]** construct and normalize the integral bounded-conductor Mellin object;
2. **[GAP-DCM]** prove defining-characteristic zero-density for its high-order twists.

Even granting the first, no cited theorem presently supplies the second.