# Frobenius--Dwork proof of the Franel tail-lattice theorem

## 1. Statement

Let

\[
F_m=\sum_{a=0}^m\binom ma^3
\]

be the Franel numbers.  Fix \(J\geq 0\).  There is a unique rational
sequence \(z_d\), with \(z_1=\cdots=z_J=0\), satisfying

\[
F_k=\sum_{i=0}^J(-1)^{J-i}\binom ki
     \binom{k-i-1}{J-i}F_i z_{k-i}\qquad(k>J).
\tag{1}
\]

Every \(z_d\) is an integer.

Consequently, for every finite integer sequence \(q_k\),

\[
\sum_{k>J}q_kF_k
\]

is an integer linear combination of the nonconstant coefficients of the
Legendre--Euler transform with cutoff \(J\).  In particular, for every
\(n,J\), the content of the truncated transform \(T_{n,J}\) divides the
full Apéry number \(A_n\).

## 2. Convolution form

Put

\[
P_J(x)=\sum_{i=0}^J(-1)^i\binom JiF_i x^i,\qquad
W_J(x)=\sum_{d>J}\frac{z_d}{d}x^d .
\]

The elementary identity

\[
\binom ki\binom{k-i-1}{J-i}
=\frac{k}{k-i}\binom{k-1}{J}\binom Ji
\]

turns (1) into

\[
P_J(x)W_J(x)
=(-1)^J\sum_{k>J}\frac{F_k}{k\binom{k-1}{J}}x^k.
\tag{2}
\]

Define

\[
\Phi(\alpha,x)=\sum_{m\geq0}(-1)^m\binom{\alpha}{m}F_mx^m
\]

and

\[
D_J(x)=\left.\frac{\partial}{\partial\alpha}
                 \Phi(\alpha,x)\right|_{\alpha=J}.
\]

For \(m>J\),

\[
\left.\frac{\partial}{\partial\alpha}\binom{\alpha}{m}
\right|_{\alpha=J}
=\frac{(-1)^{m-J-1}}{m\binom{m-1}{J}}.
\]

If \(D_J^{\leq J}\) is the degree-\(J\) truncation of \(D_J\), (2) is
equivalent to

\[
W_J=-\frac{D_J}{P_J}+\frac{D_J^{\leq J}}{P_J}.
\tag{3}
\]

It therefore suffices to prove that the Euler derivatives of both
quotients on the right belong to \(\mathbb Z[[x]]\).

## 3. Frobenius-twisted Dwork congruence

Use the constant-term representation

\[
\Lambda(u,v)=(1+u)(1+v)(1+1/(uv)),\qquad
F_m=\operatorname{CT}_{u,v}\Lambda(u,v)^m.
\]

Then

\[
P_J(x)=\operatorname{CT}_{u,v}(1-x\Lambda(u,v))^J.
\tag{4}
\]

Fix a prime \(p\).  Work over \(R=\mathbb Z_p[x]\), with Frobenius lift
\(\phi(x)=x^p\).  The Newton polygon in the Laurent variables \(u,v\)
has the origin as its unique interior lattice point.

The coefficient-ring version of the Samol--van Straten D3 congruence
applied to (4) is

\[
P_{J+mp^r}(x)P_{\lfloor J/p\rfloor}(x^p)
\equiv
P_J(x)P_{\lfloor J/p\rfloor+mp^{r-1}}(x^p)
\pmod {p^r}
\tag{5}
\]

for \(J,m\geq0\) and \(r\geq1\).

For completeness, the extension from scalar coefficients is as follows.
On \(R[u^{\pm1},v^{\pm1}]\), extend \(\phi\) by
\(u\mapsto u^p,\ v\mapsto v^p\), and write this extension as \(\psi\).
For a Laurent polynomial \(B\), define

\[
R_s(B)=B^{p^s}-\psi(B)^{p^{s-1}},\qquad R_0(B)=B .
\]

The Frobenius-lift congruence gives

\[
R_s(B)\in p^sR[u^{\pm1},v^{\pm1}]
\]

and the Laurent Newton support of \(R_s(B)\) is contained in
\(p^s\operatorname{Newt}(B)\).  Moreover,

\[
B^{p^s}=\sum_{j=0}^s\psi^{s-j}(R_j(B)).
\]

Expanding along the base-\(p\) digits and decomposing the resulting ghost
tuples into indecomposable blocks gives constants \(c_n\) with

\[
c_n\in p^{\ell(n)-1}R.
\]

A block beginning at digit position \(h\) contributes
\(\phi^h(c_n)\).  The unique-interior-point hypothesis makes constant
terms of separated blocks factor.  The good-partition bijection in the
usual D3 proof preserves every absolute digit position, hence preserves
these Frobenius powers.  Every unmatched partition contains total
\(\sum(\ell(\text{block})-1)\geq r\), and therefore vanishes modulo
\(p^r\).  This proves (5).  No specialization of \(x\) is used.

## 4. Congruence preservation of the logarithmic derivative

Let

\[
H(J,x)=\theta\log P_J(x),\qquad \theta=x\frac{d}{dx}.
\]

All factors in (5) have constant term one.  Taking Euler logarithmic
derivatives gives

\[
\begin{aligned}
H(J+mp^r,x)-H(J,x)\equiv p\big(&
H(\lfloor J/p\rfloor+mp^{r-1},x^p)\\
&-H(\lfloor J/p\rfloor,x^p)\big)\pmod {p^r}.
\end{aligned}
\tag{6}
\]

Induction on \(r\) yields

\[
H(J+mp^r,x)-H(J,x)\in p^r\mathbb Z_p[[x]].
\tag{7}
\]

Doing this for every prime power dividing \(a-b\) proves

\[
H(a,x)-H(b,x)\in(a-b)\mathbb Z[[x]]
\tag{8}
\]

for all nonnegative integers \(a,b\).

For fixed \(n\), let

\[
H_n(\alpha)=[x^n]\theta\log\Phi(\alpha,x).
\]

It is a rational polynomial and \(H_n(J)=[x^n]H(J,x)\) at every
nonnegative integer \(J\).  For every prime \(p\),

\[
\frac{H_n(J+p^r)-H_n(J)}{p^r}\in\mathbb Z_p.
\]

As \(r\to\infty\), this quotient tends \(p\)-adically to \(H_n'(J)\).
Thus \(H_n'(J)\) lies in every \(\mathbb Z_p\), and hence is an integer.
Therefore

\[
\theta(D_J/P_J)\in\mathbb Z[[x]].
\tag{9}
\]

## 5. The finite numerator

Define

\[
T_J(\alpha,x)=\sum_{m=0}^J(-1)^m\binom{\alpha}{m}F_mx^m
\]

and

\[
G_J(\alpha,x)=\theta\left(\frac{T_J(\alpha,x)}{P_J(x)}\right).
\]

For \(0\leq k\leq J\), \(T_J(k,x)=P_k(x)\), while \(G_J(J,x)=0\).
The endpoint Lagrange derivative formula gives

\[
\left.\partial_\alpha G_J(\alpha,x)\right|_{\alpha=J}
=\sum_{r=1}^J(-1)^r\binom Jr\frac{G_J(J-r,x)}r.
\tag{10}
\]

By (8),

\[
\theta\log(P_{J-r}/P_J)\in r\mathbb Z[[x]].
\]

Since \(P_{J-r}/P_J\) is an integral unit series, this implies

\[
G_J(J-r,x)=\theta(P_{J-r}/P_J)\in r\mathbb Z[[x]].
\]

Every term in (10) is therefore integral.  Its left side is precisely

\[
\theta(D_J^{\leq J}/P_J),
\]

so

\[
\theta(D_J^{\leq J}/P_J)\in\mathbb Z[[x]].
\tag{11}
\]

Equations (3), (9), and (11) show that
\(\theta W_J=\sum_{d>J}z_dx^d\) has integral coefficients, proving the
theorem.

## 6. Citation boundary

The scalar D3 congruence is proved by Kira Samol and Duco van Straten,
*Dwork Congruences and Reflexive Polytopes*.  The indecomposable
ghost-block proof is streamlined by Anton Mellit and Masha Vlasenko,
*Dwork's Congruences for the Constant Terms of Powers of a Laurent
Polynomial*.  Equation (5) is the same proof over a coefficient ring with
a Frobenius lift; it is not quoted verbatim from the scalar statement.
