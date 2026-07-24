# A degree-one signed Newton--Padé no-go theorem

## 1. Setup

Let \(A_m\) be the Apéry numbers and put

\[
c_k=\Delta^k A_0.
\]

Thus the Newton expansion of the integer-valued function \(A_x\) is

\[
A_x=\sum_{k\geq 0}c_k\binom{x}{k}.
\]

Let

\[
H=\left\lfloor\frac{n-1}{3}\right\rfloor.
\]

The degree-one signed Newton--Padé problem asks for

\[
Q(x)=q_0+q_1x,\qquad \deg P\leq H-1,
\]

such that

\[
P(j)=A_jQ(j)\qquad(0\leq j\leq H).
\tag{1.1}
\]

There is a unique primitive integral pair \((P,Q)\), up to a common sign.

## 2. Recurrence and asymptotics of the Newton coefficients

If

\[
\mathcal A(z)=\sum_{m\geq0}A_mz^m,\qquad
\mathcal C(w)=\sum_{k\geq0}c_kw^k,
\]

then binomial inversion gives

\[
\mathcal C(w)=\frac{1}{1+w}
\mathcal A\left(\frac{w}{1+w}\right).
\tag{2.1}
\]

Substitution of (2.1) into the Apéry differential equation gives, for
\(k\geq0\),

\[
\begin{aligned}
(k+4)^3c_{k+4}
={}&(2k+7)(15k^2+105k+184)c_{k+3}\\
&+(k+3)(95k^2+570k+864)c_{k+2}\\
&+48(k+2)(k+3)(2k+5)c_{k+1}\\
&+32(k+1)(k+2)(k+3)c_k.
\end{aligned}
\tag{2.2}
\]

The characteristic polynomial at infinity is

\[
\lambda^4-30\lambda^3-95\lambda^2-96\lambda-32
=(\lambda+1)^2(\lambda^2-32\lambda-32).
\tag{2.3}
\]

Its unique dominant root is

\[
\gamma=16+12\sqrt2.
\]

There is a direct way to justify the required asymptotic without invoking an
unspecified Poincaré theorem.  Edgar's Stieltjes-moment representation for
the Apéry numbers gives

\[
A_m=\int_0^{c_*}x^m\phi(x)\,dx,\qquad
c_*=17+12\sqrt2,
\tag{2.4}
\]

where \(\phi\) is positive and integrable.  At the right endpoint its
regular-singular Frobenius expansion is

\[
\phi(c_*-t)=C\,t^{1/2}(1+\beta t+O(t^2)),\qquad C>0.
\tag{2.5}
\]

See G. A. Edgar, [*The Apéry Numbers as a Stieltjes Moment
Sequence*](https://arxiv.org/abs/2005.10733), especially Proposition 25 and
Corollary 27.

Finite differences may be taken under the integral:

\[
c_k=\int_0^{c_*}(x-1)^k\phi(x)\,dx.
\tag{2.6}
\]

The part \(0\leq x\leq1\) is \(O(1)\).  Watson's endpoint lemma applied to
the remaining integral and (2.5) gives

\[
c_k=\kappa\gamma^k k^{-3/2}
\left(1+\frac{\alpha}{k}+O(k^{-2})\right),
\qquad \kappa>0.
\tag{2.7}
\]

In particular, with

\[
\rho_k=\frac{c_{k-1}}{c_k},\qquad
f(k)=k(1+\rho_k),
\]

one has

\[
\rho_k=\gamma^{-1}
\left(1+\frac{3}{2k}+O(k^{-2})\right)
\tag{2.8}
\]

and hence

\[
f(k)=\left(1+\gamma^{-1}\right)k
+\frac{3}{2\gamma}+O(k^{-1}).
\tag{2.9}
\]

Consequently \(f(k+1)-f(k)=1+\gamma^{-1}+O(k^{-1})>0\) for all
sufficiently large \(k\). Since \(f(H)\to\infty\), it follows that, for all
sufficiently large \(H\),

\[
f(H)>f(k)\qquad(0\leq k<H),
\tag{2.10}
\]

where \(f(0)=0\).

## 3. Exact solution of the degree-one Padé equations

The elementary identity

\[
x\binom{x}{k}
=k\binom{x}{k}+(k+1)\binom{x}{k+1}
\]

shows that the \(k\)-th Newton coefficient of \(A_xQ(x)\) is

\[
d_k=(q_0+kq_1)c_k+kq_1c_{k-1},
\tag{3.1}
\]

with \(c_{-1}=0\).

Condition (1.1) with \(\deg P\leq H-1\) is equivalent to \(d_H=0\).
Thus, writing

\[
g=\gcd\left(H(c_H+c_{H-1}),c_H\right),
\]

the primitive solution can be chosen as

\[
q_0=\frac{H(c_H+c_{H-1})}{g},
\qquad
q_1=-\frac{c_H}{g}.
\tag{3.2}
\]

Equations (3.1)--(3.2) give the exact factorization

\[
\boxed{
d_k=\frac{c_Hc_k}{g}\bigl(f(H)-f(k)\bigr)
}
\qquad(0\leq k<H).
\tag{3.3}
\]

By (2.7), every \(d_k\) is positive for all sufficiently large \(H\).
Therefore

\[
P(x)=\sum_{k=0}^{H-1}d_k\binom{x}{k}
\tag{3.4}
\]

has positive integral Newton coefficients.

## 4. Exponential lower bound

Take \(k=\lfloor H/2\rfloor\). Since \(d_k\) is a positive integer,

\[
|P(n)|=P(n)\geq\binom{n}{\lfloor H/2\rfloor}.
\]

As \(H/n\to1/3\), Stirling's formula yields

\[
\liminf_{n\to\infty}\frac1n\log|P(n)|
\geq
-\frac16\log\frac16-\frac56\log\frac56
=0.450561\ldots.
\tag{4.1}
\]

Hence the unique primitive signed Newton--Padé pair with
\(\deg Q=1\) cannot provide a subexponential direct-\(q=1\) certificate.

## 5. Scope

This closes the first genuinely signed denominator case left open by Q729.
It does not treat \(\deg Q\geq2\), where the high Newton coefficients impose
a higher-order signed kernel and the one-dimensional monotonicity function
\(f(k)\) is replaced by a determinant/minor system.
