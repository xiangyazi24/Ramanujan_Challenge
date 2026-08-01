PARTIAL — THE MAIN THEOREM IS NOT PROVED. SEVERAL LOAD-BEARING CLAIMS IN THE SPECIFICATION AND proof.tex ARE FALSE AS WRITTEN.

# Terminal verdict

The unconditional target

\[
G_n=\gcd(d_na_n,d_nb_n)=e^{o(n)},\qquad d_n=\operatorname{lcm}(1,\ldots,n)^3,
\]

remains open for every \(n\). The rigorous outcomes are:

1. The chart-free determinant reduction is valid after restoring a missing factor \(6\), using the physical nonwrapping strip, and distinguishing the union target \(U_B\) from the stronger incidence sum \(S_B\). A bound
   \[
   \frac1{p-1}\sum_{t\ne0}|B(t)|\ll p^{1-\kappa+o(1)}
   \]
   implies [FR_eta] for every fixed \(0<\eta<\min(1,2\kappa)\). Square-root cancellation allows every \(\eta<1/2\).
2. The claimed currently proved \(p\log p\) completion bound is not proved. It is conditional on the open complete-kernel estimate \(\max_{t\ne0,\xi}|F_t(\xi)|\ll p\). A positive convolution majorant removes the logarithm conditionally, but an \(O(p)\) error still gives no fixed \(\eta>0\) when \(L=N^{o(1)}\).
3. Exact second- and fourth-moment identities reduce §7.2 to centered determinant-value collision variance. Static point-set incidence/additive energy cannot supply it: a reflection-symmetric hyperbola clock has line intersections at most \(2\) and additive energy \(O(|A|^2)\), yet determinant mean \(\gg p\).
4. The master sum is exactly
   \[
   M(n)=\sum_{\sqrt n<p\le n}\log p\;\mathbf1_{p\mid b_{\,n\bmod p}}.
   \]
   The ledger formula using \(|R_p(n)|\) is false because it dropped \(p>r\). The top window is only quotient \(q=1\), not the full theorem. The middle contribution \((n,p,q,r)=(37,17,2,3)\) has \(17\mid G_{37}\).
5. [AVG-ZERO], \(\sum_{p\le X}|Z_p|=o(X)\), is equivalent up to \(p=5\) to the unnormalized Cesàro statement \(\sum_{n\le X}\mathrm{TOP}(n)=o(X)\). It is neither equivalent to nor sufficient for the every-\(n\) theorem, and sees only \(q=1\). No exponent improvement over \(X^{5/3}/\log X\) was obtained.
6. Contrary to the correction ledger's assertion that even an almost-all result is unknown, the existing codegree proof can be repaired. After replacing cap \(3\) by \(6\) and fixing companion height, it proves
   \[
   \#\{n\le X:\log G_n>\varepsilon n\}=O_\varepsilon((\log X)^2).
   \]
   This does not settle any exceptional \(n\).
7. A cross-prime reduction is proved below: an aggregate saving over the degree count for the growing family of gap polynomials would yield the first improvement for \(\sum_{p\le X}|Z_p|\). That uniform prime-aspect estimate remains open.

Every finite datum and exact identity used below is independently gated by CODEX_MAINTHM_verify.py; its final line is PASS.

# 1. Corrections before the main audit

## 1.1 The \(G_n\) determinant and valuation cap

Put \(A_n=d_na_n\), \(B_n=d_nb_n\). The Wronskian is

\[
a_nb_{n-1}-a_{n-1}b_n=\frac6{n^3}.
\]

The safe adjacent integer-vector determinant is

\[
A_{n+1}B_n-A_nB_{n+1}=\frac{6d_nd_{n+1}}{(n+1)^3}.
\]

Each term on the left contains a component from each vector, hence

\[
\boxed{G_nG_{n+1}\mid \frac{6d_nd_{n+1}}{(n+1)^3}.} \tag{1.1}
\]

Here \(d_n\) already is the cube of the lcm; the extra cubes in the specification are a notation mismatch. For \(p\ge5\), \(\sqrt n<p\le n\),

\[
v_p(G_n)\le6. \tag{1.2}
\]

The single-index integer determinant gives

\[
G_n\mid\frac{6d_n^2}{n^3}, \tag{1.3}
\]

which suffices for the small-prime \(O(\sqrt n)\) bound.

The proof.tex claim \(G_n\mid6d_n/n^3\) is false. At \(n=6\),

\[
G_6=270,\qquad 6d_6/6^3=6000,\qquad6000\equiv60\pmod{270}.
\]

Thus its cap \(v_p(G_n)\le3\) must be replaced by (1.2). This changes constants, not \(o(n)\) exponents.

## 1.2 Companion-height repair

Write \(a_q=A_q/C_q\) in lowest terms. The proof line

\[
|A_q|\le |a_q|C_q\le C_0(1+\sqrt2)^{4q}
\]

omits the denominator height. The measured \(\log|A_q|/q\) at \(q=10,50,100,200\) is

\[
5.213362,\ 6.314656,\ 6.200979,\ 6.463740,
\]

whereas \(\log(17+12\sqrt2)=3.525494\ldots\).

The needed result survives. Apéry integrality gives

\[
C_q\mid\operatorname{lcm}(1,\ldots,q)^3,\qquad \log C_q=O(q).
\]

Variation of parameters and the binomial formula give

\[
0<\frac{a_q}{b_q}=\sum_{k=1}^q\frac6{k^3b_kb_{k-1}}<6\zeta(3),
\qquad b_q\le(q+1)64^q.
\]

Therefore

\[
\log|\operatorname{num}(a_q)|=O(q). \tag{1.4}
\]

For \(q\le n^{1/3}\), summing logarithmic prime divisors costs \(O(n^{2/3})\). For \(q>n^{1/3}\), the associated prime has \(p<n/q<n^{2/3}\), and Chebyshev again gives \(O(n^{2/3})\).

## 1.3 Correct master reduction

For \(p>\sqrt n\), write \(n=qp+r\), \(0\le r<p\). The block congruence is

\[
p^3a_{qp+r}\equiv a_qb_r\pmod p. \tag{1.5}
\]

Since \(v_p(d_n)=3\),

\[
p\mid G_n\quad\Longleftrightarrow\quad a_qb_r\equiv0\pmod p. \tag{1.6}
\]

The \(p\mid\operatorname{num}(a_q)\) channel is covered by (1.4). Expanding the original radical gives

\[
\begin{aligned}
M(n)
&=\sum_{r<n/2}\ \sum_{\substack{p\mid b_r,\ p\mid n-r\\p>\max(\sqrt n,r)}}\log p\\
&=\boxed{\sum_{\sqrt n<p\le n}\log p\;\mathbf1_{p\mid b_{\,n\bmod p}}}. \tag{1.7}
\end{aligned}
\]

Conversely \(p>r=n\bmod p\) is automatic. Also \(r<n/2\): if \(p>n/2\), \(r=n-p\); otherwise \(r<p\le n/2\).

Equations (1.2), (1.4), and (1.6) give

\[
\log G_n\le6M(n)+O(n^{2/3}),\qquad M(n)\le\log G_n+O(1). \tag{1.8}
\]

Thus \(M(n)=o(n)\) is equivalent to the main theorem.

# 2. §7.1 — chart-free determinant reduction

## 2.1 Normalization and physical domain

Let \(p\ge5\), \(N=p-2\), \(u_r=(b_r,c_r)\), with \(c_0=0,c_1=6\). For continuants \(N_0=0,N_1=1\), Casoratian propagation gives

\[
\boxed{\det(u_r,u_{r+h})
=\frac{6N_h(r)}{\prod_{j=1}^h(r+j)^3}.} \tag{2.1}
\]

The specification omitted \(6\). Zeros and the all-\(t\ne0\) mean are unchanged because scaling permutes nonzero frequencies.

For

\[
J_D=\{d:\lfloor D/2\rfloor<d\le D\},\qquad
\Omega_D=\{(r,d):d\in J_D,\ 1\le r\le N-d\},
\]

put \(M_D=|\Omega_D|\). Exactly,

\[
M_D=
\begin{cases}
mN-m(3m+1)/2,&D=2m,\\
(m+1)N-(m+1)(3m+2)/2,&D=2m+1.
\end{cases} \tag{2.2}
\]

Thus \(M_D=(1/2+o(1))ND\) for \(D=o(N)\). The old script used \(r\le N-D\), omitting \(q(q-1)/2\) pairs, \(q=\lceil D/2\rceil\).

## 2.2 Union versus incidence

Let

\[
Z_d=\{r:(r,d)\in\Omega_D,\ \det(u_r,u_{r+d})=0\}.
\]

The actual [FR_eta] target is \(U_B=|\bigcup_{d\in J_D}Z_d|\), whereas orthogonality computes \(S_B=\sum_d|Z_d|\). Only

\[
U_B\le S_B \tag{2.3}
\]

holds. The verifier gives a strict live example \((p,D,U_B,S_B)=(11,6,3,4)\). The determinant route is sufficient but strictly stronger.

Define

\[
\nu(x)=\#\{(r,d)\in\Omega_D:\det(u_r,u_{r+d})=x\},\qquad
B(t)=\sum_x\nu(x)e_p(tx).
\]

Orthogonality gives

\[
S_B=\frac{M_D}{p}+\frac1p\sum_{t\ne0}B(t). \tag{2.4}
\]

With \(\mathcal A_1=(p-1)^{-1}\sum_{t\ne0}|B(t)|\),

\[
U_B\le S_B\le\frac{M_D}{p}+\frac{p-1}{p}\mathcal A_1. \tag{2.5}
\]

The main term is asymptotic to \(D/2\).

## 2.3 Exponent budget

The intended regime is

\[
D=\sqrt N\,L,\qquad L\to\infty,\quad L=N^{o(1)},\quad D=o(N). \tag{2.6}
\]

The target is \(D^{2-\eta}=N^{1-\eta/2}L^{2-\eta}\). A square-root strip bound is \(\sqrt{ND}=N^{3/4}L^{1/2}\), with ratio

\[
N^{-1/4+\eta/2}L^{-3/2+\eta}. \tag{2.7}
\]

This tends to zero for every fixed \(\eta<1/2\). At \(\eta=1/2\), a constant-loss bound also wins because \(L\to\infty\); a \(p^{o(1)}\) loss safely gives strict \(\eta<1/2\).

If

\[
\mathcal A_1\ll p^{1-\kappa+o(1)}, \tag{2.8}
\]

then the error/target ratio is \(N^{\eta/2-\kappa+o(1)}L^{\eta-2}\). Including the main term, (2.8) proves the stronger incidence bound and [FR_eta] for

\[
0<\eta<\min(1,2\kappa). \tag{2.9}
\]

Thus any fixed power saving closes some [FR_eta]; a logarithmic saving does not.

## 2.4 Completion: proved versus conditional

Extend the table cyclically and set

\[
A_t(d)=\sum_r e_p(t\det(u_r,u_{r+d})),\qquad
F_t(\xi)=\sum_dA_t(d)e_p(\xi d).
\]

For nonnegative \(W\ge\mathbf1_{J_D}\),

\[
B_W(t)=\frac1p\sum_\xi\widehat W(-\xi)F_t(\xi). \tag{2.10}
\]

The sharp interval has \(\sum|\widehat W|=O(p\log p)\), so

\[
\max_{t\ne0,\xi}|F_t(\xi)|\ll p
\Longrightarrow |B_W(t)|\ll p\log p. \tag{2.11}
\]

The premise is open; the repository has experiments, not a bounded-conductor theorem. The specification's “currently proved” claim is refuted.

There is a useful proved smoothing brick. Let \(J=[a,a+q-1]\), \(C=[0,q-1]\),

\[
A=[a-q+1,a+q-1],\qquad W=q^{-1}\mathbf1_A*\mathbf1_C.
\]

If \(3q-2<p\),

\[
0\le W\le1,\qquad W=1\text{ on }J,\qquad\sum W=2q-1. \tag{2.12}
\]

Parseval and Cauchy–Schwarz give

\[
\sum_\xi|\widehat W(\xi)|
\le p\sqrt{\frac{2q-1}{q}}<\sqrt2\,p. \tag{2.13}
\]

Thus the open complete bound would improve (2.11) to \(O(p)\). But \(p/D^{2-\eta}=N^{\eta/2}L^{\eta-2}\) still diverges for fixed \(\eta>0\). A complete \(p^{1-\kappa}\) bound, or equivalent average, would close (2.9).

The tempting \(\sum W^2\le(2q-1)/q\) is false: \(W=1\) at \(q\) points already gives \(\sum W^2\ge q\).

## 2.5 Corrected full-frequency experiment

The old \(148,195,462\) values used \(t=1,\ldots,59\) and a rectangle. The complete natural-strip data are:

| \(p\) | \(D\) | pairs | zeros | mean \(|B(t)|\) | \(\sqrt{M_D}\) |
|---:|---:|---:|---:|---:|---:|
| 1009 | 62 | 29,760 | 46 | 141.437464 | 172.510869 |
| 2003 | 88 | 85,118 | 68 | 227.222675 | 291.749893 |
| 4003 | 126 | 246,078 | 88 | 400.527151 | 496.062496 |
| 8009 | 178 | 700,697 | 127 | 676.916906 | 837.076460 |
| 16001 | 252 | 1,991,997 | 179 | 1125.078760 | 1411.381238 |

This extends the experiment to \(p=16001\). It is square-root-compatible but EMPIRICAL.

# 3. §7.2 — moments, geometry, and sheaves

## 3.1 Exact second and fourth moments

For the determinant histogram \(\nu\), put \(M=\sum_x\nu(x)\), \(\rho(x)=\nu(x)-M/p\). Orthogonality gives

\[
\sum_{t\ne0}|B(t)|^2=p\sum_x\nu(x)^2-M^2=p\|\rho\|_2^2. \tag{3.1}
\]

Hence

\[
\mathcal A_1\le\sqrt{\frac p{p-1}}\,\|\rho\|_2. \tag{3.2}
\]

The useful missing estimate is

\[
\sum_x\nu(x)^2=\frac{M^2}{p}+O(Mp^{o(1)}). \tag{3.3}
\]

It would give \(\mathcal A_1\ll\sqrt M\,p^{o(1)}=p^{3/4+o(1)}\) and close every fixed \(\eta<1/2\).

Let

\[
E_+(\nu)=\sum_s\left(\sum_x\nu(x)\nu(s-x)\right)^2.
\]

Then exactly

\[
\sum_{t\ne0}|B(t)|^4
=pE_+(\nu)-M^4
=p\left(E_+(\nu)-\frac{M^4}{p}\right). \tag{3.4}
\]

A fourth-moment proof of \(\mathcal A_1\ll p^{1-\kappa}\) needs

\[
E_+(\nu)-M^4/p\ll p^{4-4\kappa}. \tag{3.5}
\]

This is energy of the clock-ordered determinant-value multiset, not additive energy of the static orbit support.

## 3.2 Why gapwise Cauchy–Schwarz stalls

Let \(\nu_d\) be the histogram for one gap, \(M_d=\sum\nu_d\), \(\rho_d=\nu_d-M_d/p\). Then

\[
\rho=\sum_d\rho_d,\qquad
\|\rho\|_2^2=\sum_d\|\rho_d\|_2^2
+2\sum_{d<e}\langle\rho_d,\rho_e\rangle. \tag{3.6}
\]

Cauchy–Schwarz before controlling covariances gives

\[
\mathcal A_1
\le\left(\frac{p|J_D|}{p-1}\sum_d\|\rho_d\|_2^2\right)^{1/2}. \tag{3.7}
\]

At the random fixed-gap scale \(\|\rho_d\|_2^2\asymp p\), this is order \(p\), with no power saving. The missing statement is clockwise almost orthogonality,

\[
\left\|\sum_d\rho_d\right\|_2^2
\ll p^{o(1)}\sum_d\|\rho_d\|_2^2. \tag{3.8}
\]

No such theorem is proved.

## 3.3 Static geometry cannot imply (3.8)

Let \(p=2m+1\), choose primitive \(g\), and put

\[
j(r)=\min(r,2m-r),\qquad v_r=(g^{j(r)},g^{-j(r)}).
\]

The support lies on \(xy=1\), has reflection multiplicity two except at the center, and:

- every affine line meets it in at most two points, since substitution gives a quadratic;
- \(E_+(A)\le3|A|^2\), since a nonzero vector sum determines the unordered pair from sum and product, while zero-sum pairs form one class.

Nevertheless, on each noncrossing half,

\[
\det(v_r,v_{r+d})=g^{-d}-g^d
\]

is independent of \(r\). Parseval and \(\|B\|_\infty\le M\) give

\[
\operatorname{mean}_{t\ne0}|B(t)|\gg p
\]

for \(D=o(p)\) in a fixed positive range. The verifier checks \(p=101,211,401\), including a certified linear lower bound. Static Rudnev/Stevens–de Zeeuw/Shkredov-style support estimates alone cannot prove clock-ordered cancellation.

## 3.4 Vector Weyl sums

For

\[
S(\alpha,\beta)=\sum_{r=1}^Ne_p(\alpha b_r+\beta c_r),
\]

let \(\mu(v)\) be affine-vector multiplicity. Exact Parseval is

\[
\sum_{\alpha,\beta}|S(\alpha,\beta)|^2=p^2\sum_v\mu(v)^2. \tag{3.9}
\]

Reflection gives \(\sum\mu^2\ge2N-1\), not \(2N\), because of one fixed midpoint. Extra coincidences occur:

\[
p=73:\quad (b_r,c_r)=(5,6)\text{ at }r=1,4,68,71,
\]

and \(p=997\) also has a fourfold collision. Exactly 17 of the 166 primes \(5\le p\le1000\) have extra collisions. Thus “exactly two-to-one for every prime” is false.

Full spectra:

| \(p\) | \(\sum\mu^2\) | max nonzero \(|S|\) | ratio to \(\sqrt N\) |
|---:|---:|---:|---:|
| 1009 | 2013 | 176.033819 | 5.54730 |
| 2003 | 4001 | 244.481930 | 5.46542 |
| 3001 | 6005 | 294.307009 | 5.37418 |

The \(p=3001\) line extends the original experiment; all three lines are EMPIRICAL.

Completion actually needs

\[
S_\xi(\alpha,\beta)=\sum_u e_p(\alpha b_u+\beta c_u+\xi u),
\]

with exact identity

\[
F_t(\xi)=\sum_r e_p(-\xi r)S_\xi(-tc_r,tb_r). \tag{3.10}
\]

The measured untwisted \(S_0\) covers only \(\xi=0\). Even a hypothetical square-root bound for every \(S_\xi\), inserted by naive outer Cauchy/triangle, gives \(F_t(\xi)\ll p^{3/2}\), not \(p\). A moving-frequency correlation is missing.

## 3.5 Reflection and the sheaf claim

Strong reflection \(u_{p-1-r}=u_r\) maps

\[
(r,d)\longmapsto(p-1-r-d,d)
\]

inside the physical strip and negates the determinant. Thus

\[
\nu(x)=\nu(-x),\qquad B(t)\in\mathbf R. \tag{3.11}
\]

There is one fixed zero edge for each even \(d\), at \(2r+d=p-1\).

No fixed-conductor l-adic realization of the two-variable kernel has been constructed. Characteristic-zero nonvanishing cannot be transferred to the working prime without saturation:

\[
N_3(-3)=584=8\cdot73,\qquad
\gcd_{\mathbf F_{73}}(N_3,N_4)=X+3,
\]

and

\[
\gcd_{\mathbf F_{211}}(N_{32},N_{32}')=(X-89)^2.
\]

Thus honest-pole and characteristic-zero resultant arguments do not supply a fixed-\(p\) sheaf theorem in the mesoscopic range. The bounded-conductor route remains OPEN.

# 4. §7.3 — master sum and cross-prime routes

## 4.1 The incidence \(R_p(n)\) formula is false

The correction ledger defines

\[
R_p(n)=\{r<n/2:r\equiv n\pmod p,\ p\mid b_r\}
\]

and claims \(M(n)=\sum\log p\,|R_p(n)|\). This ignores \(p>r\). At

\[
n=101,\qquad p=31,
\]

both \(r=8,39\) lie in \(R_{31}(101)\), but only \(r=8<p\) contributes to \(M(101)\). Formula (1.7) is the corrected exact statement.

## 4.2 Digit criterion and top window

Gessel's congruence

\[
b_{mp+r}\equiv b_mb_r\pmod p
\]

iterates over the base-\(p\) expansion \(n=\sum n_ip^i\) to give

\[
p\mid b_n
\Longleftrightarrow
\text{some base-}p\text{ digit }n_i\text{ belongs to }Z_p. \tag{4.1}
\]

For \(p>n/2\), \(n=p+(n-p)\), so for \(p>5\),

\[
b_n\equiv5b_{n-p}\pmod p,\qquad
p\mid b_n\Longleftrightarrow p\mid b_{n-p}. \tag{4.2}
\]

This is the \(q=1\) channel, not a new Frobenius invariant and not the full theorem. The middle witness

\[
37=2\cdot17+3,\qquad17\mid b_3,\qquad17\mid G_{37}
\]

is absent from the top window.

The exact pointwise criterion is

\[
\sum_{\sqrt n<p\le n}\log p\;\mathbf1_{\,n\bmod p\in Z_p}=o(n), \tag{4.3}
\]

equivalently

\[
\#\{\sqrt n<p\le n:n\bmod p\in Z_p\}=o(n/\log n). \tag{4.4}
\]

The verifier constructs a reflection-symmetric row model aligning 60 quotient-\(2\) primes at \(N=10000\), while every induced quotient-\(1\) target has load at most one. It is not an Apéry counterexample; it shows row cardinality, reflection, and Lucas reindexing alone cannot promote a top-window theorem to (4.3).

## 4.3 Correct [AVG-ZERO] identity

Let

\[
\mathrm{TOP}(n)=\#\{p:n/2<p\le n,\ p\mid b_n\},\qquad
S(X)=\sum_{p\le X}|Z_p|.
\]

For \(p>5\), (4.2) and \(n=p+r\) give

\[
\sum_{n\le X}\mathrm{TOP}_{p>5}(n)
=R(X):=\#\{(p,r):p>5,\ r\in Z_p,\ p+r\le X\}. \tag{4.5}
\]

The supplied identity_check.py did not check the left side: it tested \(p\mid b_{n-p}\). At \(p=5\), \(b_1\equiv0\), so division by \(5\) is invalid. In fact \(Z_5=\{1,3\}\), also refuting the proof.tex claim \(Z(5)=0\). The endpoint correction is

\[
\sum_{n\le X}\mathrm{TOP}(n)
=R(X)+\mathbf1_{X\ge5}+\mathbf1_{X\ge7}+\mathbf1_{X\ge9}. \tag{4.6}
\]

At \(X=600\),

\[
\sum_{n\le600}\mathrm{TOP}(n)=85,\quad
R(600)=82,\quad R_{p>5}(600)=80,\quad S(600)=109.
\]

The truncation matters. Exactly,

\[
R(X)\le S(X)\le R(2X), \tag{4.7}
\]

up to the finite-prime convention. Therefore

\[
\frac1X\sum_{n\le X}\mathrm{TOP}(n)\to0
\Longleftrightarrow
\boxed{S(X)=\sum_{p\le X}|Z_p|=o(X)}. \tag{4.8}
\]

This is mean top-hit count tending to zero, not “the average version of the main theorem” without qualification. It does not imply the every-\(n\) maximum in (4.3): sparse spikes can total \(o(X)\). Conversely \(\mathrm{TOP}(n)=1\) at every \(n\) satisfies the pointwise main-theorem scale but sums to \(\asymp X\). The two norms are incomparable, and (4.8) sees no \(q\ge2\).

The measured \(S(600)=109=\pi(600)\) and mean \(|Z_p|=0.912\) below \(4000\) support \(S(X)\asymp X/\log X\), but are not theorems. Unconditionally,

\[
S(X)\ll\frac{X^{5/3}}{\log X}. \tag{4.9}
\]

This already gives the usual density-one normalized top-window conclusion by Markov; [AVG-ZERO] would give the much stronger density-one assertion \(\mathrm{TOP}(n)=0\).
