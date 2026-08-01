# The Apéry rank-two apparition law: proof, corrections, and remaining gaps

## 0. Result ledger

Write

\[
 b_n=\sum_{k=0}^n\binom nk^2\binom{n+k}k^2,
 \qquad F(t)=\sum_{n\geq0}b_nt^n,
 \qquad q(t)=1-34t+t^2,
\]

and choose the formal branches

\[
 \tau(t)=\sqrt{F(t)}=\sum_{n\geq0}\tau_nt^n,
 \qquad
 \sigma(t)=\sqrt{F(t)/q(t)}=\sum_{n\geq0}\sigma_nt^n
\]

with constant coefficient one.  For an odd prime $p$, put

\[
 A_p(t)=\sum_{n=0}^{p-1}b_nt^n,
 \qquad \chi_p=\left(\frac{-6}{p}\right).
\]

The outcome of this investigation is as follows.

1. **[PROVED] Exact relevant-branch quarter law.**  If \(\chi_p=1\), use the
   normalized polynomial \(s_p=[\tau]_{\leq(p-1)/2}\); if \(\chi_p=-1\), use
   \(s_p=[\sigma]_{\leq(p-3)/2}\).  Its reversal law is

   \[
   t^{\deg s_p}s_p(1/t)=\left(\frac{-2}{p}\right)s_p(t).
   \]

   Consequently the selected floor-quarter coefficient vanishes exactly for
   $p\equiv5,23\pmod {24}$: it is
   \(\tau_{(p-1)/4}\) in class $5$, and \(\sigma_{(p-3)/4}\) in class
   $23$.  The selected coefficients (or the middle pair when the degree is
   odd) are nonzero in the other six classes.  This closes `[GAP-BR]` from the
   working notes; no endpoint Gauss-sum evaluation is needed.

2. **[PROVED] Hypergeometric closed forms.**  The two series are algebraic
   pullbacks of one classical \({}_2F_1\), displayed in Section 1.  Exact
   finite coefficient formulae follow from Lagrange inversion.

3. **[DISPROVED] The proposed discriminant $-24$ parametrization of both
   vanishing classes.**  The form $2x^2+3y^2$ represents primes
   $p\equiv5,11\pmod {24}$, not $p\equiv23\pmod {24}$.  Primes in class
   $23$ are inert in \(\mathbf Q(\sqrt{-6})\) and have no representation by
   either reduced form of discriminant $-24$.  The quarter law is governed
   by the branch character \((-6/p)\), degree parity, and the reversal
   character \((-2/p)\), rather than by one pair \((x,y)\).

4. **[PROVED] Convolution shadow.**  Membership in $Z_p$ is an explicit
   quadratic convolution condition on the complete values of \(\tau\) or
   \(\sigma\).  **[DISPROVED]** It is not determined by their zero supports.
   No lattice-point membership law for $Z_p$ follows from the quarter zero.

5. **[PROVED] Bernoulli formula.**  For every prime $p\geq5$,

   \[
   b_p\equiv5-\frac{14}{3}p^3B_{p-3}\pmod {p^4},\qquad
   \beta_p:=\frac{b_p-5}{p^3}\equiv-\frac{14}{3}B_{p-3}\pmod p.
   \]

   This is the $p$-adic, order-three lift of the same convolution, but it is
   a Bernoulli--Wolstenholme invariant rather than a discriminant-$-24$
   coordinate.  A direct Gross--Koblitz/Jacobi-sum reformulation remains
   `[GAP-GK]`.

All machine claims below are reproduced by the three `codex_max_*.py` scripts
listed in Section 6.

## 1. The actual hypergeometric form of the two square roots

### 1.1 Differential equations and recurrences

The symmetric-square descent of the Apéry differential equation gives

\[
 (t^3-34t^2+t)y''+(2t^2-51t+1)y'+\frac{t-10}{4}y=0
\]

for \(y=\tau\).  Equivalently,

\[
 4(n+1)^2\tau_{n+1}
 =2(68n^2+34n+5)\tau_n-(2n-1)^2\tau_{n-1},
\]

with \(\tau_0=1,\tau_1=5/2\).  The companion branch satisfies

\[
 4(n+1)^2\sigma_{n+1}
 =2(68n^2+102n+39)\sigma_n-(2n+1)^2\sigma_{n-1},
\]

with \(\sigma_0=1,\sigma_1=39/2\).  The first equation and recurrence appear
explicitly in Frits Beukers, *Consequences of Apéry's work on \(\zeta(3)\)*,
pp. 1--2; the sigma equation follows by substituting
\(\tau=q^{1/2}\sigma\).  Beukers also records that the singularities of the
rank-two equation are $0,(1\pm\sqrt2)^4,\infty$.  Thus this is a Heun-type
equation in $t$, not a Gauss equation in the coefficient variable.

### 1.2 Franel/\({}_2F_1\) pullback

Let

\[
 f_m=\sum_{k=0}^m\binom mk^3,\qquad h(x)=\sum_{m\geq0}f_mx^m,
 \qquad t=\frac{x(1-8x)}{1+x}.
\]

Caruso--Fürnsinn--Vargas-Montoya--Zudilin use the exact identity

\[
 F(t(x))=(1+x)h(x)^2.
\]

The Franel series has the classical pullback

\[
 h(x)=\frac1{1-2x}
 {}_2F_1\!\left(\begin{matrix}1/3,2/3\\1\end{matrix};
                 \frac{27x^2}{(1-2x)^3}\right).
\]

Since direct algebra gives

\[
 q(t(x))=\left(\frac{1-16x-8x^2}{1+x}\right)^2,
\]

the normalized formal branches have the promised closed forms

\[
 \boxed{\quad
 \tau(t(x))=\frac{\sqrt{1+x}}{1-2x}
 {}_2F_1\!\left(\begin{matrix}1/3,2/3\\1\end{matrix};
                 \frac{27x^2}{(1-2x)^3}\right),\quad}
\]

\[
 \boxed{\quad
 \sigma(t(x))=\frac{(1+x)^{3/2}}
 {(1-16x-8x^2)(1-2x)}
 {}_2F_1\!\left(\begin{matrix}1/3,2/3\\1\end{matrix};
                 \frac{27x^2}{(1-2x)^3}\right).\quad}
\]

These identities are in characteristic zero; they are stronger and cleaner
than a fitted recurrence identity.

### 1.3 Exact coefficient formulae

Let $x=x(t)$ be the inverse series.  It satisfies

\[
 x=t\phi(x),\qquad \phi(x)=\frac{1+x}{1-8x}.
\]

Put

\[
 G_\tau(x)=\sqrt{1+x}\,h(x),\qquad
 G_\sigma(x)=\frac{(1+x)^{3/2}}{1-16x-8x^2}h(x).
\]

Lagrange--Bürmann inversion gives, for $n\geq1$,

\[
 \boxed{\quad
 \tau_n=\frac1n[x^{n-1}]G_\tau'(x)\phi(x)^n,
 \qquad
 \sigma_n=\frac1n[x^{n-1}]G_\sigma'(x)\phi(x)^n.\quad}
\]

This is already a finite hypergeometric coefficient formula.  To make it
fully explicit, write

\[
 \Phi_{n,r}=[x^r]\phi(x)^n
 =\sum_{a=0}^{\min(n,r)}\binom na
   \binom{n+r-a-1}{r-a}8^{r-a},
\]

and

\[
 g_m^{(\tau)}=[x^m]G_\tau(x)
 =\sum_{k=0}^m f_k\binom{1/2}{m-k}.
\]

Then

\[
 \tau_n=\frac1n\sum_{m=1}^n m g_m^{(\tau)}\Phi_{n,n-m}.
\]

For sigma, let $r_0=1,r_1=16$ and
\(r_j=16r_{j-1}+8r_{j-2}\).  Then

\[
 g_m^{(\sigma)}
 =\sum_{k+a+j=m} f_k\binom{3/2}{a}r_j,
 \qquad
 \sigma_n=\frac1n\sum_{m=1}^n m g_m^{(\sigma)}\Phi_{n,n-m}.
\]

The script verifies both formulae exactly through $n=39$.

### 1.4 What the half-integer recurrence does *not* imply

The spec proposed a form

\[
 c_n=\lambda^n\frac{(A)_n(B)_n}{(C)_n n!}.
\]

Its consecutive ratio must be

\[
 \frac{c_{n+1}}{c_n}
 =\lambda\frac{n^2+(A+B)n+AB}{(n+C)(n+1)}.
\]

Solving the four resulting linear equations from $n=0,1,2,3$ and testing
$n=4$ gives exact contradictions:

\[
 \begin{array}{c|c|c}
 &\text{fitted prediction at }n=4&\text{actual ratio}\\ \hline
 \tau&3812827/155270&368127/14990\\
 \sigma&14308728983/466503610&379198585/12362662.
 \end{array}
\]

Thus **[DISPROVED-HYP1]** neither sequence is one Pochhammer term.  The
displayed pullback explains why: coefficient extraction also performs
algebraic reversion.  **[GAP-HYP2]** A publication-grade exclusion of every
possible fixed two-term Pochhammer combination would require a difference
operator certificate; it is not needed for the positive closed forms above.

Most importantly, a quarter coefficient in the $t$-coordinate is not one
truncated \({}_2F_1\) value.  Applying a named Greene/Jacobsthal evaluation to
it without first handling the Lagrange reversion would be invalid.

## 2. Proof of the exact eight-class quarter law

### 2.1 Two input theorems

Caruso--Fürnsinn--Vargas-Montoya--Zudilin prove that a normalized polynomial
$B_p\in\mathbf F_p[t]$ exists with

\[
 A_p=B_p^2\quad(p\equiv1,5,7,11\bmod24),
\]

and

\[
 A_p=qB_p^2\quad(p\equiv13,17,19,23\bmod24).
\]

Since the constant coefficients are one, coefficientwise square-root
uniqueness identifies

\[
 B_p=[\tau]_{\leq(p-1)/2}\quad(\chi_p=1),\qquad
 B_p=[\sigma]_{\leq(p-3)/2}\quad(\chi_p=-1).
\]

The degrees follow because $b_{p-1}\equiv b_0=1\pmod p$.

The second input is Zhi-Wei Sun's Corollary 1.2:

\[
 A_p(1)=\sum_{n=0}^{p-1}b_n
 \equiv
 \begin{cases}
 4x^2-2p,&p=x^2+2y^2,\quad p\equiv1,3\pmod8,\\
 0,&\left(\frac{-2}{p}\right)=-1,
 \end{cases}
 \pmod p.
\]

In the first case $x\not\equiv0\pmod p$, so this says exactly

\[
 A_p(1)=0\quad\Longleftrightarrow\quad
 \left(\frac{-2}{p}\right)=-1.
\]

### 2.2 Reciprocity and the exact multiplicity at $t=1$

For $0\leq n,k<p$,

\[
 \binom{p-1-n}{k}\equiv(-1)^k\binom{n+k}{k},\qquad
 \binom{p-1-n+k}{k}\equiv(-1)^k\binom nk\pmod p.
\]

Terms with $k>n$ vanish in the second congruence.  Hence

\[
 b_{p-1-n}\equiv b_n\pmod p,
 \qquad t^{p-1}A_p(1/t)=A_p(t).
\]

Differentiating at $t=1$ shows that $A_p(1)=0$ implies $A_p'(1)=0$.
The Apéry series satisfies the third-order equation

\[
 (t^4-34t^3+t^2)y'''+(6t^3-153t^2+3t)y''
 +(7t^2-112t+1)y'+(t-5)y=0.
\]

In characteristic $p$, $F=A_pF^p$, so $A_p$ satisfies the same
differential equation.  The point $t=1$ is ordinary because the leading
coefficient is $-32\ne0\pmod p$.  If $A_p''(1)$ also vanished, uniqueness
for this third-order equation would force $A_p=0$, contradicting its
constant coefficient.  Therefore

\[
 \operatorname{ord}_{t=1}A_p=
 \begin{cases}
 2,&(-2/p)=-1,\\
 0,&(-2/p)=1.
 \end{cases}
\]

### 2.3 The branch reversal sign

Let $s_p=B_p$, $d=\deg s_p$, and
$s_p^*(t)=t^ds_p(1/t)$.  Both $1$ and $q$ are reciprocal factors, so

\[
 (s_p^*)^2=s_p^2.
\]

Since \(\mathbf F_p[t]\) is a domain and $p$ is odd,
\(s_p^*=\varepsilon_ps_p\) with \(\varepsilon_p\in\{\pm1\}\).

If $A_p(1)\ne0$, then $s_p(1)\ne0$, and evaluation at one gives
\(\varepsilon_p=1\).  If $A_p(1)=0$, the exact multiplicity result shows
that $s_p$ has a simple zero at one.  Differentiating the reversal identity
at one gives

\[
 -s_p'(1)=\varepsilon_ps_p'(1),
\]

so \(\varepsilon_p=-1\).  Thus

\[
 \boxed{\quad s_{p,d-j}=\left(\frac{-2}{p}\right)s_{p,j}
 \quad(0\leq j\leq d).\quad}
\]

This proves the previously missing leading-coefficient atom as well:
\(s_{p,d}=(-2/p)\), because $s_{p,0}=1$.

### 2.4 Forced zeros and exclusion of accidental central zeros

If $d$ is odd, reversal pairs the two middle coefficients.  If either one
were zero, both would be zero; the second-order recurrence then propagates
two consecutive zeros back to $s_{p,0}=0$, a contradiction.

Suppose the tau degree is even.  Write $p=4J+1$, so
$d=(p-1)/2=2J$.  When the reversal sign is $+1$,
\(\tau_{J+1}=\tau_{J-1}\).  If \(\tau_J=0\), the recurrence at $J$ gives

\[
 \{4(J+1)^2+(2J-1)^2\}\tau_{J-1}=0.
\]

The brace is $9/2\pmod p$, hence is nonzero for $p>3$.  This again gives
two consecutive zeros, impossible.  The sigma case has $p=4J+3$,
$d=(p-3)/2=2J$, and the corresponding brace is

\[
 4(J+1)^2+(2J+1)^2\equiv1/2\pmod p.
\]

It is also nonzero.  Therefore an even-degree middle coefficient vanishes
*if and only if* the reversal sign is $-1$.

### 2.5 Complete table

\[
\begin{array}{c|c|c|c|c}
p\bmod24&\text{selected branch}&d&(-2/p)&\text{middle result}\\ \hline
1&\tau&(p-1)/2\text{ even}&+1&\tau_{(p-1)/4}\ne0\\
5&\tau&(p-1)/2\text{ even}&-1&\tau_{(p-1)/4}=0\\
7&\tau&(p-1)/2\text{ odd}&-1&\text{middle pair nonzero}\\
11&\tau&(p-1)/2\text{ odd}&+1&\text{middle pair nonzero}\\
13&\sigma&(p-3)/2\text{ odd}&-1&\text{middle pair nonzero}\\
17&\sigma&(p-3)/2\text{ odd}&+1&\text{middle pair nonzero}\\
19&\sigma&(p-3)/2\text{ even}&+1&\sigma_{(p-3)/4}\ne0\\
23&\sigma&(p-3)/2\text{ even}&-1&\sigma_{(p-3)/4}=0.
\end{array}
\]

This is the exact relevant-branch interpretation of the eight-class law.
A stronger statement obtained by applying floor indices to *both raw series*
in every class is false: at $p=71$, both the selected
\(\sigma_{17}\) and the nonselected raw \(\tau_{17}\) vanish.

## 3. Zero sets and the failure of the proposed Jacobsthal parametrization

Define the relevant full zero set

\[
 Y_p=\{0\leq j\leq d:s_{p,j}=0\text{ in }\mathbf F_p\}.
\]

The reversal theorem proves the exact symmetry

\[
 j\in Y_p\quad\Longleftrightarrow\quad d-j\in Y_p.
\]

It also proves the central member in class $5$ and class $23$, and proves
the absence of a central member in the other six classes.  It does not
determine the remaining symmetric pairs.

### 3.1 Discriminant $-24$: the necessary correction

The reduced positive forms of discriminant $-24$ are

\[
 x^2+6y^2,\qquad 2x^2+3y^2.
\]

For primes $p>3$, genus theory gives

\[
 p=x^2+6y^2\Longleftrightarrow p\equiv1,7\pmod {24},
\]

\[
 p=2x^2+3y^2\Longleftrightarrow p\equiv5,11\pmod {24}.
\]

These four classes are exactly \((-6/p)=1\), the primes split in
\(\mathbf Q(\sqrt{-6})\).  Classes $13,17,19,23$ are inert.  In particular,
the sigma-vanishing class $23$ cannot possibly be parametrized by
$p=2x^2+3y^2$.  This is a logical obstruction, not merely a failed fit.

Only the tau class $5$ belongs to the nonprincipal-form family.  The true
uniform law uses

\[
 \text{branch }(-6/p),\quad \text{degree parity},\quad
 \text{reversal sign }(-2/p).
\]

### 3.2 Complete zero census below 3000

The full relevant-branch zero-count distributions for every prime below
$3000$ are:

\[
\begin{array}{c|c|l}
p\bmod24&\text{branch}&\#Y_p:\text{ number of primes}\\ \hline
1&\tau&0:30,\ 2:12,\ 4:4\\
5&\tau&1:47,\ 3:10\\
7&\tau&0:40,\ 2:12,\ 4:2\\
11&\tau&0:40,\ 2:12,\ 4:2\\
13&\sigma&2:42,\ 4:9,\ 6:2\\
17&\sigma&0:42,\ 2:8,\ 4:5\\
19&\sigma&0:43,\ 2:10,\ 4:1\\
23&\sigma&3:40,\ 5:12,\ 7:3.
\end{array}
\]

Thus even class $5$ is not exhausted by the forced quarter zero.  The ten
primes below $3000$ with one additional pair are

\[
101,701,1061,1181,1277,1373,1493,1949,2621,2789.
\]

For example,

\[
 p=101=2\cdot7^2+3\cdot1^2,\qquad
 Y_{101}=\{6,25,44\},
\]

where only $25=(p-1)/4$ is forced by the theorem.

There is a stronger **verified but unproved** eighth-point skeleton:

\[
 \begin{aligned}
 p\equiv13\pmod {24}:&\quad
 \{(p-5)/8,(3p-7)/8\}\subseteq Y_p,\\
 p\equiv23\pmod {24}:&\quad
 \{(p-7)/8,(p-3)/4,(3p-5)/8\}\subseteq Y_p.
 \end{aligned}
\]

It holds for every prime below $3000$.  It is marked `[GAP-EIGHTH]`; the
present reversal proof gives the symmetry and the class-$23$ middle member,
but not the outer eighth-point pair.

**[GAP-ZEROSET]** No formula for all additional pairs has been found.  The
data disprove the proposed simple linear-position/neighbor-constant laws,
and the class-$23$ inertness disproves the advertised uniform \((x,y)\)
coordinate.  What is proved is the character-forced skeleton plus reflection;
the rest is the zero set of an explicit Heun truncation.

## 4. The convolution shadow for $Z_p$

Over \(\mathbf Q\), without truncation,

\[
 \boxed{\quad b_r=\sum_{i=0}^r\tau_i\tau_{r-i}.\quad}
\]

If

\[
 c_r=\sum_{i=0}^r\sigma_i\sigma_{r-i},\qquad c_{-1}=c_{-2}=0,
\]

then $F=q\sigma^2$ gives

\[
 \boxed{\quad b_r=c_r-34c_{r-1}+c_{r-2}.\quad}
\]

Consequently, for $0\leq r<p$,

\[
 r\in Z_p
 \Longleftrightarrow
 \sum_{i=0}^r\tau_i\tau_{r-i}=0\pmod p
 \qquad(\chi_p=1),
\]

and

\[
 r\in Z_p
 \Longleftrightarrow
 c_r-34c_{r-1}+c_{r-2}=0\pmod p
 \qquad(\chi_p=-1).
\]

This is the exact convolution shadow.  Substituting the coefficient formulae
of Section 1 turns it into an explicit finite multiple hypergeometric sum.
It does **not** turn it into a condition on the support $Y_p$.

Indeed, over \(\mathbf F_5\) the sequences

\[
 a=(1,1,2),\qquad a'=(1,1,1)
\]

have the same empty zero support, but

\[
 [t^2](1+t+2t^2)^2=2\cdot2+1=0,
 \qquad
 [t^2](1+t+t^2)^2=3\ne0.
\]

Convolution vanishing depends on phases/values of all terms, not on which
terms vanish.  In particular the quarter zero deletes only one summand (two
away from the central convolution index); it does not force a zero of $b_r$.

Therefore

\[
 H(n)=\#\{p\in(n/2,n]:Q_{p,n-p}(s_p)=0\}
\]

becomes a moving quadratic-form condition in the complete universal branch
values.  **[GAP-CONV]** Without a value/phase law for those coefficients,
this is not a lattice-point count in the discriminant-$-24$ coordinates.
The negative support example shows that a parametrization of $Y_p$, even if
one existed, would still be insufficient.

## 5. The Bernoulli--Wolstenholme defect

### 5.1 A direct proof

Separate the $k=0,p$ endpoints in

\[
 b_p=\sum_{k=0}^p\binom pk^2\binom{p+k}k^2.
\]

For $1\leq k\leq p-1$, cancellation of the first-order harmonic factors
gives

\[
 \binom pk^2\binom{p+k}k^2
 \equiv\frac{p^2}{k^2}\left(1+\frac{2p}{k}\right)
 \pmod {p^4}.
\]

Hence the interior is

\[
 p^2H_{p-1}^{(2)}+2p^3H_{p-1}^{(3)}\pmod {p^4}.
\]

The classical Glaisher congruences, valid for $p\geq5$, are

\[
 H_{p-1}^{(1)}\equiv-\frac{p^2}{3}B_{p-3}\pmod {p^3},
 \qquad
 H_{p-1}^{(2)}\equiv\frac{2p}{3}B_{p-3}\pmod {p^2},
\]

and pairing $k$ with $p-k$ gives
$H_{p-1}^{(3)}\equiv0\pmod p$.  Thus the interior contributes

\[
 \frac23p^3B_{p-3}\pmod {p^4}.
\]

The $k=0$ endpoint is $1$.  Expanding

\[
 \binom{2p-1}{p-1}=\prod_{k=1}^{p-1}\left(1+\frac pk\right)
\]

and inserting the same harmonic congruences yields Glaisher's refinement

\[
 \binom{2p-1}{p-1}
 \equiv1-\frac23p^3B_{p-3}\pmod {p^4}.
\]

The $k=p$ endpoint is therefore

\[
 \binom{2p}{p}^2
 =4\binom{2p-1}{p-1}^2
 \equiv4-\frac{16}{3}p^3B_{p-3}\pmod {p^4}.
\]

Adding the three pieces proves

\[
 \boxed{\quad b_p\equiv5-\frac{14}{3}p^3B_{p-3}\pmod {p^4}.\quad}
\]

This congruence is also a specialization of Ji-Cai Liu's Theorem 1.1 in
*An extension of Gauss congruences for Apéry numbers*: in his notation take
\(A_n^{(r,s,t)}\) with \((r,s,t)=(2,2,0)\), and then \(n=m=1\).  The
resulting universal Bernoulli correction is \(-14/3\).  Thus the formula has
both the direct proof above and an independent general-theorem anchor.

The machine check covers all $93$ primes $5\leq p<500$, rather than only
the previously recorded $27/27$.

For $p\ne7$, \(\beta_p=0\) is equivalent to
$p\mid\operatorname{num}(B_{p-3})$, the Wolstenholme-prime criterion.  The
prime $p=7$ is an extra coefficient-degeneracy case because $14=0\pmod7$;
it is not a Wolstenholme prime.  This exception must be retained in any
iff statement.

### 5.2 Relation to the rank-two framework

The exact identities give the taut but useful $p$-adic rank-two expressions

\[
 \beta_p\equiv
 \frac{\sum_{i=0}^p\tau_i\tau_{p-i}-5}{p^3}\pmod p,
\]

and, with $c_r=(\sigma*\sigma)_r$,

\[
 \beta_p\equiv
 \frac{c_p-34c_{p-1}+c_{p-2}-5}{p^3}\pmod p.
\]

Thus \(\beta_p\) is in the same square-root/convolution framework, but one
level higher $p$-adically.  Its closed evaluation is controlled by harmonic
sums and $B_{p-3}$, i.e. the Wolstenholme/$p$-adic-zeta class.  The
discriminant-$-24$ coordinates only record the mod-$p$ branch; they cannot
recover this order-three lift.

**[GAP-GK]** A Gross--Koblitz expression would have to retain $p$-adic
derivatives (the harmonic terms above), not just a mod-$p$ Jacobi sum.  No
verified identity reducing those derivatives to the proposed \((x,y)\) data
was found.  Calling the displayed algebraic pullback a direct rank-two
Gross--Koblitz formula would overstate the result.

## 6. Reproducibility

Run from `problems/3.2`:

```text
python3 research/scripts/codex_max_hypergeom_quarter.py
python3 research/scripts/codex_max_zero_sets.py
python3 research/scripts/codex_max_convolution_beta.py
```

The checks are:

- `codex_max_hypergeom_quarter.py`: exact rational recurrence, square-root,
  Franel \({}_2F_1\) pullback, $q$-pullback, and Lagrange formula through
  $n=39$; exact one-Pochhammer contradictions.
- `codex_max_zero_sets.py`: every relevant prime $p<3000$; factor-branch
  recurrence, reversal sign, exact eight-class quarter table, all $219$
  ordinary double-root cases at $t=1$, full zero-set census,
  discriminant-$-24$ representations, the verified eighth skeleton, and the
  $p=71$ counterexample to the stronger raw-both-branches reading.
- `codex_max_convolution_beta.py`: both convolution identities through
  $n=59$; the beta formula for $93$ primes below $500$; $279$ harmonic
  congruences; $1032$ individual interior summand expansions; explicit
  support-only counterexample.

## 7. Precise references

1. F. Beukers, [*Consequences of Apéry's work on \(\zeta(3)\)*
   (2003)](https://webspace.science.uu.nl/~beuke106/caen.pdf), especially
   pp. 1--2 for the third-order equation, symmetric-square equation, and the
   half-integer recurrence.
2. X. Caruso, F. Fürnsinn, D. Vargas-Montoya, W. Zudilin,
   [*Galois Groups of Apéry-like Series Modulo Primes*](https://arxiv.org/abs/2510.23298),
   Bull. Aust. Math. Soc. 114 (2026), 65--78,
   doi:10.1017/S0004972725100932.  Theorem 1.2 (Theorem 2 in the preprint)
   is the square/\(q\)-square factorization; pp. 1, 3--4 give the Franel and
   \({}_2F_1\) pullbacks.
3. Z.-W. Sun, [*On sums of Apéry polynomials and related
   congruences*](https://arxiv.org/abs/1101.1946), J. Number Theory 132
   (2012), 2673--2699.  Corollary 1.2, equation (1.9), is the exact
   $A_p(1)$ evaluation used here.
4. B. K. Spearman and K. S. Williams,
   [*Representing primes by binary quadratic forms*](https://people.math.carleton.ca/~williams/papers/pdf/172.pdf),
   Amer. Math. Monthly 99 (1992), 423--426, for the discriminant-$-24$
   representation criteria.
5. S. Hong, [*Notes on Glaisher's
   Congruences*](https://doi.org/10.1007/BF02731955), Chinese Ann. Math. Ser. B
   21 (2000), 33--38, for the general inverse-power/Bernoulli congruences;
   the special cases needed for beta are also derived in Section 5 above.
6. J.-C. Liu, [*An extension of Gauss congruences for Apéry
   numbers*](https://arxiv.org/abs/2404.16636), arXiv:2404.16636 (2024),
   Theorem 1.1.  Its \((r,s,t)=(2,2,0)\), \(n=m=1\) specialization gives
   the same (b_p\bmod p^4) Bernoulli correction.

## 8. Final gap ledger

- `[GAP-HYP2]`: exclude every possible fixed two-term Pochhammer
  representation by a formal difference-operator certificate.
- `[GAP-EIGHTH]`: prove the class $13/23$ outer eighth-point zero pairs.
- `[GAP-ZEROSET]`: determine the remaining symmetric zero pairs of the Heun
  truncations.  A uniform $p=2x^2+3y^2$ parametrization is impossible as
  stated.
- `[GAP-CONV]`: find a complete value/phase law strong enough to simplify the
  quadratic convolution defining $Z_p$.  Zero support alone cannot do so.
- `[GAP-GK]`: express the order-three Bernoulli defect via verified
  $p$-adic gamma derivatives, if a useful Gross--Koblitz form exists.

The requested quarter theorem, hypergeometric pullbacks, exact convolution
shadow, and beta formula are closed.  The remaining items are genuine gaps,
not missing numerical checks.
