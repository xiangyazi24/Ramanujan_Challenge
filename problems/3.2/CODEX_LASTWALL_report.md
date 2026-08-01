NO UNBOUNDED \(L\) IS PROVED; THE \(3/2\) RECORD REMAINS UNBROKEN.

# Last-wall report

## 1. Terminal verdict

Let
\[
N=p-2,\qquad I_p=\{1,\ldots,N\},
\]
and retain the nonwrapping definitions in `CODEX_SPEC_lastwall.md`.  I did
not prove either
\[
S_D\ll N\qquad\hbox{or}\qquad Q_D\ll N
\]
at \(D=\sqrt N\,L(p)\) for any unbounded \(L(p)\).  Thus this strike does
not improve the \(3/2\) exponent.

The strongest unconditional conclusions obtained or fully re-audited here
are the following.

1. **Projective-fibre restart theorem.**  For every prime \(p\ge7\), if
   \(J\subset I_p\) has span
   \(H\), \(v\in\mathbf P^1(\mathbf F_p)\), and
   \(m_J(v)=\#\{n\in J:\pi(n)=v\}\), then for every \(1\le K\le H\),
   \[
   m_J(v)-1
   \le \left\lfloor\frac{H}{K+1}\right\rfloor
      +\frac{3K(K-1)}2.
   \]
   Consequently,
   \[
   m_J(v)-1
   \le\min_{1\le K\le H}
   \left\{
      \left\lfloor\frac{H}{K+1}\right\rfloor
      +\frac{3K(K-1)}2
   \right\}
   \le4H^{2/3},
   \tag{1.1}
   \]
   where the last inequality follows by taking
   \(K=\lceil H^{1/3}\rceil\).  Consequently, when
   \(H_r=\min(D,N-r)\ge1\),
   \[
   d_D(r)\le
   \min_{1\le K\le H_r}
   \left\{
      \left\lfloor\frac{H_r}{K+1}\right\rfloor
      +\frac{3K(K-1)}2
   \right\}
   \le4H_r^{2/3}.
   \tag{1.2}
   \]
   If \(H_r=0\), then \(d_D(r)=0\) and the assertion is trivial.

2. **Primitive-support capacity theorem.**  Put
   \(P_D=\#\{r:d_D(r)>0\}\), \(A_D=S_D-P_D\).  Then
   \[
   P_D\le N,\qquad A_D\le Q_D,
   \tag{1.3}
   \]
   and, more sharply,
   \[
   2Q_D\ge A_D+\frac{A_D^2}{P_D},
   \qquad
   S_D\le
   \frac{P_D+\sqrt{P_D^2+8P_DQ_D}}2.
   \tag{1.4}
   \]
   The case \(P_D=0\) is interpreted separately and is trivial.  Replacing
   \(P_D\) by \(N\) recovers the banked capacity inequality; primitive
   decomposition does not create a stronger unconditional contraction.

3. Combining (1.4) with the banked inverse-square bound gives
   \[
   Q_D\le66D^2(1+\log D),
   \tag{1.5}
   \]
   and hence
   \[
   S_D\le
   \min\left\{
      \frac32D(D-1),
      \frac{N+\sqrt{N^2+528ND^2(1+\log D)}}2
   \right\}.
   \tag{1.6}
   \]
   At \(D=\sqrt N\,L\), this is only
   \(S_D\ll N(1+L\sqrt{\log D})\), not \(O(N)\).

The closing obstruction is now sharply localized.  Restart controls how
many returns can follow one base, while primitive decomposition caps the
number of first-return arrows by \(N\).  Neither controls the distribution
of the nonprimitive tail strongly enough.  That missing information is a
quenched, growing-lag two-return estimate, not another scalar moment
iteration.

## 2. Exact algebra used in all three vectors

Write
\[
P(n)=34n^3+51n^2+27n+5
\]
and normalize the gap continuants by
\[
N_0(X)=0,\quad N_1(X)=1,\quad N_2(X)=P(X+1),
\]
\[
N_{h+1}(X)=P(X+h)N_h(X)-(X+h)^6N_{h-1}(X).
\tag{2.1}
\]
For \(1\le r<r+h\le N\), define the companion solution
\[
y_n^{(r)}=b_rc_n-c_rb_n.
\]
Direct recurrence propagation gives the exact normalization
\[
y_{r+h}^{(r)}
=\frac{N_h(r)}{\prod_{j=1}^h(r+j)^3}.
\tag{2.2}
\]
All denominator factors are nonzero in the physical triangle.  Thus
\[
N_h(r)=0
\quad\Longleftrightarrow\quad
\det\begin{pmatrix}b_r&c_r\\b_{r+h}&c_{r+h}\end{pmatrix}=0
\quad\Longleftrightarrow\quad
\pi(r)=\pi(r+h).
\tag{2.3}
\]

The exact addition law is
\[
N_{a+g}(r)
=N_g(r+a)N_{a+1}(r)
 -(r+a+1)^6N_{g-1}(r+a+1)N_a(r).
\tag{2.4}
\]
If \(N_a(r)=0\) in the physical triangle, then
\(N_{a+1}(r)\ne0\): otherwise \(\pi(r+a)=\pi(r+a+1)\), contradicting
\(N_1=1\).  Therefore
\[
N_{a+g}(r)=0
\quad\Longleftrightarrow\quad
N_g(r+a)=0
\qquad\text{when }N_a(r)=0.
\tag{2.5}
\]
This is the exact root-level renewal/restart equivalence.

Equations (2.1)--(2.5), including every admissible instance of (2.4) in
the three requested finite ranges, were checked independently by
`CODEX_LASTWALL_verify.py`.

## 3. Vector (7): restart for every projective fibre

### 3.1 The moving-target issue disappears after scalarization — PROVED

For a target \(v=[\alpha:\beta]\), set
\[
y_n^{(v)}=\beta b_n-\alpha c_n.
\tag{3.1}
\]
This is a nonzero homogeneous solution of the same Apéry recurrence, and
\[
y_n^{(v)}=0\quad\Longleftrightarrow\quad\pi(n)=v.
\tag{3.2}
\]
Thus an arbitrary projective fibre is an ordinary zero set of one scalar
solution.  At a regular zero \(z\), consecutive zeros are impossible, and
the Riccati state
\[
x_z=[z^3y_{z-1}:y_z]
\]
is \(\infty\).  One allowed transfer sends it to the fixed state \(0\).
The target moves only if one insists on following the original projective
orbit; it does not move in the scalar zero problem (3.1).

This statement is specifically about projective fibres.  It does not apply
to an affine level set \(b_n=A\ne0\), because \(z_n=b_n-A\) satisfies the
forced recurrence
\[
(n+1)^3z_{n+1}
=P(n)z_n-n^3z_{n-1}+4A(2n+1)^3.
\tag{3.3}
\]

### 3.2 Uniform window theorem — PROVED

First, \(N_h\bmod p\) is not the zero polynomial for
\(1\le h\le p-3\).  For \(h\ge2\), the two exact evaluations are
\[
N_h(-1)=b_{h-1}((h-1)!)^3,
\qquad
N_h(-2)=-5b_{h-2}((h-2)!)^3.
\tag{3.4}
\]
If both vanished, \(b_{h-1}\) and \(b_{h-2}\) would be consecutive
Apéry zeros, which the Casoratian/recurrence forbids.  Hence
\[
\#\{x\in\mathbf F_p:N_h(x)=0\}\le\deg N_h\le3(h-1).
\tag{3.5}
\]

Now list the occurrences of one fibre in \(J\) as
\(x_1<\cdots<x_m\).  An adjacent occurrence gap \(h=x_{i+1}-x_i\)
makes \(x_i\) a root of \(N_h\).  The number of adjacent gaps exceeding
\(K\) is at most \(\lfloor H/(K+1)\rfloor\); for each
\(1\le h\le K\), (3.5) bounds the number of gaps equal to \(h\) by
\(3(h-1)\).  Therefore
\[
m-1\le
\left\lfloor\frac H{K+1}\right\rfloor
+\sum_{h=1}^K3(h-1)
=\left\lfloor\frac H{K+1}\right\rfloor
+\frac{3K(K-1)}2.
\tag{3.6}
\]
For \(t=H^{1/3}\), choosing \(K=\lceil t\rceil\) bounds the right side
by
\[
t^2+\frac32(t+1)t\le4t^2,
\]
which proves (1.1).  Applying it to
\(J=[r,r+H_r]\) and \(v=\pi(r)\) proves (1.2).

### 3.3 Cut-edge pollution — completely classified and harmless here

A raw column can violate every \(H^{2/3}\) zero bound.  If
\(x=-m\pmod p\) and \(p\mid b_{m-1}\), endpoint factorization gives
\[
N_h(x)=0\qquad(h\ge m).
\tag{3.7}
\]
For the corresponding physical base \(r=p-m\), however,
\[
h\le N-r=(p-2)-(p-m)=m-2.
\tag{3.8}
\]
Thus the polluted ray starts two levels beyond the last admissible physical
lag.  It cannot contaminate (1.2).  The verifier found, for example, the
actual \(p=499\) ray \((x,m)=(431,68)\): it vanishes from level 68 on,
whereas the physical window ends at level 66.

### 3.4 Wall verdict — DEAD as a closing mechanism

Let \(M_D=\max_r d_D(r)\).  Equations (1.2) and the row-degree bound give
\[
M_D\le4D^{2/3},
\qquad
S_D\le\sum_{d\le D}3(d-1)=\frac32D(D-1).
\]
Consequently
\[
Q_D
\le\frac{M_D-1}{2}S_D
\le3D^{8/3}.
\tag{3.9}
\]
At \(D=\sqrt N\,L\), (3.9) is
\[
Q_D\le3N^{4/3}L^{8/3},
\]
far larger than \(N\).  Feeding \(Q_D\le(M_D-1)S_D/2\) into capacity
only gives \(S_D\le NM_D\), also worse.  Restart controls tail length per
base but supplies no bound on the intensity of first returns across bases.

**Terminal label for vector (7):** the general-fibre theorem is `PROVED`;
its use as a proof of `MESO-TOTAL` or `MESO-PAIR` is `DEAD` at the displayed
exponent.

## 4. Vector (8): primitive decomposition and bootstrap

For each base put
\[
\mathcal R_D(r)=
\{d\le D:r+d\le N,\ \pi(r+d)=\pi(r)\},
\qquad k_r=|\mathcal R_D(r)|.
\]
When nonempty, write
\(d_1(r)<\cdots<d_{k_r}(r)\).

### 4.1 Primitive collisions are exactly first returns — PROVED

A collision \((r,d)\) is primitive precisely when there is no
\(d'<d\) in \(\mathcal R_D(r)\).  Indeed, if \(d'<d\), then
\[
\pi(r+d')=\pi(r)=\pi(r+d),
\]
so \((r+d',d-d')\) is automatically the second collision in a split.
Thus each active base contributes exactly one primitive arrow and
\[
P_D=\#\{r:k_r>0\}.
\tag{4.1}
\]

The endpoint map for primitive arrows is injective.  If
\(r_1<r_2\) had the same endpoint \(s\), then
\(\pi(r_1)=\pi(r_2)=\pi(s)\) and
\(r_2-r_1<s-r_1\), contradicting that \(s\) is the first return from
\(r_1\).  Hence
\[
P_D\le N.
\tag{4.2}
\]
This is sharp in order for general reflection-symmetric words: disjoint
reflection pairs can give linearly many primitive arrows.

### 4.2 Exact split and renewal identities — PROVED

The \(j\)-th return from a base has exactly \(j-1\) earlier split points.
Therefore
\[
A_D:=\#\{\text{nonprimitive collisions}\}
=\sum_r(k_r-1)_+=S_D-P_D,
\tag{4.3}
\]
while the total split multiplicity is
\[
\sum_r\sum_{j=1}^{k_r}(j-1)
=\sum_r\binom{k_r}{2}=Q_D.
\tag{4.4}
\]
In particular, \(A_D\le Q_D\) and \(S_D\le N+Q_D\).

Equivalently, if
\[
C_p(a,g)=
\#\{r:\pi(r)=\pi(r+a)=\pi(r+a+g)\},
\]
then the bijection
\[
(d_i,d_j)\longleftrightarrow(a=d_i,g=d_j-d_i)
\]
gives
\[
Q_D=
\sum_{\substack{a,g\ge1\\a+g\le D}}C_p(a,g).
\tag{4.5}
\]
The polynomial realization of (4.5) is exactly (2.4)--(2.5).

There is also an exact first-return recursion.  If \(\tau(r)=d_1(r)\),
then
\[
k_D(r)=
\begin{cases}
0,&\tau(r)>D,\\
1+k_{D-\tau(r)}(r+\tau(r)),&\tau(r)\le D.
\end{cases}
\tag{4.6}
\]
Summing its tail merely recovers (4.3)--(4.4); it has coefficient one and
does not contract.

### 4.3 The strongest scalar bootstrap is capacity again — PROVED

On the \(P_D\) active bases put \(x_r=k_r-1\).  Then
\[
A_D=\sum x_r,
\qquad
2Q_D=\sum(x_r^2+x_r)
=A_D+\sum x_r^2.
\]
Cauchy gives
\[
2Q_D\ge A_D+\frac{A_D^2}{P_D}.
\tag{4.7}
\]
Solving this quadratic for \(A_D\), then adding \(P_D\), proves (1.4).
Using only \(P_D\le N\), (1.4) is exactly the positive-root form of
\(S_D^2\le N(S_D+2Q_D)\).  Thus primitive decomposition explains the
capacity inequality but does not improve it without a new estimate on
\(P_D\) or \(Q_D\).

Substitution of (1.5) yields (1.6), which misses the target by
\(L\sqrt{\log D}\).  The other proposed feedbacks also fail:

- \(Q_D\le(M_D-1)S_D/2\) has a coefficient growing like
  \(D^{2/3}\), so it cannot be absorbed into \(S_D\le N+Q_D\).
- Adjacent coprimality proves the nonconsecutive-zero and restart facts
  already used; it gives no saving on the number of first-return endpoints.
- The renewal triangle (4.5) is exactly \(Q_D\), not an upper bound for it.
- Passing to the third moment only rewrites the same tail: one has
  \(Q_D\le N+3\sum_r\binom{k_r}{3}\), but there is no bound for the new
  third-moment term.  The \(k_r=2\) part merely consumes the environmental
  \(O(N)\) allowance.

The exact conditional improvement exposed by (1.4) is
\[
P_DQ_D\ll N^2.
\tag{4.8}
\]
With the current \(Q_D\ll NL^2\log D\), it would suffice to prove the new
primitive-support estimate
\[
P_D\ll\frac{N}{L^2\log D}.
\tag{4.9}
\]
No banked restart, renewal, resultant, or apparition statement implies
(4.9).

### 4.4 A reflection-symmetric no-go word — PROVED abstractly, VERIFIED-N

The following construction shows that the generic inputs used above cannot
by themselves prove `MESO-PAIR`.  It is not an Apéry counterexample.

Let \(q\) be an odd prime and
\[
a_n=2qn+\bigl((n\bmod q)^2\bmod q\bigr),
\qquad0\le n<q^2.
\tag{4.10}
\]
These \(q^2\) positions are strictly increasing and have span
\[
H=2q^3-2q+1.
\tag{4.11}
\]
Let \(R_d\) be their difference multiplicity.  A difference arising from
an index gap \(k\) lies in
\([2qk-(q-1),2qk+(q-1)]\); these intervals are disjoint, so \(d\)
determines \(k\).  If \(q\nmid k\), reduction modulo \(q\) gives the
linear congruence
\[
2ki+k^2\equiv d-2qk\pmod q,
\]
so there are at most \(q\) possible indices and \(R_d\le q<d\).  If
\(q\mid k\), the residue correction is zero, \(d=2qk\ge2q^2\), and
\(R_d\le q^2\le d/2\).  Hence
\[
R_d\le d.
\tag{4.12}
\]

To impose the exact orbit-reflection pattern, let \(p_q\) be the first
prime at least \(q^5\), put \(N_q=p_q-2\), place one repeated colour at
\(1+a_n\) and all reflected positions \(N_q+1-(1+a_n)\), and give every
remaining reflection pair its own private colour.  Bertrand's bound gives
\(p_q<2q^5\).  The two special clusters are more than \(H\) apart.  For
every \(d\le H\), the resulting word therefore has
\[
C_d=2R_d+\mathbf1_{2\mid d}\le3(d-1),
\tag{4.13}
\]
including exactly one forced mirror collision at every even gap.  It has
the exact reflection \(w_n=w_{N_q+1-n}\), no adjacent equality, the usual
primitive and renewal identities of a genuine word, and
\[
\max_r d_H(r)=q^2-1\le4H^{2/3}.
\tag{4.14}
\]
Nevertheless,
\[
Q_H=2\binom{q^2}{3},
\qquad
\frac{H^2}{N_q}\ge q,
\qquad
\frac{Q_H}{N_q}\ge\frac q7\quad(q\ge5).
\tag{4.15}
\]
Thus \(H/\sqrt{N_q}\to\infty\) while \(Q_H/N_q\to\infty\), despite
(4.13), reflection, the restart-scale cap, and all word-level
primitive/renewal identities.

This construction does **not** satisfy the Apéry polynomial recurrence, so
it only proves that the listed abstract consequences cannot close the
\(Q_D\) bootstrap.  It also has \(S_H\ll N_q\), so it does not refute a
direct arithmetic proof of `MESO-TOTAL`.  Exact full-difference gates for
\(q=3,5,7,11,13,17,19\) passed.

**Terminal label for vector (8):** primitive decomposition, injection, and
renewal are `PROVED`; the scalar/bootstrap route with current inputs is
`DEAD`; (4.8) and (4.9) are the precise `CONDITIONAL` escape hatches.

## 5. Inverse-square plus small-gap arithmetic

Put
\[
\mathcal A_p(D)=\sum_{d\le D}\frac{C_d}{d^2}.
\tag{5.1}
\]
The banked localized clique inequality is
\[
Q_D\le22D^2\mathcal A_p(D).
\tag{5.2}
\]
Therefore, for (5.2) itself to certify \(Q_D\le K N\), the exact required
condition is
\[
\mathcal A_p(D)\le\frac{KN}{22D^2}.
\tag{5.3}
\]
At \(D=\sqrt N\,L\), this is
\[
\mathcal A_p(D)=O(L^{-2}).
\tag{5.4}
\]

### 5.1 The raw certificate is impossible — DEAD

Orbit reflection gives
\[
\pi(n)=\pi(p-1-n).
\]
For every even \(2\le d\le D\), the physical base
\[
r=\frac{p-1-d}{2}
\]
is consequently a gap-\(d\) collision.  Hence
\[
C_d\ge1\quad(2\mid d),
\qquad
\mathcal A_p(D)
\ge\sum_{\substack{d\le D\\2\mid d}}\frac1{d^2}
\ge\frac14.
\tag{5.5}
\]
The right side of (5.4) tends to zero for every unbounded \(L\).  Thus no
amount of upper-bound information about the remaining roots can make the
raw nonnegative certificate (5.2) prove `MESO-PAIR`.  This is a failure of
the certificate, not a lower bound on the true \(Q_D\): at \(p=997\) the
verifier found \(Q_D=0\) while the right side of (5.2) was more than
\(45N\).

This also answers the proposed polylogarithmic small-gap split.  Every
prefix ending at \(Y\ge2\) already contributes at least \(1/4\), even if
the entire tail were zero.  The exact \(h=2\) law, for \(p>13\), is
\[
C_2=1+2\,\mathbf1_{\left(\frac{-51}{p}\right)=1}.
\tag{5.6}
\]
After formally subtracting the forced mirror root, (5.6) still leaves
weight \(1/2\) on the positive-density Chebotarev class with
\((-51/p)=1\).  More generally, the fixed-\(d\) Chebotarev mean
\(C_d\sim1+\mathbf1_{2\mid d}\) is a mean over primes with \(d\) fixed;
it neither has the fixed-prime/growing-\(d\) quantifiers required here nor
makes (5.1) tend to zero.  Even a hypothetical uniform \(C_d=O(1)\) would
only give \(\mathcal A_p(D)=O(1)\), still not (5.4).

### 5.2 Leading-coefficient apparition does not supply a root bound — DEAD

The apparition law identifies rows where the nominal leading coefficient
of \(N_d\) vanishes.  It says where a degree drop can occur; it does not
bound the number of roots of the remaining polynomial.  The exact finite
audit illustrates the distinction:

- for \(p=199\), \(\rho_p=9\); row \(d=9\) has degree 22 instead of 24
  and \(C_9=0\), while row \(d=18\) has degree 49 instead of 51 and
  \(C_{18}=7\);
- for \(p=499,997\), the ranks 125 and 499 lie beyond the tested windows,
  so there is no degree-drop row there at all.

There is no implication from these facts to the moving-range weighted
estimate (5.4).

The banked axis-strip lemma can pay
\(Q_D^{\mathrm{axis}}(G)\le12D^{2/3}G^2=O(N)\) at
\(G\asymp N^{1/3}L^{-1/3}\).  This correctly removes all fixed small-gap
certificates, but it leaves the balanced two-return sum
\[
\sum_{\substack{a,g>G\\a+g\le D}}C_p(a,g),
\tag{5.7}
\]
which no one-gap estimate on \(C_d\) controls.  A tail-only inverse-square
charge cannot simply discard the short edges used by the clique-energy
proof.  Controlling (5.7) is precisely a new balanced pair-correlation
input, not a consequence of Chebotarev or apparition currently in hand.

**Terminal label for this vector:** the raw inverse-square/small-\(d\)
route is `DEAD` by (5.5); a centered or codegree-sensitive replacement
which ignores isolated mirror edges remains `CONDITIONAL` and would require
new two-gap arithmetic.

## 6. Remaining banked routes

The near-wall estimate \(K_p(H,D)\ll HD^{11/7}\) and the
\(W\)-injection are valuable anisotropic estimates, but their banked
bootstrap reaches the same global fixed point as (1.5), namely a
\(D^2\) bound up to logarithms.  At \(D=\sqrt N\,L\), this is still
\(NL^2\) times logarithms.  Neither reverses (4.5) nor supplies the
balanced quenched input (5.7).  Its verdict for the last-wall target is
therefore `DEAD with the present conversion`, not refuted as a theorem.

## 7. Exact machine audit

The persistent verifier is `CODEX_LASTWALL_verify.py`.  It computes
\(D=\lfloor p^{3/5}\rfloor\) by the integer condition
\[
D^5\le p^3<(D+1)^5,
\]
not by floating-point rounding.  All finite-field identities and all
reported rational masses use exact integer or `Fraction` arithmetic.
The audited verifier SHA-256 is
`aa2346dc9ec260ee6b3cdf5929053e61b8c04616edc24731401687cba7375630`.

Run:

```sh
PYTHONPYCACHEPREFIX=/tmp/codex-lastwall-pycache \
  python3 -m py_compile CODEX_LASTWALL_verify.py
python3 CODEX_LASTWALL_verify.py
```

The requested endpoint data are:

| \(p\) | \(N\) | \(D\) | \(S_D\) | \(Q_D\) | \(\max d_D\) | \(P_D\) | \(A_D\) | \(\mathcal A_p(D)\) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 199 | 197 | 23 | 37 | 2 | 2 | 35 | 2 | 0.622029230311 |
| 499 | 497 | 41 | 70 | 14 | 3 | 58 | 12 | 0.858421582855 |
| 997 | 995 | 62 | 85 | 0 | 1 | 85 | 0 | 0.539260467051 |

The split histograms are `{1: 2}`, `{1: 10, 2: 2}`, and `{}` in the same
order.  The exact raw-certificate failure factors
\(22D^2\mathcal A_p(D)/N\) are 36.747087, 63.875547, and 45.833346.

For every integer cutoff \(1\le H\le D\), not only at the endpoint, the
script verifies:

- dense \(\mathbf F_p[X]\) construction, degree, nonzero status, and root
  counts for every \(N_h\);
- \(N_h(r)\), companion determinant, and projective-collision equivalence;
- reflection, the forced even-gap root, and the exact \(h=2\) law;
- the full addition law and every root-level restart instance;
- every represented projective fibre as the zero set of a homogeneous
  recurrence solution;
- the exact short/long-gap window bound for every fibre;
- \(S=P+A\), primitive endpoint injection, every split witness,
  split mass \(=Q\), and first-return renewal;
- the first/second moment identities, capacity, restart cap, localized
  clique energy, and inverse-square bounds;
- every raw polluted ray arising from an Apéry zero and its exclusion from
  the physical triangle.

The aggregate audit counts are:

| gate | exact checks |
|---|---:|
| physical gap/determinant/orbit identities | 83,508 |
| addition-law identities | 2,232,496 |
| root-level restart equivalences | 4,068 |
| homogeneous fibre recurrence equations | 499,368 |
| fibre/cutoff window maxima | 33,951 |
| complete cutoffs | 126 |

The script also exhausts every difference in the abstract models for
\(q=3,5,7,11,13,17,19\), checks the exact primitive/renewal counts, and
checks the reflection-symmetric prime-length construction (4.10)--(4.15).
All gates print `PASS`.

## 8. Status ledger

| Statement or route | Status | Terminal reason |
|---|---|---|
| Companion determinant and addition/restart identities | `PROVED`, `VERIFIED-N` | Exact recurrence algebra |
| Restart for an arbitrary projective fibre | `PROVED`, `VERIFIED-N` | Scalar solution (3.1) has a fixed zero restart |
| Uniform fibre-window bound (1.1)--(1.2) | `PROVED`, `VERIFIED-N` | Short/long adjacent-gap count plus nonzero \(N_h\) |
| Cut-edge pollution exclusion in the physical triangle | `PROVED`, `VERIFIED-N` | Pollution begins at \(m\), physical lags end at \(m-2\) |
| Restart alone \(\Rightarrow\) `MESO-PAIR`/`TOTAL` | `DEAD` | Gives only \(Q_D\ll D^{8/3}\) |
| Primitive = first return; endpoint injection | `PROVED`, `VERIFIED-N` | Gives the sharp environmental cap \(P_D\le N\) |
| Split mass and renewal triangle equal \(Q_D\) | `PROVED`, `VERIFIED-N` | Exact bijection, but no upper saving |
| Primitive scalar bootstrap | `DEAD` | Its strongest inequality is capacity again |
| \(P_DQ_D\ll N^2\) or (4.9) | `CONDITIONAL` | Would close, but is a new arithmetic input |
| Raw inverse-square small-gap certificate | `DEAD` | Forced \(d=2\) mirror root gives \(\mathcal A_p(D)\ge1/4\) |
| Fixed-\(d\) Chebotarev for the needed range | `DEAD` | Wrong quantifier and nonvanishing baseline |
| Leading-coefficient apparition as a root-count bound | `DEAD` | Detects degree drops only |
| Axis deletion plus balanced two-gap estimate | `CONDITIONAL` | Requires new control of (5.7) |
| `MESO-TOTAL` or `MESO-PAIR` at unbounded \(L\) | **NOT PROVED** | The last wall remains |
