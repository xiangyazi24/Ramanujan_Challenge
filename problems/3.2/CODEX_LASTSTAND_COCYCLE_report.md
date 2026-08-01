# NOT PROVED: `[FR_eta]` remains open

## Terminal verdict

The requested fixed power saving was not proved.  All three prescribed
routes were carried to an exact stopping point, but the specification has two
foundational formulation errors which must be separated from the still-open
arithmetic problem.

1. The asserted first-order action
   \(\xi_{r+1}=M_r\xi_r\), for \(\xi_r=(b_r:c_r)\), is false.  The companion
   matrix advances the two-time state of one scalar solution, not the row
   formed by two independent solutions.  The correct collision condition is
   the fixed-Borel condition
   \[
   (T_{r,d})_{21}=0,
   \]
   not the condition that \(T_{r,d}\) fix \(\xi_r\).  Consequently the
   advertised two-run eigenvector, three-run curve, stepped-pair eigenline,
   and slowly moving eigenline all address the wrong dynamical object.
2. The displayed target in this specification counts shell incidences
   \(\sum_{d\sim D}|Z_d|\).  The target called `[FR_eta]` in
   `campaign3_questions/CTX_LASTSTAND.txt` counts the union/primitive support
   \(|\bigcup_{d\sim D}Z_d|\).  These are not equivalent; only the latter is
   bounded by the former.

After correcting those points, the exact route verdicts are:

- **R1:** horizontal dispersion gives an exact sliding-window estimate
  \(S_J\ll N+DN^{2/3}\), but it is no better than the degree budget at the
  intended scale \(D=\sqrt N L\), \(L=N^{o(1)}\).  Difference closure at a
  fixed base is false.  A private-pair construction shows that every
  pair-correlation input used by R1 can vanish while
  \(S_J=|\bigcup Z_d|=D^{2-o(1)}\).
- **R2:** the requested positive-density family of shifted correlations below
  the Bezout scale already follows from the banked triangle estimate: at
  least half of the relevant gap pairs have correlation \(O(\log D)\).
  Such a correlation can meet a shell-singleton, but an upper bound on these
  correlations gives no upper bound on the singleton first moment.  The
  private-pair obstruction makes every correlation zero while retaining
  \(D^{2-o(1)}\) singleton mass.  Characteristic-zero coprimality and
  resultants also have the wrong prime/gap quantifier order.
- **R3:** the raw companion matrices are nowhere slowly varying in the
  relevant \(p\)-adic sense.  An exact resultant shows that every adjacent
  pair has a unit-sized coefficient change and is noncommuting.  Except at
  at most twenty residues, adjacent matrices are separable and share no
  eigenline even over the algebraic closure.  Thus the proposed adiabatic
  diagonalization cannot start.  This local obstruction is not itself a
  long-product anti-concentration theorem.

The strongest clean surviving condition is a quenched **low-multiplicity
first-return bound**, stated in Section 7.  It is genuinely new arithmetic
input, not a consequence of the present degree, triangle, resultant, or local
matrix identities.  Therefore the \(3/2\) energy record remains unbroken.

All numerical values in this report are reproduced by
`CODEX_LASTSTAND_COCYCLE_verify.py`.

## 1. Foundational audit: the correct cocycle

### 1.1 Conventions

Work on the regular, nonwrapping interval
\[
I_p=\{1,\ldots,N\},\qquad N=p-2,
\]
and require \(r+d\le N\).  Put
\[
P(s)=34s^3+51s^2+27s+5,
\]
so that each scalar solution satisfies
\[
(s+1)^3y_{s+1}=P(s)y_s-s^3y_{s-1}.
\]
The verifier uses the two solutions
\[
(b_0,b_1)=(1,5),\qquad(c_0,c_1)=(0,6).
\]
Rescaling the second solution changes projective coordinates but not any
collision set.

Define
\[
\alpha_s=\frac{P(s)}{(s+1)^3},\qquad
\beta_s=\frac{s^3}{(s+1)^3},\qquad
M_s=\begin{pmatrix}0&1\\-\beta_s&\alpha_s\end{pmatrix}.
\]
For one scalar solution,
\[
\binom{y_s}{y_{s+1}}=M_s\binom{y_{s-1}}{y_s}.
\]
This is the valid companion action.

### 1.2 The asserted action on \(\xi_s\) is false

The point \(\xi_s=(b_s:c_s)\) is made from two different solutions at one
time.  It is not the two-time state of either solution.  The failure already
occurs over \(\mathbf Q\) at \(s=1\):
\[
(b_1,c_1)=(5,6),\qquad (b_2,c_2)=(73,351/4),
\]
whereas
\[
M_1(5,6)^T=(6,697/8)^T.
\]
The projective determinant is
\[
\det\begin{pmatrix}6&73\\697/8&351/4\end{pmatrix}
=-46669/8\ne0.
\]
Thus even projective equality fails.

### 1.3 `[FRAME-BOREL]`: exact corrected formulation

Set
\[
G_r=
\begin{pmatrix}
b_{r-1}&c_{r-1}\\
b_r&c_r
\end{pmatrix},
\qquad
T_{r,d}=M_{r+d-1}\cdots M_r.
\]
Then
\[
G_{r+1}=M_rG_r,qquad G_{r+d}=T_{r,d}G_r.       \tag{1.1}
\]
The Casoratian gives
\[
\det G_r=\frac6{r^3}\ne0\pmod p                \tag{1.2}
\]
throughout the physical interval.  Since \(\xi_r\) is the projective class
of the second row of \(G_r\), (1.1)--(1.2) imply
\[
\begin{aligned}
\xi_{r+d}=\xi_r
&\iff e_2^TT_{r,d}G_r\ \parallel\ e_2^TG_r\\
&\iff e_2^TT_{r,d}\ \parallel\ e_2^T\\
&\iff (T_{r,d})_{21}=0.                         \tag{1.3}
\end{aligned}
\]
This is a hit on one fixed Borel subgroup.  It is not a fixed-point condition
depending on \(\xi_r\).

For the continuants
\[
N_0(X)=0,\quad N_1(X)=1,
\]
\[
N_{d+1}(X)=P(X+d)N_d(X)-(X+d)^6N_{d-1}(X),      \tag{1.4}
\]
induction on \(d\) gives the more explicit identity
\[
(T_{r,d})_{21}
=-\frac{r^3N_d(r)}{\prod_{j=1}^d(r+j)^3}.       \tag{1.5}
\]
The denominator is nonzero on the nonwrapping interval, so (1.5) recovers the
banked continuant collision criterion exactly.  In the transposed convention
used in `chatgpt-answers/Q6520.md`, the same statement is that the upper-right
entry of the standard solution-state transfer vanishes.

### 1.4 The proposed two-run rigidity fails on the live orbit

At \(p=997\), \(d=182\), both \(r=248\) and \(r=565\) start genuine
two-runs: \(r,r+1\in Z_{182}\).  The following exact values are computed
directly from the recurrence.  Here \(v\) is the point proposed in the
specification,
\[
v(r,d)=(\alpha_r-\alpha_{r+d}:\beta_r-\beta_{r+d}).
\]

| \(r\) | \(\xi_r\) | projective \(T_{r,d}\xi_r\) | \(v(r,d)\) | \(T_{r,d}\) |
|---:|:---:|:---:|:---:|:---:|
| 248 | 758 | 344 | 798 | \(\begin{psmallmatrix}929&174\\0&297\end{psmallmatrix}\) |
| 565 | 409 | 780 | 165 | \(\begin{psmallmatrix}355&66\\0&323\end{psmallmatrix}\) |

Both lower-left entries are zero, as (1.3) requires.  Neither matrix fixes
the displayed \(\xi_r\), and neither proposed \(v(r,d)\) equals \(\xi_r\).
Thus the fresh run-rigidity statement in the supplied context, and the
dependent `CODEX_SPEC_laststand_norun.md`, are refuted as stated.

### 1.5 Two target quantities were conflated

For the integer shell
\[
B_D=\{d:\lfloor D/2\rfloor<d\le D\},
\]
write
\[
m_r=\#\{d\in B_D:r\in Z_d\},
\]
\[
S_B=\sum_r m_r=\sum_{d\in B_D}|Z_d|,
\qquad
U_B=\#\{r:m_r>0\}=\left|\bigcup_{d\in B_D}Z_d\right|.
\]
The specification's displayed count is \(S_B\); the earlier `[FR_eta]`
formulation is \(U_B\), or a primitive first-return count after dyadic
decomposition.  One only has
\[
U_B\le S_B.                                       \tag{1.6}
\]
The restart multiplicity cap gives a converse with a growing factor, not the
claimed fixed-power equivalence.

There is also an implicit scale condition.  The energy implication uses a
sufficiently slow unbounded \(L\), so below I use the intended regime
\[
D=\sqrt N L,\qquad L\to\infty,\quad L=N^{o(1)},\quad D=o(N).
\]
Without a restriction such as \(D<N\), the literal phrase "some unbounded
\(L\)" permits an empty shell beyond the observation window and is not the
intended theorem.

## 2. R1: exact horizontal dispersion and its endpoint

Let \(J\) be any interval of \(m\) consecutive gap heights and put
\[
S_J=\sum_{h\in J}|Z_h|,qquad
k_J(r)=\#\{h\in J:r\in Z_h\}.
\]

### 2.1 Exact difference charging

If \(h<k\) and \(r\in Z_h\cap Z_k\), then
\[
\xi_{r+h}=\xi_r=\xi_{r+k}.
\]
Consequently \(r+h\in Z_{k-h}\), and translation by \(h\) is injective.
Using the degree budget gives
\[
|Z_h\cap Z_k|\le C_{k-h}\le3(k-h-1).             \tag{2.1}
\]
Therefore the exact shell pair mass satisfies
\[
\begin{aligned}
I_{2,J}
&:=\sum_r\binom{k_J(r)}2\\
&=\sum_{h<k\in J}|Z_h\cap Z_k|\\
&\le3\sum_{g=1}^{m-1}(m-g)(g-1)
=3\binom m3
=\frac{m(m-1)(m-2)}2.                            \tag{2.2}
\end{aligned}
\]
In particular,
\[
\sum_r k_J(r)^2\le S_J+m(m-1)(m-2).              \tag{2.3}
\]

Truncating the difference at \(G<m\) gives the exact polynomial
\[
I_{2,J}^{(\le G)}
\le3\sum_{g=1}^{G}(m-g)(g-1)
=\frac{G(G-1)(3m-2G-2)}2.                         \tag{2.4}
\]
A dyadic layer \(U<g\le V\) is the difference of two expressions in
(2.4), hence has size \(O(mV^2)\).  Summing dyadic layers is dominated by
\(V\asymp m\) and returns \(O(m^3)\).  The locality in \(|h-k|\) creates no
logarithmic or power saving.

Global Cauchy applied to (2.3) yields only
\[
S_J^2\le N\{S_J+m(m-1)(m-2)\},
\]
and hence
\[
S_J\le
\frac{N+\sqrt{N^2+4Nm(m-1)(m-2)}}2.              \tag{2.5}
\]

### 2.2 `[DYADIC-LINNIK]`: the stronger sliding-window calculation

Fix \(1\le H\le m\).  Slide a length-\(H\) window across \(J\), including
the nonempty boundary truncations.  There are \(m+H-1\) such windows, every
height is counted \(H\) times, and a height pair at difference \(g<H\) is
counted \(H-g\) times.

For every window \(W\), the elementary inequality
\[
k\le\mathbf1_{k>0}+\binom k2
\]
gives
\[
\sum_rk_W(r)\le N+\sum_r\binom{k_W(r)}2.
\]
Summing the sliding windows and applying (2.1) gives
\[
HS_J\le(m+H-1)N
+3\sum_{g=1}^{H-1}(H-g)(m-g)(g-1).               \tag{2.6}
\]
The finite sum is exactly
\[
\sum_{g=1}^{H-1}(H-g)(m-g)(g-1)
=\frac{H(H-1)(H-2)(2m-H-1)}{12}.                 \tag{2.7}
\]
Thus
\[
S_J\le
\frac{m+H-1}{H}N
+\frac{(H-1)(H-2)(2m-H-1)}4.                     \tag{2.8}
\]
Taking \(H=\lceil N^{1/3}\rceil\), with harmless endpoint adjustment when
\(H>m\), gives
\[
S_J\ll N+mN^{2/3}.                                \tag{2.9}
\]
Combined with the degree sum,
\[
S_J\ll\min\{D^2,\,N+DN^{2/3}\}.                 \tag{2.10}
\]
At \(D=\sqrt N L\), the second branch divided by \(D^2\) contains
\(N^{1/6}/L\).  For the intended \(L=N^{o(1)}\), (2.10) therefore reverts to
the \(D^2\) branch.  The sliding calculation is exact but does not reach
`[FR_eta]`.

### 2.3 Fixed-base difference closure is false

The true cascade statement is
\[
a<b\in\operatorname{Gap}(r)
\quad\Longrightarrow\quad
b-a\in\operatorname{Gap}(r+a),                   \tag{2.11}
\]
not \(b-a\in\operatorname{Gap}(r)\).  If the positions of one projective
fibre are
\[
x_0<x_1<\cdots,
\]
then exactly
\[
\operatorname{Gap}(x_i)=\{x_j-x_i:j>i\}.         \tag{2.12}
\]
For the fibre \(\{0,6,10\}\),
\[
\operatorname{Gap}(0)=\{6,10\},
\]
so \(4\notin\operatorname{Gap}(0)\), although
\(4\in\operatorname{Gap}(6)\).  The base shift in (2.11) destroys the
claimed arithmetic-progression modulus, numerical semigroup, and power
amplification.

### 2.4 `[PRIVATE-SINGLETON-OBSTRUCTION]`

The failure is structural for the inputs used in R1.  Suppose \(D^2\ge N\)
and \(D=o(N)\).  Let
\[
t=\left\lfloor\frac{N}{2D}\right\rfloor,
\qquad d_j=D-j\quad(0\le j<t).
\]
Since \(t\le D/2\), every \(d_j\) lies in \((D/2,D]\).  Place \(t\)
disjoint blocks of lengths \(2d_j\).  Within block \(j\), pair its \(i\)-th
position with its \((d_j+i)\)-th position using a colour private to that one
pair.  Give every unused position a fresh colour.  The blocks fit because
\[
2\sum_jd_j\le2tD\le N.
\]

Then
\[
|Z_{d_j}|=d_j,qquad Z_h=\varnothing
\quad\text{for all other shell heights},
\]
the nonempty \(Z_{d_j}\) are pairwise disjoint, and
\[
S_J=U_J
=\sum_{j=0}^{t-1}(D-j)
=tD-\frac{t(t-1)}2.                               \tag{2.13}
\]
Every active base has multiplicity one, so all intersections, difference
charges, and higher cascade moments vanish.  Nevertheless
\(C_{d_j}=d_j\le3(d_j-1)\), and the collision equivalence and row-restart
constraints hold.

If \(D=\sqrt N L\), \(L\to\infty\), \(L=N^{o(1)}\), and \(D=o(N)\), then
\[
S_J=U_J\sim\frac N2
=\frac{D^2}{2L^2}=D^{2-o(1)}.                     \tag{2.14}
\]
It violates \(O(D^{2-\eta})\) for every fixed \(\eta>0\).  This is not an
Apéry counterexample.  It proves that the R1 incidence, intersection,
cascade, and restart inputs cannot imply the target.  The verifier constructs
one literal instance
\[
(N,D,t,S_J)=(100000,3641,13,47255).
\]

## 3. R2: shifted correlations are not the missing first moment

### 3.1 Exact renewal correlation

The continuants satisfy the addition law
\[
\begin{aligned}
N_{a+g}(X)
={}&N_g(X+a)N_{a+1}(X)\\
&-(X+a+1)^6N_{g-1}(X+a+1)N_a(X).                 \tag{3.1}
\end{aligned}
\]
It follows by induction in \(g\) from (1.4); the verifier expands (3.1)
exactly through a grid of independent instances.

Define
\[
\tau_{-a}Z_g=\{r:r+a\in Z_g\}.
\]
At a physical root \(N_a(r)=0\), one has \(N_{a+1}(r)\ne0\): otherwise
\(\xi_r=\xi_{r+a}=\xi_{r+a+1}\), contradicting the nonzero adjacent
Casoratian.  Reducing (3.1) at such a root gives the exact identity
\[
Z_a\cap\tau_{-a}Z_g=Z_a\cap Z_{a+g}.              \tag{3.2}
\]
Both sides count
\[
\kappa_p(a,g)
=\#\{r:\xi_r=\xi_{r+a}=\xi_{r+a+g}\}.           \tag{3.3}
\]

### 3.2 The literal positive-proportion target is already true

Assume \(2D<N\) and take \(a,g\in B_D\).  Every term in (3.3) is a triangle
of span at most \(2D\).  The banked triangle estimate therefore gives
\[
\sum_{a,g\in B_D}\kappa_p(a,g)
\le Q_{2D}
\le264D^2\{1+\log(2D)\}.                          \tag{3.4}
\]
Since \(|B_D|\ge D/2\), the average is at most
\[
1056\{1+\log(2D)\}.
\]
Markov's inequality shows that at least half of all \((a,g)\in B_D^2\)
satisfy
\[
\kappa_p(a,g)\le2112\{1+\log(2D)\}.              \tag{3.5}
\]
This is far below the \(O(D)\) Bezout scale.  Thus the stated intermediate
goal of finding a positive proportion of low-correlation gap pairs is not the
missing lemma; it already follows from a banked bound.

### 3.3 The exact singleton obstruction

For the shell multiplicities from Section 1.5, put
\[
I_B=\sum_r\binom{m_r}{2}.
\]
This pair mass has the exact correlation decomposition
\[
I_B
=\sum_{h<k\in B_D}|Z_h\cap Z_k|
=\sum_{h<k\in B_D}C_p(h,k-h).                    \tag{3.6}
\]
It counts two shell returns from the same base.  This is not the same slice
as \(\kappa_p(a,g)=C_p(a,g)\) with both \(a,g\in B_D\): the latter counts
two consecutive shell edges whose total span lies in \((D,2D]\).  Such a
triangle can meet shell-singletons.  For example, a fibre at positions
\(\{0,a,a+g\}\), with \(a,g\in B_D\), has
\(\kappa_p(a,g)=1\), while the bases \(0\) and \(a\) each have shell
multiplicity one and \(I_B=0\).

This distinction does not make an upper correlation bound control the first
moment.  In the private-pair construction of Section 2.4 every fibre has
only two positions, so every \(C_p(a,g)\) is zero while the singleton mass is
\(D^{2-o(1)}\).  For any integer \(K\ge2\), define
\[
L_B(K)=\sum_{r:1\le m_r<K}m_r.
\]
Since
\[
m\le\frac{2\binom m2}{K-1}\qquad(m\ge K),
\]
one has the exact decomposition
\[
S_B\le L_B(K)+\frac{2I_B}{K-1}.                  \tag{3.7}
\]
In particular, the \(m_r=1\) mass is completely invisible to \(I_B\).
Even the ideal conclusion \(I_B=0\) allows \(S_B=U_B=D^{2-o(1)}\), as the
construction in Section 2.4 shows.  No bound on second-order correlations
alone can prove the required first-order census: the same-base slice misses
singletons identically, while the consecutive-edge slice may see some of
them but can be zero on an arbitrarily large singleton family.

### 3.4 Characteristic-zero resultants have the wrong quantifiers

The natural translated resultant is
\[
E_{a,g}=\operatorname{Res}_X\bigl(N_a(X),N_g(X+a)\bigr).
\]
The banked separated-root argument proves \(E_{a,g}\ne0\) over
characteristic zero.  For each fixed \((a,g)\), this excludes all but finitely
many primes.  The present problem instead fixes \(p\) and lets
\(a,g\asymp D(p)\) grow.  Nothing prevents the same working prime from
dividing many different nonzero integers \(E_{a,g}\).  Resultant heights only
help after averaging over primes, which is the wrong quantifier order for a
quenched statement.

The elementary family
\[
f(X)=X+73,qquad g(X)=X+146
\]
already displays the logic: \(f,g\) are coprime over \(\mathbf Q\), their
nonzero resultant has absolute value \(73\), and both reduce to \(X\) modulo
\(73\).

Two Apéry-specific specialization failures also matter:

1. The "honest triple poles" formula is a nonzero-integer statement, not a
   nonzero-mod-every-prime statement.  Exactly
   \[
   N_3(-3)=584=8\cdot73,
   \]
   and
   \[
   \gcd_{\mathbf F_{73}}(N_3,N_4)=X+3.
   \]
   Hence the claim in `Q6730.md` that these pole evaluations remain nonzero
   modulo the working prime is false; saturation is required.
2. Squarefreeness over \(\mathbf Q\) does not imply squarefreeness at the
   working prime.  The exact live specialization
   \[
   \gcd_{\mathbf F_{211}}(N_{32},N_{32}')
   =X^2+33X+114=(X-89)^2
   \]
   has the physical multiple root \(r=89\).

The addition law is valuable, but its direct norm/resultant expansion does
not close in a fixed-dimensional scalar state: mixed polarized norms
proliferate with \(\deg N_a=3a-3\).  The full-rank tests in
`CRON_LOWRANK_REPORT.md` are consistent with this obstruction.  They are
negative evidence, not an impossibility theorem for every conceivable
arithmetic structure.

## 4. R3: exact failure of the proposed adiabatic mechanism

### 4.1 Adjacent changes are always \(p\)-adic units in some coordinate

Direct subtraction gives
\[
\alpha_{s+1}-\alpha_s
=\frac{51s^4+252s^3+435s^2+306s+77}
{(s+1)^3(s+2)^3},                                 \tag{4.1}
\]
\[
\beta_{s+1}-\beta_s
=\frac{3s^4+12s^3+15s^2+6s+1}
{(s+1)^3(s+2)^3}.                                 \tag{4.2}
\]
The two numerator quartics have exact resultant
\[
248832=2^{10}3^5.                                 \tag{4.3}
\]
Thus, for every \(p\ge5\) and every regular
\(s\not\equiv-1,-2\pmod p\), the two changes cannot both vanish modulo
\(p\).  Their common denominator is a \(p\)-adic unit, so in the entrywise
sup norm
\[
\|M_{s+1}-M_s\|_p=1.                              \tag{4.4}
\]

For \(M(a,b)=\begin{psmallmatrix}0&1\\-b&a\end{psmallmatrix}\), an exact
multiplication gives
\[
[M(a_s,b_s),M(a_t,b_t)]
=\begin{pmatrix}
b_s-b_t&a_t-a_s\\
a_tb_s-a_sb_t&b_t-b_s
\end{pmatrix}.                                    \tag{4.5}
\]
Equations (4.1)--(4.3) show that the first row in (4.5) is nonzero when
\(t=s+1\).  Adjacent physical matrices never commute.

Archimedeanly, the changes in fact have size \(O(s^{-2})\), slightly better
than the \(O(s^{-1})\) suggestion in the specification.  That observation
does not descend to \(\mathbf F_p\): distinct residues have unit-separated
Teichmuller lifts, and (4.3)--(4.4) give an exact obstruction to treating the
raw matrices as a Hensel-small perturbation of one frozen matrix.

### 4.2 `[NO-P-ADIC-ADIABATIC-EIGENFRAME]`

The characteristic polynomial of \(M_s\) is
\[
\chi_s(z)=z^2-\alpha_sz+\beta_s.
\]
Its discriminant factors as
\[
\alpha_s^2-4\beta_s
=\frac{(3s^2+3s+1)
(384s^4+768s^3+579s^2+195s+25)}{(s+1)^6}.         \tag{4.6}
\]
For \(p\ge5\), there are at most six parabolic residues.  The resultant of
two adjacent characteristic polynomials is
\[
\operatorname{Res}_z(\chi_s,\chi_{s+1})
=-\frac{24R_8(s)}{(s+1)^6(s+2)^6},                \tag{4.7}
\]
where
\[
\begin{aligned}
R_8(s)={}&108s^8+864s^7+2865s^6+5094s^5\\
&+5253s^4+3180s^3+1095s^2+198s+16.
\end{aligned}
\]
This nonzero degree-eight polynomial has at most eight roots.  After allowing
six exceptional residues for each endpoint and eight for a common root,
outside at most twenty residues the adjacent matrices are both separable and
share no eigenline over \(\overline{\mathbf F}_p\).  Equivalently, after
choosing their two eigenbases, every transition coefficient is nonzero modulo
\(p\); a zero coefficient would identify one old and one new eigenline.

Equations (4.1)--(4.7) prove that the advertised locally commuting,
near-diagonal, slowly rotating eigenline model is not a valid approximation
to this cocycle.  They do **not** prove long-product mixing.  A deterministic
sequence of locally transverse matrices may still have exceptional ordered
products.  The correct R3 object remains the quenched count
\[
\#\{(r,d):d\in B_D,\ (T_{r,d})_{21}=0\},          \tag{4.8}
\]
or an equivalent special-Borel-observable local spectral-flatness estimate.
An ordinary full-space Koopman spectral gap is impossible for the exact
time-preserving skew product, as already analyzed in
`chatgpt-answers/Q6520.md`.

## 5. Exact numerical calibration

For each prime below the verifier constructs the full orbit on
\(1,\ldots,p-2\), groups equal projective values, and enumerates all
nonwrapping pairs in
\[
D=\lceil\sqrt p\log p\rceil,qquad D/2<d\le D.
\]
Here \(S_B\) is the incidence count, \(U_B\) is the number of active bases,
"singleton" counts bases with \(m_r=1\), and
\(I_2=\sum_r\binom{m_r}{2}\).  The cutoff is computed with 80-decimal-digit
arithmetic, and every displayed \(p\) is checked prime.

| \(p\) | \(D\) | \(S_B\) | \(U_B\) | singleton | \(\max m_r\) | \(I_2\) | multiplicity histogram |
|---:|---:|---:|---:|---:|---:|---:|:---|
| 499 | 139 | 83 | 73 | 63 | 2 | 10 | `1:63, 2:10` |
| 997 | 219 | 161 | 149 | 138 | 3 | 13 | `1:138, 2:10, 3:1` |
| 1999 | 340 | 253 | 232 | 213 | 3 | 23 | `1:213, 2:17, 3:2` |
| 4001 | 525 | 371 | 350 | 332 | 3 | 24 | `1:332, 2:15, 3:3` |
| 7919 | 799 | 616 | 590 | 564 | 2 | 26 | `1:564, 2:26` |
| 16001 | 1225 | 862 | 837 | 812 | 2 | 25 | `1:812, 2:25` |
| 32003 | 1856 | 1400 | 1370 | 1340 | 2 | 30 | `1:1340, 2:30` |
| 65537 | 2840 | 2090 | 2055 | 2021 | 3 | 36 | `1:2021, 2:33, 3:1` |
| 99991 | 3641 | 2692 | 2647 | 2602 | 2 | 45 | `1:2602, 2:45` |

This finite census is consistent with a linear-size truth, but it is not an
asymptotic proof.  Its relevant diagnostic is more specific: throughout the
table, most active bases are singletons and \(I_2\) is tiny compared with
\(S_B\).  That is exactly the stratum ignored by R1/R2 correlation bounds.

## 6. Inventory of unconditional bricks

| label | status | content | proof / machine gate |
|:---|:---:|:---|:---|
| `[FRAME-BOREL]` | PROVED | \(\xi_{r+d}=\xi_r\iff(T_{r,d})_{21}=0\), with the exact continuant entry (1.5) | Section 1.3; exhaustive small-prime frame and transfer gate |
| `[FALSE-XI-ACTION]` | REFUTED | The specification's \(\xi_{r+1}=M_r\xi_r\) and dependent run eigenvector are false | Sections 1.2, 1.4; rational and live \(p=997\) witnesses |
| `[DYADIC-LINNIK]` | PROVED | Exact sliding estimate (2.8), hence \(S_J\ll N+mN^{2/3}\) | Sections 2.1--2.2; all finite-sum identities checked |
| `[FALSE-GAP-AP]` | REFUTED | Return gaps are not difference-closed at a fixed base | Section 2.3; fibre \(\{0,6,10\}\) gate |
| `[PRIVATE-SINGLETON-OBSTRUCTION]` | PROVED | R1/R2 inputs permit \(S_J=U_J=D^{2-o(1)}\) with every active row singleton | Section 2.4; literal private-block constructor |
| `[SHIFTED-RENEWAL-CORRELATION]` | PROVED | Exact equality (3.2); at least half of dyadic shifted gap pairs obey (3.5) | Sections 3.1--3.2; symbolic renewal gates and banked \(Q_{2D}\) bound |
| `[FIXED-P-SPECIALIZATION-CAVEAT]` | PROVED | Characteristic-zero coprimality, squarefreeness, and nonzero pole values do not give the needed live-prime statements | Section 3.4; exact \(p=73,211\) certificates |
| `[NO-P-ADIC-ADIABATIC-EIGENFRAME]` | PROVED | Adjacent raw matrices are unit-separated/noncommuting; outside at most twenty residues they share no eigenline | Section 4; exact symbolic resultants and finite-field gates |

None of these bricks proves `[FR_eta]`.  The first one corrects the attack
surface; the others identify precisely why the proposed reductions stop.

## 7. Sharp surviving conditional statement

The banked triangle bound controls the shell pair mass by
\[
I_B\le Q_D\le66D^2(1+\log D).                    \tag{7.1}
\]
Indeed, by (3.6),
\[
I_B=\sum_{h<k\in B_D}C_p(h,k-h)\le Q_D.          \tag{7.1}
\]
The banked estimate gives
\[
I_B\le Q_D\le66D^2(1+\log D).                    \tag{7.2}
\]
Combining (7.2) with the exact decomposition (3.7) isolates one sufficient
new input.

> **`[LOWMULT(delta,epsilon)]`.**  For some fixed
> \(\delta,\epsilon>0\), with
> \(K=\lceil D^\delta\rceil\), the actual Apéry orbit satisfies
> \[
> L_B(K)=\sum_{r:1\le m_r<K}m_r\ll D^{2-\epsilon}
> \]
> uniformly at the intended mesoscopic scale.

Indeed, (3.7) and (7.2) then give
\[
S_B
\ll D^{2-\epsilon}
+\frac{D^2(1+\log D)}{D^\delta-1}
=D^{2-\min(\delta,\epsilon)+o(1)}.                \tag{7.3}
\]
Thus `[FR_eta]` follows for every
\(\eta<\min(\delta,\epsilon)\).  Unlike another second-moment or resultant
estimate, `[LOWMULT]` directly attacks the primitive/single-return mass.  In
the corrected cocycle language it is a quenched low-multiplicity census for
zeros of the ordered coefficient \((T_{r,d})_{21}\).

The private-pair obstruction proves that no combination of the R1/R2
incidence, degree, intersection, cascade, and restart inputs can supply this
condition.  The R3 local transversality identities do not supply it either.
It must come from a new growing-gap arithmetic or special-observable spectral
estimate retaining the ordered clock.

## 8. Audit provenance and reproducibility

The independent R1 consultation correctly located the singleton obstruction;
its initially suggested fully disjoint \(\asymp d\)-sized family was resized
to the capacity-correct block construction (2.13).  The independent R3
consultation correctly rejected archimedean-to-\(p\)-adic slow variation, but
it repeated the false \(T\xi=\xi\) formulation; only the frame/Borel version
proved here was retained.  `chatgpt-answers/Q6520.md` already contained the
correct frame formulation and served as a cross-check, not as a substitute
for the live verifier.

Run

```text
python3 CODEX_LASTSTAND_COCYCLE_verify.py
```

from this directory.  The terminal line is

```text
FINAL GATE: PASS -- all claims in the terminal report were reproduced
```

The verifier uses exact integer, rational, finite-field, and SymPy polynomial
arithmetic.  It independently checks primality, the high-precision numerical
cutoffs, all displayed shell counts, the two live cocycle counterexamples,
the private-block construction, and every displayed resultant or
specialization value.

## 9. Required closing answers

1. **Least-confident step.**  The exact algebra and counting identities are
   fully checked.  The least certain part is only the strategic word
   "sharp": `[LOWMULT]` is the cleanest sufficient condition exposed by the
   exact decomposition, but the current audit cannot exclude a different new
   arithmetic invariant which bypasses that decomposition.
2. **Blind spot in the question.**  The main blind spot is the conflation of
   the solution-state cocycle with the row-projective Apéry orbit.  Secondary
   blind spots are the incidence/union conflation, the missing slow-scale
   restriction on \(L\), and the assumption that nonzero characteristic-zero
   pole/resultant data remain nonzero at the fixed working prime.
