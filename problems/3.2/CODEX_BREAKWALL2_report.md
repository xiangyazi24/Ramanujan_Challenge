NOT BROKEN, NEW RESIDUAL [MIRROR-WL1]

# Breakwall strike 2

## 1. Verdict

No win condition W1--W7 is proved for the Apéry orbit.  The wall is therefore
not broken.  The stopping condition is instead met by a new, named sufficient
condition which is strictly weaker than `[ZERO-TAIL-2]`:

\[
 [\mathrm{BDH\!\!-LAG}]
 \Longrightarrow[\mathrm{ZERO\!\!-TAIL\!\!-2}]
 \Longrightarrow_{\rm strict}[\mathrm{MIRROR\!\!-WL1}]
 \Longrightarrow_{\rm strict} W1.
\]

The new face is a weak-\(L^1\) tail condition on the non-forced reflection
pairs.  In particular it allows a single lag to attain its full degree bound;
that same profile kills `[BDH-LAG]`, `[SAME-LAG-L2]`, and
`[ZERO-TAIL-2]`.  The strictness witness below is an actual
reflection-symmetric word satisfying the degree cap, the banked window bound,
the energy bound, and \(Q_D=0\).  It is a deduction countermodel, not a
counterexample to W1 for the Apéry orbit.

The three prescribed attack directions have the following terminal status.

1. Slice disjointness proves all-pair characteristic-zero coprimality, but its
   fixed-prime intersection count is exactly the already banked shifted
   three-point quantity \(Q_D\).  Resultant height gives no quenched estimate.
2. Renewal gives an exact restart equivalence, not propagation.  An actual
   row at \(p=461\) has eleven roots and no continuation, so the proposed
   bootstrap is false.
3. Reflection gives an exact decomposition
   \(R_h=\kappa_h+2\mu_h+\varepsilon_h\).  The endpoint term
   \(\varepsilon_h\), omitted by the naive parity statement, is essential.
   This decomposition yields `[MIRROR-WL1]`.

Two elementary distributional faces, `[ZERO-EXCESS-1/2]` and
`[LAG-LORENTZ-3/2]`, are recorded in Section 8.  They are not used for the
first-line verdict.

## 2. Setup

Put \(N=p-1\), and write the nonwrapping projective orbit as

\[
 \pi(0),\ldots,\pi(N-1),
 \qquad
 \pi(r)=[b_r:c_r]\in\mathbf P^1(\mathbf F_p).
\]

For \(1\le h\le D\le N-2\), let

\[
 C_h=\{0\le r\le N-h-1:\pi(r)=\pi(r+h)\},
 \qquad R_h=|C_h|,
\]

and let \(S_D=\sum_{h\le D}R_h\).  We use the banked facts

\[
 R_h\le 3(h-1),
 \qquad \pi(N-n)=\pi(n)
\]

where the reflection identity is understood on the regular extended window
\(0\le n\le N\).  At the hatch scale

\[
 D=\sqrt N\,L,
 \qquad L\longrightarrow\infty,
 \qquad L=N^{o(1)},
 \qquad
 T_0=\left\lceil\sqrt{N/D}\right\rceil.
\]

The banked window multiplicity estimate will be used only in the form

\[
 m_J(v)-1\le 4|J|^{2/3}.                                      \tag{2.1}
\]

## 3. Exact mirror skeleton, including the endpoint

### Theorem 3.1 (physical mirror decomposition)

Define

\[
 \kappa_h=\mathbf 1_{2\mid h},
 \qquad
 \varepsilon_h=\mathbf 1_{\pi(h)=\pi(0)}.
\]

There is a nonnegative integer \(\mu_h\) such that, with no exceptional lag,

\[
 \boxed{R_h=\kappa_h+2\mu_h+\varepsilon_h.}                    \tag{3.1}
\]

Moreover

\[
 B_D:=\sum_{h\le D}\varepsilon_h
 =m_{[0,D]}(\pi(0))-1
 \le 4(D+1)^{2/3}.                                             \tag{3.2}
\]

#### Proof

Remove the possible base \(r=0\) and put

\[
 C_h^\circ=C_h\cap[1,N-h-1],
 \qquad \iota_h(r)=N-h-r.
\]

If \(r\in C_h^\circ\), strong reflection gives

\[
 \pi(\iota_h(r))
 =\pi(N-(r+h))
 =\pi(r+h)
 =\pi(r)
 =\pi(N-r)
 =\pi(\iota_h(r)+h).
\]

Also \(1\le\iota_h(r)\le N-h-1\).  Thus \(\iota_h\) is an
involution of \(C_h^\circ\).  Its fixed-point equation is
\(2r=N-h\).  Since \(N\) is even, a fixed point exists exactly when \(h\)
is even, in which case

\[
 r_h=(N-h)/2\in[1,N-h-1]
\]

and it is a collision because \(r_h+h=N-r_h\).  Hence the number of
fixed orbits is exactly \(\kappa_h\).  If \(\mu_h\) is the number of the
remaining two-element orbits, then
\(|C_h^\circ|=\kappa_h+2\mu_h\).  Adding the base \(0\) proves
(3.1).

The endpoint is unique: \(\iota_h(0)=N-h\), exactly one past the largest
allowed base \(N-h-1\).  Finally, \(\varepsilon_h=1\) precisely when the
color \(\pi(0)\) returns at index \(h\).  Equation (3.2) is therefore (2.1)
on \([0,D]\). \(\square\)

Consequently the exact physical parity law is

\[
 R_h\equiv\kappa_h+\varepsilon_h\pmod 2,                       \tag{3.3}
\]

not \(R_h\equiv\kappa_h\pmod2\) without qualification.

### Proposition 3.2 (center and pole bookkeeping)

The gap numerators satisfy

\[
 N_h(-h-1-X)=(-1)^{h-1}N_h(X).                                \tag{3.4}
\]

For even \(h\), the fixed point \(-(h+1)/2\) is not a pole and

\[
 (2X+h+1)\mid N_h(X).                                         \tag{3.5}
\]

Its residue is the physical central base \((p-1-h)/2\), giving the term
\(\kappa_h\) in (3.1).  For odd \(h\), put \(j=(h+1)/2\).  The fixed
point is the pole \(-j\), and the banked pole formula specializes to

\[
 |N_h(-j)|=((j-1)!)^6b_{j-1}^2.                               \tag{3.6}
\]

Even if \(p\mid b_{j-1}\), its residue \(p-j\) lies strictly above the
physical base bound \(p-2-h\).  Thus an odd-lag fixed pole may enter the
full polynomial root census but never changes the physical parity (3.3).
For example, at \(p=17,h=7\), the residue \(13\) is a root of \(N_7\),
but the physical interval is \([0,8]\) and \(R_7=0\).

#### Proof

Equation (3.4) is the banked continuant reflection identity (and is checked
coefficientwise by the verifier).  At an even-lag fixed point its sign is
\(-1\), so characteristic zero forces a zero; the point is a half-integer
and hence is not among the integer poles \(-1,\ldots,-h\).  This proves
(3.5).  At odd \(h\), the sign is \(+1\), so reflection forces no zero;
(3.6) is the pole-value formula with \(h-j=j-1\).  Finally

\[
 (p-j)-(p-2-h)=h+2-j=(h+3)/2>0,
\]

which proves the physical exclusion. \(\square\)

### Corollary 3.3 (the stripped variance is pure pair fluctuation)

Let

\[
 V_D=\sum_{h\le D}(R_h-\kappa_h)^2.
\]

Then

\[
 4\sum_{h\le D}\mu_h^2
 \le V_D
 \le 8\sum_{h\le D}\mu_h^2+2B_D,                            \tag{3.7}
\]

and

\[
 \sum_{h\le D}R_h^2
 \le 3\sum_{h\le D}\kappa_h
      +12\sum_{h\le D}\mu_h^2+3B_D.                         \tag{3.8}
\]

Thus `[BDH-LAG]` is, up to absolute constants, precisely

\[
 \sum_{h\le D}\mu_h^2\ll\sqrt{DN},                          \tag{3.9}
\]

and `[SAME-LAG-L2]` is equivalent to \(\sum\mu_h^2\ll N\).

#### Proof

From (3.1), \(R_h-\kappa_h=2\mu_h+\varepsilon_h\).  Since
\(\varepsilon_h\in\{0,1\}\),

\[
 4\mu_h^2\le(2\mu_h+\varepsilon_h)^2
 \le8\mu_h^2+2\varepsilon_h.
\]

Also \((x+y+z)^2\le3x^2+3y^2+3z^2\), applied to
\(\kappa_h+2\mu_h+\varepsilon_h\), gives (3.8).  The reverse estimate
\(R_h^2\ge4\mu_h^2\) is immediate.  At the hatch scale,
\(B_D=O(D^{2/3})=o(\sqrt{DN})\) and \(D+B_D=o(N)\), proving both
equivalences. \(\square\)

## 4. The new bottom face [MIRROR-WL1]

Put

\[
 A_D^\mu(t)=\#\{h\le D:\mu_h\ge t\},
 \qquad
 U_D=\left\lfloor\frac{3(D-1)}2\right\rfloor,
\]

and define

\[
 \mathcal W_D^\mu=\max_{\substack{t\in\mathbf Z\\t>T_0}}
                         tA_D^\mu(t),
 \qquad
 H_D=\sum_{t=T_0+1}^{U_D}\frac1t.                             \tag{4.1}
\]

An empty maximum or sum is zero.

> **[MIRROR-WL1].**
> \[
> \boxed{\mathcal W_D^\mu H_D\ll N.}                         \tag{4.2}
> \]

When only the weak-\(L^1\) supremum in (4.1) is used, the harmonic factor is
the sharp generic layer-cake cost.  This is not an \(L^2\) condition.

### Theorem 4.1 (weak-\(L^1\) bridge)

`[MIRROR-WL1]` implies W1.  More precisely,

\[
 S_D\le \frac D2+4(D+1)^{2/3}
             +2DT_0+2\mathcal W_D^\mu H_D.                   \tag{4.3}
\]

#### Proof

The degree bound and (3.1) give \(\mu_h\le U_D\).  Integer layer cake gives

\[
 \begin{aligned}
 \sum_{h\le D}\mu_h
 &=\sum_{t=1}^{U_D}A_D^\mu(t)\\
 &\le DT_0+
 \sum_{t=T_0+1}^{U_D}\frac{\mathcal W_D^\mu}{t}\\
 &=DT_0+\mathcal W_D^\mu H_D.
 \end{aligned}                                                \tag{4.4}
\]

Now sum (3.1), use \(\sum\kappa_h\le D/2\), and apply (3.2).  This is
(4.3).  Finally

\[
 DT_0\le\sqrt{ND}+D
 =N^{3/4}L^{1/2}+N^{1/2}L=o(N),                               \tag{4.5}
\]

while \(D+D^{2/3}=o(N)\).  Condition (4.2) therefore gives
\(S_D\ll N\), which is W1. \(\square\)

### Theorem 4.2 (`[ZERO-TAIL-2]` implies `[MIRROR-WL1]`)

Let

\[
 \mathcal Z_D=\max_{t>T_0}t^2\#\{h\le D:R_h\ge t\}.
\]

Then

\[
 \mathcal W_D^\mu H_D
 \le \frac{\mathcal Z_DH_D}{4(T_0+1)}.                       \tag{4.6}
\]

Consequently \(\mathcal Z_D\ll N\) implies (4.2), with an extra
\(o(1)\) factor.

#### Proof

If \(\mu_h\ge t\), then (3.1) gives \(R_h\ge2t\).  Hence, for integral
\(t>T_0\),

\[
 tA_D^\mu(t)
 \le t\#\{h:R_h\ge2t\}
 \le\frac{\mathcal Z_D}{4t}
 \le\frac{\mathcal Z_D}{4(T_0+1)}.
\]

Taking the maximum and multiplying by \(H_D\) proves (4.6).  Since
\(H_D=O(\log N)\) and
\(T_0=N^{1/4-o(1)}\), one has \(H_D/T_0=o(1)\). \(\square\)

## 5. Strictness and the reflection-word death certificate

### Proposition 5.1 (one full lag passes `[MIRROR-WL1]`)

There are reflection-symmetric words obeying all of the following at once:

\[
 R_h\le3(h-1),\qquad Q_D=0,\qquad
 m_J(v)-1\le4|J|^{2/3},\qquad
 \sum_v m(v)^2=O(N),                                         \tag{5.1}
\]

but with one odd \(d\asymp D\) satisfying

\[
 R_d=3(d-1),\qquad
 V_D\asymp D^2,qquad
 \mathcal Z_D\asymp D^2.                                    \tag{5.2}
\]

Nevertheless these words satisfy `[MIRROR-WL1]` with

\[
 \mathcal W_D^\mu H_D=O(D\log N)=o(N).                       \tag{5.3}
\]

Thus `[MIRROR-WL1]` is strictly weaker than `[ZERO-TAIL-2]`, as well as
strictly weaker than the two \(L^2\) faces.

#### Construction and proof

Take even \(N>9D\), choose an odd \(d\in[D-2,D]\), and put

\[
 M=\frac{3(d-1)}2.
\]

On the extended indices \(0,\ldots,N\), first give every reflection orbit
\(\{n,N-n\}\) a private color.  Let

\[
 \mathcal A_1=\{1,\ldots,d\},
 \qquad
 \mathcal A_2=
 \{2d+1,\ldots,2d+(d-3)/2\}.
\]

There are \(M\) sources in \(\mathcal A=\mathcal A_1\cup\mathcal A_2\),
and the \(2M\) points \(a,a+d\) are all distinct and smaller than \(4d\).
For every \(a\in\mathcal A\), merge the reflection orbits of \(a\) and
\(a+d\).  The resulting color class is

\[
 \{a,a+d,N-a-d,N-a\}.                                       \tag{5.4}
\]

Finally restrict the word to \(0,\ldots,N-1\).

Each class (5.4) contributes exactly two short edges of length \(d\).  All
of its other differences exceed \(D\), because
\(N-8d>D\).  The untouched reflection pairs contribute exactly the one
central edge at each even \(h\le D\).  Therefore

\[
 R_d=2M=3(d-1),\qquad
 R_h=\kappa_h\quad(h\ne d),qquad
 \varepsilon_h=0.                                            \tag{5.5}
\]

This proves the degree cap and (5.2).  The short edges have disjoint forward
bases, so every per-base short degree is at most one and \(Q_D=0\).  Every
color class has size at most four, proving the window assertion in (5.1).
Before the merges the color energy is \(2N-2\); each merge replaces two
classes of size two by one class of size four and adds eight.  Thus

\[
 \sum_vm(v)^2=2N-2+8M=2N+12(d-1)-2=O(N).                     \tag{5.6}
\]

Only \(\mu_d=M\) is nonzero.  Hence
\(\mathcal W_D^\mu=M=O(D)\) and \(H_D=O(\log N)\), proving
(5.3).  At the hatch scale, \(D\log N=N^{1/2+o(1)}=o(N)\), while

\[
 \mathcal Z_D\ge R_d^2\asymp D^2=NL^2
\]

is not \(O(N)\) because \(L\to\infty\). \(\square\)

This construction also proves that the exact mirror skeleton, the window
bound, the global energy bound, and even \(Q_D=0\) cannot imply the stripped
variance bound.  The missing input remains arithmetic anti-concentration for
the growing family \(N_h\bmod p\).

### Proposition 5.2 (`[MIRROR-WL1]` is not W1 renamed)

At the profile level, W1 does not imply `[MIRROR-WL1]`, even with the degree
cap and mirror parity.  Indeed, choose \(A=\lfloor\sqrt N\rfloor\) lags in
\([D/2,D]\), put \(\varepsilon_h=0\) and
\(\mu_h=\lfloor\sqrt N/2\rfloor\) on them, and put \(\mu_h=0\) elsewhere.
Then

\[
 S_D=\Theta(N)+O(D)=O(N),
\]

and the degree cap holds for all sufficiently large \(N\).  But at
\(t\asymp\sqrt N\),

\[
 \mathcal W_D^\mu\asymp N,
 \qquad H_D\asymp\log N,
\]

so (4.2) fails by a logarithmic factor. \(\square\)

## 6. Direction A: slice disjointness and resultants

Let

\[
 P(U)=34U^3+51U^2+27U+5,
\]

and normalize

\[
 N_0=0,\quad N_1=1,\quad
 N_{h+1}(X)=P(X+h)N_h(X)-(X+h)^6N_{h-1}(X).                  \tag{6.1}
\]

### Theorem 6.1 (all-pair characteristic-zero coprimality)

For all \(a,d\ge1\),

\[
 \begin{aligned}
 N_{a+d}(X)
 &=N_{a+1}(X)N_d(X+a)\\
 &\quad -(X+a+1)^6N_a(X)N_{d-1}(X+a+1).                     \tag{6.2}
 \end{aligned}
\]

Moreover

\[
 \gcd_{\overline{\mathbf Q}}
       \bigl(N_a(X),N_d(X+a)\bigr)=1,                        \tag{6.3}
\]

and, for all \(1\le a<c\),

\[
 \gcd_{\mathbf Q}(N_a,N_c)=1.                               \tag{6.4}
\]

In particular every same-base and shifted resultant in (6.3)--(6.4) is a
nonzero integer.

#### Proof

Equation (6.2) is the continuant concatenation identity, after the gauge
change \(N_h(X)=K_{h-1}(X+1)\); it also follows by induction from (6.1).
The banked root-strip theorem places all roots of \(N_a\) in

\[
 -a<\Re X<-1,
\]

whereas the roots of \(N_d(X+a)\) lie in

\[
 -a-d<\Re X<-a-1.
\]

The strips are disjoint, proving (6.3).  Now let \(c=a+d\), and suppose
\(\alpha\) were a common root of \(N_a,N_c\).  Evaluation of (6.2) gives

\[
 N_{a+1}(\alpha)N_d(\alpha+a)=0.
\]

Adjacent continuants are coprime, so \(N_{a+1}(\alpha)\ne0\); (6.3)
excludes \(N_d(\alpha+a)=0\).  This contradiction proves (6.4).
\(\square\)

### Proposition 6.2 (the slice double count is exactly \(Q_D\))

For the physical root sets \(Z_h=C_h\),

\[
 Z_a\cap Z_{a+d}
 =\{0\le r\le N-a-d-1:
        N_a(r)=N_d(r+a)=0\pmod p\}.                           \tag{6.5}
\]

Consequently, if \(d_D(r)=\#\{h\le D:r\in Z_h\}\), then

\[
 \begin{aligned}
 Q_D
 &=\sum_r\binom{d_D(r)}2\\
 &=\sum_{1\le a<c\le D}|Z_a\cap Z_c|\\
 &=\sum_{a+d\le D}
   \#\{r:N_a(r)=N_d(r+a)=0\}.                               \tag{6.6}
 \end{aligned}
\]

#### Proof

If \(r\in Z_a\cap Z_{a+d}\), then
\(\pi(r)=\pi(r+a)=\pi(r+a+d)\), so the shifted collision is in
\(Z_d\).  Conversely the first and shifted collisions imply the third by
transitivity.  All indices in (6.5) are physical, so the transfer matrices
are regular and the collision/\(N_h\)-root equivalence applies.  The first
equality in (6.6) counts unordered pairs of lags at each base; the remaining
equalities merely regroup those pairs. \(\square\)

Thus slice-disjointness has not created a new statistic: its pair correction
is exactly the banked W2 quantity.

### Proposition 6.3 (resultant height does not give a fixed-prime bound)

Put

\[
 \mathscr R_{a,d}=\left|\operatorname{Res}_X
       (N_a(X),N_d(X+a))\right|>0.
\]

Then

\[
 \log\mathscr R_{a,d}=O((a+d)^2\log(a+d)).                   \tag{6.7}
\]

If \(m_{a,d}(p)\) is the physical intersection in (6.5), the strongest
direct divisibility estimate is

\[
 m_{a,d}(p)\le v_p(\mathscr R_{a,d}),                         \tag{6.8}
\]

and hence

\[
 Q_D\ll\frac{D^4\log D}{\log p}.                             \tag{6.9}
\]

At \(D=\sqrt p\,L\), (6.9) is \(O(p^2L^4)\), far above the required
\(O(p)\).

#### Proof

The coefficient \(\ell^1\)-norm recurrence from (6.1) gives
\(\log\|N_h\|_1=O(h\log h)\).  Translation by \(a\) preserves the same
bound with \(h\) replaced by \(a+d\).  Hadamard's inequality on the
Sylvester matrix proves (6.7).  A physical common root gives a null vector
modulo \(p\); the number of distinct common \(\mathbf F_p\)-roots is at most
the Sylvester nullity, which is at most the \(p\)-adic valuation of its
nonzero determinant.  This proves (6.8).  Summing (6.7) over
\(a+d\le D\) and dividing by \(\log p\) proves (6.9). \(\square\)

Raw resultants require an additional warning.  At \(p=17\), the pole formula
and \(17\mid b_3\) make \(X=-4\equiv13\) a raw root of every \(N_h\) for
\(h\ge4\), but it lies outside all corresponding physical base intervals.
Also \(17\mid\operatorname{Res}(N_2,N_4)\) although the two polynomials have
no common \(\mathbf F_{17}\)-root.  Resultant divisibility is therefore only
an overcount, not an equivalence with physical slice intersection.

One exact conditional face survives.  Let \(\Gamma\) have the heavy lags
\(H=\{h:R_h>T_0\}\) as vertices, joining \(a<a+d\) when
\(p\mid\mathscr R_{a,d}\).  If

\[
 [\mathrm{HEAVY\!\!-RES\!\!-CHROM}]\qquad
 \chi(\Gamma)=O(1),                                          \tag{6.10}
\]

then W1 follows: each color class consists of pairwise disjoint \(Z_h\), so
its total mass is at most \(N\); the light mass is at most
\(DT_0=o(N)\).  Hypothesis (6.10) is open and is not logically comparable
with `[ZERO-TAIL-2]`.  Resultant height does not prove it.

## 7. Direction B: restart is not propagation

### Theorem 7.1 (two-sided restart identities)

Assume \(N_h(r)\equiv0\pmod p\) at a physical base.  Then, in
\(\mathbf F_p\),

\[
 N_{h-1}(r)N_{h+1}(r)\not\equiv0\pmod p                      \tag{7.1}
\]

and, for every \(g\ge1\),

\[
 \boxed{N_{h+g}(r)\equiv N_{h+1}(r)N_g(r+h)\pmod p.}         \tag{7.2}
\]

For \(1\le g<h\), put

\[
 D_{h,g}(X)=\prod_{j=h-g+1}^{h-1}(X+j)^6.
\]

Then the polynomial identity

\[
 \boxed{
 D_{h,g}(X)N_{h-g}(X)
 =N_{h-1}(X)N_g(X+h-g)
  -N_h(X)N_{g-1}(X+h-g)}                                    \tag{7.3}
\]

holds.  At the root \(r\), reduction modulo \(p\) gives the equivalences

\[
 \begin{aligned}
 N_{h+g}(r)\equiv0\pmod p
 &\Longleftrightarrow N_g(r+h)\equiv0\pmod p,\\
 N_{h-g}(r)\equiv0\pmod p
 &\Longleftrightarrow N_g(r+h-g)\equiv0\pmod p.               \tag{7.4}
 \end{aligned}
\]

#### Proof

Equation (6.1) at a root gives

\[
 N_{h+1}(r)\equiv-(r+h)^6N_{h-1}(r)\pmod p.
\]

All physical factors \(r+j\) are units.  If either adjacent numerator
vanished modulo \(p\), backward recurrence would force
\(N_1(r)\equiv0\pmod p\), proving (7.1).  Equation (7.2) is (6.2) with
\(a=h,d=g\), reduced modulo \(p\) after the term containing \(N_h(r)\)
vanishes.

For (7.3), fix \(h\) and view both sides as sequences in \(g\).  The bottom
recurrence shows that both satisfy

\[
 H_{g+1}=P(X+h-g)H_g-(X+h-g+1)^6H_{g-1}.
\]

At \(g=1\), both equal \(N_{h-1}\).  At \(g=2\), both equal

\[
 (X+h-1)^6N_{h-2}
 =P(X+h-1)N_{h-1}-N_h.
\]

Induction proves (7.3).  At a physical root, the carrier \(D_{h,g}(r)\)
and \(N_{h-1}(r)\) are nonzero modulo \(p\), giving (7.4). \(\square\)

The theorem is a detector, not a pump: whenever \(r+h+g\le N-1\), a longer
physical collision occurs exactly when a shifted subinterval was already a
physical collision.  In root-set notation,

\[
 |Z_h\cap Z_{h+g}|
 =\#\{r\in Z_h:r+h\in Z_g\},                                 \tag{7.5}
\]

and the right side has no lower bound in terms of \(R_h\).

This fails concretely in the Apéry orbit.  At \(p=461,h=44\),

\[
 Z_{44}=\{42,45,68,133,142,208,274,283,348,371,374\},
\]

but no \(r\in Z_{44}\) has \(r+44\in Z_{44}\), and

\[
 R_{43}=R_{45}=0,\qquad R_{88}=1=\kappa_{88}.                 \tag{7.6}
\]

Thus same-lag continuation, neighboring-lag propagation, and every exact
universal inequality
\(R_{2h}-\kappa_{2h}\ge c(R_h-\kappa_h)\) with fixed \(c>0\), no error
term, and no exceptional pair \((p,h)\), are false.  This finite example by
itself makes no claim about a suitably qualified eventual asymptotic version.
Even a hypothetical bound \(R_{2h}\ge\eta R_h\), \(\eta\le1\), would not
bootstrap: the degree allowance grows from \(3h\) to \(6h\), while the
propagated mass does not grow.  A successful pump would require a genuinely
superlinear cross-lag aggregation statement, which is another formulation of
the missing quenched arithmetic.

## 8. Two auxiliary distributional faces

These are independent elementary consequences of layer cake.  They are
included because they give additional attack coordinates, not because they
supersede the structure-aware first-line residual.

### Proposition 8.1 ([ZERO-EXCESS-1/2])

Let

\[
 E_D=\sum_{h\le D}(R_h-T_0)_+
    =\sum_{t>T_0}\#\{h:R_h\ge t\}.                            \tag{8.1}
\]

Then

\[
 [\mathrm{ZERO\!\!-EXCESS\!\!-1/2}]\qquad
 E_D\ll\frac{N}{\sqrt{T_0}}                                 \tag{8.2}
\]

implies \(S_D=o(N)\), hence W1.  Moreover `[ZERO-TAIL-2]` implies (8.2),
and the implication is strict.

#### Proof

Directly,

\[
 S_D\le DT_0+E_D.                                            \tag{8.3}
\]

At the hatch scale the two terms on the right are respectively
\(N^{3/4+o(1)}\) and \(N^{7/8+o(1)}\), proving the bridge.  If
\(\mathcal Z_D\ll N\), then

\[
 E_D\le\mathcal Z_D\sum_{t>T_0}t^{-2}
 \le\frac{\mathcal Z_D}{T_0}
 \ll\frac N{\sqrt{T_0}}.                                    \tag{8.4}
\]

The one-full-lag word of Proposition 5.1 has \(E_D=O(D)\), so it satisfies
(8.2) but fails `[ZERO-TAIL-2]`.  Thus the implication is strict. \(\square\)

The fixed power saving in (8.2) prevents this from being merely W1 with the
low tail removed.  For example, \(\asymp N/D\) upper-half lags of height
\(\asymp D\) have total mass \(\asymp N\) but excess \(\asymp N\), and so
fail (8.2).

### Proposition 8.2 ([LAG-LORENTZ-3/2])

Let

\[
 A_D(t)=\#\{h\le D:R_h\ge t\},
 \qquad
 \Lambda_D=\max_{t\ge1}t^3A_D(t)^2.                          \tag{8.5}
\]

Then

\[
 S_D^3\le27D\Lambda_D.                                      \tag{8.6}
\]

Consequently the residual

\[
 [\mathrm{LAG\!\!-LORENTZ\!\!-3/2}]\qquad
 D\Lambda_D\ll N^3                                         \tag{8.7}
\]

implies W1.  `[ZERO-TAIL-2]` implies (8.7), strictly.

#### Proof

Reorder the lag counts as
\(R_{(1)}\ge\cdots\ge R_{(D)}\).  Since
\(A_D(R_{(j)})\ge j\),

\[
 R_{(j)}\le(\sqrt{\Lambda_D}/j)^{2/3}.
\]

Therefore

\[
 S_D\le\Lambda_D^{1/3}\sum_{j=1}^D j^{-2/3}
 \le3D^{1/3}\Lambda_D^{1/3},
\]

which is (8.6).  To derive (8.7) from `[ZERO-TAIL-2]`, use
\(A_D(t)\le D\) for \(t\le T_0\), and
\(A_D(t)\le\mathcal Z_D/t^2\) for \(t>T_0\).  Since
\(T_0\le2\sqrt{N/D}\) and \(D\le N\), both ranges give
\(D\Lambda_D=O(N^3)\).  The word in Proposition 5.1 has
\(D\Lambda_D=O(D^4)=o(N^3)\), while it fails `[ZERO-TAIL-2]`, proving
strictness. \(\square\)

Equivalently, if

\[
 \sup_t t^{3/2}A_D(t)\ll\frac{N^{3/2}}{\sqrt D},              \tag{8.8}
\]

then W1 holds.  An independent pointwise sufficient condition is
\(\max_hR_h\ll N/D\): it implies (8.8) because
\(t^{3/2}A_D(t)\le(\max_hR_h)^{3/2}D\).  This is a critical pointwise
relaxation of the proposed `[PT-ANTICONC]` face, while the full condition
(8.8) is its nonuniform extension and permits sparse larger spikes.

## 9. Exact machine verification

The deliverable verifier is `CRON_breakwall2_verify.py`.  It performs the
following exact gates.

1. It rebuilds the continuants and gap numerators over \(\mathbf Z[X]\),
   checking 20 continuant renewal identities, 33 physical-gauge renewal
   identities, 28 backward restart identities, reflection/center/pole laws,
   33 shifted gcds, 36 same-base gcds, and 15 nonzero resultants.
2. It verifies the \(p=17\) cut-pole warning for \(4\le h\le8\) and the
   \(N_2,N_4\) resultant warning.  The all-height cut-pole statement follows
   from the proved pole formula, not from finite enumeration.
3. It exhaustively checks (8.6) on 19,530 small lag profiles.
4. It materializes the full reflection word of Proposition 5.1, checking
   degree bounds, \(Q_D=0\), energy, stripped variance, both strictness
   directions, and `[MIRROR-WL1]` with exact rational harmonic sums.
5. It checks the \(p=461,h=44\) restart counterexample.
6. It extracts only `orbit(p)` by AST from `CRON_b1_crosscorr.py`.  At
   \(p\in\{1009,3001,10007\}\), with
   \(L=\lceil\log\log p\rceil\), it computes every physical \(R_h\) for
   \(h\le D\), checks (3.1) root by root, verifies the displayed new
   residuals and bridge inequalities exactly, and checks the abstract
   strictness profiles.  Decimal displays are not used by any assertion.

Reproduction:

```sh
PYTHONPYCACHEPREFIX=/tmp/codex-breakwall2-pycache \
  python3 -m py_compile CRON_breakwall2_verify.py
PYTHONPYCACHEPREFIX=/tmp/codex-breakwall2-pycache \
  python3 -u CRON_breakwall2_verify.py
```

SHA-256 of the verifier:

```text
71e3dd7c02e7c931d3d83265749fb509bf41d308a6560a5ca1162589a5684726  CRON_breakwall2_verify.py
```

Output:

```text
SYMBOLIC_GATES PASS renewal=20 reflection=9 gap_renewal=33 backward=28 shifted_gcd=33 cross_gcd=36 cross_resultant=15 cut_poles=5
LORENTZ_FINITE_GATE PASS profiles=19530
MIRROR_WORD_GATE PASS N=200000 D=2000 d=1999 R_d=5994 Q_D=0 energy=423974 stripped2=35928036 Z2/N=35928036/200000 Wmu*H/N=0.084727 W1model_S=161000 W1model_Wmu*H/N=2.877057
RESTART_COUNTEREXAMPLE PASS p=461 h=44 R_h=11 continuation=0 neighbors=0,0 R_2h=1
FINITE_FIELD_GATES
p N L D T0 S maxR W2 Z2 E E^2*T0/N^2 W32sq D*W32sq/N^3 B sumMu2 Wmu WmuH/N maxR*D/N spike_old/excess/lorentz W1model_S/fails_excess/fails_lorentz
1009 1008 2 64 4 90 5 135 50 2 16/1016064 6075 388800/1024192512 0 37 0 0.000000000 320/1008 4225/True/True 992/True/True partners=90
3001 3000 3 165 5 250 8 360 98 5 125/9000000 43200 7128000/27000000000 0 126 0 0.000000000 1320/3000 27225/True/True 2952/True/True partners=250
10007 10006 3 301 6 475 8 747 147 4 96/100120036 186003 55986903/1001801080216 1 244 0 0.000000000 2408/10006 90601/True/True 9900/True/True partners=474
FINITE_FIELD_GATES PASS
ALL_BREAKWALL2_GATES_PASS
```

## 10. Terminal status

- `[MIRROR-WL1] => W1`: **PROVED**.
- `[ZERO-TAIL-2] => [MIRROR-WL1]`, with strict reverse failure:
  **PROVED**.
- Exact physical mirror skeleton including the endpoint: **PROVED**.
- Reflection/window/energy/\(Q_D\) deduction route to stripped BDH:
  **DEAD by explicit word countermodel**.
- All-pair characteristic-zero coprimality and shifted slice identity:
  **PROVED**.
- Fixed-prime resultant-height route: **DEAD at the required scale**.
- `[HEAVY-RES-CHROM] => W1`: **PROVED implication; OPEN hypothesis**.
- Renewal-only family propagation: **FALSE**, with an actual Apéry
  counterexample.
- `[ZERO-EXCESS-1/2]` and `[LAG-LORENTZ-3/2]`: **PROVED auxiliary
  residuals**.
- W1--W7 for the actual Apéry orbit: **OPEN**.
