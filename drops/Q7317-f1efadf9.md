ANSWER Q7317 f1efadf9

# Fixed-exponent top-shell Mellin route: exact lemma and exact obstruction

Repository audit pinned to public `main` commit
`1f837bc38da5fdcbba77063708d77f69ea7864e3`.  The mathematical inputs used
below are the current files

- `problems/3.2/oracleC_result.tex`,
- `problems/3.2/toric_mellin_square.tex`,
- `problems/3.2/toric_fiber_k3.tex`,
- `problems/3.2/hasse_franel_descent.tex`,
- `problems/3.2/research/scripts/q32_marked_scalar_character_sum.py`, and
- `problems/3.2/research/scripts/q32_cartier_packet_audit.py`.

The delivery branch is used only for this answer, not as mathematical source.

## Verdict

There is a strong positive exact reformulation, followed by a sharp negative
answer to the proposed Chebotarev step.

**Proved positive statement.**  Fix `n`, put `N=n-1`, and let

\[
  \frac n2<p\le n,\qquad j=n-p.
\]

Apart from the two non-hit endpoints and the finite prime `p=5`, every
 top-shell hit is exactly

\[
 p\mid b_n
 \quad\Longleftrightarrow\quad
 -\sum_{(x,y,z)\in(\mathbf F_p^\times)^3}
       \Lambda(x,y,z)^N=0\quad\hbox{in }\mathbf F_p,
\]

where the exponent `N=n-1` is independent of `p`.  Equivalently, it is the
vanishing of the `N`-th multiplicative moment of the nonzero-fiber count
function of the fixed Laurent polynomial `Lambda`.

There is also one explicit cyclotomic Frobenius trace at each `p`.  If
`omega_p` is a Teichmuller character and

\[
 \chi_{p,n}=\omega_p^N,
\]

then the algebraic integer

\[
 \widetilde M_{p,n}
 =-\sum_{u\in U(\mathbf F_p)}\chi_{p,n}(\Lambda(u))
\]

is a single compactly supported Frobenius supertrace on the fixed threefold
`U`, twisted by a rank-one Kummer sheaf, and its reduction at a chosen prime
above `p` is `b_j mod p`.

**Proved obstruction.**  This is not a single compatible system, fixed
motive, or fixed Frobenius trace as `p` varies.  The character has exact order

\[
 q_{p,n}=\frac{p-1}{\gcd(p-1,N)},
\]

so both its Kummer eigenspace and its minimal coefficient field
`Q(mu_{q_{p,n}})` move with `p`.  More fundamentally, torus orthogonality
shows that the fixed polynomial `Lambda^N` is being tested by a
**p-dependent coefficient projector**: for every interior top-shell prime,

\[
 b_j\equiv C_N(p-1)\pmod p,
\]

where `C_N(d)` is the sum of exactly 22 coefficients of `Lambda^N` at the
moving lattice positions `d mu`, `mu in P cap Z^3`.  Thus fixing the exponent
does not fix the coefficient functional.

**Chebotarev/large-sieve verdict.**  Standard Chebotarev trace-zero theorems
and the large sieve for Frobenius do not apply.  There are three independent
reasons:

1. the Kummer character/eigensummand varies with the very prime whose
   Frobenius is being evaluated;
2. the event is reduction of a cyclotomic trace modulo a prime above the
   defining characteristic `p`, not vanishing of the characteristic-zero
   trace;
3. the compatible-system large sieve reduces a fixed trace at many
   auxiliary coefficient primes, whereas here there is only one moving,
   self-characteristic congruence.

The second point is not merely formal.  The repository's explicit top-shell
example

\[
 (n,p,j)=(16,11,5)
\]

has

\[
 \widetilde M_{11,16}=-33\ne0,
 \qquad
 \widetilde M_{11,16}\equiv b_5\equiv0\pmod{11}.
\]

Therefore even a hypothetical trace-zero Chebotarev theorem would miss an
actual hit.

The route gives an exact arithmetic lemma and a clean geometric packet of
uniform local complexity, but it does **not** give

\[
 \sum_{\substack{n/2<p\le n\\p\mid b_n}}\log p=o(n).
\]

A new theorem on defining-characteristic divisibility of a moving family of
Kummer Mellin traces, or a new identity replacing the moving Cartier shell
by one fixed integral trace, would still be required.

## 1. Top-shell reduction, including the boundary exceptions

Let

\[
 b_m=\sum_{r=0}^m\binom mr^2\binom{m+r}{r}^2.
\]

For `p in (n/2,n]`, write `n=p+j`; then `0<=j<p`.  The `p`-Lucas
congruence gives

\[
 b_n=b_{p+j}\equiv b_1b_j=5b_j\pmod p.
\]

Consequently

\[
 p\mid b_n\Longleftrightarrow p\mid b_j
\]

for every `p>=7`.  The asserted equivalence is false at `p=5`, because
`b_1=5`; this is a finite exception and contributes at most `log 5` to any
asymptotic top-shell radical.

There are two endpoint cases, neither of which is a hit for `p>=7`.

- If `j=0`, then `p=n` and `b_n congruent 5 mod p`.
- If `j=p-1`, then `n=2p-1`.  Apery reflection gives
  `b_{p-1} congruent b_0=1 mod p`, so again `b_n congruent 5 mod p`.

Thus every genuine top-shell hit with `p>=7` lies in

\[
 1\le j\le p-2.
\]

On this interior range

\[
 N=n-1=(p-1)+j,
 \qquad N\equiv j\pmod{p-1},
\]

and both exponents `N` and `j` are positive.  Positivity matters at the zero
locus of `Lambda`; it is exactly why the endpoint `j=0` cannot be handled by
blindly replacing exponent zero by exponent `p-1`.

## 2. The exact fixed-exponent moment lemma

The fixed Laurent polynomial from `oracleC_result.tex` is

\[
 \Lambda(x,y,z)=
 \frac{(1+x)(1+y)(1+z)
       ((1+y)(1+z)+xyz)}{xyz}.
\]

It has the constant-term realization

\[
 \operatorname{CT}_{x,y,z}\Lambda^m=b_m
 \qquad(m\ge0).
\]

Put

\[
 T_p=(\mathbf F_p^\times)^3,
 \qquad
 \mu_p(a)=\#\{u\in T_p:\Lambda(u)=a\}
 \quad(a\in\mathbf F_p^\times).
\]

### Lemma 2.1 (fixed-exponent top-shell moment)

Let `p>=7`, `n/2<p<=n`, `j=n-p`, and suppose `1<=j<=p-2`.  Put
`N=n-1`.  Then in `F_p`,

\[
 \boxed{
 b_j
 =-\sum_{u\in T_p}\Lambda(u)^N
 =-\sum_{a\in\mathbf F_p^\times}\mu_p(a)a^N.
 }
\]

Consequently

\[
 \boxed{
 p\mid b_n
 \Longleftrightarrow
 -\sum_{u\in T_p}\Lambda(u)^N=0\quad\hbox{in }\mathbf F_p.
 }
\]

#### Proof

The marked-coordinate theorem in `oracleC_result.tex` gives, for
`1<=j<=p-2`,

\[
 b_j=-\sum_{u\in T_p}\Lambda(u)^j
 \quad\hbox{in }\mathbf F_p.
\]

Since `N-j=p-1`, Fermat gives `a^N=a^j` for every nonzero
`a in F_p`.  At a zero of `Lambda`, both powers are zero because `j,N>0`.
Hence `Lambda(u)^N=Lambda(u)^j` for every `u in T_p`, proving the first
identity.  Grouping points by the nonzero value of `Lambda` proves the
fiber-moment identity.  The last equivalence follows from
`b_n congruent 5b_j mod p` and `p not equal 5`.  QED

This is the exact positive answer to the first part of the question: the
exponent really is the fixed integer `n-1`, not a residue label chosen after
`p` is known.

### 2.2. The K3 trace form

`toric_fiber_k3.tex` proves the exact deflated fiber formula

\[
 \mu_p(a)=p^2-6p+12+\theta_p(a),
\]

where `theta_p(a)` is a weight-two Frobenius trace of dimension at most 3 and

\[
 |\theta_p(a)|\le3p.
\]

The character `a mapsto a^N` is nontrivial on `F_p^*` on the interior range,
so its sum over `F_p^*` is zero.  Lemma 2.1 therefore also gives

\[
 \boxed{
 b_j=-\sum_{a\in\mathbf F_p^\times}\theta_p(a)a^N
 \quad\hbox{in }\mathbf F_p.
 }
\]

Thus the scalar is the multiplicative Mellin transform of the fixed K3
fiber-trace function.  This is a genuine geometric compression of each
fiber, but it does not remove the moving multiplicative character.

### 2.3. What the Mellin-square formula does and does not fix

Write

\[
 A_r(t)=\sum_{k=0}^r\binom rk^2t^k.
\]

The exact identity in `toric_mellin_square.tex` is

\[
 b_j=-\sum_{t\in\mathbf F_p^\times}
          t^jA_j(t)A_{p-1-j}(t).
\]

Since `t^j=t^N`, this becomes

\[
 b_j=-\sum_{t\in\mathbf F_p^\times}
          t^N A_j(t)A_{p-1-j}(t).
\]

Only the scalar monomial `t^N` has become fixed.  The two Legendre
polynomials still have indices `j=N-(p-1)` and `p-1-j`, hence still vary with
`p`.  When `j<=(p-1)/2`, the weighted-square form

\[
 b_j=-\sum_{t\ne0}t^j(1-t)^{p-1-2j}A_j(t)^2
\]

has the same issue: both `A_j` and the weight exponent depend on `p`.
The full three-variable `Lambda^N` formulation is the only one of these
forms in which the complete algebraic integrand is fixed; Section 5 below
shows that the price is a moving coefficient extractor.

## 3. A single Frobenius trace at each prime

Let `omega_p:F_p^* -> mu_{p-1}` be a Teichmuller character.  Define

\[
 g_{p,n}=\gcd(N,p-1),
 \qquad
 q_{p,n}=\frac{p-1}{g_{p,n}},
 \qquad
 s_{p,n}=\frac{N}{g_{p,n}}.
\]

Then

\[
 \gcd(q_{p,n},s_{p,n})=1,
 \qquad
 \chi_{p,n}:=\omega_p^N=(\omega_p^{g_{p,n}})^{s_{p,n}},
\]

and `chi_{p,n}` has exact order

\[
 \boxed{q_{p,n}=\frac{p-1}{\gcd(p-1,N)}}.
\]

On the interior range `q_{p,n}>1`.  Put

\[
 U=(\mathbf G_m)^3\setminus V(\Lambda).
\]

Let `K_{chi_{p,n}}` be the rank-one Kummer sheaf with trace function
`chi_{p,n}`, and put

\[
 \mathcal L_{p,n}=\Lambda^*\mathcal K_{\chi_{p,n}}.
\]

The repository's trace formula gives the cyclotomic algebraic integer

\[
 \boxed{
 \widetilde M_{p,n}
 =-\sum_{u\in U(\mathbf F_p)}
       \chi_{p,n}(\Lambda(u))
 =-\operatorname{Tr}_{\rm alt}
   \left(\operatorname{Frob}_p\mid
         R\Gamma_c(U_{\overline{\mathbf F}_p},\mathcal L_{p,n})
   \right).
 }
\]

Its values lie in the minimal cyclotomic field

\[
 E_{p,n}=\mathbf Q(\mu_{q_{p,n}}).
\]

Because `q_{p,n}|p-1`, the rational prime `p` splits completely in
`E_{p,n}`.  Choose the prime `mathfrak p_{p,n}|p` corresponding to the
chosen Teichmuller embedding.  Reduction gives

\[
 \boxed{
 \widetilde M_{p,n}\pmod{\mathfrak p_{p,n}}
 =-\sum_{u\in T_p}\Lambda(u)^N
 =b_j\quad\hbox{in }\mathbf F_p.
 }
\]

Thus the hit event is exactly

\[
 \boxed{
 p\mid b_n
 \Longleftrightarrow
 \mathfrak p_{p,n}\mid\widetilde M_{p,n}.
 }
\]

This is a single Frobenius supertrace **for each prime**.  It is not the
statement `widetilde M_{p,n}=0`.

### 3.1. The actual top-shell counterexample to trace-zero replacement

At `(n,p,j)=(16,11,5)`, the character is quadratic:

\[
 N=15,\qquad g=5,\qquad q=2.
\]

`oracleC_result.tex` records

\[
 \widetilde M_{11,16}=-33.
\]

Hence

\[
 11\mid b_{16},
 \qquad
 11\mid\widetilde M_{11,16},
 \qquad
 \widetilde M_{11,16}\ne0.
\]

This single example rigorously rules out replacing the desired divisibility
event by the trace-zero event to which a Serre-style Chebotarev theorem
would naturally apply.

## 4. The closest integral Kummer variety, and why it is still a packet

For a fixed positive integer `q`, the natural integral cover is

\[
 Y_q:\quad v^q=\Lambda(x,y,z)
\]

over `U`.  Away from the fixed bad set of the integral compactification and
from primes dividing `q`, this is a finite etale Kummer cover.  After
adjoining `mu_q`, its compactly supported cohomology decomposes into the
`q` Kummer-character eigenspaces, and the trace in Section 3 is the
appropriate eigenspace trace.

Thus, for each **fixed** pair `(q,s mod q)`, one obtains an ordinary integral
geometric compatible system.  The top-shell problem does not remain in one
such system.  The prime itself chooses

\[
 (q,s)=\left(
 \frac{p-1}{\gcd(p-1,N)},
 \frac{N}{\gcd(p-1,N)}\bmod q
 \right).
\]

In particular, the local monodromy order of the selected rank-one Kummer
factor is exactly `q_{p,n}`.  Two factors with different `q` cannot be one
fixed rank-one Kummer local system, because their geometric monodromy groups
have different orders.

There is an exact arithmetic relation behind this movement.  If
`g=gcd(N,p-1)`, then

\[
 g\mid N,\qquad p=1+qg,\qquad s=N/g.
\]

Thus fixing `q` and the actual divisor `s=N/g` already determines `p` by

\[
 p=1+\frac{qN}{s}.
\]

Fixing only `(q,s mod q)` permits at most the divisors `s|N` in that residue
class.  The natural compatible systems therefore do not provide one system
sampled at all the top-shell primes; the system label is part of the prime
data.

### 4.1. A fixed-`n` omnibus packet

One can formally package all Kummer orders for a fixed `n`.  Define

\[
 L_n=\operatorname{lcm}_{\substack{n/2<p<n\\p\ {m prime}\\1\le n-p\le p-2}}
       q_{p,n}.
\]

Every selected character occurs in the geometric character decomposition of
an `L_n`-Kummer packet, or equivalently in the collection of covers `Y_q`
for `q|L_n`.  But the projector onto the desired eigenspace still depends on
`p`.  The ordinary Frobenius trace of the full cover is the sum of all
character traces; it is not `widetilde M_{p,n}`.

The complexity of this tautological packet is also not uniform.  Since every
`q_{p,n}<=N`,

\[
 L_n\mid\operatorname{lcm}(1,2,\ldots,N),
 \qquad
 \log L_n\le\psi(N)=N+o(N).
\]

The full `L_n`-cover has degree `L_n`; using the repository's uniform bound
`B_Lambda` for each rank-one eigensheaf gives total cohomological dimension
at most `B_Lambda L_n`.  The coefficient field needed for all projectors is
contained in `Q(mu_{L_n})`, whose degree is `phi(L_n)`.  This is a packet of
potentially exponential rank, not a bounded-rank motive attached to `n`.

This construction is useful as a precise obstruction: a fixed-`n` integral
variety can contain every desired trace, but only as a **p-dependent
isotypic component**.  Chebotarev does not turn a moving idempotent into one
fixed conjugacy condition.

## 5. The coefficient-extraction obstruction is an exact 22-term shell

The previous section describes the obstruction in sheaf language.  Torus
orthogonality makes it completely explicit over the integers.

Write

\[
 \Lambda(x,y,z)^N
 =\sum_{(u,v,w)\in\mathbf Z^3}
   c_N(u,v,w)x^uy^vz^w,
 \qquad c_N(u,v,w)\in\mathbf Z_{\ge0}.
\]

The Newton polytope of `Lambda` is

\[
 P=\{(u,v,w)\in\mathbf R^3:
 -1\le u,v,w\le1,
 \ u-v\le1,
 \ u-w\le1\}.
\]

It has exactly 22 lattice points.  The audit script gives the exact closed
coefficient formula

\[
 \boxed{
 c_N(u,v,w)
 =\sum_{k=0}^N
   \binom Nk
   \binom N{k-u}
   \binom{2N-k}{N-v}
   \binom{2N-k}{N-w},
 }
\]

with the convention that an out-of-range binomial coefficient is zero.
In particular,

\[
 c_N(0,0,0)=b_N.
\]

For a positive integer `d`, define the first Cartier shell

\[
 C_N(d)=\sum_{\mu\in P\cap\mathbf Z^3}c_N(d\mu).
\]

### Proposition 5.1 (moving-shell identity)

Under the hypotheses of Lemma 2.1, with `d=p-1`,

\[
 \boxed{
 b_j\equiv C_N(p-1)\pmod p.
 }
\]

#### Proof

For any integer `e`,

\[
 \sum_{x\in\mathbf F_p^\times}x^e
 =\begin{cases}
   -1,&p-1\mid e,\\
   0,&p-1\nmid e
  \end{cases}
 \quad\hbox{in }\mathbf F_p.
\]

Applying this in three variables gives

\[
 -\sum_{u\in T_p}\Lambda(u)^N
 =\sum_{(r,s,t)\in (p-1)\mathbf Z^3}
    c_N(r,s,t)
 \quad\hbox{in }\mathbf F_p.
\]

Now

\[
 1<\frac{N}{p-1}<2
\]

on the interior top shell.  Since the support of `Lambda^N` lies in `NP`,
a lattice point `(p-1)mu` can occur only when

\[
 \mu_i\in\{-1,0,1\},
 \qquad
 \mu_1-\mu_2\le1,
 \qquad
 \mu_1-\mu_3\le1.
\]

These are exactly the 22 points of `P cap Z^3`.  Therefore the right side is
`C_N(p-1)`.  Lemma 2.1 completes the proof.  QED

This proposition pinpoints the failure of the hoped-for fixed-motive
interpretation.  The polynomial `Lambda^N` is fixed after `n` is fixed, but
the linear functional

\[
 F\longmapsto
 \sum_{\mu\in P\cap\mathbf Z^3}[x^{d\mu_1}y^{d\mu_2}z^{d\mu_3}]F
\]

is evaluated at the moving dilation `d=p-1`.  Fixing the exponent did not
fix the cohomological or coefficient projector.

The endpoint behavior confirms that this is structural rather than a
notation issue.

- At `j=0`, exponent replacement fails at `Lambda=0`; the raw `N=p-1`
  power sum is not `b_0`.
- At `j=p-1`, one has `N=2(p-1)`, so the shell expands from the 22 points of
  `P` to the 100 lattice points of `2P`; the raw endpoint sum is again not
  the marked Apery coefficient.

### 5.2. Height of the individual moving shells

All coefficients of `Lambda` are nonnegative and

\[
 \sum_e [x^e]\Lambda=\Lambda(1,1,1)=40.
\]

Hence

\[
 0<C_N(d)\le40^N.
\]

For each hit prime, Proposition 5.1 supplies the integer divisibility

\[
 p\mid C_N(p-1).
\]

But these are different integers for different primes.  The tautological
common carrier

\[
 \mathfrak C_n=\operatorname{lcm}_{n/2<p<n}C_N(p-1)
\]

satisfies only

\[
 \log\mathfrak C_n
 \le N\log40\,\bigl(\pi(n)-\pi(n/2)\bigr)
 =O\!\left(\frac{n^2}{\log n}\right).
\]

This is far above `o(n)`.  The identity does not prove that the lcm is this
large, but it supplies no cross-`p` divisibility relation that would make it
small.  Any saving would have to be a new theorem about common factors or
linear dependence among the moving shells.

## 6. The same obstruction in Hasse-Franel coordinates

`hasse_franel_descent.tex` proves

\[
 K_p(x)^2=(1+x)^{p-1}A_p(\phi(x))
 =\sum_{m=0}^{p-1}b_m\Psi_{p,m}(x),
\]

where

\[
 \phi(x)=\frac{x(1-8x)}{1+x},
 \qquad
 \Psi_{p,m}(x)=x^m(1-8x)^m(1+x)^{p-1-m}.
\]

The polynomials `Psi_{p,m}` form a triangular, `p`-dependent basis, and

\[
 b_m=[K_p^2]_{\Psi_{p,m}}.
\]

For `m=j=N-(p-1)`, one can write formally

\[
 \Psi_{p,j}(x)
 =(1+x)^{p-1}\phi(x)^j
 =\phi(x)^N
  \left(\frac{(1+x)^2}{x(1-8x)}\right)^{p-1}.
\]

On an `F_p`-point where the displayed rational factor is defined and
nonzero, its `(p-1)`-st power is 1.  This does **not** turn the coordinate
`[K_p^2]_{Psi_{p,j}}` into a coefficient against the fixed rational function
`phi^N`:

1. the displayed equality is an equality in `F_p(x)` with a nonconstant
   `(p-1)`-st power factor, not an algebraic identity setting that factor to
   1;
2. equality of functions on an open subset of `F_p` does not preserve a
   coordinate in a polynomial basis;
3. the exceptional points `0`, `1/8`, and `-1` are precisely where zeros,
   poles, and endpoint corrections live;
4. both the truncation `K_p` and the basis `Psi_{p,m}` change with `p`.

Thus the Hasse-Franel formula reaches the same endpoint as the toric shell:
`b_j` is a defining-characteristic coordinate, not the value of one fixed
characteristic-zero period or Frobenius trace.

## 7. Exact conductor and coefficient-field ledger

There are two different notions of complexity here.  Conflating them would
hide the obstruction.

### 7.1. Geometric complexity at one prime: uniformly bounded

The compactification used in `oracleC_result.tex` has the principal divisor
of `Lambda` consisting of

- four simple zero components, and
- six simple pole components.

For every interior top-shell prime the character is nontrivial.  Its local
monodromy along each of these ten original divisors is therefore nontrivial
and tame; every Swan conductor is zero.  One fixed log resolution of this
fixed divisor works in all sufficiently large characteristics.  The
repository proves an absolute constant `B_Lambda` such that

\[
 \sum_i\dim H_c^i
 (U_{\overline{\mathbf F}_p},\mathcal L_{p,n})
 \le B_\Lambda
\]

for all `(p,n)` in the range.  Hence the rank-one geometric conductor is
`O(1)`, independent of both `p` and `n`.

This is a real gain: there is no hidden degree-`n` Laurent hypersurface in
the local sheaf.

### 7.2. Arithmetic character complexity: exactly `q_{p,n}`

The part that grows is

\[
 q_{p,n}=\frac{p-1}{\gcd(p-1,n-1)},
 \qquad
 2\le q_{p,n}\le p-1<n.
\]

The minimal coefficient field has degree

\[
 [E_{p,n}:\mathbf Q]=\varphi(q_{p,n}).
\]

Its discriminant is exactly

\[
 \boxed{
 |D_{E_{p,n}}|
 =q_{p,n}^{\varphi(q_{p,n})}
  \prod_{\ell\mid q_{p,n}}
  \ell^{-\varphi(q_{p,n})/(\ell-1)}.
 }
\]

Consequently

\[
 \log|D_{E_{p,n}}|
 \le\varphi(q_{p,n})\log q_{p,n}
 \le(n-1)\log(n-1).
\]

The full Kummer cover `Y_q` has degree exactly `q`.  It has additional bad
arithmetic reduction only at the fixed bad set for `Lambda` and at primes
dividing `q`; the latter contribution has logarithmic radical at most
`log q`.

### 7.3. Size and norm of the cyclotomic trace

The K3 deflation gives, for every complex embedding,

\[
 |\widetilde M_{p,n}|
 =\left|\sum_{a\ne0}\chi_{p,n}(a)\theta_p(a)\right|
 \le3p(p-1)<3p^2.
\]

Therefore, if the trace is nonzero,

\[
 \left|N_{E_{p,n}/\mathbf Q}(\widetilde M_{p,n})\right|
 \le(3p^2)^{\varphi(q_{p,n})}.
\]

A hit implies that this norm is divisible by `p`, but the norm itself is
again `p`-dependent, and its logarithmic height can be

\[
 O(\varphi(q_{p,n})\log p)=O(n\log n).
\]

This norm construction therefore does not produce a small common carrier.

## 8. Why Chebotarev does not yield the desired radical bound

A Chebotarev route would need a fixed global representation (for fixed `n`)
and a fixed conjugacy-stable condition on its Frobenius classes.  The natural
construction fails this in several precise ways.

### 8.1. The representation/eigenspace moves

The natural compatible system for a fixed Kummer character has a fixed order
`q` and a fixed character of `mu_q`.  Our prime `p` selects the order
`q_{p,n}` and the character exponent `s_{p,n}`.  A full Kummer packet contains
all of them only at the cost of a moving projector and rank growing with
`L_n`.  The desired scalar is not the trace of the full packet.

### 8.2. Divisibility is not trace zero

The event is

\[
 \mathfrak p_{p,n}\mid\widetilde M_{p,n},
\]

not

\[
 \widetilde M_{p,n}=0.
\]

The explicit value `-33` at `(n,p)=(16,11)` proves that these events differ
inside the actual top shell.  Weil-type size bounds do not repair this:
even the repository's rigorous `O(p^2)` bound leaves many nonzero multiples
of `p`.

### 8.3. The modulus is the Frobenius characteristic itself

The trace formula can be computed with any auxiliary `ell != p`, as in the
repository's choice `ell=2`.  The hit condition, however, reduces the
resulting algebraic integer at a prime above **the same `p` as the base
characteristic**.  This reduction is not one of the auxiliary `ell`-adic
realizations used by ordinary Chebotarev.  It is a Hasse-Witt/Cartier, or
self-characteristic, condition.

Even for a genuinely fixed motive, self-characteristic divisibility of a
weight greater than one trace is not a formal Chebotarev condition.  It
becomes trace zero only in special low-weight situations where an archimedean
bound is strictly smaller than `p`; that inequality is false here and the
`-33` example shows the failure concretely.

### 8.4. Uniformity is at the wrong scale

For each fixed `n` there are only finitely many target primes, all of size
comparable to `n`, while the natural packet and its coefficient field also
depend on `n`.  Effective Chebotarev constants depend on representation
rank, ramified primes, and field discriminant.  The ledger above gives no
uniform regime in which a Chebotarev asymptotic for a fixed object has begun
by the time the prime variable reaches `p asymp n`.

Therefore Serre's trace-zero Chebotarev method does not imply the required
`o(n)` logarithmic radical from this construction.

## 9. Why the large sieve for Frobenius does not apply

The large sieve for Frobenius in compatible systems, in the form developed
by Kowalski and used for cyclotomic trace functions by Perret-Gentil, starts
with

1. a fixed base characteristic,
2. a fixed number field of coefficients,
3. a family over auxiliary coefficient primes `lambda`,
4. characteristic polynomials independent of `lambda`, and
5. uniformly bounded conductor and controlled integral monodromy.

The top-shell diagonal has a different quantifier pattern:

- the base characteristic `p` varies;
- the coefficient field `Q(mu_{q_{p,n}})` varies;
- the Kummer character and eigensummand vary;
- the exceptional condition is reduction at `mathfrak p_{p,n}|p`, not at a
  collection of auxiliary primes `lambda` distinct from `p`.

A Frobenius large sieve can detect that one fixed algebraic trace belongs to
an exceptional set by imposing compatible restrictions modulo many
auxiliary primes.  Here the information is only one congruence

\[
 \widetilde M_{p,n}\equiv0\pmod{\mathfrak p_{p,n}}
\]

at one moving prime.  It gives no restrictions modulo the other auxiliary
coefficient primes on which the sieve relies.

A classical analytic large sieve for multiplicative characters could bound
mean squares of the sums `widetilde M_{p,n}`.  A mean-square bound does not
bound the number of exact self-characteristic zeros without an additional
integral anti-concentration statement.  Such a statement is precisely the
missing arithmetic theorem, not a consequence of the existing toric or K3
pointwise bounds.

## 10. Exact final statement

The focused route yields the following rigorous lemma and obstruction.

### Theorem 10.1 (exact status of the fixed-exponent route)

For every integer `n` and every prime `p>=7` with `n/2<p<=n`, put
`j=n-p` and `N=n-1`.

1. If `j=0` or `j=p-1`, then `p` is not a divisor of `b_n`.
2. If `1<=j<=p-2`, then
   
   \[
   p\mid b_n
   \Longleftrightarrow
   -\sum_{u\in(\mathbf F_p^\times)^3}\Lambda(u)^N=0
   \Longleftrightarrow
   C_N(p-1)=0\pmod p.
   \]
3. The middle scalar is the reduction of one cyclotomic Frobenius
   supertrace `widetilde M_{p,n}` on the fixed variety `U`, twisted by a
   rank-one Kummer sheaf of exact order
   
   \[
   q_{p,n}=(p-1)/\gcd(p-1,N).
   \]
4. The geometric boundary and Betti complexity of this rank-one twist are
   bounded absolutely, but its Kummer order, coefficient field, and selected
   eigenspace move with `p`.
5. The hit condition is divisibility of `widetilde M_{p,n}` by a prime above
   `p`, not exact vanishing; `(n,p)=(16,11)` is an explicit counterexample to
   equivalence with trace zero.
6. The coefficient-shell identity shows that the defining-characteristic
   projector itself moves with `p`.  Hence the current toric Mellin-square,
   K3-fiber, and Hasse-Franel formulas do not furnish one fixed integral
   compatible system whose fixed Frobenius trace vanishes exactly at the hit
   primes.

Accordingly, no cited Chebotarev or large-sieve theorem yields

\[
 \log\prod_{\substack{n/2<p\le n\\p\mid b_n}}p=o(n)
\]

from these identities.

This is an obstruction to the **natural constructions supplied by the
current formulas**, not a theorem excluding every unrelated future motive
that might accidentally encode the same congruences.  To make this route
work, one needs at least one genuinely new ingredient of one of the
following two forms:

- a uniform theorem for the diagonal self-characteristic divisibility
  
  \[
  \mathfrak p_{p,n}\mid\widetilde M_{p,n}
  \]
  
  across Kummer orders `q_{p,n}<=n`, strong enough to give weighted density
  `o(1)`; or
- an Apery-specific identity replacing the moving 22-coefficient shell
  `C_N(p-1)` by a fixed integral trace or by a common integer carrier of
  logarithmic height `o(n)`.

Neither ingredient is present in the current repository formulas, and
ordinary monodromy language does not supply it.

## References used for the theorem audit

- J.-P. Serre, *Quelques applications du theoreme de densite de
  Chebotarev*, Publ. Math. IHES 54 (1981), 123-201.
- E. Kowalski, *The Large Sieve and its Applications*, Cambridge Tracts in
  Mathematics 175, 2008, Chapter 8.
- E. Kowalski, *The large sieve, monodromy and zeta functions of curves*,
  J. Reine Angew. Math. 601 (2006), 29-69; arXiv:math/0503714.
- C. Perret-Gentil, *Exponential sums over finite fields and the large
  sieve*, arXiv:1703.06965, especially Definition 2.1 and Theorem 2.9.
