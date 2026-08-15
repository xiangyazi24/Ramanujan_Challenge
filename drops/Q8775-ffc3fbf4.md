ANSWER Q8775 ffc3fbf4

# LGTQ-3 extra-codimension audit: density alone is not enough

## Executive conclusion

The proposed density-only route does **not** prove

\[
\max_{h_1,h_2\le C\log n}T_{q;h_1,h_2}(n)
=o\!\left(\frac{n}{(\log n)^3}\right).
\]

There are two distinct issues.

1. **The sieve-dimension normalization in the question divides by the prime twice.** If a local forbidden set has \(\omega(\ell)=|Z_\ell|=O(1)\) residue classes modulo \(\ell\), then its local density is \(g(\ell)=\omega(\ell)/\ell\). The dimension test sums \(g(\ell)\), equivalently \(\omega(\ell)/\ell\), not \(g(\ell)/\ell\). Thus a positive mean number of zero classes gives a positive sieve dimension, not dimension zero.

2. **More importantly, LGTQ-3 contains a prime-indexed hit condition, not a fixed system of congruence exclusions.** A Selberg sieve can add the \(Z_\ell\)-classes to its dimension only when every candidate must avoid a prescribed set of classes modulo every auxiliary prime \(\ell\). In LGTQ the relevant modulus is itself the moving candidate prime, and one asks whether a target residue belongs to \(Z_p\). The mere fact that \(|Z_p|=O(1)\) gives no control over the alignment of that target with \(Z_p\).

In fact there is an exact no-go theorem: for any prescribed moving target \(\rho_p\in\mathbf F_p\), the singleton sets

\[
Z_p=\{\rho_p\}
\]

have size one and no adjacent elements, yet every target hits. Even an exact Poisson marginal law for \(|Z_p|\) does not repair this: one can conditionally place a root at \(\rho_p\) whenever \(|Z_p|\ge1\), preserving the Poisson cardinality law and the no-adjacent constraint, while retaining a positive hit probability \(1-e^{-\lambda}\), rather than the desired \(\asymp1/p\).

Therefore the missing statement is not another bound on \(|Z_p|\). It is a **spatial, centered, moving-target discrepancy estimate** (or a geometric theorem implying one). A natural exact residual is

\[
\tag{LGTQ-MTD}
\max_{q,h_1,h_2}
\left|
\sum_{p\in\mathscr P_{h_1,h_2}(X)}
\left(
\mathbf 1_{Z_p}(\rho_{q;h_1,h_2}(p))
-\frac{|Z_p|}{p}
\right)
\right|
=o\!\left(\frac{X}{(\log X)^3}\right),
\]

with the appropriate LGTQ target map \(\rho_{q;h_1,h_2}\) and prime-triple carrier \(\mathscr P_{h_1,h_2}(X)\). If LGTQ contains several Apéry hit indicators, the required residual is the corresponding mixed centered correlation estimate.

That is the actual extra-codimension problem.

---

## 1. Repository-state check

The current committed Problem 3.2 sources distinguish three facts that should not be conflated.

* `problems/3.2/proof.tex` proves the no-two-consecutive-zeros lemma by recurrence back-propagation and proves the unconditional pointwise bound
  \[
  Z(p)=\#\{0\le r<p:p\mid b_r\}=O(p^{2/3}).
  \]
* The same file records the computation \(Z(p)\le12\) through \(10^6\), mean approximately one, and a Poisson fit for the paired zero count.
* It explicitly says that the corresponding Poisson law for the **actual coefficient zeros** remains conjectural. The theorem ledger also contains an unconditional annealed Poisson theorem for a randomized/restarted Apéry chain; that is not a quenched cross-prime theorem for the deterministic sets \(Z_p\).

Thus \(|Z_p|=O(1)\) is presently a strong working hypothesis/numerical regime, not the unconditional theorem in the checked manuscript. This report grants the stronger hypothesis anyway. The negative conclusion below is therefore stronger: **even uniform boundedness and even a Poisson marginal law do not suffice without spatial information.**

The project’s own `SUBMIT/3.2/README.md` already identifies the same underlying obstruction in a different language: fiber sparsity is not yet a global counting theorem, and a classical large sieve encounters the \(Q^2\) barrier. The analysis below makes that obstruction exact for LGTQ-3.

---

## 2. Correct sieve-dimension normalization

Let \(\mathcal A\) be a sequence of total mass \(X\). For each auxiliary prime \(\ell\), let

\[
\Omega_\ell\subset\mathbf Z/\ell\mathbf Z,
\qquad
\omega(\ell)=|\Omega_\ell|,
\qquad
g(\ell)=\frac{\omega(\ell)}{\ell}.
\]

The sifted set consists of candidates avoiding \(\Omega_\ell\) for every \(\ell<z\). Its local product is

\[
V(z)=\prod_{\ell<z}\left(1-g(\ell)\right)
=
\prod_{\ell<z}\left(1-\frac{\omega(\ell)}{\ell}\right).
\]

For bounded \(\omega(\ell)\),

\[
\log V(z)
=-\sum_{\ell<z}\frac{\omega(\ell)}{\ell}+O(1).
\]

A standard dimension-\(\kappa\) hypothesis is

\[
\sum_{\ell\le z}\frac{\omega(\ell)\log\ell}{\ell}
=\kappa\log z+O(1),
\]

or, after partial summation,

\[
\sum_{\ell\le z}\frac{\omega(\ell)}{\ell}
=\kappa\log\log z+O(1).
\]

Equivalently, in the convention where \(g(\ell)\) already denotes the density \(\omega(\ell)/\ell\), the dimension is read from

\[
\sum_{\ell\le z}g(\ell)\log\ell
=\kappa\log z+O(1),
\qquad
\sum_{\ell\le z}g(\ell)
=\kappa\log\log z+O(1).
\]

The expression in the question,

\[
\sum_{\ell\le z}\frac{g(\ell)}{\ell}
\asymp\sum_{\ell\le z}\frac1{\ell^2},
\]

divides by \(\ell\) a second time. It is not the sieve-dimension sum.

### Consequences

* One forbidden residue class per prime has \(\omega(\ell)=1\), local density \(1/\ell\), and dimension one. This is exactly the local scale of primality of one linear form.
* Three generic prime forms have \(\omega(\ell)=3\) for all but finitely many \(\ell\), hence dimension three.
* A bounded zero set with weighted mean
  \[
  \frac{1}{\log z}
  \sum_{\ell\le z}\frac{|Z_\ell|\log\ell}{\ell}
  \longrightarrow\lambda>0
  \]
  would contribute dimension \(\lambda\), **provided it really appears as a new forbidden residue system**.
* True dimension zero means
  \[
  \sum_\ell\frac{\omega(\ell)}{\ell}<\infty.
  \]
  Then \(V(z)\) tends to a positive constant. Dimension zero therefore does not itself give an \(o(1)\) saving; it gives no logarithmic decay at all.

The comparison with a Kummer half-interval is also misleading in small-sieve terminology. Removing half of all classes modulo every \(\ell\) has local density \(\asymp1\), not \(1/\ell\), and lies in a large-sieve/larger-sieve regime. A constant-size Apéry fiber is sparse relative to a half-interval, but it is on exactly the same \(1/\ell\) local scale as an ordinary prime condition.

---

## 3. What a genuine extra exclusion would give

Let

\[
\mathcal H=\{0,h_1,h_1+h_2\}
\]

and consider the ordinary prime-triple carrier

\[
\mathscr P_{h_1,h_2}(X)
=
\{m\in[X,2X]:m,m+h_1,m+h_1+h_2\text{ are prime}\}.
\]

For an auxiliary prime \(\ell\), the primality sieve removes

\[
P_\ell=-\mathcal H\pmod\ell,
\qquad
\nu_\ell=|P_\ell|.
\]

For all but finitely many \(\ell\), \(\nu_\ell=3\). The corresponding product is

\[
V_{\mathcal H}(z)
=
\prod_{\ell<z}\left(1-\frac{\nu_\ell}{\ell}\right)
\asymp
\frac{\mathfrak S(\mathcal H)}{(\log z)^3},
\]

where the exceptional small primes form the singular series.

Now suppose, counterfactually but instructively, that the Apéry condition produced a **genuine additional forbidden set**

\[
\Xi_\ell\subset\mathbf F_\ell
\]

for every auxiliary prime \(\ell\), and every LGTQ candidate had to avoid \(\Xi_\ell\). Define the number of genuinely new classes

\[
e_\ell=|\Xi_\ell\setminus P_\ell|.
\]

Then the total local product is

\[
V_{\mathrm{tot}}(z)
=
\prod_{\ell<z}
\left(1-\frac{\nu_\ell+e_\ell}{\ell}\right).
\]

If \(e_\ell=O(1)\), then

\[
\log\frac{V_{\mathrm{tot}}(z)}{V_{\mathcal H}(z)}
=
-\sum_{\ell<z}\frac{e_\ell}{\ell}+O(1).
\]

This gives an exact criterion.

### Exact extra-codimension criterion

Let

\[
E(z)=\sum_{\ell<z}\frac{e_\ell}{\ell}.
\]

Then, subject to the usual sieve remainder and level-of-distribution hypotheses,

\[
S_{\mathrm{tot}}(X,z)
\ll
\frac{X\,\mathfrak S(\mathcal H)}{(\log z)^3}
\exp(-E(z)).
\]

Hence:

* If \(E(z)\) is bounded, the Apéry classes change only the constant.
* If \(E(z)\to\infty\), they give an \(o(1)\) factor relative to the dimension-three bound.
* If
  \[
  E(z)=\eta\log\log z+O(1)
  \qquad(\eta>0),
  \]
  then
  \[
  S_{\mathrm{tot}}(X,z)
  \ll
  \frac{X\,\mathfrak S(\mathcal H)}{(\log z)^{3+\eta}}.
  \]

Taking \(z=X^\theta\) with fixed \(\theta>0\) yields

\[
S_{\mathrm{tot}}(X,z)
\ll
\frac{X\,\mathfrak S(\mathcal H)}{(\log X)^{3+\eta}},
\]

which is indeed \(o(X/(\log X)^3)\). A fixed positive \(\eta\) also dominates the mild uniform growth of the singular series when \(h_i\le C\log X\).

This is the precise sense in which a bounded number of **new excluded residue classes** could create extra codimension.

But three qualifications are essential.

1. They must be exclusions for every auxiliary \(\ell\), not membership targets attached only to the moving prime.
2. They must be new after overlap with the ordinary prime-tuple classes; \(|\Xi_\ell|\) alone is not enough.
3. One still needs the sieve axioms for composite squarefree moduli and a level of distribution sufficient to take \(z=X^\theta\). Local cardinalities alone do not control the remainder terms.

LGTQ-3 does not presently supply item 1.

---

## 4. Why the Apéry hit is not an added Selberg-sieve dimension

Write the LGTQ Apéry target abstractly as

\[
\rho_{q;h_1,h_2}(p)\in\mathbf F_p.
\]

The relevant indicator is of the form

\[
H_{q;h_1,h_2}(p)
=
\mathbf 1_{Z_p}\bigl(\rho_{q;h_1,h_2}(p)\bigr),
\]

or a product of several such indicators at the three prime forms. The exact recurrence formula for \(\rho\) does not affect the present sieve audit.

The crucial point is that the modulus \(p\) varies with the candidate. This is a **diagonal moving-modulus incidence**:

\[
(p,\rho(p))\in
\mathcal Z:=\{(p,r):r\in Z_p\}.
\]

It is not a condition saying that one global integer variable avoids \(Z_\ell\) modulo every small auxiliary prime \(\ell\). Therefore there is no multiplicative local-density function \(g(d)\) obtained by Chinese remaindering the \(Z_\ell\), and no extra Selberg dimension follows from \(|Z_p|\).

There is also a hit/avoidance reversal:

* If one declares \(Z_\ell\) to be the sifted-out set, the sieve counts candidates **avoiding** Apéry zeros, the opposite of LGTQ.
* To make the surviving set consist of hits, one would have to remove the complement \(\mathbf F_\ell\setminus Z_\ell\), i.e. almost all residue classes. Gallagher’s larger sieve can exploit this only when the same global set lies in \(Z_\ell\) modulo every auxiliary \(\ell\). LGTQ again does not impose that simultaneous condition.

The heuristic \(|Z_p|/p\) is therefore an expected hit probability only after one proves that \(\rho(p)\) is spatially equidistributed relative to \(Z_p\). That equidistribution is exactly the missing theorem; it cannot be inserted as a sieve density by definition.

### Centered decomposition

For one target indicator, write

\[
T_{q;h_1,h_2}(X)
=
\sum_{p\in\mathscr P_{h_1,h_2}(X)}
\mathbf 1_{Z_p}(\rho(p)).
\]

Then identically

\[
T_{q;h_1,h_2}(X)=M_{q;h_1,h_2}(X)+D_{q;h_1,h_2}(X),
\]

where

\[
M=
\sum_{p\in\mathscr P_{h_1,h_2}(X)}
\frac{|Z_p|}{p},
\]

and

\[
D=
\sum_{p\in\mathscr P_{h_1,h_2}(X)}
\left(
\mathbf 1_{Z_p}(\rho(p))-rac{|Z_p|}{p}
\right).
\]

If \(p\asymp X\), \(|Z_p|\le K\), and the ordinary prime-triple sieve gives

\[
|\mathscr P_{h_1,h_2}(X)|
\ll
\frac{X\,\mathfrak S(\mathcal H)}{(\log X)^3},
\]

then the nominal main term is tiny:

\[
M
\ll
\frac{K\,\mathfrak S(\mathcal H)}{(\log X)^3}.
\]

So the density heuristic would provide far more than the required saving. But all of the difficulty is in \(D\). Cardinality gives only the trivial bound

\[
|D|\le |\mathscr P_{h_1,h_2}(X)|,
\]

which is exactly the unsaved dimension-three scale.

Thus the route succeeds **if and only if** one proves a centered discrepancy estimate for the actual target section. Calling \(|Z_p|/p\) a local density does not prove that estimate.

---

## 5. Exact density-only no-go theorem

### Proposition 5.1: bounded fibers and spacing do not control a moving target

Let \(\mathcal P\) be any finite set of primes and let \(\rho_p\in\mathbf F_p\) be arbitrary. There are sets \(Z_p\subset\mathbf F_p\) satisfying

\[
|Z_p|=1,
\qquad
Z_p\cap(Z_p+1)=\varnothing,
\]

for which

\[
\mathbf 1_{Z_p}(\rho_p)=1
\qquad\text{for every }p\in\mathcal P.
\]

**Proof.** Take \(Z_p=\{\rho_p\}\). A singleton has no two adjacent elements, and it contains the target. \(\square\)

Consequently no theorem whose hypotheses use only

\[
|Z_p|=O(1)
\quad\text{and}\quad
Z_p\cap(Z_p+1)=\varnothing
\]

can imply even \(o(|\mathcal P|)\) moving-target hits, let alone the LGTQ-3 saving.

The same construction applies when the carrier consists of ordinary prime triples: attach to each prime occurring in the triple the singleton containing its LGTQ target. Every ordinary prime triple survives. Hence a density-only theorem is logically compatible with the full Hardy--Littlewood-sized count.

### Proposition 5.2: even a Poisson cardinality law is insufficient

Fix \(\lambda>0\). For each large prime \(p\), choose an integer \(K_p\) with asymptotic Poisson\((\lambda)\) law, truncated only at \(p/3\); the truncation changes the law by \(o(1)\). Construct \(Z_p\) as follows.

* If \(K_p=0\), set \(Z_p=\varnothing\).
* If \(K_p\ge1\), include \(\rho_p\), then place the remaining \(K_p-1\) points with cyclic gaps at least two. This is possible for \(K_p\le p/3\).

Then

\[
|Z_p|=K_p,
\qquad
Z_p\cap(Z_p+1)=\varnothing,
\]

and

\[
\mathbf 1_{Z_p}(\rho_p)=\mathbf 1_{\{K_p\ge1\}}.
\]

Therefore

\[
\mathbb P\bigl(\rho_p\in Z_p\bigr)
\longrightarrow
1-e^{-\lambda},
\]

which is a positive constant, not \(\lambda/p\). For three independent prime coordinates, the corresponding survival probability is \((1-e^{-\lambda})^3\), still a positive constant. Thus the expected count remains a positive constant multiple of the ordinary prime-triple count.

This countermodel proves that a Poisson law for the **number** of roots contains no information about the **location** of the roots relative to the LGTQ target. One needs a joint law or a discrepancy theorem.

---

## 6. What the no-adjacent-zero theorem does and does not give

The recurrence proof is useful and exact. If

\[
b_r\equiv b_{r+1}\equiv0\pmod p,
\qquad 0\le r\le p-2,
\]

then the Apéry recurrence at \(r\) forces \(b_{r-1}\equiv0\), and back-propagation reaches \(b_0\equiv0\), contradicting \(b_0=1\). Thus

\[
Z_p\cap(Z_p+1)=\varnothing.
\]

Its sieve consequences are limited.

### 6.1 A special adjacent-hit configuration is killed

If a local LGTQ condition literally requires

\[
r\in Z_p
\quad\text{and}\quad
r+1\in Z_p
\]

in the same characteristic, its count is zero. This is a genuine exact gain for that one configuration.

### 6.2 It can certify disjointness of two translated exclusion sets

In a hypothetical exclusion sieve containing both \(Z_\ell\) and \(Z_\ell-1\), no adjacency implies

\[
|Z_\ell\cup(Z_\ell-1)|=2|Z_\ell|.
\]

Thus it prevents overlap and can double the extra local class count for exactly that unit translation. If those sets were genuine auxiliary-prime exclusions, this could increase the sieve dimension.

That is not the present LGTQ architecture: the zero condition is a moving hit, and the theorem does not manufacture simultaneous auxiliary-prime exclusions.

### 6.3 For generic shifts it changes only second-order local data

For a random bounded set of size \(K\), the expected number of pairs at a fixed nonzero gap is \(K(K-1)/(p-1)=O(1/p)\). After dividing by \(p\) to obtain a local density, the effect is \(O(1/p^2)\). Since

\[
\sum_p\frac1{p^2}<\infty,
\]

such a bounded-range pair correction normally changes an Euler-product constant, not a logarithmic exponent.

No adjacency supplies no cross-prime decorrelation. It says nothing about whether the one LGTQ target \(\rho_p\) is selected as the isolated root. The singleton counterexample already saturates this obstruction.

Finally, LGTQ asks for a maximum over many \(h_1,h_2\). Even if one unit-gap pattern is impossible, that does not control the nonunit gaps in the maximum.

---

## 7. Large-sieve and Bombieri--Vinogradov audit

### 7.1 The additive large sieve does not turn support sparsity into target cancellation

The classical inequality is

\[
\sum_{q\le Q}\ \sum_{a\bmod q}^{*}
\left|
\sum_{n\le N}c_n e(an/q)
\right|^2
\le
(N+Q^2)\sum_{n\le N}|c_n|^2.
\]

For the centered zero indicator

\[
f_p(x)=\mathbf 1_{Z_p}(x)-\frac{|Z_p|}{p},
\]

Parseval gives

\[
\sum_{a\bmod p}|\widehat f_p(a)|^2
=p|Z_p|-|Z_p|^2
\asymp p|Z_p|.
\]

This is an energy identity, not cancellation at the target. Indeed, if \(Z_p=\{\rho_p\}\), then for every nonzero frequency

\[
|\widehat f_p(a)|=1,
\]

and the phase is perfectly aligned with evaluation at \(\rho_p\). Thus an \(O(1)\)-point support can have maximal coherent Fourier behavior.

This agrees with the project’s earlier Fourier audit in `drops/Q7901-eed6cfc4.md`: cardinality, reflection, spacing, and Parseval do not prevent cross-prime frequency alignment; a genuinely mixed trilinear dispersion estimate is required.

### 7.2 The standard large-sieve exclusion bound recovers the usual dimension, but under the wrong hypothesis

If one global set of integers avoids \(\Omega_\ell\) modulo every \(\ell\), the large sieve gives

\[
S(N,Q)
\le
\frac{N+Q^2}{H(Q)},
\]

where

\[
H(Q)=
\sum_{\substack{d\le Q\\d\ \mathrm{squarefree}}}
\prod_{\ell\mid d}
\frac{\omega(\ell)}{\ell-\omega(\ell)}.
\]

For bounded \(\omega(\ell)\) of mean \(\kappa\), one has heuristically and, under standard regularity, asymptotically

\[
H(Q)\asymp(\log Q)^\kappa.
\]

Taking \(Q\asymp\sqrt N\) yields \(S\ll N/(\log N)^\kappa\). Thus this formulation again confirms that \(O(1)\) forbidden classes correspond to positive finite dimension.

But its hypothesis is simultaneous avoidance modulo all auxiliary primes. It does not apply to a candidate tested only against \(Z_p\) for its own moving modulus \(p\).

### 7.3 The \(Q^2\) barrier is real in the moving-modulus range

In LGTQ the relevant primes are on the same scale as the counting variable. Trying to average directly over moduli \(p\asymp X\) puts \(Q\asymp X\), so the large-sieve constant is \(N+X^2\). The \(X^2\) term dominates at the natural one-dimensional length. The small support of \(Z_p\) does not remove that term.

A Bombieri--Vinogradov theorem controls primes in fixed residue classes on average over moduli only up to roughly the square-root level. It does not directly cover a diagonal condition whose modulus is the candidate prime itself and whose residue set depends on that same prime.

A dispersion or reciprocity transformation could in principle move the problem to shorter auxiliary moduli, but that would be a new Apéry-specific structural theorem, not a consequence of \(|Z_p|=O(1)\).

### 7.4 Gallagher’s larger sieve also has the wrong quantifiers

Gallagher’s theorem is powerful when a **single global set** \(A\subset[1,N]\) occupies at most \(\nu(\ell)\) residue classes modulo each of many primes \(\ell\). In one standard form,

\[
|A|
\le
\frac{\sum_{\ell\in\mathcal L}\log\ell-\log N}
{\sum_{\ell\in\mathcal L}(\log\ell)/\nu(\ell)-\log N},
\]

provided the denominator is positive.

If \(A\bmod\ell\subset Z_\ell\) for every \(\ell\) and \(|Z_\ell|=O(1)\), this can be extremely strong. But LGTQ says only that the candidate associated with modulus \(p\) hits \(Z_p\). It does not put the whole candidate set into \(Z_\ell\) modulo every unrelated \(\ell\). Therefore the larger sieve cannot be invoked from the cardinality statement alone.

---

## 8. What could genuinely provide the extra codimension?

There are several mathematically valid mechanisms, but each contains information absent from \(|Z_p|\).

### 8.1 A fixed-auxiliary-prime exclusion transmutation

The cleanest purely sieve-theoretic route would be an Apéry identity converting the moving hit into a family of fixed congruence exclusions:

\[
\rho_{q;h_1,h_2}(p)\in Z_p
\quad\Longrightarrow\quad
m\notin\Xi_\ell\pmod\ell
\]

for every \(\ell\) in a positive-density auxiliary family, with

\[
\sum_{\ell<z}
\frac{|\Xi_\ell\setminus P_\ell|}{\ell}
\to\infty.
\]

Then Section 3 would apply. No such transmutation is presently in the recurrence/Lucas toolkit.

### 8.2 A centered moving-target discrepancy theorem

The minimal direct input is (LGTQ-MTD):

\[
\max_{q,h_1,h_2}
\left|
\sum_{p\in\mathscr P_{h_1,h_2}(X)}
\left(
\mathbf 1_{Z_p}(\rho(p))-rac{|Z_p|}{p}
\right)
\right|
=o\!\left(\frac{X}{(\log X)^3}\right).
\]

Together with \(|Z_p|=O(1)\), its uncentered main term is already negligible. This theorem is exactly what a successful geometric large sieve, trace-function dispersion argument, or cross-characteristic equidistribution theorem would have to deliver.

For three hit indicators, expand each as mean plus centered part. One must control every nonempty centered mixed term, especially the fully centered threefold correlation. Controlling only one-prime marginals is insufficient.

### 8.3 A true codimension-two incidence in one characteristic

Suppose recurrence algebra turns the LGTQ event into two polynomial equations

\[
F_{q,h_1,h_2}(u)=0,
\qquad
G_{q,h_1,h_2}(u)=0
\]

on a parameter space, and proves that their common locus is a proper codimension-two subvariety uniformly in the parameters. Then the local density could be \(O(1/p^2)\), rather than \(O(1/p)\), producing a genuine extra factor.

The essential word is **independent**. A zero-dimensional fiber of one equation is only codimension one. Moreover a target section can lie entirely inside that divisor. One needs a noncontainment/resultant theorem plus quantitative equidistribution of the target section.

A particularly strong version would show that every simultaneous hit forces

\[
p\mid R(q,h_1,h_2)
\]

for a nonzero integer or polynomial value \(R\) of controlled size. Then the number of possible \(p\) is bounded by the number of prime divisors of \(R\). That would bypass ordinary sieve dimension entirely.

### 8.4 Growing-group monodromy rather than fixed CM congruence classes

Q8772 correctly rules out a fixed CM/Hecke-character condition of positive constant density: it only changes the singular series. A geometric route would need a bad subset whose relative size in the relevant finite monodromy group is itself \(O(1/p)\) (or smaller), together with a large-sieve theorem uniform in the family.

For example, two independent vanishing matrix coefficients in a group of growing size can have relative density \(O(1/p^2)\). A fixed ray-class condition or fixed conjugacy subset has constant density and cannot supply the missing exponent.

### 8.5 Mixed-prime Fourier dispersion

Writing

\[
\mathbf 1_{Z_p}(r)
=\frac{|Z_p|}{p}
+\frac1p\sum_{a\ne0}\widehat{\mathbf 1_{Z_p}}(a)e_p(-ar),
\]

reduces the centered term to sums coupling the Fourier transforms for different primes with the prime-triple carrier. A saving requires cancellation in these mixed sums, not merely Parseval at each prime. This is the same structural residual identified in the repository’s earlier trilinear Fourier audit.

---

## 9. Direct answers to the five questions

### 1. Does \(|Z_p|=O(1)\) automatically give codimension saving?

**No for LGTQ-3.** The proposed dimension-zero calculation is normalized incorrectly. If \(Z_\ell\) were a genuine excluded residue set, \(|Z_\ell|\asymp1\) would contribute positive finite dimension because

\[
\sum_{\ell\le z}\frac{|Z_\ell|}{\ell}
\asymp\log\log z.
\]

But LGTQ asks for membership in \(Z_p\) at a moving modulus, so even that positive-dimension interpretation does not apply. Fiber sparsity alone does not control target alignment.

### 2. What is the bound when the local density is \(O(1/p)\)?

That is the standard finite-dimensional sieve scale. If the total number of genuinely new forbidden classes has weighted mean \(\eta>0\), then a dimension-three prime-tuple sieve becomes dimension \(3+\eta\):

\[
S(X)\ll\frac{X}{(\log X)^{3+\eta}}
\]

up to the singular series and standard sieve-function factors. This is automatically \(o(X/(\log X)^3)\).

If the extra harmonic mass merely satisfies \(E(z)\to\infty\), the saving is

\[
\exp(-E(z))=o(1).
\]

If the actual local density were \(O(1/p^2)\), then its sum would converge and it would usually change only a constant, not give a logarithmic dimension. None of these conclusions applies until the Apéry hit has been converted into a genuine auxiliary-prime exclusion.

### 3. Does minimum gap at least two in \(Z_p\) help?

Only in special aligned configurations. It kills a same-characteristic adjacent pair exactly, and it proves disjointness of \(Z_p\) and \(Z_p-1\) if both occur as exclusion sets. It does not change the first-order one-target density, does not decorrelate different primes, and does not prevent a singleton zero set from following the target. Therefore it gives no general exponent saving for the maximum over \(h_1,h_2\).

### 4. Can the classical large sieve or Bombieri--Vinogradov exploit \(|Z_p|=O(1)\)?

Not by itself. Sparse support bounds Fourier energy but not phase alignment; a singleton already has coherent Fourier coefficients at every frequency. The additive large sieve retains its \(N+Q^2\) constant, and the moving prime moduli lie beyond the Bombieri--Vinogradov square-root range. Gallagher’s larger sieve would be powerful only under the much stronger statement that one global candidate set occupies \(O(1)\) classes modulo every auxiliary prime.

A geometric or trace-function large sieve could work, but only after proving bounded complexity, nontrivial monodromy, noncontainment of the target section, and uniform cancellation. Those are precisely the missing Apéry-specific inputs.

### 5. Is there a direct density-only sieve proof of LGTQ-3?

**No.** The singleton counterexample is decisive, and the Poisson-coupled counterexample shows that even the proposed cardinality law plus no adjacency is insufficient. A direct proof must add at least one of:

* a fixed-modulus exclusion transmutation;
* a uniform centered moving-target discrepancy estimate;
* an independent codimension-two/resultant theorem;
* a mixed-prime geometric or Fourier large-sieve estimate.

The density \(|Z_p|/p\) supplies the expected main term only after such a spatial theorem is proved. It cannot serve as a substitute for it.

---

## 10. Recommended exact residual for the project

The clean next named residual is:

> **[LGTQ-MTD] Uniform moving-target discrepancy.** For every fixed \(C>0\), uniformly in the admissible \(q\) and \(1\le h_1,h_2\le C\log X\), the centered Apéry target indicator has total correlation \(o(X/(\log X)^3)\) on the corresponding prime-triple carrier.

In the one-indicator form:

\[
\max_{q,h_1,h_2}
\left|
\sum_{p\in\mathscr P_{h_1,h_2}(X)}
\left(
\mathbf 1_{Z_p}(\rho_{q;h_1,h_2}(p))
-\frac{|Z_p|}{p}
\right)
\right|
=o\!\left(\frac{X}{(\log X)^3}\right).
\]

If the exact LGTQ definition contains three zero indicators, replace this by the full centered multilinear expansion and require the same bound for every nonempty centered subproduct. The fully centered term is the genuine new obstacle.

This residual is strictly narrower than a global independence theorem: it asks only for cancellation along the specific LGTQ target graph and only on the prime-triple carrier. It is also exactly strong enough: the bounded-fiber main term is then negligible.

---

## Final verdict

The Apéry zero set being small is a **potential** codimension-one fact, not an automatically usable sieve dimension. It becomes an extra logarithmic exponent only when the zero classes are realized as genuinely new congruence exclusions for many fixed auxiliary primes. In LGTQ-3 they occur instead as a moving-target membership condition. Cardinality, Poisson marginal counts, reflection, no adjacency, and one-prime Parseval data all leave open perfect alignment with the target.

Therefore the proposed “sieve theory + Poisson size law” bypass does not close LGTQ-3. The extra codimension must come from a theorem about **where the roots are relative to the LGTQ target**, quantitatively and uniformly across the prime-triple family.

## References used for the sieve audit

* E. Bombieri, “On the large sieve,” *Mathematika* 12 (1965), 201–225.
* P. X. Gallagher, “The large sieve,” *Mathematika* 14 (1967), 14–20.
* P. X. Gallagher, “A larger sieve,” *Acta Arithmetica* 18 (1971), 77–81.
* D. R. Heath-Brown, “Lectures on sieves,” 2002 lecture notes, arXiv:math/0209360.
* Internal project state: `problems/3.2/proof.tex`, `problems/3.2/UNDERSTANDING.md`, `SUBMIT/3.2/README.md`, and `drops/Q7901-eed6cfc4.md`.
