ANSWER Q5681 faa16aeb

# P3.2 mesoscopic block audit: exact saturation and the minimal two-characteristic lemma

## Executive verdict

Write

\[
B_m=\sum_{k=0}^m\binom mk^2\binom{m+k}{k}^2
\]

for the zeta(3) Apéry numbers (the `A_j` of the question). Let
\(H=n^{1/3}\), and let `I` be an interval of length `H` on which both

\[
a=\left\lfloor\frac nq\right\rfloor
\]

and the folded branch are fixed. Put

\[
r(q)=n-aq,\qquad j(q)=\min(r(q),q-1-r(q)).
\]

I do **not** obtain the desired uniform estimate

\[
\log\prod_{\substack{q\in I\cap\mathbb P\\q\mid B_{j(q)}}}q=o(H).
\tag{0.1}
\]

The current repository does, however, determine very sharply what a proof
must add.

1. On one branch the folded index is exactly affine:
   \[
   j(q)=n-aq
   \quad\text{or}\quad
   j(q)=(a+1)q-n-1.
   \tag{0.2}
   \]
   The quotient interval has length at least `H`, so
   \(a(a+1)\le n/H=H^2\). Consequently
   \[
   a\le H,
   \qquad
   q>\frac n{a+1}\ge(1-o(1))H^2,
   \tag{0.3}
   \]
   and throughout the block \(\log q\asymp\log n\). Thus (0.1) is
   equivalent to
   \[
   K(n,I):=\#\{q\in I\cap\mathbb P:q\mid B_{j(q)}\}
   =o(H/\log n).
   \tag{0.4}
   \]

2. Every block target divides the **single fixed integer** `B_n`:
   Apéry--Lucas and reflection give
   \[
   B_n\equiv B_aB_{r(q)}\equiv B_aB_{j(q)}\pmod q.
   \tag{0.5}
   \]
   Hence
   \[
   \prod_{q\in T(n,I)}q\mid B_n.
   \tag{0.6}
   \]
   This removes the moving remainder from the characteristic-zero
   statement, but the height \(\log B_n\sim n\log(17+12\sqrt2)\) is much
   larger than `H`; ordinary size or divisor-count estimates give no
   block saving.

3. The most promising fixed-parameter Racah construction is now ruled out
   **exactly in arbitrary target rank**. For
   \(X_n=n(n+1)\), set
   \[
   \phi_m(X)=\frac{\prod_{u=0}^{m-1}(X-u(u+1))}{(m!)^2},
   \qquad
   S_M(n)=\sum_{m=0}^{M}\phi_m(X_n)^2.
   \tag{0.7}
   \]
   Since
   \(\phi_m(X_n)=\binom nm\binom{n+m}{m}\), a target prime `q`
   with \(M=(q-1)/2\) satisfies
   \[
   q\mid S_M(n).
   \tag{0.8}
   \]
   If target primes \(q_1<\cdots<q_k\) are in the block and
   \(P_i=q_1\cdots q_i\), then
   \[
   P_i^2\mid S_{M_{i+1}}(n)-S_{M_i}(n).
   \tag{0.9}
   \]
   The resulting prefix lattice has exact Smith form
   \[
   \operatorname{diag}
   (q_1,P_1^2q_2,\ldots,P_{k-1}^2q_k).
   \tag{0.10}
   \]
   After the universal node squares and the long prefix direction are
   removed, the boundary module has **unit** maximal Fitting ideal. The
   target product survives only as the index
   \(\operatorname{diag}(1,P_k)\) of two long coordinates. Thus linear
   cutoff identities, normalized differences, Smith/Fitting minors,
   resultants and subresultants cannot produce a primitive scalar carrying
   all targets.

4. Keeping the long coordinate does retain the targets, but actual-state
   positivity makes it exponentially too tall. If `Q` is the left scale
   of `I` and `m` is about `Q/2`, then
   \[
   \binom nm\ge(n/m)^m,
   \qquad
   \binom{n+m}{m}\ge((n+m)/m)^m.
   \tag{0.11}
   \]
   Since \(n\asymp aq\), one boundary square has logarithm
   \(\gg q\ge(1-o(1))H^2\). Removing that square leaves a rational
   Racah factor of height \(O(d\log n)\) for cutoff gap `d`, but the exact
   elimination calculation shows that the later target does **not** divide
   this short factor. This is the height/selectivity saturation, not merely
   a weak estimate.

5. The recurrence, reflection, companion jets, modular forms and present
   sheaf technology all stop at the same cross-characteristic wall. The
   weakest concrete new input I can isolate is a **localized adjacent
   collision estimate**. For every fixed `A>0`, prove
   \[
   \boxed{
   C_A(n,I):=
   \#\left\{\begin{array}{l}
   q<\ell\text{ in }I\cap\mathbb P,\quad
   \ell-q\le A\log n,\\
   q\mid B_{j(q)},\quad \ell\mid B_{j(\ell)}
   \end{array}\right\}
   =o_A(H/\log n)
   }
   \tag{LAC}
   \]
   uniformly for adversarial `n`, `I`, quotient `a`, and branch. This is
   strictly weaker than (0.4), but it implies it by the adjacent-gap
   packing argument proved below.

The short interval therefore gives a clean final target, but it does not
by itself create a nonzero carrier of height `o(H)`. A successful proof
must establish actual-Apéry cancellation between **different residue
characteristics**, not another one-prime zero-fibre theorem.

---

## 1. Repository state audited

This audit uses current `main` at

```text
2eb2e688d33442f32068c3d4051b21462acf5ef1
Merge remote-tracking branch 'origin/main'
```

and in particular the following current results in
`problems/3.2/research/working_notes/Q32_SEPARATION_ANALYSIS.md`:

- the exact finite Racah seam at commit `590c2a2`;
- the companion Frobenius-jet frontier at `11cd0de`;
- exact adjacent-cutoff saturation at `86be510`;
- arbitrary-rank multi-target saturation at `f3e242f`;
- the smooth-radical reduction at `33f381e`;
- global folded-boundary carriers and the corrected sparse-exception
  theorem at `ca0b248`.

The load-bearing point for the present question is that the multi-target
cutoff theorem is symbolic and valid for arbitrary target rank. It is not
based on the finite triple at `n=321`; that triple is only a regression
check.

---

## 2. Exact geometry of one quotient/branch block

Let the quotient be `a`. Its full prime-coordinate interval is

\[
\frac n{a+1}<q\le\frac na.
\tag{2.1}
\]

If this interval contains a subinterval of length `H`, then

\[
\frac n{a(a+1)}\ge H,
\qquad
 a(a+1)\le\frac nH=n^{2/3}=H^2.
\tag{2.2}
\]

This proves (0.3). In particular all candidate primes satisfy

\[
q\ge\frac n{H+1}=H^2(1+O(H^{-1}))^{-1}.
\tag{2.3}
\]

There are two affine branches.

### 2.1 Direct branch

Here `r(q)` is already folded:

\[
j(q)=r(q)=n-aq.
\tag{2.4}
\]

Therefore

\[
q\mid n-j(q).
\tag{2.5}
\]

### 2.2 Reflected branch

Here

\[
j(q)=q-1-r(q)=(a+1)q-n-1,
\tag{2.6}
\]

and hence

\[
q\mid n+j(q)+1.
\tag{2.7}
\]

These are the exact folded linear carriers already globalized in the
repository. Across an interval of `q`-length `H`, the folded index changes
by at most `aH` on the direct branch and `(a+1)H` on the reflected branch.
That can be as large as `H^2`; the folded indices are not confined to a
short characteristic-zero window even though the prime cutoffs
\((q-1)/2\) are.

Since every `q` in `I` satisfies

\[
\frac23\log n+O(1)\le\log q\le\log n,
\tag{2.8}
\]

we have

\[
\frac23K(n,I)\log n+O(K(n,I))
\le \log\prod_{q\in T(n,I)}q
\le K(n,I)\log n.
\tag{2.9}
\]

Thus the logarithmic and counting formulations are genuinely equivalent;
no prime number theorem in the short interval is being assumed.

A Brun--Selberg upper bound gives only

\[
\#(I\cap\mathbb P)=O(H/\log n),
\tag{2.10}
\]

so the desired result needs a factor tending to zero beyond the ambient
prime count.

---

## 3. The fixed coefficient reduction: useful but still critical-height

For `n=aq+r` with `a,r<q`, the p-Lucas congruence is

\[
B_n\equiv B_aB_r\pmod q.
\tag{3.1}
\]

Reflection gives

\[
B_r\equiv B_{q-1-r}=B_j\pmod q.
\tag{3.2}
\]

Consequently every target divides `B_n`, proving (0.6). Hence

\[
R(n,I):=\prod_{q\in T(n,I)}q
\mid \operatorname{rad}_{I}(B_n).
\tag{3.3}
\]

This is the cleanest characteristic-zero formulation of the block. It
uses the distinguished Apéry initial state, not a generic recurrence.
Nevertheless

\[
\log B_n=n\log(17+12\sqrt2)-\frac32\log n+O(1),
\tag{3.4}
\]

whereas `H=n^(1/3)`. Therefore

\[
\log R(n,I)\le\log B_n
\tag{3.5}
\]

is off by a factor `n/H=n^(2/3)`. Likewise

\[
K(n,I)\le\frac{\log B_n}{\log(\inf I)}\ll n/\log n
\tag{3.6}
\]

is much weaker than `H/log n`.

The known theorem of Luca--Shparlinski on prime factors of Apéry numbers
is density-one and gives lower bounds for the number and size of prime
factors. It does not upper-bound the squarefree part in a prescribed
moving interval for every `n`. No known smooth-part theorem for
P-recursive/G-function coefficients gives

\[
\log\operatorname{rad}_{I}(B_n)=o(H)
\tag{3.7}
\]

with the present uniformity.

A stronger but target-blind sufficient lemma would be:

\[
\boxed{
\sup_{\substack{I:\ |I|=n^{1/3}\\
                  I\subset(n/(a+1),n/a]}}
\log\operatorname{rad}_{I}(B_n)
\ll \frac{n^{1/3}}{(\log n)^\delta}
}
\tag{MSR}
\]

for some fixed `delta>0`. This would solve every localized block at once,
but it is currently as new as the original problem.

---

## 4. Exact Racah cutoff formulation

Put

\[
X_n=n(n+1),
\qquad
\phi_m(X)=\frac{\prod_{u=0}^{m-1}(X-u(u+1))}{(m!)^2}.
\tag{4.1}
\]

At the characteristic-zero node `X_n`,

\[
\phi_m(X_n)=\binom nm\binom{n+m}{m}.
\tag{4.2}
\]

For an odd prime `q`, let `M_q=(q-1)/2` and

\[
S_q(n)=\sum_{m=0}^{M_q}\phi_m(X_n)^2.
\tag{4.3}
\]

The Racah interpolation gives

\[
S_q(n)=H_q(X_n)\equiv B_{j(q)}\pmod q.
\tag{4.4}
\]

Thus

\[
q\in T(n,I)\quad\Longleftrightarrow\quad q\mid S_q(n).
\tag{4.5}
\]

This is the apparently favorable feature of a short prime block: all
`S_q` are nearby prefixes of one positive characteristic-zero sum.
Unfortunately, the exact prefix module is saturated.

### 4.1 Nested square divisibility

Order the target primes as

\[
q_1<\cdots<q_k,
\qquad M_i=(q_i-1)/2,
\qquad P_i=q_1\cdots q_i.
\tag{4.6}
\]

For `m>M_i`, the product defining `phi_m(X_n)` contains the factor

\[
X_n-j(q_i)(j(q_i)+1)\equiv0\pmod {q_i},
\tag{4.7}
\]

because `j(q_i)<=M_i`. Hence every later summand is divisible by
`q_i^2`, and therefore

\[
P_i^2\mid S_{q_{i+1}}(n)-S_{q_i}(n).
\tag{4.8}
\]

Together with `q_i|S_{q_i}`, these congruences have the exact integral
parametrization

\[
S_{q_i}=P_kt-\sum_{h=i}^{k-1}P_h^2c_h.
\tag{4.9}
\]

The unimodular coordinate change from prefixes to the final prefix and
successive differences gives

\[
P_k\mathbb Z\oplus P_{k-1}^2\mathbb Z
\oplus\cdots\oplus P_1^2\mathbb Z,
\tag{4.10}
\]

and sorting valuations yields the Smith form (0.10).

### 4.2 The target content disappears from the primitive boundary module

The actual cutoff graph has two long coordinates. In the localized ring
at the cutoff denominators it can be written

\[
S_{q_i}=s+WQ_i,
\tag{4.11}
\]

where `Q_1=0`. The target equations are

\[
s+WQ_i=q_it_i.
\tag{4.12}
\]

Their presentation rows are

\[
(1,Q_i,0,\ldots,-q_i,\ldots,0).
\tag{4.13}
\]

A maximal minor using the `s` column and all quotient columns except the
`i`-th is `P_k/q_i`; their gcd is one. Therefore the maximal Fitting ideal
is the unit ideal. Equivalently, after CRT parametrization the projection
onto `(s,W)` has Smith form

\[
\operatorname{diag}(1,P_k).
\tag{4.14}
\]

The complete target product is an index of the two long coordinates, not
primitive torsion in the equation module. The divided boundary variables

\[
c_i=\frac{S_{q_{i+1}}-S_{q_i}}{P_i^2}
\tag{4.15}
\]

have an identity block after the universal squares are removed.

This proves the following scoped no-go.

> **Cutoff-module no-go.** No carrier obtained by finite or growing linear
> elimination from the Racah prefixes, their contiguous differences,
> normalized boundary factors, Smith/Fitting minors, ordinary resultants or
> subresultants can retain every later target prime after the long prefix
> coordinates are removed.

This does not exclude a genuinely new nonlinear identity using the actual
initial state. It proves that such an identity is not hidden in the
existing cutoff algebra.

### 4.3 Why retaining the long coordinate is too expensive

Take `m` near `q/2`. Since `n` is approximately `aq`,

\[
\binom nm\ge(n/m)^m\gg(2a)^m,
\qquad
\binom{n+m}{m}\ge((n+m)/m)^m\gg(2a+1)^m.
\tag{4.16}
\]

Thus

\[
\log \phi_m(X_n)^2\gg q.
\tag{4.17}
\]

By (2.3), `q>=H^2(1-o(1))`; hence each genuine long prefix or boundary
square has logarithmic height `Omega(H^2)`, not `o(H)`.

For two nearby cutoffs, one can factor

\[
\frac{S_\ell-S_q}{q^2}
=\left(\frac{\phi_{M_q+1}(X_n)}q\right)^2
 Q_d(X_n),
\tag{4.18}
\]

where `d=(ell-q)/2` and

\[
h(Q_d(X_n))=O(d\log n).
\tag{4.19}
\]

But modulo `ell` the target equation is

\[
S_q+\phi_{M_q+1}(X_n)^2Q_d(X_n)\equiv0,
\tag{4.20}
\]

monic in the long prefix `S_q`. Its elimination ideal in the short
boundary coordinate is zero. In particular `ell` need not divide `Q_d`.
The exact pure-cross example `(n,q,ell)=(321,193,211)` in the repository
exhibits this nonvanishing, but (4.20) is the general proof.

Thus the only low-height Racah coordinate is not target-selective, and the
selective coordinate is exponentially tall.

---

## 5. Recurrence and reflection on the affine track

If `q<ell` are in the same branch and `h=ell-q`, then

\[
|j(\ell)-j(q)|=
\begin{cases}
a h,&\text{direct branch},\\
(a+1)h,&\text{reflected branch}.
\end{cases}
\tag{5.1}
\]

The Apéry continuant for an index gap `d` gives a low-height certificate
only when **one characteristic** sees two zeros. For example, if

\[
q\mid B_{j(q)}
\quad\text{and}\quad
q\mid B_{j(\ell)},
\tag{5.2}
\]

then `q` divides the corresponding continuant of degree `3(d-1)` and
height `O(d log n)`. The same holds with `ell` in place of `q`.

A generic block collision is instead pure cross:

\[
q\mid B_{j(q)},\qquad \ell\mid B_{j(\ell)},
\tag{5.3}
\]

while

\[
q\nmid B_{j(\ell)},\qquad \ell\nmid B_{j(q)}.
\tag{5.4}
\]

Then the same continuant is a unit modulo both moving primes. The
recurrence has not produced a common divisor; it merely transfers a free
projective coordinate separately in each characteristic.

The companion Frobenius jets make this limitation exact. At a target
`q|B_j`, the divided Apéry and companion coordinates are both multiples
of

\[
\omega_{q,k,j}=B_j/q+kW_j\pmod q.
\tag{5.5}
\]

They remain projectively rank one. Even the conjectural projective
congruence through `q^3` would propagate the same direction; its first
defect on the target locus is a unit rather than a second zero equation.
Hence adding the rational companion does not give a block selector.

Reflection likewise creates the partner zero in the same row `q`, but the
partner lies in a different outer integer column. It supplies no relation
between two different primes in the fixed block.

---

## 6. Modular and theta-cycle audit

The level-six modular parametrization gives `B_j` as the constant term of
a weakly holomorphic weight-four form whose pole order is `j`. Clearing
the pole requires weight `Omega(j)` by the valence formula. Since the
folded indices in the core can be of order `q`, this already has height
far beyond `H`.

There is a sharper characteristic-`q` Sturm formulation. For each residue
`j`, the repository constructs

\[
Q_{q,j}\in M_{2(q-1)}(\Gamma_0(6);\mathbb Z_{(q)})
\tag{6.1}
\]

such that

\[
q\mid B_j
\quad\Longleftrightarrow\quad
Q_{q,j}\equiv0\pmod q.
\tag{6.2}
\]

The selector is the optimal polynomial

\[
1-(\Theta_t-j)^{q-1};
\tag{6.3}
\]

any polynomial over `F_q` equal to one at `j` and zero at all other
residues has degree at least `q-1`. The `Gamma_0(6)` Sturm bound is exactly
`2(q-1)`.

For `K` targets in the block, multiplying the characteristic-specific
forms gives total weight at least

\[
2\sum_{q\in T(n,I)}(q-1)
\gg K H^2.
\tag{6.4}
\]

At the critical density `K~H/log n`, this is much larger than `H`.
Taking a common weight via least common multiples of the `q-1` only makes
the situation worse. Thus ordinary Sturm, theta cycles and products of
the local Hasse projectors do not yield a sub-`H` characteristic-zero
carrier.

---

## 7. The genuine sheaf simplification—and its exact obstruction

For each prime `q`, the value `B_j mod q` is a finite Mellin transform of
the fixed Apéry/Beukers--Peters K3 trace sheaf. The resulting Mellin
cohomology is a rank-two weight-three object. On the present fixed branch,
let

\[
m=n-a.
\tag{7.1}
\]

Since `q=1 mod(q-1)`, (0.2) gives

\[
j(q)\equiv
\begin{cases}
m\pmod{q-1},&\text{direct},\\
-m\pmod{q-1},&\text{reflected}.
\end{cases}
\tag{7.2}
\]

Thus all primes in the block sample

\[
\chi_q=\omega_q^{\sigma m},
\qquad \sigma\in\{1,-1\}
\tag{7.3}
\]

with one fixed integer exponent `m`. This is a real gain over the raw
moving-index formulation.

It is not, however, a fixed compatible system. The character order is

\[
d_q=\frac{q-1}{\gcd(q-1,m)},
\tag{7.4}
\]

which is generally unbounded, and the cyclotomic coefficient field and
selected prime above `q` change with `q`. The situation is

```text
one changing finite field + one changing high-order character
for each prime,
```

not

```text
one fixed finite field + many arguments,
```

nor

```text
one fixed compatible system + many Frobenius primes.
```

This is why the standard results do not apply:

- Katz's finite-field Mellin equidistribution varies characters in one
  finite-field family (or over controlled finite extensions); it is a
  vertical theorem and does not control this adversarial diagonal.
- Fouvry--Kowalski--Michel's prime trace-function estimates have a fixed
  trace function/modulus and vary prime arguments. Here the modulus is the
  prime being counted.
- Kowalski's large sieve for Frobenius requires a fixed algebraic family or
  compatible system and fixed finite quotients. Neither is present after
  the unbounded-order Kummer twist.

There is a second, independent obstruction: exact vanishing is a local
`q`-adic event, not a continuous Sato--Tate test. The reduced polynomial
representing the delta function at zero on `F_q` is

\[
\delta_{0,q}(x)=1-x^{q-1}.
\tag{7.5}
\]

It has minimal degree `q-1`: a polynomial of degree below `q-1` cannot
have `q-1` distinct nonzero roots and value one at zero. Thus any direct
exact-zero sheaf/projector has conductor or algebraic degree growing like
`q`, even though the underlying Mellin sheaf before zero detection has
bounded local conductor.

Finally, in terms of the prime scale `Q~n/a`, the block length satisfies

\[
Q^{1/3+o(1)}\le H\le Q^{1/2+o(1)}.
\tag{7.6}
\]

A uniform exact local-limit theorem over prime intervals this short would
be substantially beyond the standard bounded-conductor prime trace
estimates even if the compatible-system problem were repaired.

---

## 8. Constant term, diagonal and determinant routes

The diagonal representation

\[
B_n=\operatorname{CT}\Lambda^n
\tag{8.1}
\]

is valuable because all block targets divide this one coefficient. But
its coefficients are nonnegative, and the Racah prefix representation is
a sum of positive squares. There is no Archimedean cancellation in the
long coordinate to exploit.

A determinant made from nearby prefixes or nearby Apéry states faces the
same dichotomy:

- retaining a row/column on which a target is forced gives the full CRT
  index and exponential height;
- dividing the universal node factors or eliminating the long state
  gives a primitive determinant which is a unit at the later targets.

This is exactly what the Smith/Fitting calculations prove. Hadamard or
ordinary determinant bounds cannot reverse it.

The asymptotic approximation

\[
B_j=C\lambda^j j^{-3/2}(1+O(j^{-1})),
\qquad \lambda=17+12\sqrt2,
\tag{8.2}
\]

can make normalized real determinants small, but it does not make them
divisible by the different target primes. Cross-weighting to impose that
divisibility inserts the unknown CRT product into the coefficients. No
known irrationality or algebraic-approximation estimate for `zeta(3)`
provides the growing-order cancellation required to escape this ledger.

---

## 9. The critical fixed-gap carrier and why the block remains critical

For a fixed even prime gap `h`, the Selberg upper-bound sieve gives

\[
\#\{q\in I:q,q+h\text{ prime}\}
\ll \mathfrak S(h)\frac{H}{\log^2 n}.
\tag{9.1}
\]

Therefore the product of all ambient prime pairs at that one gap has
logarithmic height

\[
O\left(\mathfrak S(h)\frac{H}{\log n}\right)=o(H).
\tag{9.2}
\]

This is a genuine nonzero sub-`H` carrier containing every target pair at
that fixed gap. It is the closest unconditional positive result.

A dense target set, however, can use all gaps up to `A log n`. Since

\[
\sum_{h\le A\log n}\mathfrak S(h)\ll_A\log n,
\tag{9.3}
\]

multiplying the fixed-gap carriers costs

\[
O_A(H),
\tag{9.4}
\]

exactly the critical height, not `o(H)`. The short block has not saved a
logarithm. Any actual-Apéry density saving tending to infinity relative
to the ambient prime-pair sieve would close the argument.

---

## 10. The smallest concrete new lemma

Order the target primes:

\[
q_1<q_2<\cdots<q_K.
\tag{10.1}
\]

Fix `A>0` and put `L=A log n`. The total span is at most `H`, so the
number of adjacent gaps exceeding `L` is at most `H/L`. Every remaining
adjacent gap is counted by `C_A(n,I)` in `(LAC)`. Hence

\[
K-1\le\frac{H}{A\log n}+C_A(n,I).
\tag{10.2}
\]

If `(LAC)` holds, then

\[
\limsup_{n\to\infty}
\sup_I\frac{K(n,I)\log n}{H}\le\frac1A.
\tag{10.3}
\]

Since `A` is arbitrary, this proves (0.4). This quantifier order is safe:
for a prescribed epsilon choose fixed `A>2/epsilon`, then let `n` tend to
infinity in the uniform `(LAC)` estimate.

`(LAC)` is strictly weaker than full block compression. It ignores:

- all nonadjacent target pairs;
- all gaps larger than `A log n`;
- the actual number of isolated targets.

It is also weaker than the full shell pair-energy theorem. It is the
smallest pointwise two-characteristic statement currently isolated.

### 10.1 Exact analytic zero-detector form

For `x in F_q`, write

\[
\delta_q(x)=\frac1q\sum_{u\in\mathbb F_q}e_q(ux).
\tag{10.4}
\]

Then `(LAC)` is exactly

\[
\sum_{1\le h\le A\log n}
\sum_{\substack{q,q+h\in I\cap\mathbb P}}
\delta_q(B_{j(q)})
\delta_{q+h}(B_{j(q+h)})
=o_A(H/\log n).
\tag{10.5}
\]

Thus a quantitative estimate

\[
\ll_A \frac{H}{(\log n)^{1+\delta}}
\tag{10.6}
\]

for any fixed `delta>0` would suffice. This is an exact two-modulus local
limit, not a moment of normalized complex traces.

### 10.2 Geometric form

Let `M_{q,chi}` denote the rank-two weight-three Mellin Frobenius object
attached to the Apéry K3 trace sheaf. A sufficient geometric theorem is

\[
\sum_{h\le A\log n}
\sum_{\substack{q,q+h\in I\cap\mathbb P}}
\mathbf 1_{\mathrm{nonord}(M_{q,\omega_q^{\sigma m}})}
\mathbf 1_{\mathrm{nonord}(M_{q+h,\omega_{q+h}^{\sigma m}})}
=o_A(H/\log n),
\tag{10.7}
\]

uniformly for `m=n-a`, `a<=H`, and both signs. Proving (10.7) would
require a new horizontal Kummer--Mellin large sieve with:

1. two different residue characteristics in each correlation;
2. unbounded character order;
3. exact `q`-adic nonordinarity, not continuous trace equidistribution;
4. uniformity in the adversarial exponent `m`;
5. prime intervals down to the `Q^(1/3)` scale.

No standard Katz/Deligne, Chebotarev, FKM, or Frobenius-large-sieve theorem
has these quantifiers.

---

## 11. What would count as a genuine carrier breakthrough

A direct carrier solution would be an explicitly defined nonzero integer
`C_{n,I}` satisfying

\[
R(n,I)\mid C_{n,I},
\qquad
\log|C_{n,I}|=o(H)
\tag{11.1}
\]

uniformly. The current audits impose three necessary features.

1. It cannot lie in the primitive linear Racah cutoff/Fitting algebra;
   that algebra has unit target content after saturation.
2. It cannot be a product of the per-prime Sturm/Hasse projectors; their
   degree and weight are `Omega(q)`.
3. It cannot batch prime-local marked traces by ordinary CRT; the CRT
   lattice index is exactly the unknown target product.

Therefore a viable carrier must use a new **initial-state-specific
nonlinear cancellation involving the long prefix**, with coefficients
whose own height does not encode the target CRT modulus. No such identity
is presently known.

A characteristic-zero alternative is the mesoscopic radical theorem
`(MSR)`. That would avoid sheaves entirely, but it is a new every-`n`
prime-factor theorem for one Apéry coefficient, far stronger than current
density-one results.

---

## 12. Strategic conclusion

The `H=n^(1/3)` localization is real and useful:

```text
quotient a is at most H;
prime scale is at least H^2;
branch index is affine;
Mellin exponent n-a is fixed up to sign;
prime cutoffs lie in one interval of length H/2.
```

But each apparent compression stops sharply:

```text
fixed B_n:                     height Theta(n);
Racah long prefix:             height at least exp(Omega(H^2));
normalized Racah boundary:     short, but not target-divisible;
recurrence continuant:         sees two zeros only in one characteristic;
Sturm/theta projector:         degree and weight Omega(q);
Mellin sheaf:                  bounded locally, not a fixed global system;
exact zero detector:           minimal degree q-1;
all logarithmic prime gaps:    ambient carrier height O(H), exactly critical.
```

Accordingly, I do not have a complete proof of the block theorem. The
precise surviving target is `(LAC)`: an `o(1)` density saving for
short-gap **two-characteristic Apéry collisions** inside every adversarial
block. This is weaker than full block compression and strong enough to
finish it. It is also exactly the input not supplied by any current
one-prime, average-over-`n`, bounded-order, or bounded-conductor theorem.

---

## References used for the theorem-scope audit

- F. Beukers and C. A. M. Peters, *A family of K3 surfaces and zeta(3)*,
  J. reine angew. Math. 351 (1984), 42--54.
- J. Stienstra and F. Beukers, *On the Picard--Fuchs equation and the
  formal Brauer group of certain elliptic K3-surfaces*, Math. Ann. 271
  (1985), 269--304.
- N. M. Katz, *Convolution and Equidistribution: Sato--Tate Theorems for
  Finite-Field Mellin Transforms*, Annals of Mathematics Studies 180,
  Princeton University Press, 2012.
- E. Kowalski, *The large sieve, monodromy and zeta functions of curves*,
  J. reine angew. Math. 601 (2006), 29--69; arXiv `math/0503714`.
- E. Fouvry, E. Kowalski and P. Michel, *Algebraic trace functions over
  the primes*, Duke Math. J. 163 (2014), 1683--1736; arXiv `1211.6043`.
- F. Luca and I. E. Shparlinski, *Arithmetic properties of Apéry
  numbers*, J. London Math. Soc. 78 (2008), 545--562.
