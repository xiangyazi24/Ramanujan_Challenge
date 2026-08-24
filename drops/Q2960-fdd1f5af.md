ANSWER Q2960 fdd1f5af

# Verdict

I do **not** get a sublinear characteristic-zero height bound for the live Smith carrier
\[
H_q(n)=\frac{d_{n-q}(J_n)}{\gcd(d_{n-q}(J_n),d_{n-q-1}(J_{n-1}))}
\]
from the present algebra alone. I do get a new exact Apéry-specific reduction that is substantially more concrete than the full determinantal-divisor quotient and does not use generic interlacing or the false normalized-shell identity.

For fixed `q`, define one explicit principal path minor for every integer cut `c` by deleting the vertices
\[
c-1,\;2c-1,\ldots,qc-1
\]
from the Apéry Jacobi matrix. Let `M_{n,q}(c)` be its determinant and let
\[
\Gamma_q(n)=\gcd_{1\le c\le \lfloor n/q\rfloor} M_{n,q}(c),
\qquad
E_q(n)=d_{n-q-1}(J_{n-1}),
\]
\[
\boxed{\quad
\Sigma_q(n)=\frac{\Gamma_q(n)}{\gcd(\Gamma_q(n),E_q(n))}.
\quad}
\]
Then, with no unproved normalization,
\[
\boxed{\quad H_q(n)\mid \Sigma_q(n).\quad}
\]
Moreover, for every prime `p>=5` in the fixed-`q` cell
\[
q=\Big\lfloor\frac np\Big\rfloor,\qquad n=qp+r,\quad 0\le r<p,
\]
one has the exact characteristic-`p` specialization
\[
\boxed{\quad
M_{n,q}(p)\equiv (-1)^q(r!)^3 b_r\pmod p,
\quad}
\]
and consequently
\[
\boxed{\quad
p\mid\Gamma_q(n)\iff p\mid b_r.
\quad}
\]
If `p|b_r`, then additionally `p\nmid E_q(n)`, hence `p|\Sigma_q(n)` and `p|H_q(n)`.

Thus the fixed-`q` moving-prime problem is reduced to the following **single explicit path-gcd height statement**:
\[
\boxed{\log \Sigma_q(n)=o(n)\quad(q\text{ fixed}).}
\]
That statement is sufficient for the desired height theorem, and its large-prime support in the fixed-`q` cell is *exactly* the target support. What remains is a genuine characteristic-zero gcd-height problem for a one-parameter family of factored continuants, not a Smith-interlacing problem.

I do **not** prove this last height estimate below. I also explain why the obvious attempt to bound one fixed cut cannot work: the cut that exposes a target prime is `c=p`, so the useful cut itself moves with the prime.

## Source-state caveat

The GitHub connector currently exposes default branch `main` at commit `47fa0e653f52c4a71e9a8c26b31ca9f66f6bbe86`. On that branch, `problems/3.2/research/proofs/` contains only the three older files `Q32_ALL_CUTOFF_BLOCK_AND_LOSS_PROOF.md`, `Q32_FRANEL_TAIL_DWORK_PROOF.md`, and `Q32_SIGNED_PADE_DEGREE_ONE_NO_GO.md`; the four filenames named in Q2960 are not returned by the connected repository, including by exact-path fetch. I therefore do not claim to have read text the connector did not expose. The derivation below uses the live carrier stated in Q2960 and the current repository's Apéry recurrence/continuant normalization. No repository source files were edited.

# 1. Exact Jacobi normalization

Index rows and columns of `J_n` by `0,...,n-1`. The Apéry Jacobi matrix is the integral symmetric tridiagonal matrix
\[
(J_n)_{t,t}=P(t),\qquad
(J_n)_{t-1,t}=(J_n)_{t,t-1}=t^3\quad(1\le t<n),
\]
where
\[
P(t)=34t^3+51t^2+27t+5.
\]
The sign chosen for the off-diagonal is irrelevant for all principal determinants below; only its square enters the continuant recurrence.

For `a>=0`, `h>=0`, define the interval continuant
\[
K_h(a)=\det J[a,a+h-1],
\]
with the empty determinant `K_0(a)=1`. Then
\[
K_1(a)=P(a),
\]
\[
\boxed{
K_{h+1}(a)=P(a+h)K_h(a)-(a+h)^6K_{h-1}(a).
}
\]
At `a=0`, comparison with the Apéry recurrence gives
\[
\boxed{K_h(0)=(h!)^3b_h.}
\]
Indeed both sides start with `1,5` and satisfy the same renormalized recurrence
\[
B_{h+1}=P(h)B_h-h^6B_{h-1},\qquad B_h=(h!)^3b_h.
\]

There are no zero-denominator or sign issues hidden here: this is an identity in `Z`.

A useful positivity check is
\[
P(t)-t^3-(t+1)^3=4(2t+1)^3>0\qquad(t\ge0).
\]
Hence every finite interval Jacobi block is symmetric strictly diagonally dominant with positive diagonal, so every `K_h(a)>0`. In particular all gcds below are ordinary positive gcds.

# 2. The moving-cut principal minor

Fix `q>=1` and `n>=q`. Put
\[
k=n-q.
\]
For an integer `c` with
\[
1\le c\le\lfloor n/q\rfloor,
\]
delete from `J_n` the `q` vertices
\[
D_c=\{c-1,2c-1,\ldots,qc-1\}.
\]
The remaining principal submatrix has size `n-q=k`. Because the graph of a tridiagonal matrix is a path, deleting those vertices disconnects it into `q+1` intervals. Therefore its determinant factors **in characteristic zero** as
\[
\boxed{
M_{n,q}(c)
 =\left(\prod_{j=0}^{q-1}K_{c-1}(jc)\right)
   K_{n-qc}(qc).
}
\tag{2.1}
\]
This is just block-diagonal determinant factorization; no modular argument and no finite scan enter.

Since `M_{n,q}(c)` is an actual `k x k` minor of `J_n`, if
\[
D_q(n):=d_k(J_n),
\]
then
\[
\boxed{D_q(n)\mid M_{n,q}(c)\quad\text{for every admissible }c.}
\tag{2.2}
\]
Consequently, for
\[
\Gamma_q(n):=\gcd_{1\le c\le\lfloor n/q\rfloor}M_{n,q}(c),
\]
we have
\[
\boxed{D_q(n)\mid\Gamma_q(n).}
\tag{2.3}
\]
This is the first exact reduction: all minors have been replaced by one one-parameter family of factored principal path minors.

# 3. Exact specialization at the candidate cut `c=p`

Let `p` be a prime with
\[
q=\left\lfloor\frac np\right\rfloor,
\qquad n=qp+r,\qquad 0\le r<p.
\]
Then `p` is an admissible value of `c`.

For every fixed `j`, reduction of the continuant recurrence modulo `p` gives
\[
K_h(jp)\equiv K_h(0)\pmod p,
\tag{3.1}
\]
because
\[
P(jp+t)\equiv P(t),\qquad (jp+t)^6\equiv t^6\pmod p.
\]
Thus
\[
K_{p-1}(jp)
 \equiv ((p-1)!)^3b_{p-1}\pmod p.
\]
For every prime `p`,
\[
b_{p-1}\equiv1\pmod p.
\tag{3.2}
\]
One direct proof uses the binomial formula
\[
b_{p-1}=\sum_{m=0}^{p-1}
 \binom{p-1}{m}^2\binom{p-1+m}{m}^2.
\]
For `1<=m<=p-1`, the second binomial coefficient contains exactly one factor `p` in its numerator and none in `m!`, so every `m>0` summand vanishes modulo `p`; the `m=0` summand is `1`.

By Wilson,
\[
((p-1)!)^3\equiv-1\pmod p.
\]
Also
\[
K_r(qp)\equiv K_r(0)=(r!)^3b_r\pmod p.
\]
Substituting `c=p` in (2.1) therefore gives the exact sign
\[
\boxed{
M_{n,q}(p)\equiv(-1)^q(r!)^3b_r\pmod p.
}
\tag{3.3}
\]
Because `r<p`, `r!` is a `p`-unit. Hence
\[
\boxed{
p\mid M_{n,q}(p)\iff p\mid b_r.
}
\tag{3.4}
\]
This is the key Apéry-specific feature of the path family. Generic tridiagonal matrices do not have it.

# 4. Rank proof of the exact target selector

For completeness, the same statement can be seen directly from the mod-`p` Jacobi rank, and this also checks the adjacent denominator.

Assume first `r>0`. In `J_n mod p`, the off-diagonal entries at
\[
p,2p,\ldots,qp
\]
vanish. Thus `J_n mod p` splits into `q` blocks of length `p` and one tail block of length `r`.

Each length-`p` block is, after shifting indices, congruent to the same Apéry block `J_p mod p`. Its determinant is
\[
K_p(0)=(p!)^3b_p\equiv0\pmod p.
\]
All its internal off-diagonal entries are nonzero modulo `p`. Therefore its nullity is exactly one: a kernel vector is determined by its first coordinate, and if that coordinate is zero the tridiagonal recurrence forces the entire vector to vanish.

The tail block has determinant
\[
K_r(qp)\equiv(r!)^3b_r\pmod p,
\]
and, again because its internal off-diagonal entries are units, it has nullity one exactly when `p|b_r` and nullity zero otherwise.

Therefore
\[
\operatorname{nullity}_{\mathbf F_p}(J_n)=
\begin{cases}
q+1,&p\mid b_r,\\
q,&p\nmid b_r.
\end{cases}
\]
Since `k=n-q`, this is equivalent to
\[
\boxed{
p\mid d_k(J_n)\iff p\mid b_r.}
\tag{4.1}
\]
Indeed `p|d_k` iff every `k x k` minor vanishes mod `p`, i.e. iff `rank(J_n)<k`.

If `r=0`, there are exactly `q` length-`p` blocks and no tail; the nullity is `q`, so `p\nmid d_k(J_n)`. This agrees with `b_0=1` and with (3.3), whose tail factor is `K_0=1`.

Now put
\[
E_q(n):=d_{k-1}(J_{n-1}).
\]
Suppose `p>=5` and `p|b_r`. Then `r>=1`.

- If `r=1`, `J_{n-1}` consists mod `p` of exactly `q` length-`p` blocks. Its nullity is `q`, so its rank is `(n-1)-q=k-1`, hence
  \[
  p\nmid E_q(n).
  \]
- If `r>=2`, `J_{n-1}` has the same `q` length-`p` blocks and a tail of length `r-1`. The current P3.2 no-consecutive-zero lemma gives
  \[
  p\nmid b_{r-1}
  \]
  because `p|b_r`. Hence that tail is nonsingular, the total nullity is again `q`, and again
  \[
  \boxed{p\nmid E_q(n).}
  \tag{4.2}
  \]

Thus every target prime contributes to the **excess** from `J_{n-1}` to `J_n`, exactly as required by the adjacent selector.

# 5. The path-gcd upper carrier

Write
\[
D=D_q(n)=d_k(J_n),\qquad E=E_q(n)=d_{k-1}(J_{n-1}),
\]
and recall
\[
H_q(n)=\frac D{\gcd(D,E)}.
\]
From (2.3), `D|Gamma_q(n)`. Define
\[
\Sigma_q(n)=\frac{\Gamma_q(n)}{\gcd(\Gamma_q(n),E)}.
\]
For every prime `ell`,
\[
v_\ell(H_q(n))
 =\max\{v_\ell(D)-v_\ell(E),0\},
\]
whereas
\[
v_\ell(\Sigma_q(n))
 =\max\{v_\ell(\Gamma_q(n))-v_\ell(E),0\}.
\]
Since `D|Gamma_q(n)`,
\[
v_\ell(\Gamma_q(n))\ge v_\ell(D).
\]
Therefore
\[
\boxed{H_q(n)\mid\Sigma_q(n).}
\tag{5.1}
\]
This is only a divisibility. I am **not** asserting the false equality `H_q=Sigma_q`, nor any normalized-shell equality.

Now combine (3.3), (4.1), and (4.2). For every prime `p>=5` in the fixed-`q` cell,
\[
\boxed{
p\mid\Gamma_q(n)\iff p\mid b_{n-qp}.}
\tag{5.2}
\]
Proof:

- If `p|b_r`, then (4.1) gives `p|D`; since `D|Gamma`, `p|Gamma`.
- If `p|Gamma`, then in particular `p|M_{n,q}(p)`; by (3.3), `p|b_r`.

For a target prime, (4.2) says `p\nmid E`, so it survives unchanged in the adjacent normalization:
\[
\boxed{
p\mid\Sigma_q(n)\iff p\mid b_{n-qp}}
\tag{5.3}
\]
for candidate primes in the cell.

This gives the exact chain
\[
\sum_{\substack{p\text{ in fixed-}q\text{ cell}\\p\mid b_{n-qp}}}\log p
\le \log H_q(n)
\le \log\Sigma_q(n),
\tag{5.4}
\]
where the first inequality uses the fact just proved that every target prime divides `H_q(n)`.

# 6. The exact remaining height statement

Equation (5.1) shows that the following theorem would close the live quotient-height route:

> **Path-gcd height theorem (sufficient, not proved here).** For every fixed `q>=1`,
> \[
> \log\left(
> \frac{
> \gcd_{1\le c\le\lfloor n/q\rfloor}
> \left[
> \left(\prod_{j=0}^{q-1}K_{c-1}(jc)\right)K_{n-qc}(qc)
> \right]
> }{
> \gcd\!\left(
> \gcd_{1\le c\le\lfloor n/q\rfloor}
> \left[
> \left(\prod_{j=0}^{q-1}K_{c-1}(jc)\right)K_{n-qc}(qc)
> \right],
> d_{n-q-1}(J_{n-1})
> \right)
> }
> \right)=o(n).
> \]

This is a characteristic-zero statement about a gcd of explicit continuants. It is strictly more structured than the original all-minor Smith definition: every term has the closed path factorization (2.1), and at the moving specialization `c=p` its residue is exactly `(-1)^q(r!)^3b_r`.

Within this one-parameter **principal separator family**, `Gamma_q(n)` is canonical: it is the positive generator of the ideal generated by all `M_{n,q}(c)`. Thus no further Bezout optimization inside this family changes `Gamma`. I do not claim it is globally minimal among all possible minor families; the full determinantal divisor `D` is smaller in general.

# 7. Why this does not yet give a sublinear bound

The factorization (2.1) is exact but does not by itself control the characteristic-zero gcd. A single value `M_{n,q}(c)` has roughly linear-or-worse logarithmic height because its interval continuants cover a total of `n-q` vertices. Taking the gcd over all `c` may destroy almost all of that height, but proving that destruction is exactly the new arithmetic problem.

There is also a precise moving-cut obstruction to the most naive approach. For a candidate prime `p`, the cut which gives the clean unit-times-`b_r` specialization is
\[
\boxed{c=p.}
\]
A fixed finite set of cuts independent of `p` cannot contain `c=p` for all candidate primes. If one starts from a remote cut and Euclidean-descends a trailing continuant toward the zero edge at `p`, the recurrence introduces products of the intervening off-diagonal squares
\[
\prod t^6.
\]
When the distance from the cut to `p` is `Theta(n)`, the logarithm of that fraction-free certificate is also `Theta(n)` (indeed larger before normalization). This is the same height mechanism that makes a fixed-cut continuant certificate inadequate. It is a scoped obstruction to **fixed-cut path certificates**, not a global impossibility theorem for `Gamma_q`.

The new scalar therefore pinpoints what a successful positive proof has to exploit: cancellation/gcd structure **across the moving family in `c`**, not a bound on any one continuant.

# 8. Relation to the rank-one border

The exact border is
\[
J_n=
\begin{pmatrix}
J_{n-1}&(n-1)^3e_{n-2}\\
(n-1)^3e_{n-2}^{T}&P(n-1)
\end{pmatrix}.
\]
For an arbitrary symmetric tridiagonal matrix, this rank-one border does not yield a polynomial bound for the adjacent Smith quotient; the generic counterexample from the preceding audit rules that out. The reduction above uses extra Apéry structure at two essential points:

1. the off-diagonal is exactly `t^3`, so modulo a candidate `p` the path breaks at the arithmetic progression `p,2p,...,qp`;
2. after shifting by `jp`, every interval continuant reduces to the **same** Apéry continuant modulo `p`, giving (3.1) and the exact sign in (3.3).

Those are not consequences of interlacing, and they are absent from the generic counterexample.

# 9. Small-prime and quantifier caveats

The fixed-`q` application of interest has `p>sqrt(n)`, hence `n<p^2` and automatically `q<p`. For fixed `q`, all sufficiently large candidate primes therefore satisfy the block argument above.

The no-consecutive-zero step used to prove `p\nmid E_q(n)` was stated in the current P3.2 source for `p>=5`. The primes `2` and `3` are harmless separately: among the base digits one has
\[
b_0=1,\quad b_1=5,\quad b_2=73,
\]
so there is no `p=2` or `p=3` target digit. In any event, finitely many primes contribute only `O_q(log n)` to a height decomposition.

No prime-power version is asserted here. The rank argument proves exact **prime support** (`p|...`), not an equality of `p`-adic Smith exponents. The divisibility `H_q|Sigma_q`, however, is integral and includes all prime powers automatically; that is why a full characteristic-zero height bound for `Sigma_q` would be enough.

# 10. Bottom line

The positive result of this attack is the exact replacement
\[
\text{all }(n-q)\text{-minors}
\quad\leadsto\quad
\Gamma_q(n)=\gcd_c
\left(\prod_{j=0}^{q-1}K_{c-1}(jc)\right)K_{n-qc}(qc),
\]
followed by the adjacent normalization
\[
\Sigma_q(n)=\Gamma_q(n)/\gcd(\Gamma_q(n),d_{n-q-1}(J_{n-1})).
\]
It satisfies
\[
H_q(n)\mid\Sigma_q(n),
\]
and, in the fixed-`q` large-prime cell,
\[
p\mid\Sigma_q(n)\iff p\mid b_{n-qp}.
\]
The candidate-cut specialization is
\[
M_{n,q}(p)\equiv(-1)^q(r!)^3b_r\pmod p,
\]
with the sign checked.

So the remaining height problem is no longer “control an abstract adjacent Smith quotient.” It is:

\[
\boxed{
\text{prove sublinear height for the normalized gcd of the explicit moving-cut path products (2.1).}
}
\]

I do not see a valid characteristic-zero argument that proves that estimate from the recurrence alone. In particular, bounding one cut, invoking generic Smith interlacing, or identifying a naively normalized shell with `H_q` does not do it. The genuinely live algebraic direction is a cross-`c` gcd/resultant or fraction-free condensation theorem for the family `M_{n,q}(c)`.