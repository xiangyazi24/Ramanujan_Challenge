ANSWER Q415 bdadcf8f

# Executive verdict

A k-th moment strategy is the correct probabilistic framework, but Katz geometric monodromy **does not by itself give pairwise independence across distinct primes**. The obstruction is that the events

\[
E_p(n)=\{b_n\equiv0\pmod p\}=\{n\bmod p\in Z_p\}
\]

are not values of one fixed random variable over a single finite field. They are evaluations of different reductions of the same integral sequence at the same integer n. Katz equidistribution controls variation of Frobenius inside one compatible family over a parameter space; it does not automatically control correlations between the reductions modulo two unrelated primes.

The missing theorem is a genuinely two-prime (or multi-prime) Chebotarev/large-sieve statement for the Apéry sheaf.

# 1. The k-th moment reduction

Write

\[
D(n)=\sum_{\sqrt n<p\le n}1_{E_p(n)}.
\]

The second moment is

\[
\sum_{n\le N}D(n)^2
=\sum_p\#\{n\le N:E_p(n)\}
+2\sum_{p<q}\#\{n\le N:E_p(n)\cap E_q(n)\}.
\]

The diagonal term is harmless:

\[
\sum_p\frac{|Z_p|}{p}N.
\]

The hard term is

\[
R_{pq}=\#\{n\le N:n\bmod p\in Z_p,\ n\bmod q\in Z_q\}.
\]

If p q \le N, CRT gives the expected density:

\[
R_{pq}=\frac{|Z_p||Z_q|}{pq}N+O(|Z_p||Z_q|).
\]

But here p,q>\sqrt N, so generally

\[
pq>N,
\]

and there is no averaging over a full CRT box. The integer interval samples only a thin slice of the product space

\[
\mathbf F_p\times\mathbf F_q.
\]

Thus the second moment problem is exactly a thin-orbit correlation problem.

# 2. What Katz monodromy gives

The Apéry recurrence gives a hypergeometric Picard-Fuchs object. The associated l-adic sheaf has large geometric monodromy (the relevant generic group is Sp(4) in the Katz description).

This gives powerful one-parameter results:

* equidistribution of Frobenius conjugacy classes as the parameter varies;
* square-root cancellation for trace functions;
* independence of different geometric pieces inside one fixed characteristic.

For example, for a fixed prime p one can study

\[
\sum_{x\in\mathbf F_p}\psi(\mathrm{Tr}(\mathrm{Frob}_x))
\]

and obtain Deligne-type bounds.

However, the desired moment requires something different:

\[
\sum_{n\le N}1_{Z_p}(n\bmod p)1_{Z_q}(n\bmod q).
\]

This is a correlation between two reductions:

\[
\rho_p(n),\rho_q(n).
\]

Katz's theorem does not assert that

\[
(\rho_p(n),\rho_q(n))
\]

is equidistributed in

\[
\mathbf F_p\times\mathbf F_q
\]

when n ranges over a short interval.

That is a different arithmetic large sieve problem.

# 3. The right exponential sums

Introduce the Fourier detector

\[
F_p(h)=\sum_{r\in Z_p}e_p(hr).
\]

Then

\[
1_{Z_p}(n)=\frac1p\sum_hF_p(-h)e_p(hn).
\]

The pair correlation becomes

\[
R_{pq}
=\frac1{pq}\sum_{h,k}F_p(h)F_q(k)
\sum_{n\le N}e\left(n\left(\frac hp+\frac kq\right)\right).
\]

The inner sum is large precisely for near-Farey coincidences:

\[
\left\|\frac hp+\frac kq\right\|\lesssim \frac1N.
\]

Therefore the required input is not

\[
|F_p(h)|\ll \sqrt p
\]

or even a Katz bound for one prime. The needed statement is cancellation in the mixed sum

\[
\sum_{p,q}
\sum_{h,k}
F_p(h)F_q(k)
V_N\left(\frac hp-\frac kq\right).
\]

This is analogous to a large-sieve theorem for a family of compatible trace functions.

# 4. Can existing inputs prove the second moment?

No.

Known/project inputs:

## (a) Small zero sets

A bound such as

\[
|Z_p|=O(1)
\]

would give the correct expected size

\[
\lambda=\sum_p\frac{|Z_p|}{p}=O(1),
\]

but does not control pair correlations.

Sparse events can still cluster.

## (b) Fourier size

The trivial estimate

\[
|F_p(h)|\le |Z_p|
\]

is insufficient. It gives no cancellation in the dangerous near-Farey region.

## (c) Ahlgren-Ono modular form congruence

The relation

\[
b_{(p-1)/2}\equiv a_p(f_8)\pmod p
\]

controls the central zero. It gives information about one distinguished point of Z_p, not the full set Z_p or cross-prime correlations.

## (d) F_{p^2} factorization

The factorization of the local polynomial and low-degree factors over \(\mathbf F_{p^2}\) explains the geometry of individual Z_p. It does not imply that Z_p and Z_q behave independently for p\ne q.

Hence none of these inputs closes the k=2 moment.

# 5. What theorem would be enough?

A sufficient statement would be something like:

For every distinct primes p,q in the top range,

\[
R_{pq}
=\frac{|Z_p||Z_q|}{pq}N+O(N^{1-\delta})
\]

on average over p,q, or equivalently

\[
\sum_{p<q}
\left(R_{pq}-\frac{|Z_p||Z_q|}{pq}N\right)
=o(N^2).
\]

For the desired pointwise bound one needs much less than full independence, but one needs some uniform anti-clustering estimate.

A natural formulation is a two-prime Chebotarev large sieve:

\[
\sum_{p,q}\left|
\sum_{n\le N}\phi_p(n\bmod p)\phi_q(n\bmod q)
\right|
\ll N^{1+o(1)}.
\]

Here \(\phi_p\) is the centered zero indicator.

This is exactly the missing horizontal equidistribution statement.

# 6. State of the art comparison

There are several nearby theories:

## Large sieve for Frobenius

Kowalski and others developed large sieve inequalities for Frobenius trace functions. These handle many primes simultaneously, but require a family varying over a geometric parameter in a controlled way.

The Apéry problem is harder because the parameter is the integer index n and the primes are the varying characteristics.

## Chebotarev for compatible systems

Serre-style Chebotarev gives independence of Frobenius distributions in compatible Galois representations when the field and representation are fixed. Here the relevant representations live over different residue characteristics.

## Hypergeometric sheaves

Katz's results give exceptional monodromy and excellent one-prime cancellation. They do not automatically produce cross-characteristic independence.

Therefore the required result appears to be a new synthesis:

* compatible-system large sieve;
* hypergeometric monodromy;
* short-interval distribution of the index n.

# 7. Possible easier route

A full Chebotarev independence theorem may be stronger than necessary.

The empirical Poisson law suggests proving only a moment estimate:

\[
\sum_{n\le N}D(n)^2\ll N.
\]

This asks for average pair decorrelation, not pointwise independence.

Possible routes:

1. Prove a bilinear large sieve for the centered functions

\[
g_p(n)=1_{Z_p}(n\bmod p)-|Z_p|/p.
\]

2. Use the Mellin transform of the Apéry sheaf to show cancellation in

\[
\sum_{n\le N}g_p(n)g_q(n).
\]

3. Use the transfer-matrix orbit description and prove that simultaneous rank-drop conditions modulo two primes are rare.

The last route is particularly Apéry-specific: it avoids pretending that different primes are independent and instead proves that forcing two independent rank drops is arithmetically rigid.

# Final answer

Katz Sp(4) monodromy explains why each individual zero set Z_p looks random, but it does **not** prove pairwise independence of the events

\[
\{n\bmod p\in Z_p\},\quad \{n\bmod q\in Z_q\}
\]

for p\ne q beyond the CRT range.

A k=2 moment proof of

\[
D(n)=o(n/\log n)
\]

requires a new cross-prime dispersion theorem: a two-characteristic large sieve / Chebotarev-type independence result for the Apéry compatible system. This is substantially beyond current Katz equidistribution technology.

The most realistic next target is therefore not full independence, but an averaged pair-correlation bound (HM_2 style). Once that is established, higher moments could plausibly follow by the same strategy.