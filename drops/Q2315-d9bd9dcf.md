ANSWER Q2315 d9bd9dcf

# Q2315: Hensel lifting argument for the resultant

The proposed argument is **not valid as stated**. The obstruction is that reduction modulo p does not preserve the individual roots of a characteristic-zero polynomial in the way required. Hensel lifting goes in the opposite direction: a simple root modulo p lifts uniquely to a p-adic root, but a characteristic-zero algebraic root does not become the same element of \(\overline{\mathbf F}_p\) under reduction unless one has first chosen an integral model and reduced its algebraic integer.

## 1. What Hensel actually gives

Let \(F(T)\in \mathbf Z[T]\) be a primitive integral multiple of \(A_s\). Suppose \(p\nmid \operatorname{disc}(F)\). Then:

* \(F\bmod p\) is squarefree.
* Every root \(\bar\alpha\in\overline{\mathbf F}_p\) of \(F\bmod p\) lifts uniquely to a root of \(F\) in the maximal unramified extension of \(\mathbf Q_p\).

The condition is
\[
F(\bar\alpha)=0\pmod p,\qquad F'(\bar\alpha)\neq0\pmod p.
\]

Equivalently the bad primes are those dividing
\[
\operatorname{disc}(F).
\]

There is no bound of the form \(p>C(r,s)\) alone which guarantees this. A valid explicit exceptional bound is
\[
C(r,s)=\max\{\text{primes dividing the denominator},\text{primes dividing }\operatorname{disc}(\widetilde A_s),
\text{primes dividing the leading coefficient}
\},
\]
where \(\widetilde A_s\in\mathbf Z[T]\) is the chosen primitive integer polynomial.

For every
\[
p\nmid \operatorname{lc}(\widetilde A_s)\operatorname{disc}(\widetilde A_s),
\]
we only know that the mod-p roots are distinct. We do **not** know that they are reductions of the real roots in a canonical way.

## 2. Why the proposed resultant argument fails

The statement

\[
\text{roots of }A_s\bmod p = \text{characteristic-zero roots reduced mod }p
\]

is false in general. The characteristic-zero roots are algebraic numbers. They live in number fields, not in \(\mathbf Q\), and their reductions depend on the choice of prime ideal above p in the splitting field.

Even when the roots are simple and all real, the Frobenius action can permute them after reduction. Hensel gives a bijection between roots in a fixed p-adic embedding and roots modulo p, but not an identification with the original ordered real roots.

Therefore interlacing over \(\mathbf R\) does not imply
\[
A_r(\alpha)\not\equiv0\pmod p
\]
for every root \(\alpha\) of \(A_s\bmod p\).

## 3. Correct resultant criterion

For fixed r,s:

\[
p\mid \operatorname{Res}(\widetilde A_r,\widetilde A_s)
\]

if and only if the reductions have a common root in \(\overline{\mathbf F}_p\).

The correct finite exceptional set is exactly
\[
\{p:p\mid \operatorname{Res}(\widetilde A_r,\widetilde A_s)\}.
\]

The fact that the characteristic-zero resultant is nonzero only proves that this set is finite for fixed r,s.

## 4. Consequence for the P3.2 application

Here r and s depend on p. Hence the resultant is not fixed. The argument does not reduce the problem to finitely many exceptional primes because the integer
\[
\operatorname{Res}(\widetilde A_r,\widetilde A_s)
\]
changes with p.

A possible route would require a new uniform estimate such as
\[
 v_p(\operatorname{Res}(\widetilde A_r,\widetilde A_s))=0
\]
for all p with \(p>2\max(r,s)+1\), but this is precisely the missing finite-characteristic coprimality theorem. It cannot be obtained from Hensel lifting alone.

## 5. Precise Hensel condition

The strongest correct statement is:

Let \(F=\widetilde A_s\). If
\[
p\nmid \operatorname{lc}(F)\operatorname{disc}(F),
\]
then \(F\bmod p\) has exactly \(s\) distinct roots over \(\overline{\mathbf F}_p\), and each root corresponds uniquely to a root of F in an unramified extension of \(\mathbf Q_p\).

This does not imply preservation of the real-root labeling or nonvanishing of \(A_r\) on those reductions.

Hence the proposed \(C(r,s)\)-argument is invalid; the missing step is a genuine uniform mod-p resultant bound, not a Hensel lifting estimate.
