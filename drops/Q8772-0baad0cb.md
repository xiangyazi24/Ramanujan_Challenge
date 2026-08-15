ANSWER Q8772 0baad0cb

# Legendre/CM benchmark for LGTQ-3: exact Kummer saturation, not cross-prime independence

## Executive verdict

There are two different proposed benchmarks, and neither gives the desired proof of concept.

1. **Literal Legendre Picard–Fuchs coefficient model:** the LGTQ-3 anti-correlation is false, or becomes vacuous under a strict cross-unit definition of “transverse.” Kummer’s theorem gives
   \[
   Z_p=\left\{\frac{p+1}{2},\ldots,p-1\right\},
   \]
   and every zero has valuation exactly two for the squared central binomial coefficient. In a fixed quotient cell, the hit primes are therefore all primes in one explicit interval of positive proportional length. They have maximal short-gap clustering, and an elementary PNT/pigeonhole argument produces a gap pair with
   \[
   \gg_q \frac{n}{(\log n)^3}
   \]
   hit triples—the exact Q8740 obstruction scale.
2. **Fixed CM elliptic curve trace model:** for a fixed admissible gap pair, a Selberg upper-bound sieve does prove the natural estimate
   \[
   T(X;h_1,h_2)\ll \frac{X}{(\log X)^3}.
   \]
   CM/Hecke characters merely impose finite congruence or ideal-splitting conditions and alter the singular-series constant. They do not give the little-
   \(o(X/(\log X)^3)\) saving needed to contradict a Q8740 spike. In simple CM examples the three CM conditions are perfectly correlated.

There is also a basic correction: the Legendre family is **not** a CM family generically. Its \(j\)-invariant is nonconstant. Hecke characters belong to isolated CM fibers, not to the global coefficient sequence \(\binom{2d}{d}^2\).

Thus the answer is:

\[
\boxed{
\text{No CM proof of concept for the needed LGTQ-3 saving arises from Legendre.}
}
\]

The calculation does, however, give a useful negative benchmark: classical carry zeros are too rigid and too dense, while fixed-CM Frobenius conditions are too periodic. The Apéry problem needs a genuinely moving-characteristic, moving-index correlation theorem.

Repository state audited: current `main` head
[`3b32157484571a9d1997bf89c703c810e9843609`](https://github.com/xiangyazi24/Ramanujan_Challenge/commit/3b32157484571a9d1997bf89c703c810e9843609), including the central-binomial warning in [`drops/Q8763-b449ee8b.md`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/main/drops/Q8763-b449ee8b.md). I use the Q8740/LGTQ-3 formulation stated in the request because the literal symbol `LGTQ-3` is not checked into that public head.

## 1. The Legendre coefficient zero set is exactly an upper-half interval

The holomorphic Legendre period is

\[
{}_2F_1\!\left(\frac12,\frac12;1;t\right)
 =\sum_{d\ge0}a_d t^d,
\qquad
 a_d=\frac{\binom{2d}{d}^2}{16^d}.
\]

For divisibility at an odd prime \(p\), the factor \(16^d\) is a unit, so put

\[
c_d:=\binom{2d}{d}^2.
\]

### Proposition 1 — exact Kummer zero set

For an odd prime \(p\) and \(0\le d<p\),

\[
v_p(c_d)=
\begin{cases}
0,&0\le d\le (p-1)/2,\\
2,&(p+1)/2\le d\le p-1.
\end{cases}
\tag{1.1}
\]

Consequently

\[
\boxed{
Z_p:=\{0\le d<p:p\mid c_d\}
=\{(p+1)/2,\ldots,p-1\},
}
\tag{1.2}
\]

and \(|Z_p|=(p-1)/2\).

### Proof

Kummer’s theorem says that

\[
v_p\binom{2d}{d}
\]

is the number of carries when adding \(d+d\) in base \(p\). Since \(d<p\), there is no carry precisely when \(2d<p\), and exactly one carry when \(2d\ge p\). Squaring doubles the valuation. \(\square\)

This is already qualitatively opposite to the Apéry situation: the local zero density is \(1/2\), and the normalized first quotient \(c_d/p\) still vanishes modulo \(p\).

## 2. These are not coefficient zeros of the actual Legendre Hasse polynomial

Let \(p=2m+1\). For

\[
E_t:y^2=x(x-1)(x-t),
\]

the Hasse invariant is, up to the irrelevant global sign,

\[
H_p(t)=(-1)^m\sum_{d=0}^{m}\binom{m}{d}^2t^d.
\tag{2.1}
\]

Using

\[
\binom{(p-1)/2}{d}
\equiv(-1)^d4^{-d}\binom{2d}{d}\pmod p,
\]

one gets

\[
H_p(t)
\equiv(-1)^m
\sum_{d=0}^{m}rac{\binom{2d}{d}^2}{16^d}t^d
\pmod p.
\tag{2.2}
\]

Every coefficient in the natural range \(0\le d\le m\) is nonzero modulo \(p\). The Kummer zeros in (1.2) occur only in the **discarded upper half** \(m<d<p\), beyond the degree of the Hasse truncation.

Thus there are three different notions:

- the Picard–Fuchs coefficient sequence \(a_d\);
- the degree-\((p-1)/2\) Hasse polynomial \(H_p(t)\);
- the roots in the parameter \(t\) of \(H_p(t)\), which are supersingular fibers.

The proposed set \(Z_p\) uses the first notion outside the Hasse truncation. It is not the supersingular locus of the Legendre family. Adolphson–Sperber explicitly record that the Legendre Hasse invariant is the truncation of \({}_2F_1(1/2,1/2;1;t)\) at degree \((p-1)/2\): [arXiv:1209.2448](https://arxiv.org/abs/1209.2448).

## 3. Exact fixed-quotient hit interval

Fix \(q\ge1\), and for a prime in the quotient cell write

\[
\frac{n}{q+1}<p\le\frac nq,
\qquad
 d=n-qp,
\qquad 0\le d<p.
\]

Define the Legendre/Kummer hit set

\[
\mathcal P_q^{\mathrm{Leg}}(n)
:=\left\{p:\frac{n}{q+1}<p\le\frac nq,
\ p\mid c_{n-qp}\right\}.
\]

By (1.2), the hit condition is

\[
2(n-qp)\ge p+1.
\]

Therefore

\[
\boxed{
\mathcal P_q^{\mathrm{Leg}}(n)
=\left\{p\text{ prime}:
\frac{n}{q+1}<p\le\frac{2n-1}{2q+1}ight\}.
}
\tag{3.1}
\]

The interval length is exactly

\[
\frac{2n-1}{2q+1}-\frac{n}{q+1}
=\frac{n-q-1}{(q+1)(2q+1)}.
\tag{3.2}
\]

For fixed \(q\), the prime number theorem gives

\[
\#\mathcal P_q^{\mathrm{Leg}}(n)
\sim
\frac{n}{(q+1)(2q+1)\log n},
\tag{3.3}
\]

and

\[
\sum_{p\in\mathcal P_q^{\mathrm{Leg}}(n)}\log p
\sim
\frac{n}{(q+1)(2q+1)}.
\tag{3.4}
\]

So even the fixed-cell logarithmic radical is linear in \(n\), not \(o(n)\). The Kummer model realizes the maximal pointwise alignment that the Apéry program is trying to exclude.

## 4. Unconditional short-gap triple saturation

The Q8740 obstruction is not merely compatible with the Legendre model; it is forced by it.

List the primes in (3.1) as

\[
r_1<r_2<\cdots<r_M
\]

and let \(L\) be the length in (3.2). Then

\[
\sum_{i=1}^{M-2}(r_{i+2}-r_i)
=(r_{M-1}+r_M)-(r_1+r_2)
\le2L.
\tag{4.1}
\]

Since \(M\sim L/\log n\), at least \(M/2+O(1)\) of these consecutive triples satisfy

\[
r_{i+2}-r_i\le C_q\log n
\tag{4.2}
\]

for one fixed constant \(C_q\) and all large \(n\). Every such triple consists of three Kummer hits in the same quotient cell.

Write its two prime gaps as

\[
h_1=r_{i+1}-r_i,
\qquad
h_2=r_{i+2}-r_{i+1}.
\]

There are only \(O_q((\log n)^2)\) possible pairs with \(h_1+h_2\le C_q\log n\). Pigeonholing gives an exact pair \((h_1,h_2)\), depending on \(n\), for which

\[
\boxed{
\#\left\{p:
 p,p+h_1,p+h_1+h_2\in\mathcal P_q^{\mathrm{Leg}}(n)
\right\}
\gg_q\frac{n}{(\log n)^3}.
}
\tag{4.3}
\]

This proof is unconditional and uses only Kummer plus the prime number theorem.

Hence an LGTQ-3 analogue asserting

\[
\max_{h_1,h_2\ll\log n}
T_{q;h_1,h_2}(n)
=o_q\!\left(\frac{n}{(\log n)^3}\right)
\]

is false for this model. It reaches the exact critical order found in Q8740.

For a **fixed** admissible gap pair, a classical upper-bound sieve gives

\[
T_{q;h_1,h_2}^{\mathrm{Leg}}(n)
\ll_{q,h_1,h_2}
\frac{n}{(\log n)^3}.
\tag{4.4}
\]

The Hardy–Littlewood prime-tuple conjecture predicts a matching positive-constant asymptotic when the three linear forms are admissible. Thus even conjecturally there is no extra little-\(o\) saving.

## 5. If “fully transverse” includes cross-unit conditions, the model is vacuous

Suppose a strict transverse triple requires not only

\[
p_i\mid c_{d_i},
\]

but also

\[
p_i\nmid c_{d_j}\qquad(i\ne j).
\]

The Kummer model has the opposite behavior. Let the three primes have mutual gaps at most \(H\), and put

\[
d_j=n-qp_j.
\]

Trim the interval (3.1) by requiring

\[
(q+1)p_j-n>H,
\qquad
2n-(2q+1)p_j>H.
\tag{5.1}
\]

This removes only endpoint intervals of length \(O_q(H)\). For any \(|p_i-p_j|\le H\),

\[
p_i-d_j
=(q+1)p_j-n+(p_i-p_j)>0,
\tag{5.2}
\]

and

\[
2d_j-p_i
=2n-(2q+1)p_j+(p_j-p_i)>0.
\tag{5.3}
\]

Thus

\[
\frac{p_i+1}{2}\le d_j<p_i,
\]

so Kummer gives

\[
\boxed{p_i^2\mid c_{d_j}\quad\text{for every }i,j.}
\tag{5.4}
\]

The local divisibility matrix is generically all ones. Therefore:

- under the own-hit definition, the model has maximal triple spikes;
- under the cross-unit definition, the transverse count is essentially zero for the trivial reason that all nearby conditions coincide.

Neither interpretation models the Apéry phenomenon, where distinct-prime hits can remain genuinely transverse.

## 6. The Legendre family is generically non-CM

The Legendre family has

\[
j(t)=256\frac{(1-t+t^2)^3}{t^2(1-t)^2},
\tag{6.1}
\]

which is nonconstant. Hence its generic geometric endomorphism ring is \(\mathbb Z\). CM occurs only at special algebraic parameters for which \(j(t)\) is a singular modulus.

More generally, over characteristic zero a one-parameter elliptic family whose generic fiber has CM by a fixed imaginary quadratic order has constant CM \(j\)-invariant and is isotrivial after finite base change. Thus “CM elliptic curve family” cannot mean a nonisotrivial analogue of the Legendre variation with every fiber CM.

Specializing the Legendre family at a CM parameter changes the Frobenius traces of that **fixed fiber**. It does not change the global Taylor coefficients

\[
\binom{2d}{d}^2/16^d
\]

of the period at the maximally unipotent point into Hecke-character values.

## 7. What a genuine fixed-CM analogue proves

Let \(E/\mathbb Q\) have CM by an order in an imaginary quadratic field \(K\). Away from finitely many bad or ramified primes, Deuring’s criterion says

\[
a_p(E)=0
\quad\Longleftrightarrow\quad
p\text{ is inert in }K.
\tag{7.1}
\]

Equivalently, the trace-zero condition is a quadratic splitting condition

\[
\chi_K(p)=-1.
\]

For fixed gaps define

\[
T_{E;h_1,h_2}(X)
=\#\left\{X<p\le2X:
\begin{array}{l}
 p,p+h_1,p+h_1+h_2\text{ are prime},\\
 a_p(E)=a_{p+h_1}(E)=a_{p+h_1+h_2}(E)=0
\end{array}
\right\}.
\tag{7.2}
\]

The CM restrictions are a finite union of residue classes modulo the conductor/discriminant. Applying the Selberg upper-bound sieve in those progressions gives

\[
\boxed{
T_{E;h_1,h_2}(X)
\ll_{E,h_1,h_2}
\frac{X}{(\log X)^3}
}
\tag{7.3}
\]

for every fixed admissible pair.

This is the available **YES** statement. It is only the natural sieve-order bound, not the needed saving.

### CM gives correlation, not independence

Take \(K=\mathbb Q(i)\). For a standard CM curve, the good inert primes are

\[
p\equiv3\pmod4.
\]

If \(h_1\equiv h_2\equiv0\pmod4\), then whenever \(p\equiv3\pmod4\), all three shifted primes have the same inert residue class. The three CM conditions are perfectly correlated. If the shifts are incompatible modulo \(4\), the count is identically zero.

Hecke characters identify the local factor and the splitting law, but they do not make the primality of

\[
p,\quad p+h_1,\quad p+h_1+h_2
\]

independent. The Hardy–Littlewood heuristic again predicts

\[
T_{E;h_1,h_2}(X)
\sim C_{E,h_1,h_2}\frac{X}{(\log X)^3}
\]

with \(C_{E,h_1,h_2}>0\) for admissible patterns. The Hecke character modifies \(C\); it does not improve the exponent.

A modern statement of the relevant CM criterion is recorded, for example, in discussions of Deuring’s criterion for CM curves: nonsplit good primes are supersingular and have \(a_p(E)=0\). See also the original CM/Deuring theory and standard CM elliptic-curve texts.

## 8. What actually breaks in the Apéry case

The distinction is sharper than “CM versus non-CM.”

### Fixed CM curve

For a fixed CM curve:

- the motive is fixed;
- the Galois representation is abelian after restriction to the CM field;
- the trace-zero event is a fixed finite-conductor splitting condition;
- Chebotarev/Hecke theory controls its one-prime density;
- a classical sieve handles prime tuples at the natural \((\log X)^{-3}\) scale.

### Apéry coefficient hit

For Apéry:

\[
p\mid b_d,
\qquad d=n-qp,
\]

is not the trace-zero condition of one fixed CM fiber. It is a coefficient condition in the defining characteristic, with both the characteristic and the coefficient index moving together. There is no known fixed number field, finite conductor, or Hecke character \(\chi\) for which

\[
p\mid b_{n-qp}
\quad\Longleftrightarrow\quad
\chi(p)\in C
\]

for a fixed set \(C\).

The Q8740 saturation result says that the rank-two recurrence module already accounts for every scalar short-gap relation; the normalized quotient data remain free. Q8746 then shows that sparse abstract zero sets plus all current reflection/gap/sieve axioms can still spike. A successful theorem must therefore couple the **moving Frobenius/quotient-jet data across distinct residue characteristics**.

Even CM theory would not automatically solve the final short-prime-tuple step: it gives at best the critical sieve order. Hence the absence of CM is not the sole obstruction. The required input is stronger:

\[
\boxed{
\text{an Apéry-specific cross-characteristic law producing a proper locus}\
\text{and hence a saving beyond the ordinary prime-triple sieve dimension.}
}
\]

## 9. Final answer

- **Legendre coefficient/Kummer analogue:** rigorously analyzable, but it fails LGTQ-3 maximally under the own-hit interpretation and is vacuous under strict cross-transversality.
- **Actual Legendre Hasse polynomial:** its natural coefficients are all nonzero; supersingular fibers are roots in the parameter, not coefficient zeros.
- **Fixed CM trace analogue:** an \(O(X/(\log X)^3)\) upper bound is provable by a classical sieve, but a little-\(o\) saving is neither supplied by Hecke characters nor expected for admissible gaps.
- **Apéry comparison:** the missing ingredient is not merely non-CM monodromy. It is a moving-index, moving-characteristic quotient/Frobenius correlation theorem with genuine extra codimension.

So the proposed Legendre/CM example is a valuable **negative control**, not a positive proof of concept for LGTQ-3.

## References

1. E. E. Kummer, classical carry theorem for valuations of binomial coefficients.
2. A. Adolphson and S. Sperber, “Hasse invariants and mod \(p\) solutions of \(A\)-hypergeometric systems,” *J. Number Theory* 142 (2014), 183–210; [arXiv:1209.2448](https://arxiv.org/abs/1209.2448).
3. M. Deuring, foundational work on reduction of CM elliptic curves and the split/inert supersingularity criterion.
4. H. Halberstam and H.-E. Richert, *Sieve Methods*, for the upper-bound sieve for fixed prime tuples and residue-class restrictions.
5. Repository benchmark: [`drops/Q8763-b449ee8b.md`](https://github.com/xiangyazi24/Ramanujan_Challenge/blob/main/drops/Q8763-b449ee8b.md).
