# Q7703: fixed-exponent Mellin attack on the pointwise rainbow

## Verdict

**[THEOREM]** The fixed-exponent Mellin condition can be reduced exactly, in every quotient cell, to two affine prime-divisor rays for the *fixed characteristic-zero Apéry sequence*.  If

\[
n=qp+r,\qquad 1\le r\le p-2,
\]

and

\[
s=\min(r,p-1-r),
\]

then

\[
\boxed{
\mathfrak p_{p,n-q}\mid M_p(n-q)
\iff p\mid b_r
\iff p\mid b_s.}
\]

Moreover exactly one of

\[
\boxed{p=\frac{n-s}{q}}
\qquad\text{or}\qquad
\boxed{p=\frac{n+s+1}{q+1}}
\]

holds, according as `2r <= p-1` or `2r > p-1`.

So after Gross--Koblitz/Morita reduction, the moving cyclotomic coefficient field disappears from the *residual decision*.  The horizontal problem becomes a prime-divisor correlation for `b_s` along two fixed affine rays.

**[NO-GO — existing Mellin equidistribution]** The newest horizontal Mellin theorem I found, Bah--Shuddhodan, arXiv:2608.00304v2 (7 Aug 2026), really does allow the finite-field characteristic and the perverse sheaf to vary together under bounded complexity and common tannakian group.  This repairs the old “varying characteristic is unavailable” objection.  But it still averages over essentially **all multiplicative characters in each field** and produces complex/unitary Frobenius equidistribution.  It does not control the single deterministic section

\[
\chi_p=\omega_p^{\,n-q}
\]

as `p` varies, and it does not estimate divisibility at the defining prime `p`.  A single selected character has relative mass `1/(p-1)` and is invisible to any all-character equidistribution statement.

**[NO-GO — standard large sieve]** The compatible-system large sieve varies auxiliary primes `ell != p` after fixing the base characteristic.  The target event reduces the Mellin trace at the **same prime `p` that is the base characteristic**.  Setting the auxiliary `ell=p` is not an admissible specialization of the etale large sieve.  Crystalline companions exist at each fixed `p`, but do not supply an independence/sieve mechanism across the diagonal `ell=p` as `p` varies.

**[NO-GO — p-adic slope]** Gross--Koblitz gives a completely exact valuation decomposition, but every interior case has at least two slope-zero terms.  The defining-prime zero is precisely cancellation among these units.  The earlier Q7699 run exhibits genuine zero cancellations with 37, 65, and 101 minimum-slope terms.

**[CONDITIONAL — smallest closing estimate]** The all-index theorem would follow from any uniform vanishing-density estimate in the surviving quotient cells.  A particularly clean form is

\[
W_{n,q}:=\sum_{\substack{p\in I_q(n)\\p\mid b_{n-qp}}}\log p
\le \varepsilon(n)\Theta_{n,q},
\qquad \varepsilon(n)\to0,
\]

uniformly for `q < Q_n=f(n) log n`, where

\[
\Theta_{n,q}=\sum_{p\in I_q(n)}\log p.
\]

Because the quotient cells are disjoint,

\[
\sum_{q<Q_n}\Theta_{n,q}\le\vartheta(n)=O(n),
\]

so this gives `o(n)` immediately.

The stronger power-saving count requested in the question,

\[
\#B_{n,q}\ll (n/q)^{1-\delta},
\]

would give

\[
\sum_{q\le Q_n}W_{n,q}
\ll
n^{1-\delta}Q_n^\delta\log n,
\]

and therefore closes the theorem for `Q_n=n^{o(1)}`.  Controlling `q=1` alone does not close the all-index statement.

---

## 1. Geometric object and fixed-exponent Mellin trace

The repository fixes

\[
\Lambda(x,y,z)=
\frac{(1+x)(1+y)(1+z)((1+y)(1+z)+xyz)}{xyz}
\]

with

\[
\operatorname{CT}\Lambda^r=b_r.
\]

For the primitive K3 fiber trace `theta_p(a)`, `toric_fiber_k3.tex` proves the exact point-count collapse

\[
\mu_p(a)=p^2-6p+12+\theta_p(a),
\]

with `theta_p(a)` the primitive weight-two rank-three trace away from the explicitly described specializations.

The symmetric-square model in `toric_fiber_sym2.tex` gives the quadratic level-six cover

\[
8x^2+(a-1)x+a=0
\]

and the elliptic curve

\[
E_x:\quad Y^2+(1-2x)XY+x^2Y=X^3,
\]

with a rational point of exact order six.  Over the split cover the generic K3 trace is a symmetric-square trace; in the inert case it is the corresponding Asai trace.  The special corrections at `a=1` and at the two branch points are explicit multiples of `p`.  Thus for the *mod-p Mellin condition* one may replace the pointwise deflated `theta_p` by the trace of the fixed rank-three middle-extension system without changing the defining-prime zero event.

Fix `n` and a quotient `q`.  Put

\[
m=n-q.
\]

For every prime in the cell `n=qp+r`,

\[
m=q(p-1)+r,
\]

so

\[
m\equiv r\pmod{p-1}.
\]

The exact fixed-quotient Mellin theorem already recorded in the repo gives

\[
\boxed{M_p(m)\equiv-b_r\pmod{\mathfrak p_{p,m}}}
\]

for every interior `1<=r<=p-2`.  The two endpoints `r=0,p-1` are harmless because `b_0=1` and `b_{p-1}=1 mod p`.

### Uniform geometric complexity

The fixed K3/Picard--Fuchs system has rank three and four fixed singular points on the projective parameter line: `0`, `infinity`, and the two roots of `D(a)=a^2-34a+1`.  It is tame in the present range; the Kummer twist is rank one and tame, with Swan conductor zero independent of the order of the character.

After shifting the rank-three local system to a perverse sheaf on `G_m`, Grothendieck--Ogg--Shafarevich gives generic Mellin cohomological dimension

\[
3(4-2)=6
\]

when `H_c^0,H_c^2` vanish.  Resonant characters only change this by an absolute amount.  Consequently the Mellin sums have an absolute Deligne bound of the shape

\[
|M_p(\chi)|\le C p^{3/2}
\]

for an absolute `C` depending only on the fixed geometric system.

This verifies the **bounded-complexity/purity/tameness** part of the horizontal Mellin hypotheses.  Identifying the exact convolution-tannakian group for this perverse object has not been carried out in the repo; the underlying rank-three local system is the symmetric-square of the level-six elliptic system and has connected geometric local-system monodromy of `SO_3` type.  This missing tannakian identification is secondary below, because even granting the strongest common-group hypothesis does not select the required character or the defining-prime congruence.

---

## 2. Literature audit with the actual quantifiers

### 2.1 Bah--Shuddhodan 2026: varying characteristic is now allowed, but the character is averaged

Bah--Shuddhodan, **Uniform stratified vanishing and equidistribution on `G_m^d`**, arXiv:2608.00304v2, prove uniform point-count bounds for Mellin exceptional loci with constants depending only on dimension and Sawin complexity.  Their Corollary 1.4 applies to sequences

\[
(k_j,M_{j,0})
\]

with `#k_j -> infinity`, bounded complexity, fixed tannakian dimension, equality of arithmetic and geometric tannakian groups, and a common reductive tannakian group.  The characteristics may vary.

This is materially stronger than the horizontal discussion in Forey--Fresan--Kowalski, and it removes the objection that a Mellin theorem must keep the characteristic fixed.

But the output is equidistribution of

\[
\Theta_{M_{j,0},k_j}(\chi),
\qquad \chi\in X_j,
\]

where `X_j` contains all but a bounded-complexity exceptional set of characters.  The empirical measure is formed by averaging over those characters **inside each field**.

Our target is the one character

\[
\chi_{p,m}=\omega_p^m,
\qquad m=n-q\ \text{fixed in the q-cell}.
\]

The order is

\[
d_p=\frac{p-1}{\gcd(p-1,m)},
\]

which generally grows with `p`.  A sequence containing one chosen character per field has mass `1/(p-1)` inside the all-character empirical measure.  Therefore even perfect all-character equidistribution allows the chosen character to be exceptional for every `p` without changing the limit.

There is a second mismatch: Bah--Shuddhodan control complex/unitary Frobenius conjugacy classes.  The event here is

\[
\mathfrak p_{p,m}\mid M_p(m),
\]

which need not mean `M_p(m)=0` in characteristic zero.  The repo already contains the explicit example `(p,j)=(11,5)` where the lifted trace is nonzero (`-33`) although its reduction is zero modulo 11.

**Conclusion: [NO-GO].** The newest varying-characteristic Mellin equidistribution theorem does not imply a bound for this deterministic same-prime section.

### 2.2 Forey--Fresan--Kowalski

Forey--Fresan--Kowalski, arXiv:2109.11961, construct the convolution-tannakian framework and prove vertical/on-average Mellin equidistribution.  Their horizontal discussion requires quantitative stratified vanishing; Bah--Shuddhodan now supply that missing uniformity for tori.  The same selected-character/defining-prime obstruction remains.

### 2.3 Perret-Gentil large sieve: the auxiliary-prime diagonal is wrong

Perret-Gentil, **Exponential Sums Over Finite Fields and the Large Sieve**, IMRN 2020, works with a coherent family

\[
(\mathcal F_\lambda)_{\lambda\in\Lambda}
\]

over a fixed finite field of characteristic `p`.  The compatible system is reduced modulo auxiliary valuations `lambda | ell` with

\[
\ell\ne p.
\]

His Kloosterman example states this explicitly: the sieve set consists of `ell`-adic valuations with `p != ell`, and Theorem 2.7 then combines the reductions at many such `lambda`.

Our condition is reduction at

\[
\mathfrak p_{p,m}\mid p
\]

itself.  This is not one of the admissible auxiliary reductions of an etale sheaf in characteristic `p`.  Letting the base prime vary does not repair the logical problem: for each new base characteristic the forbidden auxiliary prime changes to that same characteristic.

Also, the most explicit monodromy hypotheses in Perret-Gentil are `SL` or `Sp`; the symmetric-square K3 system is naturally orthogonal.  Orthogonal variants may be obtainable by other methods, but they do not fix the diagonal `ell=p` obstruction.

**Conclusion: [NO-GO].** A standard compatible-system large sieve cannot be specialized to the defining-prime event by setting `ell=p`.

### 2.4 Crystalline companions do not supply the missing sieve

Kedlaya's **Etale and crystalline companions, II**, arXiv:2008.13053, proves that an algebraic `ell`-adic Weil sheaf on a smooth scheme over a finite field of characteristic `p` admits a crystalline `p`-adic companion with matching Frobenius characteristic polynomials.

This is exactly the right statement for moving from the fixed-characteristic etale system to a `p`-adic coefficient object at that same `p`.  It does **not** couple the companions for different rational primes, nor does it produce independent finite residual monodromy quotients across `p`.  Hence it does not turn the diagonal defining-prime condition into a horizontal large sieve.

---

## 3. Exact second moments: what they say and what they cannot say

For a fixed prime write

\[
M_p(\chi)=\sum_{a\in\mathbf F_p^\times}\Theta_p(a)\chi(a),
\]

where `Theta_p` is the clean rank-three trace, which is congruent to the manuscript's primitive trace for the defining-prime decision.

Multiplicative character orthogonality gives the exact Parseval identity

\[
\boxed{
\sum_{\chi}M_p(\chi)M_p(\chi^{-1})
=(p-1)\sum_{a\ne0}\Theta_p(a)^2.}
\tag{3.1}
\]

Under a complex embedding and unitary normalization,

\[
\boxed{
\sum_{\chi}|M_p(\chi)|^2
=(p-1)\sum_{a\ne0}|\Theta_p(a)|^2.}
\tag{3.2}
\]

Summing (3.2) over primes in one quotient cell gives a perfectly exact bilinear/second-moment identity across the cell.  Deligne gives the expected total size `O(sum p^4)` before normalization.

But (3.2) averages over all `p-1` characters.  It cannot upper-bound the number of primes for which the one selected `chi_{p,m}` has a prescribed residual property.  A deterministic one-point section can be arbitrary while contributing `O(1/p)` to the character average.

Even a hypothetical power-saving complex second moment for the selected traces

\[
\sum_{p\in I_q}|M_p(\chi_{p,m})|^2
\]

would not detect defining-prime divisibility: a nonzero algebraic integer can be divisible by the chosen prime above `p` while having generic complex size.

Thus the strongest standard second moment has the wrong observable.

---

## 4. Exact fixed-exponent arithmetic reduction

This is the main theorem produced by the audit.

### Theorem 4.1 (fixed-n prefix reduction)

Fix `n,q,p` with

\[
n=qp+r,
\qquad 1\le r\le p-2,
\]

and put

\[
s=\min(r,p-1-r).
\]

Define the fixed characteristic-zero summands

\[
A_n(k)=\left(\binom nk\binom{n+k}{k}\right)^2
\]

and prefix

\[
S_n(s)=\sum_{k=0}^s A_n(k).
\]

Then

\[
\boxed{b_r\equiv S_n(s)\equiv b_s\pmod p.}
\tag{4.1}
\]

#### Proof

For `k<p`, Lucas gives

\[
\binom nk\equiv\binom rk\pmod p.
\]

For `k<=s`, one has `r+k<=p-1`, so there is no carry and

\[
\binom{n+k}{k}\equiv\binom{r+k}{k}\pmod p.
\]

Therefore each `A_n(k)` agrees with the corresponding Apéry summand for `b_r` through `k=s`.

If `r<=(p-1)/2`, then `s=r`, so there is no omitted tail.

If `r>(p-1)/2`, then `s=p-1-r`.  For every `s<k<=r`, the factorial quotient `binom(r+k,k)` contains exactly one factor `p`; after squaring, every omitted Apéry summand is divisible by `p^2`.  This proves the first congruence in (4.1).

For the second congruence, on the left branch `s=r` trivially.  On the right branch, the candidate prime relation below gives

\[
n\equiv-s-1\pmod p.
\]

Then

\[
\binom{-s-1}{k}=(-1)^k\binom{s+k}{k},
\qquad
\binom{-s-1+k}{k}=(-1)^k\binom sk,
\]

so the squared product is exactly the `k`th Apéry summand for `b_s` modulo `p`.

Equivalently, this is the usual Apéry reflection

\[
b_r\equiv b_{p-1-r}\pmod p.
\]

QED.

### Theorem 4.2 (two affine prime-divisor rays)

With the same notation, assign the center `2r=p-1` to the left branch.

- If `2r<=p-1`, then `s=r` and
  \[
  \boxed{p=(n-s)/q,\qquad p\mid b_s.}
  \tag{4.2L}
  \]
  The allowed folded indices satisfy
  \[
  (2q+1)s\le n-q.
  \]

- If `2r>p-1`, then `s=p-1-r` and
  \[
  \boxed{p=(n+s+1)/(q+1),\qquad p\mid b_s.}
  \tag{4.2R}
  \]
  The allowed folded indices satisfy
  \[
  (2q+1)s\le n-2q-1.
  \]

Conversely, an integer `s>=1` satisfying the corresponding integrality, primality, interior, and branch inequalities reconstructs the original prime in the cell.

There is also a purely integer divisibility form:

\[
\boxed{n-s\mid q b_s}
\tag{4.3L}
\]

on the left, and

\[
\boxed{n+s+1\mid(q+1)b_s}
\tag{4.3R}
\]

on the right.

Thus a complete quotient-cell count is exactly a prime-linear-form/divisibility correlation for one fixed recurrence sequence.

### Fixed-prefix recurrence

The prefix introduced above is itself holonomic in `s` for fixed `n`.  Since

\[
\frac{A_n(s)}{A_n(s-1)}
=\left(\frac{(n-s+1)(n+s)}{s^2}\right)^2,
\]

putting `S_s=S_n(s)` gives

\[
\boxed{
s^4(S_s-S_{s-1})
=((n-s+1)(n+s))^2(S_{s-1}-S_{s-2}).}
\tag{4.4}
\]

This is an exact fixed-exponent recurrence, but because (4.1) folds it back to `b_s` modulo the candidate prime, it does not by itself create new anti-concentration.

---

## 5. Exact prime-count identity in one q-cell

Let `B_{n,q}` denote the bad primes in the interior quotient cell.  The two-ray theorem gives the exact decomposition

\[
\boxed{B_{n,q}=B^-_{n,q}\sqcup B^+_{n,q}}
\]

where

\[
B^-_{n,q}=
\left\{
\frac{n-s}{q}:\
\frac{n-s}{q}\text{ prime},\
(2q+1)s\le n-q,\
\frac{n-s}{q}\mid b_s
\right\},
\tag{5.1L}
\]

with the obvious congruence `s=n mod q`, and

\[
B^+_{n,q}=
\left\{
\frac{n+s+1}{q+1}:\
\frac{n+s+1}{q+1}\text{ prime},\
(2q+1)s\le n-2q-1,\
\frac{n+s+1}{q+1}\mid b_s
\right\},
\tag{5.1R}
\]

with `s=-n-1 mod(q+1)`.

Equivalently, if `1_P(x)` is the prime indicator,

\[
\boxed{
\begin{aligned}
\#B_{n,q}
={}&\sum_{\substack{s\ge1\\(2q+1)s\le n-q}}
1_{\mathbb P}\!\left(\frac{n-s}{q}\right)
1_{n-s\mid qb_s}\\
&+\sum_{\substack{s\ge1\\(2q+1)s\le n-2q-1}}
1_{\mathbb P}\!\left(\frac{n+s+1}{q+1}\right)
1_{n+s+1\mid(q+1)b_s}.
\end{aligned}}
\tag{5.2}
\]

This is the strongest exact horizontal identity I obtain.  A power-saving for (5.2) is already the desired theorem; it is not a restatement involving a complex Deligne bound.

One may also insert additive-character orthogonality for each candidate prime,

\[
1_{p\mid b_s}=\frac1p\sum_{h\bmod p}e^{2\pi i h b_s/p},
\tag{5.3}
\]

but this introduces a different modulus and a different single sample for each `p`.  The usual analytic large sieve has no common coefficient vector to exploit.  Equation (5.3) is exact but circular at the level of available estimates.

---

## 6. Gross--Koblitz/Morita: exact slopes, no fixed p-adic function of p

For

\[
T_{r,k}=\left(\binom rk\binom{r+k}{k}\right)^2,
\]

Q7699 proves the exact identity

\[
T_{r,k}=p^{2\mathbf 1_{r+k\ge p}}U_{p,r,k}^2,
\]

where

\[
U_{p,r,k}
=\frac{\Gamma_p(r+k+1)}
{\Gamma_p(k+1)^2\Gamma_p(r-k+1)}
\in\mathbf Z_p^\times.
\]

With

\[
s=\min(r,p-1-r),
\]

exactly `s+1` terms have valuation zero and every remaining term has valuation two.  Hence

\[
\boxed{
p\mid b_r
\iff
\sum_{k=0}^s\overline U_{p,r,k}^{\,2}=0\quad\text{in }\mathbf F_p.}
\tag{6.1}
\]

There are always at least two slope-zero terms in the interior range.  There is never a unique-minimum Newton-polygon exclusion.

For the fixed `n=321` top-half experiment in Q7699, the genuine bad primes

\[
179,\quad193,\quad211
\]

have respectively

\[
37,\quad65,\quad101
\]

minimum-slope unit contributions canceling to zero.  This verifies that the cancellation obstruction occurs at actual defining-prime zeros, not merely at generic test points.

Baldassarri--Cailotto's Dwork-family/Boyarsky theory gives `p`-adic analytic dependence of Frobenius on **hypergeometric exponents for a fixed p-adic setting**.  Here the rational prime itself changes, the fields `Q_p` change, the Morita gamma functions `Gamma_p` change, and the unit packet length `s+1` grows with the prime scale.  Thus this theory does not produce one fixed `p`-adic analytic function `F(p)` whose zeros are the horizontal bad primes.

**[NO-GO].** The p-adic reduction is exact, but it lands on the same folded Apéry divisibility problem (4.2), not a bounded-complexity analytic zero set.

---

## 7. Coefficient fields and the norm/height obstruction

For fixed `m=n-q`, the Teichmuller character has order

\[
d_p=\frac{p-1}{\gcd(p-1,m)}.
\]

The unnormalized algebraic Mellin trace lies in

\[
K_p=\mathbf Q(\mu_{d_p}).
\]

Since `d_p | p-1`, one has `p = 1 mod d_p`; hence `p` splits completely in `K_p`.  The chosen Teichmuller embedding selects one prime

\[
\mathfrak p_{p,m}\mid p
\]

with residue field `F_p`.

The target is divisibility at this **one selected prime**, not divisibility by every prime over `p` and not characteristic-zero vanishing.

Every Galois conjugate of the Mellin trace is another bounded-conductor Kummer Mellin trace, so Deligne gives

\[
|M_p(m)^\sigma|\le C p^{3/2}.
\]

Therefore, when the norm is nonzero,

\[
\log|N_{K_p/\mathbf Q}M_p(m)|
\le
\varphi(d_p)\left(\frac32\log p+O(1)\right).
\tag{7.1}
\]

Divisibility at the selected prime contributes only one factor `p` to the norm.  The degree `phi(d_p)` is unbounded and can be of order comparable to `p`.  Consequently (7.1) gives no useful zero-density estimate.

The exact integer reduction (4.2) avoids this coefficient-field inflation, but then ordinary height is still too large for a product argument.  From the Apéry sum,

\[
0<b_s
\le(s+1)2^{6s},
\]

so

\[
\log b_s\le6s\log2+\log(s+1).
\tag{7.2}
\]

In a q-cell, `s` ranges to order `P=n/q`, though only one residue class modulo `q` or `q+1` occurs on each ray.  Summing (7.2) over all candidate `s` gives total logarithmic height of order at least the scale `P^2/q`, far larger than the `P/q` scale of the prime population.  Dividing by `log P` therefore does not even improve the trivial prime count.

**[NO-GO].** Neither the cyclotomic norm nor the fixed-integer prefix/product has sufficiently small height to force a power-saving count.

---

## 8. Precise large-sieve quantifier no-go

It is useful to isolate the three independent quantifiers.

1. **Base characteristic:** `p` varies over primes in a quotient cell.
2. **Mellin character:** `chi_{p,m}=omega_p^m` varies with `p`; fixed exponent `m` does *not* mean fixed character order.
3. **Residual prime:** the event is reduction modulo the same rational prime `p` (more precisely the selected `mathfrak p_{p,m}|p`).

The newest horizontal Mellin equidistribution handles (1), and lets the collection of characters vary, but it averages over all characters and therefore does not handle the deterministic diagonal in (2).  Standard compatible-system large sieves handle many residual primes, but only auxiliary primes away from the base characteristic, so they do not handle (3).

A fixed exponent also does not define a fixed finite-order Kummer cover: the order is `d_p=(p-1)/gcd(p-1,m)`, usually growing with `p`.

Therefore there is no valid chain

```text
bounded conductor
=> existing horizontal large sieve
=> selected defining-prime zero-density.
```

The first implication is the exact place where the quantifiers change.

---

## 9. Smallest genuinely new estimate that closes the proof

Let

\[
I_q(n)=\left(\frac n{q+1},\frac nq\right]\cap\mathbb P
\]

and

\[
B_{n,q}=
\{p\in I_q(n):1\le n-qp\le p-2,\ p\mid b_{n-qp}\}.
\]

Let

\[
W_{n,q}=\sum_{p\in B_{n,q}}\log p,
\qquad
\Theta_{n,q}=\sum_{p\in I_q(n)}\log p.
\]

### Conditional estimate A: weakest clean weighted form

**[CONDITIONAL — PMLS-w]** There is `epsilon(n)->0` such that, uniformly for

\[
1\le q<Q_n=f(n)\log n,
\]

one has

\[
\boxed{W_{n,q}\le\epsilon(n)\Theta_{n,q}.}
\tag{9.1}
\]

Because the cells are disjoint,

\[
\sum_{q<Q_n}\Theta_{n,q}\le\vartheta(n)=O(n),
\]

and hence

\[
\sum_{q<Q_n}W_{n,q}=o(n).
\]

Together with the manuscript's quotient reduction, this closes the all-index theorem.

This is weaker than a power-saving count and is the smallest natural “new estimate” suggested by the present audit: **same-prime anti-concentration for the two affine divisor rays (5.1)**.

### Conditional estimate B: requested power saving

If for some fixed `delta>0`

\[
\boxed{\#B_{n,q}\ll P_q^{1-\delta},\qquad P_q=n/q}
\tag{9.2}
\]

uniformly for `q<=Q_n`, then

\[
W_{n,q}\ll(n/q)^{1-\delta}\log n.
\]

Therefore, for `0<delta<1`,

\[
\begin{aligned}
\sum_{q\le Q_n}W_{n,q}
&\ll n^{1-\delta}\log n
\sum_{q\le Q_n}q^{\delta-1}\\
&\ll_\delta
\boxed{n^{1-\delta}Q_n^\delta\log n}.
\end{aligned}
\tag{9.3}
\]

For `Q_n=n^{o(1)}` this is `o(n)`.  This is exactly the summation required by the prompt.

The q=1 cell contributes only the top-half primes.  It is one summand in (9.3), not a substitute for the full quotient range.

---

## 10. Reproducible exact experiment

New verifier:

```text
problems/3.2/research/scripts/q7703_fixed_exponent_prefix.sage
```

It checks, with exact Sage integers, for every interior prime in the requested cells:

```text
m = n-q
m mod (p-1) = r
s = min(r,p-1-r)
b_r = S_n(s) = b_s mod p
left:  p=(n-s)/q,      n-s | q*b_s
right: p=(n+s+1)/(q+1), n+s+1 | (q+1)*b_s
```

and independently verifies the fixed-prefix recurrence (4.4).

The first audit (`GitHub Actions` run `31573057876`, Sage 10.6) checked `n=321,q<=8` and `n=2000,q<=16`.  For `n=321` it reproduced the three Q7699 bad rows in the right branch:

```text
(p,r,s) = (179,142,36), (193,128,64), (211,110,100).
```

All 273 interior rows for `n=2000,q<=16` passed the exact reduction.

A second run after adding the explicit folded-`b_s` assertions was used to recheck the sharpened form.

Example command:

```bash
sage -python problems/3.2/research/scripts/q7703_fixed_exponent_prefix.sage \
  --n 321 --qmax 8
```

The experiment verifies identities only; its observed zero counts are not used as theorem evidence.

---

## 11. Final theorem/conditional/no-go ledger

- **[THEOREM]** Exact fixed-quotient Mellin congruence for every `q`, already in the repo.
- **[THEOREM]** Uniform fixed geometric complexity under the varying tame Kummer character.
- **[THEOREM]** Exact Gross--Koblitz/Morita slope profile; the zero event is equal-slope unit cancellation.
- **[THEOREM]** New fixed-exponent prefix/folding identity `b_r = S_n(s) = b_s mod p`.
- **[THEOREM]** New two-affine-ray reformulation `(n-s)/q | b_s` or `(n+s+1)/(q+1) | b_s` after the harmless scalar factor.
- **[THEOREM]** Exact vertical Mellin Parseval second moment.
- **[NO-GO]** Bah--Shuddhodan horizontal equidistribution averages over all characters; it does not control one selected `omega_p^m` per field.
- **[NO-GO]** Complex Mellin equidistribution does not detect divisibility at the defining prime.
- **[NO-GO]** Perret-Gentil/Kowalski compatible-system large sieve uses auxiliary `ell != p`; the required diagonal `ell=p` is outside its etale setup.
- **[NO-GO]** Crystalline companions are fixed-characteristic companions, not a horizontal same-prime sieve.
- **[NO-GO]** Gross--Koblitz/Dwork supplies no fixed p-adic analytic function in the varying rational prime.
- **[NO-GO]** Cyclotomic norm degree and fixed-prefix integer height are too large for a height-only density argument.
- **[CONDITIONAL]** Uniform `o(1)` log-weight density per q-cell closes the manuscript.
- **[CONDITIONAL]** The stronger `O((n/q)^(1-delta))` cell count gives `O(n^(1-delta) Q^delta log n)` and closes for `Q=n^o(1)`.

The smallest genuinely new estimate is therefore not another Deligne bound or another all-character Mellin theorem.  It is a **diagonal defining-prime anti-concentration estimate for the two affine Apéry divisor rays (5.1)**, uniformly over the slowly growing quotient range.

---

## References checked

1. A. Bah and K. V. Shuddhodan, *Uniform stratified vanishing and equidistribution on `G_m^d`*, arXiv:2608.00304v2 (2026).
2. A. Forey, J. Fresan, E. Kowalski, *Arithmetic Fourier transforms over finite fields: generic vanishing, convolution, and equidistribution*, arXiv:2109.11961, current manuscript dated 22 Mar 2026.
3. C. Perret-Gentil, *Exponential Sums Over Finite Fields and the Large Sieve*, IMRN 2020, 7139--7174, doi:10.1093/imrn/rny202.
4. K. S. Kedlaya, *Etale and crystalline companions, II*, arXiv:2008.13053.
5. W. Sawin, A. Forey, J. Fresan, E. Kowalski, *Quantitative sheaf theory*, arXiv:2101.00635.
6. F. Baldassarri, M. Cailotto, *p-adic formulas and unit root F-subcrystals of the hypergeometric system*, arXiv:math/0409207.
7. B. Zurbuchen, *Equidistribution for Tannakian monodromy groups*, arXiv:2602.21878.

Repository inputs used include `toric_fiber_k3.tex`, `toric_fiber_sym2.tex`, `oracleC_exploration.md`, `chatgpt_q32_horizontal_mellin_sieve_corrected.md`, and `chatgpt_q32_padic_horizontal_mellin.md`.
