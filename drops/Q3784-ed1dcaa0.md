ANSWER Q3784 ed1dcaa0

# Verdict

I do **not** obtain a full unconditional closure of the far reflected-gap residual.

I do obtain two genuine unconditional improvements, the first of which is a strictly stronger deletion than the one stated in the prompt:

1. **Matched lower-legality-margin deletion.** Put
   \[
   k:=p-6d>0.
   \]
   For either sign, true rows with `k<=K` have total row count
   \[
   \boxed{\#\mathcal E_{k\le K}(T)\ll_J \frac{K^3}{\log T}}
   \]
   and total logarithmic weight
   \[
   \boxed{
   \sum_{e\in\mathcal E_{k\le K}(T)}\log p_e\ll_J K^3.
   }
   \tag{A}
   \]
   Hence every
   \[
   \boxed{
   K=o\!\left(\frac{T^{2/3}}{(\log T)^{1/3}}\right)
   }
   \tag{B}
   \]
   is unconditionally deletable at cost `o_J(T^2/log T)`.

   On the **plus** branch, `k=h_+=p-6d`, so (B) improves the previous plus reflected-gap cutoff
   \[
   h_+\le T^{3/5}/\omega(T)
   \]
   all the way to
   \[
   \boxed{
   h_+\le
   \frac{T^{2/3}}{(\log T)^{1/3}\Omega(T)}
   }
   \qquad(\Omega(T)\to\infty).
   \tag{C}
   \]

   On the **minus** branch this is a genuinely different strip. In the fixed-`J` long-gap sector `d>T/J`, if `K=o(T)` then eventually `k<d`, hence
   \[
   h_-=|p-7d|=d-k\asymp_J T.
   \]
   Thus (B) removes rows lying deep inside the old *far* residual, near the lower-legality boundary `p=6d`, not near the reflection centre `p=7d`.

2. **Actual-pair energy improvement.** Summing the genuine condition
   \[
   p\mid b_x,\qquad p\mid N_h(x)
   \]
   over all reflected gaps is much cheaper than summing the marginal root envelope. It is exactly a pair of actual Apéry zeros, so for fixed `p`
   \[
   \boxed{
   \sum_d |X^+_{p,d}|\le \binom{Z(p)}2,
   \qquad
   \sum_d |X^-_{p,d}|\le 2\binom{Z(p)}2.
   }
   \tag{D}
   \]
   Therefore `Z(p)<<p^(2/3)` gives the global weighted bound
   \[
   \boxed{
   \sum_{\text{all sign rows in the shell}}\log p
   \ll_J T^{7/3}.
   }
   \tag{E}
   \]
   This strictly improves the `T^(8/3)/log T` marginal-capacity ledger, but it is still supercritical relative to `T^2/log T` by a factor of order `T^(1/3) log T`.

After (A)--(E), the first surviving theorem is a genuinely **cross-characteristic matched zero-pair dispersion estimate**. I do not find such an estimate in the pushed canonical sources. Moreover there are exact scoped obstructions to the obvious routes:

- replacing the actual `q|b_mirror` condition by a reflected formal continuant root is invalid: that formal root is the forced centre factor and occurs for every candidate index, whether or not the Apéry coefficient vanishes;
- scalar resultants lose the common root-component label, with an exact deep counterexample already proved in `gcdtail_result.tex`;
- the current row-degree/column/restart/Bezout bounds cannot force deep-pair sparsity: `gcdtail_result.tex` contains an explicit incidence-model obstruction satisfying the local `L^(2/3)` column law and even codegree one while retaining `Theta(H^2)` deep pairs;
- classical additive large sieve is already recorded in `proof.tex` as stopping one logarithm above the required cross-prime dispersion scale.

So the result of this attack is: **strictly stronger unconditional deletion, but no full far-residual closure.** No finite scan is used as proof.

---

# 0. Source boundary

The prompt-named caller-local file

```text
problems/3.2/research/proofs/Q32_Q6_REFLECTED_GAP_LOCAL_CAPACITY_DELETION.md
```

is not present on connector-visible `main`, and no connector-visible commit contains that filename. I therefore cannot honestly claim a byte-for-byte audit of that unpushed local file.

The prompt explicitly allows the exact `X^-` and `X^+` parameterizations as input. I checked them against the immediately preceding same-project delivered derivation

```text
drops/Q3779-1c5b9d1a.md
```

and audited the canonical pushed sources it uses:

```text
problems/3.2/proof.tex
problems/3.2/FABLE_SECTION_orbit_energy.tex
problems/3.2/gcdtail_result.tex
```

The relevant pushed locations are:

- `proof.tex`, around lines 620--760: pollution classification, `prop:column`, rectangular incidence, and the beginning of the continuant Bezout section;
- `proof.tex`, around lines 3160--3245: the adversarial cross-prime no-go model;
- `proof.tex`, around lines 3240--3375: Fourier/large-sieve form of the CRT error and the AP--BDH dispersion target;
- `FABLE_SECTION_orbit_energy.tex`, restart dictionary / one-step shift and `deg N_h=3(h-1)`;
- `gcdtail_result.tex`, beginning of the file: deep fixed-involution reformulation, the affine-plane incidence obstruction, and the deep scalar-resultant component-misalignment counterexample.

I use no finite checker output as an asymptotic theorem.

---

# 1. Accepted exact local parameterizations

Write a q6 row as
\[
m=6q+u=6p+r,
\qquad q=p+d,
\qquad r=u+6d<p.
\]
Thus the lower-legality margin
\[
\boxed{k:=p-6d}
\tag{1.1}
\]
satisfies
\[
\boxed{0\le u<k.}
\tag{1.2}
\]
For a true row `u=0` is impossible because `q|b_0=1`; asymptotically this endpoint is irrelevant, but I keep the nonnegative notation when convenient.

The accepted reflected-gap parameterizations are the following.

## Minus

Set
\[
h_-=|p-7d|,
\qquad
w_-:=\min(d,p-6d)=\min(d,k).
\]
Then true minus rows at fixed `(p,d)` are in bijection with
\[
X^-_{p,d}
=
\{0\le x<w_-:
 p\mid b_x,
 p\mid N_{h_-}(x),
 q\mid b_{w_--1-x}\}.
\tag{1.3}
\]

## Plus

Set
\[
h_+:=p-6d=k,
\qquad
w_+:=\min(h_+,p-h_+).
\]
Then true plus rows at fixed `(p,d)` are in bijection with
\[
X^+_{p,d}
=
\{0\le x<w_+:
 p\mid b_x,
 p\mid N_{h_+}(x),
 q\mid b_{h_+-1-x}\}.
\tag{1.4}
\]

The pushed restart convention matters here: a positional gap `h` is governed by
\[
N_h(r)=U_{h-1}(r+1),
\qquad
\deg N_h=3(h-1),
\]
not by an unshifted `U_h`. This is explicit in `FABLE_SECTION_orbit_energy.tex` in the restart dictionary and the “one-step shift” remark.

---

# 2. New matched-margin cell: both signs collapse to the same two endpoint zeros

The key observation is that the previous local-capacity proof used the inequality `u<k`, but did not globally group by the pair `(u,k)`.

Take
\[
1\le k\le K,
\qquad
0\le u<k,
\]
and define
\[
\boxed{x:=k-1-u.}
\tag{2.1}
\]
Then
\[
 u+x=k-1.
\tag{2.2}
\]

## 2.1 Plus branch

For plus, `h_+=k`. From (1.4),
\[
u=h_+-1-x=k-1-x,
\]
so (2.1) is exactly the oriented plus coordinate. Every true plus row in this cell supplies
\[
\boxed{p\mid b_x,\qquad q\mid b_u.}
\tag{2.3}
\]
The further condition `p|N_k(x)` is genuine and can only reduce the number of rows; it is not needed for the present upper bound.

Since
\[
q=7d+k,
\qquad
p=6d+k,
\tag{2.4}
\]
fixing `(u,x)` fixes `k`, and a row is then determined by `d` (equivalently by `p` or `q`). Distinct rows in the cell therefore have pairwise distinct lower primes and pairwise distinct upper primes.

## 2.2 Minus branch in the fixed-J far sector

For minus, assume `k<=K=o(T)` and the fixed-`J` long-gap condition
\[
d>T/J.
\]
For sufficiently large `T`,
\[
k<d.
\]
Hence this is the `p<7d` orientation, and
\[
w_-=k,
\qquad
h_-=d-k.
\tag{2.5}
\]
The accepted minus parameterization gives
\[
u=w_--1-x=k-1-x,
\]
so again (2.1) is the correct oriented coordinate. Every true minus row supplies
\[
\boxed{p\mid b_x,\qquad q\mid b_u.}
\tag{2.6}
\]
Here `p|b_x` is the reflected lower zero from the exact `X^-` construction; `q|b_u` is the genuine upper selected zero.

Again (2.4) holds. Distinct rows in one `(u,x)` cell have distinct `d`, hence distinct `p` and distinct `q`.

This is the useful point: the same matched cell controls plus rows near `h_+=0` **and minus rows near the lower-legality boundary `p=6d`**.

---

# 3. Exact cell multiplicity, including cross equality `p_i=q_j`

Let the rows in one fixed `(u,x)` cell be
\[
(p_i,q_i),
\qquad i=1,\ldots,R_{u,x}.
\]
By (2.3) or (2.6),
\[
\prod_i p_i\mid\operatorname{rad}(b_x),
\qquad
\prod_i q_i\mid\operatorname{rad}(b_u).
\tag{3.1}
\]
Therefore
\[
\boxed{
\left(\prod_i p_i\right)
\left(\prod_i q_i\right)
\mid
\operatorname{rad}(b_x)\operatorname{rad}(b_u).
}
\tag{3.2}
\]

No coprimality between the two endpoint products is required.

Indeed, the `p_i` are pairwise distinct and the `q_i` are pairwise distinct. If a rational prime `lambda` occurs as
\[
\lambda=p_i=q_j,
\]
then (2.3)/(2.6) imply simultaneously
\[
\lambda\mid b_x,
\qquad
\lambda\mid b_u.
\]
The product of the two radicals on the right side of (3.2) therefore contains the two valuation units needed by the left side. Thus cross equality does not collapse occurrence multiplicity.

In the current shell every endpoint satisfies
\[
p_i,q_i>T/7,
\]
so
\[
2R_{u,x}\log(T/7)
\le
\log b_x+\log b_u.
\tag{3.3}
\]
The elementary Apéry estimate
\[
b_n\le(n+1)64^n
\]
gives an absolute constant `kappa` with
\[
\log b_n\le\kappa n
\qquad(n\ge1).
\]
Hence, using `u+x=k-1`,
\[
\boxed{
R_{u,x}
\le
\frac{\kappa(k-1)}{2\log(T/7)}.
}
\tag{3.4}
\]
This is a true occurrence bound; it uses no valuation inference from a formal carrier.

---

# 4. Summation: the new `K^3` deletion

For fixed `k`, there are exactly `k` nonnegative ordered pairs
\[
(u,x),
\qquad u+x=k-1.
\]
Summing (3.4),
\[
\begin{aligned}
\#\mathcal E_{k\le K}(T)
&\le
\frac{\kappa}{2\log(T/7)}
\sum_{k=1}^K k(k-1)\\
&\ll
\boxed{\frac{K^3}{\log T}}.
\end{aligned}
\tag{4.1}
\]
Since all lower primes in the shell are `O_J(T)`,
\[
\log p\ll_J\log T.
\]
Therefore
\[
\boxed{
\sum_{e\in\mathcal E_{k\le K}(T)}\log p_e
\ll_J K^3.
}
\tag{4.2}
\]
Consequently
\[
K^3=o(T^2/\log T)
\]
is sufficient, i.e.
\[
\boxed{
K=o\!\left(T^{2/3}(\log T)^{-1/3}\right).
}
\tag{4.3}
\]

This argument is stronger than using the local `w^(2/3)` capacity plus the two-linear-form sieve. The latter gives `T K^(5/3)/log T` in log-weight and only reaches exponent `3/5`; the matched-cell height argument pays no ambient factor `T`.

---

# 5. Consequences for the prompt's far residual

Let `Omega(T)->infinity` and put
\[
K_0(T)
:=
\frac{T^{2/3}}
     {(\log T)^{1/3}\Omega(T)}.
\tag{5.1}
\]
Then (4.2) gives
\[
\sum_{k\le K_0}\log p
\ll_J
\frac{T^2}{\log T\,\Omega(T)^3}
=o_J(T^2/\log T).
\tag{5.2}
\]

## 5.1 Plus improvement

For plus,
\[
h_+=k.
\]
So (5.2) unconditionally deletes
\[
\boxed{
h_+\le K_0(T).}
\tag{5.3}
\]
This strictly improves the old
\[
h_+\le T^{3/5}/\omega(T)
\]
strip.

The new region is genuinely nonempty. For example
\[
h_+=T^{5/8}
\]
is eventually larger than every `T^(3/5)` cutoff with a slowly growing denominator, while
\[
T^{5/8}=o\!\left(T^{2/3}(\log T)^{-1/3}\right).
\]

## 5.2 New far-minus deletion

For minus with `k<=K_0`, fixed `J` gives `d>T/J` and `K_0=o(T)`, hence eventually
\[
k<d/2.
\]
Therefore
\[
h_-=d-k>\frac{T}{2J}.
\tag{5.4}
\]
In particular these rows are far outside the already-deleted centre strip
\[
h_-\le T^{1/2}/\omega(T).
\]
Thus (5.2) removes the genuinely new far-minus band
\[
\boxed{
 p-6d\le K_0(T),
 \qquad
 |p-7d|\asymp_J T.
}
\tag{5.5}
\]

An explicit nonempty scale is
\[
p-6d=T^{5/8},
\qquad d\asymp_JT,
\]
for which
\[
h_-=d-T^{5/8}\asymp_JT.
\]
This region is invisible to the old small-`h_-` theorem.

---

# 6. A second improvement: use actual zero-pair energy, not the marginal root envelope

The conditions
\[
p\mid b_x,
\qquad
p\mid N_h(x)
\]
are stronger than “`x` is a formal root of `N_h`”.

On the nonwrapping windows in (1.3)--(1.4), the restart theorem gives
\[
p\mid b_x,\ p\mid N_h(x)
\quad\Longleftrightarrow\quad
p\mid b_x,\ p\mid b_{x+h}.
\tag{6.1}
\]
So every actual row supplies an unordered pair of distinct Apéry zeros modulo `p`.

For fixed `p`, put
\[
\mathcal Z_p:=\{0\le j<p:p\mid b_j\},
\qquad
Z(p)=|\mathcal Z_p|.
\]
There are exactly
\[
\binom{Z(p)}2
\]
unordered pairs of distinct zero positions.

## Plus

A pair `(x,x+h)` determines `h`, and
\[
h=p-6d
\]
determines `d` uniquely if an integer solution exists. Hence
\[
\boxed{
\sum_d |X^+_{p,d}|
\le \binom{Z(p)}2.
}
\tag{6.2}
\]

## Minus

Here
\[
h=|p-7d|.
\]
For fixed `(p,h)` there are at most two candidate gaps
\[
d=(p-h)/7,
\qquad
 d=(p+h)/7.
\]
Thus
\[
\boxed{
\sum_d |X^-_{p,d}|
\le 2\binom{Z(p)}2.
}
\tag{6.3}
\]

The neighboring-characteristic condition `q|b_mirror` can only reduce these counts.

Using the proved pointwise zero bound
\[
Z(p)\ll p^{2/3},
\]
we get
\[
\sum_d(|X^-_{p,d}|+|X^+_{p,d}|)
\ll p^{4/3}.
\tag{6.4}
\]
The shell contains `O_J(T/log T)` lower primes and `p\asymp_JT`, so
\[
\#\{\text{all true sign rows}\}
\ll_J
\frac{T^{7/3}}{\log T}.
\tag{6.5}
\]
Weighting by `log p` gives
\[
\boxed{
\sum\log p
\ll_J T^{7/3}.
}
\tag{6.6}
\]

This is a real improvement over summing the marginal envelope
\[
\min\{h,w^{2/3}\}
\]
separately over every gap, which yields the previously recorded supercritical `T^(8/3)/log T` scale.

But (6.6) still does not close:
\[
\frac{T^{7/3}}{T^2/\log T}
=T^{1/3}\log T\to\infty.
\tag{6.7}
\]
Thus after taking full advantage of the *actual* `p`-zero-pair condition, one still needs a genuine average saving from the neighboring characteristic `q`.

---

# 7. Why the `q` condition cannot be replaced by a formal reflected-gap root

This is the first exact algebraic obstruction to a naive resultant attack.

Suppose `q` is odd and `0<=u<q`. Reflection gives
\[
b_{q-1-u}\equiv b_u\pmod q.
\]
If `q|b_u`, then the two actual zeros `u` and `q-1-u` imply a continuant root at gap
\[
g_q:=q-1-2u
\]
(up to exchanging the two endpoints when necessary).

But `g_q` is even. The continuant reflection identity
\[
N_g(-X-g-1)=(-1)^{g-1}N_g(X)
\]
therefore makes `N_g` odd about the half-integer centre, so
\[
\boxed{2X+g+1\mid N_g(X)}
\qquad(g\text{ even}).
\tag{7.1}
\]
At `X=u`,
\[
2u+g_q+1=q.
\]
Hence
\[
\boxed{q\mid N_{g_q}(u)}
\tag{7.2}
\]
holds for **every** candidate `u`, independently of whether `q|b_u`.

Therefore the implication
\[
q\mid b_u
\Longrightarrow
q\mid N_{g_q}(u)
\]
forgets all selected-zero information. A resultant or common-root construction that substitutes the formal condition (7.2) for the actual condition `q|b_u` cannot prove the required cross-characteristic saving: its `q`-side equation is tautologically satisfied by the structural centre factor.

This is the same actual-zero/formal-root distinction emphasized throughout the pushed orbit/restart proof.

---

# 8. Why current one-characteristic incidence and resultant tools do not close the remaining power

The failure here is not merely “we did not find the right summation”. There are exact scoped obstructions in the pushed source.

## 8.1 Row/column/restart bounds have a proved deep-pair model obstruction

`problems/3.2/gcdtail_result.tex` proves an explicit incidence-model proposition with the following properties:

- every column obeys the local `L^(2/3)+O(1)` interval bound;
- level `h` lies in at most `3(h-1)` columns;
- distinct levels have codegree at most one;
- nevertheless the number of well-separated deep collision triples is `(c+o(1))H^2` for an explicit `c>0`.

Thus the current row-degree bound, the full local column theorem, and even codegree one do **not** imply `o(H^2)` deep-pair sparsity. This is a theorem-level obstruction to trying to recover the missing factor in (6.7) from the existing one-characteristic incidence inequalities alone.

It is not a counterexample for the actual Apéry continuants; it proves exactly the narrower statement needed here: those inequalities, by themselves, are insufficient.

## 8.2 Scalar resultants lose the aligned root component

The same file proves a deep component-misalignment counterexample: one prime can divide two separated-block resultants through two different roots of the same first continuant, while the corresponding shifted factors have gcd one.

Consequently, a gcd of scalar resultants cannot count the aligned `x`-witnesses needed in `X^\pm`. It retains rational-prime support but loses the root-component label.

The correct invariant in `gcdtail_result.tex` is an **aligned split gcd**. But its second factorial moment misses singleton components. Our residual is precisely first-occurrence/actual-witness sensitive, so scalar resultant support cannot replace it.

## 8.3 The generic additive large sieve is already known to miss the cross-prime scale

In `proof.tex` around lines 3240--3375, the Fourier form of the CRT error is written explicitly, and the classical additive large sieve is recorded as giving
\[
O(N^2/\log N)
\]
at the top-prime scale, one logarithm larger than the corresponding required second-moment threshold. The manuscript isolates an AP--BDH-type dispersion statement as the missing cross-prime theorem.

For the present residual, (6.6) says that after using all actual `p`-pair information we still need a saving of size roughly
\[
\frac{1}{T^{1/3}\log T}
\]
from the neighboring `q`-selected condition. A generic large-sieve estimate that does not exploit additional Apéry-specific horizontal structure does not provide such a power saving.

## 8.4 Vertical bounds alone cannot control cross-prime alignment

`proof.tex` around lines 3160--3245 gives an explicit adversarial reflected zero-set model with `Z(p)=2` for every prime but with all prime columns aligned at one integer. It satisfies the vertical size/reflection constraints yet violates pointwise anti-concentration maximally.

Again, this is not an Apéry counterexample. It proves the scoped no-go relevant here: **vertical zero counts, reflection, nearest-neighbor exclusion, and gap-degree bounds do not imply horizontal decorrelation across different primes.**

The remaining `q|b_mirror` condition is exactly such horizontal data.

---

# 9. Exact residual after the stronger deletion

Let
\[
H_-(T)=T^{1/2}/\omega(T)
\]
be the already-proved minus centre cutoff, and define the stronger common margin cutoff
\[
K_0(T)=
\frac{T^{2/3}}
     {(\log T)^{1/3}\Omega(T)},
\qquad
\Omega(T)\to\infty.
\tag{9.1}
\]

After deleting (5.2), it is enough to consider:

## Minus far residual

\[
\mathcal R^-_{\rm new}(T)=
\left\{
(p,d,x):
\begin{array}{l}
p,q=p+d\text{ prime in the fixed-}J\text{ shell},\\
0<d<p/6,\\
h_-=|p-7d|>H_-(T),\\
k=p-6d>K_0(T),\\
w_-=\min(d,k),\\
0\le x<w_-,\\
p\mid b_x,\ p\mid N_{h_-}(x),\\
q\mid b_{w_--1-x}
\end{array}
\right\}.
\tag{9.2}
\]

## Plus far residual

Because `h_+=k`, the old `T^(3/5)` cutoff is superseded. The new residual is
\[
\mathcal R^+_{\rm new}(T)=
\left\{
(p,d,x):
\begin{array}{l}
p,q=p+d\text{ prime in the fixed-}J\text{ shell},\\
h_+=p-6d>K_0(T),\\
w_+=\min(h_+,p-h_+),\\
0\le x<w_+,\\
p\mid b_x,\ p\mid N_{h_+}(x),\\
q\mid b_{h_+-1-x}
\end{array}
\right\}.
\tag{9.3}
\]

Double-sign rows should still be assigned deterministically to one sign, as in the accepted `X^\pm` setup; nothing here requires counting one physical row twice.

---

# 10. Smallest remaining theorem

The exact first missing statement is now narrower than the previous marginal-capacity problem.

Define the actual `p` zero-pair set
\[
\mathcal P_p
=
\{(x,h):0\le x<x+h<p,\ p\mid b_x,\ p\mid b_{x+h}\}.
\]
The proved zero bound gives
\[
|\mathcal P_p|=\binom{Z(p)}2\ll p^{4/3}.
\]
The q6 geometry maps each member of `\mathcal P_p` to at most two minus candidate prime pairs and at most one plus candidate prime pair.

What remains is to show that the neighboring-prime zero condition selects only a power/logarithmically small fraction of these actual pairs on average.

A sufficient exact statement is:

> **Far matched-zero-pair dispersion.** For each fixed `J`, after the deletions (9.1)--(9.3),
> \[
> \sum_{(p,d,x)\in\mathcal R^-_{\rm new}(T)}1
> +
> \sum_{(p,d,x)\in\mathcal R^+_{\rm new}(T)}1
> =o_J\!\left(\frac{T^2}{\log^2T}\right).
> \tag{MZD}
> \]

Since `log p\asymp_J log T`, `(MZD)` is exactly equivalent to the required log-weight bound
\[
o_J(T^2/\log T).
\]

The key point is that `(MZD)` is no longer asking for a better one-prime zero count. It asks for cross-prime nonalignment between
\[
\{(x,h):p\mid b_x,b_{x+h}\}
\]
and
\[
\{m:q\mid b_m\},
\qquad q=p+d,
\]
under the affine q6 mirror map.

That is the first remaining theorem I cannot derive from the pushed unconditional toolkit.

---

# 11. Final assessment of the proposed routes

## Simultaneous `p|b_x`, `p|N_h(x)`

**Useful.** It improves the aggregate from the marginal `p^(5/3)` envelope to actual pair energy `p^(4/3)` per prime, giving (6.6). It does not close.

## Prime-pair geometry

**Useful in shallow affine margins.** It produced the earlier centre strips. More strongly, the matched-margin grouping above eliminates the ambient `T` factor and gives the new cubic-height deletion (4.2). In the fully macroscopic residual the prime-pair condition alone gives no further power.

## Average over `h,d,p`

**Stops at (6.6)** with current unconditional zero bounds. A further average saving would be precisely `(MZD)`.

## Large sieve

**No closure from the classical form.** The pushed proof already records the corresponding top-scale large-sieve logarithmic loss and formulates cross-prime dispersion as the missing input.

## Scalar resultants

**Invalid for aligned multiplicity.** Exact deep component-misalignment is proved in `gcdtail_result.tex`.

## Bilinear/column incidence

**Insufficient with the current axioms.** The affine-plane obstruction in `gcdtail_result.tex` proves this sharply for deep pairs.

## Formal reflected `q`-root

**Invalid substitute for selectedness.** The centre factor (7.1)--(7.2) is forced for every candidate `u`.

---

# Bottom line

There is no full unconditional far-residual closure from the currently pushed canonical tools.

There is, however, a substantial new unconditional deletion:

\[
\boxed{
\sum_{\substack{\text{true q6 sign rows}\\p-6d\le K}}
\log p
\ll_J K^3,
}
\]
so one may take
\[
\boxed{
K=o\!\left(T^{2/3}(\log T)^{-1/3}\right).
}
\]

This does two things at once:

- it improves the plus reflected-gap deletion from exponent `3/5` to the joint-height exponent `2/3`;
- it removes a new minus strip with `p-6d` this small while `|p-7d|` remains macroscopic.

After this deletion, the strongest unconditional aggregate use of the actual `p` zero-pair condition is
\[
\sum\log p\ll_JT^{7/3},
\]
and the exact remaining obstruction is cross-prime matched-zero-pair dispersion `(MZD)`.

The existing pushed sources contain theorem-level reasons why the obvious formal-root, scalar-resultant, one-characteristic incidence, and classical-large-sieve substitutions do not prove `(MZD)`. A new theorem must retain the **actual Apéry zero label in both characteristics simultaneously**.
