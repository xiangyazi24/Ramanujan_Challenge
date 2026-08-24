ANSWER Q2694 dfccaee3

# P3.2 cold audit: centered reciprocal section / TCUS

Audited against current `main` HEAD `47fa0e653f52c4a71e9a8c26b31ca9f66f6bbe86`, in particular:

- `problems/3.2/oracleA_result.tex` (exact reciprocal form, SDC normalization, doublet trace-function discussion, reflected anchored obstruction),
- `problems/3.2/new_sections.tex` (weighted CED formulation and anchored-star obstruction),
- `problems/3.2/hm3_result.tex` (strong anchored rainbow star satisfying the currently proved local consequences),
- `problems/3.2/ORACLE_COMM/ced_proof_sketch.tex` (the earlier completion heuristic and its operator-level gap).

## Verdict

The new observation contains one exact and useful identity, but it does **not** supply an established theorem route.

1. The centering
   \[
   G_p(a)=F_p(a)+\frac{Z(p)}{p-1},\qquad a\in\mathbf F_p^*,
   \]
   is normalized correctly for the actual Apéry masks, because `b_0=1`, hence `0\notin\mathcal Z_p`, and
   \[
   \sum_{a\in\mathbf F_p^*}F_p(a)=-Z(p).
   \]
   Therefore `G_p` has exact mean zero on \(\mathbf F_p^*\).

2. For distinct primes \(p,\ell\), \((k,p\ell)=1\), and `v` ranging over all CRT units,
   \[
   \sum_{v\in(\mathbf Z/p\ell\mathbf Z)^*}
   G_p(kv^{-1})G_\ell(kv^{-1})=0
   \]
   exactly. This is just CRT factorization.

3. The arithmetic point used by the reciprocal parametrization is, however,
   \[
   v=p+\ell,
   \qquad
   v^{-1}\bmod p=\ell^{-1},\quad
   v^{-1}\bmod\ell=p^{-1}.
   \]
   It is a highly structured one-point section of each \((p\ell)\)-unit group, not a sample to which the complete mean can be transferred. There is no completion theorem saying that zero complete mean controls this moving section.

4. The exponent claim `theta<1/3 closes` is **algebraically correct only after a stronger premise than the stated TCUS**: namely, one needs an \(N^\theta\)-loss estimate for the **original reciprocal `F_pF_ell` form** (or one must separately control the two linear terms created when passing from `F` to `G`). Under that stronger premise the final exponent is
   \[
   |D_N(j)|\ll N^{5/6+\theta/2+o(1)},
   \]
   so strict \(\theta<1/3\) gives `o(N)`. At \(\theta=1/3\) there is no power saving.

5. TCUS for the centered product `G_pG_ell` is **not implied** by the currently proved local Apéry-mask properties. A reflected anchored family satisfying the repo's local size, reflection, endpoint, nonadjacency, interval, gap-certificate, bounded-fiber/energy, and even bounded row-codegree properties violates TCUS already in the `K=1` block, in fact for every fixed \(\theta<1\).

6. None of the established Kloosterman-fraction / trace-function / prime-large-sieve theorems I checked applies to the actual full Apéry form. The first missing theorem hypothesis is concrete: **fixed/factorized arithmetic data across the two varying moduli**. The actual coefficient selecting a numerator is a joint incidence `1_{r in Z_p} 1_{s in Z_ell}` and therefore depends on the denominator/modulus variables themselves. On doublets the local conductor is uniformly bounded, but the slope `r_p` still moves with the characteristic; current fixed-modulus trace-function theorems do not control that horizontal motion.

So this is not a new unconditional edge yet. It is a clean reformulation of the horizontal cross-prime obstruction, plus a correct `1/3` exponent threshold conditional on controlling the right (un-centered, or fully re-expanded) reciprocal form.

---

## 1. Exact normalization of `F_p` and `G_p`

Recall
\[
F_p(a)=\sum_{r\in\mathcal Z_p}e_p(ar),
\qquad
Z(p)=|\mathcal Z_p|.
\]
Because \(b_0=1\), the actual Apéry zero set does not contain `0`. Hence for every `r in Z_p`,
\[
\sum_{a\in\mathbf F_p^*}e_p(ar)=-1,
\]
and therefore
\[
\sum_{a\in\mathbf F_p^*}F_p(a)=-Z(p).
\]
Thus
\[
G_p(a)=F_p(a)+c_p,
\qquad c_p:=\frac{Z(p)}{p-1},
\]
satisfies
\[
\sum_{a\in\mathbf F_p^*}G_p(a)=0.
\]

The exact \(L^2\) normalization is also useful. Parseval gives
\[
\sum_{a\bmod p}|F_p(a)|^2=pZ(p),
\qquad F_p(0)=Z(p).
\]
Since \(\sum_{a\ne0}F_p(a)=-Z(p)\),
\[
\begin{aligned}
\sum_{a\in\mathbf F_p^*}|G_p(a)|^2
&=pZ(p)-Z(p)^2-\frac{Z(p)^2}{p-1}\\
&=pZ(p)\left(1-\frac{Z(p)}{p-1}\right).
\end{aligned}
\]
In particular it is \(\ll pZ(p)\), and under the proved `Z(p) << p^(2/3)` bound it is asymptotic to `p Z(p)` up to a relative `O(p^{-1/3})` error.

There is one formal point about the statement of TCUS: `omega` must be normalized. If TCUS quantifies over arbitrary `omega`, it is false by scalar rescaling. The natural interpretation is the bounded positive/smooth reciprocal-tent weight coming from `oracleA_result.tex`, with `|omega| << 1`; all conclusions below use that interpretation.

---

## 2. The complete CRT-unit average really is zero

For distinct primes \(p,\ell\), and \((k,p\ell)=1\), inversion is a bijection on units, and CRT identifies
\[
(\mathbf Z/p\ell\mathbf Z)^*\simeq\mathbf F_p^*\times\mathbf F_\ell^*.
\]
Therefore
\[
\begin{aligned}
&\sum_{v\in(\mathbf Z/p\ell\mathbf Z)^*}
G_p(kv^{-1}\bmod p)G_\ell(kv^{-1}\bmod\ell)\\
&\qquad=
\left(\sum_{a\in\mathbf F_p^*}G_p(ka)\right)
\left(\sum_{b\in\mathbf F_\ell^*}G_\ell(kb)\right)=0.
\end{aligned}
\]
This uses the coprimality condition on `k`. In the exact reciprocal form that condition is present. If a later dyadic simplification drops `(k,p ell)=1`, the zero-average statement no longer applies as written.

For the arithmetic section,
\[
v=p+\ell,
\]
one has
\[
v\equiv\ell\pmod p,
\qquad v\equiv p\pmod\ell,
\]
so
\[
G_p(k/v)=G_p(k\bar\ell),
\qquad
G_\ell(k/v)=G_\ell(k\bar p).
\]
This is exactly the moving reciprocal form already isolated in `oracleA_result.tex`.

The crucial distinction is dimensional: for a fixed pair `(p,ell)`, the complete average contains about `N^2` units, while the arithmetic construction uses **one** unit, `p+ell`. As `(p,ell)` move, the modulus also moves. There is no common finite group in which the set of all `p+ell` is an equidistributed subset. In CRT coordinates the section is the graph
\[
(\ell\bmod p,\ p\bmod\ell),
\]
and after inversion the graph is
\[
(\bar\ell\bmod p,\ \bar p\bmod\ell).
\]
That graph is precisely where additive reciprocity can create coherent phases.

---

## 3. Exact relation with the repo's reciprocal SDC form

`oracleA_result.tex` proves the exact identity
\[
R_{\ne}(P,M)=\frac1M\mathcal B_{P,M},
\]
where
\[
\mathcal B_{P,M}
=\sum_{\substack{P<p,\ell\le2P\\p\ne\ell}}
\frac{(\log p)(\log\ell)}{p\ell}
\sum_{\substack{0<|k|<p\ell/M\\(k,p\ell)=1}}
\left(1-\frac{M|k|}{p\ell}\right)
F_p(k\bar\ell)F_\ell(k\bar p).
\]
For the top range `P ~ M ~ N`, the `k`-range has length `asymp N`. Split it into dyadic blocks `k ~ K`.

Let
\[
T_N:=\sum_{N<p\le2N}Z(p).
\]
The vertical Parseval quantity in the repo is
\[
Q_N=
\sum_{N<p\le2N}(\log p)^2
\left(\frac{Z(p)}p-\frac{Z(p)^2}{p^2}\right).
\]
Since `p ~ N` and `Z(p) << p^(2/3)`, eventually `Z(p)/p=o(1)`, hence
\[
Q_N\asymp \frac{(\log N)^2}{N}T_N.
\]
This comparison is important: the `2/3` bound is used here only to ensure the factor `1-Z(p)/p` is bounded away from zero. Its main exponent role comes later.

Suppose, for the moment, that a dyadic estimate of the form
\[
\operatorname{Re} C_K
\ll K N^{\theta+o(1)}T_N
\tag{*}
\]
were available for the **actual summand needed in \(\mathcal B\)**, with the outer `(log p log ell)/(p ell)` stripped off and the tent absorbed into a bounded `omega`.

Because dyadic `K` run up to `O(N)`,
\[
\sum_{K\ \mathrm{dyadic}}K\asymp N.
\]
Thus (*) gives
\[
\operatorname{Re}\mathcal B_{N,N}
\ll
\frac{(\log N)^2}{N^2}
N^{1+\theta+o(1)}T_N
\asymp N^{\theta+o(1)}Q_N.
\]
Consequently
\[
\operatorname{Re}R_{\ne}(N,N)
\ll N^{-1+\theta+o(1)}Q_N.
\]
For `theta=0` this is the SDC scale in `oracleA_result.tex`. For `theta>0` it is weaker than SDC by `N^theta`.

This point matters: TCUS(theta) with positive theta does **not** prove SDC. It can nevertheless still be enough for the final P3.2 growth target.

---

## 4. The `theta < 1/3` exponent threshold is correct, conditionally

The dilation argument in `oracleA_result.tex` uses the Fejér coefficient
\[
\gamma_{4N}(j)\gg N^{-2},\qquad N<j\le2N.
\]
If the short-arc energy is allowed an `N^theta` loss, the same argument gives
\[
\sum_{N<j\le2N}|D_N(j)|^2
\ll N^{1+\theta+o(1)}Q_N.
\]
Hence, pointwise,
\[
|D_N(j)|
\ll N^{(1+\theta)/2+o(1)}Q_N^{1/2}.
\]
The proved vertical zero bound gives
\[
Q_N\ll N^{2/3}\log N=N^{2/3+o(1)}.
\]
Therefore
\[
|D_N(j)|
\ll N^{5/6+\theta/2+o(1)}.
\]
Thus
\[
5/6+\theta/2<1
\quad\Longleftrightarrow\quad
\theta<1/3.
\]
So the advertised `1/3` threshold is not an arithmetic miracle; it is exactly the remaining `1/6` pointwise exponent margin after taking the square root of the `Z(p) << p^(2/3)` vertical energy.

At `theta=1/3` the result is only `N^(1+o(1))`; that does not prove `o(N)` without an additional logarithmic/power saving. The strict inequality is essential.

---

## 5. But TCUS on `G_p G_ell` is not the estimate needed for the original form

This is the first normalization gap in the proposal.

Write
\[
c_p=\frac{Z(p)}{p-1},\qquad F_p=G_p-c_p.
\]
Then pointwise
\[
F_pF_\ell
=G_pG_\ell-c_pG_\ell-c_\ell G_p+c_pc_\ell.
\]
Therefore a TCUS bound only for
\[
\sum\omega\,G_p(k\bar\ell)G_\ell(k\bar p)
\]
does **not** by itself bound the reciprocal form in `oracleA_result.tex`. One additionally needs the two one-sided section sums
\[
\sum\omega\,c_pG_\ell(k\bar p),
\qquad
\sum\omega\,c_\ell G_p(k\bar\ell),
\]
at a compatible scale.

The constant term `c_p c_ell` is harmless under `Z(p) << p^(2/3)`. The linear terms are not killed by complete mean-zero, because the arithmetic section again samples only `p mod ell` / `ell mod p`, not a complete unit set.

A trivial estimate for one linear dyadic block is of order
\[
K\frac{T_N^2}{N},
\]
whereas the TCUS target scale is
\[
K N^\theta T_N.
\]
Their ratio can be as large as
\[
\frac{T_N}{N^{1+\theta}}
\ll \frac{N^{2/3-\theta}}{\log N},
\]
which is not small in the desired `theta<1/3` regime. Parseval/Cauchy improves the threshold only to the same `1/3` boundary; it does not make the linear terms automatically negligible below it.

Therefore one of the following must be proved:

- a TCUS-type bound directly for the original `F_pF_ell` reciprocal block; or
- TCUS for `G_pG_ell` **plus** separate linear-section estimates.

The complete zero average of `G_pG_ell` does not provide those linear estimates.

A slightly cleaner centering of the original product is
\[
F_p(a)F_\ell(b)-c_pc_\ell,
\]
whose complete CRT-unit average is already zero. This avoids introducing the two linear terms algebraically, but it does not solve the arithmetic-section problem below.

---

## 6. A reflected anchored family violates TCUS itself

The repo already contains anchored countermodels for SDC/HM3. They can be tuned so that the new reciprocal section exhibits an explicit coherent diagonal.

Let
\[
M=N^2,
\qquad
m_0=N^2-1,
\qquad
d=2m_0+1=2N^2-1.
\]
Fix a small constant `eta>0`, say `eta=10^{-3}`, and keep only primes
\[
N<p\le(1+\eta)N,
\]
excluding the `O(1)` primes dividing
\[
E=(m_0-1)m_0(m_0+1)(m_0+2)(2m_0+1).
\]
For an active prime put
\[
a_p=m_0\bmod p,
\qquad
\mathcal Z_p^*=\{a_p,p-1-a_p\},
\]
and use the empty mask on all other primes in `(N,2N]`.

Let `L` be the number of active primes. By the PNT for the fixed proportional interval,
\[
L\asymp \frac{N}{\log N},
\qquad
T_N^*=2L\asymp\frac{N}{\log N}.
\]

### Local hypotheses

This is the same reflected anchored construction as `hm3_result.tex`, with `m_0` chosen near `M`.

- Each nonempty mask has size `2`, hence satisfies `Z(p) << p^(2/3)` and every interval bound `|Z_p cap I| <= 2 |I|^(2/3)`.
- It is reflection invariant under `r -> p-1-r`.
- The exclusions make the two residues distinct and avoid `0,1,p-2,p-1`; in particular there are no consecutive zeros.
- The reflection identity for the gap polynomial gives the necessary gap certificate for the unique pair, exactly as in `hm3_result.tex`.
- The value map can be extended with all fibers of size at most `2` and collision energy at most `2p`, again exactly as in `hm3_result.tex`.
- Since `d=2m_0+1>M`, the row-codegree proof in `hm3_result.tex` applies verbatim: distinct rows have `O(1)` common active primes.
- Every active prime nevertheless marks the same row `m_0`.

Thus this family satisfies strictly more than the local hypotheses being invoked by the proposed TCUS route.

### Reciprocal value at `k=1`

For this doublet write
\[
A=m_0=N^2-1,
\qquad
B=-m_0-1=-N^2.
\]
Then
\[
F_p^*(x)=e_p(Ax)+e_p(Bx).
\]
For distinct active `p,ell`, additive reciprocity gives
\[
\frac{\bar\ell}{p}+\frac{\bar p}{\ell}
\equiv\frac1{p\ell}\pmod1.
\]
Hence the two equal-slope terms in
\[
F_p^*(\bar\ell)F_\ell^*(\bar p)
\]
are
\[
e\left(\frac{A}{p\ell}\right)
+e\left(\frac{B}{p\ell}\right).
\]
Because `p,ell in (N,(1+eta)N]`, both `A/(p ell)` and `-B/(p ell)` lie within `O(eta)` of the integer `1`. Choosing fixed small `eta`, these two terms have uniformly positive real part, arbitrarily close to `2` as `eta -> 0`.

The two cross terms are genuine Kloosterman-fraction sums. For example
\[
\begin{aligned}
\frac{A\bar\ell}{p}+\frac{B\bar p}{\ell}
&\equiv
(A-B)\frac{\bar\ell}{p}+\frac{B}{p\ell}\\
&=d\frac{\bar\ell}{p}+\frac{B}{p\ell}\pmod1.
\end{aligned}
\]
The factor `e(B/(p ell))` is a smooth bounded weight on the fixed short proportional rectangle. Standard two-dimensional partial summation reduces it to the classical bilinear Kloosterman-fraction estimate. The Duke--Friedlander--Iwaniec bound already gives, in the balanced range with numerator `d ~ N^2`, a power saving `N^{-1/48+o(1)}` over the trivial bilinear scale. Bettin--Chandee improves the balanced saving to `1/20`; the 2026 Dong--Robles--Zeindler preprint improves the fixed-numerator bilinear saving further to `1/12`. Any one of these is more than enough here.

With prime-indicator coefficients, DFI gives the cross terms
\[
\ll N^{95/48+o(1)}=o\left(\frac{N^2}{(\log N)^2}\right)=o(L^2).
\]
The centering correction is also negligible because
\[
c_p=\frac2{p-1}=O(N^{-1}),
\qquad |F_p^*|\le2,
\]
so replacing `F_p F_ell` by `G_p G_ell` changes the total `k=1` pair sum by `O(L^2/N)`.

Therefore, for `omega=1` (and equally for the actual positive reciprocal tent at `k=1`, which is `1+O(1/N)`),
\[
\operatorname{Re} C_1^*=(2+o(1))L^2
\asymp\frac{N^2}{(\log N)^2}.
\]
But TCUS(theta) predicts
\[
\operatorname{Re}C_1^*
\ll N^{\theta+o(1)}T_N^*
\asymp\frac{N^{1+\theta+o(1)}}{\log N}.
\]
The ratio is
\[
\asymp \frac{N^{1-\theta-o(1)}}{\log N}\to\infty
\]
for every fixed `theta<1`.

So the anchored reflected family violates TCUS much more strongly than needed to rule out the proposed `theta<1/3` deduction from local hypotheses.

This counterexample is not claimed to be the actual Apéry zero family. Its role is sharper: it proves that **all of the currently available local information is compatible with failure of TCUS**. A proof for the real masks must therefore use genuinely horizontal Apéry arithmetic across different primes.

---

## 7. Why established bilinear Kloosterman-fraction theorems do not apply to the actual masks

There are two superficially relevant theorem families.

### 7.1 Duke--Friedlander--Iwaniec / Bettin--Chandee / Dong--Robles--Zeindler

The classical bilinear theorem controls forms of the shape
\[
\sum_{m\sim M}\sum_{n\sim N}
\alpha_m\beta_n\,e\left(\frac{a\bar m}{n}\right),
\]
with a fixed numerator `a` (allowed to grow within the theorem's parameter range) and factorized arbitrary coefficients `alpha_m beta_n`.

Bettin--Chandee adds a third independent numerator variable:
\[
\sum_{a,m,n}\nu_a\alpha_m\beta_n
e\left(\frac{a\bar m}{n}\right).
\]
The January 2026 Dong--Robles--Zeindler preprint `arXiv:2601.00292` strengthens the balanced bilinear saving, but retains the same essential fixed-numerator/factorized-coefficient architecture.

Expanding the actual Apéry product gives
\[
\sum_{p,\ell}
\sum_{r\in\mathcal Z_p}
\sum_{s\in\mathcal Z_\ell}
 e_p(kr\bar\ell)e_\ell(ks\bar p).
\]
By reciprocity one can rewrite one phase as
\[
e\left(k(r-s)\frac{\bar\ell}{p}\right)
\,e\left(\frac{ks}{p\ell}\right).
\]
The would-be numerator `a=k(r-s)` is **not independent of `(p,ell)`**: `r` is selected by the `p`-mask and `s` by the `ell`-mask. Equivalently, the coefficient is a joint incidence tensor
\[
\mathbf1_{r\in\mathcal Z_p}\mathbf1_{s\in\mathcal Z_\ell},
\]
not a product `alpha_ell beta_p nu_a` after grouping by `a=k(r-s)`.

This is the first missing hypothesis for applying these theorems. Treating that joint incidence tensor as if it factorized is exactly where the anchored family slips through.

The counterexample above is instructive: DFI/BC *does* control its two cross terms, because there the numerator collapses to the fixed `d=2m_0+1`; what survives is the pair of equal-slope terms, which are coherent by additive reciprocity and of size `~L^2`. That is exactly the component a generic Kloosterman-fraction cancellation theorem cannot erase.

### 7.2 KMS / Fouvry--Kowalski--Michel(-Sawin) trace-function bilinear estimates

Kowalski--Michel--Sawin (Annals 2017) proves strong bilinear bounds for (hyper-)Kloosterman trace functions modulo a **fixed prime modulus**. Fouvry--Kowalski--Michel's `Algebraic trace functions over the primes` likewise fixes the finite field/sheaf while the prime argument varies. The newer Fouvry--Kowalski--Michel--Sawin preprint `arXiv:2511.09459` extends bilinear estimates to broad bounded-conductor trace functions under monodromy hypotheses, but it still starts from a trace function on one finite field / one modulus and controls bilinear arguments inside that field.

For a full Apéry mask,
\[
a\mapsto F_p(a)=\sum_{r\in\mathcal Z_p}e_p(ar)
\]
can be represented as a direct sum of `Z(p)` rank-one Artin--Schreier traces, so its naive conductor/rank grows with `Z(p)` and can be as large as `p^(2/3)` under present theorems. This is outside the uniform bounded-conductor regime.

On the doublet sector, `oracleA_result.tex` correctly observes that this particular problem disappears: `F_p` is a rank-2 bounded-conductor trace function. But the **slope `r_p` moves with `p`**, and TCUS couples two different characteristics at once:
\[
G_p(k\bar\ell)G_\ell(k\bar p).
\]
There is no fixed sheaf over a fixed finite field whose bilinear estimate is being invoked. The anchored doublet family above has the same rank-2 bounded-conductor local description at every prime and still violates TCUS. Therefore bounded local conductor/monodromy alone cannot be the missing input.

What is missing is a horizontal compatibility theorem for the *selection of the slopes/zero fibers across characteristics*.

---

## 8. Prime large sieve does not bridge the gap either

The ordinary additive large sieve controls the collection of all fractions `a/p`, and the repo already uses it in the range `P <= 2 sqrt(N)`. At `P ~ N` its natural conductor is `P^2 ~ N^2`, which is precisely too large to deliver the `N^(1+theta)` energy needed for `theta<1/3`.

The complete-unit identity does not change that conductor. To exploit the special section `v=p+ell`, one would need a new large-sieve theorem whose hypotheses know that the modulus-dependent functions `G_p` are horizontally incompatible with anchored selections. No standard prime large sieve has such a hypothesis or conclusion.

---

## 9. The first genuinely missing arithmetic hypothesis

There are two distinct levels of missing input.

### Formal/analytic level

To deduce the original reciprocal estimate from the proposed centered TCUS, one must add either

1. a bound for the two linear centered section sums, or
2. a TCUS statement formulated directly for the original `F_pF_ell` product (or `F_pF_ell-c_pc_ell`).

Without this, the `theta<1/3` exponent calculation is being applied to a different bilinear form from the one in the exact Fourier identity.

### Arithmetic/theorem-applicability level

Even after fixing that algebra, the first missing nonlocal hypothesis is:

> **Horizontal anti-alignment of the modulus-dependent zero fibers/slopes across distinct primes.**

Concretely, an established Kloosterman-fraction theorem would need the numerator/coefficient data to factor independently of the denominator variables; an established trace-function theorem would need a fixed bounded-complexity sheaf/family to which the two varying characteristics can be compared. The actual Apéry masks presently provide neither. The anchored reflected family shows that all currently proved local consequences allow a common integer anchor, and additive reciprocity then turns that anchor into a coherent `K=1` contribution of size `~(#primes)^2`.

This is not just a renaming of SDC. It identifies the exact theorem hypothesis that fails: **the mask-selection parameter is allowed to move adversarially with the modulus**. Any valid route must derive, from Apéry-specific arithmetic, a restriction on that horizontal motion before an off-the-shelf Kloosterman/trace-function estimate can be applied.

A useful intermediate theorem would therefore be something genuinely lower-level than SDC, for example a structural statement that decomposes the actual mask Fourier family into finitely/controlled-many horizontal components whose slopes/numerators are independent of the modulus variables (or arise from a fixed compatible sheaf family), with an error small in the total `L^2` mass. Once such a decomposition exists, DFI/BC/Dong--Robles--Zeindler or modern trace-function bilinear estimates could plausibly act on the resulting factorized pieces. Nothing currently in the repo proves such a decomposition.

---

## 10. Bottom line for P3.2

- **Keep** the centered identity and the observation `v=p+ell`: they are exact and clarify the geometry of the reciprocal form.
- **Keep** the exponent calculation: an `N^theta` loss in the *correct* reciprocal short-arc estimate still closes P3.2 for every strict `theta<1/3` because `Z(p)<<p^(2/3)` leaves a `1/6` pointwise margin.
- **Do not claim** that zero complete CRT-unit mean supplies TCUS on the arithmetic section.
- **Do not claim** that TCUS for `G_pG_ell` alone closes the original `F_pF_ell` problem; the linear centering terms remain.
- **Do not cite** DFI/Bettin--Chandee/Dong--Robles--Zeindler or KMS/FKM(FS) as applying directly to the real Apéry masks. Their fixed/factorized data hypotheses are absent.
- The reflected anchored family above satisfies the stated local hypotheses and violates TCUS already at `K=1`, so any proof must use a new **cross-prime Apéry-specific horizontal selection theorem**.

### Literature checked

- Duke--Friedlander--Iwaniec, bilinear forms with Kloosterman fractions (balanced saving `1/48` in the classical bound).
- S. Bettin and V. Chandee, *Trilinear forms with Kloosterman fractions*, Adv. Math. 328 (2018), arXiv:1502.00769 (balanced saving `1/20` in the `A=1` regime).
- E. Kowalski, P. Michel, W. Sawin, *Bilinear forms with Kloosterman sums and applications*, Ann. of Math. 186 (2017), 413--500.
- E. Fouvry, E. Kowalski, P. Michel, *Algebraic trace functions over the primes*, Duke Math. J. 163 (2014), 1683--1736.
- E. Fouvry, E. Kowalski, P. Michel, W. Sawin, *Bilinear forms with trace functions*, arXiv:2511.09459 (current 2026 preprint version).
- A. Dong, N. Robles, D. Zeindler, *Bilinear forms with Kloosterman fractions and applications*, arXiv:2601.00292 (Jan. 2026; stronger balanced fixed-numerator bilinear bound, but still not the modulus-coupled Apéry mask form).

The cold-audit answer is therefore: **TCUS(theta) is not presently a theorem route; the `theta<1/3` exponent is right only after fixing the centering/algebra, and an anchored reflected family satisfying the local hypotheses violates TCUS. The first missing theorem hypothesis is horizontal/factorized control of the moving mask slopes across primes.**