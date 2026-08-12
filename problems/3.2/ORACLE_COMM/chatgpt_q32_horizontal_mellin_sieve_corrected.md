# Erratum to Q7690: horizontal Mellin sieve with the correct quotient quantifiers

## Erratum

**[ERRATUM TO Q7690]** Sections 7 and 10 of Q7690 made a false quantifier upgrade.  The `q=1` top-half radical is **not** equivalent to the full all-index conjecture.

What `proof.tex` actually proves is the following.

- The quotient class `q=1`, equivalently `n/2 < p <= n`, is exactly the top-window contribution.  A top-half Mellin sieve controls that contribution only.
- Proposition `prop:quotient-reduction` removes the small-prime range
  \[
  p\le \frac{n}{f(n)\log n}
  \]
  for any `f(n) -> infinity` with, for example, `f(n) <= log n` eventually.
- The unresolved part of the **full** theorem is therefore the union of all quotient classes
  \[
  1\le q:=\left\lfloor\frac np\right\rfloor < Q_n,
  \qquad Q_n:=f(n)\log n.
  \]
  Uniform horizontal control is required throughout this entire slowly growing quotient range.

The valid part of Q7690 is the fixed-quotient Mellin identity and the bounded-geometric-complexity interpretation.  This note states those pieces with the correct quantifiers and gives the weakest sieve hypotheses that actually imply the all-index conclusion.

---

## 1. Exact fixed-quotient Mellin identity

Fix an integer `n>=1` and a prime `p<=n`.  Write the Euclidean decomposition

\[
\boxed{n=qp+r,\qquad q=\left\lfloor\frac np\right\rfloor\ge1,
\qquad 0\le r\le p-1.}
\]

Define the Mellin exponent used by the horizontal construction by

\[
\boxed{m:=n-q.}
\]

Then identically

\[
\boxed{m=q(p-1)+r,}
\]

so

\[
\boxed{m\equiv r\pmod{p-1}.}
\]

Let `M_p(m)` denote the normalized Mellin coefficient from the fixed primitive Apéry K3 trace system, twisted by the Teichmuller/Kummer character of exponent `m`, and let `mathfrak p_p` be the chosen prime of the coefficient field above `p` induced by the Teichmuller embedding.

### Theorem 1 (fixed quotient, interior characters)

**[THEOREM]** For every quotient `q>=1`, every good prime `p>=5`, and every interior residue

\[
1\le r\le p-2,
\]

with `n=qp+r` and `m=n-q`, one has

\[
\boxed{M_p(m)\equiv -b_r\pmod{\mathfrak p_p}.}
\tag{1.1}
\]

Since the residue field at the chosen Teichmuller prime is `F_p`, this is equivalently the displayed congruence modulo `p` after the prescribed reduction.

The point is that the right side depends on the Euclidean remainder `r`, while the Mellin character depends only on

\[
m\bmod(p-1)=r.
\]

There is no restriction `q=1` in (1.1).

### Endpoint restrictions

The nontrivial Kummer/Mellin statement (1.1) is an **interior-character theorem**.  The two endpoints are trivial-character cases and should not be silently folded into it.

- **`r=0`.** Then
  \[
  n=qp,\qquad m=q(p-1).
  \]
  The character is trivial.  The Apéry residue is
  \[
  b_0=1,
  \]
  so `p` is never a bad/common prime at this endpoint.

- **`r=p-1`.** Then
  \[
  n=(q+1)p-1,\qquad m=(q+1)(p-1),
  \]
  again the trivial character.  The Apéry palindromy/Lucas endpoint gives
  \[
  b_{p-1}\equiv b_0\equiv1\pmod p,
  \]
  so this endpoint is also never bad.

Thus **every actual zero-divisibility event in the quotient cells lies in the interior range `1<=r<=p-2`, where (1.1) applies**.

For the divisor problem `p<=n`, one always has `q>=1`; the formal class `q=0` corresponds to `p>n` and is not part of the quotient reduction.

### The top-half specialization

For `q=1`,

\[
\frac n2<p\le n,
\qquad r=n-p,
\qquad m=n-1=(p-1)+r.
\]

Therefore (1.1) becomes exactly the top-half identity used in Q7690.  It is a specialization of the fixed-quotient identity, not the whole all-index statement.

---

## 2. What quotient reduction actually leaves

Let

\[
Q_n=f(n)\log n,
\]

where

\[
f(n)\to\infty
\]

and, as in `proof.tex`, one may impose `f(n)<=log n` eventually.

Proposition `prop:quotient-reduction` proves that primes

\[
p\le \frac n{Q_n}
\]

make total logarithmic radical contribution `o(n)`.  Hence the only remaining prime range is

\[
\boxed{\frac n{Q_n}<p\le n,}
\]

and every such prime lies in one and only one quotient cell

\[
I_q(n):=\left(\frac n{q+1},\frac nq\right]\cap\mathbb P,
\qquad 1\le q<Q_n.
\tag{2.1}
\]

The intervals `I_q(n)` are disjoint.  This disjointness is the key summation fact: a correct all-index sieve must sum over all quotient cells, but it must **not** pay an artificial factor `Q_n` when the estimate is normalized by the prime mass of each cell.

Let `B_{n,q}` be the set of primes in `I_q(n)` satisfying the actual common/bad condition after the harmless endpoints above have been removed, and define the logarithmic bad mass

\[
W_{n,q}:=\sum_{p\in B_{n,q}}\log p.
\tag{2.2}
\]

Also put

\[
\Theta_{n,q}:=\sum_{p\in I_q(n)}\log p.
\tag{2.3}
\]

Because the quotient cells partition the surviving primes,

\[
\sum_{1\le q<Q_n}\Theta_{n,q}
\le \vartheta(n)=O(n).
\tag{2.4}
\]

---

## 3. The weakest sieve statement that implies the full theorem

### 3.1 Literal weakest aggregate form

**[CONDITIONAL — minimal aggregate PMLS]** The weakest direct statement needed from the horizontal Mellin sieve is

\[
\boxed{
\sum_{1\le q<Q_n}W_{n,q}=o(n)
}
\tag{PMLS-agg}
\]

uniformly for every integer `n`.

For the radical target, (PMLS-agg) is already the needed high-prime estimate.  For the logarithmic common-content target, use the existing safe valuation cap

\[
v_p(G_n)\le6
\]

on the relevant common primes.  Then

\[
\sum_{\substack{p>n/Q_n\\p\mid G_n}}v_p(G_n)\log p
\le6\sum_{q<Q_n}W_{n,q}=o(n).
\tag{3.1}
\]

The quotient-reduction tail contributes another `o(n)`, so (PMLS-agg) closes the full all-index theorem.  The factor `6` changes only the absolute constant.  Where `proof.tex` has a sharper local cap (for example `3` in a narrower range), it may of course be retained; it is unnecessary for this quantifier audit.

### 3.2 Weakest clean statement uniform in each quotient cell

A convenient uniform PMLS is a relative **log-weight saving**:

**[CONDITIONAL — uniform weighted PMLS]** There is a function `epsilon(n)->0` such that for every sufficiently large `n` and every

\[
1\le q<Q_n,
\]

one has

\[
\boxed{W_{n,q}\le\epsilon(n)\Theta_{n,q}.}
\tag{PMLS-w}
\]

Then, using (2.4),

\[
\sum_{q<Q_n}W_{n,q}
\le\epsilon(n)\sum_{q<Q_n}\Theta_{n,q}
\le\epsilon(n)\vartheta(n)=o(n).
\]

**There is no factor `Q_n` in this summation.**  The cells are disjoint and the right normalization is their Chebyshev/logarithmic mass.

This is the cleanest weak uniform statement.  A fixed positive density bound

\[
W_{n,q}\le c\Theta_{n,q},\qquad 0<c<1,
\]

is **not enough**: it yields only `O(n)`, not `o(n)`.  One needs an `o(1)` relative density/log-weight saving, or something stronger such as a power-saving count.

More generally, uniformity can be weakened to a weighted-average condition

\[
\sum_{q<Q_n}\epsilon_{n,q}\Theta_{n,q}=o(n),
\qquad W_{n,q}\le\epsilon_{n,q}\Theta_{n,q};
\]

`max_q epsilon_{n,q}->0` is sufficient but not logically necessary.

### 3.3 Count power saving: sufficient but stronger than needed

Suppose an HKKS-type estimate gives a fixed `delta>0` and, uniformly in every quotient cell,

\[
\boxed{
\#B_{n,q}\ll P_q^{1-\delta},
\qquad P_q:=\frac nq.
}
\tag{HKKS-count-q}
\]

Then

\[
W_{n,q}\ll (n/q)^{1-\delta}\log n,
\]

and summing over `q<Q_n` gives

\[
\sum_{q<Q_n}W_{n,q}
\ll n^{1-\delta}\log n
\sum_{q<Q_n}q^{-(1-\delta)}
\ll n^{1-\delta}\log n\,Q_n^{\delta}.
\tag{3.2}
\]

Hence a fixed power saving closes the theorem provided

\[
\boxed{
Q_n^{\delta}\log n=o(n^{\delta}).
}
\tag{3.3}
\]

Since `Q_n=f(n)log n`, an explicit sufficient condition is

\[
\boxed{
f(n)=o\!\left(\frac{n}{(\log n)^{1+1/\delta}}\right).}
\tag{3.4}
\]

The quotient reduction in `proof.tex` is much more conservative: it permits us to choose `f(n)->infinity` **arbitrarily slowly**, with `f(n)<=log n` eventually.  Under that choice, every fixed `delta>0` in (HKKS-count-q) is vastly more than enough.

If the saving is allowed to shrink, `delta=delta(n)->0`, the exact condition corresponding to (3.2) is

\[
\boxed{
\delta(n)\log\frac{n}{Q_n}-\log\log n\longrightarrow+\infty.
}
\tag{3.5}
\]

### 3.4 Density statements

An unweighted density statement can also suffice, but it must be **vanishing density**, not merely density `<1`.

Under the paper's harmless choice `f(n)<=log n`, all surviving primes satisfy

\[
p\ge n/(\log n)^2=n^{1-o(1)},
\]

so `log p=(1+o(1))log n` uniformly.  Consequently

\[
\#B_{n,q}=o(\#I_q(n))
\]

uniformly in `q<Q_n` implies the weighted version (PMLS-w).  Direct log-weight control is preferable because it remains valid without having to compare logarithms.

---

## 4. Dyadic prime blocks: the summation audit

The surviving range can also be partitioned into dyadic prime blocks

\[
P<p\le2P,
\qquad P\ge n/Q_n.
\]

Every prime in such a block still has a unique quotient `q=floor(n/p)`; summing over quotient cells and summing over dyadic `P` are two ways of partitioning the same prime set, not independent multiplicities.

Define the aggregate bad mass in one dyadic block by

\[
W_n(P):=
\sum_{\substack{P<p\le2P\\p\text{ bad at }q=\lfloor n/p\rfloor}}
\log p.
\]

A blockwise weighted statement

\[
\boxed{W_n(P)\le\epsilon(n)P,
\qquad\epsilon(n)\to0}
\tag{PMLS-dyad}
\]

uniformly over all surviving dyadic blocks gives

\[
\sum_P W_n(P)
\ll\epsilon(n)\sum_P P
=O(\epsilon(n)n)=o(n),
\]

because the dyadic scales form a geometric series.  There is **no extraneous `log n` factor from the number of dyadic blocks** when the estimate is proportional to `P`.

Likewise a blockwise count power saving

\[
\#\{p\in(P,2P]:p\text{ bad}\}\ll P^{1-\delta}
\]

gives

\[
\sum_P W_n(P)
\ll\sum_P P^{1-\delta}\log(2P)
\ll n^{1-\delta}\log n=o(n).
\]

This blockwise statement is stronger than (HKKS-count-q) because it has already aggregated all quotient cells represented inside the block; accordingly it does not carry the factor `Q_n^delta` from (3.2).

---

## 5. What a `q=1` THMS actually proves

For `q=1`, the prime cell is exactly

\[
I_1(n)=(n/2,n]\cap\mathbb P.
\]

A top-half Mellin sieve (THMS) proving

\[
\boxed{
\sum_{\substack{n/2<p\le n\\p\text{ bad}}}\log p=o(n)
}
\tag{THMS}
\]

has the following precise consequence.

**[THEOREM — q=1 consequence only]** The entire `q=1` / top-window radical contribution is `o(n)`.  With the valuation cap `6`, the corresponding top-window logarithmic common-content contribution is also `o(n)`.

**[NON-IMPLICATION]** THMS does **not** control any cell

\[
2\le q<Q_n.
\]

It therefore does **not** imply the full all-index theorem by itself and is not equivalent to it.  The full theorem needs either (PMLS-agg), the uniform weighted PMLS, or another estimate whose sum over all surviving quotient cells is `o(n)`.

This is the decisive correction to Q7690.

---

## 6. Trace-sheaf audit: what is fixed and what moves

Q7690's useful geometric point survives, but only after separating three layers which must not be conflated.

### 6.1 Fixed primitive K3 system

**[THEOREM/STRUCTURAL INPUT]** The underlying Apéry geometry is fixed: one starts from the fixed Laurent/K3 family and its primitive middle cohomology (after removing the fixed algebraic part).  Its geometric rank, the degree of the defining maps, and the geometric singular/branch divisor are independent of `p`, `q`, and `m`.

This fixed primitive system is the source of the bounded-complexity trace function.

### 6.2 The coefficient field and Teichmuller prime move with `p`

**[CORRECTION]** There is not literally one coefficient field and one reduction prime valid for every characteristic.  To realize a character of order dividing `p-1`, one works in a coefficient field containing the relevant roots of unity and chooses the Teichmuller prime

\[
\mathfrak p_p\mid p
\]

which realizes those roots in `F_p^*`.

Thus the reduction map in (1.1) is `p`-dependent.  This arithmetic movement of the coefficient field/prime is distinct from the fixed **geometric** complexity of the K3 system.

It also means that a horizontal sieve theorem must be stated for a compatible family across varying `p`, not justified by pretending that every trace value lives in one fixed finite coefficient field.

### 6.3 The Kummer local systems move with `(p,m)`

The Mellin coefficient twists the fixed primitive system by the rank-one tame Kummer local system

\[
\mathcal L_{\omega_p^{-m}}.
\]

These local systems vary with the prime and the character exponent.  Their rank is always one; they are tamely ramified and their singular support lies on a fixed geometric boundary after the fixed pullbacks in the construction.  Therefore the twist family can have uniformly bounded geometric complexity even though the individual local system is not fixed.

### 6.4 FKM conductor convention

In the Fouvry--Kowalski--Michel geometric convention for a middle-extension sheaf on a curve, the conductor is, up to harmless convention-level additive constants,

\[
\boxed{
\operatorname{cond}_{\mathrm{geom}}(\mathcal F)
=\operatorname{rank}(\mathcal F)
+\#\operatorname{Sing}(\mathcal F)
+\sum_x\operatorname{Swan}_x(\mathcal F).
}
\tag{6.1}
\]

This is a **geometric** conductor.  It does not include the degree or discriminant of the coefficient field.

For the fixed primitive K3 system tensored with the moving Kummer characters:

- the primitive rank is fixed;
- the number of geometric singularities is fixed up to the fixed pullback/branch set;
- Kummer ramification is tame, so it adds no growing Swan conductor.

Hence the intended statement is indeed

\[
\operatorname{cond}_{\mathrm{geom}}=O(1)
\]

uniformly in the good primes and nontrivial Mellin characters.

### 6.5 Punctual corrections do not disappear

**[CORRECTION]** Primitive extraction, middle extension, and fixed finite pushforwards can leave punctual/skyscraper contributions at finitely many singular or branch fibers.  They must not be declared nonexistent.

What is valid is that their total complexity remains uniformly bounded:

- only `O(1)` geometric points are involved;
- the relevant stalk dimensions are bounded by the fixed ranks/degrees;
- trivial/exceptional characters, including the two endpoint residues above, may be split off separately.

Thus the trace formula may be written as a bounded-conductor middle-extension trace plus `O(1)` punctual correction terms, uniformly in `p` and the nontrivial character.

This bounded punctual term is harmless for a properly formulated horizontal sieve, but it must be present in the sheaf audit.

### 6.6 What bounded conductor does **not** prove

Uniformly bounded FKM conductor is only an eligibility/complexity input.  It does not by itself prove (PMLS-w) or an HKKS power saving.  A genuine horizontal theorem must still supply the needed cancellation/non-concentration uniformly across the moving Kummer twists and varying characteristics; in particular one must control exceptional geometric isomorphisms and the relevant correlation sheaves, not merely quote bounded rank.

---

## 7. Corrected theorem preserved from Q7690

### Theorem 2 (uniform horizontal Mellin sieve implies the all-index estimate)

**[CONDITIONAL THEOREM]** Let `f(n)->infinity`, with `f(n)<=log n` eventually, and put `Q_n=f(n)log n`.  Suppose the actual Apéry horizontal Mellin bad sets satisfy either the aggregate estimate

\[
\sum_{1\le q<Q_n}W_{n,q}=o(n)
\]

uniformly for all `n`, or the stronger per-cell weighted estimate (PMLS-w).

Then the surviving high-prime radical contribution is `o(n)`.  Together with Proposition `prop:quotient-reduction`, the full radical estimate is `o(n)`.  If one tracks logarithmic common content rather than only the radical, the existing valuation cap `v_p<=6` multiplies the high-prime term by at most `6`, so it is still `o(n)`.

**Proof.** The quotient reduction gives `o(n)` below `n/Q_n`.  Above that threshold every prime lies in exactly one of the disjoint cells (2.1).  Summing their log-weight is precisely the left side of (PMLS-agg); under (PMLS-w), (2.4) bounds it by `epsilon(n) theta(n)=o(n)`.  The valuation cap gives (3.1). `square`

### Corollary 3 (what q=1 proves)

**[THEOREM]** If only THMS is known, then only the `q=1` top-window contribution is `o(n)`.  No conclusion about the full all-index theorem follows without additional control of every surviving `q>=2` cell, or an aggregate substitute covering them.

---

## 8. Corrected theorem DAG

```text
n = q p + r,  m = n-q = q(p-1)+r
        |
        +--> 1 <= r <= p-2:
        |       M_p(m) == -b_r (mod mathfrak p_p)    [all q]
        |
        +--> r=0 or p-1:
                trivial Kummer endpoint; b_r != 0 mod p

prop:quotient-reduction
        |
        +--> p <= n/(f log n): o(n)                  [proved]
        |
        +--> n/(f log n) < p <= n
                |
                +--> 1 <= q < f(n) log n
                |
                +--> q=1 THMS ------------------> q=1 contribution only
                |
                +--> uniform PMLS/HKKS over ALL surviving q
                           |
                           +--> sum_q W_{n,q} = o(n)
                           |
                           +--> radical high part = o(n)
                           +--> valuation cap 6 => content high part = o(n)
                           |
                           +--> full theorem

fixed primitive K3 geometry
        + moving coefficient field / Teichmuller prime mathfrak p_p
        + moving tame Kummer L_{omega_p^{-m}}
        + uniformly bounded punctual corrections
                |
                +--> FKM geometric conductor O(1)
                |
                +--> [still needed] horizontal cancellation/non-concentration
```

---

## 9. Final status

- **[ERRATUM]** Q7690's claimed equivalence between the full all-index conjecture and the `q=1` top-half radical is false.
- **[THEOREM]** The fixed-quotient Mellin identity is valid for every `q>=1` in the interior residue range, with harmless trivial-character endpoints.
- **[THEOREM]** `q=1` is exactly the top-window class and a THMS controls exactly that contribution.
- **[THEOREM]** Proposition `prop:quotient-reduction` leaves all `1<=q<f(n)log n`, so the full theorem requires an estimate whose **sum over those cells** is `o(n)`.
- **[CONDITIONAL]** Uniform `o(1)` log-weight density in each quotient cell is sufficient and incurs no factor `f(n)log n`; a fixed positive density is not sufficient.
- **[CONDITIONAL]** A fixed count power saving `P^(1-delta)` is stronger than necessary and is sufficient for any arbitrarily slowly growing `f` used in `proof.tex`; the exact growth condition is (3.3)/(3.4).
- **[CORRECTION]** The trace-sheaf family has fixed primitive K3 geometry but `p`-dependent coefficient fields/Teichmuller primes and moving tame Kummer twists.  The appropriate FKM **geometric** conductor remains `O(1)`; bounded punctual corrections remain and must be included.
- **[REMAINING GAP]** Bounded conductor alone does not furnish the uniform horizontal PMLS/HKKS estimate.  That uniform cancellation across all surviving quotient cells is the genuine all-index input still required.

This note amends Q7690 only.  It makes no change to `proof.tex`, `DOCTRINE`, or progress files.