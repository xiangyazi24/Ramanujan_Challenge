# Q8345 — far physical one-label gateway: mixed-characteristic obstruction and surviving interface

## Verdict

I do **not** obtain an unconditional proof of GH1.  I also do not obtain a genuinely smaller positive theorem from a known result with all moving-characteristic hypotheses verified.  I obtain option **(3)** in the question: a theorem-level obstruction to three precise proposed routes, together with the first arithmetic interface that survives all three obstructions.

The decisive point is that a far physical edge is intrinsically **mixed-characteristic**.  If its endpoints are `q` and `ell`, then one has

```text
q | b_{h_q},        ell | b_{h_ell},
c_q q - c_ell ell = +/- g,
```

but one does **not** obtain

```text
q | b_{h_ell}       or       ell | b_{h_q}.
```

The existing fixed-gap continuant/resultant theorem requires two zeros in the **same** characteristic.  The integer equation joining `q` and `ell` does not change the residue field of either Apéry zero.  This is not merely a formal objection: an actual Apéry pair at

```text
X=128, m=321, q=179, ell=193
```

has folded zeros `h_q=36`, `h_ell=64`, signed lifts `s_q=-37`, `s_ell=-65`, and far gap `g=28`; both rows have Kummer order greater than the Q8345 cutoff `D` for every `0<eta<1/15`, and they are `D`-isolated.  Nevertheless

```text
179 does not divide b_64.
```

So the physical edge is **not** a gap-28 return modulo 179.  It cannot be fed into the same-characteristic continuant/resultant by forgetting which endpoint characteristic owns which zero.

There is a second, quantitative obstruction.  Even if one grants a fictitious transfer that made every gap-`g` endpoint divide one nonzero characteristic-zero resultant of logarithmic height `X^{o(1)} g^alpha`, a pointwise divisor argument can beat the GH1 vertex threshold at the first far scale only if

\[
\boxed{
\alpha<\alpha_*(\eta)
:=\frac{13/15+2\eta}{14/15+\eta}
=rac{13+30\eta}{14+15\eta}<1.
}
\tag{V1}
\]

Thus even a **linear** logarithmic-height carrier is too large for every fixed `eta<1/15`.  The banked fixed-gap continuant/resultant has roughly quadratic logarithmic height, so it misses the threshold by exactly one full power of `X` at the exponent level:

\[
2\left(\frac{14}{15}+\eta\right)
-\left(\frac{13}{15}+2\eta\right)=1.
\tag{V2}
\]

This remains fatal before one worries about constants, ordered/unordered factors, or the growth of `g` all the way to `2X`.

Third, the high-Kummer-order filter does not make an existing all-character Mellin/sheaf theorem applicable to the selected physical labels.  For a shell prime `p`, the number of multiplicative characters of order at most `D` is exactly

\[
\sum_{\substack{d\mid p-1\\d\le D}}\varphi(d)
\le \frac{D(D+1)}2=O(G),
\]

while

\[
\frac Gp
=X^{-1/15+\eta+o(1)}(\log X)^{-2}=o(1).
\]

Therefore deleting the low-order sector leaves `1-o(1)` of the character group.  An all-character equidistribution theorem can still tolerate **one exceptional high-order selected character for every prime**: that deterministic section has mass `1/(p-1)` in the character average.  Moreover GH1 asks for defining-prime divisibility at `ell=p`, whereas the standard compatible-system large sieve uses auxiliary `ell != p`.  The Kummer-order condition repairs neither quantifier mismatch.

The first surviving arithmetic object is consequently not a fixed-gap resultant.  It is the **mixed-characteristic endpoint correlation** itself.  For an oriented signed endpoint `s`, put

\[
H(s)=\begin{cases}s,&s\ge0,\\-s-1,&s<0.\end{cases}
\]

and, if `m=c_q q+s`, define the two candidate partner integers

\[
N_{q,\delta}(m,g)=c_q q-\delta g,
\qquad \delta\in\{+1,-1\}.
\]

If `s_ell=s+delta g`, then exactly

\[
N_{q,\delta}(m,g)=c_\ell\ell.
\]

There is at most one shell-prime divisor `ell` of this integer.  The actual edge event is therefore a deterministic one-candidate test coupling

\[
\boxed{
q\mid b_{H(s)},\qquad
\ell\mid b_{H(s+\delta g)},\qquad
\ell\mid N_{q,\delta}(m,g),
}
\tag{MC1}
\]

together with the already stated isolation, Kummer-order, quotient, sign, and folding masks.  `MC1` lives in two defining characteristics and preserves the physical integer relation.  None of the banked recurrence, reflection, Lucas, Hasse--Franel, or fixed-gap resultant results estimates this two-characteristic deterministic section at the GH1 scale.

A Cauchy/second-moment attempt does not make this interface smaller.  For fixed `(m,g)`, let `E_g(m)` be the number of **unordered** admissible edges and `V_g(m)` the number of incident labels.  The signed coordinates lie on a line, so the gap-`g` graph is a path forest and

\[
\boxed{
\frac{V_g(m)}2\le E_g(m)<V_g(m).
}
\tag{V3}
\]

If `T>=4` and `V_g(m)>T`, then

\[
\boxed{
V_g(m)-T
\le \frac{16}{T}\binom{E_g(m)}2.
}
\tag{V4}
\]

Pairs of distinct edges are exactly of two types: a shared-vertex three-label chain or a disjoint four-label rectangle.  Since the number of shared-edge pairs is at most `V_g(m)`, high excess forces many disjoint four-label rectangles.  Thus squaring the one-label gateway recovers the same arity as the existing top-strip four-prime rectangle gateway; it is not a genuinely smaller arithmetic premise.

So the status is:

- **PROVED:** exact signed/triangular identities, shell-partner uniqueness, path-forest geometry and ordered/unordered factors, the high-order character-count calculation, the pointwise carrier exponent barrier, and the reduction of any edge-pair second moment to three-label chains plus four-label rectangles;
- **REFUTED:** direct application of the fixed-gap continuant/resultant to a mixed-characteristic far edge; the idea that quadratic (or even linear) per-gap height is sufficient uniformly at `g>G`; the idea that `ord(omega_p^h)>D` lets all-character Mellin equidistribution control the selected deterministic labels; and the idea that the triangular carrier or the canonical affine second coordinate automatically supplies another target equation;
- **CONDITIONAL / FIRST SURVIVING INTERFACE:** a genuinely mixed-characteristic estimate for `MC1`, or an equivalent mixed-characteristic carrier/dispersion theorem retaining the physical relation `c_q q-c_ell ell=+/-g`.  This is earlier than a four-prime rectangle but is not currently supplied by a known theorem.

---

## 1. Source audit and visibility boundary

The connector-visible default branch is

```text
main@734a5a84c1e4fd8703a811aadaa2b4c7f532b20e
```

and I also inspected the current `chatgpt-drop` branch.  The literal strings `GH1`, `V_g`, and the new far-physical file names are not present in the connector-visible committed tree, so I treat the complete Q8345 statement as the authoritative interface for the newest reduction and do not invent unseen mask definitions or the value of `Y_X`.

The committed sources used in the audit include

```text
problems/3.2/meso_resultants.md
problems/3.2/FABLE_SECTION_sigma_half.tex
problems/3.2/ORACLE_COMM/chatgpt_q32_scale_sensitive_crossrow.md
problems/3.2/ORACLE_COMM/chatgpt_q7736_affine_codim2_sieve.md
problems/3.2/ORACLE_COMM/chatgpt_q32_padic_horizontal_mellin.md
problems/3.2/ORACLE_COMM/chatgpt_q32_fixed_exponent_mellin_attack.md
problems/3.2/hasse_franel_descent.tex
problems/3.2/research/scripts/q7699_padic_horizontal_mellin_n321.csv
problems/3.2/zp_million_output.txt
```

and, on `chatgpt-drop`, the immediately preceding physical-carrier audit

```text
problems/3.2/ORACLE_COMM/chatgpt_q8336_physical_racah_resultant_obstruction.md
```

No shared TeX, `DOCTRINE.md`, or `RUN_LOG_P32.md` is edited.  The owned files for Q8345 are

```text
problems/3.2/ORACLE_COMM/chatgpt_q8345_far_physical_one_label_obstruction.md
problems/3.2/ORACLE_COMM/chatgpt_q8345_far_physical_one_label_verify.py
```

---

# I. PROVED

## 2. Exact signed lift and triangular carrier

For each retained label `p`, Q8345 supplies a signed lift

\[
m=c_pp+s_p,
\qquad -X\le s_p\le X-1.
\tag{2.1}
\]

Let `h_p` be its folded Apéry-zero index.  The two possible signed representatives are

\[
s_p=h_p
\quad\text{or}\quad
s_p=-h_p-1.
\tag{2.2}
\]

In either case

\[
\boxed{\tau(s_p)=\tau(h_p).}
\tag{2.3}
\]

Therefore the triangular carrier factorizes **over the integers** as

\[
\begin{aligned}
\tau(m)-\tau(h_p)
&=\tau(m)-\tau(s_p)\\
&=(m-s_p)(m+s_p+1)\\
&=\boxed{c_pp\,(m+s_p+1)}.
\end{aligned}
\tag{2.4}
\]

This explains precisely what the physical triangular carrier contributes: it certifies the endpoint prime `p` from its own signed row.  It does not mention a second shell prime.

The Q8336 integral Racah/Newton audit sharpens this observation.  Its canonical physical-point lift satisfies

\[
\widehat R_h(\tau(m))
=(h!)^2b_h+(\tau(m)-\tau(h))Q_h(\tau(m)),
\]

so, once `p|tau(m)-tau(h)` is imposed and `p>h`, its `p`-divisibility is exactly the original event `p|b_h`.  Thus the triangular carrier plus the canonical physical interpolation row does not create a second local equation.

---

## 3. Exact partner uniqueness from physical height

Fix `(m,q,g)` and an orientation `delta in {+1,-1}`.  If

\[
s_\ell=s_q+\delta g,
\]

then subtracting the two signed lifts gives

\[
\boxed{
c_\ell\ell
=m-s_\ell
=c_qq-\delta g.}
\tag{3.1}
\]

This is the divisor statement in Q8345.

There is at most one shell prime dividing the right side.  Indeed, from
`m<X^2` and `s_ell>=-X`,

\[
0<m-s_\ell<X^2+X.
\tag{3.2}
\]

If two distinct shell primes `p_1,p_2>X` both divided it, then

\[
p_1p_2\ge(X+1)^2=X^2+2X+1>X^2+X,
\]

contradicting (3.2).  Thus for each `(m,q,g,delta)` the partner candidate is unique if it exists.

This is a strong deterministic simplification: GH1 is not summing over an uncontrolled second prime after `q,g` are fixed.  But uniqueness is not cancellation; it leaves one arithmetic primality/Apéry-zero test per oriented endpoint.

---

## 4. The two folded-root shapes

Take an edge with `|s_q-s_ell|=g`.

If the two signed lifts have the same orientation, then either

\[
s_q=h_q,\ s_\ell=h_\ell
\]

or

\[
s_q=-h_q-1,\ s_\ell=-h_\ell-1,
\]

and therefore

\[
\boxed{|h_q-h_\ell|=g.}
\tag{4.1}
\]

If the orientations are opposite, then

\[
\boxed{h_q+h_\ell+1=g.}
\tag{4.2}
\]

These are exactly the two shapes stated in Q8345.  The crucial qualifier is that the first zero is modulo `q` and the second is modulo `ell`.

---

## 5. Fixed-`g` graph geometry and exact ordered/unordered factors

For fixed `(m,g)`, put one vertex at each signed coordinate `s_q` of a retained label and join two vertices exactly when their difference in absolute value is `g` and all masks accept the pair.

Two distinct labels cannot occupy the same signed coordinate: if `s_q=s_ell`, then both shell primes divide `m-s_q`, contradicting the uniqueness argument of §3.  Hence the graph embeds in the integer line with distinct vertex positions.

Every vertex has degree at most two, because its only possible neighbours are at `s_q-g` and `s_q+g`.  Moreover a finite component cannot contain a cycle: at its smallest signed coordinate, two neighbours would both have to lie at the single point `s+g`.  Thus every nontrivial component is a path.

Let

- `E_g(m)` = number of **unordered** admissible edges;
- `V_g(m)` = number of vertices incident to at least one edge, exactly the Q8345 quantity.

If the nontrivial path components have vertex counts `v_1,...,v_C`, then

\[
V_g=\sum_i v_i,
\qquad
E_g=\sum_i(v_i-1)=V_g-C.
\tag{5.1}
\]

Since each `v_i>=2`,

\[
\boxed{\frac{V_g}{2}\le E_g<V_g.}
\tag{5.2}
\]

The number of **ordered** directed edge incidences is exactly `2E_g`; it is not `V_g` unless the graph is a matching.

Also every component lies in the signed interval of length `2X-1`, so

\[
\boxed{
v_i\le 1+\left\lfloor\frac{2X-1}{g}\right\rfloor.}
\tag{5.3}
\]

For `g>G=D^2`, this is at most

\[
1+\frac{2X}{G}
=X^{1/15-\eta+o(1)}(\log X)^2.
\]

This keeps a single chain short, but it does not bound the number of disjoint components; a matching can still have many vertices.

---

## 6. Actual Apéry falsification of same-characteristic gap transfer

The committed exact horizontal Mellin regression at `n=321` records genuine zeros at

\[
(p,r)=(179,142),(193,128),(211,110).
\]

Reflection gives folded rows

\[
(h_{179},h_{193},h_{211})=(36,64,100).
\]

Their signed physical lifts at `m=321` are

\[
321=2\cdot179-37,
\qquad
321=2\cdot193-65,
\qquad
321=2\cdot211-101.
\tag{6.1}
\]

Thus the first two endpoints have

\[
|s_{179}-s_{193}|=28,
\qquad
|h_{179}-h_{193}|=28,
\tag{6.2}
\]

and the partner integer is literally

\[
2\cdot179+28=386=2\cdot193.
\tag{6.3}
\]

Take `X=128`.  All three primes lie in `(X,2X]`.  For every
`0<eta<1/15`,

\[
\frac7{15}+\frac\eta2<\frac12,
\]

so

\[
D<\frac{\sqrt{128}}{\log128}<3,
\qquad D\le2,\qquad G\le4.
\tag{6.4}
\]

Hence `g=28` is far.  The two endpoint Kummer orders are

\[
\operatorname{ord}(\omega_{179}^{36})
=\frac{178}{\gcd(178,36)}=89,
\]

\[
\operatorname{ord}(\omega_{193}^{64})
=\frac{192}{\gcd(192,64)}=3,
\]

both greater than `D`.

The `179` zero set has exactly the reflected pair `{36,142}`: the committed million-prime scan has no prime with four zeros before `p=181`, while `142` is an exact zero and reflection forces `36`.  Therefore

\[
\boxed{179\nmid b_{64}.}
\tag{6.5}
\]

The owned verifier independently recomputes the entire row modulo 179 from the exact recurrence and checks (6.5).

The endpoints are also isolated at radius two.  No consecutive Apéry zeros are possible.  A zero at distance two from `h=64` modulo 193 would force `P(63)=0` or `P(65)=0`, but

\[
P(63)\equiv72\pmod{193},
\qquad
P(65)\equiv20\pmod{193}.
\]

The analogous residues at `h=36` modulo 179 are

\[
P(35)\equiv33\pmod{179},
\qquad
P(37)\equiv161\pmod{179}.
\]

Thus the displayed pair passes the far/high-order/isolation tests that are visible from Q8345.  Because the connector-visible tree does not contain the newest literal admissibility-mask implementation or `Y_X`, I do **not** claim that this finite pair is necessarily counted by the final `V_28(321)`.  It is used only for what a finite example can prove: it refutes the proposed universal algebraic implication

```text
mixed-characteristic physical gap g
    => same-characteristic Apéry return at gap g.
```

That implication is false on the actual Apéry sequence.

---

## 7. The fixed-gap continuant/resultant has the wrong characteristic

The banked gap theory says, schematically, that if a **single** prime `p` satisfies

\[
p\mid b_h,\qquad p\mid b_{h+g},
\tag{7.1}
\]

then `h` is a zero of the gap continuant `N_g` modulo `p`; combining two such same-characteristic conditions can be encoded by a fixed-gap resultant.  For every fixed gap, this leaves only finite prime support after the characteristic-zero resultant is known to be nonzero.

An actual GH1 edge supplies instead

\[
q\mid b_{h_q},\qquad \ell\mid b_{h_\ell},\qquad q\ne\ell.
\tag{7.2}
\]

Even when `|h_q-h_ell|=g`, (7.2) does not imply either

\[
q\mid b_{h_\ell}
\quad\text{or}\quad
\ell\mid b_{h_q}.
\tag{7.3}
\]

The linear physical relation

\[
c_qq-c_\ell\ell=\pm g
\]

only says, for example, `c_ell ell == +/- g (mod q)`.  It contains no operation that changes the reduction of the integer `b_{h_ell}` from characteristic `ell` to characteristic `q`.

Reflection, Lucas, the Apéry recurrence, and the gap continuant all act **inside a fixed residue field**.  Hasse--Franel likewise identifies `b_h` as a coefficient/Mellin coordinate in characteristic `p`; it does not identify zero coordinates in two different characteristics because their indices happen to be separated by the same integer `g`.

Therefore direct use of the fixed-gap resultant on a GH1 edge is not merely quantitatively weak; its hypothesis is absent.

---

## 8. Even a fictitious per-gap carrier has the wrong height exponent

This subsection grants more than is actually available, to isolate the quantitative wall.

Suppose that for every relevant `(m,g)` one somehow had a nonzero integer `C_{m,g}` such that every incident label divides it and

\[
\log|C_{m,g}|\le X^{o(1)}g^\alpha.
\tag{8.1}
\]

Then, since all incident labels exceed `X`,

\[
V_g(m)\log X\le\log|C_{m,g}|,
\]

so

\[
V_g(m)\le X^{o(1)}\frac{g^\alpha}{\log X}.
\tag{8.2}
\]

At the first far scale,

\[
D=X^{7/15+\eta/2+o(1)}(\log X)^{-1},
\]

\[
G=D^2=X^{14/15+\eta+o(1)}(\log X)^{-2}.
\tag{8.3}
\]

The GH1 subtraction threshold is

\[
T_X=\left(\frac14+o(1)\right)\frac{L^2}{X}
=X^{13/15+2\eta+o(1)}.
\tag{8.4}
\]

For (8.2) to be `o(T_X)` already at `g=G`, the power exponents require

\[
\alpha\left(\frac{14}{15}+\eta\right)
<\frac{13}{15}+2\eta.
\]

Equivalently,

\[
\boxed{
\alpha<\alpha_*(\eta)
=\frac{13+30\eta}{14+15\eta}<1
\quad(0<\eta<1/15).
}
\tag{8.5}
\]

So a linear carrier `alpha=1` is still too large by the power

\[
\frac1{15}-\eta>0.
\]

A quadratic carrier misses by

\[
2\left(\frac{14}{15}+\eta\right)
-\left(\frac{13}{15}+2\eta\right)=1.
\tag{8.6}
\]

Restoring the explicit `(log X)^{-2}` in `G` changes only logarithms: a roughly quadratic fixed-gap height gives a divisor-count scale whose ratio to `T_X` is `X^{1+o(1)}`.  Thus the fixed-gap resultant is not “almost enough” at `g=G`; it is a full power of `X` away.  At larger gaps the pointwise height problem only worsens.

This proves the requested comparison with the banked fixed-gap fact.  Finite prime support for each **fixed** `g` does not furnish the uniform estimate needed when `g` grows through `(G,2X)`.

---

## 9. High Kummer order does not fix the selected-character problem

For a prime `p`, multiplicative characters are indexed by exponents modulo `p-1`.  The number whose order divides exactly `d` is `phi(d)`.  Hence the number of characters of order at most `D` is

\[
\boxed{
\#\{\chi:\operatorname{ord}\chi\le D\}
=\sum_{\substack{d\mid p-1\\d\le D}}\varphi(d)
\le\sum_{d\le D}d
=\frac{D(D+1)}2.
}
\tag{9.1}
\]

For the Q8345 scale,

\[
D^2=G=X^{14/15+\eta+o(1)}(\log X)^{-2},
\]

and `p~X`, so

\[
\boxed{D^2/p=o(1).}
\tag{9.2}
\]

Thus the high-order sector has density `1-o(1)` among all characters.  Removing low-order characters does not turn an all-character theorem into a theorem for the one selected character attached to each physical label.

The repository's fixed-exponent Mellin audit checks the geometric side carefully: the rank-three K3 system has bounded complexity and the Kummer twist remains tame uniformly in its order.  Even granting the strongest common-group hypotheses for a varying-characteristic all-character equidistribution theorem, one deterministic high-order character per field contributes mass only `1/(p-1)` to that average.  It may be exceptional for every `p` without altering the limiting empirical measure.

There is a second independent obstruction.  The target is reduction at the defining prime itself.  Compatible-system large sieves work with auxiliary residual primes `ell != p` over a base field of characteristic `p`; setting the auxiliary prime equal to the moving base characteristic is not an allowed specialization.

Therefore

```text
Kummer order > D
```

is useful because it removes the explicitly countable low-order sector, but it does **not** validate a standard sheaf/large-sieve estimate for GH1.

---

## 10. The canonical affine second coordinate is not automatically present

The repository constructs a canonical inhomogeneous coordinate `kappa_r` and proves an exact affine-transfer statement.  On the ambient determinant coset, the equations

\[
b_r=0,\qquad\kappa_r=0
\]

form a genuine smooth codimension-two locus of exact density `1/(p(p+1))`.

This is important negatively: an ordinary retained Apéry zero `p|b_r` does **not** algebraically force `kappa_r=0`.  The unit-Casoratian identity makes `kappa` a legitimate transverse coordinate, not a hidden consequence of the first coordinate.  The Q7694 height audit further shows that the natural primitive rowwise numerator of this coordinate has no hidden subquadratic height compression.

Kummer order is a property of the character exponent; it does not add the equation `kappa_r=0`.  Any use of this coordinate in GH1 therefore needs a new defining-characteristic dispersion theorem.  It cannot be inserted for free as the second equation of a resultant.

---

# II. FIRST SURVIVING ARITHMETIC INTERFACE

## 11. The mixed-characteristic one-edge endpoint sieve

Define the fold map on signed integers

\[
H(s)=\begin{cases}
s,&s\ge0,\\-s-1,&s<0.\end{cases}
\tag{11.1}
\]

Fix a retained endpoint `q` with

\[
m=c_qq+s.
\]

For an orientation `delta in {+1,-1}` put

\[
\boxed{N_{q,\delta}(m,g)=c_qq-\delta g.}
\tag{11.2}
\]

If the partner signed row is `s+delta g`, then

\[
N_{q,\delta}(m,g)=m-(s+\delta g)=c_\ell\ell.
\tag{11.3}
\]

By §3 there is at most one shell-prime candidate `ell` dividing this integer.  Once that candidate is found, the genuinely arithmetic part of the edge is

\[
\boxed{
q\mid b_{H(s)},
\qquad
\ell\mid b_{H(s+\delta g)},
\qquad
\ell\mid N_{q,\delta}(m,g),
}
\tag{11.4}
\]

plus the literal retained-row and mutual-admissibility masks.

This is the first interface not killed by the audits above.  It has three useful properties.

1. It is **one-candidate**: `(m,q,g,delta)` determines at most one shell `ell`.
2. It is **actual-Apéry**: both zero tests are the real coefficients, not a synthetic model.
3. It is explicitly **mixed-characteristic**: the endpoint zero tests stay modulo their own primes while the physical relation couples the primes over `Z`.

What is missing is a theorem showing cancellation/nonalignment for these deterministic paired defining-characteristic tests as `(m,q,g)` vary.  None of the current fixed-gap continuants, all-character Mellin averages, auxiliary-prime large sieves, or one-row physical carriers provides it.

A carrier formulation of the missing input would be genuinely new.  For example, if one could construct nonzero integers `C_{m,g}` from the distinguished Apéry normalization such that

\[
\prod_{q\in\mathcal V_g(m)}q\mid C_{m,g}
\tag{11.5}
\]

and prove an aggregate **excess-height** estimate at the GH1 weight,

\[
\sum_m(k_D(m))_3
\sum_{G<g<2X}
\left(
\frac{\log|C_{m,g}|}{\log X}-T_X
\right)_+
\ll X^{o(1)}(1+X^2\Lambda^5),
\tag{11.6}
\]

then GH1 would follow immediately from

\[
V_g(m)\le\frac{\log|C_{m,g}|}{\log X}.
\]

But §8 says a pointwise construction whose height is merely a usual polynomial function of the gap must be extraordinarily small: an `X^{o(1)}g^alpha` bound needs `alpha<alpha_*(eta)<1` to win at the first far scale.  Thus (11.6) should be viewed as an **aggregate mixed-characteristic cancellation theorem**, not as an invitation to multiply the known fixed-gap resultants.

---

## 12. Why a second moment returns to a four-prime rectangle

Let

\[
T=T_X=\left(\frac14+o(1)\right)L^2/X.
\]

For all sufficiently large `X`, `T>=4`.  Fix `(m,g)` and abbreviate
`V=V_g(m)`, `E=E_g(m)`.

If `V<=T`, the GH1 positive part is zero.  If `V>T`, (5.2) gives

\[
E\ge V/2>T/2.
\]

Hence

\[
E-1>T/2-1\ge T/4,
\]

and therefore

\[
\binom E2=\frac{E(E-1)}2\ge\frac{ET}{8}.
\]

Since `V<=2E`,

\[
\boxed{
(V-T)_+\le V\le2E
\le\frac{16}{T}\binom E2.
}
\tag{12.1}
\]

This constant `16` uses **unordered** edges.  In ordered-edge notation `(E)_2=E(E-1)=2 binom(E,2)`, the corresponding constant is `8/T`.

Now classify an unordered pair of distinct unordered edges.

- If they share a vertex, their union has three labels: a length-two chain.
- If they do not share a vertex, their union has four labels: a rectangle/matching pair.

Because the gap graph has maximum degree two, the number of shared-edge pairs is exactly the number of degree-two vertices, hence at most `V<=2E`.  The number of disjoint edge pairs is at least

\[
\boxed{
\binom E2-V
\ge\frac{E(E-5)}2.
}
\tag{12.2}
\]

Thus once `E` is large, most edge pairs are disjoint four-label objects.

This answers the comparison requested in Q8345.  A naive Cauchy or edge-energy proof of GH1 does **not** reduce the arithmetic arity below the existing top-strip four-prime rectangle gateway.  It deterministically regenerates a far-physical four-label rectangle sector, with an additional three-label chain sector.  The disjoint sector is not the same literal mask as the previously studied top-strip rectangle, but it has the same four-characteristic obstruction and cannot be advertised as a smaller premise.

The genuinely earlier object is one edge, namely (11.4).  A theorem controlling that one-candidate mixed-characteristic endpoint correlation directly could beat the rectangle architecture.  No such theorem is currently banked.

---

# III. REFUTED ROUTES

## 13. REFUTED: feed every far physical edge into the fixed-gap resultant

**Claim refuted.**  The fixed-gap theorem requires two zeros modulo one prime.  GH1 supplies one zero modulo each of two different primes.  The actual `(179,193)` Apéry pair in §6 has folded gap 28 but `179 does not divide b_64`.  Therefore the required same-characteristic second zero is genuinely absent.

This refutation is scoped to the proposed transfer.  It does not refute using fixed-gap resultants elsewhere, for example to prune nonisolated rows inside a single characteristic.

## 14. REFUTED: finite prime support for each fixed gap is uniform enough

**Claim refuted by the exponent ledger.**  Even granting a common carrier, quadratic gap height has exponent deficit exactly `X^1` at `g=G`; linear height still misses by `X^{1/15-eta}`.  A pointwise gap carrier needs logarithmic-height exponent strictly below one, namely (8.5).  Therefore the existing roughly quadratic resultant is not a uniform solution for `G<g<2X`.

## 15. REFUTED: high Kummer order activates ordinary Mellin equidistribution

**Claim refuted as a quantifier implication.**  Low-order characters occupy only `o(p)` of the group.  Removing them leaves an all-character average essentially unchanged, while the target is still one selected high-order character per moving field.  Such a section can remain exceptional at every prime without violating all-character equidistribution.  Defining-prime divisibility also lies outside the auxiliary-`ell` large-sieve formalism.

## 16. REFUTED: triangular or affine local data supplies a free second equation

The triangular physical carrier is exactly one-label by (2.4) and by the Q8336 integral remainder identity.  The canonical `kappa` coordinate is genuinely transverse to `b`, not forced by it.  Neither can be counted as a second target equation without proving new arithmetic input.

---

# IV. CONDITIONAL / NEXT GATE

## 17. Minimal surviving theorem shape

The clean next theorem should address (11.4) directly.  One possible analytic statement is a **moving defining-characteristic endpoint dispersion theorem**: after imposing the exact quotient/sign/fold/isolation/order masks, prove that the deterministic candidate map

\[
(q,s)\longmapsto
N_{q,\delta}(m,g)=c_qq-\delta g
\]

hits shell primes whose own folded Apéry coefficient vanishes with enough cancellation to control the GH1 weighted excess.

Such a theorem must be uniform simultaneously in

```text
X,
m in the physical range,
g in (G,2X),
the quotient c_q,
the two orientation shapes,
and the moving defining characteristics q and ell.
```

A fixed-field sheaf theorem is insufficient.  An all-character theorem is insufficient.  A theorem at an auxiliary prime is insufficient.  A fixed-gap resultant is insufficient.  The premise has to see the actual deterministic partner map and both defining-prime coefficient-zero conditions.

A characteristic-zero carrier version is (11.5)–(11.6).  The height barrier (8.5) shows that if it is proved pointwise in `g`, it must have a genuinely sublinear gap-height exponent; more plausibly, the needed saving must occur only after the `(m,g)` aggregation, retaining signs rather than multiplying one huge resultant per gap.

I do not know a published theorem with those quantifiers, and the repository does not currently prove one.

---

# V. Status ledger

## PROVED

1. Signed fold identity `tau(s_p)=tau(h_p)` and exact triangular factorization (2.4).
2. For fixed `(m,q,g,orientation)`, at most one shell-prime partner exists.
3. The two folded shapes are exactly `|h_q-h_ell|=g` and `h_q+h_ell+1=g`.
4. For fixed `(m,g)`, the admissible signed gap graph is a path forest.
5. With unordered edge count `E` and incident-label count `V`, `V/2<=E<V`; ordered edge incidences equal `2E`.
6. Component length is at most `1+floor((2X-1)/g)`.
7. The actual Apéry data at `(X,m,q,ell,g)=(128,321,179,193,28)` refute same-characteristic gap transfer.
8. The low-Kummer-order character sector has cardinality at most `D(D+1)/2=O(G)=o(p)`.
9. A hypothetical pointwise carrier of height `X^{o(1)}g^alpha` needs `alpha<alpha_*(eta)<1`; quadratic height misses by one exact power of `X`.
10. The edge-pair inequality (12.1) has exact constant `16` for unordered edges (`8` for the ordered factorial `(E)_2`).
11. Edge pairs split exactly into three-label chains and disjoint four-label rectangles; high edge excess forces many of the latter.

## CONDITIONAL / OPEN

1. GH1 itself remains open in this report.
2. The first surviving arithmetic interface is the mixed-characteristic endpoint event (11.4).
3. A mixed-characteristic aggregate carrier/dispersion theorem of the form (11.5)–(11.6) would imply GH1.
4. A direct bilinear theorem for the deterministic partner map could also close GH1, but no existing theorem is verified to have the necessary defining-characteristic quantifiers.

## REFUTED

1. Directly applying the fixed-gap continuant/resultant to a mixed-characteristic far edge.
2. Treating roughly quadratic fixed-gap height as uniform over `G<g<2X`.
3. Treating even linear pointwise gap height as sufficient at the Q8345 exponents.
4. Inferring selected-label cancellation from the high-Kummer-order complement of an all-character average.
5. Treating the triangular physical carrier as an independent second equation.
6. Treating the canonical `kappa` coordinate as automatically zero at a retained Apéry hit.
7. Claiming that a generic second-moment attack is smaller than the existing four-prime rectangle gateway; it recreates a disjoint four-label sector.

---

## 18. Exact verifier

The owned standard-library verifier is

```text
problems/3.2/ORACLE_COMM/chatgpt_q8345_far_physical_one_label_verify.py
```

It uses only exact integer/rational arithmetic and checks:

- the actual Apéry zero rows modulo 179, 193, and 211 from the exact recurrence;
- the signed lifts at `m=321` and all triangular factorizations;
- the exact gap-28 partner integer `2*179+28=2*193`;
- `Z_179={36,142}` and `179 not dividing b_64`;
- the Kummer orders `(89,3,21)` and radius-two isolation of the relevant endpoints;
- shell-partner uniqueness on a deterministic finite regression grid;
- the path-forest and ordered/unordered edge identities by exhaustive small-set checks;
- exact low-order character counts against `sum_{d|p-1,d<=D} phi(d)`;
- the exponent ledger at a rational sample `eta=1/60`, including the exact identity `2*a_G-a_T=1`.

Its final banner is

```text
Q8345_FAR_PHYSICAL_OBSTRUCTION=PASS
```

The finite computations are used only for exact identity regression and falsification of stronger universal claims; no scan is used as an asymptotic proof.
