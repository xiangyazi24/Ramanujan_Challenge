ANSWER Q5224 ec39aebe

# Problem 3.2 — far–far active-pair attack after Q5217

## Verdict

I do **not** close the complete Q5217 far--far active-degree second moment

\[
\mathfrak D_2(T)=\sum_x w_xD_{x,T}^2
\]

at an exponent \(2\rho<145/93\).  The extra activity of both lower copies nevertheless yields two concrete advances that are not available for arbitrary Q4955 moving-hit pairs:

1. two chosen activity edges give a new exact nonzero four-prime determinant modulo the common distinguished characteristic \(p\); the attempted second CRT stage fails at an exact uncontrolled \(q\)-coordinate;
2. one **canonical** activity witness gives an unconditional three-zero cross-gap estimate.  It deletes a new part of the Q4955 far--far sector at the full canonical AD2 gate, without paying the possible \(T^{2/3}\) witness multiplicity.

The strongest new proved estimate is the following.  For each active lower copy choose one deterministic literal activity witness.  For an ordered pair of active lowers, let \(b\) be the nonwrapping distance from the chosen target of the first copy to the lower \(p\)-label of the second copy.  Then for dyadic \(B\ge2\),

\[
\boxed{
 \sum_xw_x\,#\{\text{ordered active pairs}:B\le b<2B\}
 \ll N^{o(1)}T B\min(B,T)M_0.}
\tag{V1}
\]

Consequently, for \(2\le K\le T\),

\[
\boxed{
 \mathfrak D_{2,\times\le K}(T)
 \ll N^{o(1)}TK^2M_0+N^{2/3+o(1)}M_0.}
\tag{V2}
\]

On the Q4955 far fixed-\(p\) sector the \(b=0\) term is actually

\[
\ll T^{10/9+o(1)}M_0.
\tag{V3}
\]

Thus at the **canonical** AD2 gate we may take

\[
\boxed{K=T^{26/93-\zeta}}
\]

and obtain

\[
\boxed{
 \mathfrak D_{2,\times\le K}^{\rm ff}(T)
 \ll T^{145/93-2\zeta+o(1)}M_0
      +T^{10/9+o(1)}M_0.}
\tag{V4}
\]

At the **ideal** AD2 line, whose second-moment threshold is \(23/15\), the correct choice is

\[
\boxed{K=T^{4/15-\zeta}},
\qquad
TK^2=T^{23/15-2\zeta},
\tag{V5}
\]

not \(T^{1/10}\).  This is the exact exponent corresponding to \(\rho<23/30\).

After Q4955 and (V4), the remaining canonical obstruction may therefore be restricted to the literal pair event

\[
\boxed{
 |d|>N^{8/15},\qquad
 g_p^\diamond(r_1,r_2)>N^{3/5},\qquad
 b_{1\to2}>T^{26/93-o(1)},\qquad
 b_{2\to1}>T^{26/93-o(1)}.}
\tag{V6}
\]

In this smaller event, either the two canonical activity targets coalesce, or the four physical lower/target moving primes are distinct and carry the nonzero determinant derived below.  No current Palm-weighted theorem controls either alternative with a fixed power saving.

The first exact failed implication is

\[
\boxed{
 r_i^+-r_i=\varepsilon_i h_i\pmod p
 \quad\not\Longrightarrow\quad
 m_i^+-m_i=\varepsilon_i h_i\pmod{pq}.}
\tag{V7}
\]

The exact correction is

\[
 m_i^+-m_i=\varepsilon_i h_i+pJ_i,
\]

and the integer \(J_i\) is an unconstrained \(q\)-coordinate because an unselected moving copy carries no new \(q\)-zero selector.  This is why the two activity witnesses do not generate a second independent CRT divisor.

---

# 0. Source lock

Only `xiangyazi24/Ramanujan_Challenge`, Problem 3.2, is used.

The connected canonical mathematical head remains

```text
c5d932b66ce5e4f1657b587215d290ae7a13018b
```

and I re-read through the GitHub connector:

- `problems/3.2/proof.tex`: `lem:no-consec`, `lem:gap-poly`, `lem:nonvanish`, restart/subinterval/column facts;
- `problems/3.2/energy_result.tex`: the actual three-zero count
  \[
  m_{a,b}(p)=\#\{u:N_a(u)=N_b(u+a)=0\}
  \]
  and its column/structured-energy interpretation;
- `problems/3.2/pairpalm_result.tex`: exact pair--Palm factorial hierarchy and the fact that higher Palm extension moments are genuinely new arithmetic;
- `drops/Q4760-47227b62.md`: determinant-counting reciprocal-fibre second moment and its explicit counting-measure/Palm-measure mismatch.

The connected `main` is stale for the current AT continuation.  I therefore use the already cold-audited same-project interfaces Q4955, Q5199 and Q5217 as authoritative for their later definitions, together with Q5155 for the reciprocal cubic and Q4703 for the horizontal Fourier obstruction.  No canonical source file is edited.

---

# 1. Literal two-copy Palm variables

Fix a selected Palm root

\[
x=(\mathbf b,\delta,\ell_0),
\qquad
\mathbf b=(\sigma,p,q,p'),
\qquad
P=pq,
\qquad
\Delta_x=\sigma\delta,
\]

with the unchanged nonnegative root atom

\[
\boxed{
 w_x=
 \Lambda_{\mathbf b}|u_{\ell_0}|^2
 a_{\mathbf b}(\delta,\ell_0)
 J_{\mathbf b}(\delta,\ell_0),
 \qquad M_0=\sum_xw_x.}
\tag{1.1}
\]

All four selected zero tests, phases, source masks, shell restrictions, distinctness, branch/orientation, high-centred, alias and near/reflected-near masks remain inside \(w_x\).  It appears exactly once.

Put

\[
 C_x=\langle\Delta_x\overline{p'}\rangle_P,
 \qquad
 a_x=\frac{p'C_x-\Delta_x}{P}.
\tag{1.2}
\]

For a literal unselected moving shell prime \(\ell\), define the canonical full-\(P\) completion

\[
 m_x(\ell)=\langle\Delta_x\overline{p'\ell}\rangle_P,
\tag{1.3}
\]

\[
 k_x(\ell)=\frac{\ell m_x(\ell)-C_x}{P},
 \qquad
 n_x(\ell)=a_x+p'k_x(\ell).
\tag{1.4}
\]

Then

\[
\boxed{\Delta_x=p'\ell m_x(\ell)-Pn_x(\ell).}
\tag{1.5}
\]

The two literal unselected actual-zero tests are

\[
 r_x(\ell):=m_x(\ell)\bmod p\in\mathcal Z_p,
\tag{1.6}
\]

\[
 s_x(\ell):=n_x(\ell)\bmod\ell
 =\langle-\Delta_x\overline P\rangle_\ell
 \in\mathcal Z_\ell.
\tag{1.7}
\]

There is no new \(q\)-selector.  The \(p'\)-reduction

\[
 n_x(\ell)\bmod p'
 =-\Delta_x\overline P\pmod{p'}
\tag{1.8}
\]

is the same inherited selected \(p'\)-zero residue for every unselected copy, not a new selector.

For two lower copies \(\ell_1,\ell_2\), abbreviate

\[
 m_i=m_x(\ell_i),\quad
 k_i=k_x(\ell_i),\quad
 n_i=n_x(\ell_i),\quad
 r_i=r_x(\ell_i),\quad
 s_i=s_x(\ell_i).
\]

The Q4955 reciprocal quotient displacement is

\[
\boxed{d=k_2-k_1,}
\tag{1.9}
\]

with the exact identities

\[
\boxed{\ell_2m_2-\ell_1m_1=Pd,}
\tag{1.10}
\]

\[
\boxed{n_2-n_1=p'd.}
\tag{1.11}
\]

The residual under attack already satisfies

\[
\boxed{|d|>N^{8/15}}
\tag{1.12}
\]

and

\[
\boxed{
 g_p^\diamond(r_1,r_2)
 :=\min(\|r_2-r_1\|_p,\|r_2+r_1+1\|_p)
 >N^{3/5}.}
\tag{1.13}
\]

Every moving-zero test and all lower-pair masks remain literal.

---

# 2. Attach the two activity witnesses

Fix one member of the finite ordinary/reflected orientation alphabet.  Write its sign as \(\varepsilon_i\in\{\pm1\}\).

Activity of lower copy \(i\) means that there exists a literal live target prime \(\lambda_i\) and an integer

\[
T\le h_i<2T
\]

such that the target's underlying actual distinguished-\(p\) label satisfies

\[
\boxed{
 r_i^+\equiv r_i+\varepsilon_i h_i\pmod p,}
\tag{2.1}
\]

with the chosen representatives nonwrapping.  Both \(r_i,r_i^+\) are actual \(p\)-zeros and the target also satisfies its own actual moving-characteristic zero test.

Define target completions

\[
 m_i^+=m_x(\lambda_i),\quad
 k_i^+=k_x(\lambda_i),\quad
 n_i^+=n_x(\lambda_i),
\]

and target edge quotient

\[
\boxed{e_i=k_i^+-k_i.}
\tag{2.2}
\]

Then exactly

\[
\boxed{\lambda_i m_i^+-\ell_i m_i=P e_i,}
\tag{2.3}
\]

\[
\boxed{n_i^+-n_i=p'e_i.}
\tag{2.4}
\]

The target-pair reciprocal quotient displacement is

\[
\boxed{d^+=d+e_2-e_1.}
\tag{2.5}
\]

No shortness of \(e_i\) or \(d^+\) follows from activity.

## 2.1 Witness multiplicity

For a fixed active lower \(p\)-label, possible targets lie in a consecutive interval of length \(O(T)\).  The canonical subinterval theorem gives

\[
\#\{\text{actual target }p\text{-labels}\}
\ll T^{2/3+o(1)}.
\]

A fixed target \(p\)-label has only \(O(1)\) moving-prime lifts.  Hence

\[
\boxed{
\#\{\text{literal witnesses of one active copy}\}
\ll T^{2/3+o(1)}.}
\tag{2.6}
\]

Thus attaching **all** witnesses to both lower copies can cost \(T^{4/3+o(1)}\).  The same loss appears in reverse: a fixed target can have \(T^{2/3+o(1)}\) actual predecessors in the backward interval.

To avoid this, for every active oriented copy choose deterministically the lexicographically least live pair \((h,\lambda)\).  I call it the **canonical activity witness**.  This is only a proof selection; it changes neither the active indicator nor the selected-root weight.

The important tradeoff is that the least-witness map is nonlinear, so it cannot be used as a new Fourier averaging variable.

---

# 3. Attack I: the two-witness determinant

Put

\[
A_x=C_x\bmod p.
\]

The reciprocal completion gives

\[
\ell_i r_i\equiv A_x\pmod p,
\qquad
\lambda_i r_i^+\equiv A_x\pmod p.
\]

Using (2.1),

\[
\boxed{
\varepsilon_i h_i\ell_i\lambda_i
\equiv A_x(\ell_i-\lambda_i)\pmod p.}
\tag{3.1}
\]

Equivalently, when the denominator is nonzero,

\[
\lambda_i
\equiv
\frac{A_x\ell_i}{A_x+\varepsilon_i h_i\ell_i}
\pmod p.
\tag{3.2}
\]

Eliminating \(A_x\) between the two activity edges gives

\[
\boxed{p\mid\Phi_{12},}
\tag{3.3}
\]

where

\[
\boxed{
\Phi_{12}
=
\varepsilon_1h_1\ell_1\lambda_1(\ell_2-\lambda_2)
-
\varepsilon_2h_2\ell_2\lambda_2(\ell_1-\lambda_1).}
\tag{3.4}
\]

This determinant does not exist for an arbitrary Q4955 lower pair without activity witnesses.

### Nonvanishing on the four-distinct mesoscopic sector

After the Q5217 large-\(T\) deletion, the unresolved range has \(T=o(N)\).  Assume the four physical moving primes \(\ell_1,\lambda_1,\ell_2,\lambda_2\) are pairwise distinct.  If \(\Phi_{12}=0\), then the prime \(\ell_1\) divides

\[
h_2\ell_2\lambda_2(\ell_1-\lambda_1).
\]

It divides neither \(h_2\) (because \(h_2<2T<\ell_1\)) nor either of the distinct primes \(\ell_2,\lambda_2\).  Thus \(\ell_1\mid(\ell_1-\lambda_1)\), hence \(\ell_1\mid\lambda_1\), forcing \(\ell_1=\lambda_1\), contradiction.  Therefore

\[
\boxed{\Phi_{12}\ne0.}
\tag{3.5}
\]

Also \(|\Phi_{12}|\ll TN^3\).  This gives only an **absolute** divisor capacity: for a fixed full six-variable tuple, only \(O(1)\) shell primes can divide \(\Phi_{12}\).  It does not give a relative selected-Palm saving because the mass may concentrate on one such \(p\).

If \(h_1=h_2=h\) and the orientations agree, then

\[
\Phi_{12}
=-h\left((\ell_1-\lambda_1)\ell_2\lambda_2
        -(\ell_2-\lambda_2)\ell_1\lambda_1\right),
\]

which is exactly the Q5155 reciprocal four-prime cubic after relabelling.  Thus the equal-gap specialization gives no second Pluecker carrier.

---

# 4. Exact failure of the second CRT stage

Equation (2.1) controls the target-minus-lower completed row only modulo \(p\).  Hence there is an integer \(J_i\) with

\[
\boxed{
 m_i^+-m_i=\varepsilon_i h_i+pJ_i.}
\tag{4.1}
\]

Insert this into (2.3):

\[
(\lambda_i-\ell_i)m_i
+\varepsilon_i\lambda_i h_i
+p\lambda_iJ_i
=pq e_i.
\]

The first two terms are divisible by \(p\) by (3.1).  Thus

\[
B_i:=
\frac{(\lambda_i-\ell_i)m_i
      +\varepsilon_i\lambda_i h_i}{p}
\in\mathbb Z
\]

and the exact \(q\)-stage identity is

\[
\boxed{B_i+\lambda_iJ_i=q e_i.}
\tag{4.2}
\]

The unselected copy imposes no actual \(q\)-zero condition on \(m_i\bmod q\) or \(m_i^+\bmod q\).  Therefore no current theorem forces \(J_i=0\), makes \(J_i\) short, or correlates \(J_1,J_2\) with a power saving.

This is the first exact failed implication:

\[
\boxed{
 r_i^+-r_i=\varepsilon_i h_i\pmod p
 \not\Longrightarrow
 m_i^+-m_i=\varepsilon_i h_i\pmod{pq}.}
\tag{4.3}
\]

The two activity witnesses therefore produce one common-\(p\) four-zero geometry, not a second independent \(q\)-geometry.

There is likewise no legal six-zero same-characteristic resultant:

- \(p\) sees the four actual labels \(r_1,r_1^+,r_2,r_2^+\);
- \(p'\) sees one repeated inherited residue;
- \(q\) sees no new unselected zero;
- \(\ell_1,\lambda_1,\ell_2,\lambda_2\) each see one zero in a different characteristic.

Cross-characteristic resultants would therefore be illegal.

---

# 5. Attack II: one canonical witness creates a three-zero codegree

This attack is genuinely different from the determinant argument and gives a proved deletion.

Fix one root \(x\), one orientation chart, and the canonical witness of every active lower copy.  Take an ordered pair of distinct active lowers and order their actual oriented \(p\)-labels as

\[
0\le r<u\le p-1.
\]

Let the chosen witness of the left copy be

\[
T\le h<2T,
\qquad r+h\in\mathcal Z_p,
\]

in the current nonwrapping chart.  Define the cross gap

\[
\boxed{b=|u-(r+h)|.}
\tag{5.1}
\]

Assume first \(B\le b<2B\), \(B\ge2\).

### Interlaced case

If \(r<u<r+h\), set

\[
a=u-r=h-b.
\]

Then

\[
r,\quad r+a=u,\quad r+a+b=r+h
\]

are three actual nonwrapping \(p\)-zeros, so

\[
N_a(r)=N_b(r+a)=0\pmod p.
\tag{5.2}
\]

### Noninterlaced case

If \(r<r+h<u\), then

\[
r,\quad r+h,\quad r+h+b=u
\]

are three actual zeros, so

\[
N_h(r)=N_b(r+h)=0\pmod p.
\tag{5.3}
\]

For fixed adjacent gap parameters, `lem:nonvanish` and \(\deg N_j=3(j-1)\) give directly

\[
\boxed{
 m_{a,b}(p)\le3(\min(a,b)-1),}
\tag{5.4}
\]

or the same bound with \((a,b)=(h,b)\).  This is just the smaller of the two polynomial degree bounds; it does not require the global \(H\le\sqrt p\) energy estimate.

There are \(O(T)\) choices of \(h\), \(O(B)\) choices of \(b\), and at most \(O(\min(B,T))\) starting residues for each pair.  Each of the three fixed \(p\)-labels has only \(O(1)\) literal moving-prime lifts in a fixed source box.  The canonical-witness condition and all remaining source/moving-zero/far-pair masks are deletions.

Therefore, root by root,

\[
\boxed{
\#\{\text{ordered active pairs}:B\le b<2B\}
\ll N^{o(1)}TB\min(B,T).}
\tag{5.5}
\]

Multiplication by the unchanged \(w_x\) proves (V1).

### Why there is no \(T^{2/3}\) witness loss

Every lower pair uses exactly its selected canonical target.  We never sum over all possible targets.  Once the three \(p\)-labels are fixed, all lower/target prime fibres are \(O(1)\).  Alternate activity witnesses are irrelevant to this upper bound.

---

# 6. Cross-near deletion and exponents

For \(2\le K\le T\), summing (5.4) directly over \(2\le b\le K\) yields

\[
\boxed{
\mathfrak D_{2,\times[2,K]}(T)
\ll N^{o(1)}TK^2M_0.}
\tag{6.1}
\]

The distance-one case is empty by `lem:no-consec`.

For \(b=0\), the other lower copy has the same oriented \(p\)-label as the chosen target.  A fixed \(p\)-label has \(O(1)\) moving-prime lifts, so this contributes rootwise \(O(D_{x,T})\) and globally

\[
\ll N^{2/3+o(1)}M_0.
\tag{6.2}
\]

In the Q4955 far fixed-\(p\) sector, \(b=0\) implies

\[
N^{3/5}<g_p^\diamond(r,u)\le h<2T,
\]

so \(T\gg N^{3/5}\), and hence

\[
N^{2/3}\ll T^{10/9}.
\tag{6.3}
\]

This proves (V2)--(V3).

## 6.1 Canonical gate

Take

\[
K=T^{26/93-\zeta}.
\]

Then

\[
TK^2
=T^{1+52/93-2\zeta}
=T^{145/93-2\zeta}.
\]

Thus the whole directional cross-near part satisfies the canonical AD2 target with a fixed margin.  Taking the union of both directions costs only a constant factor.

## 6.2 Ideal AD2 line

The ideal Q5217 second-moment target is \(T^{23/15-o(1)}M_0\).  Therefore set

\[
K=T^{4/15-\zeta}.
\]

Then

\[
TK^2=T^{1+8/15-2\zeta}
=T^{23/15-2\zeta}.
\]

Hence this cross-near sub-sector also satisfies AD2 with

\[
\rho=\frac{23}{30}-\zeta.
\]

---

# 7. What remains after the new deletion

At the canonical gate, fix a small \(\zeta>0\) and choose the canonical witness of each active copy.  We may now restrict to

\[
|d|>N^{8/15},
\qquad
 g_p^\diamond(r_1,r_2)>N^{3/5},
\tag{7.1}
\]

and, in both directional charts,

\[
\boxed{
\operatorname{dist}(r_1^+,r_2)>T^{26/93-\zeta},
\qquad
\operatorname{dist}(r_2^+,r_1)>T^{26/93-\zeta}.}
\tag{7.2}
\]

This is strictly smaller than the Q5217 far--far AD2 event.

There are two honest subcases.

### 7.1 Coalescent targets

The canonical targets may satisfy

\[
\lambda_1=\lambda_2.
\]

Then two distinct lower zeros return to one common target.  In a common orientation the lower gap is \(|h_1-h_2|\).  A fixed target can have as many as \(T^{2/3+o(1)}\) actual predecessors in the relevant backward interval, so the current subinterval theorem permits a substantial V-shaped star.  The existing \(T^{8/3}\) collision-energy bound, where its \(H\le\sqrt p\) range applies, is still far above the \(T^{145/93}\) AD2 gate.  Thus coalescence is not automatically harmless.

### 7.2 Four distinct physical moving primes

If the targets do not coalesce, the generic four-distinct subcase carries

\[
p\mid\Phi_{12}\ne0.
\]

The exact q-lift defect (4.2) remains, and no present Palm-relative divisor theorem converts this determinant into a power saving.

---

# 8. Attack III: disjoint-gap/phase routes and their exact failure

## 8.1 Summing all activity witnesses

Replacing the active indicator by the positive witness count exposes \((h_1,\lambda_1,h_2,\lambda_2)\), but can cost

\[
T^{2/3+o(1)}\times T^{2/3+o(1)}
=T^{4/3+o(1)}.
\]

Reindexing by targets does not remove this: each target can have the same \(T^{2/3+o(1)}\) predecessor multiplicity.  Thus a naive four-/six-zero energy over all witness choices loses too much before any gate saving can occur.

The canonical-witness choice avoids that loss, but then there is no linear witness sum left on which to apply Cauchy or a large sieve.

## 8.2 Moving-characteristic Fourier expansion

For a target characteristic \(\lambda\), the exact actual-zero identity is

\[
\mathbf1_{\mathcal Z_\lambda}(s_\lambda)-\frac{Z(\lambda)}\lambda
=
\frac1\lambda\sum_{a\ne0}F_\lambda(a)e_\lambda(-as_\lambda).
\]

Q4703 audits this trace function.  On a **successful** moving-zero hit, the summand indexed by the actual zero \(u=s_\lambda\) is the trivial Artin--Schreier trace, whose complete frequency sum is \(\lambda-1\).  Two successful activity targets therefore create two positive trivial components, not a new local orthogonality mechanism.

A useful phase theorem would have to be genuinely horizontal across changing characteristics and uniform under the selected-Palm conditioning.  No such theorem is currently banked.

---

# 9. Current-interface saturation after the cross-near deletion

This is an interface packet, **not** an actual Apéry or prime construction.

Take

\[
N\asymp T^2
\]

and choose

\[
\frac{145}{186}<\kappa<\frac45.
\]

This interval is nonempty since

\[
\frac45-\frac{145}{186}=\frac{19}{930}>0.
\]

Concentrate the formal Palm mass on one legal root and take

\[
D\asymp T^\kappa.
\]

Place the lower distinguished-\(p\) labels at spacing

\[
R=N^{3/5}=T^{6/5}.
\]

The packet fits in the \(p\)-line because \(D R\ll T^2\) when \(\kappa<4/5\).  Give each lower one target at a common gap \(h_0\asymp T\).  Since \(R\gg T\), every lower-to-foreign-target cross gap is \(\asymp R\), far above \(T^{26/93}\).  The fixed-gap root cap permits the \(D<T\) pairs at \(h_0\); a length-\(T\) interval sees only \(O(1)\) packet zeros; and the total zero count \(O(D)\) is below \(N^{2/3}=T^{4/3}\).

At the level of the current reciprocal-fibre interfaces there is also room for \(D\) lower quotient coordinates separated by

\[
N^{8/15}=T^{16/15},
\]

because \(D T^{16/15}=o(T^2)\) for \(\kappa<14/15\).

Thus the banked upper-bound interfaces still admit a formal cross-deep far--far packet with

\[
D^2\asymp T^{2\kappa}>T^{145/93}.
\]

This does **not** claim that the actual Apéry sequence realizes the packet, that all required shell primes exist in prescribed reciprocal fibres, or that the original four-copy phases align.  It shows only that the currently proved marginal/codegree interfaces do not imply the missing cross-deep AD2 bound.

---

# 10. Smallest remaining arithmetic statement

Fix \(\zeta>0\) and put

\[
K_\zeta=T^{26/93-\zeta}.
\]

Let

\[
\mathfrak D_{2,\rm ff}^{\rm crossdeep}(T;K_\zeta)
\]

denote the literal Palm-weighted ordered active-pair mass satisfying Q4955's far conditions (7.1) and both canonical-witness cross-deep conditions (7.2).  All actual moving-zero tests, reciprocal labels, source masks, distinctness conditions and the single selected-root weight are retained.

Sections 5--6 prove unconditionally that the complementary cross-near sector is already subcritical.  Therefore the smallest clean remaining input is

\[
\boxed{
\mathfrak D_{2,\rm ff}^{\rm crossdeep}(T;K_\zeta)
\ll T^{145/93-\delta+o(1)}M_0
}
\tag{FFAD2-crossdeep}
\]

for some fixed \(\delta>0\).

This is strictly narrower than Q5217's AD2 target.  On its generic four-distinct part it additionally carries the nonzero determinant \(p\mid\Phi_{12}\); its separate collision geometry is canonical-target coalescence.

For the ideal line, replace \(K_\zeta\) by

\[
T^{4/15-\zeta}
\]

and \(145/93\) by \(23/15\).

No current canonical theorem proves either remaining statement.

---

# 11. Insertion-ready LaTeX — proved statements only

```latex
\begin{lemma}[Two active witnesses and the reciprocal $p$-determinant]
\label{lem:at-active-two-witness-determinant}
Fix a selected Palm root
\[
 x=(\mathbf b,\delta,\ell_0),\qquad
 \mathbf b=(\sigma,p,q,p'),\qquad P=pq,
\]
and put
\[
 C_x=\langle\Delta_x\overline{p'}\rangle_P,
 \qquad A_x=C_x\pmod p.
\]
For $i=1,2$, let $\ell_i$ be a literal live moving-prime copy and
let $(h_i,\lambda_i)$ be a literal activity witness.  In a fixed
ordinary/reflected orientation write
\[
 r_i^+\equiv r_i+\varepsilon_i h_i\pmod p,
 \qquad \varepsilon_i\in\{\pm1\},
 \qquad T\le h_i<2T.
\]
Then
\[
 \varepsilon_i h_i\ell_i\lambda_i
 \equiv A_x(\ell_i-\lambda_i)\pmod p,
\]
and consequently
\[
 p\mid
 \Phi_{12}:=
 \varepsilon_1h_1\ell_1\lambda_1(\ell_2-\lambda_2)
 -\varepsilon_2h_2\ell_2\lambda_2(\ell_1-\lambda_1).
\]
If moreover the four physical moving primes are pairwise distinct and
$2T$ is smaller than every moving-prime shell element, then
$\Phi_{12}\ne0$.

Writing
\[
 m_i^+-m_i=\varepsilon_i h_i+pJ_i,
 \qquad e_i=k_i^+-k_i,
\]
one has the exact integer identity
\[
 \frac{(\lambda_i-\ell_i)m_i
       +\varepsilon_i\lambda_i h_i}{p}
 +\lambda_iJ_i
 =q e_i.
\]
\end{lemma}

\begin{proof}
The reciprocal completion gives
$\ell_i r_i\equiv A_x\pmod p$ and
$\lambda_i r_i^+\equiv A_x\pmod p$.  Substituting
$r_i^+=r_i+\varepsilon_i h_i$ and eliminating $r_i$ proves the first
congruence.  Eliminating $A_x$ between the two copies gives the
$\Phi_{12}$ divisibility.

If $\Phi_{12}=0$, the prime $\ell_1$ divides
$h_2\ell_2\lambda_2(\ell_1-\lambda_1)$.  Under the stated size and
four-distinct hypotheses it divides none of the first three factors,
so $\ell_1\mid(\ell_1-\lambda_1)$, forcing
$\ell_1=\lambda_1$, a contradiction.

Finally
$\lambda_i m_i^+-\ell_i m_i=pq e_i$.  Substituting
$m_i^+-m_i=\varepsilon_i h_i+pJ_i$ gives the displayed integer
identity after division by $p$.
\end{proof}

\begin{lemma}[One-witness cross-gap bound]
\label{lem:at-active-one-witness-crossgap}
Fix one selected Palm root, one literal source box and one member of the
finite orientation alphabet.  For every active lower moving-prime copy
choose deterministically one literal activity witness, for example the
lexicographically least pair $(h,\lambda)$.

Let $\mathscr P_{x,T}(B)$ be any subfamily of ordered pairs of distinct
active lower copies with the following property.  After ordering their
oriented actual $p$-labels as
\[
 0\le r<u\le p-1,
\]
use the chosen witness of the left copy,
\[
 T\le h<2T,\qquad r+h\in\mathcal Z_p,
\]
and assume the three labels are nonwrapping and distinct and
\[
 B\le |u-(r+h)|<2B,
 \qquad B\ge2.
\]
All additional literal source, moving-zero, reciprocal, prime-exclusion
and pair masks may be retained in $\mathscr P_{x,T}(B)$.
Then
\[
 \#\mathscr P_{x,T}(B)
 \ll N^{o(1)}T B\min(B,T).
\]
Consequently
\[
 \sum_xw_x\#\mathscr P_{x,T}(B)
 \ll N^{o(1)}T B\min(B,T)M_0.
\]
\end{lemma}

\begin{proof}
If $r<u<r+h$, put $b=r+h-u$ and $a=u-r=h-b$.
Then $r,r+a,r+a+b$ are three actual nonwrapping $p$-zeros, so
\[
 N_a(r)=N_b(r+a)=0\pmod p.
\]
If $r+h<u$, put $b=u-r-h$; then
$r,r+h,r+h+b$ are three actual zeros and
\[
 N_h(r)=N_b(r+h)=0\pmod p.
\]
In either case Lemma~\ref{lem:nonvanish} and
$\deg N_j=3(j-1)$ give at most
$3(\min(B,T)+O(1))$ possible starting residues for each fixed pair of
gap parameters.  There are $O(T)$ choices of $h$ and $O(B)$ choices of
$b$.

A fixed oriented $p$-label has only $O(1)$ moving-prime lifts in the
fixed shell because the reciprocal congruence fixes one residue class
modulo $p\asymp N$.  Thus restoring the two lower primes and the chosen
target costs only $N^{o(1)}$.  Every remaining literal mask is a
deletion.  This proves the result.
\end{proof}

\begin{corollary}[Cross-near deletion at the active-degree gates]
\label{cor:at-active-crossnear-gate}
In the Q4955 far--far lower-pair sector, let
$\mathfrak D_{2,\times\le K}^{\rm ff}(T)$ denote the Palm-weighted
ordered pair mass for which, in at least one directional orientation
chart, the chosen activity target lies at nonwrapping integer distance
at most $K$ from the other lower $p$-label.  Then, for $2\le K\le T$,
\[
 \mathfrak D_{2,\times\le K}^{\rm ff}(T)
 \ll N^{o(1)}TK^2M_0+N^{2/3+o(1)}M_0.
\]
The distance-one part is empty.  On the distance-zero part the far
fixed-$p$ condition forces $N^{3/5}<2T$, and hence
\[
 N^{2/3}\ll T^{10/9}.
\]
In particular, for every fixed $\zeta>0$, taking
\[
 K=T^{26/93-\zeta}
\]
gives
\[
 \mathfrak D_{2,\times\le K}^{\rm ff}(T)
 \ll
 T^{145/93-2\zeta+o(1)}M_0
 +T^{10/9+o(1)}M_0.
\]
Taking instead $K=T^{4/15-\zeta}$ gives
\[
 \mathfrak D_{2,\times\le K}^{\rm ff}(T)
 \ll
 T^{23/15-2\zeta+o(1)}M_0
 +T^{10/9+o(1)}M_0.
\]
\end{corollary}
```

---

# 12. Exponent and obstruction ledger

```text
Q4955 lower-pair residual:
  |d| > N^(8/15)
  g_p^diamond > N^(3/5)

activity witness per lower:
  actual p-target at T <= h < 2T
  actual moving-characteristic target zero
  all source/reciprocal masks retained
  possible witnesses per lower: <= T^(2/3+o(1))
  two-witness expansion loss: up to T^(4/3+o(1))

new one-edge congruence:
  eps h ell lambda = A(ell-lambda) mod p

new two-edge determinant:
  p | Phi_12
  Phi_12 != 0 on four-distinct mesoscopic sector
  equal-h specialization = Q5155 reciprocal cubic

second CRT stage:
  m^+ - m = eps h + p J
  B + lambda J = q e
  J controlled by current q-zero selector: NO

legal same-characteristic new zeros:
  p: four endpoints
  p': one repeated inherited residue
  q: no unselected zero
  moving characteristics: one zero in each changing field
  legal six-zero resultant: NO

proved one-witness cross-gap:
  b ~ B:             T B min(B,T) M0
  2 <= b <= K <= T: T K^2 M0
  b = 1:             empty
  b = 0 in far-p:    T^(10/9+o(1)) M0

canonical AD2 gate:
  target exponent: 145/93
  K = T^(26/93-zeta)
  cross-near exponent: 145/93 - 2 zeta

ideal AD2 line:
  target exponent: 23/15
  K = T^(4/15-zeta)
  cross-near exponent: 23/15 - 2 zeta

remaining cross-deep event:
  Q4955 far d
  Q4955 far p-diamond gap
  both canonical lower-to-foreign-target gaps > gate cutoff
  target coalescence OR four-distinct nonzero Phi_12

phase route:
  all-witness expansion loses T^(4/3)
  canonical witness is nonlinear
  successful moving zero has trivial trace summand
  no banked horizontal operator theorem

full far-far AD2 with rho < 145/186: NOT PROVED
new restricted far-far deletion: PROVED
smallest remaining premise: FFAD2-crossdeep
```

## Final conclusion

The activity information is genuinely useful but does not finish AD2.  It produces a new four-prime determinant in the common characteristic \(p\), yet the exact full-\(P\) identity exposes an uncontrolled \(q\)-coordinate \(J_i\), so there is no second CRT modulus and no legal six-zero same-characteristic resultant.

The strongest unconditional gain is instead combinatorial: one canonical activity edge turns the second lower label into a genuine three-zero configuration.  This proves the cross-gap bound (V1), deleting all canonical far--far pairs for which either chosen activity target approaches the other lower label within \(T^{26/93-\zeta}\) at the \(38/31\) gate, or within \(T^{4/15-\zeta}\) at the ideal AD2 line.

What remains is the strictly smaller cross-deep far--far Palm pair problem `(FFAD2-crossdeep)`.  On its generic part one has the nonzero determinant \(p\mid\Phi_{12}\); on its exceptional part the two activity witnesses coalesce at one target.  Neither is presently controlled with a power saving by the canonical Problem 3.2 interfaces.